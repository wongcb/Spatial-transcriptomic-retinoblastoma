########Spatial transcriptomic profiling of human retinoblastoma############
########created by Luozixian Wang and Raymond Wong##########################
########CERA, UNIMELB, 05/09/2025###########################################
##load the packages
library(AUCell)
library(BiocManager)
library(bmixture)
library(CARD)
library(CellChat)
library(clusterProfiler)
library(ComplexHeatmap)
library(cowplot)
library(DDRTree)
library(DESeq2)
library(devtools)
library(doParallel)
library(dplyr)
library(EnhancedVolcano)
library(enrichplot)
library(foreach)
library(future)
library(garnett)
library(ggalluvial)
library(ggplot2)
library(ggpubr)
library(GO.db)
library(grid)
library(gridGraphics)
library(gridExtra)
library(harmony)
library(infercnv)
library(kableExtra)
library(KernSmooth)
library(Matrix)
library(MuSiC)
library(NMF)
library(NNLM)
library(patchwork)
library(pathview)
library(promises)
library(RColorBrewer)
library(rdist)
library(remotes)
library(reticulate)
library(SCDC)
library(SCENIC)
library(SCopeLoomR)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(SeuratObject)
library(SeuratWrappers)
library(SingleCellExperiment)
library(spdep)
options(stringsAsFactors = FALSE)
library(STutility)
library(swne)
library(tibble)
library(tidyverse)
library(zellkonverter)

#####QC and processing#####
## filtered spots (under tissue)
infoTable <- read.csv('infoTable.csv')
se <- InputFromTable(infoTable)

##mito + ribo content
# Collect all genes coded on the mitochondrial genome
mt.genes <- grep(pattern = "^MT-", x = rownames(se), value = TRUE)
se$percent.mito <- (Matrix::colSums(se@assays$RNA@counts[mt.genes, ])/Matrix::colSums(se@assays$RNA@counts))*100

# Collect all genes coding for ribosomal proteins
rp.genes <- grep(pattern = "^RPL|^RPS", x = rownames(se), value = TRUE)
se$percent.ribo <- (Matrix::colSums(se@assays$RNA@counts[rp.genes, ])/Matrix::colSums(se@assays$RNA@counts))*100

# Keep spots with more than 20 unique genes and less than 30% mitochondrial content
se.QC <- SubsetSTData(se, expression = nFeature_RNA > 20 & percent.mito < 30)
cat("Spots removed: ", ncol(se) - ncol(se.QC), "\n")

UMIplot.se.QC <- ST.FeaturePlot(se.QC, features = "nFeature_RNA", cols = c("lightgray", "mistyrose", "red", "dark red", "black"), pt.size = 2, ncol = 2)
UMIplot.se <- ST.FeaturePlot(se, features = "nFeature_RNA", cols = c("lightgray", "mistyrose", "red", "dark red", "black"), pt.size = 2, ncol = 2)
custom_theme <- theme(legend.position = c(0.45, 0.8), # Move color legend to top
                      legend.direction = "horizontal", # Flip legend
                      legend.text = element_text(angle = 30, hjust = 1), # rotate legend axis text
                      strip.text = element_blank(), # remove strip text
                      plot.title = element_blank(), # remove plot title
                      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm")) # remove plot margins

UMIplot.se.QC & theme_bw()

## loading images
se.QC <- LoadImages(se.QC, time.resolve = FALSE, verbose = TRUE)

## manual annotation 
se.QC <- ManualAnnotation(se.QC)

## QC check counts metrics 
QC.annotated <- VlnPlot(se.QC, features = c("nFeature_RNA", "nCount_RNA"), group.by = "labels")
saveRDS(se.QC, 'se.QC.RDS')

##subset out the unnannotated
Idents(se.QC) <- 'labels'
se.QC.anno<- SubsetSTData(se.QC, idents = c('RB 1A', 'RB 1B', 'RB 1C', 'RB 1D', 'RB 2A', 'RB 2B'))

## SCTransform => update to V2
# Add a section column to your meta.data
se.QC.anno$section <- paste0("section_", GetStaffli(se.QC.anno)[[, "sample", drop = T]])
table(se.QC.anno$section)
se.QC.anno <- SCTransform(se.QC.anno, vars.to.regress = "section", vst.flavor = "v2")

##QC SCTransform: rationale of approach normalisation 
# Get raw count data 
umi_data <- GetAssayData(object = se.QC.anno, slot = "counts", assay = "RNA")
dim(umi_data)

# Calculate gene attributes
gene_attr <- data.frame(mean = rowMeans(umi_data),
                        detection_rate = rowMeans(umi_data > 0),
                        var = apply(umi_data, 1, var), 
                        row.names = rownames(umi_data)) %>%
  mutate(log_mean = log10(mean), log_var = log10(var))

# Obtain spot attributes from Seurat meta.data slot
spot_attr <- se.QC.anno[[c("nFeature_RNA", "nCount_RNA")]]

QC1 <- ggplot(gene_attr, aes(log_mean, log_var)) + 
  geom_point(alpha = 0.3, shape = 16) + 
  geom_density_2d(size = 0.3) +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  ggtitle("Mean-variance relationship")

# add the expected detection rate under Poisson model
x = seq(from = -2, to = 2, length.out = 1000)
poisson_model <- data.frame(log_mean = x, detection_rate = 1 - dpois(0, lambda = 10^x))
QC2 <- ggplot(gene_attr, aes(log_mean, detection_rate)) + 
  geom_point(alpha = 0.3, shape = 16) + 
  geom_line(data = poisson_model, color='red') +
  ggtitle("Mean-detection-rate relationship")

pdf("SCTransformQC.se.QC.anno.pdf",width=10,height=5)
QC1 - QC2
dev.off()

# https://ludvigla.github.io/STUtility_web_site/Maize_data.html
# simple workflow for Dimensionality reduction (PCA), UMAP embedding, Clustering
se.QC.anno <- se.QC.anno %>% 
  RunPCA() %>%
  RunUMAP(reduction = "pca", dims = 1:30)

se.QC.anno <- se.QC.anno %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(reduction = "pca", dims = 1:30)
se.QC.anno$seurat_clusters_pca <- se.QC.anno$seurat_clusters

#Clustering  
se.QC.anno.UMAPp1 <- DimPlot(se.QC.anno, group.by = "labels", reduction = "umap")
se.QC.anno.UMAPp2 <- DimPlot(se.QC.anno, group.by = "seurat_clusters_pca", label = TRUE, label.size = 8, reduction = "umap")
se.QC.anno.UMAPp1 - se.QC.anno.UMAPp2

se.QC.anno.clusterSTp1 <- ST.FeaturePlot(se.QC.anno, features = "seurat_clusters_pca", indices = 1, split.labels = T, pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())
se.QC.anno.clusterSTp2 <- ST.FeaturePlot(se.QC.anno, features = "seurat_clusters_pca", indices = 2, split.labels = T, pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())
se.QC.anno.clusterSTp3 <- ST.FeaturePlot(se.QC.anno, features = "seurat_clusters_pca", indices = 3, split.labels = T, pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())
se.QC.anno.clusterSTp4 <- ST.FeaturePlot(se.QC.anno, features = "seurat_clusters_pca", indices = 4, split.labels = T, pt.size = 2) & theme(plot.title = element_blank(), strip.text = element_blank())
cowplot::plot_grid(se.QC.anno.clusterSTp1, se.QC.anno.clusterSTp2, se.QC.anno.clusterSTp3, se.QC.anno.clusterSTp4, ncol = 4)

#There are also some discrepancies in the spatial distribution of clusters in the different sections.
# use harmony for integration
se.QC.anno <- RunHarmony(se.QC.anno, group.by.vars = "labels", reduction = "pca", dims.use = 1:30, assay.use = "SCT", verbose = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony") %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters()
se.QC.anno$seurat_clusters_harmony <- se.QC.anno$seurat_clusters
saveRDS(se.QC.anno, 'se.QC.anno.RDS')


#####Downstream analysis#####
#####script is organised in figure generation + module analysis structure#####
#load the RDS file
SP <- readRDS("se.QC.anno.RDS")
new.cluster.ids <- c("C0", "C1", "C2", "C3", "C4", 
                     "C5", "C6", "C7", "C8")
Idents(SP) <- SP$seurat_clusters_harmony
names(new.cluster.ids) <- levels(SP)
SP <- RenameIdents(SP, new.cluster.ids)
STeye <- readRDS("STeye.RDS")

##plot the distribution of clusters on spatial coordinates##
##figure 1B##
SP <- MaskImages(object = SP)
ImagePlot(SP, method = "raster", type = "masked", indices = 5)
ImagePlot(SP, method = "raster", type = "masked", indices = 6)
ST.FeaturePlot(object = SP, features = "seurat_clusters_harmony", pt.size = 3) + 
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 8), legend.position = "bottom")
STeye <- MaskImages(object = STeye)
ImagePlot(STeye, method = "raster", type = "masked", ncols = 2)
FeatureOverlay(STeye, features = "labels", ncols = 2)
ImagePlot(SP, method = "raster", type = "masked", indices = 4)
ImagePlot(SP, method = "raster", type = "masked", indices = 5)
ImagePlot(SP, method = "raster", type = "masked", indices = 6)
ST.FeaturePlot(object = SP, features = "seurat_clusters_harmony", pt.size = 3, indices = 4) + 
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 8), legend.position = "bottom")
ST.FeaturePlot(object = SP, features = "seurat_clusters_harmony", pt.size = 3, indices = 5) + 
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 8), legend.position = "bottom")

png("he_1.png", width = 2000, height = 1000)
ImagePlot(SP, method = "raster", type = "masked", indices = c(3,4,5,6), annotate = TRUE)
dev.off()

##plot the UMAP##
##figure 1C
DimPlot(SP, reduction = "umap.harmony", label = TRUE, pt.size = 2.0, label.size = 12) + 
  theme(legend.text = element_text(size = 16),legend.title = element_text(size = 16),legend.position = "right",axis.title.x = element_text(size = 16, face = "bold", hjust = 0.5), axis.title.y = element_text(size = 16, face = "bold", hjust = 0.5)) +
  xlab("UMAP 1") + 
  ylab("UMAP 2")
DimPlot(SP, reduction = "umap.harmony", label = TRUE, pt.size = 1.3, label.size = 4.5)
ST.FeaturePlot(object = SP, features = "seurat_clusters_harmony", pt.size = 4.5, indices = c(3, 4, 5, 6), ncol = 2) +
  theme(legend.text = element_text(size = 24), legend.title = element_text(size = 24), legend.position = "bottom")

DimPlot(SP, reduction = "umap.harmony", label = TRUE, pt.size = 2.0, label.size = 12) + 
  theme(legend.text = element_text(size = 16),legend.title = element_text(size = 16),legend.position = "right",axis.title.x = element_text(size = 16, face = "bold", hjust = 0.5), axis.title.y = element_text(size = 16, face = "bold", hjust = 0.5)) +
  xlab("UMAP 1") + 
  ylab("UMAP 2")


##inferCNV##
##figure 1D
counts_matrix = GetAssayData(SP, slot="counts")
#switch active identity to cell type
Idents(STeye) <- "labels2"
DefaultAssay(STeye) <- "RNA"
STeye@meta.data$tissue_type <- STeye@active.ident
DimPlot(STeye, reduction = "umap", label = TRUE, pt.size = 1.3, label.size = 4.5)
#subset the retina part
STretina <- subset(STeye, idents = "Retina", invert = FALSE)
## merge
RB_STretina.combined <- merge(SP, y = STretina, add.cell.ids = c("SP", "STeye"), project = "Merge")
RB_STretina.combined@meta.data[5998:6002,]
RB_STretina.combined@meta.data$tissuetype <-c(RB_STretina.combined@meta.data[1:5999,5],RB_STretina.combined@meta.data[6000:6443,22])
RB_STretina.combined@meta.data[5995:6005,]
saveRDS(RB_STretina.combined,"RB_STretina.combined.rds")
#prepare the input data for infercnv
dfcount = GetAssayData(RB_STretina.combined, slot="counts")
dfids <- rownames(RB_STretina.combined@meta.data)
tissue_types <- RB_STretina.combined@meta.data$tissuetype
tissue_info <- data.frame(CellID = dfids, CellType = tissue_types, stringsAsFactors = FALSE)
write.table(tissue_info, file = "RB_STretina_tissue_info.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#estimate CNV
#create the infercnv object for cell type compare
infercnv_obj_RB_retina = CreateInfercnvObject(raw_counts_matrix=dfcount,
                                              annotations_file="RB_STretina_tissue_info.txt",
                                              delim="\t",
                                              gene_order_file="hg38_gencode_v27.txt",
                                              ref_group_names=c("Retina"))
# perform infercnv operations to reveal cnv signal
infercnv_obj_RB_retina = infercnv::run(infercnv_obj_RB_retina,
                                       cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                       out_dir="RB_STretina_output_dir",  # dir is auto-created for storing outputs
                                       cluster_by_groups=T,   # cluster
                                       denoise=T,
                                       HMM=T, 
                                       num_threads = 16,
                                       write_expr_matrix = T)
#figure 1D was contained in the output folder#

##CNV scores calculation##
##supplementary figure 5
#The CNV score was defined as the mean squares of CNV values across the genome
cnv_score_table_v3 = data.table::fread("RB_STretina_output_dir/infercnv.observations.txt", 
                                       data.table = F) %>% 
  column_to_rownames(var = 'V1')
library(scales)
cnvScore <- function(data){
  data <- data %>% as.matrix() %>%
    t() %>% 
    scale() %>% 
    rescale(to=c(-1, 1)) %>% 
    t()
  
  cnv_score <- as.data.frame(colSums(data * data))
  return(cnv_score)
}
cnv_score_v3 <- cnvScore(cnv_score_table_v3)
cellAnnota <- subset(RB_STretina.combined@meta.data, select = c('tissuetype'))
cnv_score_v3 <- cbind(cnv_score_v3, cellAnnota[row.names(cnv_score_v3),])
names(cnv_score_v3) <- c('cnv_score_v3', 'celltype')
color <- ggsci::pal_aaas()(10)

#plotting
ggplot(cnv_score_v3, aes(celltype, cnv_score_v3, color = celltype)) +
  geom_boxplot() +
  scale_color_manual(values = color) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "NA") +
  labs(x = '', y = 'CNV Scores', title = '') +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1)) +
  stat_compare_means()

cnv_score_v3 <- cnvScore(cnv_score_table_v3)
RB_combined_subset <- subset(RB_STretina.combined, subset = labels != "Retina")
cellAnnota <- subset(RB_combined_subset@meta.data, select = c('seurat_clusters_harmony'))
cnv_score_v3 <- cbind(cnv_score_v3, cellAnnota[row.names(cnv_score_v3),])
names(cnv_score_v3) <- c('cnv_score_v3', 'Clusters')
color <- ggsci::pal_aaas()(10)
ggplot(cnv_score_v3, aes(Clusters, cnv_score_v3, color = Clusters)) +
  geom_boxplot() +
  scale_color_manual(values = color) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "NA") +
  labs(x = '', y = 'CNV Scores', title = '') +
  theme(axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1)) +
  stat_compare_means()
figureS5 <- ggplot(cnv_score_v3, aes(Clusters, cnv_score_v3, fill = Clusters)) +
  geom_boxplot() +
  scale_fill_manual(values = color) +
  theme(panel.background = element_blank()) +
  theme(axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_blank()) +
  theme(legend.position = "right", legend.title = element_text(size = 16), legend.text = element_text(size = 16), legend.key.size = unit(1.2, "cm")) +
  labs(x = '', y = 'CNV Scores', title = '') +
  theme(axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20)) +
  stat_compare_means(size = 6)
figureS5

##3D stack showing the spatial distribution of gene expression##
##figure 1E-H
#By runnig Create3DStack, we can create a z-stack of “2D point patterns” which we’ll use to interpolate expression values over and visualzie expression in 2D space. 
SP <- Create3DStack(SP)
#plot RB genes 
FeaturePlot3D(SP, features = "UBE2C", pt.size = 0.6, max.cutoff = 5)
FeaturePlot3D(SP, features = "RB1", pt.size = 0.6, max.cutoff = 5)
FeaturePlot3D(SP, features = "MDM2", pt.size = 0.6, max.cutoff = 5)
FeaturePlot3D(SP, features = "MYCN", pt.size = 0.6, max.cutoff = 5)

##check the gene expression across cell clusters##
##supplementary figure 2##
RB_check <- c("RB1", "MDM2", "E2F1", "MYCN", "UBE2C", "MDM4")
figureS2 <- DotPlot(object = SP2, features = RB_check, cols = c("cornsilk", "red2"), dot.scale = 14, 
                    group.by = "seurat_clusters_harmony", cluster.idents = F, scale.by = "size") + 
  theme(legend.position = "right", 
        axis.title = element_text(size = 20), axis.text = element_text(size = 20), 
        legend.text = element_text(size = 20), legend.title = element_text(size = 20)) + 
  labs(y = "Cluster")
figureS2

##Module score showing the cell type markers
##figure 2A-E
#Retinal progenitor cell type
RPCscore <- Seurat::AddModuleScore(object = SP, features = c("SOX2", "HES1", "MKI67", "HES5", "FZD5", "PAX6"), name = "RPC_enriched")
colnames(RPCscore@meta.data)[25] <- "RPC_enriched"
pRPC <- VlnPlot(RPCscore,features = "RPC_enriched", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
pRPC
#Pigment epithelium cell type
PECscore <- Seurat::AddModuleScore(object = SP, features = c("SERPINF1", "MITF", "BEST1", "TTR"), name = "PEC_enriched")
colnames(PECscore@meta.data)[23] <- "PEC_enriched"
pPEC <- VlnPlot(PECscore, features = "PEC_enriched", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
pPEC
#Rod cell type
RODscore <- Seurat::AddModuleScore(object = SP, features = c("RHO", "PDE6A", "CNGA1", "NRL", "GNAT1", "GNB1", "SAG", "ELOVL4", "PDE6B", "GNGT1"), name = "ROD_enriched")
colnames(RODscore@meta.data)[29] <- "ROD_enriched"
pROD <- VlnPlot(RODscore,features = "ROD_enriched", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
pROD
#Cone cell type
CONEscore <- Seurat::AddModuleScore(object = SP, features = c("ARR3", "GNGT2", "PDE6H", "GUCA1C", "GNAT2", "RXRG", "THRB", "PDC", "GNB3", "CRX"), name = "CONE_enriched")
colnames(CONEscore@meta.data)[29] <- "CONE_enriched"
pCONE <- VlnPlot(CONEscore,features = "CONE_enriched", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
pCONE
#Bipolar cell type
BCscore <- Seurat::AddModuleScore(object = SP, features = c("GRM6", "VSX2"), name = "BC_enriched")
colnames(BCscore@meta.data)[21] <- "BC_enriched"
pBC <- VlnPlot(BCscore,features = "BC_enriched", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
#Horizontal cell type
HCscore <- Seurat::AddModuleScore(object = SP, features = c("ONECUT1", "ONECUT2", "TFAP2B"), name = "HC_enriched")
colnames(HCscore@meta.data)[22] <- "HC_enriched"
pHC <- VlnPlot(HCscore,features = "HC_enriched", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
#Amacrine cell type
ACscore <- Seurat::AddModuleScore(object = SP, features = c("GAD1", "CALB1", "NRXN2", "TFAP2A", "PROX1", "GAD2"), name = "AC_enriched")
colnames(ACscore@meta.data)[25] <- "AC_enriched"
pAC <- VlnPlot(ACscore,features = "AC_enriched", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
#Retinal ganglion cell type
RGCscore <- Seurat::AddModuleScore(object = SP, features = c("POU4F2", "GAP43", "NEFL", "SNCG", "ATOH7", "EBF3", "THY1", "NRN1"), name = "RGC_enriched")
colnames(RGCscore@meta.data)[27] <- "RGC_enriched"
pRGC <- VlnPlot(RGCscore,features = "RGC_enriched", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
#Cone precursor cell type
CPscore <- Seurat::AddModuleScore(object = SP, features = c("CRX", "RXRG", "THRB"), name = "CP_enriched")
colnames(CPscore@meta.data)[22] <- "CP_enriched"
pCP <- VlnPlot(CPscore,features = "CP_enriched", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
#glial cells
Gscore <- Seurat::AddModuleScore(object = SP, features = c("CD68", "HLA-DPA1", "HLA-DPB1", "CLU"), name = "Glial_enriched")
colnames(Gscore@meta.data)[23] <- "Glial_enriched"
pG <- VlnPlot(Gscore,features = "Glial_enriched", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
#cancer-associated fibroblasts (CAF)
CAFscore <- Seurat::AddModuleScore(object = SP, features = c("ACTA2", "VIM", "FGF9"), name = "CAF_enriched")
colnames(CAFscore@meta.data)[22] <- "CAF_enriched"
pCAF <- VlnPlot(CAFscore,features = "CAF_enriched", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))
#proliferative cells
proscore <- Seurat::AddModuleScore(object = SP, features = c("MKI67", "TOP2A", "KIF14"), name = "HPCP_enriched")
colnames(proscore@meta.data)[22] <- "HPCP_enriched"
ppro <- VlnPlot(proscore,features = "HPCP_enriched", slot = "counts", log = TRUE) + scale_y_continuous(limits = c(0.5,3.5))

pRPC <- pRPC + labs(title = "Retinal Progenitor Cell Score") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
pPEC <- pPEC + labs(title = "RPE Cell Score") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
pROD <- pROD + labs(title = "Rod Cell Score") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
pCONE <- pCONE + labs(title = "Cone Cell Score") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
pHC <- pHC + labs(title = "Horizontal Cell Score") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
pAC <- pAC + labs(title = "Amacrine Cell Score") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
pRGC <- pRGC + labs(title = "RGC Score") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
pBC <- pBC + labs(title = "Bipolar Cell Score") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
plot_grid(pRPC, pPEC, pROD, pCONE, pHC, pAC, pRGC, pBC, nrow = 2)


##Spatial deconvolution##
##Conditional autoregressive-based deconvolution (CARD)##
stRB <- SP
#single cell: Wu, C. et al. Single-cell characterization of malignant phenotypes and microenvironment alteration in retinoblastoma. Cell Death Dis. 13, 438 (2022).
scRB <- readRDS("RB_data.RDS")
scRB = UpdateSeuratObject(object = scRB)
#generate sp_count and sp_location
sp_count <- stRB@assays$RNA@counts
RB4spatial <- SubsetSTData(stRB, expression = labels %in% "RB 1D")
ImagePlot(stRB, method = "raster")
ImagePlot(RB4spatial, method = "raster")
st.RB4 <- GetStaffli(RB4spatial)
st.RB4
RB4_coordinate <- st.RB4@meta.data
RB4sp_location <- data.frame(RB4_coordinate)
RB4sp_location$x <- NULL
RB4sp_location$adj_y <- NULL
RB4sp_location$pixel_x <- NULL
RB4sp_location$pixel_y <- NULL
RB4sp_location$original_x <- NULL
RB4sp_location$original_y <- NULL
RB4sp_location$warped_x <- NULL
RB4sp_location$warped_y <- NULL
RB4sp_location$sample <- NULL
colnames(RB4sp_location) <- c("x", "y")
plot(RB4sp_location, method = "raster")
RB4sp_count <- sp_count[, colnames(sp_count) %in% rownames(RB4sp_location)]
RB3spatial <- SubsetSTData(stRB, expression = labels %in% "RB 1C")
ImagePlot(stRB, method = "raster")
ImagePlot(RB3spatial, method = "raster")
st.RB3 <- GetStaffli(RB3spatial)
st.RB3
RB3_coordinate <- st.RB3@meta.data
RB3sp_location <- data.frame(RB3_coordinate)
RB3sp_location$x <- NULL
RB3sp_location$adj_y <- NULL
RB3sp_location$pixel_x <- NULL
RB3sp_location$pixel_y <- NULL
RB3sp_location$original_x <- NULL
RB3sp_location$original_y <- NULL
RB3sp_location$warped_x <- NULL
RB3sp_location$warped_y <- NULL
RB3sp_location$sample <- NULL
colnames(RB3sp_location) <- c("x", "y")
plot(RB3sp_location, method = "raster")
RB3sp_count <- sp_count[, colnames(sp_count) %in% rownames(RB3sp_location)]
RB2spatial <- SubsetSTData(stRB, expression = labels %in% "RB 1B")
ImagePlot(stRB, method = "raster")
ImagePlot(RB2spatial, method = "raster")
st.RB2 <- GetStaffli(RB2spatial)
st.RB2
RB2_coordinate <- st.RB2@meta.data
RB2sp_location <- data.frame(RB2_coordinate)
RB2sp_location$x <- NULL
RB2sp_location$adj_y <- NULL
RB2sp_location$pixel_x <- NULL
RB2sp_location$pixel_y <- NULL
RB2sp_location$original_x <- NULL
RB2sp_location$original_y <- NULL
RB2sp_location$warped_x <- NULL
RB2sp_location$warped_y <- NULL
RB2sp_location$sample <- NULL
colnames(RB2sp_location) <- c("x", "y")
plot(RB2sp_location, method = "raster")
RB2sp_count <- sp_count[, colnames(sp_count) %in% rownames(RB2sp_location)]
RB1spatial <- SubsetSTData(stRB, expression = labels %in% "RB 1A")
ImagePlot(stRB, method = "raster")
ImagePlot(RB1spatial, method = "raster")
st.RB1 <- GetStaffli(RB1spatial)
st.RB1
RB1_coordinate <- st.RB1@meta.data
RB1sp_location <- data.frame(RB1_coordinate)
RB1sp_location$x <- NULL
RB1sp_location$adj_y <- NULL
RB1sp_location$pixel_x <- NULL
RB1sp_location$pixel_y <- NULL
RB1sp_location$original_x <- NULL
RB1sp_location$original_y <- NULL
RB1sp_location$warped_x <- NULL
RB1sp_location$warped_y <- NULL
RB1sp_location$sample <- NULL
colnames(RB1sp_location) <- c("x", "y")
plot(RB1sp_location, method = "raster")
RB1sp_count <- sp_count[, colnames(sp_count) %in% rownames(RB1sp_location)]

RB2Aspatial <- SubsetSTData(stRB, expression = labels %in% "RB 2A")
ImagePlot(stRB, method = "raster")
ImagePlot(RB2Aspatial, method = "raster")
st.RB2A <- GetStaffli(RB2Aspatial)
st.RB2A
RB2A_coordinate <- st.RB2A@meta.data
RB2Asp_location <- data.frame(RB2A_coordinate)
RB2Asp_location$x <- NULL
RB2Asp_location$adj_y <- NULL
RB2Asp_location$pixel_x <- NULL
RB2Asp_location$pixel_y <- NULL
RB2Asp_location$original_x <- NULL
RB2Asp_location$original_y <- NULL
RB2Asp_location$warped_x <- NULL
RB2Asp_location$warped_y <- NULL
RB2Asp_location$sample <- NULL
colnames(RB2Asp_location) <- c("x", "y")
plot(RB2Asp_location, method = "raster")
RB2Asp_count <- sp_count[, colnames(sp_count) %in% rownames(RB2Asp_location)]
RB2Bspatial <- SubsetSTData(stRB, expression = labels %in% "RB 2B")
ImagePlot(stRB, method = "raster")
ImagePlot(RB2Bspatial, method = "raster")
st.RB2B <- GetStaffli(RB2Bspatial)
st.RB2B
RB2B_coordinate <- st.RB2B@meta.data
RB2Bsp_location <- data.frame(RB2B_coordinate)
RB2Bsp_location$x <- NULL
RB2Bsp_location$adj_y <- NULL
RB2Bsp_location$pixel_x <- NULL
RB2Bsp_location$pixel_y <- NULL
RB2Bsp_location$original_x <- NULL
RB2Bsp_location$original_y <- NULL
RB2Bsp_location$warped_x <- NULL
RB2Bsp_location$warped_y <- NULL
RB2Bsp_location$sample <- NULL
colnames(RB2Bsp_location) <- c("x", "y")
plot(RB2Bsp_location, method = "raster")
RB2Bsp_count <- sp_count[, colnames(sp_count) %in% rownames(RB2Bsp_location)]
#generate sc_count and sc_meta
scRB
scRB@meta.data$celltype
scRB@assays$RNA
scRB@active.ident
scRB@meta.data$orig.ident
count.data <- GetAssayData(object = scRB[["RNA"]], slot = "counts")
count.data
sc_count <- count.data
cellType <- scRB@meta.data$celltype
cellType
cellID <- scRB@active.ident
cellID
sampleInfo <- scRB@meta.data$orig.ident
sampleInfo
meta.data <- data.frame(cellID, cellType, sampleInfo)
meta.data$cellID <- NULL
meta.data
meta.data <- tibble::rownames_to_column(meta.data, "cellID")
rownames(meta.data) <- meta.data$cellID
sc_meta <- meta.data
print(dim(sc_count))
print(dim(sc_meta))
table(sc_meta$cellType)
##get rid of Neural Cell and Others in the sc dataset to re-do the deconvonlution
scRB <- readRDS("RB_data.RDS")
scRB_updated = UpdateSeuratObject(object = scRB)
scRB = UpdateSeuratObject(object = scRB)
scRB$celltype
celltype <- c("CP","HP-CP","Other","Neural cell","CAF","Rod-like","Glial","Cone-like")
Idents(scRB_updated) <- scRB_updated$celltype
names(celltype) <- levels(scRB_updated)
scRB_updated <- RenameIdents(scRB_updated, celltype)
DimPlot(object = scRB_updated, reduction = "umap", label = T)
scRB_updated <- subset(x = scRB_updated, idents = c("Neural cell", "Other"), invert = TRUE)
scRB_updated$celltype
scRB
scRB_updated
DimPlot(object = scRB_updated, reduction = "umap", label = T)
remove(scRB)
scRB_updated@meta.data$celltype
scRB_updated@assays$RNA
scRB_updated@active.ident
scRB_updated@meta.data$orig.ident
count.data <- GetAssayData(object = scRB_updated[["RNA"]], slot = "counts")
count.data
sc_count <- count.data
cellType <- scRB_updated@meta.data$celltype
cellType
cellID <- scRB_updated@active.ident
cellID
sampleInfo <- scRB_updated@meta.data$orig.ident
sampleInfo
meta.data <- data.frame(cellID, cellType, sampleInfo)
meta.data$cellID <- NULL
meta.data
meta.data <- tibble::rownames_to_column(meta.data, "cellID")
rownames(meta.data) <- meta.data$cellID
sc_meta <- meta.data
print(dim(sc_count))
print(dim(sc_meta))
table(sc_meta$cellType)
###RB1D
###create CARD objects
CARD_RB4obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = RB4sp_count,
  spatial_location = RB4sp_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5)
#Run deconvolution
CARD_RB4object = CARD_deconvolution(CARD_RB4obj)
##save results
saveRDS(CARD_RB4obj, "CARD_RB1Dobj_modified.RDS")
saveRDS(CARD_RB4object, "CARD_RB1Dobject_modified.RDS")
##visualization
print(CARD_RB4object@Proportion_CARD[1:2,])
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p1D <- CARD.visualize.pie(proportion = CARD_RB4object@Proportion_CARD,spatial_location = CARD_RB4object@spatial_location, colors = colors)
print(p1D)
###RB1C
#create CARD objects
CARD_RB3obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = RB3sp_count,
  spatial_location = RB3sp_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 
#Run deconvolution
CARD_RB3object = CARD_deconvolution(CARD_RB3obj)
##save results
saveRDS(CARD_RB3obj, "CARD_RB1Cobj_modified.RDS")
saveRDS(CARD_RB3object, "CARD_RB1Cobject_modified.RDS")
##visualization
print(CARD_RB3object@Proportion_CARD[1:2,])
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p1C <- CARD.visualize.pie(proportion = CARD_RB3object@Proportion_CARD,spatial_location = CARD_RB3object@spatial_location, colors = colors)
print(p1C)
##RB1B
#create CARD objects
CARD_RB2obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = RB2sp_count,
  spatial_location = RB2sp_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 
#Run deconvolution
CARD_RB2object = CARD_deconvolution(CARD_RB2obj)
##save results
saveRDS(CARD_RB2obj, "CARD_RB1Bobj_modified.RDS")
saveRDS(CARD_RB2object, "CARD_RB1Bobject_modified.RDS")
##visualization
print(CARD_RB2object@Proportion_CARD[1:2,])
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p1B <- CARD.visualize.pie(proportion = CARD_RB2object@Proportion_CARD,spatial_location = CARD_RB2object@spatial_location, colors = colors)
print(p1B)
##RB1A
#create CARD objects
CARD_RB1obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = RB1sp_count,
  spatial_location = RB1sp_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 
#Run deconvolution
CARD_RB1object = CARD_deconvolution(CARD_RB1obj)
##save results
saveRDS(CARD_RB1obj, "CARD_RB1Aobj_modified.RDS")
saveRDS(CARD_RB1object, "CARD_RB1Aobject_modified.RDS")
##visualization
print(CARD_RB1object@Proportion_CARD[1:2,])
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p1A <- CARD.visualize.pie(proportion = CARD_RB1object@Proportion_CARD,spatial_location = CARD_RB1object@spatial_location, colors = colors)
print(p1A)
##RB2A
#create CARD objects
CARD_RB2Aobj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = RB2Asp_count,
  spatial_location = RB2Asp_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 
#Run deconvolution
CARD_RB2Aobject = CARD_deconvolution(CARD_RB2Aobj)
##save results
saveRDS(CARD_RB2Aobj, "CARD_RB2Aobj_modified.RDS")
saveRDS(CARD_RB2Aobject, "CARD_RB2Aobject_modified.RDS")
##visualization
print(CARD_RB2Aobject@Proportion_CARD[1:2,])
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p2A <- CARD.visualize.pie(proportion = CARD_RB2Aobject@Proportion_CARD,spatial_location = CARD_RB2Aobject@spatial_location, colors = colors)
print(p2A)
##RB2B
#create CARD objects
CARD_RB2Bobj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = RB2Bsp_count,
  spatial_location = RB2Bsp_location,
  ct.varname = "cellType",
  ct.select = unique(sc_meta$cellType),
  sample.varname = "sampleInfo",
  minCountGene = 100,
  minCountSpot = 5) 
#Run deconvolution
CARD_RB2Bobject = CARD_deconvolution(CARD_RB2Bobj)
##save results
saveRDS(CARD_RB2Bobj, "CARD_RB2Bobj_modified.RDS")
saveRDS(CARD_RB2Bobject, "CARD_RB2Bobject_modified.RDS")
##visualization
print(CARD_RB2Bobject@Proportion_CARD[1:2,])
colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p2B <- CARD.visualize.pie(proportion = CARD_RB2Bobject@Proportion_CARD,spatial_location = CARD_RB2Bobject@spatial_location, colors = colors)
print(p2B)
#combine the pie charts
grid.arrange(p1A, p1B, p1C, p1D, p2A, p2B, nrow = 2)
#generate the results
RB1A_CARD <- data.frame(CARD_RB1object@Proportion_CARD)
RB1B_CARD <- data.frame(CARD_RB2object@Proportion_CARD)
RB1C_CARD <- data.frame(CARD_RB3object@Proportion_CARD)
RB1D_CARD <- data.frame(CARD_RB4object@Proportion_CARD)
RB2A_CARD <- data.frame(CARD_RB2Aobject@Proportion_CARD)
RB2B_CARD <- data.frame(CARD_RB2Bobject@Proportion_CARD)
write.csv(RB1A_CARD, "RB1A_CARD.CSV")
write.csv(RB1B_CARD, "RB1B_CARD.CSV")
write.csv(RB1C_CARD, "RB1C_CARD.CSV")
write.csv(RB1D_CARD, "RB1D_CARD.CSV")
write.csv(RB2A_CARD, "RB2A_CARD.CSV")
write.csv(RB2B_CARD, "RB2B_CARD.CSV")
spatial_cluster<-stRB$seurat_clusters_harmony
write.csv(spatial_cluster, "spatial_cluster.CSV")

##figure 2B-E##
RB1A <- readRDS("CARD_RB1Aobject_modified.RDS")
RB1B <- readRDS("CARD_RB1Bobject_modified.RDS")
RB1C <- readRDS("CARD_RB1Cobject_modified.RDS")
RB1D <- readRDS("CARD_RB1Dobject_modified.RDS")
RB2A <- readRDS("CARD_RB2Aobject_modified.RDS")
RB2B <- readRDS("CARD_RB2Bobject_modified.RDS")

colors = c("#FFD92F","#4DAF4A","#FCCDE5","#D9D9D9","#377EB8","#7FC97F","#BEAED4",
           "#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#1B9E77","#D95F02",
           "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D")
p1 <- CARD.visualize.pie(proportion = RB1A@Proportion_CARD,spatial_location = RB1A@spatial_location, colors = colors)
p2 <- CARD.visualize.pie(proportion = RB1B@Proportion_CARD,spatial_location = RB1B@spatial_location, colors = colors)
p3 <- CARD.visualize.pie(proportion = RB1C@Proportion_CARD,spatial_location = RB1C@spatial_location, colors = colors)
p4 <- CARD.visualize.pie(proportion = RB1D@Proportion_CARD,spatial_location = RB1D@spatial_location, colors = colors)
p5 <- CARD.visualize.pie(proportion = RB2A@Proportion_CARD,spatial_location = RB2A@spatial_location, colors = colors)
p6 <- CARD.visualize.pie(proportion = RB2B@Proportion_CARD,spatial_location = RB2B@spatial_location, colors = colors)
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)

p1 + theme(legend.text = element_text(size = 8), legend.title = element_text(size = 8)) 
p2 + theme(legend.text = element_text(size = 8), legend.title = element_text(size = 8))
p3 + theme(legend.text = element_text(size = 8), legend.title = element_text(size = 8))
p4 + theme(legend.text = element_text(size = 8), legend.title = element_text(size = 8))
p5 + theme(legend.text = element_text(size = 8), legend.title = element_text(size = 8))
p6 + theme(legend.text = element_text(size = 8), legend.title = element_text(size = 8))

print(RB1A@Proportion_CARD[1:2,])
print(RB1B@Proportion_CARD[1:2,])
print(RB1C@Proportion_CARD[1:2,])
print(RB1D@Proportion_CARD[1:2,])
print(RB2A@Proportion_CARD[1:2,])
print(RB2B@Proportion_CARD[1:2,])

ST.FeaturePlot(object = stRB, features = "seurat_clusters_harmony", pt.size = 3, indices = "RB 1A") + 
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 8), legend.position = "bottom")

ST.FeaturePlot(object = stRB, features = "seurat_clusters_harmony", pt.size = 3, indices = "RB 1B") + 
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 8), legend.position = "bottom")

ST.FeaturePlot(object = stRB, features = "seurat_clusters_harmony", pt.size = 3, indices = "RB 1C") + 
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 8), legend.position = "bottom")

ST.FeaturePlot(object = stRB, features = "seurat_clusters_harmony", pt.size = 3, indices = "RB 1D") + 
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 8), legend.position = "bottom")

ST.FeaturePlot(object = stRB, features = "seurat_clusters_harmony", pt.size = 3, indices = "RB 2A") + 
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 8), legend.position = "bottom")

ST.FeaturePlot(object = stRB, features = "seurat_clusters_harmony", pt.size = 3, indices = "RB 2B") + 
  theme(legend.text = element_text(size = 10), legend.title = element_text(size = 8), legend.position = "bottom")

## select the cell type that we are interested
ct.visualize = c("CP", "MKI67+ CP", "Glial","CAF")
## visualize the spatial distribution of the cell type proportion
p3_2 <- CARD.visualize.prop(
  proportion = RB1C@Proportion_CARD,        
  spatial_location = RB1C@spatial_location, 
  ct.visualize = ct.visualize,                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 2)                                 ### number of columns in the figure panel
print(p3_2)
p4_2 <- CARD.visualize.prop(
  proportion = RB1D@Proportion_CARD,        
  spatial_location = RB1D@spatial_location, 
  ct.visualize = ct.visualize,                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 2)                                 ### number of columns in the figure panel
print(p4_2)
p5_2 <- CARD.visualize.prop(
  proportion = RB2A@Proportion_CARD,        
  spatial_location = RB2A@spatial_location, 
  ct.visualize = ct.visualize,                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 2)                                 ### number of columns in the figure panel
print(p5_2)
p6_2 <- CARD.visualize.prop(
  proportion = RB2B@Proportion_CARD,        
  spatial_location = RB2B@spatial_location, 
  ct.visualize = ct.visualize,                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 2)                                 ### number of columns in the figure panel
print(p6_2)

pCP <- pCP + labs(title = "Cone Precursor Cell Score") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
ppro <- ppro + labs(title = "Highly-Proliferated Cone Precursor Cell Score") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
pG <- pG + labs(title = "Glial Cell Score") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
pCAF <- pCAF + labs(title = "Cancer-Associated Fibroblast Score") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
plot_grid(pCP, ppro, pG, pCAF, nrow = 4)



##cell cycle analysis##
##figure 3A##
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
SP <- CellCycleScoring(SP, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
RidgePlot(SP, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 4)

##spatial distribution of cell cycle scores
##figure 3B##
#check the python script for the spatial plot of cell cycle score#

##statistics showing the cell cycle phase
##figure 3C##
pbar <- SP@meta.data %>%
  group_by(seurat_clusters_harmony,Phase) %>%
  dplyr::count() %>%
  group_by(seurat_clusters_harmony) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters_harmony,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")
pbar

list_G2m_SP <- SP@meta.data %>%
  group_by(seurat_clusters, Phase) %>%
  dplyr::count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

##spatial distribution of enriched cell types resulted from spatial deconvolution##
##need to run the codes generating Table1 first to do the cell type enrichment first##
##supplementary figure 3##
#map back the spot identity back to the seurat object for spatial resolution
lookup_df <- enriched_spots_all[, c("barcode", "celltype")]
lookup_df <- lookup_df[!duplicated(lookup_df$barcode), ]
meta_new <- SP@meta.data %>%
  mutate(barcode = rownames(.)) %>%
  left_join(lookup_df, by = "barcode")
SP@meta.data$CARD_enriched_cell_type <- meta_new$celltype
table(SP@meta.data$CARD_enriched_cell_type, useNA = "always")

SP_CARD <- subset(SP, subset = CARD_enriched_cell_type %in% CARD_enriched_cell_type)
SP_CARD@meta.data <- SP_CARD@meta.data %>%mutate(CARD_enriched_cell_type = ifelse(CARD_enriched_cell_type == "MKI67..CP", "HP-CP", CARD_enriched_cell_type))
figureS3 <- ST.FeaturePlot(object = SP_CARD, features = "CARD_enriched_cell_type", pt.size = 3, indices = c(4, 5), ncol = 2) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20), legend.position = "bottom") + 
  guides(fill = guide_legend(override.aes = list(size = 6)))
figureS3

##RNA velocity##
##figure 3 D-I##
#check the python script of scVelo and CellRank#

##DEG analysis##
##figure 4A##
de.markers <- FindAllMarkers(SP, only.pos = TRUE)
top10 <- de.markers %>%
  dplyr::filter(p_val_adj < 0.01) %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(wt = -p_val_adj, n = 10)
SP.DEheatmap<- DoHeatmap(SP, features = top10$gene)
SP.DEheatmap

##gene set enrichment analysis with clusterprofiler##
##figure 4B-E##
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
C0M <- FindMarkers(SP, ident.1 = "C0")
C1M <- FindMarkers(SP, ident.1 = "C1")
C2M <- FindMarkers(SP, ident.1 = "C2")
C3M <- FindMarkers(SP, ident.1 = "C3")
C4M <- FindMarkers(SP, ident.1 = "C4")
C5M <- FindMarkers(SP, ident.1 = "C5")
C6M <- FindMarkers(SP, ident.1 = "C6")
C7M <- FindMarkers(SP, ident.1 = "C7")
C8M <- FindMarkers(SP, ident.1 = "C8")

#umap.harmony clusters
original_C0_list <- C0M$avg_log2FC
names(original_C0_list) <- rownames(C0M)
C0_list <- na.omit(original_C0_list)
C0_list = sort(C0_list, decreasing = TRUE)
C0_list
original_C1_list <- C1M$avg_log2FC
names(original_C1_list) <- rownames(C1M)
C1_list <- na.omit(original_C1_list)
C1_list = sort(C1_list, decreasing = TRUE)
C1_list
original_C2_list <- C2M$avg_log2FC
names(original_C2_list) <- rownames(C2M)
C2_list <- na.omit(original_C2_list)
C2_list = sort(C2_list, decreasing = TRUE)
C2_list
original_C3_list <- C3M$avg_log2FC
names(original_C3_list) <- rownames(C3M)
C3_list <- na.omit(original_C3_list)
C3_list = sort(C3_list, decreasing = TRUE)
C3_list
original_C4_list <- C4M$avg_log2FC
names(original_C4_list) <- rownames(C4M)
C4_list <- na.omit(original_C4_list)
C4_list = sort(C4_list, decreasing = TRUE)
C4_list
original_C5_list <- C5M$avg_log2FC
names(original_C5_list) <- rownames(C5M)
C5_list <- na.omit(original_C5_list)
C5_list = sort(C5_list, decreasing = TRUE)
C5_list
original_C6_list <- C6M$avg_log2FC
names(original_C6_list) <- rownames(C6M)
C6_list <- na.omit(original_C6_list)
C6_list = sort(C6_list, decreasing = TRUE)
C6_list
original_C7_list <- C7M$avg_log2FC
names(original_C7_list) <- rownames(C7M)
C7_list <- na.omit(original_C7_list)
C7_list = sort(C7_list, decreasing = TRUE)
C7_list
original_C8_list <- C8M$avg_log2FC
names(original_C8_list) <- rownames(C8M)
C8_list <- na.omit(original_C8_list)
C8_list = sort(C8_list, decreasing = TRUE)
C8_list
#BP for biological process, MF for molecular function, CC for cell compartment, ALL for all
keytypes(org.Hs.eg.db)
C0gse <- gseGO(geneList=C0_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C1gse <- gseGO(geneList=C1_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C2gse <- gseGO(geneList=C2_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C3gse <- gseGO(geneList=C3_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C4gse <- gseGO(geneList=C4_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C5gse <- gseGO(geneList=C5_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C6gse <- gseGO(geneList=C6_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C7gse <- gseGO(geneList=C7_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
C8gse <- gseGO(geneList=C8_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")

require(DOSE)
pCP0D <- dotplot(C0gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP0C <- cnetplot(C0gse, categorySize="pvalue", foldChange=C0_list, showCategory = 3)

pCP1D <- dotplot(C1gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP1C <- cnetplot(C1gse, categorySize="pvalue", foldChange=C1_list, showCategory = 3)

pCP2D <- dotplot(C2gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP2C <- cnetplot(C2gse, categorySize="pvalue", foldChange=C2_list, showCategory = 3)

pCP3D <- dotplot(C3gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP3C <- cnetplot(C3gse, categorySize="pvalue", foldChange=C3_list, showCategory = 3)

pCP4D <- dotplot(C4gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP4C <- cnetplot(C4gse, categorySize="pvalue", foldChange=C4_list, showCategory = 3)

pCP5D <- dotplot(C5gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP5C <- cnetplot(C5gse, categorySize="pvalue", foldChange=C5_list, showCategory = 3)

pCP6D <- dotplot(C6gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP6C <- cnetplot(C6gse, categorySize="pvalue", foldChange=C6_list, showCategory = 3)

pCP7D <- dotplot(C7gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP7C <- cnetplot(C7gse, categorySize="pvalue", foldChange=C7_list, showCategory = 3)

pCP8D <- dotplot(C8gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCP8C <- cnetplot(C8gse, categorySize="pvalue", foldChange=C8_list, showCategory = 3)

#figure 4B-E
pCP1D <- dotplot(C1gse, showCategory=5, split=".sign") + 
  facet_grid(.~.sign) +  
  theme(
    text = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.text.x = element_text(size = 20),
    axis.title = element_text(size = 30),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30)
  )
pCP1C <- cnetplot(C1gse, categorySize="pvalue", foldChange=C1_list, showCategory = 3) +  
  theme(
    text = element_text(size = 30),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 30)
  )

pCP4D <- dotplot(C4gse, showCategory=5, split=".sign") + facet_grid(.~.sign) + 
  facet_grid(.~.sign) +  
  theme(
    text = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.text.x = element_text(size = 20),
    axis.title = element_text(size = 30),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 30),
    legend.title = element_text(size = 30)
  )
pCP4C <- cnetplot(C4gse, categorySize="pvalue", foldChange=C4_list, showCategory = 3) +  
  theme(
    text = element_text(size = 30),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 30)
  )
pCP1D
pCP1C
pCP4D
pCP4C
figure4.2 <- plot_grid(pCP1D, pCP1C, pCP4D, pCP4C, ncol = 2)
figure4.2

##supplementary figure 6##
pCP0D <- dotplot(C0gse, showCategory=5, split=".sign") + facet_grid(.~.sign) +  
  theme(
    text = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    axis.text.x = element_text(size = 20),
    axis.title = element_text(size = 25),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )
pCP0C <- cnetplot(C0gse, categorySize="pvalue", foldChange=C0_list, showCategory = 3) +  
  theme(
    text = element_text(size = 30),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 30)
  )

pCP2D <- dotplot(C2gse, showCategory=5, split=".sign") + facet_grid(.~.sign) +  
  theme(
    text = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    axis.text.x = element_text(size = 20),
    axis.title = element_text(size = 25),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )
pCP2C <- cnetplot(C2gse, categorySize="pvalue", foldChange=C2_list, showCategory = 3) +  
  theme(
    text = element_text(size = 30),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 30)
  )

pCP3D <- dotplot(C3gse, showCategory=5, split=".sign") + facet_grid(.~.sign) +  
  theme(
    text = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    axis.text.x = element_text(size = 20),
    axis.title = element_text(size = 25),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )
pCP3C <- cnetplot(C3gse, categorySize="pvalue", foldChange=C3_list, showCategory = 3) +  
  theme(
    text = element_text(size = 30),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 30)
  )

pCP5D <- dotplot(C5gse, showCategory=5, split=".sign") + facet_grid(.~.sign) +  
  theme(
    text = element_text(size = 25),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 20),
    axis.title = element_text(size = 25),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )
pCP5C <- cnetplot(C5gse, categorySize="pvalue", foldChange=C5_list, showCategory = 3) +  
  theme(
    text = element_text(size = 30),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 30)
  )

pCP6D <- dotplot(C6gse, showCategory=5, split=".sign") + facet_grid(.~.sign) +  
  theme(
    text = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    axis.text.x = element_text(size = 20),
    axis.title = element_text(size = 25),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )
pCP6C <- cnetplot(C6gse, categorySize="pvalue", foldChange=C6_list, showCategory = 3) +  
  theme(
    text = element_text(size = 30),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 30)
  )

pCP7D <- dotplot(C7gse, showCategory=5, split=".sign") + facet_grid(.~.sign) +  
  theme(
    text = element_text(size = 25),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 20),
    axis.title = element_text(size = 25),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )
pCP7C <- cnetplot(C7gse, categorySize="pvalue", foldChange=C7_list, showCategory = 3) +  
  theme(
    text = element_text(size = 30),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 30)
  )

pCP8D <- dotplot(C8gse, showCategory=5, split=".sign") + facet_grid(.~.sign) +  
  theme(
    text = element_text(size = 25),
    axis.text.y = element_text(size = 25),
    axis.text.x = element_text(size = 20),
    axis.title = element_text(size = 25),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20)
  )
pCP8C <- cnetplot(C8gse, categorySize="pvalue", foldChange=C8_list, showCategory = 3) +  
  theme(
    text = element_text(size = 30),
    strip.text = element_text(size = 30),
    plot.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 30)
  )


pCP0_combined <- (pCP0D | pCP0C) + 
  plot_annotation(title = "C0", theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")))
pCP2_combined <- (pCP2D | pCP2C) + 
  plot_annotation(title = "C2", theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")))
pCP3_combined <- (pCP3D | pCP3C) + 
  plot_annotation(title = "C3", theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")))
pCP5_combined <- (pCP5D | pCP5C) + 
  plot_annotation(title = "C5", theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")))
pCP6_combined <- (pCP6D | pCP6C) + 
  plot_annotation(title = "C6", theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")))
pCP7_combined <- (pCP7D | pCP7C) + 
  plot_annotation(title = "C7", theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")))
pCP8_combined <- (pCP8D | pCP8C) + 
  plot_annotation(title = "C8", theme = theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")))

pCP0_combined
pCP2_combined
pCP3_combined
pCP5_combined
pCP6_combined
pCP7_combined
pCP8_combined

figureS6 <- (pCP0_combined | pCP2_combined) /
  (pCP3_combined | pCP5_combined) /
  (pCP6_combined | pCP7_combined) /
  (pCP8_combined | plot_spacer())

png("figureS6.png", width = 5000, height = 4000, res = 100)
print(figureS6)
dev.off()


##compare the RB samples with control human retina##
#supplementary figure 4##
#intergrate the RB and Eye dataset
hRetina <- subset(STeye, idents = "Retina")
Idents(STeye)
Idents(STeye) <- STeye@meta.data$labels
hChoroid <- subset(STeye, idents = "Choroid")
hSclera <- subset(STeye, idents = "Sclera")
hOptic_nerve <- subset(STeye, idents = "Optic nerve")
hRetina$type <- "Healthy_Retina"
hChoroid$type <- "Choroid"
hSclera$type <- "Sclera"
hOptic_nerve$type <- "Optic_nerve"
SP$type <- "RB"
#integration between the splitted datasets
split_retina.anchors <- FindIntegrationAnchors(object.list = list(hRetina, hChoroid, hSclera, hOptic_nerve, SP), dims = 1:20)
split_retina.combined <- IntegrateData(anchorset = split_retina.anchors, dims = 1:20)
DefaultAssay(split_retina.combined) <- "integrated"
#Run the standard workflow for visualization and clustering
split_retina.combined <- ScaleData(split_retina.combined, verbose = FALSE)
split_retina.combined <- RunPCA(split_retina.combined, npcs = 30, verbose = FALSE)
#t-SNE and Clustering
split_retina.combined <- RunUMAP(split_retina.combined, reduction = "pca", dims = 1:20)
split_retina.combined <- FindNeighbors(split_retina.combined, reduction = "pca", dims = 1:20)
split_retina.combined <- FindClusters(split_retina.combined, resolution = 0.5)
#Visualization
p1 <- DimPlot(split_retina.combined, reduction = "umap", group.by = "type", pt.size = 1, label.size = 3.5)
p1
DimPlot(split_retina.combined, reduction = "umap", group.by = "type", pt.size = 1, label.size = 3.5)
#segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
split_retina.combined_cycle <- CellCycleScoring(split_retina.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay="SCT")
DimPlot(split_retina.combined_cycle, reduction = "umap", label = TRUE, pt.size = 1, label.size = 3.5)
head(split_retina.combined_cycle[[]])
# Visualize the distribution of cell cycle markers across
RidgePlot(split_retina.combined_cycle, features = c("PCNA", "MCM6", "TOP2A", "MKI67"), ncol = 2, group.by = "type")
#barchart
pbar_cycle_split <- split_retina.combined_cycle@meta.data %>%
  group_by(type,Phase) %>%
  dplyr::count() %>%
  group_by(type) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=type,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per sample type")
pbar_cycle_split
split_retina.combined_cycle@meta.data %>%
  group_by(type,Phase) %>%
  dplyr::count() %>%
  group_by(type) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup()

#DE analysis
DefaultAssay(split_retina.combined) <- "RNA"
head(split_retina.combined[[]])
levels(split_retina.combined)
Idents(split_retina.combined) <- split_retina.combined@meta.data$type
#FindAllMarkers
Eye_RB_DEG <- FindAllMarkers(split_retina.combined, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
top10eye <- Eye_RB_DEG %>%
  dplyr::filter(p_val_adj < 0.01) %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(p_val_adj, desc(avg_log2FC)) %>%
  dplyr::slice_head(n = 10)
Eye_RB.DEheatmap <- DoHeatmap(split_retina.combined, slot = "counts", features = top10eye$gene)
Eye_RB.DEheatmap

#S4A
figureS4A <- DimPlot(split_retina.combined, reduction = "umap", group.by = "type", pt.size = 1.5) + 
  theme(legend.text = element_text(size = 20)) + 
  guides(color = guide_legend(override.aes = list(size = 5)))
figureS4A

#S4B
figureS4B <- DimPlot(split_retina.combined_cycle, reduction = "umap", pt.size = 1.5) + 
  theme(legend.text = element_text(size = 20)) + 
  guides(color = guide_legend(override.aes = list(size = 5)))
figureS4B

#S4C
figureS4C <- split_retina.combined_cycle@meta.data %>%
  group_by(type,Phase) %>%
  dplyr::count() %>%
  group_by(type) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=type,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per type")+ 
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    legend.position = "right",
    strip.text = element_text(size = 20),
    plot.title = element_text(size = 20, face = "bold")
  )
figureS4C

#S4D
figureS4D <- DoHeatmap(split_retina.combined, slot = "counts", features = top10eye$gene, angle = 45) + 
  theme(legend.position = "right", legend.title = element_text(size = 30), legend.text = element_text(size = 18)) +
  guides(fill = guide_colourbar(title = "Expression"), color = "none")
figureS4D

#S4E-F
figureS4.1 <- SP@meta.data %>%
  filter(section %in% paste0("section_", 1:4)) %>%
  group_by(seurat_clusters, Phase) %>%
  summarise(n = n(), .groups = "drop") %>% 
  group_by(seurat_clusters) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = seurat_clusters, y = percent, fill = Phase)) +
  geom_col() +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  ) +
  labs(
    x = "Cluster",  
    y = "Percent"
  )

figureS4.1

figureS4.1.2 <- SP@meta.data %>%
  filter(section %in% paste0("section_", 5:6)) %>%
  group_by(seurat_clusters, Phase) %>%
  summarise(n = n(), .groups = "drop") %>% 
  group_by(seurat_clusters) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = seurat_clusters, y = percent, fill = Phase)) +
  geom_col() +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    strip.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20)
  ) +
  labs(
    x = "Cluster",  
    y = "Percent"
  )

figureS4.1.2

stats_S4.1 <- SP@meta.data %>%
  filter(section %in% paste0("section_", 1:4)) %>%
  group_by(seurat_clusters, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

stats_S4.1

stats_S4.1_overall <- SP@meta.data %>%
  filter(section %in% paste0("section_", 1:4)) %>%
  group_by(Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(percent = 100 * n / sum(n))

stats_S4.1_overall

stats_S4.1.2 <- SP@meta.data %>%
  filter(section %in% paste0("section_", 5:6)) %>%
  group_by(seurat_clusters, Phase) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

stats_S4.1.2



##SCENIC##
RB <-  SP
DimPlot(RB, reduction = "umap.harmony", label = TRUE, pt.size = 1.8, label.size = 10)
Idents(RB)

exprMat <- as.matrix(RB@assays$RNA@data)
dim(exprMat)
exprMat[1:4,1:4]
cellInfo <- RB@meta.data[,c(14,2,3)]
cellInfo
colnames(cellInfo) = c('Harmony_Cluster', 'nGene', 'uUMI')
head(cellInfo)
table(cellInfo$Harmony_Cluster)

#initialize settings
library(SCENIC)

#input the reference model
org <- "hgnc"
dbDir <- "cisTarget_databases"
dbDir <- path.expand(dbDir)
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
dbs
lapply(dbs,function(x) file.exists(file.path(dbDir, x)))
myDatasetTitle <- "SCENIC example on Retinoblastoma"
dbs[dbs == "hg19-tss-centered-10kb-7species.mc9nr.feather"] <- "hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
dbs[dbs == "hg19-500bp-upstream-7species.mc9nr.feather"] <- "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
dbs
lapply(dbs,function(x) file.exists(file.path(dbDir, x)))
motifAnnotations_hgnc <- motifAnnotations
scenicOptions <- initializeScenic(org="hgnc", dbDir="cisTarget_databases", datasetTitle=myDatasetTitle, dbs=dbs, nCores=16)
saveRDS(scenicOptions, file = "int/scenicOptions.RDS")
scenicOptions@settings$nCores
closeAllConnections()
library(doParallel)
library(foreach)
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
dim(exprMat)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status

#use the AUC shiny app to determine the threshold (optional)
#aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
#savedSelections <- shiny::runApp(aucellApp)
#newThresholds <- savedSelections$thresholds
#scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
#saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))

#scenicOptions@settings$devType="png"
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
#tsneAUC(scenicOptions, aucType="AUC") # choose settings

#heatmap by clusters, all regulons
regulons <- loadInt(scenicOptions, "regulons")
regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$Harmony_Cluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")

#calculate the top regulons
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
write.csv(topRegulators, "Top_Regulon_clusters.csv")

#heatmap using top5 regulons per cluster
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$Harmony_Cluster), 
                                     function(cells) rowMeans(getAUC(regulonAUC)[, cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = TRUE, scale = TRUE))
top10_by_cluster <- topRegulators %>%
  group_by(CellType) %>%
  top_n(10, wt = RelativeActivity) %>%
  pull(Regulon)

top10_by_cluster <- unique(top10_by_cluster)
filtered_data <- regulonActivity_byCellType_Scaled[rownames(regulonActivity_byCellType_Scaled) %in% top10_by_cluster, ]
cluster_order <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")

ComplexHeatmap::Heatmap(filtered_data, 
                        name = "Regulon activity", 
                        row_title = "Top Regulons", 
                        column_title = "Cell Types",
                        column_order = cluster_order)

#save the results
saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo"))
export2loom(scenicOptions, exprMat)

#dimensional reduction of the regulons activity
nPcs <- c(5,15,50)
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=T, cex=.5)
par(mfrow=c(3,3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_oHC_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, cex=.5)
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 32
scenicOptions@settings$defaultTsne$perpl <- 09
scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 50
scenicOptions@settings$defaultTsne$perpl <- 9
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#save the result model
library(SCopeLoomR)
exprMat <- get_dgem(open_loom(loomPath))
dgem <- exprMat
head(colnames(dgem))
scenicOptions@fileNames$output["loomFile",] <- "output/RB_SCENIC.loom"
export2loom(scenicOptions, exprMat, addAllTsnes = T, hierarchy=c("SCENIC", "RB"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

##figure 5A##
figure5.1 <- ComplexHeatmap::Heatmap(t(filtered_data), 
                                      name = "Regulon activity", 
                                      row_title = "Cell Types", 
                                      column_title = "Top Regulons",
                                      row_order = cluster_order)
figure5.1

##neighbourhood analysis##
##figure 5B-C##
#Define neighbor region and following analysis
#C4
SP <- SetIdent(SP, value = "seurat_clusters_harmony")
SP <- RegionNeighbours(SP, id = "4", keep.within.id = T, verbose = TRUE)
#DE analysis of the cluster border
SP <- SetIdent(SP, value = "nbs_4")
nbs_4.markers <- FindMarkers(SP, ident.1 = "4", ident.2 = "nbs_4")
nbs_4.markers$gene <- rownames(nbs_4.markers)
RB.subset4 <- SubsetSTData(SP, expression = nbs_4 %in% c("4", "nbs_4"))
C4nbs_sorted.marks <- nbs_4.markers %>% arrange(-avg_log2FC) %>% top_n(n = 40, wt = abs(avg_log2FC))
#C1
FeatureOverlay(SP, features = "seurat_clusters_harmony", sampleids = c(1,2,3,4), ncols = 2, pt.size = 3)
SP <- SetIdent(SP, value = "seurat_clusters_harmony")
SP <- RegionNeighbours(SP, id = 1, keep.within.id = T, verbose = TRUE)
#DE analysis of the cluster border
SP <- SetIdent(SP, value = "nbs_1")
nbs_1.markers <- FindMarkers(SP, ident.1 = "1", ident.2 = "nbs_1")
nbs_1.markers$gene <- rownames(nbs_1.markers)
RB.subset1 <- SubsetSTData(SP, expression = nbs_1 %in% c("1", "nbs_1"))
C1nbs_sorted.marks <- nbs_1.markers %>% arrange(-avg_log2FC) %>% top_n(n = 40, wt = abs(avg_log2FC))

##figure 5B-C##
FeatureOverlay(RB.subset4, features = "nbs_4", ncols = 2, 
               sampleids = c(4, 6), cols = c("lightgray", "red"), pt.size = 3)
FeatureOverlay(RB.subset1, features = "nbs_1", ncols = 2, 
               sampleids = c(4, 6), cols = c("lightgray", "red"), pt.size = 3)

##figure 5D-E##
nbs1_list
nbs_1.markers
nbs_4.markers
EnhancedVolcano(nbs_1.markers,
                lab = rownames(nbs_1.markers), 
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Significant genes C1 neighbor vs C1',
                pCutoff = 0.1,
                FCcutoff = 0.58,
                pointSize = 4,
                labSize = 7,
                colAlpha = 1,
                legendLabels = c('Not Sig','Log2FC','padj','padj&Log2FC'),
                legendPosition = "bottom",
                legendLabSize = 16,
                legendIconSize = 5.0,
                drawConnectors = T,
                widthConnectors = 0.75)
EnhancedVolcano(nbs_4.markers,
                lab = rownames(nbs_4.markers), 
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Significant genes C4 neighbor vs C4',
                pCutoff = 0.1,
                FCcutoff = 0.58,
                pointSize = 4,
                labSize = 7,
                colAlpha = 1,
                legendLabels = c('Not Sig','Log2FC','padj','padj&Log2FC'),
                legendPosition = "bottom",
                legendLabSize = 16,
                legendIconSize = 5.0,
                drawConnectors = T,
                widthConnectors = 0.75)

#GO term enrichment of the genes
original_nbs4_list <- nbs_4.markers$avg_log2FC
names(original_nbs4_list) <- rownames(nbs_4.markers)
nbs4_list <- na.omit(original_nbs4_list)
nbs4_list = sort(nbs4_list, decreasing = TRUE)

original_nbs1_list <- nbs_1.markers$avg_log2FC
names(original_nbs1_list) <- rownames(nbs_1.markers)
nbs1_list <- na.omit(original_nbs1_list)
nbs1_list = sort(nbs1_list, decreasing = TRUE)

nbs1gse <- gseGO(geneList=nbs1_list, 
                 ont ="BP", 
                 keyType = "SYMBOL", 
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "none")
nbs4gse <- gseGO(geneList=nbs4_list, 
                 ont ="BP", 
                 keyType = "SYMBOL", 
                 nPerm = 10000, 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = organism, 
                 pAdjustMethod = "none")

pCPnbs1D <- dotplot(nbs1gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCPnbs1C <- cnetplot(nbs1gse, categorySize="pvalue", foldChange=nbs1_list, showCategory = 3)
pCPnbs4D <- dotplot(nbs4gse, showCategory=5, split=".sign") + facet_grid(.~.sign)
pCPnbs4C <- cnetplot(nbs4gse, categorySize="pvalue", foldChange=nbs4_list, showCategory = 3)

##CellChat##
#read the datasets
# Prepare input data for CellChat analysis
new.cluster.ids <- c("C0", "C1", "C2", "C3", "C4", 
                     "C5",
                     "C6", "C7", "C8")
names(new.cluster.ids) <- levels(SP)
SP <- RenameIdents(SP, new.cluster.ids)
DimPlot(SP, reduction = "umap.harmony", label = TRUE, pt.size = 1.3, label.size = 4.5)
data.input = GetAssayData(SP, slot = "data", assay = "SCT")
meta = data.frame(labels = Idents(SP), row.names = names(Idents(SP)))
###use the whole transcriptomic profile for cell chat
cellchat.all <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat.all <- addMeta(cellchat.all, meta = meta)
cellchat.all <- setIdent(cellchat.all, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat.all@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat.all@idents)) # number of cells in each cell group
##set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat.all@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat.all <- subsetData(cellchat.all)
future::plan("multisession", workers = 14)
#pre-processing of the expression data for cell-cell communication
cellchat.all <- identifyOverExpressedGenes(cellchat.all)
cellchat.all <- identifyOverExpressedInteractions(cellchat.all)
#The communication probability and infer cellular communication network
cellchat.all <- computeCommunProb(cellchat.all)
cellchat.all <- filterCommunication(cellchat.all, min.cells = 10)
cellchat.all <- computeCommunProbPathway(cellchat.all)
cellchat.all <- aggregateNet(cellchat.all)
saveRDS(cellchat.all, file = "cellchat.rds")

##figure 6A##
groupSize <- as.numeric(table(cellchat.all@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat.all@net$count, vertex.weight = rowSums(cellchat.all@net$count), 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions", edge.width.max = 14)
netVisual_circle(cellchat.all@net$weight, vertex.weight = rowSums(cellchat.all@net$weight), 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", edge.width.max = 14)

#Manifold and classification learning analysis of signaling networks
#Identify signaling groups based on their functional similarity
cellchat.all <- computeNetSimilarity(cellchat.all, type = "functional")
cellchat.all <- netEmbedding(cellchat.all, slot.name = "netP", type = "functional")

#Manifold learning of the signaling networks for a single dataset
cellchat.all <- netClustering(cellchat.all, type = "functional")
#Classification learning of the signaling networks for a single dataset
#Visualization in 2D-space
pfun <- netVisual_embedding(cellchat.all, type = "functional", label.size = 3.5)
pfunzoom <- netVisual_embeddingZoomIn(cellchat.all, type = "functional", nCol = 2)
pfun
pfunzoom

#Manifold and classification learning analysis of signaling networks
#Identify signaling groups based on their structural similarity
cellchat.all <- computeNetSimilarity(cellchat.all, type = "structural")
cellchat.all <- netEmbedding(cellchat.all, slot.name = "netP", type = "structural")

#Manifold learning of the signaling networks for a single dataset
cellchat.all <- netClustering(cellchat.all, type = "structural")
#Classification learning of the signaling networks for a single dataset
#Visualization in 2D-space
pfun.2 <- netVisual_embedding(cellchat.all, type = "structural", label.size = 3.5)
pfunzoom.2 <- netVisual_embeddingZoomIn(cellchat.all, type = "structural", nCol = 2)
pfun.2
pfunzoom.2

##figure 6B##
figure6.3 <- netVisual_embedding(cellchat.all, type = "functional", label.size = 10, dot.size = c(4, 12)) + theme(
  axis.title = element_text(size = 20), 
  axis.text = element_text(size = 20),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 20) 
)
figure6.3

##figure 6C##
#CypA
pathways.show <- c("CypA")
p1_chord <- netVisual_aggregate(cellchat.all, signaling = pathways.show, layout = "chord", vertex.label.cex = 1.5)
p1_heatmap <- grid.grabExpr(netAnalysis_signalingRole_network(cellchat.all, signaling = "CypA", width = 8, height = 2.5, font.size = 10))

#MK
pathways.show <- c("MK")
p2_chord <- netVisual_aggregate(cellchat.all, signaling = pathways.show, layout = "chord", vertex.label.cex = 1.5)
p2_heatmap <- grid.grabExpr(netAnalysis_signalingRole_network(cellchat.all, signaling = "MK", width = 8, height = 2.5, font.size = 10))

#PECAM1
pathways.show <- c("PECAM1")
p3_chord <- netVisual_aggregate(cellchat.all, signaling = pathways.show, layout = "chord", vertex.label.cex = 1.5)
p3_heatmap <- grid.grabExpr(netAnalysis_signalingRole_network(cellchat.all, signaling = "PECAM1", width = 8, height = 2.5, font.size = 10))

#SPP1
pathways.show <- c("SPP1")
p4_chord <- netVisual_aggregate(cellchat.all, signaling = pathways.show, layout = "chord", vertex.label.cex = 1.5)
p4_heatmap <- grid.grabExpr(netAnalysis_signalingRole_network(cellchat.all, signaling = "SPP1", width = 8, height = 2.5, font.size = 10))

pCypA <- plot_grid(p1_chord, p1_heatmap, ncol = 2, rel_widths = c(1, 1))
pMK <- plot_grid(p2_chord, p2_heatmap, ncol = 2, rel_widths = c(1, 1))
pPECAM1 <- plot_grid(p3_chord, p3_heatmap, ncol = 2, rel_widths = c(1, 1))
pSPP1 <- plot_grid(p4_chord, p4_heatmap, ncol = 2, rel_widths = c(1, 1))

pCypA
pMK
pPECAM1
pSPP1

# SPP1 → Group1
# MK → Group2
# CypA → Group3
# PECAM1 → Group4
figure6.5
p1_chord
p2_chord
p3_chord
p4_chord

#show the figures for all pathways
all_pathways <- c(
  "VEGF", "APP", "NECTIN", "MIF", "NOTCH", "COLLAGEN", "SPP1", "IGFBP", "FN1", "MPZ",
  "PTN", "LAMININ", "JAM", "MK", "GAP", "CADM",
  "CDH", "RBP4", "NCAM", "CypA", "CD99",
  "ANGPT", "GRN", "HSPG", "PECAM1", "AGRN", "CDH5", "ESAM"
)

plot_all_pathways <- function(pathway_list, cellchat_obj, ncol = 3) {
  plots <- list()
  for (p in pathway_list) {
    try({
      p_chord <- netVisual_aggregate(cellchat_obj, signaling = p, layout = "chord")
      p_heat <- grid.grabExpr(netAnalysis_signalingRole_network(cellchat_obj, signaling = p, width = 8, height = 2.5, font.size = 10))
      combined <- plot_grid(p_chord, p_heat, ncol = 2, rel_widths = c(1.2, 1))
      plots[[p]] <- combined
    }, silent = TRUE)
  }
  all_plot <- plot_grid(plotlist = plots, ncol = ncol)
  return(all_plot)
}

final_plot <- plot_all_pathways(all_pathways, cellchat.all, ncol = 3)

num_rows <- ceiling(length(all_pathways) / 3)

#summarize the intercellular communication resulting from the functional similarity
group_assignments <- cellchat.all@netP$similarity$functional$group$single
group_list <- split(names(group_assignments), group_assignments)

get_avg_roles <- function(pathways, centr_list) {
  cluster_names <- names(centr_list[[1]]$outdeg)
  send_mat <- matrix(0, nrow = length(pathways), ncol = length(cluster_names), dimnames = list(pathways, cluster_names))
  recv_mat <- matrix(0, nrow = length(pathways), ncol = length(cluster_names), dimnames = list(pathways, cluster_names))
  
  for (path in pathways) {
    if (!is.null(centr_list[[path]])) {
      send_mat[path, ] <- centr_list[[path]]$outdeg
      recv_mat[path, ] <- centr_list[[path]]$indeg
    }
  }
  
  avg_send <- colMeans(send_mat, na.rm = TRUE)
  avg_recv <- colMeans(recv_mat, na.rm = TRUE)
  return(list(sender = avg_send, receiver = avg_recv))
}

group_avg_roles <- lapply(group_list, get_avg_roles, centr_list = cellchat.all@netP$centr)

df_sender <- map2_df(group_avg_roles, names(group_avg_roles), ~{
  tibble(Cluster = names(.x$sender), Strength = .x$sender, Role = "Sender", Group = .y)
})

df_receiver <- map2_df(group_avg_roles, names(group_avg_roles), ~{
  tibble(Cluster = names(.x$receiver), Strength = .x$receiver, Role = "Receiver", Group = .y)
})

df_all <- bind_rows(df_sender, df_receiver)

df_top <- df_all %>%
  group_by(Group, Role) %>%
  slice_max(order_by = Strength, n = 3, with_ties = FALSE) %>%
  ungroup()
cluster_levels <- paste0("C", 0:8)

df_top <- df_top %>%
  mutate(
    Cluster_sorted = factor(Cluster, levels = rev(cluster_levels)), 
    FacetID = paste0("Group ", Group, " - ", Role)
  )

top3_functional_cellchat <- ggplot(df_top, aes(x = Cluster_sorted, y = Strength, fill = Group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ FacetID, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(title = "Top 3 Sender and Receiver Clusters per Group",
       x = "Cluster", y = "Average Strength") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )

print(top3_functional_cellchat)

df_top <- df_all %>%
  group_by(Group, Role) %>%
  slice_max(order_by = Strength, n = 5, with_ties = FALSE) %>%
  ungroup()
cluster_levels <- paste0("C", 0:8)

df_top <- df_top %>%
  mutate(
    Cluster_sorted = factor(Cluster, levels = rev(cluster_levels)), 
    FacetID = paste0("Group ", Group, " - ", Role)
  )

top5_functional_cellchat <- ggplot(df_top, aes(x = Cluster_sorted, y = Strength, fill = Group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ FacetID, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(title = "Top 5 Sender and Receiver Clusters per Group",
       x = "Cluster", y = "Average Strength") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )

print(top5_functional_cellchat)

top3_functional_cellchat
top5_functional_cellchat

#Identify signals contributing most to outgoing or incoming signaling of certain cell groups
#Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.all, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.all, pattern = "incoming")
ht1 + ht2
##figure 6D##
figure6.1.1 <- netAnalysis_signalingRole_heatmap(cellchat.all, pattern = "outgoing", font.size = 20, font.size.title = 20, width = 24, height = 20)
figure6.1.2 <- netAnalysis_signalingRole_heatmap(cellchat.all, pattern = "incoming", font.size = 20, font.size.title = 20, width = 24, height = 10)
figure6.1 <- figure6.1.1 + figure6.1.2
figure6.1

#show all the significant interactions (L-R pairs) associated with certain signaling pathways
##figure 6E-F##
figure6.2.1 <- netVisual_bubble(cellchat.all, sources.use = c("C0", "C2", "C3", "C4", "C5", "C6", "C7", "C8"), targets.use = "C1", remove.isolate = FALSE, dot.size.min = 12) + theme(
  axis.title = element_text(size = 20), 
  axis.text = element_text(size = 20),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 20) 
)
figure6.2.2 <- netVisual_bubble(cellchat.all, sources.use = "C4", targets.use = c("C0", "C1", "C2", "C3", "C5", "C6", "C7", "C8"), remove.isolate = FALSE, dot.size.min = 12) + theme(
  axis.title = element_text(size = 20), 
  axis.text = element_text(size = 20),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 20) 
)

figure6.2 <- figure6.2.1 | figure6.2.2
figure6.2

##supplementary figure 7##
#summarize the intercellular communication resulting from the functional similarity
group_assignments <- cellchat.all@netP$similarity$functional$group$single
group_list <- split(names(group_assignments), group_assignments)

get_avg_roles <- function(pathways, centr_list) {
  cluster_names <- names(centr_list[[1]]$outdeg)
  send_mat <- matrix(0, nrow = length(pathways), ncol = length(cluster_names), dimnames = list(pathways, cluster_names))
  recv_mat <- matrix(0, nrow = length(pathways), ncol = length(cluster_names), dimnames = list(pathways, cluster_names))
  
  for (path in pathways) {
    if (!is.null(centr_list[[path]])) {
      send_mat[path, ] <- centr_list[[path]]$outdeg
      recv_mat[path, ] <- centr_list[[path]]$indeg
    }
  }
  
  avg_send <- colMeans(send_mat, na.rm = TRUE)
  avg_recv <- colMeans(recv_mat, na.rm = TRUE)
  return(list(sender = avg_send, receiver = avg_recv))
}

group_avg_roles <- lapply(group_list, get_avg_roles, centr_list = cellchat.all@netP$centr)

df_sender <- map2_df(group_avg_roles, names(group_avg_roles), ~{
  tibble(Cluster = names(.x$sender), Strength = .x$sender, Role = "Sender", Group = .y)
})

df_receiver <- map2_df(group_avg_roles, names(group_avg_roles), ~{
  tibble(Cluster = names(.x$receiver), Strength = .x$receiver, Role = "Receiver", Group = .y)
})

df_all <- bind_rows(df_sender, df_receiver)

df_top <- df_all %>%
  group_by(Group, Role) %>%
  slice_max(order_by = Strength, n = 3, with_ties = FALSE) %>%
  ungroup()
cluster_levels <- paste0("C", 0:8)

df_top <- df_top %>%
  mutate(
    Cluster_sorted = factor(Cluster, levels = rev(cluster_levels)), 
    FacetID = paste0("Group ", Group, " - ", Role)
  )

figureS7.2 <- ggplot(df_top, aes(x = Cluster_sorted, y = Strength, fill = Group)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ FacetID, scales = "free_y", ncol = 2) +
  theme_bw() +
  labs(title = "Top 3 Sender and Receiver Clusters per Group",
       x = "Cluster", y = "Average Strength") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  )

figureS7.2

##supplementary figure 8##
figureS8.1 <- netVisual_bubble(cellchat.all, sources.use = "C1", targets.use = c("C0", "C2", "C3", "C4", "C5", "C6", "C7", "C8"), remove.isolate = FALSE, dot.size.min = 12) + theme(
  axis.title = element_text(size = 20), 
  axis.text = element_text(size = 20),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 20) 
)
figureS8.2 <- netVisual_bubble(cellchat.all, sources.use = c("C0", "C1", "C2", "C3", "C5", "C6", "C7", "C8"), targets.use = "C4", remove.isolate = FALSE, dot.size.min = 12) + theme(
  axis.title = element_text(size = 20), 
  axis.text = element_text(size = 20),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 20) 
)

figureS8 <- figureS8.1 | figureS8.2
figureS8

##cell type enrichment algorithm interpret spatial deconvolution results
##Table 1##
#for enriched cell type
ct_visualize <- c("CP", "MKI67..CP", "Glial", "CAF")
card_files <- c("RB1A_CARD.CSV", "RB1B_CARD.CSV", "RB1C_CARD.CSV", 
                "RB1D_CARD.CSV", "RB2A_CARD.CSV", "RB2B_CARD.CSV")
enriched_spots_list <- list()

for (file in card_files) {
  
  card_data <- read.csv(file)
  colnames(card_data)[1] <- "barcode"
  
  if (!all(ct_visualize %in% colnames(card_data))) {
    warning(paste("Missing columns in file:", file))
    next
  }
  
  section_name <- gsub("_CARD.CSV", "", file)
  card_data$section <- section_name
  
  norm_data <- card_data %>%
    rowwise() %>%
    mutate(across(all_of(ct_visualize), ~ . / sum(c_across(all_of(ct_visualize))), 
                  .names = "norm_{.col}")) %>%
    ungroup()
  
  avg_props <- norm_data %>%
    summarise(across(starts_with("norm_"), mean, na.rm = TRUE))
  
  n_spots <- nrow(norm_data)
  
  celltype_counts <- round(avg_props * n_spots)
  celltype_counts <- as.numeric(celltype_counts)
  names(celltype_counts) <- gsub("norm_", "", colnames(avg_props))
  
  selected_spots <- list()
  
  for (ct in names(celltype_counts)) {
    
    temp <- norm_data %>%
      arrange(desc(.data[[paste0("norm_", ct)]]))
    
    top_spots <- head(temp$barcode, celltype_counts[ct])
    selected_spots[[ct]] <- top_spots
  }
  
  df_selected <- data.frame(
    barcode = unlist(selected_spots),
    celltype = rep(names(selected_spots), times = sapply(selected_spots, length)),
    section = section_name,
    stringsAsFactors = FALSE
  )
  
  duplicated_barcodes <- df_selected$barcode[duplicated(df_selected$barcode)]
  
  if (length(duplicated_barcodes) > 0) {
    message("Duplicates found in section: ", section_name)
    df_selected <- df_selected %>%
      filter(!barcode %in% duplicated_barcodes)
  }
  
  enriched_spots_list[[section_name]] <- df_selected
}

enriched_spots_all <- bind_rows(enriched_spots_list)
write.csv(enriched_spots_all, "Cell_Type_Enriched_Spots.csv", row.names = FALSE)

#calculate the proportion of different cell types in each cluster
spatial_cluster <- read.csv("spatial_cluster.CSV")
colnames(spatial_cluster)[1] <- "barcode"
colnames(spatial_cluster)[2] <- "cluster"
enriched_with_cluster <- enriched_spots_all %>%
  left_join(spatial_cluster, by = "barcode")
cluster_enrichment_summary <- enriched_with_cluster %>%
  group_by(cluster, celltype) %>%
  summarise(n_spots = n()) %>%
  group_by(cluster) %>%
  mutate(total = sum(n_spots),
         proportion = n_spots / total) %>%
  arrange(cluster, desc(proportion))

print(cluster_enrichment_summary)
write.csv(cluster_enrichment_summary, "cluster_enrichment_summary.CSV", row.names = FALSE)


#table S3
#generate a list of regulons for the paper
topRegulators <- read.csv("Top_Regulon_clusters.csv")
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$Harmony_Cluster), 
                                     function(cells) rowMeans(getAUC(regulonAUC)[, cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = TRUE, scale = TRUE))
top10_by_cluster <- topRegulators %>%
  group_by(CellType) %>%
  top_n(10, wt = RelativeActivity) %>%
  pull(Regulon)

top10_by_cluster <- unique(top10_by_cluster)
filtered_data <- regulonActivity_byCellType_Scaled[rownames(regulonActivity_byCellType_Scaled) %in% top10_by_cluster, ]
cluster_order <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")

selected_regulon_paper <- data.frame(Regulons = rownames(filtered_data))
write.csv(selected_regulon_paper, "selected_regulon_paper.csv", row.names = FALSE)
print(head(selected_regulon_paper))

regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")


#load the information
selected_regulon <- read.csv("selected_regulon_paper.csv")
original_names <- selected_regulon$Regulons
selected_tf <- gsub(" \\(.*\\)", "", selected_regulon$Regulons)
selected_tf <- gsub("_extended", "", selected_tf)
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")

#extract the information
high_confident_tf <- regulonTargetsInfo %>%
  filter(TF %in% selected_tf & highConfAnnot == TRUE) %>%
  select(TF, gene, highConfAnnot, nMotifs, bestMotif, NES)

low_confident_tf <- regulonTargetsInfo %>%
  filter(TF %in% selected_tf & highConfAnnot == FALSE) %>%
  select(TF, gene, highConfAnnot, nMotifs, bestMotif, NES)

enriched_genes <- regulonTargetsInfo %>%
  filter(TF %in% selected_tf) %>%
  select(TF, gene, highConfAnnot, nMotifs)

#match with the selected list of regulons in the heatmap
result_table <- data.frame(
  Regulon_Set = character(),
  High_confident_TF = character(),
  Low_confident_TF = character(),
  Enriched_Genes = character(),
  stringsAsFactors = FALSE
)


for (i in seq_along(selected_tf)) {
  tf <- selected_tf[i]
  original_name <- original_names[i]
  
  high_confident <- regulonTargetsInfo %>%
    filter(TF == tf & highConfAnnot == TRUE) %>%
    pull(gene) %>%
    unique() %>%
    paste(collapse = "; ")
  
  low_confident <- regulonTargetsInfo %>%
    filter(TF == tf & highConfAnnot == FALSE) %>%
    pull(gene) %>%
    unique() %>%
    paste(collapse = "; ")
  
  enriched_genes <- regulonTargetsInfo %>%
    filter(TF == tf) %>%
    pull(gene) %>%
    unique() %>%
    paste(collapse = "; ")
  
  result_table <- rbind(result_table, data.frame(
    Regulon_Set = original_name,  
    High_confident_TF = high_confident,
    Low_confident_TF = low_confident,
    Enriched_Genes = enriched_genes,
    stringsAsFactors = FALSE
  ))
}

print(result_table)

write.csv(result_table, "Regulon_Information_Table.csv", row.names = FALSE)


###prapare for RNA velocity
#metadata
SP$barcode <- colnames(SP)
SP$UMAP_1 <- SP@reductions$umap.harmony@cell.embeddings[,1]
SP$UMAP_2 <- SP@reductions$umap.harmony@cell.embeddings[,2]
write.csv(SP@meta.data, file='metadata.csv', quote=F, row.names=F)
#count matrix
counts_matrix <- GetAssayData(SP, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0('counts.mtx'))
#reduction
write.csv(SP@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)
#gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)



######save session info#####
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


#####end of R script#####


