library(patchwork)
library(cluster)
library(cowplot)
library(Rtsne)
library(scatterplot3d)
library(plotly)
library(ggsci)
library(ggrepel)
library(ggiraph)
library(gplots)
library(ggridges)
library(ggstance)
library(autoplotly)
library(RColorBrewer)
library(pheatmap)
library(treemap)
library(FateID)
library(sigora)
library(colorpatch)
library(Seurat)
library(tidyverse)
library(reticulate)
library(ggbeeswarm)
library(ggrepel)
library(ggridges)
library(RColorBrewer)
library(devtools)
library(ggplot2)
library(dplyr)
library(Matrix)
library(umap)
library(rhdf5)
library(scCustomize)
library(patchwork)
library(viridis)
library(Nebulosa)



setwd("~/Desktop/lanser_nucseq/hpa_rds_files")

#Load objects

WT_KO_merged <- readRDS("WT_KO_merged_unNormalized.rds")


WT_KO_merged_unfiltered <- readRDS("WT_KO_merged_postUMAP.rds")
Idents(WT_KO_merged) <- "seurat_clusters"
WT_KO_merged <- subset(WT_KO_merged_unfiltered, idents = c("0", "1", "5", "6", "8", "9",
                                                "10", "11", 
                                                  "19", "21", "28", "30"), invert = T)

WT_KO_merged <- subset(WT_KO_merged, idents = c("0", "1", "3", "6", "27",
                                                 "26"), invert = T)

#Cluster WT & KO conditions
WT_KO_merged <- SCTransform(WT_KO_merged, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
WT_KO_merged <- RunPCA(object = WT_KO_merged)
ElbowPlot(object = WT_KO_merged, ndims = 50)
###
WT_KO_merged <- FindNeighbors(object = WT_KO_merged, dims = 1:15)
WT_KO_merged <- FindClusters(object = WT_KO_merged, resolution = 0.7)
WT_KO_merged <- RunUMAP(WT_KO_merged, dims = 1:15)
DimPlot_scCustom(WT_KO_merged, label = T, group.by = "cell_type")

DimPlot_scCustom(WT_KO_merged, label = T, group.by = "condition") + DarkTheme()
ggsave("WT_KO_Condition_Plot_res0.8_dims30.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)

WT_KO_merged_unfiltered
FeaturePlot(WT_KO_merged_unfiltered, features = c("Grik1"), label = T)

saveRDS(WT_KO_merged, "WT_KO_merged_postUMAP_filtered.rds")

tibble(
  cluster = WT_KO_merged$cell_type,
  Dataset = WT_KO_merged$condition
) %>%
  group_by(cluster, Dataset) %>%
  count() %>%
  group_by(cluster) %>%
  mutate(
    percent=(100*n)/sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    cluster=paste("Cell Type:", cluster)
  ) %>%
  ggplot(aes(x="",y=percent, fill=Dataset)) +
  geom_col(width=1) +
  coord_polar("y", start=0) +
  facet_wrap(vars(cluster)) +
  theme(axis.text.x=element_blank()) +
  xlab(NULL) +
  ylab(NULL)

WT_KO_merged@meta.data %>%
  group_by(cell_type, condition) %>%
  count() %>%
  group_by(cell_type) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=cell_type,y=percent, fill=condition)) +
  geom_col() +
  ggtitle("Percentage of KO vs WT in each cell type")

ggsave("WT_KO_PieChart_res0.8_dims30.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)

FeaturePlot_scCustom(WT_KO_merged, features = "Aldh1l1", label = T)

DimPlot_scCustom(WT_KO_merged, label = T)
#Add cell type labels 
Inhib <- subset(WT_KO_merged, idents = c("11", "4", "6", "8"))
Inhib$cell_type <- "Inhibitory Neurons"

Excit <- subset(WT_KO_merged, idents = c("5", "7", "17"))
Excit$cell_type <- "Excitatory Neurons"

Oligodendrocytes <- subset(WT_KO_merged, idents = c("0", "13"))
Oligodendrocytes$cell_type <- "Oligodendrocytes"

OPCs <- subset(WT_KO_merged, idents = c("9"))
OPCs$cell_type <- "OPCs"

Mulller <- subset(WT_KO_merged, idents = c("15"))
Mulller$cell_type <- "Muller glia?"

Astrocytes <- subset(WT_KO_merged, idents = c("2", "16"))
Astrocytes$cell_type <- "Astrocytes"

Glia <- subset(WT_KO_merged, idents = c("1", "3", "10", "12", "18"))
Glia$cell_type <- "Glia"

na <- subset(WT_KO_merged, idents = c("14"))
na$cell_type <- "NA"

WT_KO_merged <- merge(WT_KO_neurons, y = c(Oligodendrocytes, OPCs,
                                   Mulller, Astrocytes, Glia, na))


WT_KO_merged <- SCTransform(WT_KO_merged, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
WT_KO_merged <- RunPCA(object = WT_KO_merged)
ElbowPlot(object = WT_KO_merged, ndims = 50)
###
WT_KO_merged <- FindNeighbors(object = WT_KO_merged, dims = 1:15)
WT_KO_merged <- FindClusters(object = WT_KO_merged, resolution = 0.7)
WT_KO_merged <- RunUMAP(WT_KO_merged, dims = 1:15)

Idents(WT_KO_merged) <- "seurat_clusters"
DimPlot_scCustom(WT_KO_merged, label = T, group.by = "cell_type")

DimPlot_scCustom(WT_KO_merged, label = T)

saveRDS(WT_KO_merged, "WT_KO_merged_postUMAP_filtered.rds")
WT_KO_merged <- readRDS("WT_KO_merged_postUMAP_filtered.rds")


Idents(WT_KO_merged) <- "cell_type"
glia <- subset(WT_KO_merged, idents = "Glia")

glia <- SCTransform(glia, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
glia <- RunPCA(object = glia)
ElbowPlot(object = glia, ndims = 50)
###
glia <- FindNeighbors(object = glia, dims = 1:40)
glia <- FindClusters(object = glia, resolution = 0.7)
glia <- RunUMAP(glia, dims = 1:40)
DimPlot_scCustom(glia)


