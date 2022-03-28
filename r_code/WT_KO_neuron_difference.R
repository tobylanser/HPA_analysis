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
WT_KO_merged <- readRDS("WT_KO_merged_postUMAP_filtered.rds")
WT_KO_neurons <- readRDS("WT_KO_neurons_postUMAP.rds")

clust5 <- FindMarkers(WT_KO_merged, ident.1 = c("5"), only.pos = T)

Idents(WT_KO_neurons) <- "condition"
neuron_condition_Markers <- FindAllMarkers(WT_KO_neurons, only.pos = T)
neuron_condition_Markers_top20 <- neuron_condition_Markers %>%
  group_by(cluster) %>% 
  slice_max(n = 20, avg_log2FC)

write.csv(clust5_markers, "/Users/WeinerLab/Desktop/lanser_nucseq/Markers/KO_WT_difference_neuron_markers")

DimPlot_scCustom(WT_KO_merged, label = T)

FeaturePlot(WT_KO_merged, features = clust3_DEG, label = T)

neurons_diff <- subset(WT_KO_merged, idents = "5")

neurons_diff <- SCTransform(neurons_diff, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
neurons_diff <- RunPCA(object = neurons_diff)
ElbowPlot(object = neurons_diff, ndims = 50)
###  
neurons_diff <- FindNeighbors(object = neurons_diff, dims = 1:10)
neurons_diff <- FindClusters(object = neurons_diff, resolution = 0.5)
neurons_diff <- RunUMAP(neurons_diff, dims = 1:10)

DimPlot_scCustom(WT_KO_neurons, label = T, split.by = "condition")
joint1 <- Plot_Density_Joint_Only(WT_KO_neurons, features = c("Taco1", "Inpp4b"))

joint2 <- Plot_Density_Joint_Only(WT_KO_neurons, features = c("Stard9", "Nefm"))
joint1/joint2
Plot_Density_Custom(WT_KO_neurons, features = "Er81")
FeaturePlot_scCustom(WT_KO_neurons, features = "Pthlh", label = T)


DoHeatmap(WT_KO_merged, group.by = "seurat_clusters", features = neuron_condition_Markers$gene)
ggsave("/Users/WeinerLab/Desktop/lanser_nucseq/Cluster_Plots/DimPlot/WT_KO_HeatMap.jpeg", device = "jpeg", 
       width = 18, height = 12, units = "in", dpi = 600)


