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


WT_PBS_neurons <- readRDS("WT_PBS_neurons_postUMAP.rds")
#Cluster HPS & PBS conditions
WT_PBS_neurons <- SCTransform(WT_PBS_neurons, method = "glmGamPoi")
###PCA
WT_PBS_neurons <- RunPCA(object = WT_PBS_neurons)
ElbowPlot(object = WT_PBS_neurons, ndims = 50)
###
WT_PBS_neurons <- FindNeighbors(object = WT_PBS_neurons, dims = 1:35)
WT_PBS_neurons <- FindClusters(object = WT_PBS_neurons, resolution = 0.8)
WT_PBS_neurons <- RunUMAP(WT_PBS_neurons, dims = 1:35)


DimPlot_scCustom(WT_PBS_neurons, label = T, split.by = "condition")

ggsave("WT_PBS_Condition_Plot_res0.7_dims40.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)

FeaturePlot_scCustom(WT_PBS_neurons, features = "Grik1", label = T, raster = F)

saveRDS(WT_PBS_neurons, "WT_PBS_neurons_postUMAP.rds")


tibble(
  cluster = WT_PBS_neurons$seurat_clusters,
  Dataset = WT_PBS_neurons$condition
) %>%
  group_by(cluster, Dataset) %>%
  count() %>%
  group_by(cluster) %>%
  mutate(
    percent=(100*n)/sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    cluster=paste("Cluster", cluster)
  ) %>%
  ggplot(aes(x="",y=percent, fill=Dataset)) +
  geom_col(width=1) +
  coord_polar("y", start=0) +
  facet_wrap(vars(cluster)) +
  theme(axis.text.x=element_blank()) +
  xlab(NULL) +
  ylab(NULL)



