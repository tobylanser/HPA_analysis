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
WT_PBS_merged <- readRDS("WT_PBS_merged_unNormalized.rds")

WT_PBS_merged <- readRDS("WT_PBS_merged_postUMAP.rds")
#Cluster HPS & PBS conditions
WT_PBS_merged <- SCTransform(WT_PBS_merged, method = "glmGamPoi")
###PCA
WT_PBS_merged <- RunPCA(object = WT_PBS_merged)
ElbowPlot(object = WT_PBS_merged, ndims = 50)
###
WT_PBS_merged <- FindNeighbors(object = WT_PBS_merged, dims = 1:40)
WT_PBS_merged <- FindClusters(object = WT_PBS_merged, resolution = 1.2)
WT_PBS_merged <- RunUMAP(WT_PBS_merged, dims = 1:40)


DimPlot_scCustom(WT_PBS_merged, label = T)
ggsave("WT_PBS_Condition_Plot_res0.7_dims40.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)

FeaturePlot_scCustom(WT_PBS_merged, features = "Ttn", label = T, raster = F)

saveRDS(WT_PBS_merged, "WT_PBS_merged_postUMAP_filtered.rds")



WT_PBS_neurons <- subset(WT_PBS_merged, idents = c("5", "7", "8", "11", "13",
                                                  "14", "15", "16", "20", "21", "25",
                                                  "26", "31"))

tibble(
  cluster = WT_PBS_merged$seurat_clusters,
  Dataset = WT_PBS_merged$condition
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



