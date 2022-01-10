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
library(hdf5r)
library(rhdf5)
library(scCustomize)
library(patchwork)
library(viridis)
library(Nebulosa)



setwd("~/Desktop/lanser_nucseq/hpa_rds_files")

#Load objects
HPA_PBS_merged <- readRDS("HPA_PBS_merged_unNormalized.rds")
WT_KO_merged <- readRDS("WT_KO_merged_unNormalized.rds")
KO_HPA_merged <- readRDS("KO_HPA_merged_unNormalized.rds")



#Cluster HPS & PBS conditions
HPA_PBS_merged <- SCTransform(HPA_PBS_merged, method = "glmGamPoi")
###PCA
HPA_PBS_merged <- RunPCA(object = HPA_PBS_merged)
ElbowPlot(object = HPA_PBS_merged, ndims = 50)
###
HPA_PBS_merged <- FindNeighbors(object = HPA_PBS_merged, dims = 1:40)
HPA_PBS_merged <- FindClusters(object = HPA_PBS_merged, resolution = 0.7)
HPA_PBS_merged <- RunUMAP(HPA_PBS_merged, dims = 1:40)


DimPlot_scCustom(HPA_PBS_merged, label = T, split.by = "condition")
ggsave("PBS_HPA_Condition_Plot_res0.7_dims40.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)

saveRDS(HPA_PBS_merged, "HPA_PBS_merged_postUMAP.rds")


tibble(
  cluster = HPA_PBS_merged$seurat_clusters,
  Dataset = HPA_PBS_merged$condition
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



