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

WT_KO_merged <- readRDS("WT_KO_merged_unNormalized.rds")




#Cluster KO & WT conditions
WT_KO_merged <- SCTransform(WT_KO_merged, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
WT_KO_merged <- RunPCA(object = WT_KO_merged)
ElbowPlot(object = WT_KO_merged, ndims = 50)
###
WT_KO_merged <- FindNeighbors(object = WT_KO_merged, dims = 1:30)
WT_KO_merged <- FindClusters(object = WT_KO_merged, resolution = 0.8)
##
WT_KO_merged <- RunUMAP(WT_KO_merged, dims = 1:30)
DimPlot_scCustom(WT_KO_merged, split.by = "condition", label = T, raster = F, order = T)
ggsave("WT_KO_Condition_Plot_res0.8_dims30.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)

saveRDS(WT_KO_merged, "WT_KO_merged_postUMAP.rds")

tibble(
  cluster = WT_KO_merged$seurat_clusters,
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
    cluster=paste("Cluster", cluster)
  ) %>%
  ggplot(aes(x="",y=percent, fill=Dataset)) +
  geom_col(width=1) +
  coord_polar("y", start=0) +
  facet_wrap(vars(cluster)) +
  theme(axis.text.x=element_blank()) +
  xlab(NULL) +
  ylab(NULL)
ggsave("WT_KO_PieChart_res0.8_dims30.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)




