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
KO_HPA_merged <- readRDS("KO_HPA_merged_unNormalized.rds")

##Cluster KO and HPA conditions
KO_HPA_merged <- SCTransform(KO_HPA_merged, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
KO_HPA_merged <- RunPCA(object = KO_HPA_merged, verbose = F)
ElbowPlot(object = KO_HPA_merged, ndims = 50)
###
KO_HPA_merged <- FindNeighbors(object = KO_HPA_merged, dims = 1:30)
KO_HPA_merged <- FindClusters(object = KO_HPA_merged, resolution = 0.8)
KO_HPA_merged <- RunUMAP(KO_HPA_merged, dims = 1:30)

DimPlot_scCustom(KO_HPA_merged, split.by = "condition", label = T, raster = F)
ggsave("KO_HPA_Condition_Plot_res0.8_dims30.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)


FeaturePlot_scCustom(KO_HPA_merged, features = "Pdgfra", split.by = "condition", label = T, raster = F, order = T)
saveRDS(KO_HPA_merged, "KO_HPA_merged_postUMAP.rds")

tibble(
  cluster = KO_HPA_merged$seurat_clusters,
  Dataset = KO_HPA_merged$condition
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
ggsave("KO_HPA_PieChart_res0.8_dims30.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)


