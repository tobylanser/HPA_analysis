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
library(scCustomize)
library(patchwork)
library(viridis)
library(Nebulosa)
library(sctransform)


setwd("~/Desktop/lanser_nucseq/hpa_rds_files")

#Load objects
WT_KO_merged <- readRDS("WT_KO_merged_postUMAP.rds")


#Load objects




DimPlot_scCustom(WT_KO_merged, label = T)


FeaturePlot_scCustom(WT_KO_merged, features = c("Trank1", "Erbb4",
                                           "Ptprd", "Pde4d", "Slc6a1"), label = T)




FeaturePlot_scCustom(WT_KO_merged, features = "Rnf220", label = T)

Idents(WT_KO_merged) <- "seurat_clusters"
neurons <- subset(WT_KO_merged, idents = c("2", "3", "14", "18",
                                           "19", "21", "33"))

#Cluster HPS & PBS conditions
neurons <- SCTransform(neurons, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
neurons <- RunPCA(object = neurons)
ElbowPlot(object = neurons, ndims = 50)
###
neurons <- FindNeighbors(object = neurons, dims = 1:20)
neurons <- FindClusters(object = neurons, resolution = 1.2)
neurons <- RunUMAP(neurons, dims = 1:20)
DimPlot_scCustom(neurons, label = T, group.by = "orig.ident")



FeaturePlot_scCustom(neurons, features = "Xlr5a", label = T, split.by = "condition")

tibble(
  cluster = neurons$seurat_clusters,
  Dataset = neurons$condition
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


