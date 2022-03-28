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
HPA_PBS_merged <- readRDS("HPA_PBS_merged_unNormalized.rds")


#Cluster HPS & PBS conditions
HPA_PBS_merged <- SCTransform(HPA_PBS_merged, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
HPA_PBS_merged <- RunPCA(object = HPA_PBS_merged)
ElbowPlot(object = HPA_PBS_merged, ndims = 50)
###
HPA_PBS_merged <- FindNeighbors(object = HPA_PBS_merged, dims = 1:40)
HPA_PBS_merged <- FindClusters(object = HPA_PBS_merged, resolution = 1.2)
HPA_PBS_merged <- RunUMAP(HPA_PBS_merged, dims = 1:40)


DimPlot_scCustom(HPA_PBS_merged, label = T)

#ID junk clusters
junk <- subset(HPA_PBS_merged, idents = c("1", "0", "5", "4", "18", "23"))

junk <- SCTransform(junk, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
junk <- RunPCA(object = junk)
ElbowPlot(object = junk, ndims = 50)
###
junk <- FindNeighbors(object = junk, dims = 1:10)
junk <- FindClusters(object = junk, resolution = 1.2)
junk <- RunUMAP(junk, dims = 1:10)

DimPlot_scCustom(junk, label = T)

#filter out junk clusters
filtered <- subset(HPA_PBS_merged, idents = c("1", "0", "5", "4", "18", "23"), invert = T)

filtered <- SCTransform(filtered, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
filtered <- RunPCA(object = filtered)
ElbowPlot(object = filtered, ndims = 50)
###
filtered <- FindNeighbors(object = filtered, dims = 1:40)
filtered <- FindClusters(object = filtered, resolution = 1.2)
filtered <- RunUMAP(filtered, dims = 1:40)

DimPlot_scCustom(filtered, label = T)

FeaturePlot_scCustom(filtered, features = "Rbfox3", label = T, raster = F)

##More junk
more_junk <- subset(filtered, idents = c("1", "0", "4"), invert = T)

more_junk <- SCTransform(more_junk, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
more_junk <- RunPCA(object = more_junk)
ElbowPlot(object = more_junk, ndims = 50)
###
more_junk <- FindNeighbors(object = more_junk, dims = 1:15)
more_junk <- FindClusters(object = more_junk, resolution = 1.2)
more_junk <- RunUMAP(more_junk, dims = 1:15)
DimPlot_scCustom(HPA_PBS_merged, label = T, split.by = "condition")


FeaturePlot_scCustom(HPA_PBS_merged, features = "Phactr1", label = T, raster = F)

clust0 <- FindMarkers(HPA_PBS_merged, ident.1 = "0", only.pos = T)

saveRDS(filtered_2.0, "HPA_PBS_filtered.rds")

HPA_PBS_merged <- readRDS("HPA_PBS_filtered.rds")


#Funky glia
##More junk
filtered_2.0 <- subset(HPA_PBS_merged, idents = c("0", "12", "16"), invert = T)

filtered_2.0 <- SCTransform(filtered_2.0, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
filtered_2.0 <- RunPCA(object = filtered_2.0)
ElbowPlot(object = filtered_2.0, ndims = 50)
###
filtered_2.0 <- FindNeighbors(object = filtered_2.0, dims = 1:15)
filtered_2.0 <- FindClusters(object = filtered_2.0, resolution = 0.7)
filtered_2.0 <- RunUMAP(filtered_2.0, dims = 1:15)
DimPlot_scCustom(filtered_2.0, label = T, split.by = "condition")

funky_markers <- FindMarkers(filtered_2.0, ident.1 = c("6"), only.pos = T)

FeaturePlot_scCustom(filtered_2.0, features = "Zfpm2", label = T, raster = F)

tibble(
  cluster = filtered_2.0$seurat_clusters,
  Dataset = filtered_2.0$condition
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



