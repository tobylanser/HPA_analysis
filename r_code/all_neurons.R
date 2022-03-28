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
library(sctransform)

setwd("~/Desktop/lanser_nucseq/hpa_rds_files")

#Load objects
all <- readRDS("allCells_filtered.rds")

#Neurons
neuron_genes <- c("Rbfox1", "Rbfox3", "Asic2")

all <- AddModuleScore(all, features = list(neuron_genes),
                      name="neurons_enriched")
clust_24 <- FindMarkers(all, ident.1 = "24", only.pos = T)

# Plot scores
Neurons <- FeaturePlot(all,
                       features = "neurons_enriched1", label = TRUE, repel = TRUE, raster = F, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
Neurons

all_neurons <- subset(all, idents = c("2", "3", "4", "6",
                                      "7", "8", "9", "17",
                                      "22", "24"))

all_neurons <- SCTransform(all_neurons, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
all_neurons <- RunPCA(object = all_neurons)
ElbowPlot(object = all_neurons, ndims = 50)
###
all_neurons <- FindNeighbors(object = all_neurons, dims = 1:40)
all_neurons <- FindClusters(object = all_neurons, resolution = 1.2)
all_neurons <- RunUMAP(all_neurons, dims = 1:40)
DimPlot_scCustom(all_neurons, label = T, split.by = "condition")

FeaturePlot_scCustom(all_neurons, features = "Adarb2", label = T)

tibble(
  cluster = all_neuron2.0$seurat_clusters,
  Dataset = all_neuron2.0$condition
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



saveRDS(all_neurons, "all_neurons.rds")

