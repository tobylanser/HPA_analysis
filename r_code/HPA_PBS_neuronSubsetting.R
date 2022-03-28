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
HPA_PBS_merged <- readRDS("HPA_PBS_filtered.rds")
 

DimPlot_scCustom(HPA_PBS_merged, label = T, split.by = "condition")

FeaturePlot_scCustom(HPA_PBS_merged, features = "Erbb4", label = T, raster = F)
#Neurons
neuron_genes <- c("Rbfox1", "Rbfox3", "Asic2")

HPA_PBS_merged <- AddModuleScore(HPA_PBS_merged,
                               features = list(neuron_genes),
                               name="neurons_enriched")


# Plot scores
Neurons <- FeaturePlot(HPA_PBS_merged,
                       features = "neurons_enriched1", label = TRUE, repel = TRUE, raster = F, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
Neurons

HPA_PBS_neurons <- subset(HPA_PBS_merged, idents = c("1", "4", "7", "8", "6"))

HPA_PBS_neurons <- SCTransform(HPA_PBS_neurons, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
HPA_PBS_neurons <- RunPCA(object = HPA_PBS_neurons)
ElbowPlot(object = HPA_PBS_neurons, ndims = 50)
###
HPA_PBS_neurons <- FindNeighbors(object = HPA_PBS_neurons, dims = 1:40)
HPA_PBS_neurons <- FindClusters(object = HPA_PBS_neurons, resolution = 1.2)
HPA_PBS_neurons <- RunUMAP(HPA_PBS_neurons, dims = 1:40)
DimPlot_scCustom(WT_KO, label = T, split.by = "condition")

saveRDS(HPA_PBS_neurons, "Neurons_HPA_PBS.rds")

HPA_PBS_neurons <- readRDS("Neurons_HPA_PBS.rds")

test_marker <- FindMarkers(HPA_PBS_neurons, ident.1 = "12", only.pos = T)

WT_KO_neurons
FeaturePlot_scCustom(HPA_PBS, features = "Ttn", label = T, raster = F)


tibble(
  cluster = HPA_PBS_neurons$seurat_clusters,
  Dataset = HPA_PBS_neurons$condition
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
