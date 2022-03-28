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

FeaturePlot_scCustom(WT_KO_merged, features = c("Rbfox1", "Rbfox3", "Asic2"), label = T) + ncol(3)
#Neurons
neuron_genes <- c("Rbfox1", "Rbfox3", "Asic2")

WT_KO_merged <- AddModuleScore(WT_KO_merged,
                                features = list(neuron_genes),
                                name="neurons_enriched")


# Plot scores
Neurons <- FeaturePlot(WT_KO_merged,
                       features = "neurons_enriched1", label = TRUE, repel = TRUE, raster = F, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
Neurons

DimPlot_scCustom(WT_KO_merged, label = T)

FeaturePlot_scCustom(WT_KO_neurons, features = c("Ptprt"), label = T)
Plot_Density_Custom(WT_KO_merged, features = c("Trank1"))
WT_KO_merged_unfiltered
Plot_Density_Joint_Only(WT_KO_merged, features =  c("Trank1", "Dlgap1", "Erbb4", "Dpp10",
                                                "Ptprd", "Pde4d", "Adarb2"))

WT_KO_neurons

DimPlot_scCustom(WT_KO_neurons, label = T, split.by = "condition")
WT_KO_neurons <- readRDS("WT_KO_neurons.rds")
Idents(WT_KO_merged) <- "cell_type"
WT_KO_neurons <- subset(WT_KO_merged, idents = c("Excitatory Neurons", "Inhibitory Neurons"))

#Cluster WT & KO Neurons
WT_KO_neurons <- SCTransform(WT_KO_neurons, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
WT_KO_neurons <- RunPCA(object = WT_KO_neurons)
ElbowPlot(object = WT_KO_neurons, ndims = 50)
###
WT_KO_neurons <- FindNeighbors(object = WT_KO_neurons, dims = 1:10)
WT_KO_neurons <- FindClusters(object = WT_KO_neurons, resolution = 0.8)
WT_KO_neurons <- RunUMAP(WT_KO_neurons, dims = 1:10)

DimPlot_scCustom(WT_KO_neurons, label = T, group.by = "cell_type")


ggsave("WT_PBS_Condition_Plot_res0.7_dims40.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)



saveRDS(neurons_KO_WT, "WT_KO_neurons_postUMAP.rds")

WT_KO_neurons <- readRDS("WT_KO_neurons_postUMAP.rds")



tibble(
  cluster = WT_KO_neurons$seurat_clusters,
  Dataset = WT_KO_neurons$condition
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


neurons_KO_WT <- subset(WT_KO_merged, idents = c("4", "5", "6", "7", "9"))

neurons_KO_WT <- SCTransform(neurons_KO_WT, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
neurons_KO_WT <- RunPCA(object = neurons_KO_WT)
ElbowPlot(object = neurons_KO_WT, ndims = 50)
###  
neurons_KO_WT <- FindNeighbors(object = neurons_KO_WT, dims = 1:40)
neurons_KO_WT <- FindClusters(object = neurons_KO_WT, resolution = 0.5)
neurons_KO_WT <- RunUMAP(neurons_KO_WT, dims = 1:40)

DimPlot_scCustom(neurons_KO_WT, label = T, split.by = "condition")

FeaturePlot_scCustom(neurons_KO_WT, features = "Galntl6", label = T)

Neurons

clust3_DEG <- c("Gm42418", "Sox6", "Cntnap5a", "Adarb2",
                "Gpc5", "Gm26917", "Cmss1", "Dst")

inhib_markers <- c("Erbb4", "Adarb2",
                   "Robo1", "Galntl6", "Inpp4b")

excit_markers <- c("Pde4d", "Dpp10", "Ptprd",
                   "Cdh18", "Kcnip4")

merged_markers <- c("Pde4d", "Dpp10", "Kcnip4", "Ptprd",
                 "Cdh18","Erbb4",
                 "Adarb2", "Robo1", "Galntl6", "Inpp4b")

DotPlot_scCustom(neurons_KO_WT, features = merged_markers)
ggsave("/Users/WeinerLab/Desktop/lanser_nucseq/WT_KO_neuron_dotplot_markers.jpeg", device = "jpeg", 
       width = 16, height = 8, units = "in", dpi = 600)

FeaturePlot_scCustom(WT_KO_neurons, features = "Gphn", label = T)
ggsave("WT_KO_neuron_inhib_markers.jpeg", device = "jpeg", 
       width = 16, height = 8, units = "in", dpi = 600)

FeaturePlot_scCustom(neurons_KO_WT, features = "Kcnip4", label = T, label.size = 8)

DotPlot_scCustom(neurons_KO_WT, features = merged_markers)
ggsave("WT_KO_neuron_excit_markers.jpeg", device = "jpeg", 
       width = 16, height = 8, units = "in", dpi = 600)


DimPlot_scCustom(neurons_KO_WT, split.by = "condition", label = T, label.size = 8, fontface = "bold")
LabelClusters(p, id = "orig.ident")
p + legend(text.font=4)
DimPlot_scCustom(neurons_KO_WT, label = F, split.by = "condition")
LabelClusters(p, id = "ident", size = 8, fontface = "bold", color = "black")


ggsave("WT_KO_neuron_split.by_condition_noLabels.jpeg", device = "jpeg", 
       width = 16, height = 8, units = "in", dpi = 600)

Cluster_Highlight_Plot(seurat_object = WT_KO_neurons, cluster_name = c("3", "11"), split.by = "condition", highlight_color = "dodgerblue3",
                       background_color = "lightgray") 
ggsave("WT_KO_neuron_condition_clusterhighlight.jpeg", device = "jpeg", 
       width = 16, height = 8, units = "in", dpi = 600)


Meta_Highlight_Plot(seurat_object = WT_KO_neurons, meta_data_column = "cell_type",
                    meta_data_highlight = "Inhibitory Neurons", highlight_color = "firebrick", background_color = "lightgray", label = T)

Plot_Density_Joint_Only(WT_KO_neurons, features = inhib_markers)

FeaturePlot_scCustom(WT_KO_neurons, features = "Celf2")
ggsave("WT_KO_neuron_DEG_dotplot.jpeg", device = "jpeg", 
       width = 16, height = 8, units = "in", dpi = 600)
tibble(
  cluster = neurons_KO_WT$seurat_clusters,
  Dataset = neurons_KO_WT$condition
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
  ylab(NULL) +
  ggtitle("Percentage of KO vs WT in each neuronal type")
ggsave("/Users/WeinerLab/Desktop/lanser_nucseq/WT_KO_neurons_PieChart.jpeg", device = "jpeg", 
       width = 8, height = 9, units = "in", dpi = 600)


tibble(
  cluster = WT_KO_neurons$cell_type,
  Dataset = WT_KO_neurons$condition
) %>%
  group_by(cluster, Dataset) %>%
  count() %>%
  group_by(cluster) %>%
  mutate(
    percent=(100*n)/sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    cluster=paste("Cell Type:", cluster)
  ) %>%
  ggplot(aes(x="",y=percent, fill=Dataset)) +
  geom_col(width=1) +
  coord_polar("y", start=0) +
  facet_wrap(vars(cluster)) +
  theme(axis.text.x=element_blank()) +
  xlab(NULL) +
  ylab(NULL)

WT_KO_merged@meta.data %>%
  group_by(cell_type, condition) %>%
  count() %>%
  group_by(cell_type) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=cell_type,y=percent, fill=condition)) +
  geom_col() +
  ggtitle("Percentage of KO vs WT in each neuronal type")


inhibitory_clust1 <- FindMarkers(WT_KO_neurons, ident.1 = c("1"), only.pos = T)
Kcnip4
Idents(WT_KO_neurons) <- "seurat_clusters"
FeaturePlot_scCustom(WT_KO_neurons, features = "Chat", label = T)

Idents(WT_KO_neurons) <- "seurat_clusters"
Inhib <- subset(WT_KO_neurons, idents = c("2", "4", "5", "6",
                                         "8", "9", "12"))
Inhib$cell_type <- "Inhibitory Neurons"

Excit <- subset(WT_KO_neurons, idents = c("0", "1", "3",
                                          "7", "10", "11"))
Excit$cell_type <- "Excitatory Neurons"

WT_KO_neurons <- merge(Excit, Inhib)

WT_KO_neurons <- SCTransform(WT_KO_neurons, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
WT_KO_neurons <- RunPCA(object = WT_KO_neurons)
ElbowPlot(object = WT_KO_neurons, ndims = 50)
###  
WT_KO_neurons <- FindNeighbors(object = WT_KO_neurons, dims = 1:30)
WT_KO_neurons <- FindClusters(object = WT_KO_neurons, resolution = 0.7)
WT_KO_neurons <- RunUMAP(WT_KO_neurons, dims = 1:30)

DimPlot_scCustom(WT_KO_neurons, label = T, split.by = "condition")

saveRDS(WT_KO_neurons, "WT_KO_neurons_postUMAP.rds")

WT_KO_neurons <- readRDS("WT_KO_neurons_postUMAP.rds")


