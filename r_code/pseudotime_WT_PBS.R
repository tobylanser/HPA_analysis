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
library(monocle3)
library(SeuratWrappers)



setwd("~/Desktop/lanser_nucseq/hpa_rds_files")

#Load objects
WT_PBS_neurons <- readRDS("WT_KO_neurons_postUMAP.rds")

DimPlot_scCustom(WT_PBS_neurons, label = T)

cds <- as.cell_data_set(WT_PBS_neurons)
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)


cds <- learn_graph(cds)
plot_cells(cds, label_groups_by_cluster = T, label_leaves = T, label_branch_points = T)


max.avp <- which.max(unlist(FetchData(WT_PBS_neurons, "Ttn")))
max.avp <- colnames(cds)[max.avp]
cds <- order_cells(cds, root_cells = max.avp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = T, label_leaves = T, 
           label_branch_points = T)
