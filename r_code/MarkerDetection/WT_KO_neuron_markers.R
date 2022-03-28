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



setwd("~/Desktop/lanser_nucseq/hpa_rds_files")



#Load objects
WT_KO_neurons <- readRDS("WT_KO_neurons_postUMAP.rds")


##Marker Detection
#WT & KO

#Cluster markers
KO_WT_Neurons_clusterMarkers <- FindAllMarkers(neurons_KO_WT, only.pos = T)
KO_WT_Neurons_clusterMarkers_top50 <- KO_WT_Neurons_clusterMarkers %>%
  group_by(cluster) %>% 
  slice_max(n = 50, avg_log2FC)
write.csv(KO_WT_Neurons_clusterMarkers_top50, "/Users/WeinerLab/Desktop/lanser_nucseq/Markers/Top50/KO_WT_neurons_clusterMarkers_top50.csv")
write.csv(KO_WT_Neurons_clusterMarkers, "/Users/WeinerLab/Desktop/lanser_nucseq/Markers/All_markers/KO_WT_neurons_clusterMarkers.csv")
AvgExp_Ko_Wt <- AverageExpression(WT_KO_neurons)
write.csv(AvgExp_Ko_Wt$SCT, "/Users/WeinerLab/Desktop/lanser_nucseq/Markers/Average_Expression/AvgExp_Ko_Wt_neurons.csv")

#Condition Markers
Idents(WT_KO_neurons) <- "seurat_clusters"
KO_WT_Neurons_ConditionMarkers <- FindAllMarkers(WT_KO_neurons, only.pos = T)
KO_WT_Neurons_ConditionMarkers_top50 <- KO_WT_Neurons_ConditionMarkers %>%
  group_by(cluster) %>% 
  slice_max(n = 50, avg_log2FC)
write.csv(KO_WT_Neurons_ConditionMarkers_top50, "/Users/WeinerLab/Desktop/lanser_nucseq/Markers/Condition_Markers/KO_WT_Neurons_ConditionMarkers_top50.csv")
write.csv(KO_WT_Neurons_ConditionMarkers, "/Users/WeinerLab/Desktop/lanser_nucseq/Markers/Condition_Markers/KO_WT_Neurons_ConditionMarkers.csv")

top_markers <- Extract_Top_Markers(KO_WT_Neurons_ConditionMarkers)

DEG_conditions <- KO_WT_Neurons_ConditionMarkers$gene
DEG_conditions
Stacked_VlnPlot(WT_KO_neurons, features = DEG_conditions)
Plot_Density_Custom(WT_KO_neurons, features = DEG_conditions)
DimPlot_scCustom(WT_KO_neurons, label = T, raster = F)
ggsave("/Users/WeinerLab/Desktop/lanser_nucseq/Cluster_Plots/DimPlot/WT_KO_DimPlot_dims30res1.0.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)

Idents(WT_KO_merged) <- "seurat_clusters"
condition_diff <- subset(WT_KO_merged, idents = c("24"))
condition_diff <- SCTransform(condition_diff, method = "glmGamPoi", vars.to.regress = "percent.mt")
###PCA
condition_diff <- RunPCA(object = condition_diff)
ElbowPlot(object = condition_diff, ndims = 50)
###
condition_diff <- FindNeighbors(object = condition_diff, dims = 1:5)
condition_diff <- FindClusters(object = condition_diff, resolution = 0.2)
condition_diff <- RunUMAP(condition_diff, dims = 1:5)

DimPlot_scCustom(condition_diff, label = T, split.by = "condition")

clust24_markers <- FindMarkers(WT_KO_merged, ident.1 = "24", only.pos = T)

FeaturePlot_scCustom(WT_KO_merged, features = "Sox6", label = T)
FeaturePlot_scCustom(condition_diff, features = "Erbb4", label = T)

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
ggsave("/Users/WeinerLab/Desktop/lanser_nucseq/Cluster_Plots/PieCharts/WT_KO_piechart.jpeg", device = "jpeg", 
       width = 8, height = 9, units = "in", dpi = 600)


##Cell type annotation
#D1 MSNs
D1MSN_genes <- c("Tac1", "Drd1", "Isl1", "Sfxn1",
                 "Isl1")

WT_KO_neurons <- AddModuleScore(WT_KO_neurons,
                               features = list(D1MSN_genes),
                               name="D1MSN_enriched")

# Plot scores
D1_neurons <- FeaturePlot(WT_KO_neurons,
                          features = "D1MSN_enriched1", label = TRUE, repel = TRUE, raster = F, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
D1_neurons
#D2 MSNs
D2MSN_genes <- c("Drd2", "Adora2a", "Penk", "Gpr6",
                 "Gpr52", "Sp9")

WT_KO_neurons <- AddModuleScore(WT_KO_neurons,
                               features = list(D2MSN_genes),
                               name="D2MSN_enriched")

# Plot scores
D2_neurons <- FeaturePlot(WT_KO_neurons,
                          features = "D2MSN_enriched1", label = TRUE, repel = TRUE, raster = F, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
D2_neurons

#Astrocyte markers
Astro_genes <- c("Aldh1l1", "Slc1a3", "Aqp4",
                 "S100b")

WT_KO_neurons <- AddModuleScore(WT_KO_neurons,
                               features = list(Astro_genes),
                               name="Astro_enriched")

# Plot scores
Astrocyte <- FeaturePlot(WT_KO_neurons,
                         features = "Astro_enriched1", label = TRUE, repel = TRUE, raster = F, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
Astrocyte
#Oligodendrocyte markers
Oligo_genes <- c("Mog", "Mbp", "Olig2", "Sox10",
                 "Rnf220", "Plp1", "Pde1c")

WT_KO_neurons <- AddModuleScore(WT_KO_neurons,
                               features = list(Oligo_genes),
                               name="Oligo_enriched")

# Plot scores
Oligodendrocyte <- FeaturePlot(WT_KO_neurons,
                               features = "Oligo_enriched1", label = TRUE, repel = TRUE, raster = F, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
Oligodendrocyte
#Endothelial cells
Endothelial_genes <- c("Slc22a8", "Flt1", "Mylk", "Myl9",
                       "Mylk")

WT_KO_neurons <- AddModuleScore(WT_KO_neurons,
                               features = list(Endothelial_genes),
                               name="Endothelial_enriched")



# Plot scores
Endothelial <- FeaturePlot(WT_KO_neurons,
                           features = "Endothelial_enriched1", label = TRUE, repel = TRUE, raster = F, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
Endothelial
DimPlot_scCustom(WT_KO_neurons, label = T, split.by = "condition")

#Neurons
neuron_genes <- c("Trank1", "Dlgap1", "Erbb4", "Dpp10",
                  "Ptprd", "Pde4d", "Adarb2")

WT_KO_neurons <- AddModuleScore(WT_KO_neurons,
                               features = list(neuron_genes),
                               name="neurons_enriched")


# Plot scores
Neurons <- FeaturePlot(WT_KO_neurons,
                       features = "neurons_enriched1", label = TRUE, repel = TRUE, raster = F, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
Neurons



Neurons/Astrocyte|D2_neurons/Endothelial|D1_neurons/Oligodendrocyte
ggsave("/Users/WeinerLab/Desktop/lanser_nucseq/Cluster_Plots/CellType/KO_WT_CellType_Plot_res1.0_dims30.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)
