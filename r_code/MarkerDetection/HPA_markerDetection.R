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
HPA_PBS_merged <- readRDS("HPA_PBS_merged_postUMAP.rds")
WT_KO_merged <- readRDS("WT_KO_merged_postUMAP.rds")
KO_HPA_merged <- readRDS("KO_HPA_merged_postUMAP.rds")

##Marker Detection

#HPA & PBS
HPA_PBS_clusterMarkers <- FindAllMarkers(HPA_PBS_merged, only.pos = TRUE)
HPA_PBS_clusterMarkers_top50 <- HPA_PBS_clusterMarkers %>%
  group_by(cluster) %>% 
  slice_max(n = 50, avg_log2FC)
write.csv(HPA_PBS_clusterMarkers_top50, "HPA_PBS_clusterMarkers_top50.csv")
write.csv(HPA_PBS_clusterMarkers, "HPA_PBS_clusterMarkers.csv")
AvgExp_HPS_PBS <- AverageExpression(HPA_PBS_merged)
write.csv(AvgExp_HPS_PBS, "AvgExp_HPS_PBS.csv")

DimPlot_scCustom(HPA_PBS_merged, split.by = "condition", label = T, raster = F)
ggsave("HPA_PBS_DimPlot.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)


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
ggsave("HPA_PBS_piechart.jpeg", device = "jpeg", 
       width = 8, height = 9, units = "in", dpi = 600)


#KO & HPA
KO_HPA_merged_markers <- FindAllMarkers(KO_HPA_merged, only.pos = TRUE)
KO_HPA_merged_markers_top50 <- KO_HPA_merged_markers %>%
  group_by(cluster) %>% 
  slice_max(n = 50, avg_log2FC)
write.csv(KO_HPA_merged_markers_top50, "KO_HPA_merged_markers_top50.csv")
write.csv(KO_HPA_merged_markers, "KO_HPA_merged_markers.csv")
AvgExp_KO_HPA <- AverageExpression(KO_HPA_merged)
write.csv(AvgExp_KO_HPA, "AvgExp_KO_HPA.csv")

DimPlot_scCustom(KO_HPA_merged, split.by = "condition", label = T, raster = F)
ggsave("KO_HPA_DimPlot.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)

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
ggsave("KO_HPA_piechart.jpeg", device = "jpeg", 
       width = 8, height = 9, units = "in", dpi = 600)


#WT & KO
KO_WT_clusterMarkers <- FindAllMarkers(WT_KO_merged, only.pos = T)
KO_WT_clusterMarkers_top50 <- KO_WT_clusterMarkers %>%
  group_by(cluster) %>% 
  slice_max(n = 50, avg_log2FC)
write.csv(KO_WT_clusterMarkers_top50, "KO_WT_clusterMarkers_top50.csv")
write.csv(KO_WT_clusterMarkers, "KO_WT_clusterMarkers.csv")
AvgExp_Ko_Wt <- AverageExpression(WT_KO_merged)
write.csv(AvgExp_Ko_Wt, "AvgExp_Ko_Wt.csv")

DimPlot_scCustom(WT_KO_merged, split.by = "condition", label = T, raster = F)
ggsave("WT_KO_DimPlot.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)

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
ggsave("WT_KO_piechart.jpeg", device = "jpeg", 
       width = 8, height = 9, units = "in", dpi = 600)


FeaturePlot(HPA_PBS_merged, features = "mt-Co1", split.by = "condition", label = T, raster=FALSE, order = T)
Plot_Density_Custom(seurat_object = HPA_PBS_merged, features = "mt-Co1", label = T, custom_palette = PurpleAndYellow())

DimPlot_scCustom(HPA_PBS_merged, split.by = "condition", label = T, raster = F)

FeaturePlot(HPA_PBS_merged, features = D2MSN_genes)


##Cell type annotation

#D1 MSNs
D1MSN_genes <- c("Tac1", "Drd1", "Isl1", "Sfxn1",
                 "Isl1")

HPA_PBS_merged <- AddModuleScore(HPA_PBS_merged,
                                 features = list(D1MSN_genes),
                                 name="D1MSN_enriched")

# Plot scores
D1_neurons <- FeaturePlot(HPA_PBS_merged,
                          features = "D1MSN_enriched1", label = TRUE, repel = TRUE, raster = F, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
D1_neurons


#D2 MSNs
D2MSN_genes <- c("Drd2", "Adora2a", "Penk", "Gpr6",
                 "Gpr52", "Sp9")

HPA_PBS_merged <- AddModuleScore(HPA_PBS_merged,
                                 features = list(D2MSN_genes),
                                 name="D2MSN_enriched")

# Plot scores
D2_neurons <- FeaturePlot(HPA_PBS_merged,
                          features = "D2MSN_enriched1", label = TRUE, repel = TRUE, raster = F, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

D1_neurons|D2_neurons

#Astrocyte markers
Astro_genes <- c("Aldh1l1", "Slc1a3", "Aqp4", "AldoC", "Glt1",
                 "S100b")

HPA_PBS_merged <- AddModuleScore(HPA_PBS_merged,
                                 features = list(Astro_genes),
                                 name="Astro_enriched")

# Plot scores
Astrocyte <- FeaturePlot(HPA_PBS_merged,
                         features = "Astro_enriched1", label = TRUE, repel = TRUE, raster = F, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
Astrocyte
#Oligodendrocyte markers
Oligo_genes <- c("Mog", "Mbp", "Olig2", "Sox10")

HPA_PBS_merged <- AddModuleScore(HPA_PBS_merged,
                                 features = list(Oligo_genes),
                                 name="Oligo_enriched")

# Plot scores
Oligodendrocyte <- FeaturePlot(HPA_PBS_merged,
                               features = "Oligo_enriched1", label = TRUE, repel = TRUE, raster = F, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
Oligodendrocyte
#Endothelial cells
Endothelial_genes <- c("Slc22a8", "Flt1", "Mylk", "Myl9",
                       "Mylk")

HPA_PBS_merged <- AddModuleScore(HPA_PBS_merged,
                                 features = list(Endothelial_genes),
                                 name="Endothelial_enriched")



# Plot scores
Endothelial <- FeaturePlot(HPA_PBS_merged,
                           features = "Endothelial_enriched1", label = TRUE, repel = TRUE, raster = F, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
Endothelial
DimPlot_scCustom(HPA_PBS_merged, label = T, split.by = "condition")

#Neurons
neuron_genes <- c("Meis2", "Myt1l", "Mylk", "Bcl11b",
                  "Rarb", "Rxrg", "Foxp2", "Isl1", "Camta1", 
                  "Dpf1", "Bcl11a", "Tubb3", "Peg3", "Hivep2", 
                  "Trerf1", "Klf16", "Zfhx2", "Npas2", "Zfp941",
                  "Pou3f1", "Sp9", "Ebf1", "Gatad1")

HPA_PBS_merged <- AddModuleScore(HPA_PBS_merged,
                                 features = list(neuron_genes),
                                 name="neurons_enriched")


# Plot scores
Neurons <- FeaturePlot(HPA_PBS_merged,
                       features = "neurons_enriched1", label = TRUE, repel = TRUE, raster = F, order = T) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
Neurons



Neurons/Astrocyte|D2_neurons/Endothelial|D1_neurons/Oligodendrocyte
ggsave("PBS_HPA_CellType_Plot.jpeg", device = "jpeg", 
       width = 16, height = 9, units = "in", dpi = 600)

