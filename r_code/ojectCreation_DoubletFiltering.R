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
library(DoubletFinder)
library(parallel)
library(glmGamPoi)



setwd("~/Desktop/lanser_nucseq/hpa_rds_files")

##GEM1
gem1_data_dir <- '~/Desktop/lanser_nucseq/all_filtered_feature_bc_matrix/gem1_filtered_feature_bc_matrix'
list.files(gem1_data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
gem1_data_files <- Read10X(data.dir = gem1_data_dir)
gem1_seurat_object = CreateSeuratObject(gem1_data_files, min.cells = 3, min.genes = 350)
gem1_seurat_object@meta.data$condition <- "HPA"
gem1_seurat_object@meta.data$Sex <- "M"
gem1_seurat_object@meta.data$orig.ident <- "GEM1"
gem1_seurat_object <- gem1_seurat_object[!grepl("Malat1", row.names(gem1_seurat_object)), ]
gem1_seurat_object[["percent.mt"]] <- PercentageFeatureSet(gem1_seurat_object, pattern = "^mt-")
plot1 <- FeatureScatter(gem1_seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gem1_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- VlnPlot(object = gem1_seurat_object, features = c("nFeature_RNA"), pt.size = 0.1) +NoLegend()
plot1 + plot2 + plot3
gem1_seurat_object <- subset(gem1_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 700 & percent.mt < 7)
gem1_seurat_object <- SCTransform(gem1_seurat_object, method = "glmGamPoi", vars.to.regress = "percent.mt")
gem1_seurat_object <- RunPCA(object = gem1_seurat_object, verbose = F)
ElbowPlot(object = gem1_seurat_object, ndims = 30)
gem1_seurat_object <- FindNeighbors(object = gem1_seurat_object, dims = 1:20)
gem1_seurat_object <- FindClusters(object = gem1_seurat_object, resolution = 1.0)
##
gem1_seurat_object <- RunUMAP(gem1_seurat_object, dims = 1:20)

DimPlot_scCustom(gem1_seurat_object, label = T, split.by = "condition")

sweep.res.list_gem1 <- paramSweep_v3(gem1_seurat_object, PCs = 1:20, sct = T, num.cores = 2)
sweep.stats_gem1 <- summarizeSweep(sweep.res.list_gem1, GT = FALSE)
bcmv_gem1 <- find.pK(sweep.stats_gem1)
homotypic.prop <- modelHomotypic(gem1_seurat_object$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.035*nrow(gem1_seurat_object@meta.data)) #tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
gem1_seurat_object <- doubletFinder_v3(gem1_seurat_object, PCs = 1:20, pN = 0.25, pK = 0.04, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
gem1_seurat_object@meta.data[,"CellTypes_DF"] <- gem1_seurat_object$DF.classifications_0.25_0.04_293
gem1_seurat_object@meta.data$CellTypes_DF[which(gem1_seurat_object$pANN_0.25_0.04_293 == "Doublet")] <- "Doublet"
DimPlot(gem1_seurat_object, group.by="CellTypes_DF", reduction="umap", pt.size=0.5, order=c("Coll.Duct.TC","Doublet"), cols=c("#66C2A5","#FFD92F","#8DA0CB","#A6D854","#E78AC3","#B3B3B3","#E5C494","black","#FC8D62"))
saveRDS(gem1_seurat_object, "gem1_seurat_withDoublets.rds")

##GEM2
gem2_data_dir <- '~/Desktop/lanser_nucseq/all_filtered_feature_bc_matrix/gem2_filtered_feature_bc_matrix'
list.files(gem2_data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
gem2_data <- Read10X(data.dir = gem2_data_dir)
gem2_seurat_object = CreateSeuratObject(gem2_data, min.cells = 3, min.genes = 350)
gem2_seurat_object@meta.data$condition <- "HPA"
gem2_seurat_object@meta.data$Sex <- "M"
gem2_seurat_object@meta.data$orig.ident <- "GEM2"
gem2_seurat_object <- gem2_seurat_object[!grepl("Malat1", row.names(gem2_seurat_object)), ]
gem2_seurat_object[["percent.mt"]] <- PercentageFeatureSet(gem2_seurat_object, pattern = "^mt-")
plot1 <- FeatureScatter(gem2_seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gem2_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- VlnPlot(object = gem2_seurat_object, features = c("nFeature_RNA"), pt.size = 0.1) +NoLegend()
plot1 + plot2 + plot3
gem2_seurat_object <- subset(gem2_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 500 & percent.mt < 10)
gem2_seurat_object <- SCTransform(gem2_seurat_object, method = "glmGamPoi", vars.to.regress = "percent.mt")
gem2_seurat_object <- RunPCA(object = gem2_seurat_object, verbose = F)
ElbowPlot(object = gem2_seurat_object, ndims = 50)
gem2_seurat_object <- FindNeighbors(object = gem2_seurat_object, dims = 1:20)
gem2_seurat_object <- FindClusters(object = gem2_seurat_object, resolution = 1.0)
##
gem2_seurat_object <- RunUMAP(gem2_seurat_object, dims = 1:20)
DimPlot_scCustom(gem2_seurat_object, label = T, split.by = "condition")

sweep.res.list_gem2 <- paramSweep_v3(gem2_seurat_object, PCs = 1:20, sct = T, num.cores = 2)
sweep.stats_gem2 <- summarizeSweep(sweep.res.list_gem2, GT = FALSE)
bcmv_gem2 <- find.pK(sweep.stats_gem2)
homotypic.prop <- modelHomotypic(gem2_seurat_object$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.035*nrow(gem2_seurat_object@meta.data)) #tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
gem2_seurat_object <- doubletFinder_v3(gem2_seurat_object, PCs = 1:20, pN = 0.25, pK = 0.02, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
gem2_seurat_object@meta.data[,"CellTypes_DF"] <- gem2_seurat_object$DF.classifications_0.25_0.02_290
gem2_seurat_object@meta.data$CellTypes_DF[which(gem2_seurat_object$pANN_0.25_0.02_290 == "Doublet")] <- "Doublet"
DimPlot(gem2_seurat_object, group.by="CellTypes_DF", reduction="umap", pt.size=0.5, order=c("Coll.Duct.TC","Doublet"), cols=c("#66C2A5","#FFD92F","#8DA0CB","#A6D854","#E78AC3","#B3B3B3","#E5C494","black","#FC8D62"))
saveRDS(gem2_seurat_object, "gem2_seurat_withDoublets.rds")

##GEM3
gem3_data_dir <- '~/Desktop/lanser_nucseq/all_filtered_feature_bc_matrix/gem3_filtered_feature_bc_matrix'
list.files(gem3_data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
gem3_data <- Read10X(data.dir = gem3_data_dir)
gem3_seurat_object = CreateSeuratObject(gem3_data, min.cells = 3, min.genes = 350)
gem3_seurat_object@meta.data$condition <- "PBS"
gem3_seurat_object@meta.data$Sex <- "M"
gem3_seurat_object@meta.data$orig.ident <- "GEM3"
gem3_seurat_object <- gem3_seurat_object[!grepl("Malat1", row.names(gem3_seurat_object)), ]
gem3_seurat_object[["percent.mt"]] <- PercentageFeatureSet(gem3_seurat_object, pattern = "^mt-")
plot1 <- FeatureScatter(gem3_seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gem3_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- VlnPlot(object = gem3_seurat_object, features = c("nFeature_RNA"), pt.size = 0.1) +NoLegend()
plot1 + plot2 + plot3
gem3_seurat_object <- subset(gem3_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 750 & percent.mt < 5)
gem3_seurat_object <- SCTransform(gem3_seurat_object, method = "glmGamPoi", vars.to.regress = "percent.mt")
gem3_seurat_object <- RunPCA(object = gem3_seurat_object, verbose = F)
ElbowPlot(object = gem3_seurat_object, ndims = 20)
gem3_seurat_object <- FindNeighbors(object = gem3_seurat_object, dims = 1:10)
gem3_seurat_object <- FindClusters(object = gem3_seurat_object, resolution = 1.2)
gem3_seurat_object <- RunUMAP(gem3_seurat_object, dims = 1:10)
DimPlot_scCustom(gem3_seurat_object, label = T, split.by = "condition")

sweep.res.list_gem3 <- paramSweep_v3(gem3_seurat_object, PCs = 1:10, sct = T, num.cores = 2)
sweep.stats_gem3 <- summarizeSweep(sweep.res.list_gem3, GT = FALSE)
bcmv_gem3 <- find.pK(sweep.stats_gem3)
homotypic.prop <- modelHomotypic(gem3_seurat_object$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.069*nrow(gem3_seurat_object@meta.data)) #tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
gem3_seurat_object <- doubletFinder_v3(gem3_seurat_object, PCs = 1:10, pN = 0.25, pK = 0.13, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
gem3_seurat_object@meta.data[,"CellTypes_DF"] <- gem3_seurat_object$DF.classifications_0.25_0.13_1099
gem3_seurat_object@meta.data$CellTypes_DF[which(gem3_seurat_object$pANN_0.25_0.13_1099 == "Doublet")] <- "Doublet"
DimPlot(gem3_seurat_object, group.by="CellTypes_DF", reduction="umap", pt.size=0.5, order=c("Coll.Duct.TC","Doublet"), cols=c("#66C2A5","#FFD92F","#8DA0CB","#A6D854","#E78AC3","#B3B3B3","#E5C494","black","#FC8D62"))
saveRDS(gem3_seurat_object, "gem3_seurat_withDoublets.rds")


##GEM4
gem4_data_dir <- '~/Desktop/lanser_nucseq/all_filtered_feature_bc_matrix/gem4_filtered_feature_bc_matrix'
list.files(gem4_data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
gem4_data <- Read10X(data.dir = gem4_data_dir)
gem4_seurat_object = CreateSeuratObject(gem4_data, min.cells = 3, min.genes = 350)
gem4_seurat_object@meta.data$condition <- "PBS"
gem4_seurat_object@meta.data$Sex <- "F"
gem4_seurat_object@meta.data$orig.ident <- "GEM4"
gem4_seurat_object <- gem4_seurat_object[!grepl("Malat1", row.names(gem4_seurat_object)), ]
gem4_seurat_object[["percent.mt"]] <- PercentageFeatureSet(gem4_seurat_object, pattern = "^mt-")
plot1 <- FeatureScatter(gem4_seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gem4_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- VlnPlot(object = gem4_seurat_object, features = c("nFeature_RNA"), pt.size = 0.1) +NoLegend()
plot1 + plot2 + plot3
gem4_seurat_object <- subset(gem4_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 1000 & percent.mt < 10)
gem4_seurat_object <- SCTransform(gem4_seurat_object, method = "glmGamPoi", vars.to.regress = "percent.mt")
gem4_seurat_object <- RunPCA(object = gem4_seurat_object, verbose = F)
ElbowPlot(object = gem4_seurat_object, ndims = 20)
gem4_seurat_object <- FindNeighbors(object = gem4_seurat_object, dims = 1:15)
gem4_seurat_object <- FindClusters(object = gem4_seurat_object, resolution = 1.2)
gem4_seurat_object <- RunUMAP(gem4_seurat_object, dims = 1:15)
DimPlot_scCustom(gem4_seurat_object, label = T, split.by = "condition")

sweep.res.list_gem4 <- paramSweep_v3(gem4_seurat_object, PCs = 1:15, sct = T, num.cores = 2)
sweep.stats_gem4 <- summarizeSweep(sweep.res.list_gem4, GT = FALSE)
bcmv_gem4 <- find.pK(sweep.stats_gem4)
homotypic.prop <- modelHomotypic(gem4_seurat_object$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.09*nrow(gem4_seurat_object@meta.data)) #tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
gem4_seurat_object <- doubletFinder_v3(gem4_seurat_object, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
gem4_seurat_object@meta.data[,"CellTypes_DF"] <- gem4_seurat_object$DF.classifications_0.25_0.09_1940
gem4_seurat_object@meta.data$CellTypes_DF[which(gem4_seurat_object$pANN_0.25_0.09_1940 == "Doublet")] <- "Doublet"
DimPlot(gem4_seurat_object, group.by="CellTypes_DF", reduction="umap", pt.size=0.5, order=c("Coll.Duct.TC","Doublet"), cols=c("#66C2A5","#FFD92F","#8DA0CB","#A6D854","#E78AC3","#B3B3B3","#E5C494","black","#FC8D62"))
saveRDS(gem4_seurat_object, "gem4_seurat_withDoublets.rds")



tibble(
  cluster = gem4_seurat_object$seurat_clusters,
  Dataset = gem4_seurat_object$CellTypes_DF
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



##GEM5
gem5_data_dir <- '~/Desktop/lanser_nucseq/all_filtered_feature_bc_matrix/gem5_filtered_feature_bc_matrix'
list.files(gem5_data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
gem5_data <- Read10X(data.dir = gem5_data_dir)
gem5_seurat_object = CreateSeuratObject(gem5_data, min.cells = 3, min.genes = 350)
gem5_seurat_object@meta.data$condition <- "WT"
gem5_seurat_object@meta.data$Sex <- "M"
gem5_seurat_object@meta.data$orig.ident <- "GEM5"
gem5_seurat_object <- gem5_seurat_object[!grepl("Malat1", row.names(gem5_seurat_object)), ]
gem5_seurat_object[["percent.mt"]] <- PercentageFeatureSet(gem5_seurat_object, pattern = "^mt-")
plot1 <- FeatureScatter(gem5_seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gem5_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- VlnPlot(object = gem5_seurat_object, features = c("nFeature_RNA"), pt.size = 0.1) +NoLegend()
plot1 + plot2 + plot3
gem5_seurat_object <- subset(gem5_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 750 & percent.mt < 7)
gem5_seurat_object <- SCTransform(gem5_seurat_object, method = "glmGamPoi", vars.to.regress = "percent.mt")
gem5_seurat_object <- RunPCA(object = gem5_seurat_object, verbose = F)
ElbowPlot(object = gem5_seurat_object, ndims = 20)
gem5_seurat_object <- FindNeighbors(object = gem5_seurat_object, dims = 1:10)
gem5_seurat_object <- FindClusters(object = gem5_seurat_object, resolution = 1.5)
gem5_seurat_object <- RunUMAP(gem5_seurat_object, dims = 1:10)
DimPlot_scCustom(gem5_seurat_object, label = T, split.by = "condition")

sweep.res.list_gem5 <- paramSweep_v3(gem5_seurat_object, PCs = 1:10, sct = T, num.cores = 2)
sweep.stats_gem5 <- summarizeSweep(sweep.res.list_gem5, GT = FALSE)
bcmv_gem5 <- find.pK(sweep.stats_gem5)
homotypic.prop <- modelHomotypic(gem5_seurat_object$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.15*nrow(gem5_seurat_object@meta.data)) #tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
gem5_seurat_object <- doubletFinder_v3(gem5_seurat_object, PCs = 1:10, pN = 0.25, pK = 0.16, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
gem5_seurat_object@meta.data[,"CellTypes_DF"] <- gem5_seurat_object$DF.classifications_0.25_0.16_5429
gem5_seurat_object@meta.data$CellTypes_DF[which(gem5_seurat_object$pANN_0.25_0.16_5429 == "Doublet")] <- "Doublet"
DimPlot(gem5_seurat_object, group.by="CellTypes_DF", reduction="umap", pt.size=0.5, order=c("Coll.Duct.TC","Doublet"), cols=c("#66C2A5","#FFD92F","#8DA0CB","#A6D854","#E78AC3","#B3B3B3","#E5C494","black","#FC8D62"))
saveRDS(gem5_seurat_object, "gem5_seurat_withDoublets.rds")

tibble(
  cluster = gem5_seurat_object$seurat_clusters,
  Dataset = gem5_seurat_object$DF.classifications_0.25_0.16_5429
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

##GEM6
gem6_data_dir <- '~/Desktop/lanser_nucseq/all_filtered_feature_bc_matrix/gem6_filtered_feature_bc_matrix'
list.files(gem6_data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
gem6_data <- Read10X(data.dir = gem6_data_dir)
gem6_seurat_object = CreateSeuratObject(gem6_data, min.cells = 3, min.genes = 350)
gem6_seurat_object@meta.data$condition <- "WT"
gem6_seurat_object@meta.data$Sex <- "F"
gem6_seurat_object@meta.data$orig.ident <- "GEM6"
gem6_seurat_object <- gem6_seurat_object[!grepl("Malat1", row.names(gem6_seurat_object)), ]
gem6_seurat_object[["percent.mt"]] <- PercentageFeatureSet(gem6_seurat_object, pattern = "^mt-")
plot1 <- FeatureScatter(gem6_seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gem6_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- VlnPlot(object = gem6_seurat_object, features = c("nFeature_RNA"), pt.size = 0.1) +NoLegend()
plot1 + plot2 + plot3
gem6_seurat_object <- subset(gem6_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 1000 & percent.mt < 5)
gem6_seurat_object <- SCTransform(gem6_seurat_object, method = "glmGamPoi", vars.to.regress = "percent.mt")
gem6_seurat_object <- RunPCA(object = gem6_seurat_object, verbose = F)
ElbowPlot(object = gem6_seurat_object, ndims = 50)
gem6_seurat_object <- FindNeighbors(object = gem6_seurat_object, dims = 1:20)
gem6_seurat_object <- FindClusters(object = gem6_seurat_object, resolution = 1.2)
gem6_seurat_object <- RunUMAP(gem6_seurat_object, dims = 1:20)
DimPlot_scCustom(gem6_seurat_object, label = T, split.by = "condition")

sweep.res.list_gem6 <- paramSweep_v3(gem6_seurat_object, PCs = 1:20, sct = T, num.cores = 2)
sweep.stats_gem6 <- summarizeSweep(sweep.res.list_gem6, GT = FALSE)
bcmv_gem6 <- find.pK(sweep.stats_gem6)
homotypic.prop <- modelHomotypic(gem6_seurat_object$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.20*nrow(gem6_seurat_object@meta.data)) #tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
gem6_seurat_object <- doubletFinder_v3(gem6_seurat_object, PCs = 1:20, pN = 0.25, pK = 0.25, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
gem6_seurat_object@meta.data[,"CellTypes_DF"] <- gem6_seurat_object$DF.classifications_0.25_0.25_8686
gem6_seurat_object@meta.data$CellTypes_DF[which(gem6_seurat_object$pANN_0.25_0.25_8686 == "Doublet")] <- "Doublet"
DimPlot(gem6_seurat_object, group.by="CellTypes_DF", reduction="umap", pt.size=0.5, order=c("Coll.Duct.TC","Doublet"), cols=c("#66C2A5","#FFD92F","#8DA0CB","#A6D854","#E78AC3","#B3B3B3","#E5C494","black","#FC8D62"))
saveRDS(gem6_seurat_object, "gem6_seurat_withDoublets.rds")


tibble(
  cluster = gem6_seurat_object$seurat_clusters,
  Dataset = gem6_seurat_object$CellTypes_DF
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


##GEM7
gem7_data_dir <- '~/Desktop/lanser_nucseq/all_filtered_feature_bc_matrix/gem7_filtered_feature_bc_matrix'
list.files(gem7_data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
gem7_data <- Read10X(data.dir = gem7_data_dir)
gem7_seurat_object = CreateSeuratObject(gem7_data, min.cells = 3, min.genes = 350)
gem7_seurat_object@meta.data$condition <- "KO"
gem7_seurat_object@meta.data$Sex <- "F"
gem7_seurat_object@meta.data$orig.ident <- "GEM7"
gem7_seurat_object <- gem7_seurat_object[!grepl("Malat1", row.names(gem7_seurat_object)), ]
gem7_seurat_object[["percent.mt"]] <- PercentageFeatureSet(gem7_seurat_object, pattern = "^mt-")
plot1 <- FeatureScatter(gem7_seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gem7_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- VlnPlot(object = gem7_seurat_object, features = c("nFeature_RNA"), pt.size = 0.1) +NoLegend()
plot1 + plot2 + plot3
gem7_seurat_object <- subset(gem7_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 750 & percent.mt < 7)
gem7_seurat_object <- SCTransform(gem7_seurat_object, method = "glmGamPoi", vars.to.regress = "percent.mt")
gem7_seurat_object <- RunPCA(object = gem7_seurat_object, verbose = F)
ElbowPlot(object = gem7_seurat_object, ndims = 20)
gem7_seurat_object <- FindNeighbors(object = gem7_seurat_object, dims = 1:15)
gem7_seurat_object <- FindClusters(object = gem7_seurat_object, resolution = 1.2)
gem7_seurat_object <- RunUMAP(gem7_seurat_object, dims = 1:15)
DimPlot_scCustom(gem7_seurat_object, label = T, split.by = "condition")

sweep.res.list_gem7 <- paramSweep_v3(gem7_seurat_object, PCs = 1:15, sct = T, num.cores = 2)
sweep.stats_gem7 <- summarizeSweep(sweep.res.list_gem7, GT = FALSE)
bcmv_gem7 <- find.pK(sweep.stats_gem7)
homotypic.prop <- modelHomotypic(gem7_seurat_object$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.046*nrow(gem7_seurat_object@meta.data)) #tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
gem7_seurat_object <- doubletFinder_v3(gem7_seurat_object, PCs = 1:15, pN = 0.25, pK = 0.11, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
gem7_seurat_object@meta.data[,"CellTypes_DF"] <- gem7_seurat_object$DF.classifications_0.25_0.11_465
gem7_seurat_object@meta.data$CellTypes_DF[which(gem7_seurat_object$pANN_0.25_0.11_465 == "Doublet")] <- "Doublet"
DimPlot(gem7_seurat_object, group.by="CellTypes_DF", reduction="umap", pt.size=0.5, order=c("Coll.Duct.TC","Doublet"), cols=c("#66C2A5","#FFD92F","#8DA0CB","#A6D854","#E78AC3","#B3B3B3","#E5C494","black","#FC8D62"))
saveRDS(gem7_seurat_object, "gem7_seurat_withDoublets.rds")



##GEM8
gem8_data_dir <- '~/Desktop/lanser_nucseq/all_filtered_feature_bc_matrix/gem8_filtered_feature_bc_matrix'
list.files(gem8_data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
gem8_data <- Read10X(data.dir = gem8_data_dir)
gem8_seurat_object = CreateSeuratObject(gem8_data, min.cells = 3, min.genes = 350)
gem8_seurat_object@meta.data$condition <- "KO"
gem8_seurat_object@meta.data$Sex <- "M"
gem8_seurat_object@meta.data$orig.ident <- "GEM8"
gem8_seurat_object <- gem8_seurat_object[!grepl("Malat1", row.names(gem8_seurat_object)), ]
gem8_seurat_object[["percent.mt"]] <- PercentageFeatureSet(gem8_seurat_object, pattern = "^mt-")
plot1 <- FeatureScatter(gem8_seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gem8_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- VlnPlot(object = gem8_seurat_object, features = c("nFeature_RNA"), pt.size = 0.1) +NoLegend()
plot1 + plot2 + plot3
gem8_seurat_object <- subset(gem8_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 750 & percent.mt < 5)
gem8_seurat_object <- SCTransform(gem8_seurat_object, method = "glmGamPoi", vars.to.regress = "percent.mt")
gem8_seurat_object <- RunPCA(object = gem8_seurat_object, verbose = F)
ElbowPlot(object = gem8_seurat_object, ndims = 20)
gem8_seurat_object <- FindNeighbors(object = gem8_seurat_object, dims = 1:15)
gem8_seurat_object <- FindClusters(object = gem8_seurat_object, resolution = 1.2)
gem8_seurat_object <- RunUMAP(gem8_seurat_object, dims = 1:15)
DimPlot_scCustom(gem8_seurat_object, label = T, split.by = "condition")

sweep.res.list_gem8 <- paramSweep_v3(gem8_seurat_object, PCs = 1:15, sct = T, num.cores = 2)
sweep.stats_gem8 <- summarizeSweep(sweep.res.list_gem8, GT = FALSE)
bcmv_gem8 <- find.pK(sweep.stats_gem8)
homotypic.prop <- modelHomotypic(gem8_seurat_object$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.05*nrow(gem8_seurat_object@meta.data)) #tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
gem8_seurat_object <- doubletFinder_v3(gem8_seurat_object, PCs = 1:15, pN = 0.25, pK = 0.26, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
gem8_seurat_object@meta.data[,"CellTypes_DF"] <- gem8_seurat_object$DF.classifications_0.25_0.26_576
gem8_seurat_object@meta.data$CellTypes_DF[which(gem8_seurat_object$pANN_0.25_0.26_576 == "Doublet")] <- "Doublet"
DimPlot(gem8_seurat_object, group.by="CellTypes_DF", reduction="umap", pt.size=0.5, order=c("Coll.Duct.TC","Doublet"), cols=c("#66C2A5","#FFD92F","#8DA0CB","#A6D854","#E78AC3","#B3B3B3","#E5C494","black","#FC8D62"))
saveRDS(gem8_seurat_object, "gem8_seurat_withDoublets.rds")



##GEM9
gem9_data_dir <- '~/Desktop/lanser_nucseq/all_filtered_feature_bc_matrix/gem9_filtered_feature_bc_matrix'
list.files(gem9_data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
gem9_data <- Read10X(data.dir = gem9_data_dir)
gem9_seurat_object = CreateSeuratObject(gem9_data, min.cells = 3, min.genes = 350)
gem9_seurat_object@meta.data$condition <- "HPA"
gem9_seurat_object@meta.data$Sex <- "F"
gem9_seurat_object@meta.data$orig.ident <- "GEM9"
gem9_seurat_object <- gem9_seurat_object[!grepl("Malat1", row.names(gem9_seurat_object)), ]
gem9_seurat_object[["percent.mt"]] <- PercentageFeatureSet(gem9_seurat_object, pattern = "^mt-")
plot1 <- FeatureScatter(gem9_seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gem9_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- VlnPlot(object = gem9_seurat_object, features = c("nFeature_RNA"), pt.size = 0.1) +NoLegend()
plot1 + plot2 + plot3
gem9_seurat_object <- subset(gem9_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 1000 & percent.mt < 5)
gem9_seurat_object <- SCTransform(gem9_seurat_object, method = "glmGamPoi", vars.to.regress = "percent.mt")
gem9_seurat_object <- RunPCA(object = gem9_seurat_object, verbose = F)
ElbowPlot(object = gem9_seurat_object, ndims = 50)
gem9_seurat_object <- FindNeighbors(object = gem9_seurat_object, dims = 1:15)
gem9_seurat_object <- FindClusters(object = gem9_seurat_object, resolution = 2.0)
gem9_seurat_object <- RunUMAP(gem9_seurat_object, dims = 1:15)
DimPlot_scCustom(gem9_seurat_object, label = T, split.by = "condition")

sweep.res.list_gem9 <- paramSweep_v3(gem9_seurat_object, PCs = 1:15, sct = T, num.cores = 2)
sweep.stats_gem9 <- summarizeSweep(sweep.res.list_gem9, GT = FALSE)
bcmv_gem9 <- find.pK(sweep.stats_gem9)
homotypic.prop <- modelHomotypic(gem9_seurat_object$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.20*nrow(gem9_seurat_object@meta.data)) #tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
gem9_seurat_object <- doubletFinder_v3(gem9_seurat_object, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
gem9_seurat_object@meta.data[,"CellTypes_DF"] <- gem9_seurat_object$DF.classifications_0.25_0.09_9677
gem9_seurat_object@meta.data$CellTypes_DF[which(gem9_seurat_object$pANN_0.25_0.09_9677 == "Doublet")] <- "Doublet"
DimPlot(gem9_seurat_object, group.by="CellTypes_DF", reduction="umap", pt.size=0.5, order=c("Coll.Duct.TC","Doublet"), cols=c("#66C2A5","#FFD92F","#8DA0CB","#A6D854","#E78AC3","#B3B3B3","#E5C494","black","#FC8D62"))
saveRDS(gem9_seurat_object, "gem9_seurat_withDoublets.rds")

tibble(
  cluster = gem9_seurat_object$seurat_clusters,
  Dataset = gem9_seurat_object$CellTypes_DF
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


##GEM10
gem10_data_dir <- '~/Desktop/lanser_nucseq/all_filtered_feature_bc_matrix/gem10_filtered_feature_bc_matrix'
list.files(gem10_data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
gem10_data <- Read10X(data.dir = gem10_data_dir)
gem10_seurat_object = CreateSeuratObject(gem10_data, min.cells = 3, min.genes = 350)
gem10_seurat_object@meta.data$condition <- "PBS"
gem10_seurat_object@meta.data$Sex <- "M"
gem10_seurat_object@meta.data$orig.ident <- "GEM10"
gem10_seurat_object <- gem10_seurat_object[!grepl("Malat1", row.names(gem10_seurat_object)), ]
gem10_seurat_object[["percent.mt"]] <- PercentageFeatureSet(gem10_seurat_object, pattern = "^mt-")
plot1 <- FeatureScatter(gem10_seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gem10_seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- VlnPlot(object = gem10_seurat_object, features = c("nFeature_RNA"), pt.size = 0.1) +NoLegend()
plot1 + plot2 + plot3
gem10_seurat_object <- subset(gem10_seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 750 & percent.mt < 10)
gem10_seurat_object <- SCTransform(gem10_seurat_object, method = "glmGamPoi", vars.to.regress = "percent.mt")
gem10_seurat_object <- RunPCA(object = gem10_seurat_object, verbose = F)
ElbowPlot(object = gem10_seurat_object, ndims = 20)
gem10_seurat_object <- FindNeighbors(object = gem10_seurat_object, dims = 1:15)
gem10_seurat_object <- FindClusters(object = gem10_seurat_object, resolution = 1.2)
gem10_seurat_object <- RunUMAP(gem10_seurat_object, dims = 1:15)
DimPlot_scCustom(gem10_seurat_object, label = T, split.by = "condition")

sweep.res.list_gem10 <- paramSweep_v3(gem10_seurat_object, PCs = 1:15, sct = T, num.cores = 2)
sweep.stats_gem10 <- summarizeSweep(sweep.res.list_gem10, GT = FALSE)
bcmv_gem10 <- find.pK(sweep.stats_gem10)
homotypic.prop <- modelHomotypic(gem10_seurat_object$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.065*nrow(gem10_seurat_object@meta.data)) #tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
gem10_seurat_object <- doubletFinder_v3(gem10_seurat_object, PCs = 1:15, pN = 0.25, pK = 0.1, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
gem10_seurat_object@meta.data[,"CellTypes_DF"] <- gem10_seurat_object$DF.classifications_0.25_0.1_963
gem10_seurat_object@meta.data$CellTypes_DF[which(gem10_seurat_object$pANN_0.25_0.1_963 == "Doublet")] <- "Doublet"
DimPlot(gem10_seurat_object, group.by="CellTypes_DF", reduction="umap", pt.size=0.5, order=c("Coll.Duct.TC","Doublet"), cols=c("#66C2A5","#FFD92F","#8DA0CB","#A6D854","#E78AC3","#B3B3B3","#E5C494","black","#FC8D62"))
saveRDS(gem10_seurat_object, "gem10_seurat_withDoublets.rds")

gem1_Wdoublets <- readRDS("gem1_seurat_withDoublets.rds")
gem2_Wdoublets <- readRDS("gem2_seurat_withDoublets.rds")
gem3_Wdoublets <- readRDS("gem3_seurat_withDoublets.rds")
gem4_Wdoublets <- readRDS("gem4_seurat_withDoublets.rds")
gem5_Wdoublets <- readRDS("gem5_seurat_withDoublets.rds")
gem6_Wdoublets <- readRDS("gem6_seurat_withDoublets.rds")
gem7_Wdoublets <- readRDS("gem7_seurat_withDoublets.rds")
gem8_Wdoublets <- readRDS("gem8_seurat_withDoublets.rds")
gem9_Wdoublets <- readRDS("gem9_seurat_withDoublets.rds")
gem10_Wdoublets <- readRDS("gem10_seurat_withDoublets.rds")

Idents(object = gem1_Wdoublets) <- "CellTypes_DF"
Idents(object = gem2_Wdoublets) <- "CellTypes_DF"
Idents(object = gem3_Wdoublets) <- "CellTypes_DF"
Idents(object = gem4_Wdoublets) <- "CellTypes_DF"
Idents(object = gem5_Wdoublets) <- "CellTypes_DF"
Idents(object = gem6_Wdoublets) <- "CellTypes_DF"
Idents(object = gem7_Wdoublets) <- "CellTypes_DF"
Idents(object = gem8_Wdoublets) <- "CellTypes_DF"
Idents(object = gem9_Wdoublets) <- "CellTypes_DF"
Idents(object = gem10_Wdoublets) <- "CellTypes_DF"

gem1_singlets <- subset(gem1_Wdoublets, idents ="Singlet")
gem2_singlets <- subset(gem2_Wdoublets, idents ="Singlet")
gem3_singlets <- subset(gem3_Wdoublets, idents ="Singlet")
gem4_singlets <- subset(gem4_Wdoublets, idents ="Singlet")
gem5_singlets <- subset(gem5_Wdoublets, idents ="Singlet")
gem6_singlets <- subset(gem6_Wdoublets, idents ="Singlet")
gem7_singlets <- subset(gem7_Wdoublets, idents ="Singlet")
gem8_singlets <- subset(gem8_Wdoublets, idents ="Singlet")
gem9_singlets <- subset(gem9_Wdoublets, idents ="Singlet")
gem10_singlets <- subset(gem10_Wdoublets, idents ="Singlet")

saveRDS(gem1_singlets, "gem1_singlets.rds")
saveRDS(gem2_singlets, "gem2_singlets.rds")
saveRDS(gem3_singlets, "gem3_singlets.rds")
saveRDS(gem4_singlets, "gem4_singlets.rds")
saveRDS(gem5_singlets, "gem5_singlets.rds")
saveRDS(gem6_singlets, "gem6_singlets.rds")
saveRDS(gem7_singlets, "gem7_singlets.rds")
saveRDS(gem8_singlets, "gem8_singlets.rds")
saveRDS(gem9_singlets, "gem9_singlets.rds")
saveRDS(gem10_singlets, "gem10_singlets.rds")


HPA_PBS_merged <- merge(gem1_singlets, y = c(gem2_singlets, gem3_singlets, gem4_singlets,
                                             gem9_singlets, gem10_singlets),
                        add.cell.ids = c("GEM1", "GEM2", "GEM3", "GEM4",
                                         "GEM9", "GEM10"))

WT_KO_merged <- merge(gem5_singlets, y = c(gem6_singlets, gem7_singlets, gem8_singlets),
                      add.cell.ids = c("GEM5", "GEM6", "GEM7", "GEM8"))


KO_HPA_merged <- merge(gem1_singlets, y = c(gem2_singlets, gem7_singlets, gem8_singlets, gem9_singlets),
                       add.cell.ids = c("GEM1", "GEM2", "GEM7", "GEM8", "GEM9"))

WT_PBS_merged <- merge(gem3_singlets, y = c(gem4_singlets, gem5_singlets, gem6_singlets, gem10_singlets))

saveRDS(HPA_PBS_merged, "HPA_PBS_merged_unNormalized.rds")
saveRDS(WT_KO_merged, "WT_KO_merged_unNormalized.rds")
saveRDS(KO_HPA_merged, "KO_HPA_merged_unNormalized.rds")
saveRDS(WT_PBS_merged, "WT_PBS_merged_unNormalized.rds")