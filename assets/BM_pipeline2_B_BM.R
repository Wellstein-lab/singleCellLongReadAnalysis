#compare BM samples

#Integration and Label Transfer
#https://satijalab.org/seurat/v3.1/integration.html

#Tutorial: Integrating stimulated vs. control PBMC datasets to learn cell-type specific responses
#https://satijalab.org/seurat/v3.1/immune_alignment.html


library("dplyr", lib.loc="~/Library/R/3.6/library")
library("Seurat")

## import data
Blin_neg_data = Read10X(data.dir = "~/dropbox/Single_cell_analysis/Cellranger_output/Blin_neg/outs/raw_feature_bc_matrix")
B_BM_tot_data = Read10X(data.dir = "~/dropbox/Single_cell_analysis/Cellranger_output/B_BM_tot/outs/raw_feature_bc_matrix")

##make Seurat
Blin_neg_data = CreateSeuratObject(counts = Blin_neg_data, project = "Blin_neg", min.cells = 5, min.features = 200)
B_BM_tot_data = CreateSeuratObject(counts = B_BM_tot_data, project = "B_BM_tot", min.cells = 5, min.features = 200)

## add metadata Mitochondria and ribosomal RNA
B_BM_tot_data[["percent.MT"]] <- PercentageFeatureSet(B_BM_tot_data, pattern = "^MT")
B_BM_tot_data[["percent.RPS"]] <- PercentageFeatureSet(B_BM_tot_data, pattern = "^RPS")
B_BM_tot_data[["percent.RPL"]] <- PercentageFeatureSet(B_BM_tot_data, pattern = "^RPL")
VlnPlot(B_BM_tot_data, features = c("percent.RPS", "percent.RPL", "percent.MT"), ncol = 3)
VlnPlot(B_BM_tot_data, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
B_BM_tot_data_nCount_nFeature <- FeatureScatter(B_BM_tot_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(Alin_neg_data_nCount_nFeature, Blin_neg_data_nCount_nFeature,B_BM_tot_data_nCount_nFeature))
FeatureScatter(B_BM_tot_data, feature1 = "nCount_RNA", feature2 = "percent.MT")

## filter out cells with more than 25% MT and less than 200 counts (UMIs)
Blin_neg_subset_200_25 <- subset(Blin_neg_data, subset = nFeature_RNA > 200  & percent.mt < 25)
FeatureScatter(Blin_neg_subset_200_25, feature1 = "nCount_RNA", feature2 = "percent.mt")
VlnPlot(Blin_neg_subset_200_25, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


## print out metadata
Blin_neg_metadata = Blin_neg_data@meta.data
Blin_neg_subset_metadata = Blin_neg_subset_200_25@meta.data

## normalize
Blin_neg_subset_200_25<- NormalizeData(Blin_neg_subset_200_25, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features (feature selection) select 2000
Blin_neg_subset_200_25 <- FindVariableFeatures(Blin_neg_subset_200_25, selection.method = "vst", nfeatures = 2000)

### Identify the 25 most highly variable genes
Blin_neg_subset_200_25_top25 <- head(VariableFeatures(Blin_neg_subset_200_25), 25)

# plot variable features with and without labels
plot1 = VariableFeaturePlot(Blin_neg_subset_200_25)
plot2 <- LabelPoints(plot = plot1, points = Blin_neg_subset_200_25_top25, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

#Scaling the data
all.genes <- rownames(Blin_neg_subset_200_25)
Blin_neg_subset_200_25 <- ScaleData(Blin_neg_subset_200_25, features = all.genes, vars.to.regress = "percent.mt")

#Perform linear dimensional reduction

Blin_neg_subset_200_25 <- RunPCA(Blin_neg_subset_200_25, features = VariableFeatures(object = Blin_neg_subset_200_25))

# Visualizations:
# graphs with genes vs PC for each PC
VizDimLoadings(Blin_neg_subset_200_25, dims = 1:15, reduction = "pca")

#scatter blot with 2 dimensions
DimPlot(Blin_neg_subset_200_25, reduction = "pca")

#heatmaps
DimHeatmap(Blin_neg_subset_200_25, dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)

#Determine the ‘dimensionality’ of the dataset, We identify ‘significant’ PCs as those 
#who have a strong enrichment of low p-value features.

ElbowPlot(Blin_neg_subset_200_25)

#Determine the 'dimensionality' of the dataset, We identify 'significant' PCs as those 
#who have a strong enrichment of low p-value features.

Blin_neg_subset_200_25 <- JackStraw(Blin_neg_subset_200_25, num.replicate = 100)

ScoreJackStraw(Blin_neg_subset_200_25)
#An object of class Seurat 
#12947 features across 489 samples within 1 assay 
#Active assay: RNA (12947 features)
# 3 dimensional reductions calculated: pca, umap, tsne
 
 JackStrawPlot(Blin_neg_subset_200_25, dims = 1:15)
#Error in JackStrawPlot(Blin_neg_subset_200_25, dims = 1:3) : 
# Jackstraw procedure not scored for all the provided dims. Please run ScoreJackStraw.


#Cluster the cells
## PCA
Blin_neg_subset_200_25 <- FindNeighbors(Blin_neg_subset_200_25, dims = 1:10)
Blin_neg_subset_200_25 <- FindClusters(Blin_neg_subset_200_25, resolution = 0.5)

#Number of nodes: 489
#Number of edges: 12634

##UMAP
Blin_neg_subset_200_25 <- RunUMAP(Blin_neg_subset_200_25, dims = 1:10)
DimPlot(Blin_neg_subset_200_25, reduction = "umap")

##tSNE
Blin_neg_subset_200_25 <- RunTSNE(Blin_neg_subset_200_25, dims = 1:10)
DimPlot(Blin_neg_subset_200_25, reduction = "tsne")


# find markers for every cluster compared to all remaining cells, report only the positive ones
Blin_neg.markers <- FindAllMarkers(Blin_neg_subset_200_25, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
Blin_neg.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# heatmap of clusters
top10 <- Blin_neg.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
View(top10)
DoHeatmap(Blin_neg_subset_200_25, features = top10$gene) + NoLegend()

##Violinplot

VlnPlot(Blin_neg_subset_200_25, features = c("LYZ", "SRGN", "CD14"))

##Cluster graph

FeaturePlot(Blin_neg_subset_200_25, features = c("LYZ", "SRGN"), blend = TRUE, reduction = "tsne")

FeaturePlot(Blin_neg_subset_200_25, features = c("LYZ"), reduction = "tsne")

FeaturePlot(Blin_neg_subset_200_25, features = c("LYZ","FCN1", "HLA-DRB5", "SPINK2", "HBB", "SLC25A4"), reduction = "tsne")

# extract gene and UMI numbers:
quantile(B_BM_tot_subset_200_25_upd$nCount_RNA)
quantile(B_BM_tot_subset_200_25_upd$nFeature_RNA)
quantile(Blin_neg_subset_200_25_upd$nCount_RNA)
quantile(Blin_neg_subset_200_25_upd$nFeature_RNA)
quantile(BM.B.combined_upd$nCount_RNA)
quantile(BM.B.combined_upd$nFeature_RNA)
