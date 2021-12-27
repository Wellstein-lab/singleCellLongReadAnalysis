#compare BM samples

#Integration and Label Transfer
#https://satijalab.org/seurat/v3.1/integration.html

#Tutorial: Integrating stimulated vs. control PBMC datasets to learn cell-type specific responses
#https://satijalab.org/seurat/v3.1/immune_alignment.html


library("dplyr", lib.loc="~/Library/R/3.6/library")
library("Seurat")

## import data
Alin_neg_data = Read10X(data.dir = "~/dropbox/Single_cell_analysis/Cellranger_output/Alin_neg/outs/raw_feature_bc_matrix")
Blin_neg_data = Read10X(data.dir = "~/dropbox/Single_cell_analysis/Cellranger_output/Blin_neg/outs/raw_feature_bc_matrix")
B_BM_tot_data = Read10X(data.dir = "~/dropbox/Single_cell_analysis/Cellranger_output/B_BM_tot/outs/raw_feature_bc_matrix")

##make Seurat
Alin_neg_data = CreateSeuratObject(counts = Alin_neg_data, project = "Alin_neg", min.cells = 5, min.features = 200)
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
Alin_neg_subset_200_25 <- subset(Alin_neg_data, subset = nFeature_RNA > 200  & percent.mt < 25)
FeatureScatter(Alin_neg_subset_200_25, feature1 = "nCount_RNA", feature2 = "percent.mt")
VlnPlot(Alin_neg_subset_200_25, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
head(Alin_neg_subset_200_25@meta.data)
#                 orig.ident nCount_RNA nFeature_RNA percent.mt percent.RPS percent.RPL
#AAACGAAAGAGGGCGA   Alin_neg       1326          595  16.289593   14.555053   19.834087
#AAAGAACTCTGCGTCT   Alin_neg        690          458  10.000000    7.536232    9.130435
#AAAGGGCCAAATGATG   Alin_neg       1298          327   4.699538    2.927581    5.161787
#AAAGGTAAGATCCCGC   Alin_neg        383          252  23.759791    6.005222    5.744125
#AAAGGTAGTCCGAAGA   Alin_neg        538          337  11.710037    8.921933   14.126394
#AACAACCAGAGAATCT   Alin_neg        604          363  17.384106    5.629139    9.933775
head(Alin_neg_data@meta.data)
#                 orig.ident nCount_RNA nFeature_RNA percent.mt percent.RPS percent.RPL
#AAACCCATCCGGGACT   Alin_neg        777          345  44.015444    1.158301    4.375804
#AAACGAAAGAGGGCGA   Alin_neg       1326          595  16.289593   14.555053   19.834087
#AAACGCTTCTCCCATG   Alin_neg        400          255  28.000000    7.000000    8.250000
#AAAGAACTCTGCGTCT   Alin_neg        690          458  10.000000    7.536232    9.130435
#AAAGGGCCAAATGATG   Alin_neg       1298          327   4.699538    2.927581    5.161787
#AAAGGTAAGATCCCGC   Alin_neg        383          252  23.759791    6.005222    5.744125

## print out metadata
Alin_neg_metadata = Alin_neg_data@meta.data
Alin_neg_subset_metadata = Alin_neg_subset_200_25@meta.data

## normalize
Alin_neg_subset_200_25<- NormalizeData(Alin_neg_subset_200_25, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features (feature selection) select 2000
Alin_neg_subset_200_25 <- FindVariableFeatures(Alin_neg_subset_200_25, selection.method = "vst", nfeatures = 2000)

### Identify the 25 most highly variable genes
Alin_neg_subset_200_25_top25 <- head(VariableFeatures(Alin_neg_subset_200_25), 25)

# plot variable features with and without labels
plot1 = VariableFeaturePlot(Alin_neg_subset_200_25)
plot2 <- LabelPoints(plot = plot1, points = Alin_neg_subset_200_25_top25, repel = TRUE, xnudge = 0, ynudge = 0)
CombinePlots(plots = list(plot1, plot2))

#Scaling the data
all.genes <- rownames(Alin_neg_subset_200_25)
Alin_neg_subset_200_25 <- ScaleData(Alin_neg_subset_200_25, features = all.genes, vars.to.regress = "percent.mt")

#Perform linear dimensional reduction

Alin_neg_subset_200_25 <- RunPCA(Alin_neg_subset_200_25, features = VariableFeatures(object = Alin_neg_subset_200_25))

# Visualizations:
# graphs with genes vs PC for each PC
VizDimLoadings(Alin_neg_subset_200_25, dims = 1:15, reduction = "pca")

#scatter blot with 2 dimensions
DimPlot(Alin_neg_subset_200_25, reduction = "pca")

#heatmaps
DimHeatmap(Alin_neg_subset_200_25, dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)

#Determine the ‘dimensionality’ of the dataset, We identify ‘significant’ PCs as those 
#who have a strong enrichment of low p-value features.

ElbowPlot(Alin_neg_subset_200_25)

#Determine the 'dimensionality' of the dataset, We identify 'significant' PCs as those 
#who have a strong enrichment of low p-value features.

Alin_neg_subset_200_25 <- JackStraw(Alin_neg_subset_200_25, num.replicate = 100)

ScoreJackStraw(Alin_neg_subset_200_25)
#An object of class Seurat 
#12947 features across 489 samples within 1 assay 
#Active assay: RNA (12947 features)
# 3 dimensional reductions calculated: pca, umap, tsne
 
 JackStrawPlot(Alin_neg_subset_200_25, dims = 1:15)

#Cluster the cells
## PCA
Alin_neg_subset_200_25 <- FindNeighbors(Alin_neg_subset_200_25, dims = 1:10)
Alin_neg_subset_200_25 <- FindClusters(Alin_neg_subset_200_25, resolution = 0.5)

#Number of nodes: 489
#Number of edges: 12634

##UMAP
Alin_neg_subset_200_25 <- RunUMAP(Alin_neg_subset_200_25, dims = 1:10)
DimPlot(Alin_neg_subset_200_25, reduction = "umap")

##tSNE
Alin_neg_subset_200_25 <- RunTSNE(Alin_neg_subset_200_25, dims = 1:10)
DimPlot(Alin_neg_subset_200_25, reduction = "tsne")


# find markers for every cluster compared to all remaining cells, report only the positive ones
Alin_neg.markers <- FindAllMarkers(Alin_neg_subset_200_25, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
Alin_neg.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
B_BM_tot.markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_logFC) -> B_BM_tot.top4_FC_markers
write.csv(B_BM_tot.top4_FC_markers, "B_BM_tot.top4_FC_markers.csv", row.names = FALSE)

# heatmap of clusters
top10 <- Blin_neg.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
View(top10)
DoHeatmap(Alin_neg_subset_200_25, features = top10$gene) + NoLegend()

##Violinplot

VlnPlot(Alin_neg_subset_200_25, features = c("LYZ", "SRGN", "CD14"))

##Cluster graph

FeaturePlot(Alin_neg_subset_200_25, features = c("LYZ", "SRGN"), blend = TRUE, reduction = "tsne")

FeaturePlot(Alin_neg_subset_200_25, features = c("LYZ"), reduction = "tsne")

FeaturePlot(Alin_neg_subset_200_25, features = c("LYZ","FCN1", "HLA-DRB5", "SPINK2", "HBB", "SLC25A4"), reduction = "tsne")

# extract gene and UMI numbers:
quantile(B_BM_tot_subset_200_25_upd$nCount_RNA)
quantile(B_BM_tot_subset_200_25_upd$nFeature_RNA)
quantile(Blin_neg_subset_200_25_upd$nCount_RNA)
quantile(Blin_neg_subset_200_25_upd$nFeature_RNA)
quantile(Alin_neg_subset_200_25_upd$nFeature_RNA)
quantile(Alin_neg_subset_200_25_upd$nCount_RNA)
quantile(BM.AB.combined_upd$nCount_RNA)
quantile(BM.AB.combined_upd$nFeature_RNA)
