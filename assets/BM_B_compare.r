# add selection metadata
## do the same with Alin_neg and B_BM_tot
#rename everything to BM.AB.merge...

Blin.neg[["selection"]] = "lin.neg"
B.BM.tot[["selection"]] = "total"

#merge both seurat ojects

BM.B.merge = merge(B.BM.tot,Blin.neg)

#split the objects

BM.B.split <- SplitObject(BM.B.merge, split.by = "selection")

# Perform integration
# We then identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData.
BM.B.anchors <- FindIntegrationAnchors(object.list = BM.B.split, dims = 1:20)
BM.B.combined <- IntegrateData(anchorset = BM.B.anchors, dims = 1:20)

# Perform an integrated analysis
# Now we can run a single integrated analysis on all cells!

DefaultAssay(BM.B.combined) <- "integrated"
BM.B.combined <- ScaleData(BM.B.combined)
BM.B.combined <- RunPCA(BM.B.combined)
BM.B.combined <- RunUMAP(BM.B.combined, , dims = 1:20)
BM.B.combined <- RunTSNE(BM.B.combined, , dims = 1:20)
BM.B.combined = FindNeighbors(BM.B.combined, dims = 1:20)
BM.B.combined = FindClusters(BM.B.combined, resolution = 0.5)

# visualisation cluster plot by sample
p1 <- DimPlot(BM.B.combined, reduction = "tsne", group.by = "selection")
p2 <- DimPlot(BM.B.combined, reduction = "tsne", label = TRUE)
plot_grid(p1, p2)

# by cluster
DimPlot(BM.B.combined, reduction = "tsne", split.by = "selection")

#Identify conserved cell type markers

BM_all.markers <- FindAllMarkers(BM.B.combined, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, grouping.var = "selection")
##doesn't split samples

BM.B.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

topgenes = c("HIST3H2A","PRTN3", "PCNA", "S100A9", "FTH1", "AHSP","TOP2A","ID2","CCDC50","IGHA1")
top2genes = c("SLC25A4","CTSG", "LYZ", "CXCL8", "FCN1", "HBB","CENPE","HLA-DPB1","TCF4","IGKC")

DotPlot(BM.B.combined, features = rev(topgenes), cols = c("blue", "red"), dot.scale = 8, split.by = "selection") + RotatedAxis()

#Identify differential expressed genes across conditions
#find sample-specific genes in each cluster

theme_set(theme_cowplot())
BM.B.combined_cluster0 <- subset(BM.B.combined, idents = 0)
Idents(BM.B.combined_cluster0) <- "selection"
avg.BM.B.combined_cluster0 <- log1p(AverageExpression(BM.B.combined_cluster0)$RNA)
avg.BM.B.combined_cluster0$gene <- rownames(BM.B.combined_cluster0)

BM.B.combined_cluster3 <- subset(BM.B.combined, idents = 3)
Idents(BM.B.combined_cluster3) <- "selection"
avg.BM.B.combined_cluster3 <- log1p(AverageExpression(BM.B.combined_cluster3)$RNA)
avg.BM.B.combined_cluster3$gene <- rownames(avg.BM.B.combined_cluster3)

genes.to.label = c("HIST3H2A", "CD14", "SLC25A4", "PCNA", "LYZ")
p1 <- ggplot(avg.BM.B.combined_cluster0, aes(total, lin.neg)) + geom_point() + ggtitle("BM.B.combined_cluster0")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)
p2 <- ggplot(avg.BM.B.combined_cluster3, aes(total, lin.neg)) + geom_point() + ggtitle("BM.B.combined_cluster3")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)
plot_grid(p1, p2)

# find significant different genes between total and lin_neg

BM.B.combined$celltype.selection <- paste(Idents(BM.B.combined), BM.B.combined$selection, sep = "_")
BM.B.combined$celltype <- Idents(BM.B.combined)
Idents(BM.B.combined) <- "celltype.selection"

BM.B_cluster0_compare <- FindMarkers(BM.B.combined, ident.1 = "0_lin.neg", ident.2 = "0_total")
head(BM.B_cluster0_compare, n = 15)

#pick the top 6 genes
cluster_0_genes = c("DEFA3","DEFA4","IGKC","CTSG","IGLC2","GCC2")
# cluster plot
FeaturePlot(BM.B.combined, features = cluster_0_genes, reduction = "tsne",split.by = "selection", max.cutoff = 3, cols = c("grey", "red"))
FeaturePlot(BM.B.combined, features = markers.to.plot, split.by = "selection", reduction = "tsne", do.hover = TRUE, data.hover = c("ident", "PC1", "nFeature"), max.cutoff = 3, cols = c("grey", "red"))

##violin plot
plots = VlnPlot(BM.B.combined, features = cluster_0_6genes, split.by = "selection", group.by = "seurat_clusters" ,pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 2,legend = "none")

#Identify differential expressed genes across conditions
#find sample-specific genes in each cluster but I changed to idents to 0_lin.neg ...

theme_set(theme_cowplot())
BM.B.combined_cluster1 <- subset(BM.B.combined, idents = c("1_lin.neg","1_total"))
Idents(BM.B.combined_cluster1) <- "selection"
BM.B.combined_cluster1 <- log1p(AverageExpression(BM.B.combined_cluster1)$RNA)
BM.B.combined_cluster1$gene <- rownames(BM.B.combined_cluster1)
write.table(BM.B.combined_cluster4,"BM.B.combined_cluster4.txt")

## cigar plot
ggplot(BM.B.combined_cluster1, aes(total, lin.neg)) + geom_point() + ggtitle("BM.B.combined_cluster1")
# once I know the outliers:
p1 <- ggplot(BM.B.combined_cluster1, aes(total, lin.neg)) + geom_point() + ggtitle("BM.B.combined_cluster1")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)