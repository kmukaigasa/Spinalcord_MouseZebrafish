library(Seurat)
#Read 10X
mSC_noNorm.data <- Read10X(data.dir = "./filtered_feature_bc_matrix")
mSC <- CreateSeuratObject(counts = mSC_noNorm.data, project = "mSC_scRNA",
                          min.cells = 0, min.features = 0)
#Filtering
VlnPlot(mSC, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(mSC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
mSC <- subset(mSC, subset = nFeature_RNA > 200 & nCount_RNA < 90000)
#Add sample metadata
gemgroup <- sapply(strsplit(rownames(mSC@meta.data), split="-"), "[[", 2)
names(gemgroup) <- rownames(mSC@meta.data)
gemgroup <- plyr::mapvalues(
  x = gemgroup,
  from = 1:12,
  to = c(
    rep("E9.5", length=2),
    rep("E10.5", length=2),
    rep("E11.5", length=2),
    rep("E12.5", length=3),
    rep("E13.5", length=3))
)
mSC <- AddMetaData(mSC, metadata = gemgroup, col.name = 'Emb_day')
mSC@meta.data$Emb_day <- factor(mSC@meta.data$Emb_day,
                                   levels = c("E9.5", "E10.5", "E11.5", "E12.5", "E13.5"))
#Normalization
mSC <- NormalizeData(mSC)
mSC <- FindVariableFeatures(mSC, selection.method = "vst", nfeatures = 2000)
mSC <- ScaleData(mSC, features = rownames(mSC))
#Dimentinality reduction & clustering
mSC <- RunPCA(mSC, npcs = 100)
mSC <- RunTSNE(mSC, dims = 1:50)
mSC <- FindNeighbors(mSC, dims = 1:50)
mSC <- FindClusters(mSC, resolution = 1)
#Plot
DimPlot(mSC, reduction = "tsne", group.by = "Emb_day")
DimPlot(mSC, reduction = "tsne", label = TRUE) + NoLegend()
FeaturePlot(mSC, features = "Sim1")

#Subset V3
mSC.V3 <- subset(mSC, subset = (Tubb3 > 0 | Sox2 > 0) & (`Nkx2-2` > 0 | Sim1 > 0))
#Normalization
mSC.V3 <- NormalizeData(mSC.V3)
mSC.V3 <- FindVariableFeatures(mSC.V3)
mSC.V3 <- ScaleData(mSC.V3, features = rownames(mSC.V3))
#Dimensionality reduction & clustering
mSC.V3 <- RunPCA(mSC.V3)
DimHeatmap(mSC.V3, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(mSC.V3, dims = 16:30, cells = 500, balanced = TRUE)
DimHeatmap(mSC.V3, dims = 31:45, cells = 500, balanced = TRUE)
ElbowPlot(mSC.V3, ndims = 50)
mSC.V3 <- RunTSNE(mSC.V3, dims = 1:30)
mSC.V3 <- FindNeighbors(mSC.V3, dims = 1:30)
mSC.V3 <- FindClusters(mSC.V3, resolution = 0.8)
#Plot
DimPlot(mSC.V3, reduction = "tsne", group.by = "Emb_day", pt.size = 2)
DimPlot(mSC.V3, reduction = "tsne", label = TRUE, pt.size = 2, label.size = 6) + NoLegend()
FeaturePlot(mSC.V3, features = "Robo3", pt.size = 2)

#Find markers
marker.all <- FindAllMarkers(mSC.V3, min.pct = 0.25, only.pos = TRUE)
write.csv(marker.all, file = "marker.all.csv", quote = FALSE)
