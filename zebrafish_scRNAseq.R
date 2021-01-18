library(Seurat)
zWE <- readRDS("Farnsworth_2019.rds")
zWE <- UpdateSeuratObject(zWE)
#Add metadara
md <- zWE[["sample_name"]]
md$sample_name <- plyr::mapvalues(
  x = md$sample_name,
  from = c("1_olig2_24a", "2_olig2_24b", "3_olig2_48a", "4_elav_48b", "5_olig2120a", "6_olig2120b"),
  to = c("1dpf", "1dpf", "2dpf", "2dpf", "5dpf", "5dpf")
)
zWE <- AddMetaData(zWE, metadata = md, col.name = 'Age')

#Subset 1-2dpf whole embryo
zWE.1_2dpf <- subset(zWE, Age == "1dpf" | Age == "2dpf")
#Scaling
zWE.1_2dpf <- FindVariableFeatures(zWE.1_2dpf)
zWE.1_2dpf <- ScaleData(zWE.1_2dpf, features = rownames(zWE.1_2dpf))
#Dimentionality reduction
zWE.1_2dpf <- RunPCA(zWE.1_2dpf, npcs = 100)
DimHeatmap(zWE.1_2dpf, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(zWE.1_2dpf, dims = 16:30, cells = 500, balanced = TRUE)
DimHeatmap(zWE.1_2dpf, dims = 31:45, cells = 500, balanced = TRUE)
DimHeatmap(zWE.1_2dpf, dims = 46:60, cells = 500, balanced = TRUE)
DimHeatmap(zWE.1_2dpf, dims = 61:75, cells = 500, balanced = TRUE)
DimHeatmap(zWE.1_2dpf, dims = 76:90, cells = 500, balanced = TRUE)
ElbowPlot(zWE.1_2dpf, ndims = 60)
zWE.1_2dpf <- RunTSNE(zWE.1_2dpf, dims = 1:54)
zWE.1_2dpf <- FindNeighbors(zWE.1_2dpf, dims = 1:54)
zWE.1_2dpf <- FindClusters(zWE.1_2dpf, resolution = 1.2)

#Plot
DimPlot(zWE.1_2dpf, reduction = "tsne", label = TRUE) + NoLegend()
DimPlot(zWE.1_2dpf, reduction = "tsne", group.by = "Age")
FeaturePlot(zWE.1_2dpf, reduction = "tsne", features = "sox10", order = TRUE)

#Highlight spinal cord cells
DimPlot(zWE.1_2dpf, reduction = "tsne",
        cells.highlight = WhichCells(zWE.1_2dpf, idents = c("2","4", "42")),
        sizes.highlight = 0.3) + NoLegend()

#Find markers of zWE.1_2dpf
marker_zWE.1_2dpf <- FindAllMarkers(zWE.1_2dpf, min.pct = 0.25, only.pos = TRUE)
write.csv(marker_zWE.1_2dpf, file = "marker_zWE.1_2dpf.csv", quote = FALSE)

#Subset spinal cord from zWE.1_2dpf
zWE.1_2SC <- subset(zWE.1_2dpf, idents = c("2", "4", "42"))

#Scaling
zWE.1_2SC <- FindVariableFeatures(zWE.1_2SC)
zWE.1_2SC <- ScaleData(zWE.1_2SC, features = rownames(zWE.1_2SC))
#Dimentionality reduction
zWE.1_2SC <- RunPCA(zWE.1_2SC)
DimHeatmap(zWE.1_2SC, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(zWE.1_2SC, dims = 16:30, cells = 500, balanced = TRUE)
DimHeatmap(zWE.1_2SC, dims = 31:45, cells = 500, balanced = TRUE)
ElbowPlot(zWE.1_2SC, ndims = 50)
zWE.1_2SC <- RunTSNE(zWE.1_2SC, dims = 1:28)
zWE.1_2SC <- FindNeighbors(zWE.1_2SC, dims = 1:28)
zWE.1_2SC <- FindClusters(zWE.1_2SC, resolution = 1.2)

#Plot
DimPlot(zWE.1_2SC, reduction = "tsne", label = TRUE, label.size = 8, pt.size = 1.2) + NoLegend()
DimPlot(zWE.1_2SC, reduction = "tsne", group.by = "Age", label = TRUE,
        label.size = 8, pt.size = 1.2) + NoLegend()
FeaturePlot(zWE.1_2SC, reduction = "tsne", features = "elavl3", order = TRUE, pt.size = 1.2)

#Highlight progenitor and post-mitotic cells
prog <- WhichCells(zWE.1_2SC, idents = c("0", "3", "4", "7", "11", "15", "17"))
pm <- WhichCells(zWE.1_2SC, idents = c("0", "3", "4", "7", "11", "15", "17"), invert = T)
DimPlot(zWE.1_2SC, reduction = "tsne",
        cells.highlight = list(prog, pm),
        cols.highlight = "pink",
        cols = "lightblue") + NoLegend()
DimPlot(zWE.1_2SC, reduction = "tsne", cells.highlight = prog, cols.highlight = "green") + NoLegend()

#Find markers of zWE.1_2SC
marker_zWE.1_2SC <- FindAllMarkers(zWE.1_2SC, min.pct = 0.25, only.pos = TRUE)
write.csv(marker_zWE.1_2SC, file = "marker_zWE.1_2SC.csv", quote = FALSE)

#Subset progenitor cells
zWE.1_2SCp <- subset(zWE.1_2SC, idents = c("0", "3", "4", "7", "11", "15", "17"))
FeaturePlot(zWE.1_2SCp, reduction = "tsne", features = "pax6b", order = TRUE, pt.size = 1.3)
