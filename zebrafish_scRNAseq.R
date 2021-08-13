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

#subset V2
zWE.V2 <- subset(zWE.1_2SC, subset = (foxn4 > 0 | vsx1 > 0 | vsx2 > 0 | gata2a > 0 | gata3 > 0) &
                   otx1 == 0 & otx2a == 0 & otx2b == 0 & six3a == 0 & six3b == 0 & emx2 == 0)
zWE.V2 <- FindVariableFeatures(zWE.V2)
zWE.V2 <- ScaleData(zWE.V2, features = rownames(zWE.V2))
zWE.V2 <- RunPCA(zWE.V2)
DimHeatmap(zWE.V2, dims = 1:15, balanced = TRUE)
DimHeatmap(zWE.V2, dims = 16:30, balanced = TRUE)
ElbowPlot(zWE.V2, ndims = 50)
zWE.V2 <- RunTSNE(zWE.V2, dims = 1:15)
zWE.V2 <- FindNeighbors(zWE.V2, dims = 1:15)
zWE.V2 <- FindClusters(zWE.V2, resolution = 0.6)
#Plot
DimPlot(zWE.V2, reduction = "tsne", label = T, label.size = 8, pt.size = 2) + NoLegend()
DimPlot(zWE.V2, reduction = "tsne", group.by = "Age", label = TRUE,
        label.size = 8, pt.size = 2) + NoLegend()
FeaturePlot(zWE.V2, reduction = "tsne", features = "pkd2l1", order = T, pt.size = 4)
#Find markers of V2
write.csv(FindAllMarkers(zWE.V2, min.pct = 0.25, only.pos = T),
          file = "marker_zV2nonotx.csv", quote = F)

#Subset subtypes
id_sim <- c(rep("posi", length(WhichCells(zWE.1_2SC, expression = sim1a > 0))),
            rep("nega", length(WhichCells(zWE.1_2SC, expression = sim1a > 0, invert = T))))
names(id_sim) <- c(WhichCells(zWE.1_2SC, expression = sim1a > 0),
                   WhichCells(zWE.1_2SC, expression = sim1a > 0, invert = T))
zWE.1_2SC <- AddMetaData(zWE.1_2SC, metadata = id_sim, col.name = "exp_sim1a")
DimPlot(zWE.1_2SC, reduction = "tsne", group.by = "exp_sim1a")
Idents(zWE.1_2SC) <- "exp_sim1a"
VlnPlot(zWE.1_2SC, features = c("sim1a", "nkx2.2a", "nkx2.2b", "robo3", "cntn2",
                                "lhx1a", "olig3", "slc17a6b"),
        idents = "posi", cols = "yellow2")

id_en1 <- c(rep("posi", length(WhichCells(zWE.1_2SC, expression = en1b > 0))),
            rep("nega", length(WhichCells(zWE.1_2SC, expression = en1b > 0, invert = T))))
names(id_en1) <- c(WhichCells(zWE.1_2SC, expression = en1b > 0),
                   WhichCells(zWE.1_2SC, expression = en1b > 0, invert = T))
zWE.1_2SC <- AddMetaData(zWE.1_2SC, metadata = id_en1, col.name = "exp_en1b")
DimPlot(zWE.1_2SC, reduction = "tsne", group.by = "exp_en1b")
Idents(zWE.1_2SC) <- "exp_en1b"
marker_en1b <- FindMarkers(zWE.1_2SC, ident.1 = "posi", ident.2 = "nega",
                           min.pct = 0.25, only.pos = T)
write.csv(marker_en1b, file = "marker_en1b.csv", quote = F)
DoHeatmap(zWE.1_2SC, cells = WhichCells(zWE.1_2SC, expression = en1b > 0),
          features = c("en1b", "foxd3", "lhx1a", "lhx5", "pax2a", "pax5", "pax8",
                       "sp9", "foxp2", "foxd3", "onecut2", "otpa", "otpb", "gbx2",
                       "pou2f2a.1", "neurod1", "zfhx3", 
                       "sox11a", "bhlhe22", "mab21l1", "mab21l2",
                       "lmo3", "slc32a1", "slc6a5", "gad2", "nova2", "rtn1a",
                       "cnr1", "nxph1", "hmx2", "cacna2d3", "pnoca",
                       "adarb2", "myt1la", "nova1", "nrxn3b",
                       "pdzrn4", "rab3c", "celf4",
                       "olfm1a", "olfm1b", "vstm2l", "rcan2", "olfm2a", "chl1a",
                       "ank3b",
                       "cirbpb", "amer2", "ywhah", "zc4h2", "jagn1a", "tnrc6c1",
                       "rnasekb", "fscn1a", "celf3a", "ckbb", "syt11a", "hes6",
                       "nrxn1a", "ssbp3b"),
          label = F)

id_evx <- c(rep("posi", length(WhichCells(zWE.1_2SC, expression = evx2 > 0))),
            rep("nega", length(WhichCells(zWE.1_2SC, expression = evx2 > 0, invert = T))))
names(id_evx) <- c(WhichCells(zWE.1_2SC, expression = evx2 > 0),
                   WhichCells(zWE.1_2SC, expression = evx2 > 0, invert = T))
zWE.1_2SC <- AddMetaData(zWE.1_2SC, metadata = id_evx, col.name = "exp_evx2")
DimPlot(zWE.1_2SC, reduction = "tsne", group.by = "exp_evx2")
Idents(zWE.1_2SC) <- "exp_evx2"
write.csv(FindMarkers(zWE.1_2SC, ident.1 = "posi", ident.2 = "nega",
                      min.pct = 0.25, only.pos = T),
          file = "marker_evx2.csv", quote = F)
DoHeatmap(zWE.1_2SC, cells = WhichCells(zWE.1_2SC, expression = evx2 > 0),
          features = c("evx2", "nrn1a", "adcyap1b", "arl4d", "elavl4", "pou3f1",
                       "robo3", "nhlh2", "slc17a6b", "sort1a", "pou2f2a.1",
                       "map1b", "stmn2b", "rgmb", "rab3c",
                       "gng3", "lhx1a", "inab", "mab21l2", "epb41a", "tubb5",
                       "hmgb3a", "gpm6ab", "nsg2", "ebf3a", "nova2", "gng2",
                       "elavl3", "lzts1", "mllt11", "rtn1a",
                       "scrt2", "tuba1a", "uncx", "ebf1a", "meis2a", "tmeff1b",
                       "cadm4", "zc4h2", "dpysl3", "pak1", "ptp4a1",
                       "cbln2b", "neflb", "pdzrn4", "tmem163a", "olfm2a",
                       "lhx5", "pax5", "synpr", "mapk10", "lnpk", "mab21l1",
                       "lhx9", "lhx2b", "fez1", "ywhah",
                       "dclk2a", "onecut1", "oaz1b",
                       "gabarapb", "ptmab", "otpa"),
          label = F)
