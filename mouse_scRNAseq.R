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

#Subset V3 postmitotic
mSC.V3post <- subset(mSC.V3, idents = c("1", "2", "3", "4", "6", "8", "11"))
Idents(mSC.V3post) <- "orig.ident"
VlnPlot(mSC.V3post, features = c("Sim1", "Nkx2-2", "Robo3", "Cntn2", "Olig3", "Lhx1", "Slc17a6"),
        cols = "royalblue2", ncol = 4)
FeatureScatter(mSC.V3post, feature1 = "Olig3", feature2 = "Lhx1", pt.size = 2,
               group.by = "orig.ident") + NoLegend()

#Subset V2
mSC.V2 <- subset(mSC,
                 subset = (Foxn4 > 0 | Vsx1 > 0 | Vsx2 > 0 | Gata2 > 0 | Gata3 > 0) & 
                   ( Tubb3 > 0 | Sox2 > 0))
mSC.V2 <- FindVariableFeatures(mSC.V2)
mSC.V2 <- ScaleData(mSC.V2, features = rownames(mSC.V2))
mSC.V2 <- RunPCA(mSC.V2)
DimHeatmap(mSC.V2, dims = 1:15, cells = 500, balanced = TRUE)
DimHeatmap(mSC.V2, dims = 16:30, cells = 500, balanced = TRUE)
DimHeatmap(mSC.V2, dims = 31:45, cells = 500, balanced = TRUE)
ElbowPlot(mSC.V2, ndims = 50)
mSC.V2 <- RunTSNE(mSC.V2, dims = 1:20)
mSC.V2 <- FindNeighbors(mSC.V2, dims = 1:20)
mSC.V2 <- FindClusters(mSC.V2, resolution = 0.6)
#Plot
DimPlot(mSC.V2, reduction = "tsne", group.by = "Emb_day", pt.size = 2)
DimPlot(mSC.V2, reduction = "tsne", label = T, label.size = 7, pt.size = 2) + NoLegend()
FeaturePlot(mSC.V2, features = "Vsx2", order = T, pt.size = 2.5)
#Find markers
write.csv(FindAllMarkers(mSC.V2, min.pct = 0.25, only.pos = T),
          file = "marker_V2all.csv", quote = F)

#Subset subtype
id_En1 <- c(rep("posi", length(WhichCells(mSC, expression = (En1 > 0 & Tubb3 > 0)))),
            rep("nega", length(WhichCells(mSC, expression = (En1 > 0 & Tubb3 > 0), invert = T))))
names(id_En1) <- c(WhichCells(mSC, expression = (En1 > 0 & Tubb3 > 0)),
                   WhichCells(mSC, expression = (En1 > 0 & Tubb3 > 0), invert = T))
mSC <- AddMetaData(mSC, metadata = id_En1, col.name = "exp_En1")
DimPlot(mSC, reduction = "tsne", group.by = "exp_En1")
Idents(mSC) <- "exp_En1"
marker_En1 <- FindMarkers(mSC, ident.1 = "posi", ident.2 = "nega", min.pct = 0.25, only.pos = T)
write.csv(marker_En1, file = "marker_En1.csv", quote = F)
DoHeatmap(mSC, cells = WhichCells(mSC, expression = En1 > 0 & Tubb3 > 0),
          features = c("En1", "Foxd3", "Lhx1", "Lhx5", "Pax2", "Pax5", "Pax8", "Sp9", "Foxp2",
                       "Foxd3", "Onecut2", "Otp", "Gbx2", "Pou2f2", "Neurod1", "Zfhx3",
                       "Sox11", "Bhlhe22", "Mab21l1",
                       "Mab21l2", "Lmo3", "Slc32a1", "Slc6a5", "Gad2",
                       "Nova2", "Rtn1",
                       "Cnr1", "Nxph1", "Hmx2", "Cacna2d3", "Pnoc",
                       "Adarb2", "Myt1l", "Nova1", "Nrxn3", "Pdzrn4", "Rab3c",
                       "Celf4", "Olfm1", "Vstm2l", "Rcan2", "Olfm2", "Chl1", "Ank3",
                       "Cirbp", "Amer2", "Ywhah", "Zc4h2", "Jagn1", "Tnrc6c",
                       "Rnasek", "Fscn1", "Celf3", "Ckb", "Syt11", "Hes6",
                       "Nrxn1", "Ssbp3"),
          label = F)

id_Evx <- c(rep("posi", length(WhichCells(mSC, expression = (Evx1 > 0 | Evx2 > 0) & Tubb3 > 0))),
            rep("nega", length(WhichCells(mSC, expression = (Evx1 > 0 | Evx2 > 0) & Tubb3 > 0, 
                                          invert = T))))
names(id_Evx) <- c(WhichCells(mSC, expression = (Evx1 > 0 | Evx2 > 0) & Tubb3 > 0),
                   WhichCells(mSC, expression = (Evx1 > 0 | Evx2 > 0) & Tubb3 > 0, invert = T))
mSC <- AddMetaData(mSC, metadata = id_Evx, col.name = "exp_Evx1_2")
Idents(mSC) <- "exp_Evx1_2"
DimPlot(mSC, reduction = "tsne")
write.csv(FindMarkers(mSC, ident.1 = "posi", ident.2 = "nega",
                      min.pct = 0.25, only.pos = T),
          file = "marker_Evx1_2.csv", quote = F)
DoHeatmap(mSC, cells = WhichCells(mSC, expression = (Evx1 > 0 | Evx2 > 0) & Tubb3 > 0),
          features = c("Evx1", "Nrn1", "Adcyap1", "Arl4d", "Elavl4", "Pou3f1",
                       "Robo3", "Nhlh2", "Slc17a6", "Sort1", "Pou2f2", "Map1b",
                       "Stmn1", "Rgmb", "Rab3c", "Gng3", "Lhx1", "Ina",
                       "Mab21l2", "Epb41", "Tubb5", "Hmgb3", "Gpm6a",
                       "Nsg2", "Ebf3", "Nova2", "Gng2", "Elavl3", "Lzts1",
                       "Mllt11", "Rtn1", "Scrt2", "Tuba1a", "Uncx",
                       "Ebf1", "Meis2", "Tmeff1", "Cadm4", "Zc4h2", "Dpysl3",
                       "Pak1", "Ptp4a1",
                       "Cbln2", "Nefl", "Pdzrn4", "Tmem163", "Olfm2",
                       "Lhx5", "Pax5", "Synpr", "Mapk10", "Lnpk",
                       "Mab21l1", "Lhx9", "Lhx2", "Fez1", "Ywhah", "Dclk2",
                       "Onecut1", "Oaz1", "Gabarap", "Ptma", "Otp"),
          label = F)
