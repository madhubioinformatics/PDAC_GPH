
#Autophagy genes in PDAC and Fibroblast (FB)
#Clustering and non-linear dimensionality reduction
seurat_obj <- readRDS("processed_seurat_object.rds")

#Subset of PDAC from Seurat object
PDAC <-subset(seurat_obj, subset = CellType == "PDAC")
head(PDAC)

#Subset of Fibroblasts(FB) from Seurat object
FB<- subset(seurat_obj, subset = CellType == "Fibroblast")
head(FB)


#Subset of Fibroblasts(FB) from Seurat object
MG<- subset(seurat_obj, subset = CellType == "Macrophage")
head(MG)

#Subset of Fibroblasts(FB) from Seurat object
T<- subset(seurat_obj, subset = CellType == "T_cell")
head(T)

#Subset of Fibroblasts(FB) from Seurat object
NK<- subset(seurat_obj, subset = CellType == "NK_cell")
head(NK)

#Subset of Fibroblasts(FB) from Seurat object
Endo<- subset(seurat_obj, subset = CellType == "Endothelial")
head(Endo)


DimPlot(PDAC, reduction = "tsne", group.by = "condition", label = 'T')
DimPlot(PDAC, reduction = "umap", group.by = "condition", label = 'T')

DimPlot(FB, reduction = "tsne", group.by = "condition", label = 'T')
DimPlot(FB, reduction = "umap", group.by = "condition", label = 'T')

DimPlot(MG, reduction = "tsne", group.by = "condition", label = 'T')
DimPlot(MG, reduction = "umap", group.by = "condition", label = 'T')

DimPlot(T, reduction = "tsne", group.by = "condition", label = 'T')
DimPlot(T, reduction = "umap", group.by = "condition", label = 'T')

DimPlot(NK, reduction = "tsne", group.by = "condition", label = 'T')
DimPlot(NK, reduction = "umap", group.by = "condition", label = 'T')

DimPlot(Endo, reduction = "tsne", group.by = "condition", label = 'T')
DimPlot(Endo, reduction = "umap", group.by = "condition", label = 'T')

#subset of celltypes 
T <- FindVariableFeatures(T, selection.method = "vst", nfeatures = 2000)
T <- NormalizeData(T, normalization.method = "LogNormalize", scale.factor = 10000)
T <- ScaleData(T, features = all.genes)
T <- RunPCA(T, features = VariableFeatures(object = T))
T <- RunUMAP(T, reduction = "pca", dims = 1:30)
T <- RunTSNE(T, dims = 1:30)
T <- FindNeighbors(T, reduction = "pca", dims = 1:30)
T <- FindSubCluster(T, "T_cell", "RNA_nn", subcluster.name = "T_subcluster",  resolution = 0.3, algorithm = 1)
T <- SetIdent(T, value = T@meta.data$T_subcluster)

DotPlot(T, features=c("Cd3d", "Cd3e", "Cd3g", "Cd4", "Lef1", "Igfbp4", "Ctla4", 
                      "Tnfrsf4", "Pdcd1", "Foxp3", "Cd8a", "Nkg7", "Gzmb"), dot.scale=9, cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()


T <- RenameIdents(T,
                  'T_cell_0' = 'CD4 T cell1',
                  'T_cell_1' = 'CD4 T cell2',
                  'T_cell_2' = 'CD8 Treg1',
                  'T_cell_3' = 'CD8 T cell',
                  'T_cell_4' = 'CD4 Treg1',
                  'T_cell_5' = 'CD4 Treg2',
                  'T_cell_6' = 'NK cell',
                  'T_cell_7' = 'CD8 Treg2'
)

DimPlot(T, reduction = "umap", label = TRUE, label.size = 5, pt.size = 1.5) + NoLegend()
DimPlot(T, reduction = "tsne", label = TRUE, label.size = 5, pt.size = 1.5) + NoLegend()

#Percent of cells and number of cells in cluster wise
dittoBarPlot(T, "Sample.ID", group.by = "T_subcluster")

dittoBarPlot(T, "Sample.ID", group.by = "T_subcluster", scale = "count")

saveRDS(T, file='figures/T_obj.rds')


# Macrophage cells sub clustering
MG<- subset(seurat_obj, idents = "Macrophage")
head(MG)

MG <- NormalizeData(MG, normalization.method = "LogNormalize", scale.factor = 10000)
MG <- FindVariableFeatures(MG, selection.method = "vst", nfeatures = 2000)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(MG)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2

MG <- ScaleData(MG, features = all.genes)
MG <- RunPCA(MG, features = VariableFeatures(object = MG))
MG <- RunUMAP(MG, reduction = "pca", dims = 1:30)
MG <- RunTSNE(MG, dims = 1:30)
MG <- FindNeighbors(MG, reduction = "pca", dims = 1:30)
#MG <- FindClusters(MG, resolution = 0.2)
MG <- FindSubCluster(MG, "Macrophage", "RNA_nn", subcluster.name = "MG_subcluster",  resolution = 0.2, algorithm = 1)
MG <- SetIdent(MG, value = MG@meta.data$MG_subcluster)

DimPlot(object = MG, reduction = "umap", split.by = "condition")
DimPlot(object = MG, reduction = "tsne", split.by = "condition", ncol = 4)


DotPlot(MG, features = c("Cd86", "Cd80", "Cd68", "Il1r", "Tlr2", "Socs3", "Nos1", "Nos3", "Ptgs2", 
                         "Cd163", "Cd206", "Cd200r", "Tgm2", "Fizz1", "Arg1"), dot.scale=9, cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()

MG <- RenameIdents(MG,
                  'Macrophage_0' = 'Monocyte',
                  'Macrophage_1' = 'Granulocyte1',
                  'Macrophage_2' = 'cDC',
                  'Macrophage_3' = 'M1',
                  'Macrophage_4' = 'M2-1',
                  'Macrophage_5' = 'M2-2',
                  'Macrophage_6' = 'pDC',
                  'Macrophage_7' = 'Granulocyte2'
)

DimPlot(MG, reduction = "umap", group.by = "condition", label.size = 5, pt.size = 1.5)
DimPlot(MG, reduction = "tsne", label = FALSE, label.size = 5, pt.size = 1.5)

MG$subcluster <- levels(MG)
saveRDS(MG, file='figures/MG_obj.rds')

#Percent of cells and number of cells in cluster wise
dittoBarPlot(MG, "Sample.ID", group.by = "MG_subcluster")
dittoBarPlot(MG, "Sample.ID", group.by = "MG_subcluster", scale = "count")



# Fibroblast cells sub clustering
FB<- subset(seurat_obj, idents = "Fibroblast")
head(FB)

FB <- NormalizeData(FB, normalization.method = "LogNormalize", scale.factor = 10000)
FB <- FindVariableFeatures(FB, selection.method = "vst", nfeatures = 2000)
plot1 <- VariableFeaturePlot(FB)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

FB <- ScaleData(FB, features = all.genes)
FB <- RunPCA(FB, features = VariableFeatures(object = FB))
FB <- RunUMAP(FB, reduction = "pca", dims = 1:30)
FB <- RunTSNE(FB, dims = 1:30)
FB <- FindNeighbors(FB, reduction = "pca", dims = 1:30)
#FB <- FindClusters(FB, resolution = 1.2)
FB <- FindSubCluster(FB, "Fibroblast", "RNA_nn", subcluster.name = "FB_subcluster",  resolution = 1.2, algorithm = 1)
FB <- SetIdent(FB, value = FB@meta.data$FB_subcluster)

DimPlot(object = FB, split.by = "condition")
DimPlot(object = FB, split.by = "condition", ncol = 4)


DotPlot(FB, features = c("Acta2", "Tagln", "Postn", "Tpm1", "Tpm2", "Pdgfra", "Apod", "Cxcl12", "Lmna", "Clec3b",
                         "Nusap1", "Cenpf", "Pttg1", "Stmn1"), dot.scale=9, cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + RotatedAxis()


FB <- RenameIdents(FB,
                   'Fibroblast_0' = 'qCAFs',
                   'Fibroblast_1' = 'apCAFs',
                   'Fibroblast_2' = 'pCAFs',
                   'Fibroblast_3' = "myCAFs",
                   'Fibroblast_4' = "iCAFs"
)

DimPlot(FB, reduction = "umap", label = FALSE, label.size = 5, pt.size = 5.5)
DimPlot(FB, reduction = "tsne", label = FALSE, label.size = 5, pt.size = 5.5)

saveRDS(FB, file='figures/FB_obj.rds')

#iCAFs markers
DotPlot(FB, features = c("Il6", "C3", "Cfd", "Clec3b", "Has1", "Igf1", "Cxcl12", "Ly6c1", "Itm2a"), dot.scale=9, cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + RotatedAxis()

#myCAFs markers
DotPlot(FB, features = c("Postn1", "Thy1", "Clo11a1", "Cthrc1", "Fap", "Tagln", "Thbs2", "Cxcl14", "Mmp11",
                          "Acta2", "Col1a1", "Col10a1"), dot.scale=9, cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + RotatedAxis()

#meCAFs markers
DotPlot(FB, features = c("Bnip3", "Slc2a1",
                         "Hilpda", "Vegfa", "Ndrg1"), dot.scale=9, cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + RotatedAxis()

#pCAFs markers
DotPlot(FB, features = c("Pdpn", "Tenpf", "Pttg1", "Stmn1", "Top2a", "Nusap1", "Cenpf"), dot.scale=9, cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + RotatedAxis()


#apCAFs markers
DotPlot(FB, features = c("Cd74", "H2-ab1",
                         "Saa3", "Slpi", "Nkam4", "Comp"), dot.scale=9, cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + RotatedAxis()

#panCAFs markers
DotPlot(FB, features = c("Pdpn", "Dcn", "Col1a1", "Fap", "Vim"), dot.scale=6, cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + RotatedAxis()

#Percent of cells and number of cells in cluster wise
dittoBarPlot(FB, "Sample.ID", group.by = "FB_subcluster")
dittoBarPlot(FB, "Sample.ID", group.by = "FB_subcluster", scale = "count")


#Differential markers within sub clusters
cluster_markers <- FindAllMarkers(
  FB,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 1.0,
  method='MAST'
)
head(cluster_markers)
#cluster_markers$subcluster <- paste0(unlist(cluster_markers$subcluster), '-', cluster_markers$subcluster)
write.csv(cluster_markers, file='subcluster_markers_FB.csv', quote=FALSE, row.names=FALSE)

# plot the number of DEGs per cluster:
df <- as.data.frame(rev(table(cluster_markers$cluster)))
head(df)
colnames(df) <- c('cluster', 'n_DEGs')

# bar plot of the number of cells in each sample
p <- ggplot(df, aes(y=n_DEGs, x=reorder(cluster, -n_DEGs), fill=cluster)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  NoLegend() + RotatedAxis() +
  ylab(expression(italic(N)[DEGs])) + xlab('') +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.y=element_line(colour="lightgray", size=0.5),
  )

pdf('FB_sub_cluster/DEGs_barplot_FB.pdf', width=4, height=3)
print(p)
dev.off()
cluster_markers <- read.csv("cluster_markers_FB.csv")
head(cluster_markers)
df <- as.data.frame(rev(table(cluster_markers$cluster)))

# plot the top 3 DEGs per cluster as a heatmap:
top_DEGs <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(10, wt=avg_log2FC) %>%
  .$gene

head(top_DEGs)

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
DoHeatmap(FB, features = top_DEGs, size = 5.1, angle = 35, label=TRUE)+ scale_fill_gradientn(colours = rev(mapal)) + theme(panel.spacing.x = unit(2.5, "cm")) + theme(text = element_text(size = 15))


#Autophagy genes in Fibroblast cells
genes <- c("Ambra1", "Ulk1", "Atg13", "Atg101",
           "Pik3r4" , "Becn1", "Atg14", "Nrbf2",
           "Atg16l1", "Atg12", "Atg5", "Gabarap", "Nbr1", "Sqstm1", 
           "Atg4c", "Atg4a", "Casp3", "Gabarapl1", "Ins1", "Ins2", 
           "Casp8", "Prkaa1", "Fadd", "Atg7",  "Ulk3", "Pik3cg", 
           "Gabarapl2", "Atg3", "Ptpn6", "Mcoln3", "Rxra", "Rxrb", "Pdia3", "Prkaa2")

DotPlot(FB, features=genes, dot.scale=8, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()

DotPlot(FB, features=genes, dot.scale=7, cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()

#subset of sub clusters
table(FB$FB_subcluster)
qCAFs<- subset(FB, subset=FB_subcluster == "Fibroblast_0")
apCAFs<- subset(FB, subset=FB_subcluster == "Fibroblast_1")
pCAFs<- subset(FB, subset=FB_subcluster == "Fibroblast_2")
myCAFs<- subset(FB, subset=FB_subcluster == "Fibroblast_3")
iCAFs<- subset(FB, subset=FB_subcluster == "Fibroblast_4")

DotPlot(qCAFs, features=genes, dot.scale=8, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()

DotPlot(apCAFs, features=genes, dot.scale=8, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()

DotPlot(pCAFs, features=genes, dot.scale=8, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()

DotPlot(myCAFs, features=genes, dot.scale=8, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()

DotPlot(iCAFs, features=genes, dot.scale=6, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()

#Percent of cells and number of cells in cluster wise
dittoBarPlot(FB, "Sample.ID", group.by = "subcluster")
dittoBarPlot(FB, "Sample.ID", group.by = "subcluster", scale = "count")

#ER stress marker genes
DotPlot(FB, features = c("Xbp1", "Ddit3", "Hspa5", "Atf4", "Hyou1", "Nfya", "Ire1", "Atf6", "Cebpa", "Ern1", "Mapk1", "Casp3"), dot.scale=9, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()

DotPlot(FB, features = c("Xbp1", "Ddit3", "Hspa5", "Atf4", "Hyou1", "Nfya", "Ire1", "Atf6", "Cebpa", "Ern1", "Mapk1", "Casp3"), dot.scale=9, group.by="FB_subcluster", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()




