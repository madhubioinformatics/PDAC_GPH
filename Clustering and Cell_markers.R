#Clustering and non-linear dimensionality reduction
seurat_obj <- readRDS("filter_seurat_obj1.rds")
# KNN and clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)
head(seurat_obj)

# non-linear reductions (UMAP & t-SNE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
seurat_obj <- RunTSNE(seurat_obj, dims = 1:30)


DimPlot(seurat_obj, reduction = "umap", group.by='seurat_clusters', label.size = 5, label=TRUE) +
  umap_theme + NoLegend() + ggtitle('UMAP colored by seurat clusters')
DimPlot(seurat_obj, reduction = "umap", group.by='seurat_clusters', label.size = 4, label=TRUE, split.by = "Sample.ID") +
  umap_theme + NoLegend() + ggtitle('UMAP colored by seurat clusters')

DimPlot(seurat_obj, reduction = "tsne", group.by='seurat_clusters', label.size = 4, label=TRUE) +
  umap_theme + NoLegend() + ggtitle('t-SNE colored by seurat clusters')

plot_grid(ncol = 1, DimPlot(seurat_obj, reduction = "umap", group.by = "Sample.ID"))
plot_grid(ncol = 1, DimPlot(seurat_obj, reduction = "tsne", group.by = "Sample.ID"))

saveRDS(seurat_obj, file="normalized_seurat_obj.rds")


#Crude cell-type annotation
seurat_obj <- readRDS("normalized_seurat_obj.rds")

# set up list of canonical cell type markers
T_cell <- c('Cd3d', 'Cd3e', 'Cd3g')
NK_cell <- c('Nkg7', 'Il2rb')
B_cell <- c('Cd79b', 'Ms4a1', 'Mzb1')
Macrophage <- c('Cd14', 'Ms4a7', 'Cd68', 'C1qc')
Endothelium <- c('Cdh5', 'Cldn5', 'Pecam1')
PDAC <- c('Krt18', 'Mmp7', 'Krt19')#PDAC
Fibroblasts <- c('Col1a1', 'Lum', 'Pdgfra')


features <- list("T cell" = T_cell, "NK" = NK_cell, "B cell" = B_cell, "Macro" = Macrophage, "Endo" = Endothelium, "Fibro" = Fibroblasts, "PDAC" = PDAC)
head(features)

#Dotplot
DotPlot(object = seurat_obj, features=features, dot.scale=7, cols="RdBu", cluster.idents=T) + theme(axis.text.x = element_text(angle = 90))

#Vlnplot
VlnPlot(seurat_obj, features = T_cell, ncol = 1, group.by = "seurat_clusters", pt.size =0)
VlnPlot(seurat_obj, features = NK_cell, ncol = 1, group.by = "seurat_clusters", pt.size =0)
VlnPlot(seurat_obj, features = B_cell, ncol = 1, group.by = "seurat_clusters", pt.size =0)
VlnPlot(seurat_obj, features = Macrophage, ncol = 1, group.by = "seurat_clusters", pt.size =0)
VlnPlot(seurat_obj, features = Endothelium, ncol = 1, group.by = "seurat_clusters", pt.size =0)
VlnPlot(seurat_obj, features = Fibroblasts, ncol = 1, group.by = "seurat_clusters", pt.size =0)
VlnPlot(seurat_obj, features = PDAC, ncol = 1, group.by = "seurat_clusters", pt.size =0)


#Stack vlnplot
features <- c('Nkg7', 'Cd3d', 'Mzb1', 'Cd19', 'Cd68', 'C1qc', 'Cdh5', 'Pecam1', 'Lum', 'Pdgfra', 'S100a8', 'S100a9', 'Krt18', 'Krt19', 'Rgs5', 'Prom1')
a <- VlnPlot(seurat_obj, features, stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Identity on x-axis")
a


#makers finding between clusters
head(seurat_obj)
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
head(markers)
markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

write.csv(markers, file='Cluster_markers.csv', quote=FALSE, row.names=FALSE)

markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) -> top3
DoHeatmap(seurat_obj, features = top3$gene) + NoLegend() + theme(text = element_text(size = 0))


#Features of data after filters
VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by='seurat_clusters',
  ncol = 1, pt.size=0)


cluster_annotations <- list(
  '0' = 'PDAC',
  '1' = 'PDAC',
  '2' = 'PDAC',
  '3' = 'PDAC',
  '4' = 'PDAC',
  '5' = 'Macrophage',
  '6' = 'T_cell',
  '7' = 'PDAC',
  '8' = 'PDAC',
  '9' = 'PDAC',
  '10' = 'NK_cell',
  '11' = 'T_cell',
  '12' = 'T_cell',
  '13' = 'Macrophage',
  '14' = 'PDAC',
  '15' = 'Macrophage',
  '16' = 'Macrophage',
  '17' = 'PDAC',
  '18' = 'T_cell',
  '19' = 'PDAC',
  '20' = 'NK_cell',
  '21' = 'T_cell',
  '22' = 'B_cell',
  '23' = 'Macrophage',
  '24' = 'Fibroblast',
  '25' = 'Endothelial',
  '26' = 'PDAC'
)


DimPlot(seurat_obj, label = TRUE)
# add CellType to seurat metadata
seurat_obj$CellType <- unlist(cluster_annotations[seurat_obj$seurat_clusters])
seurat_obj$CellType_cluster <- paste0(seurat_obj$CellType, '-', seurat_obj$seurat_clusters)


DimPlot(seurat_obj, reduction = "umap", group.by='CellType') +
  umap_theme + ggtitle('UMAP colored by cell type annotations')


DimPlot(seurat_obj, reduction = "umap", group.by='CellType_cluster', label=FALSE) +
  umap_theme + ggtitle('UMAP colored by cell type + cluster') + NoLegend()


#TSNE MAP visualization
DimPlot(seurat_obj, reduction = "tsne", group.by='CellType', label=FALSE) +
  umap_theme + ggtitle('TSNE colored by cell type annotations')


DimPlot(seurat_obj, reduction = "tsne", group.by='CellType_cluster', label=TRUE) +
  umap_theme + ggtitle('TSNE colored by cell type + cluster') + NoLegend()


DimPlot(seurat_obj, reduction = "tsne", split.by='Sample.ID', label=FALSE) +
  umap_theme + ggtitle('TSNE colored by cell type annotations')


DimPlot(seurat_obj, reduction = "tsne", group.by='CellType', label=TRUE) +
  umap_theme + ggtitle('TSNE colored by cell type + cluster') + NoLegend()

#Saving and loading Seurat objects

# add barcode and UMAP to metadata
seurat_obj$barcode <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]

# save seurat object
saveRDS(seurat_obj, file='processed_seurat_object.rds')

# Plot ideal resolution on tSNE and UMAP
plotTheme <- theme_classic(base_size = 18)
nClust <- uniqueN(Idents(seurat_obj))         # Setup color palette
colCls <- colorRampPalette(brewer.pal(n = 10, name = "Paired"))(nClust)

p1 <- DimPlot(seurat_obj, reduction = "tsne", pt.size = 0.1, label = FALSE, 
              label.size = 3, cols = colCls) + plotTheme + coord_fixed()
p2 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.1, label = FALSE, 
              label.size = 3, cols = colCls) + plotTheme + coord_fixed()
ggsave(p1 + p2 + plot_layout(guides = "collect"), 
       width = 10, height = 4, filename = "clustDimPlot.png")

#Split Dimplot by samples
DimPlot(seurat_obj, reduction = "umap", split.by = "Sample.ID", pt.size = 0.1, label = FALSE, 
              label.size = 3, cols = colCls) + plotTheme + coord_fixed()

colGEX = c("grey85", brewer.pal(7, "Reds"))      # color for gene expression
# Proportion / cell number composition per library
ggData <- as.matrix(prop.table(table(seurat_obj$CellType, seurat_obj$Sample.ID), margin = 2))
pheatmap(ggData, color = colorRampPalette(colGEX)(100), 
         cluster_rows = FALSE, cutree_cols = 2, 
         display_numbers = TRUE, number_format = "%.3f", angle_col = 315,
         width = 4, height = 6, filename = "clustComLibH.png")

ggData = data.frame(prop.table(table(seurat_obj$CellType, seurat_obj$Sample.ID), margin = 2))
colnames(ggData) = c("CellType", "Sample.ID", "value")
p1 <- ggplot(ggData, aes(Sample.ID, value, fill = CellType)) +
  geom_col() + xlab("Sample.ID") + ylab("Proportion of Cells (%)") +
  scale_fill_manual(values = colCls) + plotTheme + coord_flip() + theme(axis.title = element_text(size=8), axis.text.y = element_text(size=8),
                                                                        axis.text.x = element_text(size=8))

ggData = data.frame(table(seurat_obj$CellType, seurat_obj$Sample.ID))
colnames(ggData) = c("CellType", "Sample.ID", "value")
p2 <- ggplot(ggData, aes(Sample.ID, value, fill = CellType)) +
  geom_col() + xlab("Sample.ID") + ylab("Cell Number") +
  scale_fill_manual(values = colCls) + plotTheme + coord_flip() + theme(axis.title = element_text(size=8), axis.text.y = element_text(size=8),
                                                                        axis.text.x = element_text(size=8))

ggsave(p1 + p2 + plot_layout(guides = "collect"), 
       width = 8, height = 3, filename = "clustComLib.png")


colLib = brewer.pal(7, "Paired")                 # color for libraries
# Proportion / cell number composition per CellType
ggData = data.frame(prop.table(table(seurat_obj$Sample.ID, seurat_obj$CellType), margin = 2))
colnames(ggData) = c("Sample.ID", "CellType", "value")
p1 <- ggplot(ggData, aes(CellType, value, fill = Sample.ID)) +
  geom_col() + xlab("CellType") + ylab("Proportion of Cells (%)") +
  scale_fill_manual(values = colLib) + plotTheme + coord_flip() + theme(axis.title = element_text(size=8), axis.text.y = element_text(size=8),
                                                                        axis.text.x = element_text(size=8))

ggData = data.frame(table(seurat_obj$Sample.ID, seurat_obj$CellType))
colnames(ggData) = c("Sample.ID", "CellType", "value")
p2 <- ggplot(ggData, aes(CellType, value, fill = Sample.ID)) +
  geom_col() + xlab("CellType") + ylab("Cell Number") +
  scale_fill_manual(values = colLib) + plotTheme + coord_flip() + theme(axis.title = element_text(size=8), axis.text.y = element_text(size=8),
                                                                        axis.text.x = element_text(size=8))
ggsave(p1 + p2 + plot_layout(guides = "collect"), 
       width = 7, height = 2, filename = "clustComClust.png")


# save seurat object
saveRDS(seurat_obj, file='processed_seurat_object.rds')

