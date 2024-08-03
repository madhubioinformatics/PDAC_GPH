

#Normalization
seurat_obj <- readRDS("filter_seurat_obj.rds")
head(seurat_obj)

VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by='Sample.ID',
  ncol = 1, pt.size=0)

# plot the number of cells in each sample post filtering
df <- as.data.frame(rev(table(seurat_obj$Sample.ID)))
colnames(df) <- c('SampleID', 'n_cells')
p <- ggplot(df, aes(y=n_cells, x=reorder(SampleID, -n_cells), fill=SampleID)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  NoLegend() + RotatedAxis() +
  ylab(expression(italic(N)[cells])) + xlab('Sample ID') +
  ggtitle(paste('Total cells post-filtering:', sum(df$n_cells))) +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.y=element_line(colour="lightgray", size=0.5),
  )

print(p)

# log normalize data
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# scale data:
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

#Feature Selection

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 4000
)

p <- LabelPoints(
  VariableFeaturePlot(seurat_obj),
  points = head(VariableFeatures(seurat_obj),10),
  repel = TRUE
) + theme(legend.position="bottom")

print(p)

#Linear Dimensionality Reduction

seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(object = seurat_obj),
  npcs=50
)

# plot the top genes contributing to the first 3 PCs
p <- VizDimLoadings(seurat_obj, dims = 1:3, reduction = "pca", ncol=3)
p

p <- DimPlot(seurat_obj, reduction = "pca", group.by='Sample.ID')
p

DimHeatmap(seurat_obj, dims = 1:16, cells = 500, balanced = TRUE, ncol=4)

ElbowPlot(seurat_obj, ndims = 50)

saveRDS(seurat_obj, file="filter_seurat_obj1.rds")
