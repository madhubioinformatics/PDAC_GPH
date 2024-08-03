
#scRNA-Seq data analysis and visualization
install.packages("BiocManager")
BiocManager::install("Seurat")
BiocManager::install("scCustomize")
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
BiocManager::install("Matrix", force = TRUE)
BiocManager::install("viridis")

library(Seurat)
library(tidyverse)
library(Matrix)
library(cowplot)
theme_set(theme_cowplot())
library(scCustomize)
library(ggplot2)
library(viridis)
library(patchwork)
library(DoubletFinder)

#Load CellRange ouput data
cellranger_out_dir <- "/dfs7/swaruplab/msaddala/test/ganji/raw/"
set_dir <- cellranger_out_dir
# create seurat object:
samples.Set <- dir(cellranger_out_dir) # getting all name in the directory in the dataset

seurat_listSet <- lapply(samples.Set, function(x){
  print(x)
  data <- Read10X(paste0(set_dir,x,''))
  
  cur_seurat <- CreateSeuratObject(
    counts =  data,
    project='Ganji',
    min.cells=3, min.features=200
  )
  cur_seurat$Sample.ID <- x
  return(cur_seurat)
})

# merge seurat object
seurat_obj <- merge(x=seurat_listSet[[1]], y=seurat_listSet[2:length(seurat_listSet)])
head(seurat_obj)


#Add metadata to the seurat object

meta <- read.csv("/dfs7/swaruplab/msaddala/test/ganji/metadata_ganji.csv")
head(meta)
rownames(meta) <- meta$Sample.ID
metadata <- meta[seurat_obj$Sample.ID,]
head(metadata)

for(meta in names(metadata)){
  seurat_obj@meta.data[[meta]] <- metadata[match(as.character(seurat_obj@meta.data$Sample.ID), metadata$Sample.ID), meta]
}

saveRDS(seurat_obj, file="/dfs7/swaruplab/msaddala/test/ganji/meta_seurat_obj.rds")

# use table function to get the number of cells in each Sample as a dataframe
df <- as.data.frame(rev(table(seurat_obj$Sample.ID)))
colnames(df) <- c('SampleID', 'n_cells')

# bar plot of the number of cells in each sample
p <- ggplot(df, aes(y=n_cells, x=reorder(SampleID, -n_cells), fill=SampleID)) +
  geom_bar(stat='identity') +
  scale_y_continuous(expand = c(0,0)) +
  NoLegend() + RotatedAxis() +
  ylab(expression(italic(N)[cells])) + xlab('Sample ID') +
  ggtitle(paste('Total cells:', sum(df$n_cells))) +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.y=element_line(colour="lightgray", size=0.5),
  )
print(p)


# calculate the percentage of mitochondrial reads per cell
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-") #mt for mouse and rat and MT for humans only.

# plot distributions of QC metrics, grouped by SampleID
png('figures/basic_qc.png', width=10, height=12, res=200, units='in')
VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by='Sample.ID',
  ncol = 1, pt.size=0)
dev.off()

# apply filter
seurat_obj <- subset(seurat_obj, nCount_RNA >= 150 & nCount_RNA <= 30000 & percent.mt <= 10)

# split aggregated data by sample
table(seurat_obj$Sample.ID)

# plot distributions of QC metrics, grouped by SampleID after filtered
png('figures/basic_qc_filtered.png', width=10, height=12, res=200, units='in')
VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by='Sample.ID',
  ncol = 1, pt.size=0)
dev.off()

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

saveRDS(seurat_obj, file="/dfs7/swaruplab/msaddala/test/ganji/filter_seurat_obj.rds")
