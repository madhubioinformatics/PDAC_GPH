
#Autophagy genes in PDAC and Fibroblast (FB)
#Clustering and non-linear dimensionality reduction
seurat_obj <- readRDS("processed_seurat_object.rds")

#Subset of PDAC from Seurat object
PDAC <-subset(seurat_obj, subset = CellType == "PDAC")
head(PDAC)

#Subset of Fibroblasts(FB) from Seurat object
FB<- subset(seurat_obj, subset = CellType == "Fibroblast")
head(FB)

DimPlot(PDAC, reduction = "tsne", group.by = "CellType_cluster", label = 'T')
DimPlot(PDAC, reduction = "umap", group.by = "CellType_cluster", label = 'T')
DimPlot(seurat_obj, reduction = "umap", split.by = "condition", label = 'T')


#Autophagy genes
genes <- c("Ambra1", "Ulk1", "Rb1cc1", "Atg13", "Atg101", "Pik3c3",
          "Pik3r4" , "Becn1", "Atg14", "Uvrag", "Nrbf2",
           "Atg16l1", "Atg12", "Atg5", "Gabarap", "Nbr1", "Sqstm1", 
           "Atg4c", "Atg4a", "Casp3", "Gabarapl1", "Ins1", "Ins2", 
           "Casp8", "Prkaa1", "Fadd", "Atg7",  "Ulk3", "Pik3cg", 
          "Gabarapl2", "Atg3", "Ptpn6", "Mcoln3", "Rxra", "Rxrb", "Stat3", "Pdia3", "Prkaa2")


DotPlot(seurat_obj, features=genes, dot.scale=7, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + RotatedAxis()

DotPlot(PDAC, features=genes, dot.scale=7, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + RotatedAxis()



PDAC_8 <-subset(PDAC, subset=CellType_cluster == "PDAC-8")
head(PDAC_8)
table(PDAC_8$condition)

DotPlot(PDAC_8, features=genes, dot.scale=8, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()

PDAC_9 <-subset(PDAC, subset=CellType_cluster == "PDAC-9")
head(PDAC_9)
table(PDAC_9$condition)


DotPlot(PDAC_9, features=genes, dot.scale=7, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + RotatedAxis()


DotPlot(FB, features=genes, dot.scale=8, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + RotatedAxis()



#Immuno markers in PDAC and Fibroblast
Immuno <- c("Cd47", "Cd9", "Ctla2a", "Gprc5b", "Cd177", "Cd44",
           "Gpr160" , "Cd151", "Cd276", "Pdcd4", "Pdcd7",
           "Pdcd6", "Gpr180", "Cd74", "Cd38", "Cd81", "Cd68", 
           "Gpr39", "Gpr89a", "Alcam", "Pdcd10", "Mme", "Cd14", 
           "Gpr108", "Pdcd2", "Pdcd5", "Cd5", "Cd3e", "Cd4", "Cd8b", 
           "Cd82", "Cd99", "Cd36", "Cd109", "Cd180", "Cd48", "Cd34", 
           "Senp3", "Clip1", "Ctla4", "Pdcd1", "Cd274", "Foxp3", 
           "Il10", "Il6", "Il2", "Il12", "Il4", "Il13", "Ifn1", "Tigit")#Immune

DotPlot(seurat_obj, features=Immuno, dot.scale=6, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()


DotPlot(PDAC, features=Immuno, dot.scale=6, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + RotatedAxis()


DotPlot(FB, features=Immuno, dot.scale=6, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + RotatedAxis()


#ER stress marker genes in conditions
DotPlot(seurat_obj, features = c("Xbp1", "Ddit3", "Hspa5", "Atf4", "Hyou1", "Nfya", "Ire1", "Atf6", "Cebpa", "Ern1", "Mapk1", "Casp3"), dot.scale=9, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()

DotPlot(PDAC, features = c("Xbp1", "Ddit3", "Hspa5", "Atf4", "Hyou1", "Nfya", "Ire1", "Atf6", "Cebpa", "Ern1", "Mapk1", "Map1lc3a", "Map1lc3b"), dot.scale=9, group.by="condition", cols="RdBu") + theme(axis.text.x = element_text(angle = 90)) +
  coord_flip() + RotatedAxis()




