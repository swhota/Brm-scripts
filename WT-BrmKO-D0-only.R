##this script was run using Seurat 3.1.4
#remotes::install_version("Seurat", version = "3.1.4")
library(Seurat)

#load data
brm.data <- Read10X(data.dir = "/path/to/directory/Cellranger_aggr/")
brmd0d4d6d10 <- CreateSeuratObject(counts = brm.data, project = "Brm")

#add gemgroup and names as metadata columns
classification.vec <- as.numeric(gsub(".*-","", (colnames(brmd0d4d6d10@assays$RNA@data))))
names(classification.vec) <- colnames(brmd0d4d6d10@assays$RNA@data)
brmd0d4d6d10 <- AddMetaData(brmd0d4d6d10, classification.vec, "gem.group")
head(brmd0d4d6d10@meta.data)

#create a column with actual genotype names 
classification.vec1 <- as.numeric(gsub(".*-","", (colnames(brmd0d4d6d10@assays$RNA@data))))
names(classification.vec1) <- colnames(brmd0d4d6d10@assays$RNA@data)
# rename gem group assignments as something different 
classification.vec1[classification.vec1==1] <- "D0-WT-1"
classification.vec1[classification.vec1==2] <- "D0-WT-2"
classification.vec1[classification.vec1==3] <- "D0-BrmKO-1"
classification.vec1[classification.vec1==4] <- "D0-BrmKO-2"
classification.vec1[classification.vec1==5] <- "D4-WT"
classification.vec1[classification.vec1==6] <- "D4-BrmKO"
classification.vec1[classification.vec1==7] <- "D6-WT"
classification.vec1[classification.vec1==8] <- "D6-BrmKO"
classification.vec1[classification.vec1==9] <- "D10-WT-2"
classification.vec1[classification.vec1==10] <- "D10-KO-2"
classification.vec1[classification.vec1==11] <- "D10-WT-1"
classification.vec1[classification.vec1==12] <- "D10-KO-1"
brmd0d4d6d10 <- AddMetaData(brmd0d4d6d10, classification.vec1, "genotype.timepoint")
head(brmd0d4d6d10@meta.data)

#create a column with actual genotype names 
classification.vec2 <- as.numeric(gsub(".*-","", (colnames(brmd0d4d6d10@assays$RNA@data))))
names(classification.vec2) <- colnames(brmd0d4d6d10@assays$RNA@data)
# rename gem group assignments as something different 
classification.vec2[classification.vec2==1] <- "D0"
classification.vec2[classification.vec2==2] <- "D0"
classification.vec2[classification.vec2==3] <- "D0"
classification.vec2[classification.vec2==4] <- "D0"
classification.vec2[classification.vec2==5] <- "D4"
classification.vec2[classification.vec2==6] <- "D4"
classification.vec2[classification.vec2==7] <- "D6"
classification.vec2[classification.vec2==8] <- "D6"
classification.vec2[classification.vec2==9] <- "D10"
classification.vec2[classification.vec2==10] <- "D10"
classification.vec2[classification.vec2==11] <- "D10"
classification.vec2[classification.vec2==12] <- "D10"
brmd0d4d6d10 <- AddMetaData(brmd0d4d6d10, classification.vec2, "timepoint")
head(brmd0d4d6d10@meta.data)

#create a column with actual genotype names 
classification.vec3 <- as.numeric(gsub(".*-","", (colnames(brmd0d4d6d10@assays$RNA@data))))
names(classification.vec3) <- colnames(brmd0d4d6d10@assays$RNA@data)
# rename gem group assignments as something different 
classification.vec3[classification.vec3==1] <- "D0-WT"
classification.vec3[classification.vec3==2] <- "D0-WT"
classification.vec3[classification.vec3==3] <- "D0-KO"
classification.vec3[classification.vec3==4] <- "D0-KO"
classification.vec3[classification.vec3==5] <- "D4-WT"
classification.vec3[classification.vec3==6] <- "D4-KO"
classification.vec3[classification.vec3==7] <- "D6-WT"
classification.vec3[classification.vec3==8] <- "D6-KO"
classification.vec3[classification.vec3==9] <- "D10-WT"
classification.vec3[classification.vec3==10] <- "D10-KO"
classification.vec3[classification.vec3==11] <- "D10-WT"
classification.vec3[classification.vec3==12] <- "D10-KO"
brmd0d4d6d10 <- AddMetaData(brmd0d4d6d10, classification.vec3, "niceorder")
head(brmd0d4d6d10@meta.data)

#QC
#calculate percent mitochondrial reads
brmd0d4d6d10[["percent.mt"]] <- PercentageFeatureSet(brmd0d4d6d10, pattern = "^mt-")
# Visualize QC metrics
VlnPlot(brmd0d4d6d10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(brmd0d4d6d10, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(brmd0d4d6d10, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
#apply cutoffs to remove dying cells and potential doublets
brmd0d4d6d10 <- subset(brmd0d4d6d10, subset = nFeature_RNA > 750 & nFeature_RNA < 4000 & nCount_RNA > 2500 & nCount_RNA < 20000)

##subset day 0 datasets
Idents(brmd0d4d6d10) <- "timepoint"
d0only <- subset(brmd0d4d6d10, idents = "D0")

##run SCTransform
d0only <- SCTransform(d0only, vars.to.regress = c("percent.mt"), verbose = FALSE)

#dimensionality reduction
d0only <- RunPCA(d0only)
ElbowPlot(d0only, ndims = 50)
#cluster cells
d0only <- FindNeighbors(d0only, dims = 1:40)
d0only <- FindClusters(d0only, resolution = c(0.1,0.2,0.3,0.4,0.5))
#umap
d0only <- RunUMAP(d0only, dims = 1:40)
Idents(d0only) <- "SCT_snn_res.0.1"
DimPlot(d0only, reduction = "umap", label = T, label.size = 8) + NoLegend()
DimPlot(object = d0only, reduction = "umap", group.by = "niceorder", cols = c("plum","skyblue"), pt.size = 1.0)

#save object
save.image(file = "WT_BrmKO_D0.Rdata")

##differential gene test
Idents(d0only) <- "niceorder"
D0_WTvsKO <- FindMarkers(d0only, ident.1 = "D0-WT", ident.2 = "D0-KO", logfc.threshold = 0.1)
write.csv(D0_WTvsKO, file = "D0_WTvsKO_differentialgenes.csv")




