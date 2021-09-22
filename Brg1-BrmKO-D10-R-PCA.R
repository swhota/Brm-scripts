##this script was run using Seurat 4.0.1
#install.packages('Seurat')
library(Seurat)

#load data
brgbrm.data <- Read10X(data.dir = "~/path/to/directory/Cellranger_aggr/")
brgbrm <- CreateSeuratObject(counts = brgbrm.data, project = "Brg-Brm")

#add gemgroup and names as metadata columns
classification.vec <- as.numeric(gsub(".*-","", (colnames(brgbrm@assays$RNA@data))))
names(classification.vec) <- colnames(brgbrm@assays$RNA@data)
brgbrm <- AddMetaData(brgbrm, classification.vec, "gem.group")
head(brgbrm@meta.data)

#create a column with actual genotype names 
classification.vec1 <- as.numeric(gsub(".*-","", (colnames(brgbrm@assays$RNA@data))))
names(classification.vec1) <- colnames(brgbrm@assays$RNA@data)
# rename gem group assignments as something different
classification.vec1[classification.vec1==1] <- "D10-Brg-WT"
classification.vec1[classification.vec1==2] <- "D10-Brg-KO"
classification.vec1[classification.vec1==3] <- "D10-Brm-WT-2"
classification.vec1[classification.vec1==4] <- "D10-Brm-KO-2"
classification.vec1[classification.vec1==5] <- "D10-Brm-WT-1"
classification.vec1[classification.vec1==6] <- "D10-Brm-KO-1"
brgbrm <- AddMetaData(brgbrm, classification.vec1, "genotype.timepoint")
head(brgbrm@meta.data)

#create a column with actual genotype names 
classification.vec2 <- as.numeric(gsub(".*-","", (colnames(brgbrm@assays$RNA@data))))
names(classification.vec2) <- colnames(brgbrm@assays$RNA@data)
# rename gem group assignments as something different
classification.vec2[classification.vec2==1] <- "D10-Brg-WT"
classification.vec2[classification.vec2==2] <- "D10-Brg-KO"
classification.vec2[classification.vec2==3] <- "D10-Brm-WT-2"
classification.vec2[classification.vec2==4] <- "D10-Brm-KO-2"
classification.vec2[classification.vec2==5] <- "D10-Brm-WT-1"
classification.vec2[classification.vec2==6] <- "D10-Brm-KO-1"
brgbrm <- AddMetaData(brgbrm, classification.vec2, "batch")
head(brgbrm@meta.data)

#create a column with actual genotype names 
classification.vec3 <- as.numeric(gsub(".*-","", (colnames(brgbrm@assays$RNA@data))))
names(classification.vec3) <- colnames(brgbrm@assays$RNA@data)
# rename gem group assignments as something different
classification.vec3[classification.vec3==1] <- "exp1"
classification.vec3[classification.vec3==2] <- "exp1"
classification.vec3[classification.vec3==3] <- "exp2"
classification.vec3[classification.vec3==4] <- "exp2"
classification.vec3[classification.vec3==5] <- "exp2"
classification.vec3[classification.vec3==6] <- "exp2"
brgbrm <- AddMetaData(brgbrm, classification.vec3, "experiment")
head(brgbrm@meta.data)

#QC
#calculate percent mitochondrial reads
brgbrm[["percent.mt"]] <- PercentageFeatureSet(brgbrm, pattern = "^mt-")
# Visualize QC metrics
VlnPlot(brgbrm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(brgbrm, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(brgbrm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#apply cutoffs to remove dying cells and potential doublets
brgbrm <- subset(brgbrm, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & nCount_RNA > 700 & nCount_RNA < 6000)

##run SCTransform and reciprocal PCA
list <- SplitObject(brgbrm, split.by = "experiment")
list <- lapply(X = list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = list, nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list, anchor.features = features)
list <- lapply(X = list, FUN = RunPCA, features = features)

anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT", 
                                  anchor.features = features, dims = 1:40, reduction = "rpca", k.anchor = 1)
brgbrmint <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:40)
brgbrmint <- RunPCA(brgbrmint, verbose = FALSE)
brgbrmint <- RunUMAP(brgbrmint, reduction = "pca", dims = 1:40)

#Visualize UMAP
DimPlot(brgbrmint, reduction = "umap", group.by = "genotype.timepoint", label = F, repel = TRUE, 
        cols = c("#006637","#62337F","#d94801","#2171b5"))

#save object
save.image(file = "Brg1_BrmKO_D10_RPCA.Rdata")
