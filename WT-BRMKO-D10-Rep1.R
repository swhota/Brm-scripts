##this script was run using Seurat v3.1.4
#remotes::install_version("Seurat", version = "3.1.4")
library(Seurat)

#load data
brm.data <- Read10X(data.dir = "/path/to/directory/Cellranger_aggr/")
WT_BrmKO_D10_rep1 <- CreateSeuratObject(counts = brm.data, project = "Brm")

#add gemgroup and names as metadata columns
classification.vec <- as.numeric(gsub(".*-","", (colnames(WT_BrmKO_D10_rep1@assays$RNA@data))))
names(classification.vec) <- colnames(WT_BrmKO_D10_rep1@assays$RNA@data)
WT_BrmKO_D10_rep1 <- AddMetaData(WT_BrmKO_D10_rep1, classification.vec, "gem.group")
head(WT_BrmKO_D10_rep1@meta.data)
tail(WT_BrmKO_D10_rep1@meta.data)

#create a column with actual genotype names 
classification.vec1 <- as.numeric(gsub(".*-","", (colnames(WT_BrmKO_D10_rep1@assays$RNA@data))))
names(classification.vec1) <- colnames(WT_BrmKO_D10_rep1@assays$RNA@data)
# rename gem group assignments as something different 
classification.vec1[classification.vec1==1] <- "WT_BrmKO_D10_rep1-WT-1"
classification.vec1[classification.vec1==2] <- "WT_BrmKO_D10_rep1-KO-1"
WT_BrmKO_D10_rep1 <- AddMetaData(WT_BrmKO_D10_rep1, classification.vec1, "genotype.timepoint")
head(WT_BrmKO_D10_rep1@meta.data)

#QC
#calculate percent mitochondrial reads
WT_BrmKO_D10_rep1[["percent.mt"]] <- PercentageFeatureSet(WT_BrmKO_D10_rep1, pattern = "^mt-")
# Visualize QC metrics
VlnPlot(WT_BrmKO_D10_rep1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#apply cutoffs to remove dying cells and potential doublets
plot1 <- FeatureScatter(WT_BrmKO_D10_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WT_BrmKO_D10_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
WT_BrmKO_D10_rep1 <- subset(WT_BrmKO_D10_rep1, subset = nFeature_RNA > 250 & nFeature_RNA < 3000)

##Cell cycle scoring- load cellcycle genes
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
#Assign cell cycle scores
WT_BrmKO_D10_rep1 <- CellCycleScoring(object = WT_BrmKO_D10_rep1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(x = WT_BrmKO_D10_rep1@meta.data)
WT_BrmKO_D10_rep1@meta.data$CC.Difference <- WT_BrmKO_D10_rep1@meta.data$S.Score - WT_BrmKO_D10_rep1@meta.data$G2M.Score
WT_BrmKO_D10_rep1 <- NormalizeData(WT_BrmKO_D10_rep1)
WT_BrmKO_D10_rep1 <- FindVariableFeatures(WT_BrmKO_D10_rep1, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(WT_BrmKO_D10_rep1), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(WT_BrmKO_D10_rep1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
#scale data and regress variables
WT_BrmKO_D10_rep1 <- ScaleData(WT_BrmKO_D10_rep1, vars.to.regress = c("percent.mito", "nUMI", "CC.difference", "gem.group"))

#dimensionality reduction
WT_BrmKO_D10_rep1 <- RunPCA(WT_BrmKO_D10_rep1)
ElbowPlot(WT_BrmKO_D10_rep1, ndims = 30)
#cluster cells
WT_BrmKO_D10_rep1 <- FindNeighbors(WT_BrmKO_D10_rep1, dims = 1:15)
WT_BrmKO_D10_rep1 <- FindClusters(WT_BrmKO_D10_rep1, resolution = c(0.1,0.2,0.3,0.4,0.5))
WT_BrmKO_D10_rep1 <- FindClusters(WT_BrmKO_D10_rep1, resolution = c(1.0,1.5,2.0))
WT_BrmKO_D10_rep1 <- FindClusters(WT_BrmKO_D10_rep1, resolution = c(0.6))
#umap
WT_BrmKO_D10_rep1 <- RunUMAP(WT_BrmKO_D10_rep1, dims = 1:25, min.dist=0.1, n.neighbors = 5)
Idents(WT_BrmKO_D10_rep1) <- "RNA_snn_res.0.3"
DimPlot(WT_BrmKO_D10_rep1, reduction = "umap", label = T, label.size = 8) + NoLegend()
Idents(WT_BrmKO_D10_rep1) <- "genotype.timepoint"
DimPlot(WT_BrmKO_D10_rep1, reduction = "umap", label = T, label.size = 8, cols = c("red", "blue")) 

#save
save.image(file = "WT_BrmKO_D10_rep1.Rdata")

