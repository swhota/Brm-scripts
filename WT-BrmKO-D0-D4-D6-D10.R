##this script was run using Seurat 3.1.4
#remotes::install_version("Seurat", version = "3.1.4")
library(Seurat)

#load data
brm.data <- Read10X(data.dir = "~/path/to/directory/Cellranger_aggr")
brmd0d4d6d10 <- CreateSeuratObject(counts = brm.data, project = "Brm")
rm(brm.data)

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
#apply cutoffs to remove dying cells and potential doublets
plot1 <- FeatureScatter(brmd0d4d6d10, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(brmd0d4d6d10, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
brmd0d4d6d10 <- subset(brmd0d4d6d10, subset = nFeature_RNA > 750 & nFeature_RNA < 4000 & nCount_RNA > 2500 & nCount_RNA < 20000)

#cellcycle scoring
#A list of cell cycle markers, from Tirosh et al, 2015 is read and segregated into markers of G2/M phase and markers of S phase
cc.genes <- readLines(con = "~/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
#Assign cell cycle scores
brmd0d4d6d10 <- CellCycleScoring(object = brmd0d4d6d10, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
head(x = brmd0d4d6d10@meta.data)
brmd0d4d6d10@meta.data$CC.Difference <- brmd0d4d6d10@meta.data$S.Score - brmd0d4d6d10@meta.data$G2M.Score

#Run SCTransform
brmd0d4d6d10 <- SCTransform(brmd0d4d6d10, vars.to.regress = c("percent.mt", "CC.Difference"), verbose = FALSE)

#dimensionality reduction
brmd0d4d6d10 <- RunPCA(brmd0d4d6d10)
ElbowPlot(brmd0d4d6d10, ndims = 50)
#cluster cells
brmd0d4d6d10 <- FindNeighbors(brmd0d4d6d10, dims = 1:40)
brmd0d4d6d10 <- FindClusters(brmd0d4d6d10, resolution = c(0.1,0.2,0.3,0.4,0.5))
brmd0d4d6d10 <- FindClusters(brmd0d4d6d10, resolution = c(1.0,1.5,2.0))

#umap
brmd0d4d6d10 <- RunUMAP(brmd0d4d6d10, dims = 1:40)
Idents(brmd0d4d6d10) <- "SCT_snn_res.1"
DimPlot(brmd0d4d6d10, reduction = "umap", label = T, label.size = 8) + NoLegend()
DimPlot(object = brmd0d4d6d10, reduction = "umap", group.by = "niceorder",
        cols = c("#ffe6da","#e5eef8", "#d94801","#2171b5","#fdd0a2","#c6dbef","#fd8d3c","#6baed6")) + NoLegend()

#feature plots
FeaturePlot(brmd0d4d6d10, features = c("Mesp1"), cols = c("grey","purple","navy","darkblue"))
FeaturePlot(brmd0d4d6d10, features = c("Smarcd3"), cols = c("grey","purple","navy","darkblue"))
FeaturePlot(brmd0d4d6d10, features = c("Tnnt2"), cols = c("grey","purple","navy","darkblue"))
FeaturePlot(brmd0d4d6d10, features = c("Sox2"), cols = c("grey","purple","navy","darkblue"))
FeaturePlot(brmd0d4d6d10, features = c("Gbx2"), cols = c("grey","purple","navy","darkblue"))
FeaturePlot(brmd0d4d6d10, features = c("Crabp2"), cols = c("grey","purple","navy","darkblue"))
FeaturePlot(brmd0d4d6d10, features = c("Wt1"))
FeaturePlot(brmd0d4d6d10, features = c("Kdr"))
FeaturePlot(brmd0d4d6d10, features = c("Hbb-b2"))

##find all markers for cell types
Idents(brmd0d4d6d10) <- "SCT_snn_res.1"
brmd0d4d6d10_markers <- FindAllMarkers(brmd0d4d6d10)
write.csv(brmd0d4d6d10_markers, file = "brmd0d4d6d10_markers_D101st2nd_res1.csv")

#save object
save.image(file = "WT_BrmKO_D0_D4_D6_D10.Rdata")


