library("URD")
library(Seurat)

#load seurat object and then convert to URD object
load("D4_D101st2nd_lowBMP_downsampled.Robj")
D4_D101st2nd_lowBMP_downsampled <- SetAllIdent(object = D4_D101st2nd_lowBMP_downsampled, id = "res.2.0")

D4_D101st2nd_lowBMP_downsampled_URD <- seuratToURD(D4_D101st2nd_lowBMP_downsampled)
vector <- WhichCells(D4_D101st2nd_lowBMP_downsampled)
rm(vector)

#calculate diffusion map
#using square root 131 of total cells 17205 and lower than auto sigma of 32
D4_D101st2nd_lowBMP_downsampled_URD <- calcDM(D4_D101st2nd_lowBMP_downsampled_URD, knn=131, sigma.use=26)
save(D4_D101st2nd_lowBMP_downsampled_URD, file = "D4_D101st2nd_lowBMP_downsampled_URD.Robj")
#plot diffusion maps using pairs of dimensions, for small datasets, structure of diff may already be apparent
plotDimArray(D4_D101st2nd_lowBMP_downsampled_URD, reduction.use = "dm", dims.to.plot = 1:8, outer.title = "Diffusion Map (Sigma 26, 131 NNs): Stage", label="niceorder", plot.title="", legend=T)


#define roots 
root.cells <- cellsInCluster(D4_D101st2nd_lowBMP_downsampled_URD, "res.2", cluster=c(2,16))

# Then we run 'flood' simulations
combined.floods <- floodPseudotime(D4_D101st2nd_lowBMP_downsampled_URD, root.cells = root.cells, n=80, minimum.cells.flooded = 10, verbose=T)

# The we process the simulations into a pseudotime
D4_D101st2nd_lowBMP_downsampled_URD <- floodPseudotimeProcess(D4_D101st2nd_lowBMP_downsampled_URD, combined.floods, floods.name="pseudotime")

#We can make sure that enough simulations have been performed by looking at the change in cell pseudotime as more simulations are added. 
#Here, we can see that an asymptote was reached around 80 simulations, so 80 was just enough.
pseudotimePlotStabilityOverall(D4_D101st2nd_lowBMP_downsampled_URD)

#We can also plot pseudotime on the tSNE (to confirm that it makes sense)
plotDim(D4_D101st2nd_lowBMP_downsampled_URD,reduction.use = "tsne", "pseudotime")


#now plot distribution of pseudotime for each timepoint. 
#are timepoints in the correct order? is there overlap between neighboring timepoints (as expected), 
#but they do not completely collapse on top of each other (which often indicates that sigma is too large in the diffusion map).
plotDists(D4_D101st2nd_lowBMP_downsampled_URD, "pseudotime", "niceorder", plot.title="Pseudotime by timepoint")


save(D4_D101st2nd_lowBMP_downsampled_URD, file = "D4_D101st2nd_lowBMP_downsampled_URD.Robj")

#define tip cells; here we chose all WT clusters from day 23
D4_D101st2nd_lowBMP_downsampled <- SetAllIdent(object = D4_D101st2nd_lowBMP_downsampled, id = "res.2")
tipcells <- SubsetData(D4_D101st2nd_lowBMP_downsampled, ident.use = c(6,17,8,0,24,10,13,19,5,15,7,11,18,21,25,28))
#convert tipcells object to URD
tipcellsURD <- seuratToURD(tipcells)
# Copy cluster identities from tipcells URD object to a new clustering ("tip.clusters") in the full URD object.
D4_D101st2nd_lowBMP_downsampled_URD@group.ids[rownames(tipcellsURD@group.ids), "tip.clusters"] <- tipcellsURD@group.ids$res.2

# Determine the parameters of the logistic used to bias the transition probabilities. The procedure
# is relatively robust to this parameter, but the cell numbers may need to be modified for larger
# or smaller data sets.
D4_D101st2nd_lowBMP_downsampled_URD.ptlogistic <- pseudotimeDetermineLogistic(D4_D101st2nd_lowBMP_downsampled_URD, "pseudotime", optimal.cells.forward=50, max.cells.back=100, do.plot = T)

save(combined.floods, file="D4_D101st2nd_lowBMP_combinedfloods.Robj")

# Bias the transition matrix acording to pseudotime
D4_D101st2nd_lowBMP_downsampled_URD.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(D4_D101st2nd_lowBMP_downsampled_URD, "pseudotime", logistic.params=D4_D101st2nd_lowBMP_downsampled_URD.ptlogistic))

# Simulate the biased random walks from each tip
D4_D101st2nd_lowBMP_downsampled_URD.walks <- simulateRandomWalksFromTips(D4_D101st2nd_lowBMP_downsampled_URD, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = D4_D101st2nd_lowBMP_downsampled_URD.biased.tm, 
                                                                        n.per.tip = 25000, root.visits = 1, max.steps = 2500, verbose = F)

#3906 cells not visited
# Process the biased random walks into visitation frequencies
D4_D101st2nd_lowBMP_downsampled_URD <- processRandomWalksFromTips(D4_D101st2nd_lowBMP_downsampled_URD, D4_D101st2nd_lowBMP_downsampled_URD.walks, verbose = F)

# Load the cells used for each tip into the URD object
lowBMP.tree <- loadTipCells(D4_D101st2nd_lowBMP_downsampled_URD, "tip.clusters")

# Build tree1
lowBMP.tree <- buildTree(lowBMP.tree, pseudotime = "pseudotime", tips.use = c(6,17,8,0,24,10,13,19,5,15,7,11,18,21,25,28), divergence.method = "preference", 
                           cells.per.pseudotime.bin = 80, bins.per.pseudotime.window = 8, min.cells.per.segment = 10,
                           p.thresh=0.05)
#plotDim(D4_D101st2nd_lowBMP_downsampled_URD, label="res.0.7", na.rm=F)

save(lowBMP.tree, file = "lowBMP_D4_D6_D101st2nd.Robj")
save(D4_D101st2nd_lowBMP_downsampled_URD, file = "D4_D101st2nd_lowBMP_downsampled_URD.Robj")

# Name the segments based on our previous determination of the identity of tips 1 and 2.
#combined.tree <- nameSegments(combined.tree, segments=c("14","1"), segment.names = c("Endoderm", "Cardiomyocytes"), short.names = c("Endo", "CMs"))

#plot 2D tree
#plot only wt
plotTree(lowBMP.tree, "niceorder", title="Day4_Day6_Day10_lowBMP", cell.alpha = 1.0,cell.size = 0.6,tree.alpha = 0.3,
         discrete.colors = c("#d94801", "#2171b5", "#d94801", "#2171b5", "#fdd0a2", "#c6dbef",  "#fd8d3c", "#6baed6"))


#plot all timepoints and genotypes
plotTree(lowBMP.tree, label.segments = T)


save(D4_D101st2nd_lowBMP_downsampled_URD, file="D4_D101st2nd_lowBMP_downsampled_URD.Robj")
save(lowBMP.tree, file="lowBMP_tree.Robj")

#to draw feature plots of gene expression on tree
genestoplot <- c("Foxp2")
for (gene in genestoplot) {
  plot(plotTree(lowBMP.tree, gene))
}

#get cells sorted by niceorder in each segment
lowBMP1_cellnumbers <- table(lowBMP.tree@group.ids$niceorder, lowBMP.tree@group.ids$segment)
write.csv(lowBMP1_cellnumbers, file = "lowBMP1_cellnumbers.csv")


#code from Daniel Carlin, to extract cell barcodes from a particular segment 
segment18<-list()
segment18[['intermediates']]<-c(combined.tree@tree$`cells.in.segment`$`18`)
#after this switch to differential gene test in Seurat


#set identity of the extracted cells in your seurat object
segment43vec <- segment43[["seg43"]]
D4_D101st2nd_lowBMP_downsampled <- SetIdent(D4_D101st2nd_lowBMP_downsampled, cells.use = segment43vec, ident.use = "segment43")
segment43_markers <- FindMarkers(D4_D101st2nd_lowBMP_downsampled, ident.1 = "segment43", max.cells.per.ident = 500)
