library("URD")
library(Seurat)

#load seurat object and then convert to URD object
load("D4_D101st2nd_downsampled.Robj")
D4_D101st2nd_downsampled <- SetAllIdent(object = D4_D101st2nd_downsampled, id = "res.2")

D4_D101st2nd_highlow_downsampled_URD <- seuratToURD(D4_D101st2nd_downsampled)
vector <- WhichCells(D4_D101st2nd_downsampled)
rm(vector)

#calculate diffusion map
#using square root 157 of 24643 total cells  and lower than auto sigma of 28.7
D4_D101st2nd_highlow_downsampled_URD <- calcDM(D4_D101st2nd_highlow_downsampled_URD, knn=157, sigma.use=22)
save(D4_D101st2nd_highlow_downsampled_URD, file = "D4_D101st2nd_highlow_downsampled_URD.Robj")
#plot diffusion maps using pairs of dimensions, for small datasets, structure of diff may already be apparent
plotDimArray(D4_D101st2nd_highlow_downsampled_URD, reduction.use = "dm", dims.to.plot = 1:8, outer.title = "Diffusion Map (Sigma 22, 157 NNs): Stage", label="niceorder", plot.title="", legend=T)


#define roots 
root.cells <- cellsInCluster(D4_D101st2nd_highlow_downsampled_URD, "res.2", cluster=3)

# Then we run 'flood' simulations
combined.floods <- floodPseudotime(D4_D101st2nd_highlow_downsampled_URD, root.cells = root.cells, n=80, minimum.cells.flooded = 10, verbose=T)

# The we process the simulations into a pseudotime
D4_D101st2nd_highlow_downsampled_URD <- floodPseudotimeProcess(D4_D101st2nd_highlow_downsampled_URD, combined.floods, floods.name="pseudotime")

#We can make sure that enough simulations have been performed by looking at the change in cell pseudotime as more simulations are added. 
#Here, we can see that an asymptote was reached around 80 simulations, so 80 was just enough.
pseudotimePlotStabilityOverall(D4_D101st2nd_highlow_downsampled_URD)

#We can also plot pseudotime on the tSNE (to confirm that it makes sense)
plotDim(D4_D101st2nd_highlow_downsampled_URD,reduction.use = "tsne", "pseudotime")


#now plot distribution of pseudotime for each timepoint.
plotDists(D4_D101st2nd_highlow_downsampled_URD, "pseudotime", "niceorder", plot.title="Pseudotime by timepoint")


save(D4_D101st2nd_highlow_downsampled_URD, file = "D4_D101st2nd_highlow_downsampled_URD.Robj")

#define tip cells; here we choose all clusters from D10
D4_D101st2nd_downsampled <- SetAllIdent(object = D4_D101st2nd_downsampled, id = "res.2")
tipcells <- SubsetData(D4_D101st2nd_downsampled, ident.use = c(4,12,10,18,13,25,23,28,1,19,30,11,26,22,9))
#convert tipcells object to URD
tipcellsURD <- seuratToURD(tipcells)
# Copy cluster identities from tipcells URD object to a new clustering ("tip.clusters") in the full URD object.
D4_D101st2nd_highlow_downsampled_URD@group.ids[rownames(tipcellsURD@group.ids), "tip.clusters"] <- tipcellsURD@group.ids$res.2

# Determine the parameters of the logistic used to bias the transition probabilities. The procedure
# is relatively robust to this parameter, but the cell numbers may need to be modified for larger
# or smaller data sets.
D4_D101st2nd_highlow_downsampled_URD.ptlogistic <- pseudotimeDetermineLogistic(D4_D101st2nd_highlow_downsampled_URD, "pseudotime", optimal.cells.forward=50, max.cells.back=100, do.plot = T)

#save(combined.floods, file="D4_D101st2nd_highlowBMP_combinedfloods.Robj")

# Bias the transition matrix acording to pseudotime
D4_D101st2nd_highlow_downsampled_URD.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(D4_D101st2nd_highlow_downsampled_URD, "pseudotime", logistic.params=D4_D101st2nd_highlow_downsampled_URD.ptlogistic))

# Simulate the biased random walks from each tip
D4_D101st2nd_highlow_downsampled_URD.walks <- simulateRandomWalksFromTips(D4_D101st2nd_highlow_downsampled_URD, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = D4_D101st2nd_highlow_downsampled_URD.biased.tm, 
                                                                         n.per.tip = 25000, root.visits = 1, max.steps = 2500, verbose = F)

#3906 cells not visited
# Process the biased random walks into visitation frequencies
D4_D101st2nd_highlow_downsampled_URD <- processRandomWalksFromTips(D4_D101st2nd_highlow_downsampled_URD, D4_D101st2nd_highlow_downsampled_URD.walks, verbose = F)

# Load the cells used for each tip into the URD object
highlowBMP.tree <- loadTipCells(D4_D101st2nd_highlow_downsampled_URD, "tip.clusters")

# Build tree1, 6777 cells not visited
highlowBMP.tree <- buildTree(highlowBMP.tree, pseudotime = "pseudotime", tips.use = c(4,12,10,18,13,25,23,28,1,19,30,11,26,22,9), divergence.method = "preference", 
                         cells.per.pseudotime.bin = 80, bins.per.pseudotime.window = 8, min.cells.per.segment = 10,
                         p.thresh=0.05)
#plotDim(D4_D101st2nd_highlow_downsampled_URD, label="res.0.7", na.rm=F)

save(highlowBMP.tree, file = "highlowBMP_D4_D6_D101st2nd.Robj")
save(D4_D101st2nd_highlow_downsampled_URD, file = "D4_D101st2nd_highlow_downsampled_URD.Robj")


#plot 2D tree
#plot all
plotTree(highlowBMP.tree, "niceorder", cell.alpha = 1.0,cell.size = 0.6,tree.alpha = 0.3,
discrete.colors = c("#006837", "#980043", "#DE2D26", "#006837","#2C7FB8","#980043","#FC9272","#78C679","#7FCDBB","#DF65B0","#FB6A4A","#31A354","#41B6C4","#DD1C77"))

#plot all timepoints and genotypes
plotTree(highlowBMP.tree, label.segments = T)
plotTree(highlowBMP.tree, "segment")

save(D4_D101st2nd_highlow_downsampled_URD, file="D4_D101st2nd_highlow_downsampled_URD.Robj")
save(highlowBMP.tree, file="combineddeepfull_tree.Robj")

#to draw feature plots of gene expression on tree
genestoplot <- c("Utf1", "Pou5f1", "T", "Tdgf1", "Mesp1", "Mixl1")
for (gene in genestoplot) {
  plot(plotTree(highlowBMP.tree, gene))
}

genestoplot <- c("Tgfb2")
for (gene in genestoplot) {
  plot(plotTree(highlowBMP.tree, gene))
}

#name the tips
new.seg.names <- c("cm-4","cm-12","cm-1","epi-10","fib-13","fib-25","fib-18","cm-19","neuro-26","hema-23","hema-28","neuro-9","tuj1-30","retin11","neuro-22") 
segs.to.name <- c("4","12","1","10","13","25","18","19","26","23","28","9","30","11","22")
highlowBMP.tree <- nameSegments(highlowBMP.tree, segments = segs.to.name, segment.names = new.seg.names)

#to draw feature plots of gene expression on tree
genestoplot <- c("Wt1","Hba-x","Tcf21","Sox2","Tubb3","Rax","Col1a1")
for (gene in genestoplot) {
  plot(plotTree(highlowBMP.tree, gene))
}

#get cells sorted niceorder in each segment
highlowBMP_URD_cellnumbers <- table(highlowBMP.tree@group.ids$niceorder, highlowBMP.tree@group.ids$segment)
write.csv(highlowBMP_URD_cellnumbers, file = "highlowBMP_URD_cellnumbers.csv")


#code from Daniel Carlin, to extract cell barcodes from a particular segment 
segment18<-list()
segment18[['intermediates']]<-c(highlowBMP.tree@tree$`cells.in.segment`$`18`)
#after this switch to differential gene test in Seurat
