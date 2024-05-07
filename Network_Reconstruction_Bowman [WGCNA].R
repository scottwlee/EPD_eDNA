source("http://bioconductor.org/biocLite.R")
source("https://bioconductor.org/biocLite.R")


##########################################################################################
########## Preparation ##########
##########################################################################################
### Package installation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("GO.db", "preprocessCore", "impute"))

install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "reshape", "fastcluster", "dynamicTreeCut", "survival"))
install.packages("WGCNA")
install.packages("vegan")

### Working directory setup
setwd("D://")

### Package loading
library("WGCNA")
library("vegan")


##########################################################################################
########## Input data ##########
##########################################################################################
data<-read.table("XXX.txt",header=T,na.strings="NA",row.names=1)
HellingerData<-decostand(data,method = "hellinger")
write.table(t(HellingerData), file = "YYY.txt", sep="\t")

rt=read.table("YYY.txt",sep="\t",header=T)
dim(rt)

### WGCNA, "Code chunk 2"
datExpr0 = as.data.frame(t(rt[,-1]))
names(datExpr0) = rt[,1] # If the first line is not named by ID, write it as "rt[,1]" 
rownames(datExpr0) = names(rt[,-1])


##########################################################################################
########## Check missing value and filter ##########
##########################################################################################
### WGCNA, "Code chunk 3"/ Check missing value
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

# if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


##########################################################################################
########## sample cluster ##########
##########################################################################################
sampleTree = hclust(dist(datExpr0), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small

# SizeGrWindow (12,9)
pdf(file = "1_sample_cluster.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

### Plot a line to show the cut
## abline(h = 10000, col = "red")
dev.off()

### Determine cluster under the line 
## clust = cutreeStatic(sampleTree, cutHeight = 10000, minSize = 10)
## table(clust)
## keepSamples = (clust==1)
## datExpr0 = datExpr0[keepSamples, ]
### (In general, not used unless the figure is quite complex and simplification is necessary)


##########################################################################################
########## Input trait data/ Environmental Factors ##########
##########################################################################################
### Loading environmental factors data
traitData = read.table("environmental_factors.txt",row.names=1,header=T,comment.char = "",check.names=F)
dim(traitData)
names(traitData)
# Remove columns with information we don't need

### Form a data frame analogous to expression data that will hold the environmental factors
fpkmSamples = rownames(datExpr0)
traitSamples =rownames(traitData)
traitRows = match(fpkmSamples, traitSamples)
datTraits = traitData[traitRows,]
rownames(datTraits) 
collectGarbage()

### Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath

# SizeGrWindow (12,12)
pdf(file="2_sample_heatmap.pdf",width=12,height=12)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()


##########################################################################################
########## Network construction ##########
##########################################################################################
enableWGCNAThreads()

# Choose a set of soft-thresholding powers
powers = c(1:20)

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# Plot the results
# SizeGrWindow (9,5)
pdf(file="3_scale_independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.75

### Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.75,col="red") # This line corresponds to using the R^2 cut-off of h

### Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()



### Choose the soft power
# Creates connection strengths between pairs of samples
# Calculates (either correlation or distance) network adjacency from given expression data 

sft # Based on above plot we chose
softPower = # e.g., 1, 2, 3, ...
adjacency = adjacency(datExpr0, power = softPower) # Power means to specify soft thresholding power
# Input sample expression data (datExpr0)
# In correlation network, the adjacency is constructed from the correlations (values between -1 and 1, high numbers meaning high similarity)

### Turn adjacency into topological overlap
# This is a measurement of the true interconnectedness
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

### Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
# SizeGrWindow (12,9)
pdf(file="4_gene_clustering.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.4)
dev.off()



### Set the minimum module size/ We like large modules, so we set the minimum module size relatively high
minModuleSize = # e.g., 1, 2, 3, ...

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

### The following command gives the module labels and the size of each module/ Label 0 is reserved for unassigned genes
table(dynamicMods)

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# Plot the dendrogram and colors underneath
pdf(file="5_Dynamic_Tree.pdf",width=8,height=6) # SizeGrWindow (8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()



### Merging of modules whose expression profiles are very similar

# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
pdf(file="6_Clustering_module.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.75 # We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge/ Height setting can be modified
abline(h=MEDissThres, col = "red")
dev.off()



# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors = merge$colors

# Eigengenes of the new merged modules
mergedMEs = merge$newMEs

# To see what the merging did to our module colors, we plot the gene dendrogram again, with the original and merged module colors underneath
pdf(file="7_merged_dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()



moduleColors = mergedColors # Rename to moduleColors

# Construct numerical labels corresponding to the colors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs


##########################################################################################
########## Relate modules to external clinical traits ##########
##########################################################################################
# Define numbers of genes and samples
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

pdf(file="8_Module_trait.pdf",width=10,height=6) # SizeGrWindow (10,6)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               setStdMargins = FALSE,
               cex.text = 0.4,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


##########################################################################################
########## Define variable weight containing all column of datTraits ##########
##########################################################################################
### GS and MM
modNames = substring(names(MEs), 3) # Names (colors) of the modules

geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

traitNames=names(datTraits)

geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

# Each row corresponds to a module eigengene, column to a trait



### Plot GS vs MM for each trait vs each module

## Example; royalblue and CK
# module="royalblue"
# column = match(module, modNames)
# moduleGenes = moduleColors==module

# trait="CK"
# traitColumn=match(trait,traitNames)

# SizeGrWindow (7,7)

# par(mfrow = c(1,1))
## verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                     abs(geneTraitSignificance[moduleGenes, traitColumn]),
#                     xlab = paste("Module Membership in", module, "module"),
#                     ylab = paste("Gene significance for ",trait),
#                     main = paste("Module membership vs. gene significance\n"),
#                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
##

picDir="module_trait"
dir.create(picDir)
for (trait in traitNames){
  traitColumn=match(trait,traitNames)  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      pdfFile=paste("9_", trait, "_", module,".pdf",sep="")
      outPdf=paste(picDir,pdfFile,sep="\\")
      pdf(file=outPdf,width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}



names(datExpr0)
probes = names(datExpr0)


##########################################################################################
########## Export GS and MM ##########
##########################################################################################
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.xls",sep="\t",row.names=F)


##########################################################################################
########## visualizing the gene network ##########
##########################################################################################
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7
diag(plotTOM) = NA

### Call the plot function
pdf(file="10_allgene_heatmap.pdf",width=9, height=9) # SizeGrWindow (9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()



nSelect = 400
set.seed(10) ### For reproducibility, we set the random seed
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster

selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]

# Open a graphical window
# SizeGrWindow (9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# The color palette; Setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7
diag(plotDiss) = NA

pdf(file="11_selectgene_heatmap.pdf",width=9, height=9)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
dev.off()


##########################################################################################
########## visualizing the gene network of eigengenes ##########
##########################################################################################
# SizeGrWindow (5,7.5)
pdf(file="12_module_dendrogram.pdf",width=6, height=6)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()



### Or devide into two parts
# Plot the dendrogram
# SizeGrWindow (6,6)
pdf(file="13_module_heatmap.pdf",width=6, height=6)
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()



pdf(file="14_dendrogram_heatmap.pdf", width=5, height=7.5) # Plot the heatmap matrix (Note; This plot will overwrite the dendrogram plot)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
dev.off()


##########################################################################################
########## Exporting to Cytoscape ##########
##########################################################################################
cytoDir="CytoscapeInput"
dir.create(cytoDir)
for (mod in 1:nrow(table(moduleColors)))
{  
  modules = names(table(moduleColors))[mod]  
  probes = names(datExpr0)
  inModule = (moduleColors == modules)
  modProbes = probes[inModule]
  modGenes = modProbes  
  modTOM = TOM[inModule, inModule] 
  dimnames(modTOM) = list(modProbes, modProbes)
  edges_File = paste("CytoscapeInput-edges-", modules , ".txt", sep="")
  nodes_File = paste("CytoscapeInput-nodes-", modules, ".txt", sep="")
  outEdge=paste(cytoDir,edges_File,sep="\\")
  outNode=paste(cytoDir,nodes_File,sep="\\")
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = outEdge,
                                 nodeFile = outNode,
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])
}