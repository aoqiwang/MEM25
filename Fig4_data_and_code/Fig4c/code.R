library(MASS)
library(class)
library(cluster)
library(impute)
library(Hmisc)
library(WGCNA)
library(flashClust)
library(sva)
library(preprocessCore)
library(car)
library(pheatmap)
library(preprocessCore)
library(sva)
library(cluster)
library(gplots)
library(ggplot2)
library(edgeR)
library(tidyverse)
library(RColorBrewer)
options(stringsAsFactors = F)
dat0 <- read.table("log2_gene_expression.txt", header = T, sep = "\t", row.names = 1)
cdata <- as.matrix(dat0)
cdata2 <-data.frame(cdata)
datExpr0 = as.data.frame(t(cdata))
datExprFemale=datExpr0
sif <- read.table("sif.txt", header = T, sep = "\t", row.names = 1)
datExpr=datExpr0
sample=sif
set.seed(13)
powers = c(c(1:10), seq(from = 12, to=20, by=2)) #设置软阈值
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5) #网络拓扑分析
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1],
-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",
ylab="Scale Free Topology Model Fit,signed R^2",
type="n",
main = paste("Scale independence"))
text(sft$fitIndices[,1],
-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
plot(sft$fitIndices[,1],
sft$fitIndices[,5],
xlab="Soft Threshold (power)",
ylab="Mean Connectivity",
type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1],
sft$fitIndices[,5],
labels=powers,
cex=cex1,col="red")
datTraits=sample
weight = as.data.frame(datTraits$Stress_severity)
names(weight) = "weight"
GS.weight = as.numeric(cor(datExprFemale, weight, use = "p"))
mergingThresh=0.15
# previous analysis
moduleLabelsManual2 = matchLabels(moduleLabelsManual1, moduleLabelsAutomatic)
# Convert labels to colors for plotting
moduleColorsManual2 = labels2colors(moduleLabelsManual2)
###Clustering methods may identify modules whose expression profiles are very similar. More specifically, a module detection method may result in modules whose eigengenes are highly correlated. Since highly correlated modules are not distinct, it may be advisable to merge them. The following code shows how to create a cluster tree of module eigengenes and how to merge them (if their pairwise correlation is larger than 0.75).
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
adjacency = adjacency(datExpr, power =7)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = 'mcquitty')
minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2,pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = 'mcquitty')
MEDissThres = 0
plot(METree, main = 'Clustering of module eigengenes', xlab = '', sub = '')
abline(h=MEDissThres, col = 'red')
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c('Dynamic Tree Cut', 'Merged dynamic'), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
person=cor(datExpr,use = 'p')
corr<-TOM
Colors<-mergedColors
colnames(corr)<-colnames(datExpr)
rownames(corr)<-colnames(datExpr)
names(Colors)<-colnames(datExpr)
colnames(person)<-colnames(datExpr)
rownames(person)<-colnames(datExpr)
umc = unique(mergedColors)
lumc = length(umc)
for (i in c(1:lumc)){
if(umc[i]== 'grey'){
next
}
ME=MEs[, paste('ME',umc[i], sep='')]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,Colors==umc[i]])),nrgcols=30,rlabels=F,rcols=umc[i], main=umc[i], cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=umc[i], main='', cex.main=2,ylab='eigengene expression',xlab='array sample')
}
kME=signedKME(datExpr, mergedMEs, outputColumnName = 'kME', corFnc = 'cor', corOptions = "use='p',method='spearman'")
if (dim(datExpr)[2]>=1500) nSelect=1500 else nSelect=dim(datExpr)[2]
set.seed(1)



select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = 'average')
selectColors = moduleColors[select]
plotDiss = selectTOM^7
#TOMplot(plotDiss, selectTree, selectColors, main = 'Network heatmap plot')
TOMplot(plotDiss, selectTree,  main = 'Network heatmap plot')



merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
kME=signedKME(datExpr, mergedMEs, outputColumnName = 'kME', corFnc = 'cor', corOptions = "use='p',method='spearman'")
if (dim(datExpr)[2]>=1500) nSelect=1500 else nSelect=dim(datExpr)[2]
set.seed(1)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
selectTree = hclust(as.dist(selectTOM), method = 'average')
selectColors = moduleColors[select]
plotDiss = selectTOM^7
TOMplot(plotDiss, selectTree, selectColors, main = 'Network heatmap plot')



MEs = moduleEigengenes(datExpr, Colors)$eigengenes
MET = orderMEs(MEs)
plotEigengeneNetworks(MET, 'Eigengene adjacency heatmap', marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)


moduleTraitCor = cor(MET, sample, use = 'p')
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), '\n(',signif(moduleTraitPvalue, 1), ')', sep = '')
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(2, 4, 2, 1.5))
labeledHeatmap(Matrix = moduleTraitCor,
xLabelsAngle = 0,
cex.lab = 0.5,
xLabels = colnames(sample),
yLabels = names(MET),
ySymbols = names(MET),
colorLabels = FALSE,
colors = blueWhiteRed(100),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
yColorWidth=0.02,
xColorWidth = 0.05,
main = paste('Module-trait relationships'))


moduleTraitCor = cor(MET, sample, use = 'p')
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), '\n(',signif(moduleTraitPvalue, 1), ')', sep = '')
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(2, 4, 2, 1.5))
labeledHeatmap(Matrix = moduleTraitCor,
xLabelsAngle = 0,
cex.lab = 0.5,
xLabels = colnames(sample),
yLabels = names(MET),
ySymbols = names(MET),
colorLabels = FALSE,
colors = blueWhiteRed(6),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
yColorWidth=0.02,
xColorWidth = 0.05,
main = paste('Module-trait relationships'))


moduleTraitCor = cor(MET, sample, use = 'p')
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), '\n(',signif(moduleTraitPvalue, 1), ')', sep = '')
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(2, 4, 2, 1.5))
labeledHeatmap(Matrix = moduleTraitCor,
xLabelsAngle = 0,
cex.lab = 0.5,
xLabels = colnames(sample),
yLabels = names(MET),
ySymbols = names(MET),
colorLabels = FALSE,
colors = blueWhiteRed(6),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1.5,1.5),
yColorWidth=0.02,
xColorWidth = 0.05,
main = paste('Module-trait relationships'))
