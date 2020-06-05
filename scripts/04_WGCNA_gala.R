##################### WGCNA on GALA II #################

library(DESeq2)
library(edgeR)
library(WGCNA)
library(dendextend)
library(heatmap3)
library(ppcor)

#####################
### PRELIMINARIES ###
#####################

### Load raw GALA II count matrix
raw_counts <- read.table("..data/GALA_raw_counts.txt", header=T, sep="\t", row.names=1)
raw_counts <- as.matrix(raw_counts)

### Vst normalize data
design <- data.frame(row.names = colnames(raw_counts), subject = colnames(raw_counts))
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = design,
  design = ~1)
vst <- varianceStabilizingTransformation(dds)
vstMat <- assay(vst)

### Size factor normalize data
expr_norm<-estimateSizeFactors(dds)
expr_norm<-counts(expr_norm, normalized=T)






### Remove lowly expressed genes
### 17473 genes
good_genes <- rowSums(raw_counts >= 10) >= (ncol(raw_counts) * .15)
vstMat_rm <- vstMat[good_genes,]

### Load phenotype data
phen <- read.table("./TABLES/fixed_phenotype_age_gender_bmi_asthma.csv", header=T, sep=",")
rownames(phen) <- phen$SubjectID

### 695 samples
int.samples <- intersect(rownames(phen), colnames(vstMat_rm))
traits <- phen[int.samples,]

### compute residual matrix (regressing out age and gender)
vstMat_res <- vstMat_rm
for(gene in rownames(vstMat_rm)) {
    fit <- lm(vstMat_rm[gene,]~ traits$age + as.factor(traits$Male))
    res <- residuals(fit)
    vstMat_res[gene,] <- res
}

write.table(vstMat_res, file="./TABLES/695_expr_vst_res.txt", sep='\t', quote=F)







#############
### WGCNA ###
#############

####check genes
expMat <- t(vstMat_res)
gsg <- goodSamplesGenes(expMat, verbose=3)

if (!gsg$allOK)
{
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", paste(colnames(expMat)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:", paste(rownames(expMat)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    expMat = expMat[gsg$goodSamples, gsg$goodGenes]
}

#Now cluster donors (in contrast to clustering genes later on...)
sampleTree = hclust(dist(expMat), method = "average")

# Plot the sample tree; tips will be ordered to some extent by library size,
# Look for outliers along this continuum
pdf("./IMAGES/695.res.SampleTree_outliers.wgcna.pdf", height=48, width=64)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
dev.off()

enableWGCNAThreads()

#First, choose a set of candidate powers to look at
powers = c(1:20)

# Call the network topology analysis function
sft = pickSoftThreshold(expMat, powerVector = powers, verbose = 5, networkType="signed")

# Plot the results:
pdf("./IMAGES/695.res.pickSoftThresholdPlots.wgcna.pdf")
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# pick soft threshold = 9
softPower <- 9

####step by step WGCNA

#1. Create Similarity Matrix
pearson <- WGCNA::cor(as.matrix(expMat),method="pearson")

#2. Convert Similarity Matrix to Adjacency Matrix using Adjacency Function
adjacency.p <- adjacency.fromSimilarity(pearson,type = "signed",power=softPower)

#3. Convert Adjacency to TOM dissimilarity
TOMdissim.p <- 1 - TOMsimilarity(adjacency.p,TOMType = "signed",TOMDenom = "min")

#4. Perform hierarchical clustering on the dissimilarity matrix
geneTree = hclust(as.dist(TOMdissim.p), method = "average")

deepsplit <- 0.82

#5. cut tree based on hierarchical clustering
modules = cutreeDynamic(dendro = geneTree,method='hybrid', distM = TOMdissim.p,pamStage=F, pamRespectsDendro = F,
maxCoreScatter=min(geneTree$height)+deepsplit*(max(geneTree$height)-min(geneTree$height)),
minGap=(1-deepsplit)*(max(geneTree$height)-min(geneTree$height)),
cutHeight=quantile(geneTree$height,.99), minClusterSize=30)

modcolors = labels2colors(modules)

table(modcolors)

# dendogram
pdf("./IMAGES/695.res.dendrogram_signed_wgcna.pdf",width=14,height=7)
#plotDendroAndColors(dendro=geneTree,colors=cbind(modcolors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE,main=' ')
plotDendroAndColors(dendro=geneTree,colors=modcolors, dendroLabels = FALSE,main=' ')
dev.off()

#######
MEs <- moduleEigengenes(expMat, modcolors)$eigengenes
rownames(MEs) <- rownames(expMat)
MEs_no_grey <- MEs[,-which(colnames(MEs) %in% c("MEgrey"))]

write.table(MEs_no_grey,  file="./TABLES/695.res.WGCNA.ME.txt", sep='\t', quote=F)

# create table for genes and their corresponding module
gene2module = data.frame(gene=colnames(expMat), module=modcolors)
gene2module<-na.omit(gene2module)
write.table(gene2module, file="./TABLES/695.res.WGCNA.gene2module.txt",sep="\t",quote=F,row.names=F)

datKME <- signedKME(expMat, MEs)
write.table(datKME, file="./TABLES/695.res.WGCNA.KME.txt",sep="\t",quote=F)

gene2module_with_cor <- gene2module
gene2module_with_cor$cor <- NA

for(i in unique(gene2module_with_cor$module)) {
    kME_name <- paste0("kME",i)
    idx <- which(gene2module_with_cor$module==i)
    gene.idx <- as.character(gene2module_with_cor[idx,"gene"])
    gene.idx <- gsub("-",".",gene.idx)
    gene2module_with_cor$cor[idx] <- datKME[gene.idx,kME_name]
    #print(kME_name)
}
write.table(gene2module_with_cor, file="./TABLES/695.res.WGCNA.gene2module.with.cor.txt",sep="\t",quote=F,row.names=F)

###hub genes
ADJ1 <- abs(cor(expMat, use="p"))^softPower
Alldegrees1=intramodularConnectivity(ADJ1, gene2module$module)
hub_genes<-data.frame()
for(i in 1:length(unique(gene2module$module))){
    Alldegrees1$Module = gene2module$module
    tmp = Alldegrees1[Alldegrees1$Module == unique(gene2module$module)[i], ]
    hub_genes<-rbind(hub_genes, head(tmp[order(tmp$kWithin, decreasing=T),], n=nrow(tmp)))
}

write.table(hub_genes, file="./TABLES/695.res.WGCNA.hub_genes.txt",sep="\t",quote=F)

#####module dendogram
moduleTree <- hclust(dist(t(MEs_no_grey)), method = "average")

pdf("./IMAGES/695.res.ModuleTree.wgcna.pdf")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(moduleTree, main = "Module clustering", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
dev.off()
