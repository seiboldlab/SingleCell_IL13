#########################################
#############Seurat analysis#############
#########################################

library(Seurat)

#################################### First bring in the dataset and summary stats

#Bring in the expression matrix that has been combined across runs
IL13<-read.table("../data/chip97570_expressionMatrix.txt",sep="\t",row.names=1,header=T,stringsAsFactors=F)
sum1<-read.table("../data/IL13_11d_summaryStats.txt",sep="\t",row.names=1,header=T,stringsAsFactors=F)

#Toss controls
IL13<-IL13[,-grep("Ctrl",colnames(IL13))]
sum1<-sum1[-grep("Ctrl",rownames(sum1)),]

#Remove barcodes and X from sample names
colnames(IL13)<-sapply(strsplit(sapply(strsplit(colnames(IL13),"X"),function(x)x[2]),"_"),function(x)x[1])
all(rownames(sum1) == colnames(IL13))





###################################### Filter

#Remove outliers (13 cells removed; 789 cells remaining)
outliers<-names(colSums(IL13))[which(colSums(IL13) > 60000)] #removes 10 cells
outliers<-append(outliers,names(colSums(IL13))[which(colSums(IL13) < 2000)]) #Removes 3 cells
sum1f<-sum1[-which(colnames(IL13) %in% outliers),]
IL13f<-IL13[,-which(colnames(IL13) %in% outliers)]

##Toss mtDNA pseudogenes, ribosomal genes, and mt ribosomal genes
IL13f<-IL13f[-c(grep("MTAT|MT-|MTCO|MTCY|MTERF|MTND|MTRF|MTRN|MRPL|MRPS|RPL|RPS",rownames(IL13f))),]

#Filter the genes based on expression (all genes expressed in at least 0.5 percent of the cells)
IL13f<-IL13f[which(rowSums(IL13f > 0) / ncol(IL13f) >= 0.005),]








###################################### SEURAT analysis on filtered dataset - aligning the two donors

###Now import data into Seurat
IL13s<-CreateSeuratObject(raw.data = IL13f, min.cells = 0, min.genes = 0, project = "IL13_expression")

#Normalize
IL13s<-NormalizeData(object = IL13s, normalization.method = "LogNormalize", scale.factor = 10000)

#Add metadata
#1 = T71 BSA
#2 = T71 IL13
#3 = T72 BSA
#4 = T72 IL13

###Now we want to break up the dataset into T71 and T72 and then import and scale these separately
IL13_T71<-CreateSeuratObject(raw.data = IL13f[,grep("97570.1|97570.2",colnames(IL13f))], 
	min.cells = 1, min.genes = 1, project = "IL13_T71_expression")
IL13_T71<-NormalizeData(object = IL13_T71, normalization.method = "LogNormalize", scale.factor = 10000)
IL13_T71 <- ScaleData(object = IL13_T71)

IL13_T72<-CreateSeuratObject(raw.data = IL13f[,grep("97570.3|97570.4",colnames(IL13f))], 
	min.cells = 1, min.genes = 1, project = "IL13_T72_expression")
IL13_T72<-NormalizeData(object = IL13_T72, normalization.method = "LogNormalize", scale.factor = 10000)
IL13_T72 <- ScaleData(object = IL13_T72)

#Add metadata column distinguishing IL13 from BSA control
treat.col<-data.frame(rep("BSA",nrow(IL13_T71@meta.data)))
treat.col[,1]<-as.character(treat.col[,1])
treat.col[grep("97570.2",rownames(IL13_T71@meta.data)),]<-"IL13"
rownames(treat.col)<-rownames(IL13_T71@meta.data)
colnames(treat.col)<-"treatment" 
IL13_T71<-AddMetaData(IL13_T71,metadata=treat.col,col.name="treatment")

treat.col<-data.frame(rep("BSA",nrow(IL13_T72@meta.data)))
treat.col[,1]<-as.character(treat.col[,1])
treat.col[grep("97570.4",rownames(IL13_T72@meta.data)),]<-"IL13"
rownames(treat.col)<-rownames(IL13_T72@meta.data) 
colnames(treat.col)<-"treatment" 
IL13_T72<-AddMetaData(IL13_T72,metadata=treat.col,col.name="treatment")


#Now find variable genes
IL13_T71 <- FindVariableGenes(object = IL13_T71, do.plot = F)
IL13_T72 <- FindVariableGenes(object = IL13_T72, do.plot = F)

#And take the union of the top 2000 most variable genes for the two datasets
hvg.T71 <- rownames(x = head(x = IL13_T71@hvg.info, n = 2000))
hvg.T72 <- rownames(x = head(x = IL13_T72@hvg.info, n = 2000))
hvg.union <- union(x = hvg.T71, y = hvg.T72)

#Lastly, we set the 'protocol' in each dataset for easy identification
IL13_T71@meta.data[, "protocol"] <- "T71"
IL13_T72@meta.data[, "protocol"] <- "T72"








###################################### Run canonical correlation analysis
#Identify common sources of variation between the two datasets.
IL13s<-RunCCA(object = IL13_T71, object2 = IL13_T72, genes.use = hvg.union)

#Visualize results of CCA plot CC1 versus CC2
p1 <- DimPlot(object = IL13s, reduction.use = "cca", group.by = "protocol", pt.size = 2, 
    do.return = TRUE)
p2 <- VlnPlot(object = IL13s, features.plot = "CC1", group.by = "protocol", do.return = TRUE)

pdf("IL13s_CCA_T71_T72_seurat.pdf")
plot_grid(p1, p2)
dev.off()

#Look at the extreme loading genes for each CC dimension
PrintDim(object = IL13s, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

#Choosing the number of CCs for downstream analysis and aligning. 
DimHeatmap(object = IL13s, reduction.type = "cca", cells.use = 500, dim.use = 1:9, 
    do.balanced = TRUE)

#We should search for cells whose expression level cannot be well-explained by low-dimensional CCA,
#compared to low-dimensional PCA
IL13s <- CalcVarExpRatio(object = IL13s, reduction.type = "pca", grouping.var = "protocol", 
    dims.use = 1:3)

#Discard cells where the variance explained by CCA is < 2-fold (ratio < 0.5) compared to PCA
IL13s.all.save <- IL13s
IL13s <- SubsetData(object = IL13s, subset.name = "var.ratio.pca", accept.low = 0.5) #33 cells tossed

#Now align the CCA subspaces
IL13s <- AlignSubspace(object = IL13s, reduction.type = "cca", grouping.var = "protocol", 
    dims.align = 1:3)

#Make sure the combined dataset is scaled to regress out nUMI
IL13s<-ScaleData(IL13s,vars.to.regress=c("nUMI"))

#Visualize the aligned CCA
p1 <- VlnPlot(object = IL13s, features.plot = "ACC1", group.by = "protocol", do.return = TRUE)
p2 <- VlnPlot(object = IL13s, features.plot = "ACC2", group.by = "protocol",  do.return = TRUE)
plot_grid(p1, p2)

#Run a single integrated analysis on all cells
IL13s <- RunTSNE(object = IL13s, reduction.use = "cca.aligned", dims.use = 1:3, do.fast = TRUE, perplexity=30)
IL13s <- FindClusters(object = IL13s, reduction.type = "cca.aligned", dims.use = 1:3, 
	resolution = 0.6, save.SNN = TRUE)

pdf("IL13s_10day_aligned_withDonor_tSNE_6clusters.pdf",height=3,width=11)
#dev.new(height=3,width=11)
p1 <- TSNEPlot(object = IL13s, do.return = TRUE, pt.size = 				
	0.8,colors.use=c("yellow2","darkgreen","saddlebrown","blue","seagreen2","purple"))
p2 <- TSNEPlot(object = IL13s, group.by = "protocol", colors.use=c("grey","black"),do.return = TRUE, pt.size = 0.8)
p3<-TSNEPlot(object = IL13s, colors.use=c("grey","black"),group.by = "treatment", 
	do.return = TRUE, pt.size = 0.8)

plot_grid(p1,p2,p3,nrow=1)
dev.off()

pdf("IL13s_10day_aligned_tSNE_6clusters.pdf",height=2.5,width=6)
dev.new(height=2.5,width=6)
p1 <- TSNEPlot(object = IL13s, do.return = TRUE, pt.size = 				
	0.8,colors.use=c("yellow2","darkgreen","saddlebrown","blue","seagreen2","purple"),no.legend=F)
p3<-TSNEPlot(object = IL13s, colors.use=c("grey","black"),group.by = "treatment", 
	do.return = TRUE, pt.size = 0.8,no.legend=F)
plot_grid(p1,p3,nrow=1)
dev.off()











######################################################Subclustering ciliated cells

#Create ciliated data subset
cilCells<-rownames(IL13s@meta.data)[which(IL13s@meta.data$res.0.6 == 1)]
IL13s_cil<-SubsetData(IL13s,cells.use=cilCells)

#########Cluster 1 (cilated)
#Select genes
IL13s_cil<-FindVariableGenes(IL13s_cil ,mean.function = ExpMean, dispersion.function = LogVMR, 
	x.low.cutoff = 0.4, x.high.cutoff = 8, y.cutoff = 0.6, do.contour=F)


###########Performing linear dimensional reduction
#Perform PCA on the scaled data 
IL13s_cil<-RunPCA(IL13s_cil, pc.genes = IL13s_cil@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)

#Clustering
IL13s_cil<-FindClusters(IL13s_cil,dims.use = 1:3,resolution = 0.7,print.output = 0,save.SNN = T,force.recalc=T)
IL13s_cil<-RunTSNE(IL13s_cil,dims.use = 1:3, do.fast = T,max_iter=1000,perplexity=23)

#Plot these groups
#pdf('IL13s_cil_tSNE_byCluster_seurat.pdf')
TSNEPlot(IL13s_cil,colors.use=c("tomato","cornflowerblue","green3","purple","yellow","darkgreen","red"),pt.size=2.5)
dev.off()




######## Now create a new meta.data column with assignments to these four ciliated clusters merged with the 
######## secretory clusters from the entire dataset

#Start new cluster vector
clusterVec<-data.frame("clusters_8"=IL13s@meta.data$res.0.6,row.names=rownames(IL13s@meta.data),stringsAsFactors=F)
#Change current cluster names
clusterVec[which(IL13s@meta.data$res.0.6 == 5),]<- 1
clusterVec[which(IL13s@meta.data$res.0.6 == 1),]<- 3
clusterVec[rownames(IL13s_cil@meta.data[which(IL13s_cil@meta.data$res.0.7 == 1),]),1]<-4
clusterVec[rownames(IL13s_cil@meta.data[which(IL13s_cil@meta.data$res.0.7 == 0),]),1]<-5
clusterVec[rownames(IL13s_cil@meta.data[which(IL13s_cil@meta.data$res.0.7 == 3),]),1]<-6
clusterVec[rownames(IL13s_cil@meta.data[which(IL13s_cil@meta.data$res.0.7 == 2),]),1]<-7
#Add to metadata
IL13s<-AddMetaData(IL13s,metadata=clusterVec,col.name="clusters_8")

#Plot the new composite clusters
pdf("IL13s_10day_aligned_tSNE_8clusters.pdf",height=2.5,width=6)
dev.new(height=5,width=5.5)
p1 <- TSNEPlot(object = IL13s, do.return = TRUE, group.by="clusters_8",pt.size = 4,
	colors.use=c("yellow2","purple","saddlebrown","darkgreen","cornflowerblue","pink","blue","red3"),no.legend=T)
p3<-TSNEPlot(object = IL13s, colors.use=c("grey","black"),group.by = "treatment", 
	do.return = TRUE, pt.size = 4,no.legend=T)
plot_grid(p1,p3,nrow=1)
dev.off()


pdf("IL13s_10day_aligned_withDonor_tSNE_8clusters.pdf",height=3,width=11)
dev.new(height=2.5,width=10)
p1 <- TSNEPlot(object = IL13s, do.return = TRUE, group.by="clusters_8",pt.size = 1.2,
	colors.use=c("yellow2","darkorange","yellow4","saddlebrown","springgreen","#00c100","springgreen4","#003d00"),no.legend=F)
p2 <- TSNEPlot(object = IL13s, group.by = "protocol", colors.use=c("grey","black"),do.return = TRUE, pt.size = 1.2)
p3<-TSNEPlot(object = IL13s, colors.use=c("grey","black"),group.by = "treatment", 
	do.return = TRUE, pt.size = 1.2)
plot_grid(p1,p2,p3,nrow=1)
dev.off()





