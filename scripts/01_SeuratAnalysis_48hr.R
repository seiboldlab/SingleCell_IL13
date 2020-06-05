#########################################
#############Seurat analysis#############
#########################################
library(Seurat)

############### PRELIMINARY ANALYSIS

### Bring in raw counts for the two 48 hour chips
merged_1<-read.table("..data/chip91497_expressionMatrix.txt")
merged_2<-read.table("../data/chip94576_expressionMatrix.txt")

#Toss controls
merged_1<-merged_1[,-grep("Ctrl",colnames(merged_1))]
merged_2<-merged_2[,-grep("Ctrl",colnames(merged_2))]

#Remove barcodes and X from sample names
colnames(merged_1)<-sapply(strsplit(sapply(strsplit(colnames(merged_1),"X"),function(x)x[2]),"_"),function(x)x[1])
colnames(merged_2)<-sapply(strsplit(sapply(strsplit(colnames(merged_2),"X"),function(x)x[2]),"_"),function(x)x[1])

#Combine the chips
IL13<-cbind(merged_1,merged_2)

#Combine mapping summaries across chips and toss controls
sum_all<-read.table("../data/IL13_48hr_summaryStats.txt",header=T,stringsAsFactors=F)
sum_all<-sum_all[-grep("Ctrl",rownames(sum_all)),]
all(rownames(sum_all) == colnames(IL13))

#Do cell filtering
IL13f<-IL13
sum_allf<-sum_all

# Remove outliers
outliers<-names(colSums(merged_1))[which(colSums(merged_1) > 46000)]
outliers<-append(outliers,names(colSums(merged_2))[which(colSums(merged_2) > 28500)])
sum_allf<-sum_allf[-which(colnames(IL13f) %in% outliers),]
IL13f<-IL13f[,-which(colnames(IL13f) %in% outliers)]

# Remove all remaining cells with % mapped and % mapped to genes values < 50
IL13f<-IL13f[,-unique(c(which(sum_allf$pct_mapped_reads < 50),which(sum_allf$pct_reads_in_gene < 50)))]
sum_allf<-sum_allf[-unique(c(which(sum_allf$pct_mapped_reads < 50),which(sum_allf$pct_reads_in_gene < 50))),]

##Toss mtDNA pseudogenes, ribosomal genes, and mt ribosomal genes
IL13f<-IL13f[-c(grep("^MTAT|^MT-|^MTCO|^MTCY|^MTERF|^MTND|^MTRF|^MTRN|^MRPL|^MRPS|^RPL|^RPS",rownames(IL13f))),]

#Filter the genes based on expression (all genes expressed in at least 0.5 percent of the cells)
IL13f<-IL13f[which(rowSums(IL13f > 0) / ncol(IL13f) >= 0.005),]

#Finally, remove cells with fewer than 2,500 UMIs
sum_allf<-sum_allf[-which(colSums(IL13f) < 2500),]
IL13f<-IL13f[,-which(colSums(IL13f) < 2500)]












########################### SEURAT ANALYSIS

###Now import data into Seurat
IL13s<-CreateSeuratObject(raw.data = IL13f, min.cells = 0, min.genes = 0, project = "IL13s_expression")

#Normalize
IL13s<-NormalizeData(object = IL13s, normalization.method = "LogNormalize", scale.factor = 10000)

##### Add metadata
#1 = T71 BSA
#2 = T72 BSA
#3 = T71 IL13
#4 = T72 IL13

#Add metadata column distinguishing IL13 from BSA control
treat.col<-data.frame(rep("BSA",nrow(IL13s@meta.data)))
treat.col[,1]<-as.character(treat.col[,1])
treat.col[grep("91497.3",rownames(IL13s@meta.data)),]<-"IL13"
treat.col[grep("91497.4",rownames(IL13s@meta.data)),]<-"IL13"
treat.col[grep("94576.3",rownames(IL13s@meta.data)),]<-"IL13"
treat.col[grep("94576.4",rownames(IL13s@meta.data)),]<-"IL13"
rownames(treat.col)<-rownames(IL13s@meta.data) #Need to have same row names as the object before adding
colnames(treat.col)<-"treatment" #Need to pre-add the colname (make it same as specified in function call below)
IL13s<-AddMetaData(IL13s,metadata=treat.col,col.name="treatment")

#Add metadata column distinguishing the two donors
treat.col<-data.frame(rep("T71",nrow(IL13s@meta.data)))
treat.col[,1]<-as.character(treat.col[,1])
treat.col[grep("91497.2",rownames(IL13s@meta.data)),]<-"T72"
treat.col[grep("91497.4",rownames(IL13s@meta.data)),]<-"T72"
treat.col[grep("94576.2",rownames(IL13s@meta.data)),]<-"T72"
treat.col[grep("94576.4",rownames(IL13s@meta.data)),]<-"T72"
rownames(treat.col)<-rownames(IL13s@meta.data) #Need to have same row names as the object before adding
colnames(treat.col)<-"donor" #Need to pre-add the colname (make it same as specified in function call below)
IL13s<-AddMetaData(IL13s,metadata=treat.col,col.name="donor")

#Add metadata column distinguishing the two chips
treat.col<-data.frame(rep("91497",nrow(IL13s@meta.data)))
treat.col[,1]<-as.character(treat.col[,1])
treat.col[grep("94576",rownames(IL13s@meta.data)),]<-"94576"
rownames(treat.col)<-rownames(IL13s@meta.data) #Need to have same row names as the object before adding
colnames(treat.col)<-"chip" #Need to pre-add the colname (make it same as specified in function call below)
IL13s<-AddMetaData(IL13s,metadata=treat.col,col.name="chip")

#Scale data
##For the purpose of visualization, want to regress out unwanted sources of variation
IL13s <- ScaleData(object = IL13s, vars.to.regress = c("nUMI","donor","chip"))

#Select genes
IL13s<-FindVariableGenes(IL13s ,mean.function = ExpMean, dispersion.function = LogVMR, 
	x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.1, do.contour=F)

#Perform PCA on the scaled data 
IL13s<-RunPCA(IL13s, pc.genes = IL13s@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
PCElbowPlot(IL13s)

#Clustering
IL13s<-FindClusters(IL13s,dims.use = 1:8,resolution = 1.1,print.output = 0,save.SNN = T,force.recalc=T)
IL13s<-RunTSNE(IL13s,dims.use = 1:8, do.fast = T,max_iter=1000,perplexity=45)

#Order and name clusters
treat.col<-data.frame("clusters_9_newOrder"=IL13s@meta.data$res.1.1,row.names=rownames(IL13s@meta.data))
treat.col[,1]<-as.character(treat.col[,1])
treat.col$clusters_9_newOrder[which(IL13s@meta.data$res.1.1 == "0")]<-"c6"
treat.col$clusters_9_newOrder[which(IL13s@meta.data$res.1.1 == "1")]<-"c5"
treat.col$clusters_9_newOrder[which(IL13s@meta.data$res.1.1 == "2")]<-"c3"
treat.col$clusters_9_newOrder[which(IL13s@meta.data$res.1.1 == "3")]<-"c4"
treat.col$clusters_9_newOrder[which(IL13s@meta.data$res.1.1 == "4")]<-"c8"
treat.col$clusters_9_newOrder[which(IL13s@meta.data$res.1.1 == "5")]<-"c7"
treat.col$clusters_9_newOrder[which(IL13s@meta.data$res.1.1 == "6")]<-"c9"
treat.col$clusters_9_newOrder[which(IL13s@meta.data$res.1.1 == "7")]<-"c2"
treat.col$clusters_9_newOrder[which(IL13s@meta.data$res.1.1 == "8")]<-"c1"
IL13s<-AddMetaData(IL13s,metadata=treat.col,col.name="clusters_9_newOrder")

#Plot color based on clusters
IL13s<-SetIdent(IL13s,ident.use = IL13s@meta.data$clusters_9_newOrder)
pdf('IL13s_tSNE_byCluster_seurat.pdf',height=5,width=4)
dev.new(height=5,width=4)
TSNEPlot(IL13s,colors.use=c("darkturquoise","palegreen","springgreen3","darkgreen","yellow2","yellow4","darkorange",
	"saddlebrown","midnightblue"),pt.size=2.5,no.legend=T)
dev.off()

#Plot color based on the treatment
pdf('IL13s_tSNE_byTreatment_seurat.pdf',height=5,width=4)
dev.new(height=5,width=4)
TSNEPlot(IL13s,colors.use=c("grey","black"),group.by="treatment",pt.size=2.5,no.legend=T)
dev.off()

#Plot color based on the donor
pdf('IL13s_tSNE_byDonor_seurat_noLegend.pdf',height=5,width=4)
TSNEPlot(IL13s,colors.use=c("grey","black"),group.by="donor",pt.size=2.5,no.legend=T)
dev.off()

#Plot color based on the chip)
pdf('IL13s_tSNE_byChip_seurat_noLegend.pdf',height=5,width=4)
TSNEPlot(IL13s,colors.use=c("grey","black"),group.by="chip",pt.size=2.5,no.legend=T)
dev.off()




