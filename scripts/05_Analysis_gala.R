####################GALA FIGURE 7
library(beeswarm)
library(Seurat)
library(plotrix)
library(heatmap3)
library(DESeq2)

#####################
### PRELIMINARIES ###
#####################

### Load raw GALA II count matrix
GALA_raw <- read.table("..data/GALA_raw_counts.txt", header=T, sep="\t", row.names=1)
GALA_raw <- as.matrix(GALA_raw)

### Vst normalize data
design <- data.frame(row.names = colnames(GALA_raw), subject = colnames(GALA_raw))
dds <- DESeqDataSetFromMatrix(
  countData = GALA_raw,
  colData = design,
  design = ~1)
vst <- varianceStabilizingTransformation(dds)
GALA_vst <- assay(vst)

### Size factor normalize data
GALA_expr_norm<-estimateSizeFactors(dds)
GALA_expr_norm<-counts(GALA_expr_norm, normalized=T)



### Read in phenotype data for GALA 
GALA_phenotype <- read.table("../data/695.GALA.final.phen.with.add.txt", row.names = 1, header = T, sep='\t')
rownames(GALA_phenotype) <- GALA_phenotype$SubjectID

### Read in WGCNA module eigengenes 
GALA_eig <- read.table("GALA_datasets/695.res.WGCNA.ME.txt", header = T)

#Order GALA datasets by the phenotype table
GALA_raw<-GALA_raw[,rownames(GALA_phenotype)]
GALA_expr_norm<-GALA_expr_norm[,rownames(GALA_phenotype)]
GALA_vst<-GALA_vst[,rownames(GALA_phenotype)]
GALA_eig<-GALA_eig[rownames(GALA_phenotype),]

#Add in the module eigengenes to GALA_phenotype (MEskyblue is the type 2 module)
GALA_phenotype<-cbind(GALA_phenotype,GALA_eig)





#####################
####### PLOTs #######
#####################


###1. Boxplots for type 2 low vs. high individuals, genes: MUC2, MUC5AC, MUC5B, MUC5AC/MUC5B, SCGB1A1, FOXJ1
#genes to plot
z <- c("MUC2", "MUC5AC", "MUC5B", "SCGB1A1", "FOXJ1")

MUC5AC_5B_ratio <- t(log(GALA_expr_norm["MUC5AC",rownames(GALA_phenotype)]/GALA_expr_norm["MUC5B",rownames(GALA_phenotype)]))

ratioVec<-c("less","less","greater")

pdf("GALA_mucins_boxplots.pdf", width = 6, height = 5)
dev.new(width = 6, height = 5)
par(mfrow=c(2,3),bty='n')
for(i in 1:3){
	boxplot(log(as.matrix(GALA_expr_norm)[z[i],rownames(GALA_phenotype)])~GALA_phenotype$type2_status, 
		names = c("Type 2 High", "Type 2 Low"), plot = TRUE, main = z[i], ylab = "Normalized Expression", 
		cex.lab = 0.75, cex.axis = 0.5, boxlwd = 0.25, lwd = 0.5, lwd.axis = 0.5, outline = FALSE, las = 1, cex.main = 1)
	beeswarm(log(as.matrix(GALA_expr_norm)[z[i],rownames(GALA_phenotype)])~GALA_phenotype$type2_status, 
		col = c("red", "blue"), pch = 16, cex = 0.35, add = TRUE, outline = FALSE)
	mtext(paste("p = ",formatC(wilcox.test(log(as.matrix(GALA_expr_norm)[z[i],rownames(GALA_phenotype)])[which(GALA_phenotype$type2_status == "type2_low")],
		log(as.matrix(GALA_expr_norm)[z[i],rownames(GALA_phenotype)])[which(GALA_phenotype$type2_status == "type2_high")],alternative=ratioVec[i])$p.value,
		format="e",digits=2),sep=""),side=3,line=0.2,at=2.5,adj=1,cex=1)
	}

boxplot(MUC5AC_5B_ratio~GALA_phenotype$type2_status, names = c("Type 2 High", "Type 2 Low"), plot = TRUE, main = "MUC5AC/MUC5B", 
	ylab = "Normalized Expression", cex.lab = 0.75, cex.axis = 0.5, boxlwd = 0.25, lwd = 0.5, lwd.axis = 0.5, outline = FALSE, las = 1, cex.main = 1)
beeswarm(MUC5AC_5B_ratio~GALA_phenotype$type2_status, col = c("red", "blue"), pch = 16, cex = 0.35, add = TRUE, outline = FALSE)
mtext(paste("p = ",formatC(wilcox.test(MUC5AC_5B_ratio[which(GALA_phenotype$type2_status == "type2_low")],
	MUC5AC_5B_ratio[which(GALA_phenotype$type2_status == "type2_high")],alternative="less")$p.value,
	format="e",digits=2),sep=""),side=3,line=0.2,at=2.5,adj=1,cex=1)

for(i in 4:5){
	boxplot(log(as.matrix(GALA_expr_norm)[z[i],rownames(GALA_phenotype)])~GALA_phenotype$type2_status, 
		names = c("Type 2 High", "Type 2 Low"), plot = TRUE, main = z[i], ylab = "Normalized Expression", cex.lab = 0.75, 
		cex.axis = 0.5, boxlwd = 0.25, lwd = 0.5, lwd.axis = 0.5, outline = FALSE, las = 1, cex.main = 1)
	beeswarm(log(as.matrix(GALA_expr_norm)[z[i],rownames(GALA_phenotype)])~GALA_phenotype$type2_status, 
		col = c("red", "blue"), pch = 16, cex = 0.35, add = TRUE, outline = FALSE)
	mtext(paste("p = ",formatC(wilcox.test(log(as.matrix(GALA_expr_norm)[z[i],rownames(GALA_phenotype)])[which(GALA_phenotype$type2_status == "type2_low")],
		log(as.matrix(GALA_expr_norm)[z[i],rownames(GALA_phenotype)])[which(GALA_phenotype$type2_status == "type2_high")],alternative="greater")$p.value,
		format="e",digits=2),sep=""),side=3,line=0.2,at=2.5,adj=1,cex=1)
	}
dev.off()










###2. Boxplots for mean expression of soluble apical secretome IL-13 DEPs comparing type 2 high vs. low individuals (will end up with two plots, 
## one for up and one for downregulated DEPs)

upregulated_apical_secretome_DEPs <- c("CLCA1", "ITLN1", "SERPINB2", "MUC2", "ITLN2", "CST1", "SERPINB10", "FCGBP", "ST6GAL1", "PRB2", 
	"HSPA1L", "FASN", "CORO1B", "PSME2", "KRT16", "ACTA1", "CDH6", "ALOX15", "CAP1", "SERPINB13", "ALPL", "ACTA2", "LYZ", "TFF3", "LGALS7", 
	"THBS1", "TIMP1", "TNC", "MSLN", "HSPG2", "LAMB3", "UBA1", "CEACAM5", "MYH9", "PDIA3", "CTSC", "SERPINB4", "MUC5AC", "GSN")

GALA_expr_norm_transposed_up <- t(GALA_expr_norm)
GALA_expr_norm_transposed_up <- GALA_expr_norm_transposed_up[order(match(rownames(GALA_expr_norm_transposed_up), rownames(GALA_phenotype))), ]
GALA_expr_norm_transposed_up <- log(as.matrix(GALA_expr_norm_transposed_up)+1)
GALA_expr_norm_transposed_up <- as.data.frame(GALA_expr_norm_transposed_up)
GALA_expr_norm_transposed_up <- GALA_expr_norm_transposed_up[, upregulated_apical_secretome_DEPs]
GALA_expr_norm_transposed_up$means <- rowMeans(GALA_expr_norm_transposed_up)


downregulated_apical_secretome_DEPs <- c("TMC5", "PROM1", "LTF", "ALDH3B1", "BPIFA1", "MYO1B", "TMC4", "ACTG2", "KRT13", "KLK11", "SLC34A2", 
	"SLC44A4", "MUC16", "CLIC6", "HIST1H2BH", "MSMB", "GNB2", "MYOF", "YWHAH", "FABP5", "H2AFY", "C4A", "EHD4", "RARRES1", "CFH", "BAIAP2", 
	"BASP1", "CP", "YWHAQ", "CIB1", "HIST1H2AH", "ALDH3A1", "LRG1", "GOLM1", "MUC5B", "C3", "TPPP3", "CAPZA1", "CLU", "CD59", "SLC9A3R1", "KRT2", 
	"FN1", "PPIB", "ACTC1", "SERPINA3", "MUC4", "SCGB1A1", "SERPINF1", "KRT1", "STOM", "SERPINA1", "IGFBP2", "LGALS3BP", "CTSB", "TAGLN2", "CTSD", "EZR")

GALA_expr_norm_transposed_down <- t(GALA_expr_norm)
GALA_expr_norm_transposed_down <- GALA_expr_norm_transposed_down[ order(match(rownames(GALA_expr_norm_transposed_down), rownames(GALA_phenotype))), ]
GALA_expr_norm_transposed_down <- log(as.matrix(GALA_expr_norm_transposed_down)+1)
GALA_expr_norm_transposed_down <- as.data.frame(GALA_expr_norm_transposed_down)
GALA_expr_norm_transposed_down <- GALA_expr_norm_transposed_down[, downregulated_apical_secretome_DEPs]
GALA_expr_norm_transposed_down$means <- rowMeans(GALA_expr_norm_transposed_down)

#pdf("GALA_apical_secretome_DEPs_boxplots.pdf")
dev.new(width = 6, height = 5)
par(mfrow=c(2,3),bty='n')
boxplot(GALA_expr_norm_transposed_up$means~GALA_phenotype$type2_status, names = c("Type 2 High", "Type 2 Low"), plot = TRUE, 
	main = "Upregulated Apical Secretome DEPs", ylab = "Mean Normalized Expression", cex.lab = 0.75, cex.axis = 0.5, 
	boxlwd = 0.25, lwd = 0.5, lwd.axis = 0.5, outline = FALSE, las = 1, cex.main = 1)
beeswarm(GALA_expr_norm_transposed_up$means~GALA_phenotype$type2_status, col = c("red", "blue"), 
	pch = 16, cex = 0.35, add = TRUE, outline = FALSE)
mtext(paste("p = ",formatC(wilcox.test(GALA_expr_norm_transposed_up$means[which(GALA_phenotype$type2_status == "type2_low")],
	GALA_expr_norm_transposed_up$means[which(GALA_phenotype$type2_status == "type2_high")],alternative="less")$p.value,
	format="e",digits=2),sep=""),side=3,line=0.2,at=2.5,adj=1,cex=1)

boxplot(GALA_expr_norm_transposed_down$means~GALA_phenotype$type2_status, names = c("Type 2 High", "Type 2 Low"), plot = TRUE, 
	main = "Downregulated Apical Secretome DEPs", ylab = "Mean Normalized Expression", cex.lab = 0.75, cex.axis = 0.5, 
	boxlwd = 0.25, lwd = 0.5, lwd.axis = 0.5, outline = FALSE, las = 1, cex.main = 1)
beeswarm(GALA_expr_norm_transposed_down$means~GALA_phenotype$type2_status, col = c("red", "blue"), pch = 16, cex = 0.35, 
	add = TRUE, outline = FALSE)
mtext(paste("p = ",formatC(wilcox.test(GALA_expr_norm_transposed_down$means[which(GALA_phenotype$type2_status == "type2_low")],
	GALA_expr_norm_transposed_down$means[which(GALA_phenotype$type2_status == "type2_high")],alternative="greater")$p.value,
	format="e",digits=2),sep=""),side=3,line=0.2,at=2.5,adj=1,cex=1)
#dev.off()











### 3. GALA expression heat map

#First do differential expression between Type 2 high and Type 2 low samples
all(colnames(GALA_raw) == rownames(GALA_phenotype))

#Specify design
subject <- colnames(GALA_raw)
treatment <- GALA_phenotype$type2_status
design <- data.frame(row.names = colnames(GALA_raw), subject = subject, treatment = treatment)

########Now do DESeq2
dds <- DESeqDataSetFromMatrix(
	countData = GALA_raw,
	colData = design,
	design = ~treatment) 

#Set control to reference
dds$treatment <- relevel(dds$treatment, "type2_low")

#Do DE
dds <- DESeq(dds)
res <- results(dds)

#Prepare table for printing
res <- cbind("gene"=rownames(res),as.data.frame(res))
res <- res[order(res$padj),]
res <- res[which(!is.na(res$padj)),] #toss genes with mean norm counts below default specified threshold

#Write DE results to text file
write.table(res,file = "GALA_DEGs_DESeq2.txt",quote = F, sep = "\t",row.names=F)

#Significant GALA DEGs
DEGs_gala_up<-res[which(res$log2FoldChange > 0 & res$padj < 0.05),]
DEGs_gala_down<-res[which(res$log2FoldChange < 0 & res$padj < 0.05),]


#Now get gene lists from transcriptomics to plot
type_2_inflam_genes <- c("CPA3", "MS4A2", "TPSB2", "TPSAB1", "IL1RL1", "POSTN", "CLCA1", "CST1", "DPP4", "AGR2")
transporters_atp<-c("ATP1A1","ATP1B1","ATP13A5","ATP6V0E1","ATP7B")
transporters_slc<-c("SLC5A8","SLC9A9","SLC9B2","SLC12A2","SLC22A23","SLC31A1","SLC37A1","SLC38A6","SLC39A7","SLC39A8","SLCO3A1","SLCO4C1")
endopeptidase_inhibitors<-c("CAPN14","CSTB","FETUB","SERPINB2","SERPINB6","SERPINB10","SPINK5","CSTA","TIMP1","SERPINB4","SERPINB3")
mucin_associated<-c("FCGBP","GSN","ITLN1","MUC5AC","SCIN")
glycosylation<-c("ADAMTSL3","B3GNT6","B3GNT2","CHST6","GALK2","GALNT1","GCNT3")
pro_inflammatory<-c("ALOX15","ANXA3","CCL26","CD36","DDX58","LGALS7","NOS2","RIPK2","SH2D1B")
d_reductases<-c("AKR1C2","AKR1C3","DHRS3","DHRS9","GLUD1","RDH10","STEAP1B","STEAP4")
d_innate_detox<-c("ABCA13","CP","CYP2B7P","GLUL","GSTA1","GSTA2","LYN","SCGB1A1","SCGB3A1","WFDC2")
d_innate_recruit<-c("CD74","CXCL1","CXCL17","IL18","S100A2","S100A4","S100A6","S100A8","S100A9")
cilia_genes_use<-c("ARMC4","CC2D2A","CFAP46","CEP126","DNAH5","DNAH7","DYNC2H1","IFT57","IFT81","KIF3A","KIF21A","RFX3","RFX2","TP73","TUBA1A","WDR35")

#Isolate genes that are significant in GALA
type_2_inflam_genes<-type_2_inflam_genes[which(type_2_inflam_genes %in% DEGs_gala_up$gene)]
transporters_atp<-transporters_atp[which(transporters_atp %in% DEGs_gala_up$gene)]
transporters_slc<-transporters_slc[which(transporters_slc %in% DEGs_gala_up$gene)]
endopeptidase_inhibitors<-endopeptidase_inhibitors[which(endopeptidase_inhibitors %in% DEGs_gala_up$gene)]
mucin_associated<-mucin_associated[which(mucin_associated %in% DEGs_gala_up$gene)]
glycosylation<-glycosylation[which(glycosylation %in% DEGs_gala_up$gene)]
pro_inflammatory<-pro_inflammatory[which(pro_inflammatory %in% DEGs_gala_up$gene)]
d_reductases<-d_reductases[which(d_reductases %in% DEGs_gala_down$gene)]
d_innate_detox<-d_innate_detox[which(d_innate_detox %in% DEGs_gala_down$gene)]
d_innate_recruit<-d_innate_recruit[which(d_innate_recruit %in% DEGs_gala_down$gene)]
cilia_genes_use<-cilia_genes_use[which(cilia_genes_use %in% DEGs_gala_down$gene)]

##getting the colors and ordering for the heatmap ready
rownames(GALA_phenotype) <- GALA_phenotype$SubjectID
GALA_phenotype <- GALA_phenotype[colnames(GALA_expr_norm),]

asthma_type_cols <- c("red", "blue") #red for type 2 high subtype, blue for type 2 low subtype
gene_cols <- c("orange", "pink", "yellow3", "purple")
gene_cols <- rev(gene_cols)

##ordering the genes (rows of input matrix)
gene_order <- c(type_2_inflam_genes, transporters_atp, transporters_slc, endopeptidase_inhibitors,mucin_associated,
	glycosylation,pro_inflammatory,d_reductases,d_innate_detox,d_innate_recruit,cilia_genes_use)
gene_order <- rev(gene_order)

gene_groups <-c(rep("orange",length(type_2_inflam_genes)),
rep("pink",length(transporters_atp)),
rep("yellow",length(transporters_slc)),
rep("purple",length(endopeptidase_inhibitors)),
rep("orange",length(mucin_associated)),
rep("pink",length(glycosylation)),
rep("yellow",length(pro_inflammatory)),
rep("purple",length(d_reductases)),
rep("orange",length(d_innate_detox)),
rep("pink",length(d_innate_recruit)),
rep("yellow",length(cilia_genes_use)))
gene_groups <- rev(gene_groups)

##ordering the cells ("columns" of input matrix) - get type 2 low and sort by eigengene and then get type 2 high and sort by eigengenes
cells_order <- c(which(GALA_phenotype[order(GALA_phenotype$MEskyblue),]$type2_status == "type2_low"), 
	which(GALA_phenotype[order(GALA_phenotype$MEskyblue),]$type2_status == "type2_high"))
#cells_order <- rownames(GALA_phenotype[order(GALA_phenotype$MEskyblue),]) #If just want to sort by eigengene, but not bin into T2L and H first

##making the heatmap matrix
mat <- GALA_vst[,order(GALA_phenotype$MEskyblue)][gene_order, cells_order]

##getting colors for type2_high vs. type2_low designations
type <- color.scale(as.numeric(as.factor(GALA_phenotype[order(GALA_phenotype$MEskyblue),][cells_order,]$type2_status)), 
	extremes = asthma_type_cols, color.spec = "rgb")
#note that the continuous line is merely conceptual as it doesn't match the true order of eigengenes (which should be ordered by
#GALA_phenotype[order(GALA_phenotype$MEskyblue),]
type_continuous <- color.scale(as.numeric(as.factor(GALA_phenotype[order(GALA_phenotype$MEskyblue),][cells_order,]$MEskyblue)), 
	extremes = asthma_type_cols[c(2,1)], color.spec = "rgb")
ColSideColors<-data.frame(type,type_continuous)

##getting colors for gene groups
gene_group <- color.scale(as.numeric(as.factor(gene_groups)), extremes=sort(unique(gene_groups)), color.spec="rgb")
RowSideColors <- data.frame(gene_group=gene_group)

###make the heatmap
mycols = colorRampPalette(c("blue","white","red"))(1000)
breakscale = c(-8,seq(-1.3,1.3, length.out=length(mycols)-1), 8)

pdf("GALA_heatmap_figure7.pdf",height=24,width=18)
heatmap3(mat, scale="row",Rowv=NA,Colv=NA,col=mycols, ColSideColors=as.matrix(ColSideColors), 
	RowSideColors=as.matrix(RowSideColors), cexRow=0.7,cexCol=2,breaks=breakscale,labCol="",useRaster=F)
dev.off()









