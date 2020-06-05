#########################################
#############Monocle2 analysis#############
#########################################

library(data.table)
library(heatmap3)
library(plotrix)
library(devtools)
library(monocle)
library(Seurat)
library(dplyr)
library(Matrix)

####################### Get data
#Bring in Seurat object for the acute dataset
attach("../data/IL13_acute_seurat.rda")

#Create datasets with raw expression
#Mucus secretory populations
exp_muc<-as.data.frame(IL13s@raw.data[,which(IL13s@meta.data$clusters_9 == "c5"| IL13s@meta.data$clusters_9 == "c4")])
#Defense secretory populations
exp_def<-as.data.frame(IL13s@raw.data[,which(IL13s@meta.data$clusters_9 == "c1"| IL13s@meta.data$clusters_9 == "c0")])



######################## Get pseudotime trajectories for each population
 
################################################
################### exp_muc ####################
################################################

#######################Read in attributes table for the cells (based on Seurat metadata) and read dataset into Monocle
currMetadata<-IL13s@meta.data[which(IL13s@meta.data$clusters_9 == "c5"| IL13s@meta.data$clusters_9 == "c4"),]
sample_sheet = data.frame("nGene"=currMetadata$nGene,"nUMI"=currMetadata$nUMI,"treatment"=currMetadata$treatment,
	"clusters_9"=as.character(currMetadata$clusters_9),stringsAsFactors=F,row.names=rownames(currMetadata))
gene_names = data.frame("gene_short_name"=rownames(exp_muc),stringsAsFactors=F,row.names=rownames(exp_muc))
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_names) 
HSMM_muc<-newCellDataSet(as.matrix(exp_muc),phenoData = pd,featureData = fd, expressionFamily=negbinomial())

##Means and dispersions	
HSMM_muc <- estimateSizeFactors(HSMM_muc)
HSMM_muc <- estimateDispersions(HSMM_muc)

##Toss genes not expressed in at least 1% of cells
HSMM_muc <- detectGenes(HSMM_muc, min_expr = 0.1)
expressed_genes<-row.names(subset(fData(HSMM_muc), num_cells_expressed >= 0.01 * ncol(exp_muc)))
HSMM_muc<-HSMM_muc[expressed_genes,]




####################### Get genes for ordering
#Get differentially expressed genes between the two clusters
clustering_DEG_genes_muc <- differentialGeneTest(HSMM_muc,fullModelFormulaStr = '~clusters_9',cores = detectCores())




#######################Pseudotime plotting
#Make trajectory using the top 500 DEGs
ordering_genes_muc <-row.names(clustering_DEG_genes_muc)[order(clustering_DEG_genes_muc$qval)][1:500]
HSMM_muc <- setOrderingFilter(HSMM_muc, ordering_genes_muc)
HSMM_muc <- reduceDimension(HSMM_muc, max_components=2,reduction_method="DDRTree")
HSMM_muc <- orderCells(HSMM_muc)

#MeanVarPlot
#pdf("IL13s_meanVarPlotForOrdering_monocle.pdf")
plot_ordering_genes(HSMM_muc)
#dev.off()

#Plot trajectories
pdf('cellTrajectory_muc_500DEGs_monocle.pdf')
p1<-plot_cell_trajectory(HSMM_muc, color_by = "clusters_9",show_branch_points = F,cell_size=2)+ 
	scale_colour_manual(values=c("darkgreen","purple"))
	
p2<-plot_cell_trajectory(HSMM_muc, color_by="Pseudotime",show_branch_points = F,cell_size=2) + 
	scale_colour_gradient(low="black",high="lightblue")

p3<-plot_cell_trajectory(HSMM_muc, color_by="treatment",show_branch_points = FALSE,cell_size=2) + 
	scale_colour_manual(values=c("grey","black"))

p4<-plot_cell_trajectory(HSMM_muc, color_by="State",show_branch_points = FALSE,cell_size=2) + 
	scale_colour_manual(values=c("red3","cornflowerblue","orange","turquoise"))

plot_grid(p1,p2,p3,p4)
dev.off()




#Test and plot for pseudotime dependence of the top variable genes in the dataset
#First, subset genes based on expression mean and dispersion for each #2,116 genes
disp_table_muc <- dispersionTable(HSMM_muc)
unsup_clustering_genes_muc <- subset(disp_table_muc, mean_expression >= 0.01 & 
	dispersion_empirical >= 0.01 * dispersion_fit)$gene_id
HSMM_muc <- setOrderingFilter(HSMM_muc, unsup_clustering_genes_muc)
pdf("plot_ordering_genes_var_muc.pdf")
plot_ordering_genes(HSMM_muc)
dev.off()
length(unsup_clustering_genes_muc)

to_be_tested <- row.names(subset(fData(HSMM_muc),gene_short_name %in% unsup_clustering_genes_muc))
to_be_tested <- to_be_tested[-grep("^ENSG",to_be_tested)]
cds_subset <- HSMM_muc[to_be_tested,]

diff_test_var_muc <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = detectCores())
diff_test_var_muc<-diff_test_var_muc[order(diff_test_var_muc$qval),]

diff_test_var_muc[,c("gene_short_name", "pval", "qval")]

pdf("Genes_bestPseudotimeGenes_topVarGenes_muc_sig.pdf")
ncount<-1
for(i in 1:(length(rownames(diff_test_var_muc)[1:450])) / 15){
	cds_subset <- HSMM_muc[rownames(diff_test_var_muc)[1:450][ncount:(ncount + 14)],]
	print(plot_genes_in_pseudotime(cds_subset, color_by = "treatment",nrow=5,ncol=3))
	ncount<-ncount + 15
}
dev.off()



###### Now make heat map of these genes
#Function for calculating binned values for monocle
getBinnedValues<-function(pData,metadata_col_index,quantile=F,k=5){
	#Sort metadata by pseudotime
	pData_byPseud<-pData[order(pData$Pseudotime),]
	#Establish dataframe
	metadataTab<-data.frame(matrix(nrow=100,ncol=2))
	colnames(metadataTab)<-c(colnames(pData_byPseud)[metadata_col_index],"pseudotime")
	
	#Create bins for categorical variable
	d<-pData_byPseud[,metadata_col_index]
	p<-pData_byPseud$Pseudotime
	bins<-round(seq(1,length(d),by=(length(d) / 100)))
	bins_p<-seq(min(p),max(p),length.out = 100)
	for(i in 1:100){
		#For each bin, get pseudotime midpoint
		if(quantile){ #if you want to create bins based on the data (equal number of cells per bin)
			metadataTab[i,2]<-mean(range(p[bins[i]:(round(bins[i] + (length(d) / 100)) - 1)]))
		}else{ #if you want to bin the way Monocle does, i.e., bins at equal intervals across pseudotime and takes the midpoint
		metadataTab[i,2]<-mean(c(bins_p[i],(bins_p[i] + (bins_p[2] - bins_p[1]))))
		}
		#For each bin, get dominant categorial value
		if(quantile){#if you want to create bins based on the data (equal number of cells per bin)
			cellTab<-sort(table(d[bins[i]:(round(bins[i] + (length(d) / 100)) - 1)]),decreasing = T)
			metadataTab[i,1]<-sample(names(cellTab)[which(cellTab == max(cellTab))],1)
			#if you want to bin the way Monocle does, i.e., bins at equal intervals across pseudotime and takes the 
			#most representative of k nearest neighbor cells
		}else{
			cellTab<-sort(table(pData_byPseud[order(abs(pData_byPseud$Pseudotime - metadataTab[i,2])),metadata_col_index][1:k]))
			metadataTab[i,1]<-sample(names(cellTab)[which(cellTab == max(cellTab))],1)
		}
	}
	return(metadataTab)
}

#Plot a heat map, where rows are genes and where columns are cells, ordered
#by pseudotime
genes_subset<-rownames(diff_test_var_muc)[which(diff_test_var_muc$qval < 0.05)]
annotations<-getBinnedValues(pData=pData(HSMM_muc),metadata_col_index=4)
pdf("Heatmap_bestPseudotimeGenes_topVarGenes_muc_sig.pdf")
plot_pseudotime_heatmap(HSMM_muc[genes_subset,],
	num_clusters = 6,
	cores = detectCores(),
	show_rownames = T,
	add_annotation_col=annotations)
dev.off()



##Export table of significant genes
write.table(diff_test_var_muc[which(diff_test_var_muc$qval < 0.05),c(5,3:4,6)],row.names=F,sep="\t",quote=F,file="BestPseudotimeGenes_topVarGenes_muc_sig.txt")



















################################################
################### exp_def ####################
################################################

#######################Read in attributes table for the cells (based on Seurat metadata) and read dataset into Monocle
currMetadata<-IL13s@meta.data[which(IL13s@meta.data$clusters_9 == "c1"| IL13s@meta.data$clusters_9 == "c0"),]
sample_sheet = data.frame("nGene"=currMetadata$nGene,"nUMI"=currMetadata$nUMI,"treatment"=currMetadata$treatment,
	"clusters_9"=as.character(currMetadata$clusters_9),stringsAsFactors=F,row.names=rownames(currMetadata))
gene_names = data.frame("gene_short_name"=rownames(expression_mat),stringsAsFactors=F,row.names=rownames(exp_def))
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_names) 
HSMM_def<-newCellDataSet(as.matrix(exp_def),phenoData = pd,featureData = fd,
	expressionFamily=negbinomial())

##Means and dispersions	
HSMM_def <- estimateSizeFactors(HSMM_def)
HSMM_def <- estimateDispersions(HSMM_def)

##Toss genes not expressed in at least 1% of cells
HSMM_def <- detectGenes(HSMM_def, min_expr = 0.1)
expressed_genes<-row.names(subset(fData(HSMM_def), num_cells_expressed >= 0.01 * ncol(exp_def)))
HSMM_def<-HSMM_def[expressed_genes,]




####################### Get genes for ordering
#Get differentially expressed genes between the two clusters
clustering_DEG_genes_def <- differentialGeneTest(HSMM_def,fullModelFormulaStr = '~clusters_9',cores = detectCores())




#######################Pseudotime plotting
#Make trajectory using the top 500 DEGs
ordering_genes_def <-row.names(clustering_DEG_genes_def)[order(clustering_DEG_genes_def$qval)][1:500]
HSMM_def <- setOrderingFilter(HSMM_def, ordering_genes_def)
HSMM_def <- reduceDimension(HSMM_def, max_components=2,reduction_method="DDRTree")
HSMM_def <- orderCells(HSMM_def,reverse=T)

#MeanVarPlot
#pdf("IL13s_meanVarPlotForOrdering_monocle.pdf")
plot_ordering_genes(HSMM_def)
#dev.off()

#Plot trajectories
pdf('cellTrajectory_def_500DEGs_monocle.pdf')
p1<-plot_cell_trajectory(HSMM_def, color_by = "clusters_9",show_branch_points = F,cell_size=2)+ 
	scale_colour_manual(values=c("saddlebrown","yellow2"))
	
p2<-plot_cell_trajectory(HSMM_def, color_by="Pseudotime",show_branch_points = F,cell_size=2) + 
	scale_colour_gradient(low="black",high="lightblue")

p3<-plot_cell_trajectory(HSMM_def, color_by="treatment",show_branch_points = FALSE,cell_size=2) + 
	scale_colour_manual(values=c("grey","black"))

p4<-plot_cell_trajectory(HSMM_def, color_by="State",show_branch_points = FALSE,cell_size=2) + 
	scale_colour_manual(values=c("red3","cornflowerblue","orange","turquoise","midnightblue","violet","green3"))

plot_grid(p1,p2,p3,p4)
dev.off()




#Test and plot for pseudotime dependence of the top variable genes in the dataset
#First, subset genes based on expression mean and dispersion
disp_table_def <- dispersionTable(HSMM_def)
unsup_clustering_genes_def <- subset(disp_table_def, mean_expression >= 0.01 & 
	dispersion_empirical >= 0.01 * dispersion_fit)$gene_id
HSMM_def <- setOrderingFilter(HSMM_def, unsup_clustering_genes_def)
pdf("plot_ordering_genes_var_def.pdf")
plot_ordering_genes(HSMM_def)
dev.off()
length(unsup_clustering_genes_def)

to_be_tested <- row.names(subset(fData(HSMM_def),gene_short_name %in% unsup_clustering_genes_def))
to_be_tested <- to_be_tested[-grep("^ENSG",to_be_tested)]
cds_subset <- HSMM_def[to_be_tested,]

diff_test_var_def <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = detectCores())
diff_test_var_def<-diff_test_var_def[order(diff_test_var_def$qval),]

diff_test_var_def[,c("gene_short_name", "pval", "qval")]

pdf("Genes_bestPseudotimeGenes_topVarGenes_def_sig.pdf")
ncount<-1
for(i in 1:(length(rownames(diff_test_var_def)[1:420]) / 15)){
	cds_subset <- HSMM_def[rownames(diff_test_var_def)[1:420][ncount:(ncount + 14)],]
	print(plot_genes_in_pseudotime(cds_subset, color_by = "treatment",nrow=5,ncol=3))
	ncount<-ncount + 15
}
dev.off()



#Now make heat map of these genes
genes_subset<-rownames(diff_test_var_def)[which(diff_test_var_def$qval < 0.05)]
annotations<-getBinnedValues(pData=pData(HSMM_def),metadata_col_index=4)

#Plot a heat map, where rows are genes and where columns are cells, ordered
#by pseudotime
pdf("Heatmap_bestPseudotimeGenes_topVarGenes_def_sig.pdf")
plot_pseudotime_heatmap(HSMM_def[genes_subset,],
	num_clusters = 6,
	cores = detectCores(),
	show_rownames = T,
	add_annotation_col=annotations)
dev.off()



##Export table of significant genes
write.table(diff_test_var_def[which(diff_test_var_def$qval < 0.05),c(5,3:4,6)],
	row.names=F,sep="\t",quote=F,file="BestPseudotimeGenes_topVarGenes_def_sig.txt")


















############################ Plot pseudotime for transcription factors identified using IPA for both mucus and defense datasets
TFs_muc<-read.table("../data/TranscriptionFactors_pseudotimeGenes_muc.txt",header=T,stringsAsFactors=F,sep="\t")
TFs_def<-read.table("../data/TranscriptionFactors/IPA/TranscriptionFactors_pseudotimeGenes_def.txt",header=T,stringsAsFactors=F,sep="\t")

TFs_muc<-diff_test_var_muc[which(rownames(diff_test_var_muc[which(diff_test_var_muc$qval < 0.05),]) %in% 
	TFs_muc$ID[which(TFs_muc$Type == "transcription regulator")]),c(4,6,5)]
TFs_def<-diff_test_var_def[which(rownames(diff_test_var_def[which(diff_test_var_def$qval < 0.05),]) %in% 
	TFs_def$ID[which(TFs_def$Type == "transcription regulator")]),c(4,6,5)]

###### DEGs shared by the two secretory cell types
TFs_mucANDdef<-TFs_muc[which(rownames(TFs_muc) %in% rownames(TFs_def)),]

#pdf("Genes_TFs_shared_def.pdf")
dev.new(height=5,width=7)
cds_subset <- HSMM_def[rownames(TFs_mucANDdef),]
print(plot_genes_in_pseudotime(cds_subset, color_by = "treatment",nrow=4,ncol=4))

#pdf("Genes_TFs_shared_muc.pdf")
dev.new(height=5,width=7)
cds_subset <- HSMM_muc[rownames(TFs_mucANDdef),]
print(plot_genes_in_pseudotime(cds_subset, color_by = "treatment",nrow=4,ncol=4))


######## Unique DEGs in one or the other cell type
TFs_mucUnique<-TFs_def[-which(rownames(TFs_muc) %in% rownames(TFs_def)),]
TFs_defUnique<-TFs_def[-which(rownames(TFs_def) %in% rownames(TFs_muc)),]

#pdf("Genes_TFs_unique_def.pdf")
dev.new(height=5,width=7)
cds_subset <- HSMM_def[rownames(TFs_defUnique),]
print(plot_genes_in_pseudotime(cds_subset, color_by = "treatment",nrow=4,ncol=4))

#pdf("Genes_TFs_unique2_muc.pdf")
dev.new(height=5,width=7)
cds_subset <- HSMM_muc[rownames(TFs_mucUnique)[17:31],]
print(plot_genes_in_pseudotime(cds_subset, color_by = "treatment",nrow=4,ncol=4))


















######################## Finally, plot smooth curved expression across pseudotime for overlaid TFs

#This function takes the plot_genes_in_pseudotime and gives it more options that I want
#Most importantly, it takes a list of datasets (cds_subsets) so that points and curves can be overlaid onto the same plot.
#I can also take a set of colors for points from each dataset and category 
#	(there need to be as many colors as there are number of categories to color x number of datasets.
#Print_points specifies whether or not you want points printed at all
#Curve_cols gives colors used for each curve (in order of the datasets in the list).
#Fix_y gives whether or not the y-axis should be fixed across all genes or allowed to adjust to the expression of each gene
#overlay gives whether you want to overlay the same genes from different datasets ("datasets") or different genes from the same dataset ("genes") - in the case
	#of different genes from the same dataset, you can only plot one dataset at a time. So cds_subsets should just be a list of the same dataset with different genes subsetted
	#that will all be overlaid onto a single plot. Also, I don't really have print_points = TRUE working for this case. Finally, note that right now, because there
	#are such vast differences in average expression across genes, all curves are scaled to be between 0 and 1.
prepPseudotimePlotDatasets<-function(cds_subsets,color_by = NULL,color_vec = NULL,nrow = 4,ncol = 4,print_points=TRUE,
	curve_cols = NULL,fix_y = T,overlay="datasets",min_expr = NULL,cell_size = 0.75,panel_order = NULL,trend_formula = "~ sm.ns(Pseudotime, df=3)",
	label_by_short_name = TRUE,relative_expr = FALSE,vertical_jitter = NULL, horizontal_jitter = NULL){
	
	#Set of up list of processed datasets
	cds_expr_final<-list()
	
	#For each dataset to overlay, get ready for plotting
	for(i in 1:length(cds_subsets)){
		cds_subset<-cds_subsets[[i]]
    	f_id <- NA
    	Cell <- NA
    	if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
    	    "negbinomial.size")) {
    	    integer_expression <- TRUE
    	}else{
    	    integer_expression <- FALSE
    	    relative_expr <- TRUE
    	}
    	
    	#Get size factor adjusted expression
    	if (integer_expression) {
    	    cds_exprs <- exprs(cds_subset)
    	    if (relative_expr) {
    	        if (is.null(sizeFactors(cds_subset))) {
    	            stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
    	        }
    	        cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    	    }
    	    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
    	}else{
    	    cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
    	}
    	
    	#Set min expression to the lowest detection limit
    	if (is.null(min_expr)) {
    	    min_expr <- cds_subset@lowerDetectionLimit
    	}
    	
    	#Add in meta.data
    	colnames(cds_exprs) <- c("f_id", "Cell", "expression")
    	cds_pData <- pData(cds_subset)
    	cds_fData <- fData(cds_subset)
    	cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
    	cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
    	#Set adjusted_expression
    	if (integer_expression) {
    	    cds_exprs$adjusted_expression <- cds_exprs$expression
    	}else{
    	    cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
    	}
    	#Set feature_label
    	if (label_by_short_name == TRUE) {
    	    if (is.null(cds_exprs$gene_short_name) == FALSE) {
    	        cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
    	        cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    	    }else{
    	        cds_exprs$feature_label <- cds_exprs$f_id
    	    }
    	}else{
    	    cds_exprs$feature_label <- cds_exprs$f_id
    	}
    	cds_expr_final[[length(cds_expr_final) + 1]]<-cds_exprs
    }
    	    	
    #Get smooth curve for expression of each gene across pseudotime (but note cells not yet ordered by psudotime)
    #and merge these with the expression + metadata table	
    for(i in 1:length(cds_subsets)){
    	cds_subset<-cds_subsets[[i]]
    	cds_expr_final[[i]]$f_id <- as.character(cds_expr_final[[i]]$f_id) 
    	cds_expr_final[[i]]$feature_label <- factor(cds_expr_final[[i]]$feature_label)
    	new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
    	model_expectation <- genSmoothCurves(cds_subset, cores = 1, 
    	    trend_formula = trend_formula, relative_expr = relative_expr, new_data = new_data)
    	colnames(model_expectation) <- colnames(cds_subset)
    	expectation <- ddply(cds_expr_final[[i]], .(f_id, Cell), function(x) data.frame(expectation = model_expectation[x$f_id, 
    	    x$Cell]))
    	
    	#If overlaying genes, move the curves to be between 0 and 1
    	if(overlay == "genes"){
    		expectation$expectation<-(expectation$expectation - min(expectation$expectation)) / 
    			(max(expectation$expectation) - min(expectation$expectation))
    	}
    	cds_expr_final[[i]] <- merge(cds_expr_final[[i]], expectation)
    	
    	#Set expression below the specified minimum expression value to the minimum expression value
    	cds_expr_final[[i]]$expression[cds_expr_final[[i]]$expression < min_expr] <- min_expr
    	
    	#Do the same for the expected expression value
    	cds_expr_final[[i]]$expectation[cds_expr_final[[i]]$expectation < min_expr] <- min_expr
    	if (is.null(panel_order) == FALSE) {
    	    cds_expr_final[[i]]$feature_label <- factor(cds_expr_final[[i]]$feature_label, 
    	        levels = panel_order)
    	}
    }
 

    #########Now that we have the list of overlaid datasets prepped, plot each of them using this:
    #If there is more than one dataset, make the categorical variable specific to each dataset
    for(i in 1:length(cds_subsets)){
    	categorical_index<-grep(color_by,colnames(cds_expr_final[[i]]))
    	cds_expr_final[[i]][,categorical_index]<-paste(cds_expr_final[[i]][,categorical_index],"_",i,sep="")
    }
    
    #Then establish an empty plot with y axis fitting the global expression range (across all genes) and with the
    #x axis fitting the pseudotime range (across all datasets)
    appendedData<-data.frame()
    for(i in 1:length(cds_subsets)){
    	appendedData<-rbind(appendedData,cds_expr_final[[i]])
    }
 
    #Now make a ggplot
    q <- ggplot(aes(Pseudotime, expression), data = appendedData)

	#Make scatter plot of expression against pseudotime across all genes 
	#cell_size gives the size of the points
	#color_by gives the categories to color by
	#I've added point color specification here
	if(print_points){
	    if(is.null(color_by) == FALSE){
		     q <- q + geom_point(aes_string(color = color_by), size = I(cell_size), alpha = 0.3, shape = 16,
		        position = position_jitter(horizontal_jitter, vertical_jitter),data = appendedData) +
		        scale_color_manual(values=color_vec[order(unique(appendedData[,categorical_index]))])
		}else{
		    q <- q + geom_point(size = I(cell_size), position = position_jitter(horizontal_jitter, 
	    	  	vertical_jitter),data = appendedData)
	    }
	}
    
    #Plot the expression expectation curves across the combined plot, one for each dataset
 	for(i in 1:length(cds_expr_final)){
	    q <- q + geom_line(aes(x = Pseudotime, y = expectation), color = curve_cols[i], #size = 1,
        	data = cds_expr_final[[i]])
    }
        
    #Now break up the plots based on gene (if overlaying onto different datasets)
    #scale_y_log10 places expression on a log scale
    #scales = "free_y" (default) lets the y-axis shift to fit the data
    #scales = "fixed" keeps the x and y axes the same across plots
    if(fix_y){
    	if(overlay == "datasets"){
 		   	q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow, 
        		ncol = ncol, scales = "fixed")
        }
    }else{
    	if(overlay == "datasets"){
 		   	q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow, 
        		ncol = ncol, scales = "free_y")
        }
    }
    
    if (min_expr < 1) {
        q <- q + expand_limits(y = c(min_expr, 1))
    }
    
    #Add labels
    if (relative_expr) {
        q <- q + ylab("Relative Expression")
    }else{
        q <- q + ylab("Absolute Expression")
    }
    if(overlay == "genes"){
    	q <- q + ylab("Scaled Expression")
    }
    q <- q + xlab("Pseudotime")
 
    #Get rid of ggplot background stuff
    q <- q + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		panel.background = element_blank(), axis.line = element_line(colour = "black"))
	q
}


######### Make separate plots for each gene in the paper
library(plyr)

#Upregulated genes
cds_subsets<-list(HSMM_def[rownames(TFs_mucANDdef)[c(15,16,7,12,9)],],HSMM_muc[rownames(TFs_mucANDdef)[c(15,16,7,12,9)],])
#pdf("Genes_TFs_shared_up.pdf")
dev.new(height = 1.5, width = 6)
prepPseudotimePlotDatasets(cds_subsets=cds_subsets,color_by = "treatment",
	color_vec = c("yellow3","saddlebrown","purple","darkgreen"),nrow = 1,ncol = 5,
	curve_cols=c("black","red"),fix_y = T,cell_size = 1.5)

#Downregulated genes
cds_subsets<-list(HSMM_def[rownames(TFs_mucANDdef)[c(8,1,6,5,3)],],HSMM_muc[rownames(TFs_mucANDdef)[c(8,1,6,5,3)],])
#pdf("Genes_TFs_shared_down.pdf")
dev.new(height = 1.5, width = 6)
prepPseudotimePlotDatasets(cds_subsets=cds_subsets,color_by = "treatment",
	color_vec = c("yellow2","saddlebrown","purple","darkgreen"),nrow = 1,ncol = 5,
	curve_cols=c("black","red"),fix_y = T,cell_size = 1.5)






