############# integrate all HC/PFC samples

library(Seurat)
library(harmony)
library(knitr)
library(kableExtra)
library(ggplot2)
library(cowplot)
library(plyr)
library(dplyr)
library(scater)

### 1.1 Load the Cell Ranger Matrix Data and create the base Seurat object.
	options(stringsAsFactors = FALSE)
	setwd("/projects/compsci/zhangjg/scRNA_test")
	dataset_loc <- "./data/"
	result_loc <- "./results/"
	tissue <- "HC"		## tissue: "HC" , "PFC"
	sample_info <- "sample_metrics_final.csv"
	
	All_Ids = read.csv(sample_info)
	Sample_ID = subset(All_Ids, Sample == tissue); 
	ids <- as.vector(Sample_ID$SCBID)
	rownames(Sample_ID)= 1:length(ids)	

	d10x.data <- sapply(ids, function(i){print(i)
	  d10x <- Read10X(file.path(dataset_loc,i,"/filtered_feature_bc_matrix"))
	  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")			### add the subject name to column names
	  d10x
	})

	experiment.data <- do.call("cbind", d10x.data)
	experiment.aggregate <- CreateSeuratObject(
	  experiment.data,
	  project = "scRNA workshop",
	  min.cells = 10,
	  min.features = 200,
	  names.field = 2,
	  names.delim = "\\-")
	  
	  experiment.aggregate$R_S <- experiment.aggregate$orig.ident
	  experiment.aggregate$geno <- experiment.aggregate$orig.ident
	  experiment.aggregate$stain <- experiment.aggregate$orig.ident
	  experiment.aggregate$age <- experiment.aggregate$orig.ident
	  
	  for (i in 1:length(ids)){print(i)	  
		experiment.aggregate$R_S = mapvalues(experiment.aggregate$R_S, c(as.vector(Sample_ID$SCBID[i])), c(as.vector(Sample_ID$R_S[i])))
		experiment.aggregate$geno = mapvalues(experiment.aggregate$geno,  c(as.vector(Sample_ID$SCBID[i])), c(as.vector(Sample_ID$Genotype[i])))
		experiment.aggregate$stain = mapvalues(experiment.aggregate$stain,  c(as.vector(Sample_ID$SCBID[i])), c(as.vector(Sample_ID$Line_Short[i])))
		experiment.aggregate$age = mapvalues(experiment.aggregate$age,  c(as.vector(Sample_ID$SCBID[i])), c(as.vector(Sample_ID$Age[i])))
	  }
	  
### 1.2 The percentage of reads that map to the mitochondrial genome and ribosomal genes
			
	experiment.aggregate$percent.mito <- PercentageFeatureSet(experiment.aggregate, pattern = "^mt-")
	
	mt_gene <- grep(pattern = "^mt-", x = rownames(experiment.aggregate), value = TRUE)
	rpS_gene <- grep(pattern = "^Rps", x = rownames(experiment.aggregate), value = TRUE)
	rpL_gene <- grep(pattern = "^Rpl", x = rownames(experiment.aggregate), value = TRUE)
	ps_gene <- grep(pattern = "-ps", x = rownames(experiment.aggregate), value = TRUE)
	
	rpS_gene <- setdiff(rpS_gene, ps_gene)
	rpL_gene <- setdiff(rpL_gene, ps_gene)

	experiment.aggregate$percent.rpS <- PercentageFeatureSet(experiment.aggregate, features = rpS_gene)  		## "^Rp[sl][[:digit:]]"
	experiment.aggregate$percent.rpL <- PercentageFeatureSet(experiment.aggregate, features = rpL_gene)
	experiment.aggregate$percent.ps <- PercentageFeatureSet(experiment.aggregate, features = ps_gene)
		
			
### 1.3 Cell filtering  (cell level filter)
	## Here we choose those that have percent mitochondrial genes max of 5% 
	## and unique UMI counts under 20,000 or greater than 500
	
	experiment.aggregate <- subset(experiment.aggregate, percent.mito <= 5)
	experiment.aggregate <- subset(experiment.aggregate, nCount_RNA >= 500 & nCount_RNA <= 20000)
	experiment.aggregate <- subset(experiment.aggregate, percent.rpL <= 5)
	experiment.aggregate <- subset(experiment.aggregate, percent.rpS <= 5)
	experiment.aggregate <- subset(experiment.aggregate, percent.ps <= 5)
	
	## and remove the MT , 

	mt_gene <- grep(pattern = "^mt-", x = rownames(experiment.aggregate), value = TRUE)
	Nomt_gene <- setdiff(rownames(experiment.aggregate), mt_gene)
	experiment.aggregate <- experiment.aggregate[Nomt_gene,]
	
	#Nomt_gene <- setdiff(rownames(experiment.aggregate), union(rpS_gene,rpL_gene,ps_gene))
	#experiment.aggregate <- experiment.aggregate[Nomt_gene,]
	
#### 2. “integrate”	
	
	### 2.1  Harmony

	experiment.aggregate <- NormalizeData(experiment.aggregate) %>% 
	  FindVariableFeatures() %>% 
	  ScaleData(vars.to.regress =c("nCount_RNA", "percent.rpL","percent.rpS","percent.ps")) %>% 
	  RunPCA(verbose = FALSE)
	experiment.aggregate <- RunHarmony(experiment.aggregate, group.by.vars = "orig.ident",plot_convergence = TRUE)
	#experiment.aggregate <- RunHarmony(experiment.aggregate, group.by.vars = "orig.ident",plot_convergence = TRUE,max.iter.harmony = 100)
	experiment.aggregate <- RunUMAP(experiment.aggregate, reduction = "harmony", dims = 1:30)
	experiment.aggregate <- RunTSNE(experiment.aggregate, reduction = "harmony", dims = 1:30)
	
    ###### Option!	
		#pdf(file = paste0(result_loc,"plot/","ElbowPlot.pdf"), width = 20, height = 20)
		#ElbowPlot(experiment.aggregate)	
		#dev.off()
		
	experiment.aggregate <- FindNeighbors(experiment.aggregate, reduction = "harmony", dims = 1:30) 
	experiment.aggregate <- FindClusters(object = experiment.aggregate, resolution = seq(0,2,0.1), verbose = FALSE)
	
	saveRDS(experiment.aggregate, file=paste0(result_loc,"rds/","HC_all_harmony_rm_BD.rds"))

	
### Part II All subjects on Local	
	library(clustree)
	library(Seurat)
	library(knitr)
	library(ggplot2)
	library(cowplot)
	library(dplyr)
	library(clusterExperiment)
	########
	setwd("C:/Users/zhangj/Documents/Zhang/Hippocampus/Catherine/Data/WGCNA/byStrain/code/scRNA/")
	setwd("C:/Users/zhangj/Documents/Zhang/Hippocampus/Catherine/Data/WGCNA/byStrain/code/scRNA/All_samples/results_HC/")

	comRS = readRDS("./rds/HC_all_harmony_rm_BD.rds")

	DimPlot(object = comRS, reduction = "tsne", group.by = "orig.ident", label = F,repel = TRUE)
	DimPlot(object = comRS, reduction = "umap", group.by = "orig.ident", label = F,repel = TRUE)

	DimPlot(object = comRS, reduction = "tsne", group.by = "age", label = F,repel = TRUE)
	DimPlot(object = comRS, reduction = "umap", group.by = "age", label = F,repel = TRUE)

	DimPlot(object = comRS, reduction = "tsne", group.by = "stain", label = F,repel = TRUE)
	DimPlot(object = comRS, reduction = "umap", group.by = "stain", label = F,repel = TRUE)

	DimPlot(object = comRS, reduction = "tsne", group.by = "geno", label = F,repel = TRUE)
	DimPlot(object = comRS, reduction = "umap", group.by = "geno", label = F,repel = TRUE)

	DimPlot(object = subset(comRS,R_S %in% c("Sus","Res")), reduction = "tsne", group.by = "R_S", label = F,repel = TRUE)
	DimPlot(object = subset(comRS,R_S %in% c("Sus","Res")), reduction = "umap", group.by = "R_S", label = F,repel = TRUE)

	DimPlot(object = comRS, reduction = "tsne", group.by = "R_S", label = F,repel = TRUE)
	DimPlot(object = comRS, reduction = "umap", group.by = "R_S", label = F,repel = TRUE)


	DimPlot(object = comRS, group.by="RNA_snn_res.0.5", pt.size=0.5, do.label = TRUE, reduction = "tsne", label = T)
	DimPlot(object = comRS, group.by="RNA_snn_res.0.5", pt.size=0.5, do.label = TRUE, reduction = "umap", label = T)


	Idents(comRS) <- "RNA_snn_res.0.5"
	write.csv(table(Idents(comRS),comRS$orig.ident),"Cluster_sample.csv")
	
###0.5	
	
	
	DefaultAssay(comRS) <- "RNA"
	genes = c('Rbfox3','Gad1','Gad2','Slc17a6','Slc17a7','C1qb','P2ry12','Inpp5d','Cyth4','Abca9','Csf1r','Aqp4','Gja1','S1pr1','Mbp','Plp1','Mog','Tnr','Cspg4','Pdgfra','Flt1','Dcn','Pecam1')
	DotPlot(comRS, features = genes, cols = c("white", "blue")) + RotatedAxis()
	
	my_level = c(1:2,5:7,9,11,13,16,18,21:23,10,14,15,17,19,20,3,24,4,31,0,12,29,26,8,25,27,28,30)
	levels(comRS)=my_level
	DotPlot(comRS, features = genes, cols = c("white", "blue")) + RotatedAxis()
	

