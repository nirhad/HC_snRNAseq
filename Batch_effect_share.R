############ integrate all HC/PFC samples

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
	  experiment.aggregate$strain <- experiment.aggregate$orig.ident
	  experiment.aggregate$age <- experiment.aggregate$orig.ident
	  experiment.aggregate$tissue <- experiment.aggregate$orig.ident
	  experiment.aggregate$mouse <- experiment.aggregate$orig.ident
	  experiment.aggregate$Run <- experiment.aggregate$orig.ident
	  experiment.aggregate$havest <- experiment.aggregate$orig.ident
	  experiment.aggregate$XRun <- experiment.aggregate$orig.ident
	  experiment.aggregate$i7Index <- experiment.aggregate$orig.ident
	  
	  for (i in 1:length(ids)){print(i)	  
		experiment.aggregate$R_S = mapvalues(experiment.aggregate$R_S, c(as.vector(Sample_ID$SCBID[i])), c(as.vector(Sample_ID$R_S[i])))
		experiment.aggregate$geno = mapvalues(experiment.aggregate$geno,  c(as.vector(Sample_ID$SCBID[i])), c(as.vector(Sample_ID$Genotype[i])))
		experiment.aggregate$strain = mapvalues(experiment.aggregate$strain,  c(as.vector(Sample_ID$SCBID[i])), c(as.vector(Sample_ID$Line_Short[i])))
		experiment.aggregate$age = mapvalues(experiment.aggregate$age,  c(as.vector(Sample_ID$SCBID[i])), c(as.vector(Sample_ID$Age[i])))
		experiment.aggregate$tissue = mapvalues(experiment.aggregate$tissue,  c(as.vector(Sample_ID$SCBID[i])), c(as.vector(Sample_ID$Sample[i])))
		experiment.aggregate$mouse = mapvalues(experiment.aggregate$mouse,  c(as.vector(Sample_ID$SCBID[i])), c(as.vector(Sample_ID$MouseID[i])))
		experiment.aggregate$Run = mapvalues(experiment.aggregate$Run,  c(as.vector(Sample_ID$SCBID[i])), c(as.vector(Sample_ID$Run[i])))
		experiment.aggregate$havest = mapvalues(experiment.aggregate$havest,  c(as.vector(Sample_ID$SCBID[i])), c(as.vector(Sample_ID$Harvest[i])))
		experiment.aggregate$XRun = mapvalues(experiment.aggregate$XRun,  c(as.vector(Sample_ID$SCBID[i])), c(as.vector(Sample_ID$XRun[i])))
		experiment.aggregate$i7Index = mapvalues(experiment.aggregate$i7Index,  c(as.vector(Sample_ID$SCBID[i])), c(as.vector(Sample_ID$i7Index[i])))
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
								
		
### 1.4 Cell filtering  (cell level filter)
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
	
#### 2. “pre-processing”	
	

	experiment.aggregate <- NormalizeData(experiment.aggregate) %>% FindVariableFeatures() %>% ScaleData(vars.to.regress =c("nCount_RNA", "percent.rpL","percent.rpS","percent.ps")) %>% RunPCA(verbose = FALSE)
	saveRDS(experiment.aggregate, file=paste0(result_loc,"rds/","HC_harmony_all_forBatch.rds"))
	
	### no batch correction
	pcs.use <- 30
	experiment.aggregate <- RunTSNE(object = experiment.aggregate, dims = 1:pcs.use)
	experiment.aggregate <- RunUMAP(object = experiment.aggregate, dims = 1:pcs.use)
	saveRDS(experiment.aggregate, file=paste0(result_loc,"rds/","HC_harmony_all_forBatch_1.rds"))
	
	### Batch correction using Harmony
	comRS = readRDS(paste0(result_loc,"rds/","HC_harmony_all_forBatch.rds"))
	batchid <- c("Run")
	comRS <- RunHarmony(comRS, group.by.vars = batchid,plot_convergence = TRUE)
	comRS <- RunUMAP(comRS, reduction = "harmony", dims = 1:30)
	comRS <- RunTSNE(comRS, reduction = "harmony", dims = 1:30)
	saveRDS(comRS, file=paste0(result_loc,"rds/","HC_harmony_all_forBatch_",batchid,".rds"))
	
	##### Clustering
	comRS = readRDS(paste0(result_loc,"rds/","HC_harmony_all_forBatch_Run.rds"))
	comRS <- FindNeighbors(comRS, reduction = "harmony", dims = 1:30) 
	comRS <- FindClusters(object = comRS, resolution = seq(0,2,0.1), verbose = FALSE)
	saveRDS(comRS, file=paste0(result_loc,"rds/","HC_harmony_all_forBatch_Run_clustering.rds"))
	
		
	
#### localcheck the plots for batch effect
	
	library(clustree)
	library(Seurat)
	library(knitr)
	library(ggplot2)
	library(cowplot)
	library(dplyr)

	########
	setwd("C:/Users/zhangj/Documents/Zhang/Hippocampus/Catherine/Data/WGCNA/byStrain/code/scRNA/All_samples/results_HC/")

	##### Step 1 check the batch befor correction
	comRS = readRDS("./rds/HC_harmony_all_forBatch_1.rds")
	#ElbowPlot(comRS)
	batchid <- "Run"   ### "havest" "mouse" "tissue" "Run" "geno" "age" "strain" "R_S" "i7Index" "10XRun"
	
	DimPlot(object = comRS, reduction = "pca" ,group.by = batchid)
	DimPlot(object = comRS, reduction = "tsne" ,group.by = batchid)
	DimPlot(object = comRS, reduction = "umap" ,group.by = batchid)
	
	##### Step 2 batch correction
	comRS = readRDS("./rds/HC_harmony_all_forBatch_Run.rds")  #### revising
	batchid <- c("Run")
	
	DimPlot(object = comRS, reduction = "pca" ,group.by = batchid)
	DimPlot(object = comRS, reduction = "tsne" ,group.by = batchid)
	DimPlot(object = comRS, reduction = "umap" ,group.by = batchid)

	
########### BATCH EFFECT	
### Check correlation varibles (batch)
library(vcd)
setwd("C:/Users/zhangj/Documents/Zhang/Hippocampus/Catherine/Data/WGCNA/byStrain/code/scRNA/All_samples/results_HC/")
options(stringsAsFactors = FALSE)	
a=read.csv("sample_metrics_final.csv")
a=subset(a,Sample=="HC")
xx=xtabs(~Harvest+XRun, data = a)
summary(assocstats(xx))


#################################################
#### impact of batch effect factors ANOVA Results
#################################################
	
	library(clustree)
	library(Seurat)
	library(knitr)
	library(ggplot2)
	library(cowplot)
	library(dplyr)

	########
	setwd("C:/Users/zhangj/Documents/Zhang/Hippocampus/Catherine/Data/WGCNA/byStrain/code/scRNA/All_samples/results_HC/")

	##### Step 1 check the batch befor correction
	comRS = readRDS("./rds/HC_harmony_all_forBatch_1.rds")

	b1=comRS@meta.data[,c("mouse", "Run","i7Index", "XRun","havest","age","geno","strain")]		### "havest" "mouse" "tissue" "Run" "geno" "age" "strain" "R_S" "i7Index" "10XRun"
	b2=comRS@reductions$pca@cell.embeddings			### comRS@reductions$harmony@cell.embeddings
	
#########  plot	

	library(ggplot2)
	library(reshape2)
	  
	  ####  pc_var_1 is varicae explaied; pc_var_2 is significance
		pc_var_1 <- matrix(NA,ncol=dim(b2)[2],nrow=dim(b1)[2])
		pc_var_2 <- pc_var_1
		
		for(i in 1:dim(b2)[2])
		{
		  print(i)
		  for(j in 1:dim(b1)[2])
		  {
			m1 = aov(b2[,i]~b1[,j])
			pc_var_1[j,i] <- round(prop.table(anova(m1)[,"Sum Sq"])*100,digits = 2)[1]
			pc_var_2[j,i] <- summary(m1)[[1]]$`Pr(>F)`[1]
		  }
		}
		
		pc_var_2[pc_var_2>0.05] = NA;pc_var_2[!is.na(pc_var_2)] = 1;
		pc_var_1 = pc_var_1*pc_var_2
		colnames(pc_var_1)=colnames(b2);rownames(pc_var_1)=colnames(b1)
	  
	  pc_var_1 = pc_var_1[,1:30]
	  #### option 1 for plot
	  melted_cordata <- melt(pc_var_1)
	  colnames(melted_cordata)= c("Var1","Var2","value")
	  
	  plot <- ggplot(data = melted_cordata, aes(x=Var1, y=Var2, fill=value, label= value))
	  plot_tile <- plot + geom_tile()
	  plot_fill_color <- plot_tile + scale_fill_gradient2(limits=c(0,10))
	  plot_label <- plot_fill_color + geom_text()
	  
	  ##### option 2 for plot
	  library("gplots")
	  my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 30)
	  heatmap.2(data.matrix(pc_var_1),cellnote=round(pc_var_1,3),notecol="black",density.info="none",trace="none",col = my_palette, margins=c(8,10),dendrogram="row",Colv="NA",scale="row",main="PC snRNAseq",cexRow=1,cexCol=1)
	  mtext("Association between PCs and Phenotypes")

	


			### Harmony correction (Run)
			#################################################
			#### impact of batch effect factors ANOVA Results
			#################################################

			library(clustree)
			library(Seurat)
			library(knitr)
			library(ggplot2)
			library(cowplot)
			library(dplyr)

			########
			setwd("C:/Users/zhangj/Documents/Zhang/Hippocampus/Catherine/Data/WGCNA/byStrain/code/scRNA/All_samples/results_HC/")

			##### Step 1 check the batch befor correction
			comRS = readRDS("./rds/HC_harmony_all_forBatch_Run_clustering.rds")
			b1=comRS@meta.data[,c("mouse", "Run","i7Index", "XRun","havest","age","geno","strain","RNA_snn_res.0.5")]		### "havest" "mouse" "tissue" "Run" "geno" "age" "strain" "R_S" "i7Index" "10XRun"
			b2=comRS@reductions$harmony@cell.embeddings

			#########  plot	

			library(ggplot2)
			library(reshape2)

			####  pc_var_1 is varicae explaied; pc_var_2 is significance
			pc_var_1 <- matrix(NA,ncol=dim(b2)[2],nrow=dim(b1)[2])
			pc_var_2 <- pc_var_1

			for(i in 1:dim(b2)[2]){
			  print(i)
			  for(j in 1:dim(b1)[2]){
				m1 = aov(b2[,i]~b1[,j])
				pc_var_1[j,i] <- round(prop.table(anova(m1)[,"Sum Sq"])*100,digits = 2)[1]
				pc_var_2[j,i] <- summary(m1)[[1]]$`Pr(>F)`[1]
			  }
			}

			pc_var_2[pc_var_2>0.05] = NA;pc_var_2[!is.na(pc_var_2)] = 1;
			pc_var_1 = pc_var_1*pc_var_2
			colnames(pc_var_1)=colnames(b2);rownames(pc_var_1)=colnames(b1)
			v=pc_var_1
			#### option 1 for plot
			pc_var_1 = pc_var_1[,1:30]
			melted_cordata <- melt(pc_var_1)
			colnames(melted_cordata)= c("Var1","Var2","value")

			plot <- ggplot(data = melted_cordata, aes(x=Var1, y=Var2, fill=value, label= value))
			plot_tile <- plot + geom_tile()
			plot_fill_color <- plot_tile + scale_fill_gradient2(limits=c(0,10))
			plot_label <- plot_fill_color + geom_text()

			##### option 2 for plot
			library("gplots")
			my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 100)
			heatmap.2(data.matrix(pc_var_1),cellnote=round(pc_var_1,3),notecol="black",density.info="none",trace="none",col = my_palette, margins=c(8,10),dendrogram="row",Colv="NA",scale="row",main="Harmony snRNAseq",cexRow=1,cexCol=1)
			mtext("Association between PCs and Phenotypes")

