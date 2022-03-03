library(ggplot2)
library(Seurat)

setwd("~/Desktop/Projects/scRNAseq/aov/")

'%!in%' <- function(x,y)!('%in%'(x,y))
## -----------------------------------------------------------------------------------------
hc <- readRDS("HC_harmony_all_Run_doublets_5%_BXD161.rds")

# rerun doublet finder
library(DoubletFinder)
ids <- as.vector(unique(hc$orig.ident))
Drate = 0.05
Idents(hc) <- 'RNA_snn_res.2'
tm0 = NULL
for (i in ids){
  print(i)
  seu_kidney = subset(hc, orig.ident== i)
  seu_kidney = doublet_D(seu_kidney, Drate)
  colnames(seu_kidney)=c(1:ncol(seu_kidney))
  tm0 = rbind(tm0, seu_kidney)
}

tm1=tm0[rownames(hc@meta.data),]
hc$Doublets = tm1[,dim(tm1)[2]]

DimPlot(object = hc, reduction = "umap" ,group.by = 'Doublets')
DimPlot(object = hc, reduction = "tsne" ,group.by = 'Doublets')

saveRDS(hc, "HC_harmony_all_Run_doublets_5%_BXD161_doublt2.rds")

res2.0 <- hc@meta.data$Doublets
res0.5 <- hc@meta.data$Doublet

hc <- readRDS("HC_harmony_all_Run_doublets_5%_BXD161_doublt2.rds")
df <- hc@meta.data
doubltes <- c()
clusters <- as.character(unique(hc$RNA_snn_res.2))
for (cluster in clusters){
  dat <- subset(df, RNA_snn_res.2 == cluster)
  d <- sum(dat$Doublets == "Doublet")
  s <- sum(dat$Doublets == 'Singlet')
  t <- nrow(dat)
  doublets.per.cluster <- data.frame(Cluster = cluster, n.Singlet = s, n.Doublets = d, n.Cells = t)
  doubltes <- rbind(doubltes, doublets.per.cluster)
} 
doubltes$pct.doublets <-(doubltes$n.Doublets/doubltes$n.Cells) * 100
DimPlot(object = hc, reduction = "umap" ,group.by = 'RNA_snn_res.2', , label = T,repel = TRUE, pt.size=0.25)

Idents(hc) <- 'RNA_snn_res.0.5'
DefaultAssay(hc) <- "RNA"
new.cluster.ids <- c('Oligodendrocytes' ,'Glutamatergic Neurons', 'Glutamatergic Neurons' , 'Glutamatergic Neurons', 'Astrocytes',
                     'Microglia', 'Glutamatergic Neurons', 'Glutamatergic Neurons', 'Glutamatergic Neurons', 'GABAergic Neurons',
                     'Glutamatergic Neurons', 'Glutamatergic Neurons', 'OPCs', 'GABAergic Neurons', 'Glutamatergic Neurons',
                     'GABAergic Neurons', 'Glutamatergic Neurons', 'GABAergic Neurons',
                     'Glutamatergic Neurons' ,'Glutamatergic Neurons', 'Differentiating Cells', 'Glutamatergic Neurons',
                     'Glutamatergic Neurons', 'Doublet', 'Endothelial Cells', 'Cajal-Retzius cells', 'Glutamatergic Neurons',
                     'Doublet', 'Doublet', 'OPCs', 'Endothelial Cells', 'Glutamatergic Neurons', 'Glutamatergic Neurons',
                     'Doublet', 'Doublet', 'Glutamatergic Neurons', 'Doublet')
names(new.cluster.ids) <- c('0' ,'1', '2' , '3', '4' , '5',
                            '6', '7', '8', '9', '10', '11', 
                            '12', '13', '14', '15', '16', '17',
                            '18' ,'19', '20', '21', '22', '23', 
                            '24', '25', '26', '27', '28', '29',
                            '30', '31', '32', '33', '34', '35', '36')
hc <- RenameIdents(hc, new.cluster.ids)
hc$Cell.Type <- Idents(hc)

doublet.clusters <- subset(doubltes, pct.doublets >= 30)$Cluster %>% as.character()
for (cluster in doublet.clusters) {
  print(cluster)
  hc <- subset(hc, RNA_snn_res.2 != cluster)
  gc()
}
DimPlot(object = hc, reduction = "umap" ,group.by = 'RNA_snn_res.0.5', label = T,repel = TRUE, pt.size=0.25)
DimPlot(object = hc, reduction = "umap" ,group.by = 'Cell.Type', label = T,repel = TRUE, pt.size=0.25)

