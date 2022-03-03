library(dplyr)
library(Seurat)
"%!in%" <- function(x,y)!("%in%"(x,y))

nCores <- 67 # set number of cores for de
resolution <- "RNA_snn_res.0.3" # desired cluster resolution
path <- getwd() # path to seurat object directory
sfname <- "New_Microglia_10PCs_1000Features_Iter1_2021-07-14.rds" # name of seurat object 

# read seurat object
hc <- readRDS(file.path(path, sfname))
Idents(hc) <- resolution
DefaultAssay(hc) <- "RNA"

clusters <- hc[[resolution]] %>% 
  unlist %>% 
  unique

# prepare count data for each cluster

# define groups
variable1 <- c("0", "5")
variable2 <- c("6", "14") # if not second variable set to 1; if 3 or more variable needed adjust code

threshold <- 20 # define threshold 
i <- 0

count_data <- list()
for (cluster in clusters){ 
  i <- i + 1
  
  df <- hc[,hc[[resolution]] == cluster]
  df.dat <- GetAssayData(df, slot = "data") %>% as.data.frame()
  df.meta <- df@meta.data
  
  genes.pct <- matrix(nrow = nrow(df.dat))
  
  # use only genes expressed in at least 20 percent of the cells in at least one experimental group
  for (f1 in variable1){ ## change to appropriate variable ##
    for (f2 in variable2){ ## change to appropriate variable ##
      
      df.subset <- df.dat[,rownames(subset(df.meta,  geno == f1 & age == f2))]
      df.subset.pct <- (1-(rowSums(df.subset == 0)/ncol(df.subset))) * 100
      res <- data.frame(df.subset.pct)
      colnames(res) <- paste(f1, f2, "pct", sep = "_")
      print(colnames(res) <- paste(f1, f2, "pct", sep = "_"))
      genes.pct <- cbind(genes.pct, res)
    }
  }
  
  pct.threshold <- 20 # define threshold 
  genes.pct <- genes.pct[,-1]
  genes.filtered <- apply(genes.pct, 1, function(x) x > threshold ) %>% colMeans() 
  genes.filtered <- genes.filtered[genes.filtered != 0] %>% names
  
  # save normalized data
  normdat <- GetAssayData(df, slot = "data") %>% as.data.frame()
  normdat <- normdat[,rownames(df.meta)]
  normdat <- normdat[genes.filtered,]
  # save count data
  countdat <- GetAssayData(df, slot = "counts") %>% as.data.frame()
  countdat <- countdat[,rownames(df.meta)]
  countdat <- countdat[genes.filtered,]
  
  combined <- list(Counts = countdat, NormData = normdat, meta.data = df.meta)
  count_data[[i]] <- combined
  names(count_data)[i] <- paste0('Cluster', cluster)
}
saveRDS(count_data, "ClusterCountDat.RDS")


###########################
# differential expression #
###########################
count_data <- readRDS("ClusterCountDat.RDS")
tryCatch({
  results <- lapply(count_data, function(dat){
    # read count matrix
    #dat <- x
    df.meta <- dat$meta.data
    df.meta <- within(df.meta, {
      age <- factor(age, levels=c("6", "14"))
      geno <- factor(geno, levels=c("0", "5"))
      orig.ident <- factor(orig.ident)
      Run <- factor(Run)
    })
    genes <- rownames(dat[["Counts"]])[1:20]
    
    library(doParallel)
    registerDoParallel(cores=nCores)
    
    x <- foreach (i=1:length(genes), .combine = rbind) %dopar% {
      tryCatch({
        
        require(glmmTMB); require("bbmle"); require("dplyr")
        g <- genes[i]
        lmm.res <- c()
        df.meta$expression <- dat[["Counts"]][g,] %>% unname() %>% unlist()
        df.meta$normexpression <- dat[["NormData"]][g,] %>% unname() %>% unlist()
        
        # differential expression using negative binomial mixed models
        # original identity is used as random variable
        # Run was added to the model as batch factor
        mod.tbm <- glmmTMB::glmmTMB(expression ~  age + geno + Run + (1|orig.ident), data = df.meta, family=glmmTMB::nbinom2(link = "log"), ziformula= ~0) ## adjut model as needed
        
        mod.lmm <- car::Anova(mod.tbm) %>% as.data.frame() # get main effect
        mod.lmm$gene <- g
        mod.lmm$contrast <- rownames(car::Anova(mod.tbm) %>% as.data.frame())
        lmm.res <- rbind(lmm.res, mod.lmm)
        
        remove(mod.tbm, mod.lmm)
        return(lmm.res)
      }, error = function(e){c(print(paste0("ERROR: ", g, " did not converge")))}) # model convergence problem are frequent for clusters with low count data likely due to poor fit
    }
    stopImplicitCluster()
    log <- warnings()# log to store all warnings
    
    return(as.data.frame(x))
    
    #saveRDS(log, paste0(ct, '_Cluster',  cluster, '_WarningLog_ResSus.rds'))
  })
  }, error = function(e){})

saveRDS(results, 'diffexp.RDS')


# remove genes that did not converge
# add false discovery
diff_exp <- readRDS("H:/diffexp.RDS")
corrected_de <- 
  lapply(diff_exp, function(i){
  
  dt <- i[which(!is.na(i$Chisq)),] %>% 
    filter(contrast != 'Run') %>% 
    rename(pvalue = `Pr(>Chisq)` )
  
  dt <- split(dt, dt$contrast)
  dt <- lapply(dt, function(j){
    j$fdr <-  with(j, p.adjust(pvalue, 'fdr'))
    return(j)
  })
  dt <- do.call(rbind, dt)
  return(dt)
})


# bind all clusters
merged_de <- c()
for(i in seq_along(corrected_de)){
  n <- names(corrected_de)[i]
  tmp <- corrected_de[[i]]
  tmp$Cluster <- n %>% gsub('Cluster', '', .)
  merged_de <- rbind(merged_de, tmp) 
}

rownames(merged_de) <- NULL
