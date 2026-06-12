library(here)



sce.clus <- readRDS(here("output", "sce.clus.RDS"))


### selecting "bad" clusters

### by umi count per cluster
pct_cells<- 50
min_umi  <- lapply(det_sce_qc, function(x){summary(colData(x)$sum)}) #3000 # using a the 1st quartile for each dataset

### Select clusters with `r pct_cells` % cells having less than `r min_umi`umi counts
umi_df <- low_umi_clusters <- list()
for (n in names(all.sce.norm)){
  umi_df[[n]] <- as.data.frame(colData(all.sce.norm[[n]])[,c("sum", k)])
  colnames(umi_df[[n]]) <- c("umi", "cluster_names")
  
  # calculate the percentages of cells that are "low umi"
  umi_df[[n]] <-
    umi_df[[n]] %>%
    mutate(low_umi = as.numeric(umi) < min_umi[[n]][2]) %>%
    group_by(cluster_names) %>%
    summarise(n_low_umi = sum(low_umi), n_cells = n()) %>%
    mutate(pct_cell_low_umi = (n_low_umi/n_cells)*100) %>%
    arrange(desc(pct_cell_low_umi))
  
  # Filter the clusters that have a high proportion of "low umi"
  low_umi_clusters[[n]] <-
    umi_df[[n]] %>%
    dplyr::filter(pct_cell_low_umi > pct_cells) #%>% 
  #.$cluster_names
  print(low_umi_clusters[[n]])
  print(n)
  print(umi_df[[n]])
}


### by mito percent per cluster
pct_cells<- 50
pct_mt    <- 10 #lapply(det_sce_qc, function(x){summary(colData(x)$subsets_Mito_percent)})

# Select clusters with % of cells with higher mito % than threshold
mt_df <- high_mito_clusters <- list()
for (n in names(all.sce.norm)){
  mt_df[[n]] <- as.data.frame(colData(all.sce.norm[[n]])[,c("subsets_Mito_percent", k)])
  colnames(mt_df[[n]]) <- c("subsets_Mito_percent", "cluster_names")
  
  # calculate the percentages of cells that are "high mito"
  mt_df[[n]] <-
    mt_df[[n]] %>%
    mutate(high_pct = as.numeric(subsets_Mito_percent) > pct_mt) %>%
    group_by(cluster_names) %>%
    summarise(n_high_mt = sum(high_pct), n_cells = n()) %>%
    mutate(pct_cell_high_mt = (n_high_mt/n_cells)*100) %>% 
    arrange(desc(pct_cell_high_mt))
  
  # Filter the clusters that have a high proportion of "high mito"
  high_mito_clusters[[n]] <-
    mt_df[[n]] %>%
    dplyr::filter(pct_cell_high_mt > pct_cells) #%>% 
  #.$cluster_names
  print(n)
  print(mt_df[[n]])
}

### select clusters with high ribosomal genes
pct_cells <- 50
pct_ribo    <- 30 #lapply(det_sce_qc, function(x){summary(colData(x)$subsets_Ribo_percent)})
# Select clusters with pct_cells % of cells with higher that pct_ribo % of ribo genes
ribo_df <- high_ribo_clusters <- list()
for (n in names(all.sce.norm)){
  ribo_df[[n]] <- as.data.frame(colData(all.sce.norm[[n]])[,c("subsets_Ribo_percent", k)])
  colnames(ribo_df[[n]]) <- c("subsets_Ribo_percent", "cluster_names")
  
  # calculate the percentages of cells that are "high mito"
  ribo_df[[n]] <-
    ribo_df[[n]] %>%
    mutate(high_pct = as.numeric(subsets_Ribo_percent) > pct_ribo) %>%
    group_by(cluster_names) %>%
    summarise(n_high_ribo = sum(high_pct), n_cells = n()) %>%
    mutate(pct_cell_high_ribo = (n_high_ribo/n_cells)*100) %>% 
    arrange(desc(pct_cell_high_ribo))
  
  # Filter the clusters that have a high proportion of "high mito"
  high_ribo_clusters[[n]] <-
    ribo_df[[n]] %>%
    dplyr::filter(pct_cell_high_ribo > pct_cells) #%>% 
  #.$cluster_names
  print(n)
  print(high_ribo_clusters[[n]])
  print(ribo_df[[n]])
  
}

### CHECK: what happens when removing ribosomal genes
# removing ribosomal genes to see if it makes a difference to pca #####
ribo.test <- all.sce.norm
for (a in names(ribo.test)){
  rownames(ribo.test[[a]]) <- rowData(ribo.test[[a]])$Symbol
}
ribo_genes_list <- all.hvgs.symbols <- hvg.no.ribo <- plot.list <- list()
for (n in names(ribo.test)) {
  ribo_genes_list[[n]] <- grep("^RP[SL]", rownames(ribo.test[[n]]), value = TRUE)
  ribo.test[[n]] <- ribo.test[[n]][!(rownames(ribo.test[[n]]) %in% ribo_genes_list[[n]]),]
  
  all.hvgs.symbols[[n]] <- rowData(ribo.test[[n]])$Symbol[which(all.hvgs[[n]] %in% rowData(ribo.test[[n]])$ID)]
  hvg.no.ribo[[n]] <- all.hvgs.symbols[[n]][!(all.hvgs.symbols[[n]] %in% ribo_genes_list[[n]])]
  # PCA
  ribo.test[[n]] <- runPCA(ribo.test[[n]], subset_row=hvg.no.ribo[[n]])
  
  
  # plotting TSNE
  set.seed(00101001101)
  # When dimred is specified, no additional feature selection or standardization is performed. This means that any settings of ntop, subset_row and scale are ignored.
  ribo.test[[n]] <- runTSNE(ribo.test[[n]], dimred="PCA", perplexity=90)
  plot.list[[1]] <- plotReducedDim(ribo.test[[n]], dimred="TSNE", colour_by = "subsets_Ribo_percent") + ggtitle(paste0("TSNE - without ribo genes for ", n))
  plot.list[[2]] <- plotReducedDim(all.sce.norm[[n]], colour_by="subsets_Ribo_percent", dimred = "TSNE") + ggtitle(paste0("TSNE - with ribo genes for ", n))
  
  print(patchwork::wrap_plots(plot.list))
}
### not sure - will have to see properly later ###


# sce.doublet <- all.sce.norm
# ### doublets
# dbl.dens <- list()
# for (i in names(sce.doublet)){
#   dbl.dens[[i]] <- computeDoubletDensity(sce.doublet[[i]], subset.row=rowData(sce.doublet[[i]])$subset, d=ncol(reducedDim(sce.doublet[[i]], "pca_denoise")))
#   sce.doublet[[i]]$DoubletScore <- dbl.dens[[i]]
#   print(plotTSNE(sce.doublet[[i]], colour_by="DoubletScore") + ggtitle(paste0("doublet scores for ", i)))
# }
# 
# # how many cells are singlets and doublets?
# dbl.calls <- list()
# for (i in names(dbl.dens)){
#   dbl.calls[[i]] <- doubletThresholding(data.frame(score=dbl.dens[[i]]),
#                                    method="griffiths", returnType="call")
#   sce.doublet[[i]]$DoubletCalls <- dbl.calls[[i]]
#   print(summary(dbl.calls[[i]]))
#   print(plotColData(sce.doublet[[i]], x="cluster_k10", y="DoubletScore", colour_by=I(dbl.calls[[i]])) + ggtitle(paste0("doublet scores for ", i)))
#   
# }

### doublet finder

# All cells from a cluster with a large average doublet score should be considered suspect, and close 
# neighbors of problematic clusters should be treated with caution. A cluster containing only a small 
# proportion of high-scoring cells is safer, though this prognosis comes with the caveat that true doublets often lie 
# immediately adjacent to their source populations and end up being assigned to the same cluster. It is worth confirming 
# that any interesting results of downstream analyses are not being driven by those cells, e.g., by checking that DE in 
# an interesting gene is not driven solely by cells with high doublet scores. While clustering is still required for 
# interpretation, the simulation-based strategy is more robust than findDoubletClusters() to the quality of the clustering 
# as the scores are computed on a per-cell basis.

sce.doublet <- all.sce.norm
dbl.finder <- list()
for (i in names(sce.doublet)){
  dbl.finder[[i]] <- scDblFinder(sce.doublet[[i]], clusters=sce.doublet[[i]]$cluster_k10)
  print(plotTSNE(dbl.finder[[i]], colour_by="scDblFinder.score") + ggtitle(paste0("doublet finder scores for ", i)))
  print(table(dbl.finder[[i]]$scDblFinder.class))
}


### filtering out doublets
sce_fil <- lapply(dbl.finder, function(a){a[,!(a$scDblFinder.class == "doublet")]})