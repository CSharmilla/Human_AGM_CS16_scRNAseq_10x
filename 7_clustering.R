library(here)
library(bluster)
library(clustree)
library(scater)
library(scuttle)


sce.hvg.dimred <- readRDS(file = here("output", "sce.hvg.dimred.RDS"))

# clustering
if (!file.exists(here("output", "clusters_list.RDS"))){
  
  ks <- c(5, 10, 20, 30, 40)
  names <- paste0("cluster_k", as.character(ks)) 
  clusters_list <- list()
  
  
  clusters_list <- mapply(function(cluster_name, k) {clusterRows(reducedDim(sce.hvg.dimred, "pca_denoise"), NNGraphParam(k = k, cluster.fun="louvain"))},
                          names,
                          ks,
                          USE.NAMES = TRUE, # use names of first argument
                          SIMPLIFY = FALSE)
  saveRDS(clusters_list, here("output", "clusters_list.RDS"))
} else{
  clusters_list <- readRDS(here("output", "clusters_list.RDS"))
}


# plotting and saving clusters in sce object
if (!file.exists(here("output", "sce.clus.RDS"))){
  # creating new to keep this neat
 sce.clus <- sce.hvg.dimred
  ### A higher k generally corresponds to a lower resolution clustering
  
  idx <- c(1:length(names))
  new_names <- paste0("cluster_name_", rev(idx))
  
  for (i in idx){
    k <- ks[i]
    cluster_name <- names[i]
    colData(sce.clus)[names[i]] <- clusters_list[cluster_name]
    colData(sce.clus)[new_names[i]] <- colData(sce.clus)[names[i]]
    # and plot the result
    # plotDim <-
    #   plotReducedDim(
    #     sce.clus,
    #     "TSNE",
    #     colour_by = cluster_name,
    #     point_size = 0.6,
    #     point_alpha = 0.3,
    #     text_by = cluster_name,
    #     text_size = 3
    #   )  + scale_color_hue() + ggtitle(paste("Clustering with k =", k))
    # print(plotDim)
  }
  
  # changing rownames from ensembl to gene names
  # sanity check
  all(rownames(sce.clus) == rowData(sce.clus)[, "ID"])
  new_names <- uniquifyFeatureNames(rownames(sce.clus), rowData(sce.clus)[, "Symbol"])
  rownames(sce.clus) <- new_names
  rowData(sce.clus)$new_names <- new_names
  
  # saving
  saveRDS(sce.clus, here("output", "sce.clus.RDS"))
  
}else{
  sce.clus <- readRDS(here("output", "sce.clus.RDS"))
}

# cluster QC - k = 10 looks good

plot_list_clustree <- clustree(sce.clus, prefix = "cluster_name_", edge_arrow = FALSE) +
  scale_color_discrete(name="cluster_k", labels=c("40","30", "20", "10", "5"))  #+ ggtitle(n)#labels=c("8", "7", "6", "4","3", "2", "1"))
print(plot_list_clustree)

k <- "cluster_k10" 

plotUMAP(sce.clus, colour_by = k, point_size=1, text_by = k, text_size =5)
plotUMAP(sce.clus, colour_by = "subsets_Mito_percent", point_size=1, text_by = k, text_size =5)


