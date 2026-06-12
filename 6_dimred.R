
library(here)
library(scran)
library(PCAtools)


sce.hvg <- readRDS(file = here("output", "sce.hvg.RDS"))
poisson_genevar_sce <- readRDS(file = here("output", "poisson_genevar_sce.RDS"))

if (!file.exists(here("output", "sce.hvg.dimred.RDS"))){
  # PCA using the same HVG data from above
  # sce.hvg <- fixedPCA(sce.hvg)
  # reducedDimNames(sce.hvg)
  # dim(reducedDim(sce.hvg, "PCA"))
  # 
  # percent.var <- attr(reducedDim(sce.hvg), "percentVar")
  # 
  # # find elbow
  # chosen.elbow <- findElbowPoint(percent.var)
  # chosen.elbow
  # plot(percent.var, xlab="PC", ylab="Variance explained (%)")
  # abline(v=30, col="red") # picked manually because findelbowpoint function doesnt work
  # 
  # plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")
  # plotReducedDim(sce_poisson_hvg, dimred="PCA")
  
  # using technical noise to identify threshold for informative PCs instead of above - this method performs best when men-variance trend reflects the actual technical noise therefore using poisson
  set.seed(111001001)
  sce.hvg.denoise.pca <- denoisePCA(sce.hvg, technical=poisson_genevar_sce, subset.row = hvg, name = "pca_denoise") # make sure to use unordered original sce object otherwise it throws an error.
  ncol(reducedDim(sce.hvg.denoise.pca))
  
  # creating a new entry with only the interesting pcs
  #reducedDim(sce.hvg, "PCA.denoised") <- reducedDim(sce.hvg)[,1:ncol(reducedDim(denoised))]
  #reducedDimNames(sce_poisson_hvg)
  
  # tsne
  set.seed(00101001101)
  # tSNE is very computationally intensive. To mitigate this, using dimred="PCA" will instruct the function to use the top PCs in the sce to do tSNE calculations.
  sce.hvg.denoise.pca <- runTSNE(sce.hvg.denoise.pca, dimred="pca_denoise", perplexity=100)
  plotReducedDim(sce.hvg.denoise.pca, dimred="TSNE")
  
  # UMAP
  set.seed(00101101)
  sce.hvg.denoise.pca <- runUMAP(sce.hvg.denoise.pca, dimred="pca_denoise", min_dist = 0.3)#, n_neighbors = 15)
  plotReducedDim(sce.hvg.denoise.pca, dimred="UMAP")
  
  #sce.filtered.sparse <- runTSNE(sce.filtered.sparse, dimred="PCA", perplexity=80)
  #plotReducedDim(sce.filtered.sparse, dimred="TSNE")
  
  saveRDS(sce.hvg.denoise.pca, file = here("output", "sce.hvg.dimred.RDS"))
} else{
  sce.hvg.dimred <- readRDS(file = here("output", "sce.hvg.dimred.RDS"))
  
}
