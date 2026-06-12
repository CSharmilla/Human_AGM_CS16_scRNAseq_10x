
library(here)
library(scran)

sce.filtered2 <- readRDS(file = here("output", "sce.filtered2.RDS"))

## Normalisation

### Library size normalisation

# library size is the total sum of counts across all genes for each cell
# library size factor is The “library size factor” for each cell is then directly proportional to its library size where the proportionality constant is defined such that the mean size factor across all cells is equal to 1
# from scran package 
# normalisation by deconvolution
set.seed(100)
clust.sce <- quickCluster(sce.filtered2)
sce.filtered2 <- computeSumFactors(sce.filtered2, cluster=clust.sce, min.mean=0.1)
sce.filtered2 <- logNormCounts(sce.filtered2)

summary(sizeFactors(sce.filtered2))

# plotting sf normalised data
hist(log10(sce.filtered2$sizeFactor), xlab="Log10[Size factor]", col='grey80')

saveRDS(sce.filtered2, file = here("output", "sce.filtered2.norm.RDS"))
