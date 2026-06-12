
library(here)
library(scran)


sce.filtered2.norm <- readRDS(file = here("output", "sce.filtered2.norm.RDS"))

## Feature Selection

### per-gene variation
# simplest way is to compute the variance of the log-normalized expression values (i.e., “log-counts” ) for each gene across all cells

# modelGeneVar() fits a trend to the variance with respect to abundance across all genes
genevar_sce <- modelGeneVar(sce.filtered2.norm)
# now to visualise the fit
# this assumes that for any abundance, variation in expression for most genes in driven by uninteresting processes like sampling
# this is not always right since sometimes many genes at a particular abundance
# are affected by a biological process
# this would inflate the fitted trend in that abundance interval and compromise the detection of relevent genes
fit <- metadata(genevar_sce)
plot(fit$mean, fit$var, xlab="mean of log-expression", ylab="variance of log-expression")
curve(fit$trend(x), col="blue", add = TRUE, lwd = 3)


# quantifying technical noise - in the absence of spike in data, you can create a trend by making distributional assumptions about the noise. assuming a poisson distribution for technical noise. 
# leads to increased residual (higher error) compared to the previous plot. this can be interpreted as there being more biological variance. 
poisson_genevar_sce <- modelGeneVarByPoisson(sce.filtered2.norm)
poisson_genevar_sce_order <- poisson_genevar_sce[order(poisson_genevar_sce$bio, decreasing=TRUE),]
poisson_fit <- metadata(poisson_genevar_sce_order)
plot(poisson_genevar_sce_order$mean, poisson_genevar_sce_order$total, xlab="mean of log-expression", ylab="variance of log-expression")
curve(poisson_fit$trend(x), col="blue", add = TRUE, lwd = 3)

# selecting the most highly variable genes
hvg <- getTopHVGs(poisson_genevar_sce, prop = 0.3)

sce_poisson_hvg <- sce.filtered2.norm[hvg,] 
sce_poisson_hvg
sce.filtered2.norm
#poisson_hvg_df <- data.frame(rowData(sce_poisson_hvg))

## adding hvgs sce and the original sce to the same structure
#altExp(sce.filtered2.norm, "poisson_hvg") <- sce_poisson_hvg
#altExpNames(sce_poisson_hvg)

# adding row subset
rowSubset(sce.filtered2.norm, "hvg") <- hvg

saveRDS(sce.filtered2.norm, file = here("output", "sce.hvg.RDS"))
saveRDS(poisson_genevar_sce, file = here("output", "poisson_genevar_sce.RDS"))


