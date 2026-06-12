library(here)
library(scater)
library(edgeR)

sce.filtered <- readRDS(file = here("output", "sce.filtered.RDS"))
sce_qc <- readRDS(file = here("output", "sce_qc.RDS"))

### per feature quality control

sce.filtered <- addPerFeatureQC(sce.filtered) # detected - the percentage of cells in which the gene was detected; 
# and mean - the mean UMI count for the gene across all cells
rowData(sce.filtered)

# cell sparsity: for each cell, what is the proportion of genes undetected?
#colData(sce.filtered)$cell_sparsity <- 1 - (colData(sce.filtered)$detected / nrow(sce.filtered)) # 1 - the number of genes detected per cell by the total number of genes

# gene sparsity: for each gene, what is the proportion of cells in which it is not detected
rowData(sce.filtered)$gene_sparsity <- (100 - rowData(sce.filtered)$detected) / 100

#hist(colData(sce.filtered)$cell_sparsity, breaks=50, col="grey80", xlab="Cell sparsity", main="")
hist(rowData(sce.filtered)$gene_sparsity, breaks=50, col="grey80", xlab="Gene sparsity", main="")

# removing genes present in less than x cells
min.cells <- 1 - (5 / ncol(sce.filtered))
sparse.genes <- rowData(sce.filtered)$gene_sparsity > min.cells
table(sparse.genes)

rowData(sce.filtered)$sparse.genes <- sparse.genes

# filter out sparse genes - not filtering these out since sometimes they seem to filter out important genes to certain cell populations
#sce.filtered.sparse <- sce.filtered[!rowData(sce.filtered)$sparse.genes, ]
#sce.filtered.sparse
# recalculate QC metrics
#colData(sce.filtered.sparse) <- colData(sce.filtered.sparse)[,1:2]
#sce.filtered.sparse <- addPerCellQC(sce.filtered.sparse)
#colData(sce.filtered.sparse)

a <- data.frame(colData(sce.filtered))

### diagnosing cell type loss

discard <- sce_qc$discard
lost <- calculateAverage(counts(sce_qc)[,discard])
kept <- calculateAverage(counts(sce_qc)[,!discard])

logged <- cpm(cbind(lost, kept), log= TRUE, prior.count=2)
logFC <- logged[,1] - logged[,2]
abundance <- rowMeans(logged)

## this plot shows Log-fold change in expression in the discarded cells compared to the retained cells
plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16)
is.mito <- grep("MT-", rowData(sce_det)$Symbol)
points(abundance[is.mito], logFC[is.mito], col="dodgerblue", pch=16)


saveRDS(sce.filtered, file = here("output", "sce.filtered2.RDS"))
