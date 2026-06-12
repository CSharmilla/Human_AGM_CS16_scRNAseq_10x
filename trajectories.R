library(here)
library(scater)
library(slingshot)


sce.clus <- readRDS(here("output", "sce.clus.RDS"))
rowdata <- data.frame(rowData(sce.clus))

k <- "cluster_k5" 
plotUMAP(sce.clus, colour_by = k, point_size=1, text_by = k, text_size =5)

# with slingshot

# fit a single principal curve to yields a pseudotime ordering of cells based on their relative positions when projected onto the curve
sce.sling <- slingshot(sce.clus, reducedDim='pca_denoise', cluster = k)
# creating matrix of pseudotimes
pseudo.paths <- slingPseudotime(sce.sling)
head(pseudo.paths)

# Taking the rowMeans just gives us a single pseudo-time for all cells. Cells
# in segments that are shared across paths have similar pseudo-time values in 
# all paths anyway, so taking the rowMeans is not particularly controversial.
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

# Need to loop over the paths and add each one separately.
gg <- plotUMAP(sce.sling, colour_by=I(shared.pseudo))
embedded <- embedCurves(sce.sling, "UMAP")
embedded <- slingCurves(embedded)
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg <- gg + geom_path(data=embedded, aes(x=UMAP1, y=UMAP2), size=1.2)
}

gg
