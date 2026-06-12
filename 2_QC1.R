library(here)
library(scater)




sce <- readRDS(here("output", "sce"))



# Find undetected genes
detected_genes <- rowSums(counts(sce)) > 0
table(detected_genes)
sce_det <- sce[detected_genes,]

rowData(sce_det)
colData(sce_det)


# use "subsets" to get metrics for the subset of mitochondrial genes
sce_qc <- addPerCellQC(sce_det, subsets=list(Mito=grep("MT-", rowData(sce_det)$Symbol), Ribo=grep("^RP[SL]", rowData(sce_det)$Symbol)))      #grep("MT-", rowData(sce)$Symbol)))
colData(sce_qc)

summary(colData(sce_qc)$sum)  # library size (total sum of UMI counts across all genes)
summary(colData(sce_qc)$detected)  # number of genes detected 
summary(colData(sce_qc)$subsets_Mito_percent) # mitochondrial proportions

# create plots to check data
scater::plotColData(sce_qc, x="Sample", y="sum") + scale_y_log10() + ggtitle("Total count")
plotColData(sce_qc, x="Sample", y="detected") + scale_y_log10() + ggtitle("Detected features")
plotColData(sce_qc, x="Sample", y="subsets_Mito_percent") + scale_y_log10() + ggtitle("Mito percent")


# adaptive thresholds
# library size
low_lib_size <- isOutlier(sce_qc$sum, log=TRUE, type="lower") # default nmads = 3
table(low_lib_size)
attr(low_lib_size, "thresholds")
colData(sce_qc)$low_lib_size <- low_lib_size
plotColData(sce_qc, 
            x="Sample", 
            y="sum",
            colour_by = "low_lib_size") + 
  scale_y_log10() + 
  labs(y = "Total count", title = "Total count") +
  guides(colour=guide_legend(title="Discarded"))

#reasons <- perCellQCFilters(sce_qc, 
#    sub.fields=c("subsets_Mito_percent"))
#reasons_df <- data.frame(reasons)
#colSums(as.matrix(reasons))

# number of genes
low_n_features <- isOutlier(sce_qc$detected, log=TRUE, type="lower")
table(low_n_features)
attr(low_n_features, "thresholds")
colData(sce_qc)$low_n_features <- low_n_features
plotColData(sce_qc, 
            x="Sample", 
            y="detected",
            colour_by = "low_n_features") + 
  scale_y_log10() + 
  labs(y = "Genes detected", title = "Genes detected") +
  guides(colour=guide_legend(title="Discarded"))

# mito genes
high_Mito_percent  <- isOutlier(sce_qc$subsets_Mito_percent, type="higher") #For the mitochondrial content the exclusion zone is in the 
# higher part of the distribution. For this reason we do not need to worry about log transforming 
# the data as want to remove the long right hand tail anyway
table(high_Mito_percent)
attr(high_Mito_percent, "thresholds")
colData(sce_qc)$high_Mito_percent <- high_Mito_percent
plotColData(sce_qc, 
            x="Sample", 
            y="subsets_Mito_percent",
            colour_by = "high_Mito_percent") + 
  labs(y = "Percentage mitochondrial UMIs", title = "Mitochondrial UMIs") +
  guides(colour=guide_legend(title="Discarded"))

#df <- data.frame(colData(sce_qc))
# summary of discarded cells
data.frame(`Library Size` = sum(low_lib_size),
           `Genes detected` = sum(low_n_features),
           `Mitochondrial UMIs` = sum(high_Mito_percent),
           Total = sum(low_lib_size | low_n_features | high_Mito_percent))



# adding discard column

colData(sce_qc)$discard <- with(colData(sce_qc), ifelse(low_lib_size == TRUE | low_n_features == TRUE | sum < 5000 & high_Mito_percent == TRUE, TRUE, FALSE))
#df <- data.frame(colData(sce_qc))

# mitochondrial content vs library size
plotColData(sce_qc, 
            x="sum", 
            y="subsets_Mito_percent", 
            colour_by="discard")

# detected genes vs mito content
plotColData(sce_qc, 
            x="detected", 
            y="subsets_Mito_percent", 
            colour_by="discard")

# det genes vs lib size
plotColData(sce_qc, 
            x="sum", 
            y="detected", 
            colour_by="discard")


# filter out poor quality cells
sce.filtered <- sce_qc[,!sce_qc$discard]

saveRDS(sce.filtered, file = here("output", "sce.filtered.RDS"))
saveRDS(sce_qc, file = here("output", "sce_qc.RDS"))

