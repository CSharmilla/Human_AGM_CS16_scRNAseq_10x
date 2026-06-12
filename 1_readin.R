library(here)
library(DropletUtils)

if (!file.exists(here("output", "sce"))){
  sample <- "Y:/Sharmilla/Edie_AGM_CS16_17/scRNAseq_GSE151876/GSE151876_CellRanger_CD43enriched/outs/filtered_feature_bc_matrix"
  sce <- read10xCounts(samples = sample, col.names = TRUE)
  saveRDS(sce, file = here("output", "sce"))
} else {
  sce <- readRDS(here("output", "sce"))
}

