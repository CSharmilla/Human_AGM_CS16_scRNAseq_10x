library(here)
library(scater)
library(scran)
library(patchwork)
library(dplyr)


sce.clus <- readRDS(here("output", "sce.clus.RDS"))

k <- "cluster_k5" 


rowdata <- data.frame(rowData(sce.clus))

# PTPRC = CD45, SPN = CD43
endothelial_markers <- c("CDH5", "CD34", "SPN", "PTPRC")
plot_list <- lapply(endothelial_markers,
                    function(x)plotUMAP(sce.clus, colour_by = x, text_by = k, text_colour = "red") + ggtitle(x))

patchwork::wrap_plots(plot_list) +  patchwork::plot_annotation(title = "endothelial")

plotUMAP(sce.clus, colour_by = k, point_size=1, text_by = k, text_size =5)
plotUMAP(sce.clus, colour_by = "SPN", point_size=1, text_by = k, text_size =5)
plotUMAP(sce.clus, colour_by = "PTPRC", point_size=1, text_by = k, text_size =5)
plotUMAP(sce.clus, colour_by = "NR2F2", point_size=1, text_by = k, text_size =5)

hsc_markers <- c("RUNX1", "MLLT3", "HOXA9", "MECOM", "HLF", "SPINK2")

plot_list <- lapply(hsc_markers,
                    function(x)plotUMAP(sce.clus, colour_by = x, text_by = k, text_colour = "red") + ggtitle(x))

patchwork::wrap_plots(plot_list) +  patchwork::plot_annotation(title = "hsc")


# get markers
# using scoremarkers function
marker.info <- scoreMarkers(sce.clus, sce.clus$cluster_k5)
marker.info

for (i in 1:16){
  chosen <- marker.info[[i]]
  name <- paste0('cluster', i)
  assign(name, data.frame(chosen))#[,grepl("AUC", colnames(chosen))]))
  name$gene <- row.names(name)
}

# using findmarkers
ribo <- grepl("^RP[SL]", rownames(sce.clus))
mito <- grepl("MT-", rownames(sce.clus))
keep <- !(ribo|mito)
genes.to.use <- sce.clus[keep,]
find.markers.1 <- findMarkers(genes.to.use, groups = genes.to.use$cluster_k5, direction = "up", lfc = 1)
for (i in 1:16){
  chosen <- find.markers.1[[i]]
  name <- paste0('lfc1_fm_cluster', i)
  assign(name, data.frame(chosen))#[,grepl("AUC", colnames(chosen))]))
  name$gene <- row.names(name)
}

plot.list <- lapply(find.markers.1, function(a){a %>% 
                                                as_tibble(rownames = "rownames") %>%
                                                dplyr::filter(FDR<0.05) %>%
                                                arrange(FDR)}) # without as_tibble it throws an error
howmany <- 50
plots <- lapply(names(plot.list), function(a){x <- ifelse(nrow(plot.list[[a]])<howmany, nrow(plot.list[[a]]), howmany)
                                      plotDots(sce.clus, plot.list[[a]]$rownames[1:x], group = "cluster_k5") + 
                                        ggtitle(paste0("Top ", howmany, " genes (lfc>1) for cluster ", a))})
print(plots)

find.markers <- findMarkers(genes.to.use, groups = genes.to.use$cluster_k5, direction = "up")
for (i in 1:16){
  chosen <- find.markers[[i]]
  name <- paste0('fm_cluster', i)
  assign(name, data.frame(chosen))#[,grepl("AUC", colnames(chosen))]))
  name$gene <- row.names(name)
}

