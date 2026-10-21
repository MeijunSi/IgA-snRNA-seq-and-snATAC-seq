# ArchR UMAP and C7 GeneScore visualization.
source("R/common.R")
suppressPackageStartupMessages({library(ArchR);library(patchwork)})
stopifnot(inherits(proj1,"ArchRProject"))
meta <- as.data.frame(getCellColData(proj1))
stopifnot(all(c("group","celltype") %in% names(meta)))
proj1$plot_group <- as.character(standard_group(meta$group))
dir <- output_dir("snATAC_Visualization")
types <- sort(unique(as.character(meta$celltype)))
pal <- setNames(grDevices::hcl.colors(length(types),"Dynamic"),types)
p1 <- plotEmbedding(proj1,colorBy="cellColData",name="plot_group",
                    embedding="UMAP",pal=group_colors)
p2 <- plotEmbedding(proj1,colorBy="cellColData",name="celltype",
                    embedding="UMAP",pal=pal)
save_figure(p1,"UMAP_ByCondition",dir,7,6)
save_figure(p2,"UMAP_ByCellType",dir,8,6)
save_figure(p1+p2,"UMAP_Combined",dir,14,6)
if ("GeneScoreMatrix" %in% getAvailableMatrices(proj1)) {
  tryCatch({
    p <- plotEmbedding(proj1,colorBy="GeneScoreMatrix",name="C7",
                       embedding="UMAP",continuousSet="solarExtra")
    save_figure(p,"C7_GeneScore_UMAP",dir,7,6)
  },error=function(e) message("C7 GeneScore skipped: ",conditionMessage(e)))
}
# GeneScore is accessibility-based; do not label it RNA expression.
