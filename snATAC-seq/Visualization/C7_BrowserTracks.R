# Genome browser: C7 accessibility grouped by condition.
source("R/common.R")
suppressPackageStartupMessages(library(ArchR))
stopifnot(inherits(proj1,"ArchRProject"))
meta <- as.data.frame(getCellColData(proj1))
stopifnot("group" %in% colnames(meta))
proj1$plot_group <- as.character(standard_group(meta$group))
dir <- output_dir("snATAC_Visualization")
tracks <- plotBrowserTrack(ArchRProj=proj1,groupBy="plot_group",
                           geneSymbol="C7",upstream=50000,downstream=50000,
                           plotSummary=c("bulkTrack","geneTrack"))
if (!is.list(tracks)) tracks <- list(tracks)
for (i in seq_along(tracks))
  save_figure(tracks[[i]],paste0("C7_BrowserTrack_",i),dir,12,5)
