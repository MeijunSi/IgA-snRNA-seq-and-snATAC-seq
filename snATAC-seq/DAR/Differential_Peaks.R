# Differential accessible regions by group (ArchR Wilcoxon with bias correction).
# Requires PeakMatrix, group, TSSEnrichment and nFrags.
source("R/common.R")
suppressPackageStartupMessages({library(ArchR);library(ggplot2)})
stopifnot(inherits(proj1,"ArchRProject"))
stopifnot("PeakMatrix" %in% getAvailableMatrices(proj1))
meta <- as.data.frame(getCellColData(proj1))
stopifnot(all(c("group","TSSEnrichment","nFrags") %in% names(meta)))
proj1$plot_group <- as.character(standard_group(meta$group))
dir <- output_dir("snATAC_DAR")
marker_peaks <- getMarkerFeatures(ArchRProj=proj1,useMatrix="PeakMatrix",
                                  groupBy="plot_group",bias=c("TSSEnrichment","log10(nFrags)"),
                                  testMethod="wilcoxon")
# Signed log2FC direction depends on marker group. Keep original ArchR objects
# in memory to avoid discarding the associated peak annotations.
plot_marker <- tryCatch(
  plotMarkers(seMarker=marker_peaks,name="IgA_Nephropathy",cutOff="FDR <= 0.05 & Log2FC >= 0.5"),
  error=function(e) {message(conditionMessage(e));NULL}
)
if (inherits(plot_marker,"ggplot")) save_figure(plot_marker,"IgAN_Differential_Peaks",dir,7,6)
message("Result object: marker_peaks. No CSV generated.")
