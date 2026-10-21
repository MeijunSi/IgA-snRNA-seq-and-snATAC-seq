# chromVAR motif deviation UMAPs, not motif enrichment or TF expression.
# Requires MotifMatrix, with deviation scores in the 'z' assay.
source("R/common.R")
suppressPackageStartupMessages({library(ArchR);library(ggplot2)})
stopifnot(inherits(proj1,"ArchRProject"))
stopifnot("MotifMatrix" %in% getAvailableMatrices(proj1))
dir <- output_dir("snATAC_Motif")
motif_features <- getFeatures(proj1,useMatrix="MotifMatrix")
targets <- c("ELK4","FLI1","CREM","MEF2A","NR2F1","NR3C1","TCF21")
motif_selected <- character()
for (tf in targets) {
  hits <- grep(paste0("(^|:)",tf,"(_|$)"),motif_features,value=TRUE,ignore.case=TRUE)
  if (!length(hits)) hits <- grep(paste0("^",tf,"_"),motif_features,value=TRUE,ignore.case=TRUE)
  if (length(hits)) motif_selected <- c(motif_selected,hits[1]) else message("Motif not found: ",tf)
}
for (motif in unique(motif_selected)) {
  tryCatch({
    p <- plotEmbedding(ArchRProj=proj1,colorBy="MotifMatrix",name=motif,
                       embedding="UMAP",imputeWeights=NULL,continuousSet="solarExtra")
    save_figure(p,paste0("chromVAR_",gsub("[^A-Za-z0-9_-]","_",motif)),dir,7,6)
  },error=function(e)message("Motif plot skipped for ",motif,": ",conditionMessage(e)))
}
# Interpretation: chromVAR deviation scores quantify motif accessibility.
