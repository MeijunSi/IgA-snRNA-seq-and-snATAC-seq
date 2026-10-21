# Tn5 footprint profiles for selected FLI1 and NR2F1 motifs.
# Requires existing motif positions; motifs are not interchangeable across databases.
source("R/common.R")
suppressPackageStartupMessages(library(ArchR))
stopifnot(inherits(proj1,"ArchRProject"))
meta <- as.data.frame(getCellColData(proj1))
proj1$plot_group <- as.character(standard_group(meta$group))
positions <- getPositions(proj1,name="Motif")
matches <- grep("^(FLI1|NR2F1)(_|$)",names(positions),value=TRUE,ignore.case=TRUE)
if (!length(matches)) stop("FLI1/NR2F1 motif identifiers not found in this project.")
message("Selected motif IDs: ",paste(matches,collapse=", "))
seFoot <- getFootprints(ArchRProj=proj1,positions=positions[matches],
                         groupBy="plot_group",flank=250,minCells=25)
# ArchR writes its native footprint PDF under the project's Plots/ folder.
plotFootprints(seFoot=seFoot,names=matches,normMethod="Subtract",
               ArchRProj=proj1,plotName="FLI1_NR2F1_Footprints",addDOC=FALSE)
