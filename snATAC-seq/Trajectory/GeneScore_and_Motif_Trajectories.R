# Existing ArchR pseudotime trajectories; no trajectory or direction is fabricated.
source("R/common.R")
suppressPackageStartupMessages({library(ArchR);library(ComplexHeatmap);library(SummarizedExperiment)})
stopifnot(inherits(proj1,"ArchRProject"))
trajectory_name <- NULL  # Replace with an actual trajectory column created via addTrajectory().
if (is.null(trajectory_name)) stop("Set trajectory_name to your existing ArchR trajectory name.")
dir <- output_dir("snATAC_Trajectory")
save_h <- function(h,name) {
  pdf(file.path(dir,paste0(name,".pdf")),width=10,height=12)
  ComplexHeatmap::draw(h);dev.off()
  png(file.path(dir,paste0(name,".png")),width=10,height=12,units="in",res=600)
  ComplexHeatmap::draw(h);dev.off()
}
available <- getAvailableMatrices(proj1)
if ("GeneScoreMatrix" %in% available) {
  trajGene <- getTrajectory(ArchRProj=proj1,name=trajectory_name,
                             useMatrix="GeneScoreMatrix",log2Norm=TRUE)
  h <- plotTrajectoryHeatmap(trajGene,varCutOff=.9,maxFeatures=1000)
  save_h(h,"GeneScore_Trajectory_Heatmap")
  idx <- which(toupper(rownames(trajGene))=="C7")
  if (length(idx)) {
    values <- as.numeric(SummarizedExperiment::assay(trajGene)[idx[1],])
    d <- data.frame(relative_position=seq(0,100,length.out=length(values)),gene_score=values)
    p <- ggplot2::ggplot(d,ggplot2::aes(relative_position,gene_score))+
      ggplot2::geom_line(color="#D95F02")+
      ggplot2::labs(x="Relative trajectory position",y="GeneScore",title="C7 accessibility trajectory")+
      ggplot2::theme_classic()
    save_figure(p,"C7_Accessibility_Trajectory",dir,7,5)
  }
}
if ("MotifMatrix" %in% available) {
  trajMotif <- getTrajectory(ArchRProj=proj1,name=trajectory_name,
                              useMatrix="MotifMatrix",log2Norm=FALSE)
  h <- plotTrajectoryHeatmap(trajMotif,varCutOff=.9,maxFeatures=500)
  save_h(h,"Motif_Trajectory_Heatmap")
}
