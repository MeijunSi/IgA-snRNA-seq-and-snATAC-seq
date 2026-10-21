############################################################
# IgAN snATAC-seq visualization
# Input: ArchRProject 'proj1'; metadata 'group', 'celltype'
# Output: PDF and 600-dpi PNG only (no CSV)
# Optional: set trajectory_name to an existing trajectory
############################################################
suppressPackageStartupMessages({
  library(ArchR)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(SummarizedExperiment)
})
stopifnot(inherits(proj1, "ArchRProject"))
outdir <- file.path(getwd(), "IgAN_snATAC_Figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
trajectory_name <- NULL
meta <- as.data.frame(getCellColData(proj1))
stopifnot(all(c("group", "celltype") %in% colnames(meta)))
cleaned <- tolower(gsub("[^a-z0-9]", "", as.character(meta$group)))
groups <- ifelse(cleaned %in% c("hc", "ctl", "control", "healthy", "healthycontrol", "normal"),
                 "Healthy_Control",
                 ifelse(cleaned %in% c("iga", "igan", "iganephropathy", "nephritis"),
                        "IgA_Nephropathy", NA_character_))
if (anyNA(groups)) stop("Unrecognized group labels: ",
                        paste(unique(as.character(meta$group)[is.na(groups)]), collapse = ", "))
proj1$plot_group <- groups
group_colors <- c(Healthy_Control = "#238B8D", IgA_Nephropathy = "#E66B64")
types <- sort(unique(as.character(meta$celltype)))
type_colors <- setNames(grDevices::hcl.colors(length(types), "Dynamic"), types)
save_plot <- function(p, name, width = 8, height = 6) {
  ggsave(file.path(outdir, paste0(name, ".pdf")), p, width = width, height = height,
         device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf", limitsize = FALSE)
  ggsave(file.path(outdir, paste0(name, ".png")), p, width = width, height = height,
         dpi = 600, limitsize = FALSE)
}
save_heatmap <- function(h, name, width = 10, height = 9) {
  grDevices::pdf(file.path(outdir, paste0(name, ".pdf")), width = width, height = height)
  ComplexHeatmap::draw(h); grDevices::dev.off()
  grDevices::png(file.path(outdir, paste0(name, ".png")),
                 width = width, height = height, units = "in", res = 600)
  ComplexHeatmap::draw(h); grDevices::dev.off()
}
try_analysis <- function(label, expr) {
  tryCatch(force(expr), error = function(e) {
    message("Skipped ", label, ": ", conditionMessage(e))
    NULL
  })
}
matrices <- getAvailableMatrices(proj1)
message("Available matrices: ", paste(matrices, collapse = ", "))
p_group <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData",
                         name = "plot_group", embedding = "UMAP",
                         pal = group_colors, size = 0.15) + ggtitle("Condition")
p_type <- plotEmbedding(ArchRProj = proj1, colorBy = "cellColData",
                        name = "celltype", embedding = "UMAP",
                        pal = type_colors, size = 0.15) + ggtitle("Cell type")
save_plot(p_group, "UMAP_ByCondition", 7, 6)
save_plot(p_type, "UMAP_ByCellType", 9, 7)
save_plot(p_group + p_type + plot_layout(ncol = 2), "UMAP_Condition_CellType", 14, 6)
if ("GeneScoreMatrix" %in% matrices) {
  try_analysis("C7 GeneScore UMAP", {
    p <- plotEmbedding(ArchRProj = proj1, colorBy = "GeneScoreMatrix",
                       name = "C7", embedding = "UMAP", continuousSet = "solarExtra",
                       size = 0.15) + ggtitle("C7 chromatin accessibility (GeneScore)")
    save_plot(p, "C7_GeneScore_UMAP", 7, 6)
  })
}
# ArchR browser tracks are generated from a real project and correct reference genome.
try_analysis("C7 browser track", {
  plots <- plotBrowserTrack(ArchRProj = proj1, groupBy = "plot_group",
                            geneSymbol = "C7", upstream = 50000, downstream = 50000,
                            plotSummary = c("bulkTrack", "geneTrack"))
  if (!is.list(plots)) plots <- list(plots)
  for (i in seq_along(plots)) {
    save_plot(plots[[i]], paste0("C7_BrowserTrack_", i), 12, 5)
  }
})
has_motifs <- !is.null(tryCatch(getPositions(proj1, name = "Motif"),
                               error = function(e) NULL))
if ("PeakMatrix" %in% matrices && has_motifs) {
  try_analysis("differential motif enrichment", {
    markers <- getMarkerFeatures(ArchRProj = proj1, useMatrix = "PeakMatrix",
                                 groupBy = "plot_group",
                                 bias = c("TSSEnrichment", "log10(nFrags)"),
                                 testMethod = "wilcoxon")
    enrichment <- peakAnnoEnrichment(seMarker = markers, ArchRProj = proj1,
                                     peakAnnotation = "Motif",
                                     cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
    p <- plotEnrichHeatmap(seEnrich = enrichment, n = 40, transpose = TRUE)
    save_heatmap(p, "Differential_Motif_Enrichment", 10, 10)
  })
}
if (has_motifs) {
  try_analysis("FLI1/NR2F1 footprints", {
    positions <- getPositions(proj1, name = "Motif")
    matches <- grep("^(FLI1|NR2F1)(_|$)", names(positions),
                    ignore.case = TRUE, value = TRUE)
    if (!length(matches)) stop("No FLI1 or NR2F1 motif positions found.")
    se_foot <- getFootprints(ArchRProj = proj1, positions = positions[matches],
                             groupBy = "plot_group", flank = 250, minCells = 25)
    # With ArchRProj supplied, plotFootprints writes a PDF to the ArchR project plots directory.
    plotFootprints(seFoot = se_foot, names = matches, normMethod = "Subtract",
                   ArchRProj = proj1, plotName = "FLI1_NR2F1_Footprints",
                   addDOC = FALSE)
  })
}
# No trajectory is inferred: supply a trajectory column generated using addTrajectory().
if (!is.null(trajectory_name)) {
  if ("GeneScoreMatrix" %in% matrices) {
    try_analysis("gene accessibility trajectory", {
      traj <- getTrajectory(ArchRProj = proj1, name = trajectory_name,
                            useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
      h <- plotTrajectoryHeatmap(traj, varCutOff = 0.9, maxFeatures = 1000)
      save_heatmap(h, "GeneScore_Pseudotime_Heatmap", 10, 12)
      idx <- which(toupper(rownames(traj)) == "C7")
      if (length(idx)) {
        values <- as.numeric(SummarizedExperiment::assay(traj)[idx[1], ])
        d <- data.frame(Position = seq(0, 100, length.out = length(values)),
                        Accessibility = values)
        p <- ggplot(d, aes(Position, Accessibility)) +
          geom_line(color = "#D95F02", linewidth = 1) +
          labs(title = "C7 accessibility trajectory",
               x = "Relative trajectory position", y = "GeneScore")
        save_plot(p, "C7_Accessibility_Trajectory", 7, 5)
      }
    })
  }
  if ("MotifMatrix" %in% matrices) {
    try_analysis("motif trajectory", {
      traj <- getTrajectory(ArchRProj = proj1, name = trajectory_name,
                            useMatrix = "MotifMatrix", log2Norm = FALSE)
      h <- plotTrajectoryHeatmap(traj, varCutOff = 0.9, maxFeatures = 500)
      save_heatmap(h, "Motif_Pseudotime_Heatmap", 10, 12)
    })
  }
} else {
  message("Trajectories skipped: set 'trajectory_name' to an existing ArchR trajectory column.")
}
message("snATAC visualization complete: ", outdir)
