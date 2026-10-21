############################################################
# IgAN snRNA-seq visualization
# Input: Seurat object 'obj'; metadata 'group', 'celltype'
# Output: PDF and 600-dpi PNG only (no CSV)
############################################################
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})
stopifnot(inherits(obj, "Seurat"))
stopifnot(all(c("group", "celltype") %in% colnames(obj[[]])))
stopifnot("RNA" %in% Assays(obj))
outdir <- file.path(getwd(), "IgAN_snRNA_Figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj <- Seurat::JoinLayers(obj, assay = "RNA")
}
standardize_group <- function(x) {
  cleaned <- tolower(gsub("[^a-z0-9]", "", as.character(x)))
  result <- ifelse(cleaned %in% c("hc", "ctl", "control", "healthy", "healthycontrol", "normal"), "Healthy_Control",
                   ifelse(cleaned %in% c("iga", "igan", "iganephropathy", "nephritis"), "IgA_Nephropathy", NA_character_))
  if (anyNA(result)) stop("Unknown group labels: ", paste(unique(as.character(x)[is.na(result)]), collapse = ", "))
  factor(result, levels = c("Healthy_Control", "IgA_Nephropathy"))
}
obj$plot_group <- standardize_group(obj$group)
obj$plot_celltype <- factor(as.character(obj$celltype))
group_colors <- c(Healthy_Control = "#238B8D", IgA_Nephropathy = "#E66B64")
base_colors <- c("#8DA0CB", "#FC8D62", "#66C2A5", "#8CC63F", "#D8C57A", "#C9B458",
                 "#7B2CBF", "#B9B8D8", "#5DA5DA", "#FFA500", "#DDA0DD", "#F4A3A8",
                 "#FF69B4", "#2ED9C3", "#A6761D", "#7570B3", "#E7298A", "#1B9E77")
types <- levels(obj$plot_celltype)
type_colors <- setNames(grDevices::colorRampPalette(base_colors)(length(types)), types)
theme_publication <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text = element_text(color = "black"), legend.title = element_text(face = "bold"))
save_plot <- function(p, name, width, height) {
  pdf_path <- file.path(outdir, paste0(name, ".pdf"))
  png_path <- file.path(outdir, paste0(name, ".png"))
  ggsave(pdf_path, p, width = width, height = height, limitsize = FALSE,
         device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf")
  ggsave(png_path, p, width = width, height = height, dpi = 600, limitsize = FALSE)
}
reduction <- if ("umap.harmony" %in% Reductions(obj)) "umap.harmony" else "umap"
if (!reduction %in% Reductions(obj)) stop("Missing UMAP reduction.")
p_types <- DimPlot(obj, reduction = reduction, group.by = "plot_celltype",
                   split.by = "plot_group", cols = type_colors, pt.size = 0.05,
                   raster = TRUE, ncol = 2) &
  theme_publication & theme(axis.text = element_blank(), axis.ticks = element_blank(),
                            axis.title = element_blank(), legend.title = element_blank())
print(p_types)
save_plot(p_types, "UMAP_CellTypes_ByCondition", 12, 5.5)
p_conditions <- DimPlot(obj, reduction = reduction, group.by = "plot_group",
                        cols = group_colors, pt.size = 0.05, raster = TRUE) +
  ggtitle("Condition") + theme_publication
print(p_conditions)
save_plot(p_conditions, "UMAP_ByCondition", 7, 6)
composition <- obj[[]] %>%
  dplyr::count(plot_group, plot_celltype, name = "Cells") %>%
  dplyr::group_by(plot_group) %>%
  dplyr::mutate(Proportion = Cells / sum(Cells)) %>%
  dplyr::ungroup()
p_composition <- ggplot(composition, aes(x = Proportion, y = plot_group, fill = plot_celltype)) +
  geom_col(width = 0.55, color = "black", linewidth = 0.15) +
  scale_fill_manual(values = type_colors) +
  scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "Cell proportion", y = NULL, fill = "Cell type") + theme_publication
print(p_composition)
save_plot(p_composition, "CellType_Composition", 10, 4.5)
if ("C7" %in% rownames(obj[["RNA"]])) {
  p_c7 <- FeaturePlot(obj, features = "C7", reduction = reduction, split.by = "plot_group",
                      order = TRUE, min.cutoff = "q05", max.cutoff = "q99",
                      cols = c("#EEEEEE", "#FDAE91", "#CB181D"),
                      pt.size = 0.08, raster = TRUE, ncol = 2) &
    theme_publication & theme(axis.text = element_blank(), axis.ticks = element_blank(),
                              axis.title = element_blank())
  print(p_c7)
  save_plot(p_c7, "C7_FeaturePlot_ByCondition", 11, 5.5)
  # split.plot = FALSE displays side-by-side violins for both conditions.
  p_violin <- VlnPlot(obj, features = "C7", group.by = "plot_celltype",
                      split.by = "plot_group", split.plot = FALSE,
                      cols = group_colors, pt.size = 0) +
    labs(title = "C7", x = "Cell type", y = "Normalized RNA expression") +
    theme_publication +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p_violin)
  save_plot(p_violin, "C7_Expression_Violin", 14, 5)
} else {
  message("C7 is absent from the RNA assay; C7 RNA plots skipped.")
}
complement_genes <- c("C4BPA", "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9",
                      "CFHR1", "CFHR2", "CFHR5", "CFI", "CFH", "C3", "C1QA", "C1QB")
available_genes <- intersect(complement_genes, rownames(obj[["RNA"]]))
if (length(available_genes) > 0) {
  p_dot <- DotPlot(obj, features = available_genes, group.by = "plot_celltype",
                   assay = "RNA", scale = TRUE, dot.scale = 6) +
    scale_color_gradient2(low = "#377EB8", mid = "white", high = "#F46D43", midpoint = 0) +
    labs(title = "Complement-related gene expression", x = "Gene", y = "Cell type",
         color = "Scaled average expression", size = "Percent expressed") +
    theme_publication + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p_dot)
  save_plot(p_dot, "Complement_Genes_DotPlot", 12, 7)
}
p_overview <- p_types / p_composition + plot_layout(heights = c(2, 1))
print(p_overview)
save_plot(p_overview, "snRNA_Overview", 13, 10)
message("snRNA visualization complete: ", outdir)
