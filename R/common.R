# Shared utilities: source from repository root.
suppressPackageStartupMessages({library(ggplot2)})
output_dir <- function(section) {
  path <- file.path(getwd(), "figures", section)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}
save_figure <- function(plot, name, dir, width = 8, height = 6) {
  grDevices::pdf(file.path(dir, paste0(name, ".pdf")), width = width, height = height, useDingbats = FALSE)
  print(plot); grDevices::dev.off()
  ggplot2::ggsave(file.path(dir, paste0(name, ".png")), plot = plot,
                  width = width, height = height, dpi = 600, limitsize = FALSE)
}
standard_group <- function(x) {
  x0 <- tolower(gsub("[^a-z0-9]", "", as.character(x)))
  result <- ifelse(x0 %in% c("hc","ctl","control","normal","healthy","healthycontrol"), "Healthy_Control",
                   ifelse(x0 %in% c("iga","igan","iganephropathy","nephritis"), "IgA_Nephropathy", NA_character_))
  if (anyNA(result)) stop("Unknown group labels: ", paste(unique(as.character(x)[is.na(result)]), collapse = ", "))
  factor(result, levels = c("Healthy_Control","IgA_Nephropathy"))
}
group_colors <- c(Healthy_Control = "#238B8D", IgA_Nephropathy = "#E66B64")
