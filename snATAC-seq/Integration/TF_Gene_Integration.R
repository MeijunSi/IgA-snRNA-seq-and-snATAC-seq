# ArchR TF-to-gene correlation using motif deviation and GeneScore.
# Requires MotifMatrix, GeneScoreMatrix, and aligned cells in the ArchR project.
# Correlations are exploratory and do not imply direct TF binding.
source("R/common.R")
suppressPackageStartupMessages({library(ArchR);library(ggplot2);library(SummarizedExperiment)})
stopifnot(inherits(proj1,"ArchRProject"))
stopifnot(all(c("MotifMatrix","GeneScoreMatrix") %in% getAvailableMatrices(proj1)))
dir <- output_dir("snATAC_Integration")
tfs <- c("ELK4","FLI1","NR2F1","TCF21")
genes <- c("C7","CFH","C3","VCAM1")
motif_names <- getFeatures(proj1,useMatrix="MotifMatrix")
gene_names <- getFeatures(proj1,useMatrix="GeneScoreMatrix")
matched_motifs <- unique(unlist(lapply(tfs,function(tf)head(grep(paste0("^",tf,"_"),motif_names,value=TRUE),1))))
matched_genes <- intersect(genes,gene_names)
if (!length(matched_motifs)||!length(matched_genes)) stop("Requested motif/gene features were not found.")
# Use pseudobulk by cell type for robust exploratory correlations,
# rather than treating per-cell sparse deviations as independent replicates.
meta <- as.data.frame(getCellColData(proj1))
stopifnot("celltype" %in% names(meta))
motif_se <- getMatrixFromProject(proj1,useMatrix="MotifMatrix")
gene_se <- getMatrixFromProject(proj1,useMatrix="GeneScoreMatrix")
m_assay <- if ("z" %in% SummarizedExperiment::assayNames(motif_se)) "z" else
  stop("MotifMatrix deviation z assay is missing; calculate motif deviations first.")
m <- SummarizedExperiment::assay(motif_se,m_assay)
g <- SummarizedExperiment::assay(gene_se)
m_idx <- match(matched_motifs,rownames(m))
g_idx <- match(matched_genes,rownames(g))
common <- Reduce(intersect,list(colnames(m),colnames(g),rownames(meta)))
if (length(common)<20) stop("Insufficient shared cells.")
ct <- as.character(meta[common,"celltype"])
summarize_by_type <- function(mat,idx) {
  sapply(split(seq_along(common),ct),function(ii) Matrix::rowMeans(mat[idx,common[ii],drop=FALSE]))
}
ms <- summarize_by_type(m,m_idx)
gs <- summarize_by_type(g,g_idx)
if (ncol(ms)<4) stop("Too few cell types for robust descriptive correlations.")
d <- expand.grid(TF=matched_motifs,Gene=matched_genes,stringsAsFactors=FALSE)
d$rho <- mapply(function(tf,gene) {
  cor(as.numeric(ms[tf,]),as.numeric(gs[gene,]),method="spearman")
},d$TF,d$Gene)
p <- ggplot(d,aes(Gene,TF,fill=rho))+geom_tile(color="white")+
  geom_text(aes(label=sprintf("%.2f",rho)),size=3)+
  scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",limits=c(-1,1))+
  labs(title="Exploratory motif–gene association",x="GeneScore",y="Motif deviation",fill="Spearman rho")+
  theme_minimal(base_size=11)
save_figure(p,"TF_Gene_Correlation_Heatmap",dir,8,6)
# NOT peak-to-gene linkage and NOT TF binding inference.
