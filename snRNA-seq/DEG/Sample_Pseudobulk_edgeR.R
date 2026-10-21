# Sample-level pseudobulk differential expression per cell type.
# REQUIRED: obj metadata columns group, celltype, sample_id; raw RNA counts.
# Biological replicates must be independent samples. Results remain in R memory;
# plots are PDF/PNG only, not CSV.
source("R/common.R")
suppressPackageStartupMessages({library(Seurat);library(edgeR);library(ggplot2);library(Matrix)})
stopifnot(inherits(obj,"Seurat"),all(c("group","celltype","sample_id") %in% colnames(obj[[]])))
DefaultAssay(obj) <- "RNA"
meta <- obj[[]]
meta$analysis_group <- standard_group(meta$group)
meta$sample_id <- as.character(meta$sample_id)
# Avoid silently assigning one sample to both conditions.
stopifnot(all(tapply(as.character(meta$analysis_group),meta$sample_id,
                       function(x) length(unique(x))) == 1))
counts <- Seurat::GetAssayData(obj, assay="RNA", layer="counts")
if (ncol(counts) != nrow(meta)) stop("Count matrix and metadata are not aligned.")
types <- unique(as.character(meta$celltype))
dir <- output_dir("snRNA_DEG")
deg_results <- list()
for (ct in types) {
  idx <- which(as.character(meta$celltype) == ct)
  if (length(idx) < 20) next
  md <- meta[idx,,drop=FALSE]
  groups <- paste(md$sample_id,as.character(md$analysis_group),sep="___")
  f <- factor(groups)
  # Aggregate raw counts by sample and condition, not by individual cell.
  design_sparse <- Matrix::sparse.model.matrix(~0+f)
  pb <- counts[,idx,drop=FALSE] %*% design_sparse
  sample_names <- sub("___.*$","",levels(f))
  condition <- sub("^.*___","",levels(f))
  condition <- factor(condition,levels=c("Healthy_Control","IgA_Nephropathy"))
  if (any(table(condition)<2)) {message("Skip ",ct,": fewer than two samples in a group");next}
  y <- edgeR::DGEList(counts=pb,group=condition)
  keep <- edgeR::filterByExpr(y,group=condition)
  y <- y[keep,,keep.lib.sizes=FALSE]
  if (nrow(y)<10) next
  y <- edgeR::calcNormFactors(y)
  design <- model.matrix(~condition)
  y <- edgeR::estimateDisp(y,design)
  fit <- edgeR::glmQLFit(y,design,robust=TRUE)
  test <- edgeR::glmQLFTest(fit,coef=2)
  tab <- edgeR::topTags(test,n=Inf)$table
  tab$gene <- rownames(tab)
  deg_results[[ct]] <- tab
  d <- data.frame(gene=tab$gene,logFC=tab$logFC,neglogFDR=-log10(pmax(tab$FDR,1e-300)),
                  significant=tab$FDR<0.05 & abs(tab$logFC)>=0.25)
  p <- ggplot(d,aes(logFC,neglogFDR,color=significant))+
    geom_point(size=.8,alpha=.7)+scale_color_manual(values=c("FALSE"="grey75","TRUE"="#D95F02"))+
    labs(title=paste("Differential RNA expression:",ct),x="log2 fold change (IgAN vs HC)",
         y="-log10 FDR",color="FDR < 0.05 and |log2FC| >= 0.25")+theme_classic()
  safe <- gsub("[^A-Za-z0-9_-]","_",ct)
  save_figure(p,paste0("DEG_Volcano_",safe),dir,7,5)
}
message("Pseudobulk DEG results are held in 'deg_results' (not written to disk).")
