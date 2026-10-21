# GO Biological Process and KEGG enrichment from sample-level pseudobulk DEG
# Prerequisite: run snRNA-seq/DEG/Sample_Pseudobulk_edgeR.R
# Human gene SYMBOL annotation. No CSV output.
source("R/common.R")
suppressPackageStartupMessages({
  library(clusterProfiler);library(org.Hs.eg.db);library(enrichplot);library(ggplot2)
})
stopifnot(exists("deg_results"),length(deg_results)>0)
dir <- output_dir("snRNA_Enrichment")
for (ct in names(deg_results)) {
  dat <- deg_results[[ct]]
  req <- c("gene","FDR","logFC")
  if (!all(req %in% names(dat))) next
  genes <- unique(dat$gene[!is.na(dat$FDR)&dat$FDR<0.05&abs(dat$logFC)>=0.25])
  universe <- unique(dat$gene[!is.na(dat$FDR)])
  if (length(genes)<5) next
  target <- suppressMessages(bitr(genes,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db))
  bg <- suppressMessages(bitr(universe,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db))
  if (is.null(target)||nrow(target)<5||is.null(bg)||nrow(bg)<5) next
  safe <- gsub("[^A-Za-z0-9_-]","_",ct)
  run_one <- function(expr,name) {
    tryCatch({
      x <- expr
      if (!is.null(x)&&nrow(as.data.frame(x))>0) {
        p <- enrichplot::dotplot(x,showCategory=15)+ggtitle(paste(name,ct))
        save_figure(p,paste0(safe,"_",name),dir,10,7)
      }
      invisible(x)
    },error=function(e) message("Skipped ",ct," ",name,": ",conditionMessage(e)))
  }
  run_one(enrichGO(gene=unique(target$ENTREZID),universe=unique(bg$ENTREZID),
                   OrgDb=org.Hs.eg.db,ont="BP",pAdjustMethod="BH",readable=TRUE),"GO_BP")
  run_one(enrichKEGG(gene=unique(target$ENTREZID),universe=unique(bg$ENTREZID),
                     organism="hsa",pAdjustMethod="BH"),"KEGG")
}
