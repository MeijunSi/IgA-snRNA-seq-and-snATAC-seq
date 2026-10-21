# Preranked Hallmark GSEA: edgeR pseudobulk results, per cell type.
# Requires msigdbr, fgsea and the pseudobulk object 'deg_results'.
source("R/common.R")
suppressPackageStartupMessages({library(msigdbr);library(fgsea);library(ggplot2)})
stopifnot(exists("deg_results"),length(deg_results)>0)
dir <- output_dir("snRNA_GSEA")
hallmark <- msigdbr::msigdbr(species="Homo sapiens",collection="H")
pathways <- split(hallmark$gene_symbol,hallmark$gs_name)
gsea_results <- list()
for (ct in names(deg_results)) {
  d <- deg_results[[ct]]
  if (!all(c("gene","logFC","PValue") %in% names(d))) next
  d <- d[is.finite(d$logFC)&is.finite(d$PValue)&!is.na(d$gene),]
  # Signed Wald-like ranking; not interchangeable with the edgeR QL F statistic.
  d$rank <- sign(d$logFC)*(-log10(pmax(d$PValue,1e-300)))
  d <- d[order(-abs(d$rank)),]
  d <- d[!duplicated(d$gene),]
  ranks <- setNames(d$rank,d$gene)
  ranks <- sort(ranks,decreasing=TRUE)
  if (length(ranks)<100) next
  res <- fgsea::fgseaMultilevel(pathways=pathways,stats=ranks,minSize=15,maxSize=500)
  gsea_results[[ct]] <- res
  res <- as.data.frame(res)
  res <- res[!is.na(res$padj)&res$padj<0.05,]
  if (!nrow(res)) next
  res <- res[order(res$padj),]
  res <- head(res,15)
  res$pathway <- factor(res$pathway,levels=rev(res$pathway))
  p <- ggplot(res,aes(x=NES,y=pathway,color=padj))+
    geom_point(size=3)+scale_color_gradient(low="#C43C39",high="#3B80AC")+
    labs(title=paste("Hallmark GSEA:",ct),x="Normalized enrichment score",y=NULL,color="FDR")+
    theme_classic(base_size=11)
  save_figure(p,paste0("GSEA_Hallmark_",gsub("[^A-Za-z0-9_-]","_",ct)),dir,10,6)
}
message("Full GSEA results are retained in 'gsea_results'; no CSV output.")
