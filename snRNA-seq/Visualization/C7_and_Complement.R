# C7 gene expression by condition and cell type.
source("R/common.R")
suppressPackageStartupMessages({library(Seurat);library(ggplot2);library(patchwork)})
stopifnot(inherits(obj,"Seurat"),all(c("group","celltype") %in% colnames(obj[[]])))
stopifnot("RNA" %in% Assays(obj))
dir <- output_dir("snRNA_Visualization")
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]],"Assay5")) obj <- JoinLayers(obj,assay="RNA")
obj$plot_group <- standard_group(obj$group)
obj$plot_celltype <- factor(as.character(obj$celltype))
red <- if ("umap.harmony" %in% Reductions(obj)) "umap.harmony" else "umap"
stopifnot(red %in% Reductions(obj))
if ("C7" %in% rownames(obj[["RNA"]])) {
  p <- FeaturePlot(obj,features="C7",reduction=red,split.by="plot_group",
                   cols=c("#EEEEEE","#FDAE91","#CB181D"),order=TRUE,
                   raster=TRUE,pt.size=.08) & theme_classic()
  save_figure(p,"C7_FeaturePlot",dir,11,5)
  v <- VlnPlot(obj,features="C7",group.by="plot_celltype",split.by="plot_group",
               split.plot=FALSE,cols=group_colors,pt.size=0)+
    labs(x="Cell type",y="Normalized RNA expression")+theme_classic()+
    theme(axis.text.x=element_text(angle=45,hjust=1))
  save_figure(v,"C7_Violin",dir,13,5)
} else message("C7 not present; expression figures skipped.")
genes <- intersect(c("C3","C4BPA","C5","C6","C7","C8A","C8B","C8G","C9",
                     "CFH","CFI","CFHR1","CFHR2","CFHR5","C1QA","C1QB"),rownames(obj[["RNA"]]))
if (length(genes)) {
  p <- DotPlot(obj,features=genes,group.by="plot_celltype",assay="RNA",
               scale=TRUE,dot.scale=6)+
    labs(x="Gene",y="Cell type",color="Scaled average expression",size="Percent expressed")+
    theme_classic()+theme(axis.text.x=element_text(angle=45,hjust=1))
  save_figure(p,"Complement_Genes_DotPlot",dir,12,7)
}
# Cell-level violins are descriptive; use sample pseudobulk for inference.
