# UMAP and pooled cell composition. Run from repository root with 'obj' in memory.
source("R/common.R")
suppressPackageStartupMessages({library(Seurat);library(dplyr);library(patchwork);library(scales)})
stopifnot(inherits(obj,"Seurat"), all(c("group","celltype") %in% colnames(obj[[]])))
dir <- output_dir("snRNA_Visualization")
obj$plot_group <- standard_group(obj$group)
obj$plot_celltype <- factor(as.character(obj$celltype))
reduction <- if ("umap.harmony" %in% Reductions(obj)) "umap.harmony" else "umap"
stopifnot(reduction %in% Reductions(obj))
types <- levels(obj$plot_celltype)
cols <- setNames(grDevices::hcl.colors(length(types),"Dynamic"),types)
p <- DimPlot(obj,reduction=reduction,group.by="plot_celltype",split.by="plot_group",cols=cols,
             raster=TRUE,pt.size=0.05,ncol=2) & theme_classic(base_size=11)
save_figure(p,"UMAP_CellTypes_ByCondition",dir,12,5)
p_group <- DimPlot(obj,reduction=reduction,group.by="plot_group",
                   cols=group_colors,pt.size=0.05,raster=TRUE) + theme_classic()
save_figure(p_group,"UMAP_ByCondition",dir,7,6)
d <- obj[[]] %>% dplyr::count(plot_group,plot_celltype,name="n") %>%
  dplyr::group_by(plot_group) %>% dplyr::mutate(proportion=n/sum(n)) %>% dplyr::ungroup()
p_comp <- ggplot(d,aes(x=proportion,y=plot_group,fill=plot_celltype))+
  geom_col(width=.6)+scale_fill_manual(values=cols)+
  scale_x_continuous(labels=scales::percent,expand=c(0,0))+
  labs(x="Cell proportion",y=NULL,fill="Cell type")+theme_classic()
save_figure(p_comp,"Pooled_Cell_Composition",dir,10,4)
message("Pooled proportions are descriptive, not sample-level differential abundance.")
