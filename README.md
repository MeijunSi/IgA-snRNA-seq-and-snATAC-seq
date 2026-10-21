# IgA Nephropathy: snRNA-seq and snATAC-seq Visualization

Reproducible R visualization scripts for human IgA nephropathy (IgAN) kidney multiomics.

## Scripts

- [snRNA-seq/IgAN_snRNA_Visualization.R](snRNA-seq/IgAN_snRNA_Visualization.R): Seurat cell-type and condition UMAPs, group composition, C7 feature/violin plots, and complement-associated expression dot plot.
- [snATAC-seq/IgAN_snATAC_Visualization.R](snATAC-seq/IgAN_snATAC_Visualization.R): ArchR UMAPs, C7 GeneScore embedding, browser tracks, differential motif enrichment, and optional motif footprints and trajectories.

## Input objects

The scripts expect objects **already loaded into R**:

- `obj`: Seurat object with metadata columns `group` and `celltype`, RNA assay, and `umap` or `umap.harmony`.
- `proj1`: ArchRProject with metadata columns `group` and `celltype`, UMAP, and Arrow files accessible at their original paths.

The group labels `HC`, `CTL`, `Healthy_Control`, `IgA`, `IgAN`, and `IgA_Nephropathy` are standardized for presentation. Modify the label map if your dataset uses other names. Input matrices are **not included**.

## Running

Install compatible Bioconductor and CRAN packages, then run within the R session where the corresponding object is loaded:

```r
source("snRNA-seq/IgAN_snRNA_Visualization.R")
source("snATAC-seq/IgAN_snATAC_Visualization.R")
```

The scripts produce PDF (vector where feasible) and 600-dpi PNG plots under `IgAN_snRNA_Figures/` or `IgAN_snATAC_Figures/` in the current working directory. **No CSV files are written.**

## Dependencies

- R, Seurat (v5 recommended), ggplot2, dplyr, patchwork, scales
- ArchR, ComplexHeatmap and SummarizedExperiment for snATAC
- A functional Cairo PDF graphics device for some systems (base PDF is used as fallback)

## Interpretation and requirements

- Composition bars summarize pooled cell counts by group, **not** sample-level differential-abundance tests.
- Violin plots show cell-level RNA expression and should not be treated as independent biological replicates. For statistical comparisons, analyze sample-level pseudobulk.
- C7 GeneScore reflects chromatin accessibility, **not** C7 RNA abundance.
- Browser plots require correct genome annotations, genomic coordinates and readable Arrow files.
- Motif enrichment requires PeakMatrix, peak annotations, and bias covariates. Footprints require the matching motif positions. Trajectory figures require an existing ArchR trajectory annotation; no pseudotime is invented.
- Some plots may be skipped with an explanatory message when prerequisites are absent. No analysis outcome is asserted without data.

## Reproducibility and sharing

Do not commit patient/sample identifying metadata, raw sequencing reads, large ArchR Arrow files, serialized Seurat/ArchR objects, or credentials. The included `.gitignore` excludes common generated and data files. Record package versions using `sessionInfo()` in your local analysis record.

## License

No license is specified here. Add one after choosing the terms under which you want the code reused.
