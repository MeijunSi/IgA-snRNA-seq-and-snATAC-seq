# IgA Nephropathy | snRNA-seq and snATAC-seq

Modular **English-language R analysis scripts**, organized by task rather than figure number.

## Quick navigation

### snRNA-seq | Seurat `obj`

| Analysis | Script | Requirements |
|---|---|---|
| Visualization: UMAP and composition | [UMAP_and_Composition.R](snRNA-seq/Visualization/UMAP_and_Composition.R) | `obj`, `group`, `celltype`, UMAP |
| Visualization: C7 and complement genes | [C7_and_Complement.R](snRNA-seq/Visualization/C7_and_Complement.R) | RNA assay with normalized expression |
| DEG: pseudobulk edgeR | [Sample_Pseudobulk_edgeR.R](snRNA-seq/DEG/Sample_Pseudobulk_edgeR.R) | **Raw counts**, `sample_id`, `group`, `celltype`, biological replicates |

### snATAC-seq | ArchR `proj1`

| Analysis | Script | Requirements |
|---|---|---|
| Visualization: UMAP and C7 GeneScore | [UMAP_and_C7.R](snATAC-seq/Visualization/UMAP_and_C7.R) | `proj1`, UMAP, GeneScoreMatrix for C7 |
| Visualization: C7 browser tracks | [C7_BrowserTracks.R](snATAC-seq/Visualization/C7_BrowserTracks.R) | ArchR Arrow files, correct reference genome |
| DAR: differential accessible regions | [Differential_Peaks.R](snATAC-seq/DAR/Differential_Peaks.R) | PeakMatrix, TSS enrichment, fragment counts |
| Motif: differential enrichment | [Differential_Motif_Enrichment.R](snATAC-seq/Motif/Differential_Motif_Enrichment.R) | PeakMatrix, motif peak annotations |
| Motif: FLI1 and NR2F1 footprints | [FLI1_NR2F1_Footprints.R](snATAC-seq/Motif/FLI1_NR2F1_Footprints.R) | motif positions and readable insertion files |
| Trajectory: GeneScore and motif activity | [GeneScore_and_Motif_Trajectories.R](snATAC-seq/Trajectory/GeneScore_and_Motif_Trajectories.R) | existing trajectory from `addTrajectory()` |

All scripts reuse [R/common.R](R/common.R) for English group names, colors and figure exports.

## Repository structure

```text
R/
  common.R
snRNA-seq/
  Visualization/
    UMAP_and_Composition.R
    C7_and_Complement.R
  DEG/
    Sample_Pseudobulk_edgeR.R
snATAC-seq/
  Visualization/
    UMAP_and_C7.R
    C7_BrowserTracks.R
  DAR/
    Differential_Peaks.R
  Motif/
    Differential_Motif_Enrichment.R
    FLI1_NR2F1_Footprints.R
  Trajectory/
    GeneScore_and_Motif_Trajectories.R
```

## How to run

From the **repository root** in an R session with `obj` and/or `proj1` already loaded:

```r
source("snRNA-seq/Visualization/UMAP_and_Composition.R")
source("snRNA-seq/Visualization/C7_and_Complement.R")

# Additional metadata required: sample_id
source("snRNA-seq/DEG/Sample_Pseudobulk_edgeR.R")

source("snATAC-seq/Visualization/UMAP_and_C7.R")
source("snATAC-seq/Visualization/C7_BrowserTracks.R")

# Run only after the corresponding ArchR analyses have been prepared:
source("snATAC-seq/DAR/Differential_Peaks.R")
source("snATAC-seq/Motif/Differential_Motif_Enrichment.R")
source("snATAC-seq/Motif/FLI1_NR2F1_Footprints.R")
# First edit trajectory_name in the trajectory script:
source("snATAC-seq/Trajectory/GeneScore_and_Motif_Trajectories.R")
```

Run modules individually, not blindly all at once: motif, DAR and trajectory analyses can be computationally expensive. The trajectory script deliberately stops until you supply an existing trajectory name.

### Metadata and groups

- Group column: `group`
- Cell annotation column: `celltype`
- Biological replicate column for pseudobulk DEG: `sample_id` (edit the script if yours differs).
- Recognized control labels: HC, CTL, Healthy_Control, control.
- Recognized IgAN labels: IgA, IgAN, IgA_Nephropathy, nephritis.

## Figures and analysis objects

Figures are saved only as **PDF and 600-dpi PNG**, under `figures/<analysis>/`. No CSV files are written. DEG results are kept as R object `deg_results`; ArchR results are kept in the current session. Native ArchR footprint PDFs may be stored in the project's own `Plots/` directory.

## Biological and statistical notes

- UMAP and C7 violin plots are descriptive. Do not use individual cells as biological replicates.
- Pooled composition plots are not differential-abundance tests.
- Pseudobulk DEG sums raw RNA counts **by sample within cell type** and uses edgeR. Independent samples, not cells, are the inferential units. Replication, batch design and count availability must be checked before use.
- ATAC GeneScore is a measure of estimated gene accessibility, not RNA expression.
- Motif enrichment is not direct evidence of transcription-factor binding or activity.
- Footprints require motif-specific positions, and trajectory plots require a biologically justified trajectory previously calculated with ArchR.

## Packages

Seurat, ggplot2, dplyr, patchwork, scales, Matrix, edgeR, ArchR, ComplexHeatmap, SummarizedExperiment.

## Earlier monolithic scripts (kept for reference)

- [IgAN_snRNA_Visualization.R](snRNA-seq/IgAN_snRNA_Visualization.R)
- [IgAN_snATAC_Visualization.R](snATAC-seq/IgAN_snATAC_Visualization.R)

**Use the task-oriented scripts above for new analysis.** The earlier combined files are retained so no prior work is lost.

## Data and privacy

Do not commit patient-identifying records, raw reads, Seurat/ArchR serialized projects, Arrow files, credentials or restricted raw data. No license has yet been selected.

## Extended downstream analyses

### snRNA-seq: functional analysis (requires pseudobulk DEG results)

- [GO and KEGG enrichment](snRNA-seq/DEG/GO_KEGG_Enrichment.R) — overrepresentation analysis using a tested-gene background and human gene annotation (`clusterProfiler`, `org.Hs.eg.db`).
- [Hallmark GSEA](snRNA-seq/DEG/GSEA_Hallmark.R) — preranked MSigDB Hallmark enrichment using `fgsea` and `msigdbr`.

Run [Sample_Pseudobulk_edgeR.R](snRNA-seq/DEG/Sample_Pseudobulk_edgeR.R) first in the **same R session**. Both enrichment scripts depend on the in-memory `deg_results` list; neither invents results or exports CSV.

### snATAC-seq: motif activity and TF–gene exploration

- [chromVAR motif activity](snATAC-seq/Motif/chromVAR_Motif_Activity.R) — motif deviation UMAPs, requiring a precomputed `MotifMatrix`.
- [TF–gene integration](snATAC-seq/Integration/TF_Gene_Integration.R) — an **exploratory** correlation of motif deviations against GeneScore aggregated by cell type. Requires matching features, cell barcodes, sufficient cell types, and the motif deviation `z` assay. **This is not peak-to-gene linkage, TF occupancy, or causal evidence.**

### Suggested execution order

```r
# Run from the repository root with obj already loaded.
source("snRNA-seq/DEG/Sample_Pseudobulk_edgeR.R")
source("snRNA-seq/DEG/GO_KEGG_Enrichment.R")
source("snRNA-seq/DEG/GSEA_Hallmark.R")

# Run with proj1 already loaded; only after needed ArchR matrices exist.
source("snATAC-seq/Motif/chromVAR_Motif_Activity.R")
source("snATAC-seq/Integration/TF_Gene_Integration.R")
```

Extra R dependencies: `clusterProfiler`, `org.Hs.eg.db`, `enrichplot`, `msigdbr`, `fgsea`.

**Validation status:** These are reproducible templates committed to the repository, not outputs confirmed by execution on the user's dataset. Check metadata, gene IDs, ArchR feature identifiers and package versions before publication.
