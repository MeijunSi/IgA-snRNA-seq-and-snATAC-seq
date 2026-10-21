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
