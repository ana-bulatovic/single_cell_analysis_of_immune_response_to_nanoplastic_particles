# Single-Cell Analysis of Immune Response to Nanoplastic Particles

Reproducible single-cell RNA-seq (scRNA-seq) workflow for studying how human peripheral blood mononuclear cells (PBMC) respond to carboxylated polystyrene nanoparticles (PSNP) of different sizes.

**Data:** [Zenodo 15866724](https://zenodo.org/records/15866724) · **Reference atlas:** [Zenodo 4546839](https://zenodo.org/records/4546839) (`ref.Rds`)

---

## Overview

This repository implements an end-to-end analysis of a four-condition PBMC experiment:

| Sample | Condition |
|--------|-----------|
| Sample 1 | **40 nm** PSNP |
| Sample 2 | **200 nm** PSNP |
| Sample 3 | **40 + 200 nm** mixture |
| Sample 4 | **Untreated control** |

The pipeline covers quality control, batch integration (Combat), clustering, cell-type annotation, composition analysis, differential expression (DE), pathway enrichment, **size-specific effect classification**, and supplementary analyses (module scores, pseudobulk export, annotation validation).

---

## Features

- **QC & integration** — Scanpy-based filtering, normalization, HVG selection, Combat integration, UMAP, Leiden clustering
- **Multi-method annotation** — Seurat/Azimuth (`ref.Rds`), CoDi labels, optional Tangram / cell2location, literature marker validation
- **Differential expression** — Wilcoxon tests per cell type (each exposure vs control)
- **Pathway enrichment** — Enrichr (GO, KEGG, Reactome) on **UP- and DOWN-regulated** DE genes vs control, with summary tables and biological interpretation
- **Annotation cross-validation** — CoDi vs literature marker contingency matrix (counts, percentages, heatmaps)
- **Size-specific effects** — Classifies DE genes into unique 40 nm, unique 200 nm, shared solo, shared all three, and mixture-only emergent modules
- **Reproducibility** — All parameters in `config/config.yaml`; timestamped run logs in `results/run_logs/`

---

## Requirements

- **Python** 3.10+ (recommended)
- **R** 4.x (optional, for Seurat/Azimuth annotation with `ref.Rds`)
- ~2 GB disk space after download (raw data ≈ 800 MB)

---

## Quick start

```bash
git clone https://github.com/YOUR_USERNAME/single_cell_analysis_of_immune_response_to_nanoplastic_particles.git
cd single_cell_analysis_of_immune_response_to_nanoplastic_particles

# Option A — pip (recommended on Windows)
python -m pip install -r requirements.txt

# Option B — Conda + R (for full Azimuth workflow)
conda env create -f environment.yml
conda activate nanoplastic-scRNA
```

> On Windows, the pipeline uses **Combat** integration by default (`integration_method: combat` in `config/config.yaml`) because `harmonypy` is often unavailable.

### Download data

```bash
python scripts/download_data.py
```

Downloads `.h5ad` matrices, CoDi annotation CSV, and the Azimuth `ref.Rds` reference into `data/raw/`.

### Run the main pipeline

```bash
python scripts/run_pipeline.py
```

### Typical runtime (consumer laptop)

| Step | Duration |
|------|----------|
| `pip install -r requirements.txt` | 5–15 min (once) |
| `python scripts/download_data.py` | 10–40 min |
| `python scripts/run_pipeline.py` | 30–90 min |
| `python scripts/run_additional_analyses.py` | 1–3 min |
| `python scripts/plot_supplementary_figures.py` | 1–3 min |
| `Rscript scripts/azimuth_annotation.R` (optional) | 15–45 min |

---

## Pipeline steps

| Step | Description |
|------|-------------|
| 1 — QC | Cell filtering (`min_genes`, `max_genes`, `min_counts`, `max_mt_percent`) |
| 2 — Integration | Normalization, module scores on full transcriptome, 3000 HVGs, Combat, UMAP, Leiden |
| 3 — Annotation | `ref.Rds` atlas: Seurat, CoDi, Tangram, cell2location + literature marker validation |
| 4 — Figures | UMAP, dotplots, module score maps |
| 5 — Composition | Cell-type proportions per condition |
| 6 — DE | Wilcoxon per cell type, each exposure vs control |
| 7 — Pathway | Enrichr on **UP** and **DOWN** DE genes; heatmaps, summary CSV, EN interpretation report |
| 8 — Size-specific | Gene classification by particle size (40 nm / 200 nm / mix) |
| 9 — Additional | Module scores, pseudobulk, **CoDi vs marker contingency matrix**, interpretation |
| 10 — Save | `data/processed/integrated_annotated.h5ad` |

---

## Cell-type annotation

The default workflow maps query cells to a local PBMC reference atlas (`ref.Rds`) using **Seurat/Azimuth**, then validates labels with literature markers and external CoDi annotations.

```bash
python scripts/download_data.py              # ref.Rds + idx.annoy
python scripts/run_pipeline.py               # integration
python scripts/prepare_azimuth_h5ad.py
Rscript scripts/azimuth_annotation.R         # Seurat + ref.Rds
python scripts/run_reference_annotation.py   # Tangram + cell2location (optional)
python scripts/run_pipeline.py               # full run with all labels
```

Configure in `config/config.yaml`:

```yaml
annotation:
  primary_method: seurat    # seurat | codi | tangram | cell2location | markers
  run_seurat: true
  run_tangram: false
  run_cell2location: false
```

| Column | Purpose |
|--------|---------|
| `cell_type_primary` | Drives DE and composition (default: Seurat/ref.Rds) |
| `cell_type_marker` | Literature marker scoring (validation) |
| `cell_type_codi_norm` | External CoDi labels from the dataset |

Validation dotplots: `results/figures/annotation_validation/`

### CoDi vs marker gene contingency matrix

Compares cell types assigned by **literature marker panels** (`cell_type_marker`) with **CoDi reference labels** (`cell_type_codi_norm`). Generated in pipeline step 9 (or standalone):

```bash
python scripts/run_additional_analyses.py
```

| Output | Description |
|--------|-------------|
| `results/tables/annotation_crosstab_marker_codi.csv` | Cell counts per marker × CoDi pair (contingency matrix) |
| `results/tables/annotation_crosstab_marker_codi_row_pct.csv` | Row-normalized % (dominant CoDi label per marker type) |
| `results/tables/annotation_codi_marker_mapping.csv` | Long-format mapping with counts and percentages |
| `results/figures/additional_analyses/annotation_confusion_marker_codi.png` | Heatmap (cell counts) |
| `results/figures/additional_analyses/annotation_confusion_marker_codi_normalized.png` | Heatmap (row %) |

### ref.Rds / Seurat vs CoDi contingency matrix

Compares **primary annotation** from `ref.Rds` (Seurat/Azimuth, `cell_type_ref`) with **CoDi** labels. Requires `azimuth_annotations.csv` from `Rscript scripts/azimuth_annotation.R`.

| Output | Description |
|--------|-------------|
| `results/tables/annotation_crosstab_ref_codi.csv` | Cell counts per ref.Rds × CoDi pair |
| `results/tables/annotation_crosstab_ref_codi_row_pct.csv` | Row-normalized % |
| `results/tables/annotation_codi_ref_mapping.csv` | Long-format mapping |
| `results/figures/additional_analyses/annotation_confusion_ref_codi.png` | Heatmap (cell counts) |
| `results/figures/additional_analyses/annotation_confusion_ref_codi_normalized.png` | Heatmap (row %) — **recommended for presentation** |

---

## Pathway enrichment (UP vs DOWN vs control)

Step 7 runs Enrichr separately on genes **up-regulated** and **down-regulated** in each exposure compared to control (Wilcoxon DE, padj < 0.05, \|log2FC\| > 0.25).

```bash
# Full pipeline (includes UP + DOWN enrichment)
python scripts/run_pipeline.py

# Regenerate figures + summary from saved CSV (add --enrich-missing if only UP was run)
python scripts/plot_pathway_figures.py --enrich-missing
```

| Output | Description |
|--------|-------------|
| `results/tables/pathway_enrichment_all.csv` | Full Enrichr output; **`direction`** column = `UP` or `DOWN` |
| `results/tables/pathway_enrichment_summary.csv` | Top pathways per cell type × exposure × direction with biological notes |
| `results/figures/pathway_enrichment/pathways_{cell_type}_{UP\|DOWN}.png` | Heatmap per cell type and regulation direction |
| `deliverables/Pathway_Enrichment_Interpretation_EN.md` | English report linking enriched pathways to tissue/organism-level effects |

**How to read direction:** `UP` = pathway genes are higher in the exposed condition than control; `DOWN` = lower in exposure than control.

---

## Optional: Azimuth PBMC reference (R)

```bash
python scripts/check_rscript.py
Rscript scripts/install_r_packages.R
python scripts/prepare_azimuth_h5ad.py
Rscript scripts/azimuth_annotation.R
python scripts/extract_azimuth_markers.py
```

Outputs:

- `results/tables/azimuth_annotations.csv`
- `results/tables/azimuth_marker_panels_l1.yaml`
- `results/figures/azimuth_marker_dotplot_l1.png`

---

## Additional analyses

Run after the main pipeline (requires `data/processed/integrated_annotated.h5ad`):

```bash
python scripts/run_additional_analyses.py
```

| Analysis | Method | Key output |
|----------|--------|------------|
| Cell-cycle module | `scanpy.tl.score_genes_cell_cycle` | `cell_cycle_scores_by_condition.csv` |
| IFN signature | `score_genes` on ISG panel | `ifn_scores_by_condition.csv` |
| Antigen presentation | `score_genes` on HLA/MHC panel | `antigen_presentation_scores.csv` |
| Pseudobulk | Sum UMI per condition × cell type | `pseudobulk_counts_condition_celltype.csv` |
| Annotation validation | Marker vs CoDi vs Azimuth agreement | `annotation_agreement_metrics.csv` |
| **CoDi vs marker matrix** | Contingency tables + heatmaps | `annotation_crosstab_marker_codi*.csv`, `annotation_confusion_marker_codi*.png` |
| **ref.Rds vs CoDi matrix** | Primary Seurat annotation vs CoDi | `annotation_crosstab_ref_codi*.csv`, `annotation_confusion_ref_codi*.png` |
| Figures | Bar / violin / heatmap / confusion | `results/figures/additional_analyses/` |

Module scores are computed during pipeline step 2 (before HVG subsetting). The standalone script reads them from the `.h5ad` object.

---

## Figure scripts

```bash
python scripts/plot_figures.py                 # UMAP from .h5ad
python scripts/plot_pathway_figures.py         # pathway heatmaps + UP/DOWN summary report
python scripts/plot_pathway_figures.py --enrich-missing   # re-run DOWN enrichment if missing
python scripts/plot_supplementary_figures.py   # volcano, size-specific, DE summary
python scripts/make_slides.py                  # → deliverables/nanoplastic_scRNA_results.pptx
```

---

## Project structure

```
.
├── config/
│   └── config.yaml              # All analysis parameters
├── data/
│   ├── raw/                     # Downloaded Zenodo files (gitignored)
│   └── processed/               # integrated_annotated.h5ad (gitignored)
├── scripts/
│   ├── run_pipeline.py          # Main analysis pipeline
│   ├── download_data.py         # Zenodo data fetcher
│   ├── run_additional_analyses.py
│   ├── azimuth_annotation.R     # Seurat/Azimuth mapping
│   └── ...
├── results/
│   ├── figures/                 # Generated plots (gitignored)
│   ├── tables/                  # CSV outputs (gitignored)
│   └── run_logs/                # Timestamped pipeline logs
└── deliverables/                # Reports and interpretation documents
```

Large data files and generated results are excluded via `.gitignore`. Clone the repo, install dependencies, and run the scripts to reproduce outputs locally.

---

## Key outputs

### Processed data

- `data/processed/integrated_annotated.h5ad`

### Figures

- `results/figures/` — UMAP, composition, marker dotplots
- `results/figures/pathway_enrichment/` — pathway heatmaps per cell type (**UP** and **DOWN** vs control)
- `results/figures/supplementary/` — DE volcano, size-specific plots
- `results/figures/additional_analyses/` — module scores, annotation agreement, **CoDi vs marker heatmaps**

### Tables

- `differential_expression_all.csv`
- `pathway_enrichment_all.csv` — includes `direction` (`UP` / `DOWN`)
- `pathway_enrichment_summary.csv` — top pathways with regulation direction and biological notes
- `annotation_crosstab_marker_codi.csv` — CoDi vs marker contingency matrix
- `annotation_crosstab_ref_codi.csv` — ref.Rds/Seurat vs CoDi contingency matrix
- `annotation_codi_ref_mapping.csv` — long-format ref.Rds ↔ CoDi mapping
- `size_specific_effects_summary.csv`
- `size_specific_interpretation.csv`
- `module_scores_by_condition.csv`
- `annotation_agreement_metrics.csv`

---

## Documentation

| Document | Description |
|----------|-------------|
| [`deliverables/Complete_Results_Report_EN.md`](deliverables/Complete_Results_Report_EN.md) | Full results report |
| [`deliverables/Pathway_Enrichment_Interpretation_EN.md`](deliverables/Pathway_Enrichment_Interpretation_EN.md) | UP/DOWN pathway enrichment vs control with biological interpretation (generated by pipeline) |
| [`deliverables/Size_Specific_Effects_Interpretation_EN.md`](deliverables/Size_Specific_Effects_Interpretation_EN.md) | Biological interpretation of size-specific gene modules |
| [`deliverables/Analysis_Results_Report.md`](deliverables/Analysis_Results_Report.md) | Analysis summary (Serbian) |
| [`deliverables/marker_gene_selection_report.md`](deliverables/marker_gene_selection_report.md) | Marker gene selection rationale |

---

## Limitations

- **Single donor** — no biological replicates across donors; treat statistics as within-sample descriptive effects
- **DE on 3000 HVGs** — some markers may be excluded from DE tests (module scores use the full gene space)
- **Annotation circularity** — Azimuth panels and marker scoring can partially validate the same labels

---

## Citation

If you use this code or reproduce the analysis, please cite the original dataset:

> Zenodo record **15866724** — single-cell RNA-seq of human PBMC exposed to polystyrene nanoparticles.  
> DOI: [10.5281/zenodo.15866724](https://doi.org/10.5281/zenodo.15866724)

Reference atlas (Azimuth PBMC):

> Zenodo record **4546839**.  
> DOI: [10.5281/zenodo.4546839](https://doi.org/10.5281/zenodo.4546839)

---

## License

Analysis code in this repository is provided for academic and research use. Check the Zenodo dataset license for data usage terms.
