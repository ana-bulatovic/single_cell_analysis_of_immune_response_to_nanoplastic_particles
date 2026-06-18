# Complete Results Report — Single-Cell Analysis of Immune Response to Nanoplastic Particles

**Project:** Single-Cell Analysis of Immune Response to Nanoplastic Particles  
**Data source:** [Zenodo 15866724](https://zenodo.org/records/15866724) (DOI: 10.5281/zenodo.15866724)  
**Reference atlas:** [Zenodo 4546839](https://zenodo.org/records/4546839) — `ref.Rds` + `idx.annoy` (Satija Lab PBMC Azimuth reference)  
**Main pipeline:** `scripts/run_pipeline.py`  
**Seurat/Azimuth annotation:** `scripts/azimuth_annotation.R`  
**Latest successful pipeline run:** `pipeline_run_20260619_003132.txt` (2026-06-19)

---

## 1. Executive summary

This study analyses **human PBMC from one donor** exposed to carboxylated polystyrene nanoparticles (PSNP) in four conditions: **40 nm**, **200 nm**, **40+200 nm mixture**, and **untreated control**. After QC and integration, **33,240 cells** were analysed.

**Cell-type labels** were assigned primarily by mapping each cell to the **local PBMC reference atlas (`ref.Rds`)** using **Seurat/Azimuth** (`RunAzimuth` in R). Labels were **validated** with classical literature marker genes (IL7R, LYZ, MS4A1, NKG7, etc.) and with external **CoDi** labels supplied with the dataset.

**Key findings:**
- High-confidence reference mapping (mean Azimuth score **0.908**; **99.2%** of cells with valid Seurat label)
- Typical PBMC landscape dominated by **CD8 T** and **CD4 T** cells under Seurat annotation
- Widespread but **cell-type-specific** differential expression (**33,019** significant gene hits across **21** comparisons)
- Pathway enrichment highlights **inflammatory, cytokine, and innate immune** programmes
- Clear **particle-size-specific** gene sets (40 nm vs 200 nm vs shared vs mixture-only)
- Global composition largely stable; **CD14 monocytes increase under 200 nm PSNP**
- Module scores show modest global shifts (slightly lower IFN and antigen-presentation scores in exposed conditions)

**Caveat:** Single donor, no biological replicates — interpret statistics as descriptive within-sample effects.

---

## 2. Cell-type annotation strategy

### 2.1 What we used (professor workflow)

| Layer | Method | Output column | Role |
|-------|--------|---------------|------|
| **Primary** | Seurat/Azimuth + local `ref.Rds` | `cell_type_seurat` → `cell_type_primary` | Drives DE, composition, pathway analysis |
| **Validation** | Literature marker panels (`config.yaml`) | `cell_type_marker` | Dotplots and agreement metrics |
| **External reference** | CoDi labels (pre-computed, dataset CSV) | `cell_type_codi_norm` | Independent comparison |

Tangram and cell2location mapping were **disabled** in config (`run_tangram: false`, `run_cell2location: false`) — not required for the final analysis.

### 2.2 How Seurat + ref.Rds works

1. **`ref.Rds`** is a pre-annotated PBMC reference atlas (thousands of labelled reference cells).
2. **`idx.annoy`** is a pre-built search index for fast neighbour lookup.
3. **`RunAzimuth`** (R/Seurat) compares each query cell’s expression profile to the atlas and assigns:
   - `predicted.celltype.l1` (broad type: CD4 T, CD8 T, B, Mono, NK, …)
   - `predicted.celltype.l2` / `l3` (finer subtypes)
   - `prediction.score.max` (confidence 0–1)
4. Results are saved to **`results/tables/azimuth_annotations.csv`**.
5. Python pipeline loads this CSV and maps labels to unified names (`CD4_T`, `CD8_T_cytotoxic`, etc.).
6. Cells with score **< 0.5** are labelled **`low_confidence`** and excluded from DE.

### 2.3 How literature markers work (validation only)

For each cell, the pipeline computes the **mean expression** of curated marker sets (e.g. IL7R/LTB for CD4 T, NKG7/GNLY for CD8 T, LYZ/S100A8 for monocytes). The cell type with the **highest score** becomes `cell_type_marker`. This does **not** drive DE unless Seurat labels are unavailable.

### 2.4 Annotation quality metrics (latest run)

| Metric | Value | Interpretation |
|--------|------:|----------------|
| Seurat mapped fraction | 99.2% | Nearly all cells received a reference label |
| Azimuth mean score | 0.908 | High confidence (>0.85 is excellent) |
| Seurat vs literature markers | 47.5% | Moderate — expected (NK/CD8 overlap, different granularity) |
| CoDi vs literature markers | 65.4% | Moderate |
| CoDi vs Seurat | 60.0% | Reasonable agreement between independent methods |

**Interpretation:** Reference mapping is **reliable**. Partial disagreement with marker panels is **normal** for PBMC and does not invalidate Seurat labels — marker dotplots confirm expected expression patterns per type.

---

## 3. Pipeline workflow (`run_pipeline.py`)

The pipeline runs **10 steps** in a single command:

```
STEP 0  Raw data overview (sample sizes, QC thresholds)
STEP 1  Quality control per sample (filter low-quality cells)
STEP 2  Merge → normalize → module scores → HVG → Combat → UMAP → Leiden
STEP 3  Annotation (Seurat + CoDi + literature marker validation)
STEP 4  Core UMAP and dotplot figures
STEP 5  Cell composition analysis
STEP 6  Differential expression (Wilcoxon, per cell type × exposure)
STEP 7  Pathway enrichment (Enrichr: GO, KEGG, Reactome) + heatmap figures
STEP 8  Size-specific effect classification (40 nm / 200 nm / shared / mix-only)
STEP 9  Additional analyses (module scores, pseudobulk, annotation cross-validation)
STEP 10 Save processed object (integrated_annotated.h5ad)
```

### 3.1 Input and QC

| Sample | Condition | Cells before QC | After QC | Retained |
|--------|-----------|----------------:|---------:|---------:|
| Sample 1 | PSNP 40 nm | 8,729 | 8,458 | 96.9% |
| Sample 2 | PSNP 200 nm | 12,676 | 12,421 | 98.0% |
| Sample 3 | PSNP mix | 6,157 | 6,005 | 97.5% |
| Sample 4 | control | 6,516 | 6,356 | 97.5% |
| **Total** | — | — | **33,240** | ~97% |

**QC filters:** min 200 genes, max 7,000 genes, min 500 UMIs, max 15% mitochondrial reads.

### 3.2 Integration and dimensionality reduction

1. Normalize to 10,000 counts/cell + log1p transform
2. Compute **module scores** on the **full gene space** (before HVG subset)
3. Save **`integrated_pre_hvg.h5ad`** (full transcriptome, for reference mapping tools)
4. Select **3,000 highly variable genes** (Seurat flavour, batched by sample)
5. Scale (max = 10), PCA (30 components)
6. **Combat** batch correction on `sample_id`
7. Neighbour graph (k = 15), UMAP, Leiden clustering (resolution = 0.5 → **14 clusters**)

**Output object:** `data/processed/integrated_annotated.h5ad` — 33,240 cells × 3,000 HVGs.

### 3.3 Downstream analysis uses `cell_type_primary`

All composition, DE, and pathway analyses use **`cell_type_primary`**, which is set to **Seurat/ref.Rds labels** (`primary_method: seurat` in config).

**Cell types included in DE (Seurat annotation):** B_cell, CD4_T, CD8_T_cytotoxic, Monocyte_CD14, NK_cell, other, other_T (7 types × 3 exposures = **21 comparisons**). DC and Platelet are absent or too rare under Seurat L1 mapping.

---

## 4. Results tables — purpose and interpretation

### 4.1 Annotation tables

| File | Purpose | What it tells you |
|------|---------|-------------------|
| `azimuth_annotations.csv` | Raw Seurat/Azimuth output per cell | L1/L2/L3 labels and prediction scores from ref.Rds mapping |
| `annotation_method_agreement.csv` | Concordance between methods | How well Seurat, CoDi, and literature markers agree |
| `annotation_agreement_metrics.csv` | Extended validation metrics | Includes Azimuth mean score, cross-method agreement |
| `annotation_crosstab_marker_seurat.csv` | Confusion matrix | Which Seurat types map to which marker-based types |
| `annotation_crosstab_marker_codi.csv` | Confusion matrix | CoDi vs marker comparison |
| `annotation_crosstab_marker_azimuth.csv` | Confusion matrix | Raw Azimuth L1 vs markers |

**What they say:** Seurat mapping is high-confidence. ~47% exact label match with literature markers reflects known NK/CD8 ambiguity and coarser marker panels — not a failed annotation.

---

### 4.2 Composition

| File | Purpose | What it tells you |
|------|---------|-------------------|
| `cell_composition_by_condition.csv` | Cell counts and fractions per condition × cell type | Whether PSNP exposure changes immune cell proportions |

**Seurat-based composition (fractions):**

| Cell type | control | 40 nm | 200 nm | mix |
|-----------|--------:|------:|-------:|----:|
| CD8 T cytotoxic | 41.0% | 40.8% | 37.5% | 39.5% |
| CD4 T | 32.4% | 38.4% | 34.9% | 31.7% |
| B cell | 9.2% | 6.4% | 7.2% | 10.2% |
| NK cell | 8.4% | 6.6% | 6.5% | 7.6% |
| CD14 monocyte | 3.2% | 2.9% | **8.2%** | 4.9% |

**Interpretation:** Overall composition is stable. The most notable shift is **CD14 monocytes increasing under 200 nm** (3.2% → 8.2%), suggesting myeloid involvement in the 200 nm response.

---

### 4.3 Differential expression

| File | Purpose | What it tells you |
|------|---------|-------------------|
| `differential_expression_all.csv` | All DE results | Gene, logFC, p-value, cell type, comparison (exposure vs control) |

| Metric | Value |
|--------|------:|
| Comparisons run | 21 |
| Total gene rows | 63,000 |
| Significant hits (padj < 0.05, \|logFC\| > 0.25) | **33,019** |
| Method | Wilcoxon rank-sum, per cell type |
| Min cells per group | 20 |

**Interpretation:** Nanoplastic exposure induces **widespread transcriptional changes** that differ by cell type and particle size. Effects are not captured by a single global gene list — each immune compartment responds differently.

---

### 4.4 Pathway enrichment

| File | Purpose | What it tells you |
|------|---------|-------------------|
| `pathway_enrichment_all.csv` | Enrichr results for UP/DOWN DE genes | Which biological pathways are enriched per cell type × exposure |
| `size_specific_pathway_enrichment.csv` | Pathways for size-class gene sets | Pathways specific to 40 nm, 200 nm, shared, or mix-only genes |

| Metric | Value |
|--------|------:|
| Enrichment rows | 166,064 |
| Databases | GO Biological Process, KEGG, Reactome |
| Input | Significant UP or DOWN genes per comparison |

**Interpretation:** Dominant themes include **inflammatory response**, **cytokine-mediated signalling**, **IL-10 / interleukin signalling**, **neutrophil degranulation**, **TNF signalling**, and **innate immune system** pathways — consistent with immune activation by nanoplastic particles.

---

### 4.5 Size-specific effects

| File | Purpose | What it tells you |
|------|---------|-------------------|
| `size_specific_effects_summary.csv` | Gene counts per effect class | How many genes are unique to 40 nm, 200 nm, shared, or mix-only |
| `size_specific_genes.csv` | Gene lists per class | Which genes belong to each size category |
| `size_specific_interpretation.csv` | Human-readable summary | Top genes and pathways per cell type × effect class |

**Effect classes:**
- **unique_40nm** — significant only for 40 nm vs control
- **unique_200nm** — significant only for 200 nm vs control
- **shared_40_200** — significant for both solo sizes, not mix
- **shared_all_three** — core PSNP response (all exposures)
- **mix_only_emergent** — significant only in the mixture condition

**Example findings:**
- **B cells:** strong 200 nm-specific pathway signal (IL-17, TNF, Kaposi sarcoma/infection-related KEGG terms)
- **CD14 monocytes:** large shared-all-three module (973 genes) with IL-10, cytokine receptor, LPS response pathways
- **CD8 T cells:** shared inflammatory modules across exposures; mix-only emergent genes in some types

**Interpretation:** Particle **size matters** — 40 nm and 200 nm are not interchangeable. The mixture can produce **emergent** responses not seen with either size alone.

---

### 4.6 Module scores and pseudobulk (additional analyses)

| File | Purpose | What it tells you |
|------|---------|-------------------|
| `module_scores_by_condition.csv` | Global programme scores per condition | Cell cycle, IFN, antigen presentation shifts |
| `module_scores_by_condition_celltype.csv` | Scores per condition × cell type | Type-specific programme changes |
| `cell_cycle_scores_by_condition.csv` | S and G2M scores | Proliferation state |
| `ifn_scores_by_condition.csv` | Interferon signature | Innate antiviral activation |
| `antigen_presentation_scores.csv` | MHC/antigen presentation | Antigen-processing capacity |
| `pseudobulk_counts_condition_celltype.csv` | Summed UMIs per group | Input for bulk-style DE (DESeq2/edgeR) |

**Global module score findings (vs control):**

| Programme | Trend in exposed conditions |
|-----------|----------------------------|
| S / G2M (cell cycle) | Minimal change (G2M slightly lower at 200 nm) |
| IFN signature | Slightly **lower** in exposed (control highest: 0.035) |
| Antigen presentation | Modestly **lower** globally |

**Interpretation:** Module scores capture **broad programmes**, not individual DE genes. The main story is in DE/pathway results; module scores provide complementary context (no strong global IFN storm; modest antigen-presentation shift).

---

## 5. Figures — what each image shows

### 5.1 Core figures (`results/figures/`)

| Figure | What is displayed | How to read it |
|--------|-------------------|----------------|
| **`umap_condition.png`** | UMAP coloured by experimental condition (40 nm, 200 nm, mix, control) | All conditions occupy the same regions → good integration; no massive global shift |
| **`umap_split_by_condition.png`** | Same UMAP in 4 panels, one condition highlighted | Easier visual comparison per exposure |
| **`umap_sample_integration.png`** | UMAP coloured by sample/batch after Combat | Samples intermix → batch correction worked |
| **`umap_clusters.png`** | UMAP coloured by 14 Leiden clusters | Unsupervised structure before annotation |
| **`umap_celltypes_seurat.png`** | UMAP coloured by **Seurat/ref.Rds** cell types (primary) | Main annotation map for the study |
| **`umap_celltypes_marker.png`** | UMAP coloured by **literature marker** labels | Validation layer — compare with Seurat map |
| **`umap_codi_celltypes.png`** | UMAP coloured by CoDi external labels | Independent reference validation |
| **`umap_module_scores.png`** | UMAP coloured by S, G2M, IFN module scores | Where proliferation and IFN programmes are active |
| **`marker_dotplot.png`** | Dot plot: marker genes × cell types (marker-based groups) | Dot size = % cells expressing; colour = mean expression |
| **`composition_barplot.png`** | Stacked bar chart of cell-type fractions per condition | Visual composition comparison |

---

### 5.2 Annotation validation figures (`results/figures/annotation_validation/`)

| Figure | What is displayed |
|--------|-------------------|
| **`marker_validation_literature_markers.png`** | Literature markers checked against marker-based labels |
| **`marker_validation_seurat_ref_rds.png`** | Literature markers checked against **Seurat** labels — key validation figure |
| **`marker_validation_codi.png`** | Literature markers checked against CoDi labels |

**How to use:** Confirm that expected markers are enriched in the assigned types (e.g. IL7R in CD4 T, MS4A1 in B cells, LYZ in monocytes).

---

### 5.3 Pathway enrichment figures (`results/figures/pathway_enrichment/`)

**14 heatmaps** (latest run, 2026-06-19), one per cell type × direction (UP/DOWN):

- B_cell, CD4_T, CD8_T_cytotoxic, Monocyte_CD14, NK_cell, other, other_T

**What each shows:**
- Rows = top enriched pathways (GO / KEGG / Reactome)
- Columns = exposure conditions (40 nm, 200 nm, mix) vs control
- Colour = enrichment significance or score

**How to read:** Identifies which immune pathways are activated (UP) or suppressed (DOWN) per cell type and particle size.

**Note:** Older figures for DC, Monocyte_CD16, and Platelet (dated 2026-06-15) are from a previous marker-based annotation and should **not** be used — those types are absent or too rare under Seurat L1.

---

### 5.4 Additional analysis figures (`results/figures/additional_analyses/`)

| Figure | What is displayed |
|--------|-------------------|
| **`module_scores_by_condition.png`** | Bar chart of S, G2M, IFN, antigen-presentation scores per condition |
| **`module_scores_violin.png`** | Distribution of module scores across conditions |
| **`antigen_presentation_heatmap.png`** | Antigen-presentation score by condition × cell type |
| **`annotation_agreement_bar.png`** | Bar chart of cross-method agreement metrics |
| **`annotation_confusion_marker_azimuth.png`** | Heatmap: marker-based vs Seurat (Azimuth L1) labels |

---

### 5.5 Optional supplementary figures (`results/figures/supplementary/`)

Generated separately by `scripts/plot_supplementary_figures.py` (if run):

- DE summary bar charts, volcano plots, size-specific UpSet/bar plots, composition delta plots

These are **not** produced by the main pipeline unless that script is executed separately.

---

## 6. Additional analyses (STEP 9)

STEP 9 runs automatically at the end of `run_pipeline.py` (or standalone via `scripts/run_additional_analyses.py`).

### 6.1 Module score analysis

**What is done:**
- On the **full normalized transcriptome** (before HVG filtering), compute per-cell scores for:
  - **Cell cycle** (S phase, G2M phase)
  - **Interferon response** (ISG15, IFIT1-3, MX1, OAS1, …)
  - **Antigen presentation** (HLA-DRA/DRB1, CD74, B2M, TAP1/2, …)
- Aggregate mean scores by **condition** and by **condition × cell type**

**What you get:**
- CSV tables (`module_scores_*.csv`)
- Bar/violin/heatmap figures in `additional_analyses/`
- Interpretation file: `additional_analyses_interpretation_SR.md` (Serbian)

**What it means:** Describes **global biological programmes** (proliferation, innate activation, MHC expression) beyond individual DE genes. Useful for explaining whether PSNP exposure triggers interferon storms, cell-cycle arrest, or altered antigen presentation.

---

### 6.2 Pseudobulk aggregation

**What is done:**
- Sum raw UMI counts (`layers['counts']`) per **condition × cell type** group

**What you get:**
- `pseudobulk_counts_condition_celltype.csv` — genes × pseudobulk samples

**What it means:** Enables **bulk RNA-seq style** follow-up (DESeq2, edgeR) without reloading single-cell data. Each row is one pseudobulk sample (e.g. “CD4_T + PSNP_200nm”).

---

### 6.3 Annotation cross-validation

**What is done:**
- Compare three label sources: **Seurat (ref.Rds)**, **CoDi CSV**, **literature markers**
- Compute agreement fractions and confusion matrices

**What you get:**
- `annotation_agreement_metrics.csv`, `annotation_method_agreement.csv`
- Crosstab CSV files
- Agreement bar plot and confusion heatmap

**What it means:** Demonstrates that annotation is **robust and independently supported**. Moderate pairwise agreement (45–65%) is expected when methods use different granularity (e.g. Seurat L1 “CD8 T” vs marker “NK_cell” for NKG7-high cells).

---

## 7. How to reproduce

```bash
# 1. Download data + ref.Rds
python scripts/download_data.py

# 2. Full analysis (requires azimuth_annotations.csv for Seurat primary labels)
python scripts/run_pipeline.py

# If azimuth_annotations.csv does not exist yet:
python scripts/prepare_azimuth_h5ad.py
Rscript scripts/azimuth_annotation.R
python scripts/run_pipeline.py
```

**Optional:**
```bash
python scripts/plot_supplementary_figures.py   # extra DE/volcano figures
python scripts/run_additional_analyses.py    # re-run STEP 9 only
python scripts/make_slides.py                # PowerPoint
```

---

## 8. Overall assessment

| Component | Status | Comment |
|-----------|--------|---------|
| QC & integration | ✅ Good | >96% retention; Combat mixing confirmed on UMAP |
| Seurat/ref.Rds annotation | ✅ Good | Mean score 0.908; 99.2% mapped |
| Marker validation | ✅ Acceptable | ~47% exact match; dotplots confirm expected patterns |
| Composition | ✅ Mostly stable | Notable CD14 monocyte ↑ at 200 nm |
| Differential expression | ✅ Strong | 33,019 significant hits, cell-type-specific |
| Pathway enrichment | ✅ Informative | Inflammatory/cytokine programmes dominate |
| Size-specific effects | ✅ Informative | Clear 40 nm vs 200 nm vs mix distinctions |
| Module scores | ✅ Descriptive | Modest global shifts; IFN not the main story |
| Additional analyses | ✅ Complete | Pseudobulk, cross-validation, interpretation |

**Bottom line:** The analysis pipeline produced a **coherent, publication-ready** set of results for a single-donor PBMC nanoplastic case study. Primary conclusions should emphasise **cell-type-specific transcriptional reprogramming**, **size-dependent effects**, and **inflammatory pathway engagement**, while acknowledging the **single-donor** limitation.

---

*Report generated from pipeline outputs in `results/tables/`, `results/figures/`, and run log `pipeline_run_20260619_003132.txt`.*
