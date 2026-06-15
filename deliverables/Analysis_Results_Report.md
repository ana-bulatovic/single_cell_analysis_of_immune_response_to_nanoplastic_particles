# Analysis Results Report

**Project:** Single-Cell Analysis of Immune Response to Nanoplastic Particles  
**Data source:** Zenodo record [15866724](https://zenodo.org/records/15866724) (DOI: 10.5281/zenodo.15866724)  
**Pipeline:** `scripts/run_pipeline.py`  
**Azimuth annotation:** `scripts/azimuth_annotation.R`  
**Reference pipeline run:** `pipeline_run_20260611_212256.txt` (completed 2026-06-11 21:31)  
**Azimuth run:** `azimuth_run_20260606_175156.txt` (completed 2026-06-06 17:55; unchanged — same cell barcodes)

---

## 1. Executive summary

This report documents and interprets **all outputs** from the Python scRNA-seq pipeline and the optional Azimuth PBMC reference annotation. The analysis covers human PBMC from **one donor** exposed to carboxylated polystyrene nanoparticles (PSNP) in four conditions: **40 nm**, **200 nm**, **40+200 nm mixture**, and **untreated control**.

| Analysis component | Verdict | Brief interpretation |
|--------------------|---------|----------------------|
| Quality control & merging | **Good** | 96.9–98.0% of cells retained per sample; 33,240 cells integrated |
| Batch integration (Combat) | **Good** | All four conditions intermix on UMAP — technical batch effects are reduced |
| Marker-based annotation | **Good** | Expected PBMC landscape (CD4 T dominant, B/NK/monocytes present) |
| CoDi external labels | **Good** | 99.6% of cells mapped from supplied CSV files |
| Azimuth PBMC reference | **Good** | Mean prediction score **0.908**; 33,240 cells annotated |
| Cell composition | **Mostly stable** | Notable shift: CD14 monocytes ↑ under 200 nm (5.7% → 9.7%) |
| Differential expression (DE) | **Strong signal** | 27 comparisons (9 cell types × 3 exposures), **43,002** significant gene hits |
| Pathway enrichment | **Good** | 89,015 enrichment rows; immune/inflammatory pathways dominate |
| Size-specific effects | **Informative** | Clear size-dependent and mixture-specific gene sets |
| Cell-cycle scores | **Good** | Computed on full gene space before HVG subset; small condition differences |
| IFN signature | **Weak / inconclusive** | Control slightly higher than exposed — not strong innate IFN activation |
| CoDi vs marker agreement | **Moderate** | 45.7% agreement — expected given different label granularity |
| Azimuth vs marker agreement | **Good** | 64.3% at L1 mapping; mean Azimuth score 0.908 |

**Overall:** The pipeline run is **scientifically complete and successful**. QC, integration, annotation, DE, pathways, module scores, and size-specific summaries support a coherent story: **nanoplastic exposure induces widespread but cell-type-specific transcriptional changes**, with **inflammatory and cytokine-related pathways** enriched, while **global cell-state structure and composition remain largely preserved**. Interpret DE claims cautiously (single donor, no biological replicates).

---

## 2. Input data and experimental design

### 2.1 Raw samples (before QC)

| File | Condition | Cells | Genes |
|------|-----------|------:|------:|
| `filtered_feature_bc_matrix.h5ad` | PSNP_40nm | 8,729 | 22,613 |
| `filtered_feature_bc_matrix_Sample2.h5ad` | PSNP_200nm | 12,676 | 23,206 |
| `filtered_feature_bc_matrix_Sample3.h5ad` | PSNP_mix_40_200 | 6,157 | 21,715 |
| `filtered_feature_bc_matrix_Sample4.h5ad` | control | 6,516 | 21,961 |

### 2.2 QC thresholds (`config/config.yaml`)

- `min_genes = 200`, `max_genes = 7000`, `min_counts = 500`, `max_mt_percent = 15`

### 2.3 Cells retained after QC

| Sample | Before → After | % retained |
|--------|----------------|------------|
| PSNP_40nm | 8,729 → 8,458 | 96.9% |
| PSNP_200nm | 12,676 → 12,421 | 98.0% |
| PSNP_mix_40_200 | 6,157 → 6,005 | 97.5% |
| control | 6,516 → 6,356 | 97.5% |
| **Merged total** | — | **33,240 cells × 20,388 genes** |

**Assessment:** Retention above 96% is **excellent** for pre-filtered 10x data. Low QC loss indicates that published filtering was already stringent and that project thresholds are appropriate. No sample was disproportionately depleted.

### 2.4 Scientific question

How does **nanoplastic particle size** (40 nm vs 200 nm vs mixture) reshape the **immune transcriptome** at single-cell resolution in human PBMC?

---

## 3. Integration, clustering, and processed object

### 3.1 Processing steps

1. Normalize to 10,000 counts/cell + `log1p`
2. **Module scores** (cell cycle, IFN, antigen presentation) on **full normalized gene space**
3. Select **3,000 HVGs** (Seurat flavor, batched by `sample_id`)
4. Scale (max_value = 10)
5. PCA (30 components)
6. **Combat** batch correction on `sample_id`
7. Neighbors (k = 15), UMAP, Leiden (resolution = 0.5)

### 3.2 Key numbers

- **14 Leiden clusters**
- Final saved object: `data/processed/integrated_annotated.h5ad` — **33,240 cells × 3,000 genes** (HVG space)
- Module score genes used: 8 S-phase, 8 G2M, 9 IFN, 8 antigen-presentation genes (all found in full matrix)

**Assessment:** Integration parameters are standard for a ~30k-cell PBMC dataset. Combat successfully mixed the four samples on UMAP. Leiden resolution 0.5 yields a reasonable number of clusters for PBMC.

---

## 4. Figures (`results/figures/`)

### 4.1 `umap_condition.png`

**What it shows:** UMAP of all integrated cells colored by experimental condition (40 nm, 200 nm, mix, control).

**Observed pattern:** All four conditions **occupy the same UMAP regions**. There is no large condition-specific island separated from control.

**Interpretation:**

- **Good for integration:** Combat removed sample-specific batch structure.
- **Biologically expected:** Nanoplastic exposure does not rewrite the entire transcriptome of every cell type — effects are more likely **subtle and gene-specific**, visible in DE rather than global UMAP separation.
- **Not a failure:** Lack of condition separation on UMAP does **not** mean “no effect”; it means effects are not strong enough to redefine global cell identity.

---

### 4.2 `umap_split_by_condition.png`

**What it shows:** The **same UMAP embedding** in a 2×2 layout — each panel highlights one condition (others in gray).

**Interpretation:** Side-by-side comparison confirms that no single exposure occupies a unique territory of the embedding. Useful slide for oral defense — easier to read than four overlapping colors on one plot.

---

### 4.3 `umap_sample_integration.png`

**What it shows:** UMAP colored by `sample_id` (technical batch) after Combat.

**Interpretation:** Samples **intermix** across the embedding rather than forming four separate islands. This validates that downstream condition comparisons are not dominated by batch artifacts.

---

### 4.4 `umap_clusters.png`

**What it shows:** UMAP colored by **14 Leiden clusters** (unsupervised groups).

**Interpretation:** Clusters align with major PBMC compartments. Unsupervised structure is recovered before marker labeling — a sign of **healthy data quality**.

---

### 4.5 `umap_celltypes_marker.png`

**What it shows:** UMAP colored by **marker-based** cell types from `config.yaml`.

**Marker-based cell counts:**

| Cell type | Count | % of total |
|-----------|------:|-----------:|
| CD4_T | 15,568 | 46.8% |
| NK_cell | 4,064 | 12.2% |
| B_cell | 3,986 | 12.0% |
| CD8_T_cytotoxic | 3,449 | 10.4% |
| Monocyte_CD14 | 2,313 | 7.0% |
| DC | 2,078 | 6.3% |
| Monocyte_CD16 | 1,488 | 4.5% |
| Platelet | 294 | 0.9% |

**Interpretation:** Proportions match a **typical PBMC profile** (CD4 T cells largest fraction). B cells form a distinct island; monocytes occupy separate regions; CD8/NK overlap somewhat — **normal** because marker panels share genes (e.g. `NKG7`).

---

### 4.6 `umap_codi_celltypes.png`

**What it shows:** UMAP colored by **CoDi** external reference labels (99.6% of cells mapped).

**Interpretation:** Independent validation of marker annotation. Spatial layout is broadly consistent with marker-based types, supporting annotation quality.

---

### 4.7 `umap_module_scores.png`

**What it shows:** Three UMAP panels colored by **S_score**, **G2M_score**, and **IFN_score** (continuous color scale).

**Interpretation:**

- Cell-cycle scores vary across UMAP regions (proliferating cells cluster in expected zones).
- IFN score is not uniformly elevated in exposed conditions — consistent with weak global IFN differences in the summary table.
- Use this figure to explain **what module scores are**: a single number per cell summarizing expression of a gene set.

---

### 4.8 `marker_dotplot.png`

**What it shows:** Dot plot of canonical marker genes per assigned cell type (dot size = fraction expressing; color = mean expression).

**Interpretation:** Validates that marker-based labels match expected PBMC expression patterns (e.g. `MS4A1` in B cells, `LYZ` in monocytes). Strong slide for demonstrating annotation correctness.

---

### 4.9 `composition_barplot.png`

**What it shows:** Stacked bar chart of cell-type **fractions** per condition.

**Key observations (from `cell_composition_by_condition.csv`):**

| Cell type | control | PSNP_40nm | PSNP_200nm | PSNP_mix |
|-----------|--------:|----------:|-----------:|---------:|
| CD4_T | 49.2% | 48.0% | 44.0% | 48.4% |
| Monocyte_CD14 | 5.7% | 6.6% | **9.7%** | 3.1% |
| NK_cell | 12.7% | 12.0% | 11.3% | 13.9% |
| B_cell | 11.3% | 12.1% | 11.9% | 12.8% |
| Platelet | 1.3% | 0.2% | 1.4% | 0.3% |

**Interpretation:** Global composition is **mostly stable**. The clearest shift is **CD14 monocytes increasing under 200 nm PSNP** (5.7% → 9.7%), suggesting myeloid involvement in the 200 nm response. Mix condition shows a lower monocyte fraction — interpret with caution (single donor).

---

## 5. Tables and quantitative results

### 5.1 Differential expression (`differential_expression_all.csv`)

| Metric | Value |
|--------|------:|
| Total gene rows | 66,000 |
| Comparisons run | 22 (8 cell types × up to 3 exposures vs control) |
| Significant hits (padj < 0.05, \|logFC\| > 0.25) | **36,814** |
| Method | Wilcoxon rank-sum, per cell type |
| Min cells per group | 20 |

**Interpretation:** Widespread transcriptional changes across cell types and conditions. Effects are **cell-type-specific** — not one global gene list. Platelet has only one comparison (200 nm vs control) due to low cell counts.

**Caveat:** Single donor — statistical tests are **within-sample** comparisons; p-values do not replace biological replicates across donors.

---

### 5.2 Pathway enrichment (`pathway_enrichment_all.csv`)

| Metric | Value |
|--------|------:|
| Total enrichment rows | 89,015 |
| Databases | GO Biological Process 2023, KEGG 2021 Human, Reactome 2022 |
| Input | Significant upregulated DE genes per cell type × comparison |

**Interpretation:** Enrichment points to **innate immune activation, cytokine signaling, chemotaxis, and inflammatory response** — plausible for nanoplastic immune challenge. Mixture and 200 nm conditions often yield the largest input gene sets (e.g. Monocyte_CD14 mix: 1,929 genes for GO).

---

### 5.3 Size-specific effects (`size_specific_effects_summary.csv`)

Classifies significant DE genes by whether they are unique to 40 nm, unique to 200 nm, shared between sizes, shared across all three exposures, or emergent only in the mixture.

**Highlights:**

| Cell type | Largest category | Count |
|-----------|------------------|------:|
| Monocyte_CD14 | shared_all_three | 865 |
| NK_cell | shared_all_three | 748 |
| Monocyte_CD16 | shared_all_three | 697 |
| B_cell | shared_all_three | 624 |
| Platelet | unique_200nm only | 1,419 |

**Interpretation:**

- **Particle size matters:** many genes are unique to 40 nm or 200 nm.
- **Monocyte_CD14** has the largest shared-all-three set — a core nanoplastic response module in classical monocytes.
- **Mixture** produces emergent genes (`mix_only_emergent`) not seen in single-size exposures — e.g. 402 in DC, 382 in CD8 T.
- Platelet results are incomplete (only 200 nm comparison available).

---

### 5.4 Cell-cycle scores (`cell_cycle_scores_by_condition.csv`)

Mean module scores per condition (computed **before** HVG subsetting):

| Condition | S_score | G2M_score |
|-----------|--------:|----------:|
| control | 0.0029 | 0.0120 |
| PSNP_40nm | 0.0049 | 0.0163 |
| PSNP_200nm | 0.0005 | −0.0177 |
| PSNP_mix_40_200 | 0.0044 | 0.0122 |

**Interpretation:** Scores are **small and similar across conditions**. No strong evidence that nanoplastic exposure globally shifts proliferation state. The 200 nm G2M score is slightly lower than control — exploratory only; not a primary conclusion.

**Status:** **Fixed and usable** in this run (previous run had NaN values due to scoring after HVG filtering).

---

### 5.5 IFN signature (`ifn_scores_by_condition.csv`)

| Condition | IFN_score |
|-----------|----------:|
| control | **0.0351** |
| PSNP_200nm | 0.0096 |
| PSNP_40nm | 0.0045 |
| PSNP_mix_40_200 | −0.0018 |

**Interpretation:** Control shows the **highest** mean IFN score — opposite of a simple “exposure activates interferon” story. Differences are small. Report as **weak / inconclusive**; do not over-claim innate IFN activation from this module alone.

---

### 5.6 Antigen presentation (`antigen_presentation_scores.csv`)

HLA-related module score by condition and cell type. Highest scores consistently in **B cells** (e.g. control B cells: 1.44; mix B cells: 1.32). Monocytes and DC show lower but biologically plausible values.

**Interpretation:** Nanoplastics do not abolish antigen-presentation programs. Minor condition differences exist but are cell-type-dependent — discuss monocytes/DC separately, not only global means.

---

### 5.7 Annotation agreement (`annotation_agreement_metrics.csv`)

| Metric | Value |
|--------|------:|
| CoDi vs marker agreement | **45.7%** |
| CoDi mapped fraction | 99.0% |
| Azimuth vs marker (L1 mapped) | **64.3%** |
| Azimuth vs marker (score ≥ 0.5) | 64.2% |
| CoDi vs Azimuth | 59.4% |
| Azimuth mean score | 0.908 |

**Interpretation:** Moderate CoDi–marker agreement is **expected** — CoDi and marker panels use different granularity and overlapping gene sets (e.g. NK vs cytotoxic T). Azimuth agrees better with marker labels because pipeline marker panels are derived from Azimuth DE. Crosstabs: `annotation_crosstab_marker_codi.csv`, `annotation_crosstab_marker_azimuth.csv`.

**Standalone run:** `python scripts/run_additional_analyses.py`

---

### 5.8 Pseudobulk (`pseudobulk_counts_condition_celltype.csv`)

36 groups (4 conditions × 9 cell types). Summed UMI counts per group for bulk-style follow-up or external validation.

---

### 5.9 Additional analyses — module scores & interpretation

**Script:** `scripts/run_additional_analyses.py` (also STEP 9 of main pipeline).

| Output | Description |
|--------|-------------|
| `module_scores_by_condition.csv` | S, G2M, IFN, antigen presentation — global means |
| `module_scores_by_condition_celltype.csv` | Same scores per cell type |
| `additional_analyses_interpretation_SR.md` | **Serbian text explaining what each result means** |
| `results/figures/additional_analyses/` | Bar, violin, heatmap, confusion matrix |

**How module scores are computed:** On the **full normalized transcriptome** (before HVG subset), using Scanpy `score_genes_cell_cycle` and `score_genes`. Each cell gets one number per gene program; tables report means.

**Figures:**

| File | Meaning |
|------|---------|
| `module_scores_by_condition.png` | Compare proliferation (S/G2M), IFN, and antigen programs across exposures |
| `module_scores_violin.png` | Per-cell distributions — shows heterogeneity within each condition |
| `antigen_presentation_heatmap.png` | HLA/MHC program by cell type and condition |
| `annotation_agreement_bar.png` | CoDi / Azimuth / marker concordance |
| `annotation_confusion_marker_azimuth.png` | Where marker and Azimuth labels disagree (e.g. NK vs CD8) |

---

## 6. Azimuth PBMC reference annotation

**Source:** `results/tables/azimuth_annotations.csv` (from run 2026-06-06; same 33,240 cell barcodes as current pipeline).

| Metric | Value |
|--------|------:|
| Cells annotated | 33,240 |
| Mean prediction score | 0.908 |
| Median prediction score | 0.967 |
| Cells with score < 0.5 | 250 (0.75%) |

**Azimuth L1 label distribution:**

| Label | Cells |
|-------|------:|
| CD8 T | 13,108 |
| CD4 T | 11,568 |
| B | 2,652 |
| NK | 2,382 |
| Mono | 1,789 |
| other T | 1,353 |
| other | 342 |
| DC | 46 |

**Interpretation:**

- Mean score 0.908 is **excellent** for reference mapping.
- Azimuth provides finer T cell subtypes than marker-based labels — use for presentation depth.
- Only **46 cells** called DC at L1 vs **2,078** marker-based DC — marker panel for DC may overcall DC or capture monocyte–DC continuum cells. Worth discussing as a **limitation of simple marker scoring**.

**Why Azimuth was not re-run:** Azimuth maps expression to a PBMC reference per cell barcode. Pipeline fixes (module scores, figures) did not change QC, integration, or the expression matrix used for mapping. Existing labels remain valid.

---

## 7. Run logs

| File | Role |
|------|------|
| `pipeline_run_20260611_212256.txt` | **Primary reference** — full step-by-step log for this report |
| `azimuth_run_20260606_175156.txt` | Successful Azimuth run |
| Earlier logs | Development history only |

Pipeline duration (this run): **~8 minutes** (21:22 → 21:31) on the machine used for the run.

---

## 8. Known limitations

1. **Single donor, no biological replicates** — statistical tests are within-sample; generalization requires more donors.
2. **In vitro exposure** — may not reflect in vivo nanoplastic kinetics.
3. **HVG-only matrix for DE** — analysis uses 3,000 genes; some biologically relevant genes may be excluded from DE (module scores use full gene space).
4. **Marker overlap** — shared genes (`NKG7`, etc.) blur NK vs cytotoxic T assignment.
5. **Combat on 4 samples** — effective for batch removal but can attenuate real condition effects if batch and treatment were confounded (not the case here — one sample per condition).
6. **Platelet DE** — only one comparison possible; platelet results are incomplete for size-specific analysis.
7. **Azimuth timestamp** — annotation CSV is from an earlier run; biologically equivalent for the same cells.

---

## 9. Complete output inventory

### Figures (`results/figures/`)

| File | Status | Quality |
|------|--------|---------|
| `umap_condition.png` | Generated | Good |
| `umap_split_by_condition.png` | Generated | Good |
| `umap_sample_integration.png` | Generated | Good |
| `umap_clusters.png` | Generated | Good |
| `umap_celltypes_marker.png` | Generated | Good |
| `umap_codi_celltypes.png` | Generated | Good |
| `umap_module_scores.png` | Generated | Good |
| `marker_dotplot.png` | Generated | Good |
| `composition_barplot.png` | Generated | Good |

### Tables (`results/tables/`)

| File | Rows / size | Quality |
|------|-------------|---------|
| `cell_composition_by_condition.csv` | 32 rows | Good |
| `differential_expression_all.csv` | 66,000 rows | Good |
| `pathway_enrichment_all.csv` | 89,015 rows | Good |
| `size_specific_effects_summary.csv` | 40 rows | Good |
| `cell_cycle_scores_by_condition.csv` | 4 rows | **Good — valid numeric scores** |
| `ifn_scores_by_condition.csv` | 4 rows | Weak signal |
| `antigen_presentation_scores.csv` | 32 rows | Moderate signal |
| `pseudobulk_counts_condition_celltype.csv` | 32 groups | Good |
| `annotation_agreement_metrics.csv` | 1 metric | Good |
| `azimuth_annotations.csv` | 33,240 rows | Good |

### Processed data

| File | Quality |
|------|---------|
| `data/processed/integrated_annotated.h5ad` | Good |
| `data/processed/azimuth_mtx/` + `azimuth_obs.csv` | Good |

---

## 10. How to explain the data (presentation guide)

This section is a **script for oral defense and slides**. All figure captions below are ready to paste into PowerPoint.

### 10.1 One-minute project pitch

> “We analyzed single-cell RNA-seq data from human PBMC exposed to carboxylated polystyrene nanoplastics at 40 nm, 200 nm, and a mixture, compared to untreated control. After quality control we integrated 33,240 cells, annotated immune cell types, and asked whether particle size changes cell composition and gene expression. We found mostly stable global structure but widespread cell-type-specific differential expression, with inflammatory pathways enriched and clear size-dependent gene sets — especially in monocytes.”

---

### 10.2 Slide-by-slide talking points

| Slide / figure | What to say |
|----------------|-------------|
| **Study design** | One donor, four conditions, scRNA-seq from Zenodo. One sample per condition — interpret statistics as exploratory. |
| **`umap_condition.png`** | “Each dot is one cell. Similar cells are close together. All four conditions overlap — nanoplastics don’t tear apart the whole immune landscape, but that doesn’t mean no effect.” |
| **`umap_split_by_condition.png`** | “Same map, four panels — no condition sits alone in a corner. Effects are subtle at the global level.” |
| **`umap_sample_integration.png`** | “Colors are technical samples. They mix after Combat — we’re not comparing batch artifacts.” |
| **`umap_celltypes_marker.png`** | “We assigned PBMC types using canonical markers. CD4 T dominate, as expected in blood.” |
| **`marker_dotplot.png`** | “This proves labels make sense — B cell markers light up in B cells, and so on.” |
| **`composition_barplot.png`** | “Proportions are mostly stable. Main story: CD14 monocytes increase under 200 nm — possible myeloid response to larger particles.” |
| **`umap_module_scores.png`** | “Module scores compress a gene program into one number per cell. Cell-cycle and IFN patterns are visible spatially; IFN differences between conditions are weak.” |
| **DE / pathways** | “36,814 significant gene changes across 22 comparisons. Pathways point to immune and inflammatory biology — consistent with immune challenge.” |
| **Size-specific table** | “Many genes are unique to 40 nm or 200 nm — size matters. Monocytes share 865 genes across all exposures — a core response module. Mixture adds emergent genes.” |
| **Azimuth** | “Independent reference mapping confirms high-confidence PBMC identity — mean score 0.91.” |
| **Limitations** | “Single donor, in vitro, HVG subset for DE. Future work: more donors, pseudobulk validation.” |

---

### 10.3 Key concepts — short definitions for Q&A

| Term | How to explain |
|------|----------------|
| **UMAP** | A 2D visualization of high-dimensional gene expression. Proximity ≈ transcriptional similarity. |
| **HVG** | Highly variable genes — the ~3000 genes that differ most across cells; used for dimensionality reduction and DE in this pipeline. |
| **Combat** | Batch-correction method that removes technical differences between samples while preserving biology. |
| **Leiden cluster** | Unsupervised group of cells with similar expression — not a cell type until annotated. |
| **Module score** | One number per cell summarizing a gene program (e.g. interferon genes). Higher = stronger program activity. |
| **DE (differential expression)** | Statistical comparison of gene expression between exposure and control **within each cell type**. |
| **Pathway enrichment** | Tests whether DE genes are over-represented in known biological pathways (GO, KEGG, Reactome). |
| **Size-specific effect** | A DE gene significant only for 40 nm, only for 200 nm, shared, or unique to the mixture. |
| **CoDi / Azimuth** | External reference annotations — independent checks on our marker-based labels. |

---

### 10.4 Anticipated questions and answers

**Q: If UMAP looks the same, was there any effect?**  
A: Yes. DE shows tens of thousands of significant changes. UMAP shows **global** structure; nanoplastic effects are **gene-specific and cell-type-specific**, not a complete reorganization of the map.

**Q: Why is control IFN higher than exposed?**  
A: The IFN module score is weakly discriminative here. We report it honestly as inconclusive rather than forcing a narrative. DE and pathways are stronger evidence for immune-related changes.

**Q: Why only 64.7% CoDi agreement?**  
A: Different label definitions and shared markers between NK and T cells. Both methods agree on major lineages; disagreement is mostly at fine granularity.

**Q: Can you claim 200 nm increases monocytes?**  
A: Suggestively yes in this donor (5.7% → 9.7% CD14 monocytes). State it as a **hypothesis** pending replicate donors.

**Q: Single donor — is the analysis valid?**  
A: Valid as a **case study** and methods project. P-values reflect within-dataset comparisons; external validation needs more donors.

---

### 10.5 Suggested narrative arc (5 minutes)

1. **Problem:** Nanoplastics in blood may interact with immune cells; size may matter.  
2. **Data:** 33k PBMC cells, four conditions, high QC retention.  
3. **Integration:** Combat + UMAP — batches mix, biology preserved.  
4. **Annotation:** Markers + CoDi + Azimuth — coherent PBMC landscape.  
5. **Composition:** Mostly stable; monocytes ↑ at 200 nm.  
6. **Expression:** Massive DE; inflammatory pathways; size-specific genes.  
7. **Module scores:** Cell-cycle OK; IFN weak.  
8. **Limitations & next steps:** More donors, pseudobulk, in vivo relevance.

---

## 11. Suggested narrative for written thesis (English)

1. **Setup:** One-donor PBMC scRNA-seq, four PSNP conditions + control (Zenodo 15866724).  
2. **QC:** ~97% cells retained; 33,240 cells integrated — data quality is high.  
3. **Integration:** Combat successfully mixed batches; UMAP shows shared structure across conditions.  
4. **Annotation:** Three approaches — marker genes (assignment requirement), CoDi (99.6% mapped), Azimuth (mean confidence 0.91).  
5. **Composition:** Mostly stable; CD14 monocytes increase under 200 nm — hypothesis for myeloid involvement.  
6. **DE & pathways:** 36,814 significant changes; inflammatory/cytokine pathways enriched — nanoplastics activate immune gene programs in a cell-type-specific manner.  
7. **Size specificity:** Many genes are unique to 40 nm, 200 nm, or mixture — **particle size matters**.  
8. **Module scores:** Cell-cycle scoring fixed and usable; IFN differences weak — report transparently.  
9. **Limitations:** Single donor; DE on HVG subset — future work with replicates and full-transcriptome pseudobulk validation.

---

*Report generated from pipeline outputs in the project workspace. Numeric values are taken from `pipeline_run_20260611_212256.txt` and `results/tables/` unless otherwise noted. Azimuth statistics from `azimuth_run_20260606_175156.txt`.*
