# Marker Gene Selection Report

**Project:** Single-Cell Analysis of Immune Response to Nanoplastic Particles  
**Data:** Human PBMC from one donor, four conditions (40 nm PSNP, 200 nm PSNP, mixture, control)  
**Source:** Zenodo record [15866724](https://zenodo.org/records/15866724) (DOI: 10.5281/zenodo.15866724)

---

## 1. Purpose of marker genes in this project

The project assignment requires **cell type annotation** of immune populations (T cells, B cells, NK cells, monocytes, etc.) using **marker genes**, alongside validation with the **Azimuth PBMC reference** and comparison to **CoDi** labels supplied with the dataset.

The genes listed in `config/config.yaml` under `markers:` are **not learned by machine learning from the four samples**. They are **curated canonical markers** widely used in human peripheral blood mononuclear cell (PBMC) single-cell studies. Each gene (or small set) is expressed preferentially in one lineage or subtype, so average expression across a marker set acts as a simple score for assigning cell types.

In the pipeline (`scripts/run_pipeline.py`, function `marker_based_annotation`):

1. For each cell and each defined cell type, the mean expression of that type’s marker genes (present in the data matrix) is computed.
2. The cell is assigned the label with the **highest** mean score.
3. The result is stored in `adata.obs["cell_type_marker"]`.

That label is then used for:

- **Composition analysis** — cell type proportions across the four exposure conditions  
- **Differential expression (DE)** — comparisons of each exposed condition vs control **within** each major cell type  
- **Pathway enrichment** — GO / KEGG / Reactome on DE genes per cell type  
- **Size-specific interpretation** — which responses are unique to 40 nm, 200 nm, mixture, or shared  
- **Figures** — e.g. UMAP coloured by `cell_type_marker`  
- **Validation** — agreement with CoDi CSV labels and (optionally) Azimuth PBMC annotations  

Marker-based labels are one of **three** labelling approaches in this repository; they fulfil the assignment’s “marker genes” requirement directly, while Azimuth and CoDi provide independent checks.

---

## 2. Why curated markers instead of learning from this dataset?

| Reason | Explanation |
|--------|-------------|
| **Biology is known** | PBMC lineages are well characterised; standard marker panels exist in atlases and tutorials (e.g. Seurat PBMC, Azimuth). |
| **Small experimental design** | Four samples from one donor are insufficient to reliably discover de novo cell-type markers by unsupervised learning alone. |
| **Assignment wording** | The task asks for annotation *using marker genes* and the *Azimuth PBMC reference*, not training a new classifier on these `.h5ad` files. |
| **Interpretability** | Curated markers are easy to justify in a report, slides, and viva. |
| **CoDi as external reference** | Precomputed CoDi labels allow comparison without replacing the marker-based method required in the brief. |

Unsupervised **Leiden clusters** identify groups by expression similarity; they do not assign biological names. Marker genes bridge clusters to interpretable immune types for downstream DE and composition.

---

## 3. General selection principles

1. **Lineage specificity** — genes enriched in one PBMC compartment (T, B, NK, myeloid, platelet).  
2. **Detection in scRNA-seq** — reasonably abundant in PBMC datasets (not only protein-level markers).  
3. **Compatibility with CoDi / Azimuth categories** — labels map to coarse types used in `load_codi_annotations` (e.g. CD14+ monocyte → `Monocyte_CD14`).  
4. **Small panels** — multiple genes per type reduce noise from a single dropout-prone gene.  
5. **Overlap handled by scoring** — e.g. `NKG7` appears in both cytotoxic T and NK sets; assignment uses **highest mean score** across types (see limitations in Section 6).

**Note:** After integration, the working matrix is restricted to ~3000 highly variable genes (HVGs). A marker counts only if it remains in `adata.var_names`. Rare dropout of a marker gene slightly weakens that type’s score but does not change the biological rationale for including it in the config.

---

## 4. Marker sets: gene-by-gene rationale

### 4.1 T_cell — pan–T cell markers

| Gene | Role | Why included |
|------|------|----------------|
| **CD3D** | T cell receptor (TCR) complex subunit | Canonical surface/signalling complex for mature T cells; strong pan-T signal in blood. |
| **CD3E** | TCR complex subunit | Paired with CD3D; redundant markers stabilise the score if one gene is lowly detected. |
| **TRBC1** | TCR beta constant region | Transcript for TCR β chain; high in T cells, low in B/NK/myeloid. |

**Project role:** Defines the broad **T cell** compartment before splitting CD4 vs CD8 cytotoxic signatures.

---

### 4.2 CD4_T — CD4+ helper / memory T cells

| Gene | Role | Why included |
|------|------|----------------|
| **IL7R** | CD127, IL-7 receptor | Classic CD4+ T / memory T marker in PBMC; low on many cytotoxic CD8+ cells. |
| **LTB** | Lymphotoxin β | Expressed in CD4+ T and related lymphoid cells; supports CD4-associated identity. |
| **MALAT1** | Long non-coding RNA, nuclear | Highly expressed in many lymphocytes; boosts scoring in activated/structured T populations (used here as part of a CD4-associated panel, not lineage-specific alone). |

**Project role:** Separates **CD4+ T** from CD8 cytotoxic and other lineages for composition and DE within helper T cells.

---

### 4.3 CD8_T_cytotoxic — cytotoxic T cells

| Gene | Role | Why included |
|------|------|----------------|
| **NKG7** | Natural killer granule protein 7 | Cytotoxic lymphocyte effector programme (shared with NK at transcript level). |
| **GNLY** | Granulysin | Cytotoxic granule content; effector T / NK-associated. |
| **PRF1** | Perforin-1 | Key cytotoxic effector molecule. |
| **GZMB** | Granzyme B | Cytotoxic serine protease; strong in CD8+ effector cells. |

**Project role:** Identifies **cytotoxic T** responses to nanoplastics (e.g. activation, cytotoxicity pathways in DE). Overlap with NK markers is expected; scoring distinguishes the winning label per cell.

---

### 4.4 B_cell — B lymphocytes

| Gene | Role | Why included |
|------|------|----------------|
| **MS4A1** | CD20 | Standard B cell surface marker gene in blood. |
| **CD79A** | B cell receptor component | BCR signalling; B-lineage specific. |
| **CD74** | MHC class II invariant chain | High in B cells and antigen-presenting cells; here supports B cell identity in a B-focused panel. |

**Project role:** Quantifies **B cell** abundance per condition and enables B-cell-specific DE vs control.

---

### 4.5 NK_cell — natural killer cells

| Gene | Role | Why included |
|------|------|----------------|
| **NKG7** | NK granule protein | Core NK / cytotoxic lymphocyte marker. |
| **KLRD1** | CD94 | NK receptor complex; NK-enriched vs classical T. |
| **TRAC** | TCR α constant | Low in true NK (lack rearranged TCR programme); included to contrast with T cells—NK typically score high on NKG7/KLRD1 and relatively lower on pan-TCR programme when T markers are separated in other sets. |

**Project role:** Captures **innate cytotoxic** compartment, relevant for inflammatory responses to particles.

---

### 4.6 Monocyte_CD14 — classical CD14+ monocytes

| Gene | Role | Why included |
|------|------|----------------|
| **LYZ** | Lysozyme | Strong myeloid/monocyte signal. |
| **S100A8** | Calprotectin subunit | Classical monocyte / inflammatory myeloid. |
| **S100A9** | Calprotectin subunit | Paired with S100A8; acute inflammatory monocytes. |
| **CTSS** | Cathepsin S | Lysosomal protease enriched in monocytes / APCs. |

**Project role:** Central for the biological question (e.g. monocyte activation by 40 nm vs 200 nm PSNP); DE and pathway steps often focus on this population.

---

### 4.7 Monocyte_CD16 — non-classical CD16+ monocytes

| Gene | Role | Why included |
|------|------|----------------|
| **FCGR3A** | CD16 | Defining marker for non-classical / patrolling monocytes. |
| **LST1** | Leukocyte-specific transcript 1 | Monocyte subset marker. |
| **MS4A7** | Membrane protein | Enriched in CD16+ monocyte cluster in standard PBMC references. |

**Project role:** Splits **CD14+ vs CD16+** monocyte responses; size-specific effects may differ between subsets.

---

### 4.8 DC — dendritic cells

| Gene | Role | Why included |
|------|------|----------------|
| **FCER1A** | FcεRI α | Conventional dendritic cell (cDC) marker in blood. |
| **CST3** | Cystatin C | Expressed in dendritic cells and some monocytes; paired with FCER1A for DC scoring. |

**Project role:** Rare but biologically important **antigen-presenting** population for nanoplastic–immune interaction; included for completeness of PBMC annotation.

---

### 4.9 Platelet — platelets / platelet contaminants

| Gene | Role | Why included |
|------|------|----------------|
| **PPBP** | Pro-platelet basic protein | Highly specific platelet transcript in scRNA-seq. |
| **PF4** | Platelet factor 4 | Classic platelet chemokine; identifies platelet RNA in PBMC preps. |

**Project role:** Flags **platelet** signal (often ambient RNA or doublets in PBMC); useful for QC interpretation and excluding mis-assigned lymphocytes from platelet-rich droplets.

---

## 5. How markers connect to project tasks

| Task | Marker role |
|------|-------------|
| **Cell type annotation** | Primary implementation via `cell_type_marker` from config panels. |
| **Composition analysis** | Counts and fractions by `cell_type_marker` × `condition`. |
| **Differential expression** | Subsets cells by `cell_type_marker` before exposure vs control tests. |
| **Pathway enrichment** | Applied to DE gene lists **per annotated cell type**. |
| **Size-specific effects** | Interpreted within cell types defined by markers (and cross-checked with DE tables). |
| **Comparison to CoDi** | `codi_to_marker` mapping aligns external labels to the same coarse types. |
| **Azimuth (optional)** | Finer PBMC labels from `ref.Rds`; report should compare Azimuth vs `cell_type_marker`. |

---

## 6. Limitations (to state in the thesis / slides)

1. **Winner-takes-all scoring** — each cell gets exactly one label; mixed or transitional states are not modelled.  
2. **Marker overlap** — `NKG7` in both CD8 cytotoxic and NK sets; biology overlaps; highest score decides.  
3. **HVG subset** — annotation runs on genes retained after HVG selection; some markers may be absent from the matrix.  
4. **Processed expression** — scores use normalized/log/scaled data in the integrated object, not raw counts.  
5. **Not a substitute for Azimuth** — assignment requires both marker-based and reference-based annotation; markers alone do not replace the RDS workflow.

---

## 7. References and further reading (for citation in report)

- Hao et al., *Integrated analysis of multimodal single-cell data*, Cell (2021) — Azimuth / PBMC reference framework.  
- Satija et al. / Seurat PBMC tutorials — standard marker genes for human blood.  
- Dataset CoDi labels (Zenodo 15866724) — external cell-type predictions for validation.  
- Classic PBMC markers: CD3 (T), MS4A1/CD79A (B), NKG7/KLRD1 (NK), LYZ/S100A8/S100A9 (CD14 monocytes), FCGR3A (CD16 monocytes), FCER1A (DC), PPBP/PF4 (platelets).

---

## 8. Summary

The marker genes in `config.yaml` were chosen because they are **established, interpretable indicators** of major human PBMC lineages relevant to **blood–nanoplastic immune interactions**. In this project they **assign each cell a type** for composition and DE, implement the coursework requirement for **marker-based annotation**, and complement **Azimuth** and **CoDi** as independent references—not because they were statistically discovered de novo from the four PSNP samples alone.

For reproduction: edit markers in `config/config.yaml`, run `python scripts/run_pipeline.py`, and inspect `adata.obs["cell_type_marker"]` and outputs under `results/tables/` and `results/figures/`.
