# Biological Interpretation of Size-Specific Transcriptional Effects

**Project:** Single-Cell Analysis of Immune Response to Nanoplastic Particles  
**Analysis step:** STEP 8 — size-specific effect classification (`scripts/run_pipeline.py`)  
**Primary annotation:** Seurat/Azimuth L1 (`cell_type_primary`)  
**Source tables:** `results/tables/size_specific_effects_summary.csv`, `size_specific_interpretation.csv`, `size_specific_pathway_enrichment.csv`  
**Pipeline run:** `pipeline_run_20260619_003132.txt`

---

## 1. Purpose and rationale

Polystyrene nanoparticles (PSNP) of different diameters differ in **surface-area-to-volume ratio**, **cellular uptake route**, and **interaction with pattern-recognition receptors (PRRs)** and phagocytic machinery. To move beyond “PSNP exposure changes gene expression” and ask whether **particle size shapes the immune response**, we classified significantly differentially expressed (DE) genes into five mutually exclusive effect classes per cell type:

| Effect class | Definition |
|--------------|------------|
| **unique_40nm** | Significant in 40 nm vs control, but not in 200 nm or mix |
| **unique_200nm** | Significant in 200 nm vs control, but not in 40 nm or mix |
| **shared_40_200** | Significant in both solo exposures, but **not** in the mixture |
| **shared_all_three** | Significant in 40 nm, 200 nm, **and** mix — a core, size-independent PSNP response |
| **mix_only_emergent** | Significant only in the 40+200 nm mixture — not explained by either solo exposure |

Genes were taken from per-cell-type Wilcoxon DE (adjusted *p* < 0.05, |log₂FC| > 0.25) across four conditions: control, 40 nm, 200 nm, and 40+200 nm mix. Pathway enrichment (Enrichr: GO, KEGG, Reactome) was run on UP-regulated genes within each class.

**Biological question addressed:** Do 40 nm and 200 nm PSNP trigger the same transcriptional programmes in human PBMC, and does co-exposure to both sizes produce additive, redundant, or **emergent** effects?

---

## 2. Overview of results

Across **seven Seurat L1 cell types** (B cells, CD4 T, CD8 cytotoxic T, CD14 monocytes, NK cells, *other*, *other_T*), every effect class contained non-zero gene sets, confirming that **particle size is not a neutral experimental detail** — it partitions the transcriptomic response.

Three global patterns emerge:

1. **A conserved “core PSNP module” (`shared_all_three`)** is present in all cell types, with the largest gene set in **CD14 monocytes** (973 genes). Pathway enrichment consistently points to **cytokine–cytokine receptor interaction**, **IL-10 / interleukin signalling**, **TNF signalling**, and **cellular response to LPS** — a coherent picture of **innate immune and inflammatory reprogramming** triggered by nanoplastic exposure regardless of particle diameter.

2. **Size-specific modules (`unique_40nm`, `unique_200nm`)** are substantial (typically 150–400 genes per class per type), indicating that **40 nm and 200 nm are not interchangeable** at the transcriptional level. Smaller particles enrich pathways linked to **high surface-area signalling** (e.g. MAPK, AGE–RAGE, calcium transport, VEGF/platelet activation in NK cells), whereas 200 nm-specific sets more often align with **myeloid-type engagement** (phagocytosis, FCGR3A-mediated uptake, chemokine and IL-17/TNF cascades).

3. **The mixture is not a simple sum of solo exposures.** Genes significant in both solo conditions but absent in the mix (`shared_40_200`) suggest **masking or redirection** when both sizes are present. Conversely, **`mix_only_emergent`** genes (223–287 per type) reveal **non-additive mixture biology** — programmes activated only when both particle sizes co-occur.

---

## 3. Cell-type-specific interpretation

### 3.1 CD14 monocytes — primary innate sentinels

Monocytes carry the **largest core response module** (973 genes in `shared_all_three`), with the strongest pathway statistics in the entire analysis:

- **Reactome Interleukin-10 signalling** (padj = 1.37×10⁻¹⁰)
- **KEGG Cytokine–cytokine receptor interaction** (padj = 4.14×10⁻⁹)
- **Reactome Signalling by interleukins** (padj = 7.02×10⁻⁸)
- **GO Cellular response to lipopolysaccharide** (padj = 3.49×10⁻⁷)

**Interpretation:** CD14⁺ monocytes are the cell type most comprehensively reprogrammed by PSNP exposure. The enrichment profile resembles **trained innate immunity / endotoxin-response-like transcriptional states**, consistent with monocytes acting as first-line sensors of particulate stress. This aligns with composition analysis showing a **relative increase in CD14 monocytes under 200 nm PSNP**.

The **`mix_only_emergent`** monocyte module (234 genes) enriches for **Th17 cell differentiation**, **graft-versus-host disease**, **cell adhesion molecules**, and **intestinal immune network for IgA production** — suggesting that the mixture may push monocytes toward a **T-cell polarisation and tissue-immune crosstalk** phenotype not seen with either particle size alone.

Size-unique monocyte gene sets (196 for 40 nm; 157 for 200 nm) did not yield strong standalone pathway terms at FDR < 0.05, but their existence confirms **quantitatively distinct monocyte responses** to each diameter even when pathway-level signal is diffuse.

---

### 3.2 CD4 T cells — cytokine-driven adaptive reshaping

CD4 T cells show a balanced distribution across effect classes (316 unique 40 nm; 244 unique 200 nm; 511 shared solo; 582 core; 264 mix-only).

**`shared_40_200` (511 genes)** — significant under both solo exposures but lost in the mix — is strongly enriched for:

- **Reactome Neutrophil degranulation** (padj = 1.71×10⁻⁷)
- **GO Cytokine-mediated signalling** (padj = 1.81×10⁻⁶)
- **Reactome Innate immune system** (padj = 8.51×10⁻⁵)

**Interpretation:** Solo PSNP exposures activate in CD4 T cells a **cytokine/inflammatory programme** that may be driven indirectly by myeloid-derived signals (neutrophil degranulation is a hallmark of innate–adaptive crosstalk). The **absence of this module in the mixture** implies that co-exposure **suppresses or replaces** this solo-size programme — an important non-additive finding for mixture toxicity assessment.

**`shared_all_three` (582 genes)** enriches for **viral protein–cytokine interaction**, **IL-10 signalling**, and **TNF / rheumatoid arthritis** KEGG terms — a **conserved stress/inflammatory axis** maintained across all exposure arms.

**`unique_200nm`** links to **AGE–RAGE signalling**, **chemokine signalling**, and **fluid shear stress / atherosclerosis** — pathways associated with **endothelial–immune interface stress and myeloid recruitment**, consistent with larger particles engaging phagocytic and danger-sensing routes.

**`mix_only_emergent`** enriches **ECM–receptor interaction**, **focal adhesion**, and **PI3K–Akt signalling**, suggesting mixture-specific **adhesion and survival signalling** in CD4 T cells.

---

### 3.3 CD8 cytotoxic T cells — innate-like activation in cytotoxic compartment

CD8 T cells mirror CD4 patterns in **`shared_40_200`** (537 genes): **neutrophil degranulation**, **interleukin signalling**, and **innate immune system** — again pointing to **indirect innate immune imprinting** on cytotoxic T cells under solo PSNP exposure, attenuated in the mixture.

**`shared_all_three` (526 genes)** enriches:

- **GO Cellular response to LPS** (padj = 3.95×10⁻³)
- **GO Regulation of macrophage activation** (padj = 4.58×10⁻³)
- **GO Positive regulation of IL-2 production** (padj = 1.22×10⁻²)

**Interpretation:** Even in cytotoxic T cells — not classical phagocytes — PSNP exposure induces a **macrophage/monocyte-coupled inflammatory state**. This supports a model where **myeloid sensing of particles propagates transcriptional changes across lymphocyte compartments**, rather than each T cell detecting particles autonomously.

**`unique_200nm`** shows **cytokine–cytokine receptor interaction**, reinforcing size-dependent immune signalling in the CD8 compartment.

---

### 3.4 B cells — humoral arm with strong 200 nm-specific inflammation

B cells exhibit one of the largest **`shared_all_three`** modules (689 genes), indicating broad PSNP sensitivity in the humoral lineage.

**`unique_200nm` (296 genes)** enriches:

- **KEGG IL-17 signalling** (padj = 1.13×10⁻²)
- **KEGG TNF signalling** (padj = 1.28×10⁻²)
- **KEGG Legionellosis** (padj = 1.24×10⁻²)

**Interpretation:** 200 nm PSNP specifically engages **T helper 17 / TNF-axis inflammatory pathways** in B cells — a pattern typical of **T-dependent B-cell activation in inflammatory microenvironments**, possibly reflecting stronger myeloid co-stimulation with larger particles.

**`shared_40_200` (385 genes)** enriches **inflammatory response**, **cytokine-mediated signalling**, and **IL-10 signalling** — solo-size inflammatory B-cell programmes again **absent in the mixture**.

**`mix_only_emergent` (228 genes)** includes **ITGAM** (CD11b) among top genes and calcitonin-receptor-related Reactome terms, hinting at **myeloid–B-cell interface molecules** emerging only under mixed exposure.

---

### 3.5 NK cells — cytotoxic innate lymphocytes with distinct size fingerprints

**`unique_40nm` (208 genes)** enriches **VEGF signalling**, **platelet activation**, **Fc epsilon RI signalling**, and **osteoclast differentiation** — pathways tied to **degranulation, Fc-receptor signalling, and tissue-remodelling cues**, consistent with small particles triggering **high surface-area receptor clustering**.

**`unique_200nm` (272 genes)** enriches **cytokine–cytokine receptor interaction**, **JAK–STAT signalling**, and **positive regulation of NF-κB** — a more classical **cytokine-driven NK activation** profile associated with larger-particle / myeloid-context exposure.

**`shared_40_200` (391 genes)** shows **inflammatory response** among top GO terms, again a solo-exposure programme not retained in the mix.

The large **`shared_all_three`** set (672 genes) lacks strong FDR-passing pathway terms in the summary table, but its size indicates a **broad conserved NK transcriptional shift** under any PSNP condition.

---

### 3.6 *Other* and *other_T* — residual compartments

The Seurat L1 **“other”** cluster (likely heterogeneous myeloid/dendritic/progenitor-like cells) shows:

- **`unique_40nm`:** **AGE–RAGE**, **MAPK**, **regulation of macrophage cytokine production** — direct evidence of **small-particle danger signalling** in a myeloid-enriched compartment.
- **`unique_200nm`:** **FCGR3A-mediated phagocytosis**, **RAC1/CDC42 GTPase cycles**, **lipid and atherosclerosis** — textbook **phagocytic uptake and cytoskeletal remodelling** for larger particles.
- **`mix_only_emergent`:** **Axon guidance** and **MECP2-related transcriptional regulation** — less canonical for immunity; may reflect **stress-induced developmental/neuronal gene co-option** or heterogeneous cell-state mixing.

***other_T*** (miscellaneous T-lineage cells) shows **`shared_40_200`** enrichment for **cytokine-mediated signalling** and **IL-10**, and **`unique_200nm`** enrichment for **IL-17**, **collagen degradation**, and **ECM organisation** — reinforcing the **Th17 / matrix-remodelling** theme under 200 nm exposure.

---

## 4. Integrative biological model

Based on the size-specific classification, we propose the following model for PBMC response to carboxylated PSNP:

```
                    ┌─────────────────────────────────────┐
                    │   Core PSNP module (all exposures)   │
                    │  cytokine/IL-10/TNF/LPS-like innate  │
                    └─────────────────┬───────────────────┘
                                      │
              ┌───────────────────────┼───────────────────────┐
              ▼                       ▼                       ▼
       ┌─────────────┐         ┌─────────────┐         ┌─────────────┐
       │  40 nm solo │         │ 200 nm solo │         │ 40+200 mix  │
       │  high SA/V  │         │ phagocytosis│         │ non-additive│
       │ MAPK/RAGE/  │         │ FCGR3A/IL-17│         │ emergent +  │
       │ VEGF/FcR    │         │ chemokine   │         │ masking of  │
       └──────┬──────┘         └──────┬──────┘         │ solo modules│
              │                       │                └──────┬──────┘
              └───────────┬───────────┘                       │
                          ▼                                   ▼
              shared_40_200: innate/cytokine         mix_only: adhesion,
              programmes in T/B/NK                   Th17, ECM/PI3K
              (lost in mixture)                      (not in solo arms)
```

**Key biological conclusions for report/defense:**

1. **Nanoplastic immune effects are partly size-independent** — a shared cytokine/inflammatory module operates across 40 nm, 200 nm, and mixture in all major PBMC lineages, with **CD14 monocytes as the epicentre**.

2. **40 nm and 200 nm PSNP are biologically distinct exposures**, not merely dose-scaled versions of the same stimulus. Smaller particles favour **surface-area-dependent signalling** (MAPK, RAGE, VEGF/FcR); larger particles favour **phagocytic and IL-17/TNF myeloid programmes**.

3. **Mixture toxicity cannot be predicted from solo exposures alone.** The presence of large `shared_40_200` (solo-only) and `mix_only_emergent` (mixture-only) gene sets demonstrates **masking, redirection, and emergent transcriptional programmes** — critical for environmental risk assessment where humans encounter **polydisperse** nanoplastic mixtures.

4. **Innate–adaptive crosstalk is central:** neutrophil degranulation, macrophage activation, and IL-10/TNF/IL-17 pathways appear across lymphocyte types, suggesting **particle sensing is initiated in myeloid cells** and propagated as transcriptional state changes in T, B, and NK compartments.

---

## 5. Suggested wording for oral defense (short form)

> “We classified differentially expressed genes by exposure geometry — unique to 40 nm, unique to 200 nm, shared between solo sizes, shared across all conditions, or emergent only in the mixture. Monocytes showed the largest conserved inflammatory module, enriched for IL-10, cytokine receptors, and LPS-response pathways. Smaller particles preferentially activated surface-area-related signalling such as MAPK and AGE–RAGE, whereas 200 nm particles enriched phagocytosis and IL-17/TNF programmes. Importantly, the mixture was not additive: solo-exposure inflammatory modules disappeared in the mix, and new mixture-specific pathways — including Th17 differentiation and ECM/focal adhesion — appeared only when both particle sizes were present. This supports the conclusion that particle size and co-exposure pattern jointly shape the human immune transcriptional response to nanoplastics.”

---

## 6. Limitations

- **Single donor, no biological replicates:** effect classes are descriptive within-sample partitions; independent cohort validation is required for generalisation.
- **Carboxylated polystyrene model particles** may not fully represent environmental nanoplastics (mixed polymers, coatings, protein corona in vivo).
- **Pathway enrichment** reflects gene-set overlap; some top genes are uncharacterised lncRNAs or pseudogenes — functional follow-up (flow cytometry, cytokine panels, phagocytosis assays) would strengthen causal claims.
- **Seurat L1 annotation** collapses fine subtypes (e.g. no separate DC or CD16 monocyte class in this reference level), which may pool distinct biology into *other* / *other_T* clusters.

---

## 7. Supporting files and figures

| Output | Location |
|--------|----------|
| Gene counts per class | `results/tables/size_specific_effects_summary.csv` |
| Full gene lists | `results/tables/size_specific_genes.csv` |
| Pathway enrichment | `results/tables/size_specific_pathway_enrichment.csv` |
| Machine-readable summary | `results/tables/size_specific_interpretation.csv` |
| Bar / UpSet plots | `results/figures/supplementary/size_specific_*.png` |
| Pathway dot/bar plots per class | `results/figures/supplementary/size_pathways_*.png` |

---

*Document prepared for inclusion in the project report and oral defense. Generated from pipeline outputs (2026-06-19).*
