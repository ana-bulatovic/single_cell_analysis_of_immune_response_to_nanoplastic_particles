# Dodatne analize — interpretacija rezultata

Ovaj fajl objašnjava šta znače tabele i figure iz `run_additional_analyses.py` / STEP 9 pipeline-a.

## 1. Module score analize

**Šta se radi:** Na punom normalizovanom transkriptomu (pre HVG filtriranja) računaju se
skorovi genetskih programa pomoću `scanpy.tl.score_genes` / `score_genes_cell_cycle`.
Svaka ćelija dobija jedan broj po programu; ovde se računaju proseci po uslovu i tipu ćelije.

**Geni u panelima:**
- S / G2M faza: klasični cell-cycle geni (MCM5, PCNA, MKI67, ...)
- IFN: ISG15, IFIT1-3, MX1, OAS1, OASL, IFI6, RSAD2
- Antigen presentation: HLA-DRA/DRB1, CD74, B2M, TAP1/2, HLA-DPA1/DPB1

### Globalni nalazi po uslovu

**S_score**
  - S_score | PSNP 40 nm: 0.0049 (Δ vs control +0.0020) — bez značajne promene
  - S_score | PSNP 200 nm: 0.0005 (Δ vs control -0.0024) — bez značajne promene
  - S_score | PSNP mix: 0.0044 (Δ vs control +0.0015) — bez značajne promene

**G2M_score**
  - G2M_score | PSNP 40 nm: 0.0163 (Δ vs control +0.0044) — bez značajne promene
  - G2M_score | PSNP 200 nm: -0.0177 (Δ vs control -0.0296) — snižen u odnosu na control
  - G2M_score | PSNP mix: 0.0122 (Δ vs control +0.0003) — bez značajne promene

**IFN_score**
  - IFN_score | PSNP 40 nm: 0.0045 (Δ vs control -0.0306) — snižen u odnosu na control
  - IFN_score | PSNP 200 nm: 0.0096 (Δ vs control -0.0255) — snižen u odnosu na control
  - IFN_score | PSNP mix: -0.0018 (Δ vs control -0.0369) — snižen u odnosu na control

**antigen_presentation_score**
  - antigen_presentation_score | PSNP 40 nm: 0.2979 (Δ vs control -0.0676) — snižen u odnosu na control
  - antigen_presentation_score | PSNP 200 nm: 0.3360 (Δ vs control -0.0294) — snižen u odnosu na control
  - antigen_presentation_score | PSNP mix: 0.3197 (Δ vs control -0.0457) — snižen u odnosu na control

### Antigen presentation po tipu ćelije

Najviši skorovi su očekivani kod B ćelija i antigen-prezentujućih tipova (DC, monociti).
Poređenje uslova unutar istog tipa pokazuje da li izloženost menja MHC/HLA program.

**B cell**
  - PSNP 40 nm: 1.451 (Δ -0.220 vs control)
  - PSNP 200 nm: 1.618 (Δ -0.054 vs control)
  - PSNP mix: 1.523 (Δ -0.148 vs control)

**Monocyte CD14**
  - PSNP 40 nm: 0.059 (Δ -0.180 vs control)
  - PSNP 200 nm: 0.104 (Δ -0.134 vs control)
  - PSNP mix: 0.039 (Δ -0.200 vs control)

**DC**
  - PSNP 40 nm: 1.647 (Δ -0.037 vs control)
  - PSNP 200 nm: 1.764 (Δ +0.081 vs control)
  - PSNP mix: 1.916 (Δ +0.232 vs control)

## 2. Pseudobulk matrica

**Šta se radi:** Za svaku kombinaciju `condition × cell_type` sabiraju se sirovi UMI
brojevi po genu (iz `layers['counts']`).

**Šta znači:** Omogućava bulk RNA-seq stil analizu (DESeq2, edgeR) bez ponovnog
učitavanja pojedinačnih ćelija. Broj grupa = broj uslova × broj tipova ćelija.

## 3. Validacija anotacije

**Šta se radi:** Porede se tri nezavisna izvora tipova ćelija:
marker panel (pipeline), CoDi (Zenodo CSV), Azimuth PBMC referenca (opciono).

| Metrika | Vrednost | Značenje |
|---------|----------|----------|
| ref_marker_agreement | 0.454 | Fraction of cells where ref.Rds primary label matches literature markers |
| codi_marker_agreement | 0.630 | Fraction of CoDi-mapped cells with same label as literature markers |
| codi_ref_agreement | 0.599 | Fraction of CoDi-mapped cells with same label as ref.Rds primary |
| codi_mapped_fraction | 0.990 | Fraction of cells with a CoDi reference label |
| azimuth_marker_agreement | 0.456 | Fraction of Azimuth-mapped cells matching marker annotation (L1 mapped) |
| azimuth_marker_agreement_score_ge_0.5 | 0.458 | Same as above, only cells with Azimuth score >= 0.5 |
| codi_azimuth_agreement | 0.600 | Fraction where CoDi normalized label matches Azimuth L1 mapped label |
| azimuth_mean_score | 0.908 | Mean Azimuth prediction score across all cells |

**Kako čitati:**
- 100% slaganje je retko — NK/CD8 i granularnost tipova razlikuju metode.
- ~45–65% je tipično za PBMC sa više anotacionih šema.
- Azimuth mean score > 0.85 = visoko poverenje u referentno mapiranje.

**Contingency matrix (CoDi vs markeri):**
- `annotation_crosstab_marker_codi.csv` — broj ćelija po paru tipova
- `annotation_crosstab_marker_codi_row_pct.csv` — % unutar svakog marker tipa
- `annotation_codi_marker_mapping.csv` — long format (marker → CoDi mapiranje)

**Contingency matrix (ref.Rds/Seurat vs CoDi):**
- `annotation_crosstab_ref_codi.csv` — broj ćelija po paru tipova
- `annotation_crosstab_ref_codi_row_pct.csv` — % unutar svakog ref.Rds tipa
- `annotation_codi_ref_mapping.csv` — long format (ref.Rds → CoDi mapiranje)

## 4. Figure

- `results/figures/additional_analyses/module_scores_by_condition.png` — bar chart modula
- `results/figures/additional_analyses/module_scores_violin.png` — raspodela po uslovu
- `results/figures/additional_analyses/antigen_presentation_heatmap.png`
- `results/figures/additional_analyses/annotation_agreement_bar.png`
- `results/figures/additional_analyses/annotation_confusion_marker_codi.png` — CoDi vs markeri
- `results/figures/additional_analyses/annotation_confusion_marker_codi_normalized.png`
- `results/figures/additional_analyses/annotation_confusion_ref_codi.png` — ref.Rds vs CoDi
- `results/figures/additional_analyses/annotation_confusion_ref_codi_normalized.png`
- `results/figures/additional_analyses/annotation_confusion_marker_azimuth.png` (ako postoji Azimuth)
