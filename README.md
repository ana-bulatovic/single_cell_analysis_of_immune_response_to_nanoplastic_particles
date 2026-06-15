# Single-Cell Analysis of Immune Response to Nanoplastic Particles

Ovaj repozitorijum implementira kompletan scRNA-seq workflow za dataset sa Zenodo zapisa [15866724](https://zenodo.org/records/15866724):
- Sample 1: **40 nm** PSNP
- Sample 2: **200 nm** PSNP
- Sample 3: **40 + 200 nm** (mix)
- Sample 4: **control**

Analiza pokriva: QC, integraciju (Combat), klasterizaciju, anotaciju, kompoziciju, DE, pathway obogaćivanje, size-specific efekte i **dodatne analize** (module score, pseudobulk, validacija anotacije).

---

## 1) Podešavanje okruženja

**Opcija A (preporučeno na Windows-u):**

```bash
cd single_cell_analysis_of_immune_response_to_nanoplastic_particles
python -m pip install -r requirements.txt
```

**Opcija B (Conda + R za Azimuth):**

```bash
conda env create -f environment.yml
conda activate nanoplastic-scRNA
```

> Na Windows-u pipeline koristi **Combat** integraciju (`integration_method: combat` u `config/config.yaml`), jer `harmonypy` često nije dostupan.

### Procena trajanja (tipičan laptop)

| Korak | Trajanje |
|-------|----------|
| `pip install -r requirements.txt` | 5–15 min (jednom) |
| `python scripts/download_data.py` (~800 MB) | 10–40 min |
| `python scripts/run_pipeline.py` | 30–90 min |
| `python scripts/run_additional_analyses.py` | 1–3 min |
| `python scripts/plot_supplementary_figures.py` | 1–3 min |
| `python scripts/make_slides.py` | < 1 min |
| `Rscript scripts/azimuth_annotation.R` (opciono) | 15–45 min |

---

## 2) Preuzimanje podataka

```bash
python scripts/download_data.py
```

Preuzima `.h5ad`, CoDi `.csv` i Azimuth `ref.Rds` u `data/raw/`.

---

## 3) Glavni pipeline

```bash
python scripts/run_pipeline.py
```

### Koraci pipeline-a

| Korak | Šta radi |
|-------|----------|
| 1 QC | Filtriranje ćelija (`min_genes`, `max_genes`, `min_counts`, `max_mt_percent`) |
| 2 Integracija | Normalizacija, **module score na punom transkriptomu**, HVG (3000), Combat, UMAP, Leiden |
| 3 Anotacija | Marker paneli (Azimuth-izvedeni ili `config.yaml`) + CoDi CSV |
| 4 Figure | UMAP, dotplot, module score mape |
| 5 Kompozicija | Proporcije tipova po uslovu |
| 6 DE | Wilcoxon po tipu ćelije, svaki exposure vs control |
| 7 Pathway | Enrichr (GO / KEGG / Reactome) + figure |
| 8 Size-specific | Klasifikacija gena po veličini čestice (40 / 200 / mix) |
| 9 **Dodatne analize** | Module score tabele, pseudobulk, CoDi/Azimuth validacija, figure, interpretacija |
| 10 Save | `data/processed/integrated_annotated.h5ad` |

### Anotacija ćelija

U `config/config.yaml`:
- `azimuth.use_panels_for_pipeline: true` — koristi `results/tables/azimuth_marker_panels_l1.yaml` (preporučeno posle Azimuth run-a)
- `false` — koristi klasične markere iz `config.yaml` (`IL7R`, `LYZ`, `MS4A1`, …)

---

## 4) Azimuth PBMC referenca (opciono, R)

```bash
python scripts/check_rscript.py
Rscript scripts/install_r_packages.R
python scripts/prepare_azimuth_h5ad.py
Rscript scripts/azimuth_annotation.R
python scripts/extract_azimuth_markers.py
```

Izlazi:
- `results/tables/azimuth_annotations.csv`
- `results/tables/azimuth_marker_panels_l1.yaml`
- `results/figures/azimuth_marker_dotplot_l1.png`

---

## 5) Dodatne analize (odvojeno pokretanje)

```bash
python scripts/run_additional_analyses.py
```

**Zahteva:** `data/processed/integrated_annotated.h5ad` (iz glavnog pipeline-a).

### Šta konkretno radi svaka analiza

| # | Analiza | Metoda | Output | Šta znači rezultat |
|---|---------|--------|--------|-------------------|
| 1 | **Cell-cycle module** | `scanpy.tl.score_genes_cell_cycle` na punom transkriptomu (S i G2M geni) | `cell_cycle_scores_by_condition.csv`, `module_scores_by_condition.csv` | Da li izloženost menja proliferaciju/stres ćelija. Skorovi blizu 0 = većina PBMC nije u aktivnom deljenju (očekivano za krv). |
| 2 | **IFN signature** | `score_genes` na ISG15, IFIT1-3, MX1, OAS1, … | `ifn_scores_by_condition.csv` | Innate interferon program. Viši skor = jača IFN aktivacija. Ako je control viši od izloženih, IFN nije glavni dokaz efekta. |
| 3 | **Antigen presentation** | `score_genes` na HLA-DRA, CD74, TAP1, … | `antigen_presentation_scores.csv`, `module_scores_by_condition_celltype.csv` | MHC/HLA program po tipu i uslovu. Relevantno za B, DC, monocite. |
| 4 | **Pseudobulk** | Zbir UMI po `condition × cell_type` | `pseudobulk_counts_condition_celltype.csv` | Bulk RNA-seq matrica za DESeq2/edgeR validaciju bez ponovnog učitavanja ćelija. |
| 5 | **Validacija anotacije** | Poređenje marker vs CoDi vs Azimuth | `annotation_agreement_metrics.csv`, `annotation_crosstab_*.csv` | Koliko se slažu nezavisne metode. 45–65% je tipično (NK/CD8 granica, različita granularnost). |
| 6 | **Interpretacija** | Automatski tekst na srpskom | `additional_analyses_interpretation_SR.md` | Objašnjenje svake tabele i figure uz brojeve iz trenutnog run-a. |
| 7 | **Figure** | bar / violin / heatmap / confusion | `results/figures/additional_analyses/` | Vizuelni prikaz modula i slažnosti anotacija za prezentaciju. |

**Napomena:** Module score-ovi se računaju u koraku 2 pipeline-a (pre HVG). Standalone skripta ih **čita iz** `.h5ad` — ne ponovo računa pun transkriptom. Ako fale, pokreni ponovo `run_pipeline.py`.

---

## 6) Dodatne skripte za figure

```bash
python scripts/plot_figures.py              # UMAP iz .h5ad
python scripts/plot_pathway_figures.py      # pathway heatmap iz CSV
python scripts/plot_supplementary_figures.py  # volcano, size-specific, DE summary iz CSV
```

---

## 7) Prezentacija i dokumentacija

```bash
python scripts/make_slides.py
```

Izlaz: `deliverables/nanoplastic_scRNA_results.pptx`

Detaljni izveštaj: `deliverables/Analysis_Results_Report.md`  
Scenario za video: `deliverables/video_script.md`  
Word dokumentacija: `deliverables/Project_Documentation.docx`

---

## 8) Glavni output fajlovi

### Podaci
- `data/processed/integrated_annotated.h5ad`

### Figure
- `results/figures/` — UMAP, composition, marker dotplot
- `results/figures/pathway_enrichment/` — pathway heatmap po tipu
- `results/figures/supplementary/` — DE volcano, size-specific (iz CSV)
- `results/figures/additional_analyses/` — module score i anotacija

### Tabele
- `differential_expression_all.csv`
- `pathway_enrichment_all.csv`
- `size_specific_effects_summary.csv`
- `module_scores_by_condition.csv`
- `annotation_agreement_metrics.csv`
- `additional_analyses_interpretation_SR.md`

---

## 9) Reproducibilnost

Svi parametri su u `config/config.yaml`. Logovi svakog pokretanja: `results/run_logs/`.

---

## 10) Ograničenja (za rad / odbranu)

- **Jedan donor** — nema bioloških replikata između donorova
- **DE na 3000 HVG** — neki markeri mogu biti isključeni iz DE (module score koristi pun gen prostor)
- **Anotacija** — Azimuth paneli i marker scoring mogu delimično kružno validirati iste labele
