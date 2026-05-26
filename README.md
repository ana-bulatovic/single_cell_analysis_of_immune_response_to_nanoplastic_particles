# Single-Cell Analysis of Immune Response to Nanoplastic Particles

Ovaj repozitorijum implementira kompletan scRNA-seq workflow za dataset sa Zenodo zapisa `15866724`:
- Sample 1: 40 nm PSNP
- Sample 2: 200 nm PSNP
- Sample 3: 40 + 200 nm (mix)
- Sample 4: control

Analiza je organizovana tako da direktno pokriva sve stavke zadatka: QC, integraciju, klasterizaciju, anotaciju, kompoziciju, DE i size-specific interpretaciju.

## 1) Podešavanje okruženja

**Opcija A (preporučeno na Windows-u, bez Conda):**

```bash
cd single_cell_analysis_of_immune_response_to_nanoplastic_particles
python -m pip install -r requirements.txt
```

**Opcija B (Conda + R za Azimuth):**

```bash
conda env create -f environment.yml
conda activate nanoplastic-scRNA
```

> Napomena: `harmonypy` često ne radi na Windows-u. Pipeline automatski koristi **Combat** integraciju (`integration_method: combat` u `config/config.yaml`).

### Procena trajanja (tipičan laptop)

| Korak | Trajanje |
|-------|----------|
| `pip install -r requirements.txt` | 5–15 min (jednom) |
| `python scripts/download_data.py` (~800 MB) | 10–40 min (zavisi od interneta) |
| `python scripts/run_pipeline.py` | 30–90 min |
| `python scripts/make_slides.py` | < 1 min |
| `Rscript scripts/azimuth_annotation.R` (opciono) | 15–45 min |

## 2) Preuzimanje podataka

```bash
python scripts/download_data.py
```

Skripta preuzima `.h5ad`, `.csv` i `.rds` fajlove iz Zenodo zapisa i čuva ih u `data/raw/`.

## 3) Pokretanje glavne analize

```bash
python scripts/run_pipeline.py
```

Pipeline radi sledeće:
1. **QC & preprocessing**
   - Filter: `min_genes=200`, `max_genes=7000`, `min_counts=500`, `pct_mt<=15`.
   - Justifikacija:
     - `<200` gena je često low-quality/debris.
     - veoma visok broj gena može ukazati na doublete.
     - visoki mt% ukazuje na oštećene/stresirane ćelije.
2. **Integracija & clustering**
   - Normalizacija + HVG selekcija (3000 gena).
   - Harmony batch correction po `sample_id`.
   - PCA, UMAP, Leiden klasteri.
3. **Cell type annotation**
   - Marker-based anotacija (T, B, NK, monocyte, DC, platelet).
   - Učitavanje CoDi CSV anotacija za poređenje.
4. **Composition analysis**
   - Proporcije cell tipova po uslovu.
5. **Differential expression & pathways**
   - Za svaki cell tip: svaki exposure vs control.
   - GO/KEGG/Reactome enrichment (gseapy/Enrichr).
6. **Size-specific effects**
   - `unique_40nm`, `unique_200nm`, `shared_40_200`, `shared_all_three`, `mix_only_emergent`.

## 4) Azimuth PBMC anotacija (R)

Ako želiš da dodaš PBMC referencu iz Azimuth modela, prvo proveri da li je R instaliran:

```bash
python scripts/check_rscript.py
```

Ako ti `Rscript` nije pronađen, proveri da li je instaliran i dostupna komanda:

```bash
where Rscript
```

Ako `where Rscript` ne vrati putanju, koristi punu lokaciju Rscript-a. Na primer:

```powershell
"C:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" scripts/install_r_packages.R
"C:\Program Files\R\R-4.6.0\bin\x64\Rscript.exe" scripts/azimuth_annotation.R
```

Ako `Rscript` radi normalno, onda pokreni:

```bash
Rscript scripts/install_r_packages.R
python scripts/download_data.py
Rscript scripts/azimuth_annotation.R
```

Ako se pojavi greška da `C:/Program Files/R/R-4.6.0/library` nije writable, koristi skriptu `scripts/install_r_packages.R`. Ona automatski instalira pakete u korisničku biblioteku.

Pipeline sada preuzima lokalni Azimuth PBMC reference model iz Zenodo zapisa `4546839` i čuva ga u `data/raw/ref.Rds`.

Izlaz: `results/tables/azimuth_annotations.csv`

Ako želiš da koristiš drugu lokalnu `.rds` referencu, zameni `reference_file` u `scripts/azimuth_annotation.R`.

## 5) Dodatne analize (3–5)

Implementirano u `scripts/run_pipeline.py`:
1. Cell-cycle score po uslovu.
2. IFN (interferon) signature score.
3. Antigen presentation signature po uslovu/cell tipu.
4. Pseudobulk matrica (condition × cell_type).
5. Sličnost marker anotacije i CoDi anotacije.

## 6) Generisanje PowerPoint-a

```bash
python scripts/make_slides.py
```

Izlaz: `deliverables/nanoplastic_scRNA_results.pptx`

## 7) Word dokumentacija (objašnjenje koda i genomike)

```bash
python scripts/generate_word_documentation.py
```

Izlaz: `deliverables/Project_Documentation.docx`

## 8) Predlog strukture za video (5–10 min)

Kratak scenario je u fajlu `deliverables/video_script.md`.

## 9) Glavni output fajlovi

- `data/processed/integrated_annotated.h5ad`
- `results/figures/umap_condition.png`
- `results/figures/umap_clusters.png`
- `results/figures/umap_celltypes_marker.png`
- `results/figures/composition_barplot.png`
- `results/tables/differential_expression_all.csv`
- `results/tables/pathway_enrichment_all.csv`
- `results/tables/size_specific_effects_summary.csv`

## Napomena o reproducibilnosti

Sve parametre možete menjati centralno u `config/config.yaml` bez izmene koda.
