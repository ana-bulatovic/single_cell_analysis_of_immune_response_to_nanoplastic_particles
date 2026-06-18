# Anotacija sa ref.Rds — uputstvo (prema zahtevu profesora)

## Šta profesor traži

| Uloga | Metoda | Gde u projektu |
|-------|--------|----------------|
| **Primarna anotacija** | Lokalni `ref.Rds` + Azimuth (R) | `cell_type_ref` → `cell_type_primary` |
| **Provera / validacija** | Klasični markeri iz literature | `cell_type_marker` (config.yaml) |
| **Azimuth (opciono)** | Ti biraš za šta ćeš | L1/L2/L3, score, marker paneli, UMAP |

---

## Fajlovi koje treba imati

```
data/raw/ref.Rds          # referentni PBMC atlas (Zenodo 4546839)
data/raw/idx.annoy        # indeks za Azimuth (isti Zenodo record)
data/processed/integrated_annotated.h5ad   # nakon prvog pipeline run-a
```

Preuzimanje: `python scripts/download_data.py` (ako već nisu tu).

---

## Redosled pokretanja

### 1. Prvi Python pipeline (integracija)

```bash
python scripts/run_pipeline.py
```

Ovo radi QC, integraciju, UMAP i **literaturne markere** (`cell_type_marker`).
Ako još nema `azimuth_annotations.csv`, DE/kompozicija privremeno koriste markere
(`primary_annotation: ref_rds` u config-u, ali CSV još ne postoji).

### 2. Priprema matrice za R

```bash
python scripts/prepare_azimuth_h5ad.py
```

Piše `data/processed/azimuth_mtx/` i `azimuth_obs.csv`.

### 3. Azimuth sa lokalnim ref.Rds (R)

```bash
Rscript scripts/azimuth_annotation.R
```

Skripta automatski koristi `data/raw/` (ref.Rds + idx.annoy) ako postoje.
Izlaz: **`results/tables/azimuth_annotations.csv`**

Kolone u CSV-u:
- `predicted.celltype.l1` — gruba PBMC klasa (CD4 T, CD8 T, B, Mono, …)
- `predicted.celltype.l2`, `predicted.celltype.l3` — finije podklase (opciono)
- `prediction.score.max` — pouzdanost mape (0–1)

### 4. Ponovni pipeline (primarna anotacija iz ref.Rds)

```bash
python scripts/run_pipeline.py
```

Sada STEP 3 učitava Azimuth CSV i postavlja:
- `cell_type_ref` — mapirani L1 tipovi (npr. `CD4_T`, `Monocyte_CD14`)
- `cell_type_primary` — **ovo ide u DE, kompoziciju i pathway analizu**

Ćelije sa score &lt; 0.5 dobijaju labelu `low_confidence` i **ne ulaze u DE**.

---

## Validacija klasičnim markerima

U `config/config.yaml` sekcija `markers:` (IL7R, LYZ, MS4A1, NKG7, …).

Pipeline računa `cell_type_marker` nezavisno od ref.Rds — **samo za proveru**.

Pogledaj:
- `results/figures/marker_dotplot.png` — da li markeri „sedu“ na očekivane tipove
- `results/figures/umap_celltypes_marker.png` vs `umap_celltypes_ref.png`
- `results/tables/annotation_agreement_metrics.csv` — % slaganja ref ↔ marker ↔ CoDi
- `results/tables/annotation_crosstab_ref_marker.csv`

Dodatno (bez punog pipeline-a):

```bash
python scripts/run_additional_analyses.py
```

---

## Azimuth — šta TI možeš da koristiš (opciono)

Azimuth nije obavezan za DE ako već imaš `cell_type_primary` iz ref.Rds.
Koristi ga dodatno po potrebi:

| Šta | Kako | Kada ima smisla |
|-----|------|-----------------|
| **L2/L3 finije tipove** | kolone `azimuth_l2`, `azimuth_l3` u h5ad / CSV | podanaliza T ćelija, monocita |
| **Marker paneli iz referentnog atlasa** | `python scripts/extract_azimuth_markers.py` | uporedi sa literature markerima |
| **UMAP po Azimuth L1** | ručno u Scanpy / dodatne figure | prezentacija, validacija |
| **Score filter** | `min_prediction_score` u config | isključi neizvesne ćelije |
| **Samo validacija** | `primary_annotation: markers` u config | ako želiš DE po markerima, ref samo za proveru |

Ekstrakcija Azimuth markera (opciono):

```bash
python scripts/extract_azimuth_markers.py
# → results/tables/azimuth_marker_panels_l1.yaml
# NE uključuj use_panels_for_pipeline: true osim ako namerno želiš DE po Azimuth markerima
```

---

## Kolone u `adata.obs` (rezime)

| Kolona | Izvor | Uloga |
|--------|-------|-------|
| `cell_type_primary` | ref.Rds ili markeri | **DE, kompozicija, pathway** |
| `cell_type_ref` | Azimuth + ref.Rds | primarna ref anotacija |
| `cell_type_marker` | config.yaml markeri | **validacija** |
| `cell_type_codi` | CoDi CSV | spoljna validacija |
| `azimuth_l1`, `azimuth_l2`, `azimuth_l3` | Azimuth CSV | opciono, finije labele |
| `azimuth_score` | Azimuth CSV | pouzdanost |

---

## Za izveštaj / prezentaciju (jedna rečenica)

> Tipove ćelija dodelili smo mapiranjem na lokalni PBMC referentni atlas (`ref.Rds`, Azimuth);
> dodelu smo proverili klasičnim markerima iz literature (IL7R, LYZ, MS4A1, …) i CoDi referencom;
> Azimuth L2/L3 i score koristili smo dodatno za [navedi šta si ti izabrala].

---

## Česta pitanja

**Zašto dva puta `run_pipeline.py`?**  
Prvi put pravi integrisani objekat; Azimuth u R-u radi na njemu; drugi put učitava CSV i menja primarnu anotaciju za DE.

**Šta ako nema ref.Rds?**  
R skripta pada na `pbmcref` (ugrađena referenca). Za zahtev profesora preuzmi Zenodo 4546839.

**Mono CD14 vs CD16?**  
Azimuth L1 ima samo „Mono“ → mapira se na `Monocyte_CD14`. Finiju podelu vidi u `azimuth_l2`/`l3`.

**Platelet?**  
Azimuth često ne daje trombocite; marker PPBP/PF4 u validacionom sloju ih i dalje hvata.
