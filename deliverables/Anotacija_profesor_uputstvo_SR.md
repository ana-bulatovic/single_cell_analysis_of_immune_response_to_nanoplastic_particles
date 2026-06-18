# Anotacija prema uputstvu profesora

## Zahtev

1. **ref.Rds** kao referentni anotirani dataset (Zenodo 4546839)
2. Anotacija alatima: **Seurat**, **CoDi**, **Tangram**, **cell2location**
3. **Provera** klasičnim marker genima iz literature (`config.yaml`)

---

## Metode u projektu

| Alat | Kolona u podacima | Kako se pokreće |
|------|-------------------|-----------------|
| **Seurat** (Azimuth + ref.Rds) | `cell_type_seurat` | `Rscript scripts/azimuth_annotation.R` |
| **CoDi** | `cell_type_codi_norm` | automatski iz `data/raw/*_CoDi_KLD.csv` |
| **Tangram** | `cell_type_tangram` | Python (ref.Rds → h5ad) |
| **cell2location** | `cell_type_cell2location` | Python signature mapping iz ref.Rds |
| **Literaturni markeri** | `cell_type_marker` | validacija (IL7R, LYZ, MS4A1, …) |

**Primarna anotacija za DE/kompoziciju:** `cell_type_primary` (default: `seurat`)

Podesi u `config/config.yaml`:
```yaml
annotation:
  primary_method: seurat   # seurat | codi | tangram | cell2location | markers
```

---

## Redosled pokretanja

```bash
# 1. Preuzmi podatke + ref.Rds
python scripts/download_data.py

# 2. Integracija (pravi integrated_pre_hvg.h5ad)
python scripts/run_pipeline.py

# 3. Seurat/Azimuth sa ref.Rds
python scripts/prepare_azimuth_h5ad.py
Rscript scripts/azimuth_annotation.R

# 4. Tangram + cell2location (ref.Rds atlas)
python scripts/run_reference_annotation.py

# 5. Pun pipeline sa svim labelama
python scripts/run_pipeline.py
```

Korak 4 može i bez R-a ako već imaš `azimuth_annotations.csv`:
```bash
python scripts/run_reference_annotation.py --skip-r
```

---

## Validacija markerima

Za svaku metodu generiše se dotplot literaturnih markera:

```
results/figures/annotation_validation/
  marker_validation_literature_markers.png
  marker_validation_seurat_ref_rds.png
  marker_validation_codi.png
  marker_validation_tangram.png
  marker_validation_cell2location.png
```

Tabele slaganja metoda vs markeri:
```
results/tables/annotation_method_agreement.csv
results/tables/annotation_crosstab_marker_*.csv
```

---

## Napomene

- **Seurat/Azimuth** je standardni način korišćenja ref.Rds u R-u (`RunAzimuth`).
- **CoDi** labele dolaze iz dataseta (prethodno pokrenut CoDi algoritam na istim ćelijama).
- **cell2location** u punoj verziji zahteva spatial transcriptomics; ovde se koristi **signature mapping** iz istog ref.Rds atlasa (referentni potpisi tipova ćelija).
- **Tangram** je opciono: `pip install tangram-sc`

---

## Za izveštaj (1 rečenica)

> Tipove ćelija dodelili smo mapiranjem na referentni PBMC atlas (ref.Rds) pomoću Seurat/Azimuth, Tangram i cell2location-style signature mappinga, uporedili sa CoDi referencom, i validirali klasičnim marker genima iz literature za svaki tip.
