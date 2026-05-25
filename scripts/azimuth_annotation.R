# =============================================================================
# Azimuth PBMC anotacija (R)
# =============================================================================
# Ova skripta mapira integrisane ćelije na Satija lab PBMC referencu
# koristeći Azimuth paket u Seurat-u.
#
# Preduslov: prvo pokrenuti Python pipeline (run_pipeline.py) da se
#             kreira data/processed/integrated_annotated.h5ad
#
# Izlaz: results/tables/azimuth_annotations.csv
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)  # konverzija h5ad <-> h5seurat
  library(Azimuth)     # referentna mapa tipova ćelija
  library(readr)
  library(dplyr)
})

processed_file <- "data/processed/integrated_annotated.h5ad"
output_csv <- "results/tables/azimuth_annotations.csv"

if (!file.exists(processed_file)) {
  stop("Nedostaje integrisani fajl. Prvo pokreni: python scripts/run_pipeline.py")
}

# Konvertuje AnnData (Python) u Seurat objekat (R)
obj <- LoadH5Seurat(Convert(processed_file, dest = "h5seurat", overwrite = TRUE))

reference_file <- "data/raw/ref.Rds"
if (!file.exists(reference_file)) {
  stop(
    "Nedostaje Azimuth referentni model. Pokreni `python scripts/download_data.py` da preuzmeš azimuth reference iz Zenodo record 4546839."
  )
}

# RunAzimuth: za svaku ćeliju predviđa tip na 3 nivoa (gruba → fina anotacija)
# koristi lokalni Azimuth PBMC reference model iz `data/raw/ref.Rds`
obj <- RunAzimuth(obj, reference = reference_file)

# Izvlači predikcije i score poverenja u CSV
az <- obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  select(
    cell_id,
    predicted.celltype.l1,  # nivo 1: npr. lymphocyte, myeloid
    predicted.celltype.l2,  # nivo 2: npr. CD4 T, CD14 monocyte
    predicted.celltype.l3,  # nivo 3: najfinija anotacija
    prediction.score.max    # koliko je model siguran (0–1)
  )

write_csv(az, output_csv)
message("Sačuvano: ", output_csv)
