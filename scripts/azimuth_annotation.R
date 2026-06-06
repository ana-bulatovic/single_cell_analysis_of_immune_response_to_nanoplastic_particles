# =============================================================================
# Azimuth PBMC anotacija (R)
# =============================================================================
# Mapira integrisane celije na Satija lab PBMC referencu (Azimuth + Seurat).
#
# Preduslov: python scripts/run_pipeline.py  -> integrated_annotated.h5ad
#            python scripts/download_data.py -> ref.Rds (optional; pbmcref used by default)
#
# Izlaz:  results/tables/azimuth_annotations.csv
#         results/run_logs/azimuth_run_YYYYMMDD_HHMMSS.txt
# =============================================================================

# Run from repo root; if invoked via Rscript scripts/..., chdir to parent of scripts/
args_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(args_file) == 1L) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", args_file)))
  if (basename(script_dir) == "scripts") {
    setwd(dirname(script_dir))
  }
}

log_dir <- "results/run_logs"
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_path <- file.path(log_dir, paste0("azimuth_run_", stamp, ".txt"))
log_con <- file(log_path, open = "wt", encoding = "UTF-8")
on.exit(close(log_con), add = TRUE)

log_msg <- function(...) {
  text <- paste(..., collapse = " ")
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", text)
  message(line)
  writeLines(line, log_con, useBytes = FALSE)
  flush(log_con)
}

log_msg("AZIMUTH ANNOTATION - START")
log_msg("Working directory:", getwd())
log_msg("Log file:", normalizePath(log_path, winslash = "/"))

processed_mtx <- "data/processed/azimuth_mtx"
processed_obs <- "data/processed/azimuth_obs.csv"
source_h5ad <- "data/processed/integrated_annotated.h5ad"
output_csv <- "results/tables/azimuth_annotations.csv"
# RunAzimuth expects a built-in reference name (e.g. "pbmcref") or a directory with ref.Rds + idx.annoy
reference <- Sys.getenv("AZIMUTH_REFERENCE", unset = "pbmcref")
reference_dir <- "data/raw"

if (!file.exists(source_h5ad)) {
  log_msg("ERROR: Missing", source_h5ad, "- run: python scripts/run_pipeline.py")
  stop("Missing integrated AnnData file.")
}
if (!dir.exists(processed_mtx) || !file.exists(processed_obs)) {
  log_msg("Preparing Azimuth 10x export via Python...")
  prep_status <- system2("python", c("scripts/prepare_azimuth_h5ad.py"))
  if (prep_status != 0 || !dir.exists(processed_mtx)) {
    log_msg("ERROR: Failed to create", processed_mtx)
    stop("Run manually: python scripts/prepare_azimuth_h5ad.py")
  }
}
if (!file.exists(reference_dir) || !file.exists(file.path(reference_dir, "ref.Rds"))) {
  log_msg("Note: local ref.Rds not found; using Azimuth reference:", reference)
}

log_msg("Loading R packages (Seurat, Azimuth)...")
suppressPackageStartupMessages({
  library(Seurat)
  library(Azimuth)
  library(readr)
  library(dplyr)
})
log_msg("Packages loaded.")

log_msg("Loading 10x matrix into Seurat object...")
log_msg("Input:", processed_mtx)
t0 <- Sys.time()
counts <- Read10X(processed_mtx, gene.column = 2)
obj <- CreateSeuratObject(counts = counts, project = "psnp")
meta <- read_csv(processed_obs, show_col_types = FALSE)
meta <- as.data.frame(meta)
rownames(meta) <- meta$cell_id
meta$cell_id <- NULL
common <- intersect(colnames(obj), rownames(meta))
obj <- subset(obj, cells = common)
obj <- AddMetaData(obj, meta[common, , drop = FALSE])
log_msg(
  "Load done in",
  round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
  "min;",
  "cells:", ncol(obj),
  "features:", nrow(obj)
)

log_msg("Running RunAzimuth with reference:", reference)
log_msg("(expect ~15-45 min on a typical laptop for ~33k cells)")
t1 <- Sys.time()
obj <- RunAzimuth(obj, reference = reference)
log_msg(
  "RunAzimuth done in",
  round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1),
  "min"
)

log_msg("Extracting predictions to CSV:", output_csv)
score_col <- if ("prediction.score.max" %in% colnames(obj@meta.data)) {
  "prediction.score.max"
} else {
  "mapping.score"
}

az <- obj@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  select(
    cell_id,
    predicted.celltype.l1,
    predicted.celltype.l2,
    predicted.celltype.l3,
    all_of(score_col)
  ) %>%
  rename(prediction.score.max = all_of(score_col))

write_csv(az, output_csv)
log_msg("Saved", nrow(az), "cell rows to", output_csv)
log_msg(
  "Mean prediction score:",
  round(mean(az$prediction.score.max, na.rm = TRUE), 3)
)
log_msg(
  "AZIMUTH ANNOTATION - COMPLETE (total",
  round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
  "min)"
)
