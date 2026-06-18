# Export local ref.Rds (Zenodo 4546839) to 10x MTX + metadata for Python mapping tools.
# Output: data/processed/ref_mtx/  and  data/processed/ref_obs.csv

args_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(args_file) == 1L) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", args_file)))
  if (basename(script_dir) == "scripts") {
    setwd(dirname(script_dir))
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

ref_rds <- "data/raw/ref.Rds"
out_dir <- "data/processed/ref_mtx"
out_obs <- "data/processed/ref_obs.csv"
max_cells <- as.integer(Sys.getenv("REF_EXPORT_MAX_CELLS", unset = "30000"))

if (!file.exists(ref_rds)) {
  stop("Missing ", ref_rds, " — run: python scripts/download_data.py")
}

message("Loading ", ref_rds, " ...")
ref <- readRDS(ref_rds)
if (!inherits(ref, "Seurat")) {
  stop("ref.Rds is not a Seurat object.")
}

label_candidates <- c(
  "celltype.l1", "predicted.celltype.l1", "celltype.l2",
  "cell_type", "CellType", "annotation"
)
label_col <- label_candidates[label_candidates %in% colnames(ref@meta.data)][1]
if (is.na(label_col)) {
  stop("No cell-type column found in ref.Rds metadata.")
}
message("Using reference label column: ", label_col)

DefaultAssay(ref) <- DefaultAssay(ref)
counts <- GetAssayData(ref, assay = DefaultAssay(ref), layer = "counts")
if (is.null(counts)) {
  counts <- GetAssayData(ref, assay = DefaultAssay(ref), slot = "counts")
}

meta <- ref@meta.data
meta$cell_id <- rownames(meta)
meta$ref_celltype <- as.character(meta[[label_col]])

if (ncol(counts) > max_cells) {
  message("Subsampling reference to ", max_cells, " cells (stratified by ", label_col, ") ...")
  set.seed(42)
  idx <- unlist(lapply(split(seq_len(ncol(counts)), meta$ref_celltype), function(i) {
    if (length(i) <= max(1L, floor(max_cells * length(i) / ncol(counts)))) {
      return(i)
    }
    sample(i, size = max(1L, floor(max_cells * length(i) / ncol(counts))))
  }))
  idx <- unique(idx)
  counts <- counts[, idx, drop = FALSE]
  meta <- meta[idx, , drop = FALSE]
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Writing 10x matrix to ", out_dir, " (", ncol(counts), " cells x ", nrow(counts), " genes)")

writeMM(counts, file = gzfile(file.path(out_dir, "matrix.mtx.gz"), "w"))

writeLines(colnames(counts), con = gzfile(file.path(out_dir, "barcodes.tsv.gz"), "w"))
features <- data.frame(
  gene_id = rownames(counts),
  gene_symbol = rownames(counts),
  feature_type = "Gene Expression"
)
write.table(
  features,
  file = gzfile(file.path(out_dir, "features.tsv.gz"), "w"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

meta_out <- meta[, c("cell_id", "ref_celltype"), drop = FALSE]
write.csv(meta_out, out_obs, row.names = FALSE)
message("Saved metadata: ", out_obs)
message("REF EXPORT COMPLETE")
