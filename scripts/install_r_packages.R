# Seurat 5.x needs spatstat packages explicitly on some Windows installs
spatstat_deps <- c(
  "spatstat.data",
  "spatstat.geom",
  "spatstat.random",
  "spatstat.univar",
  "spatstat.utils",
  "spatstat.explore"
)
required_cran <- c(spatstat_deps, "Seurat", "SeuratDisk", "remotes", "readr", "dplyr", "Matrix")
user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib) && !dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
if (!nzchar(user_lib)) {
  user_lib <- file.path(Sys.getenv("HOME"), "R", paste0("library", R.version$major, ".", R.version$minor))
  if (!dir.exists(user_lib)) {
    dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  }
}
.libPaths(c(user_lib, .libPaths()))

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com", lib = user_lib)
  }
}

verify_load <- function(pkg) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  message(if (ok) paste0("[OK] ", pkg) else paste0("[MISSING] ", pkg))
  ok
}

message("Using R library path: ", user_lib)
for (pkg in required_cran) {
  install_if_missing(pkg)
}

# presto (Azimuth dependency) needs Rtools g++ on Windows — install from:
# https://cran.r-project.org/bin/windows/Rtools/
if (!requireNamespace("presto", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cran.rstudio.com", lib = user_lib)
  }
  message("Installing presto from GitHub (requires Rtools on Windows)...")
  tryCatch(
    remotes::install_github("immunogenomics/presto", lib = user_lib),
    error = function(e) {
      message(
        "presto install failed: ", conditionMessage(e),
        "\nInstall Rtools 4.5, restart terminal, then re-run this script."
      )
    }
  )
}

if (!requireNamespace("Azimuth", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cran.rstudio.com", lib = user_lib)
  }
  remotes::install_github("satijalab/azimuth", lib = user_lib)
}

message("\nPackage check:")
ok_seurat <- verify_load("Seurat")
ok_disk <- verify_load("SeuratDisk")
ok_presto <- verify_load("presto")
ok_az <- verify_load("Azimuth")

if (ok_seurat && ok_disk && ok_az) {
  message("\nReady. Run: Rscript scripts/azimuth_annotation.R")
} else if (ok_seurat && ok_disk && !ok_az) {
  message(
    "\nSeurat is ready but Azimuth is missing (usually presto/Rtools on Windows).",
    "\n1) Install Rtools: https://cran.r-project.org/bin/windows/Rtools/",
    "\n2) Re-run: Rscript scripts/install_r_packages.R"
  )
} else {
  message("\nSome packages still missing — see messages above.")
}
