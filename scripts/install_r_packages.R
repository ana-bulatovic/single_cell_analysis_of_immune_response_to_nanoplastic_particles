required_cran <- c("Seurat", "SeuratDisk", "remotes", "readr", "dplyr")
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

message("Using R library path: ", user_lib)
for (pkg in required_cran) {
  install_if_missing(pkg)
}

if (!requireNamespace("Azimuth", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cran.rstudio.com", lib = user_lib)
  }
  remotes::install_github("satijalab/azimuth", lib = user_lib)
}

message("R packages installed or already present. Now run: Rscript scripts/azimuth_annotation.R")
