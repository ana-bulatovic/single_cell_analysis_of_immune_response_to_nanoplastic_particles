"""
Export minimal inputs for Azimuth in R (avoids SeuratDisk HDF5 / Seurat 5 issues).

Reads:  data/processed/integrated_annotated.h5ad
Writes: data/processed/azimuth_mtx/          (10x matrix for Read10X)
        data/processed/azimuth_obs.csv       (cell metadata)
"""

from pathlib import Path

import gzip
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmwrite


def _write_10x_mtx(mtx_dir: Path, adata: sc.AnnData) -> None:
    """Write genes x cells matrix in 10x MTX layout for Seurat Read10X."""
    mtx_dir.mkdir(parents=True, exist_ok=True)
    x = adata.X
    if not sparse.issparse(x):
        x = sparse.csr_matrix(x)
    x = sparse.csc_matrix(x.T)  # genes x cells

    tmp_mtx = mtx_dir / "matrix.mtx"
    mmwrite(tmp_mtx, x)
    with open(tmp_mtx, "rb") as src, gzip.open(mtx_dir / "matrix.mtx.gz", "wb") as dst:
        dst.write(src.read())
    tmp_mtx.unlink()

    barcodes = pd.Series(adata.obs_names.astype(str), name="barcode")
    barcodes.to_csv(mtx_dir / "barcodes.tsv.gz", index=False, header=False, compression="gzip")

    features = pd.DataFrame(
        {
            "gene_id": adata.var_names.astype(str),
            "gene_symbol": adata.var_names.astype(str),
            "feature_type": "Gene Expression",
        }
    )
    features.to_csv(mtx_dir / "features.tsv.gz", index=False, header=False, sep="\t", compression="gzip")


def main() -> None:
    src = Path("data/processed/integrated_annotated.h5ad")
    mtx_dir = Path("data/processed/azimuth_mtx")
    obs_csv = Path("data/processed/azimuth_obs.csv")

    if not src.exists():
        raise FileNotFoundError(f"Missing {src}. Run: python scripts/run_pipeline.py")

    print(f"Reading {src} ...")
    adata = sc.read_h5ad(src)

    if "counts" in adata.layers:
        x = adata.layers["counts"]
        print("Using layers['counts'] as expression matrix.")
    else:
        x = adata.X
        print("Using adata.X as expression matrix (no counts layer).")

    if sparse.issparse(x):
        x = sparse.csr_matrix(x, dtype=np.float32)
    else:
        x = np.asarray(x, dtype=np.float32)

    out = sc.AnnData(X=x, obs=adata.obs.copy(), var=adata.var.copy())
    out.obs_names_make_unique()
    out.var_names_make_unique()

    mtx_dir.mkdir(parents=True, exist_ok=True)
    _write_10x_mtx(mtx_dir, out)
    print(f"Saved 10x matrix: {mtx_dir} ({out.n_obs:,} cells x {out.n_vars:,} genes)")

    obs = out.obs.copy()
    for col in obs.columns:
        obs[col] = obs[col].astype(str)
    obs.insert(0, "cell_id", obs.index.astype(str))
    obs.to_csv(obs_csv, index=False)
    print(f"Saved metadata: {obs_csv}")


if __name__ == "__main__":
    main()
