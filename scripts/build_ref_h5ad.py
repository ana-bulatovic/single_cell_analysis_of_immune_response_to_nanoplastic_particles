"""Build AnnData reference object from ref.Rds export (10x MTX + metadata)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import scanpy as sc


def build_ref_h5ad(
    mtx_dir: Path | None = None,
    obs_csv: Path | None = None,
    out_h5ad: Path | None = None,
) -> sc.AnnData:
    mtx_dir = mtx_dir or Path("data/processed/ref_mtx")
    obs_csv = obs_csv or Path("data/processed/ref_obs.csv")
    out_h5ad = out_h5ad or Path("data/processed/ref_pbmc.h5ad")

    if not mtx_dir.exists() or not obs_csv.exists():
        raise FileNotFoundError(
            f"Missing reference export. Run: Rscript scripts/export_ref_h5ad.R\n"
            f"  Expected: {mtx_dir} and {obs_csv}"
        )

    adata = sc.read_10x_mtx(mtx_dir, gex_only=False, make_unique=True)
    meta = pd.read_csv(obs_csv)
    meta = meta.set_index("cell_id")
    common = adata.obs_names.intersection(meta.index)
    adata = adata[common].copy()
    adata.obs["ref_celltype"] = meta.loc[common, "ref_celltype"].astype(str).values
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write(out_h5ad)
    print(f"Saved reference AnnData: {out_h5ad} ({adata.n_obs:,} cells x {adata.n_vars:,} genes)")
    return adata


if __name__ == "__main__":
    build_ref_h5ad()
