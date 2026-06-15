"""
Extract marker genes per cell type using Azimuth reference labels.

Uses FULL gene expression from raw .h5ad files (post-QC), not the 3000-gene HVG
matrix in integrated_annotated.h5ad — canonical markers like CD3D are often
missing from HVG and produce nonsense panels.

Prerequisites (in order):
  1. python scripts/run_pipeline.py
  2. python scripts/prepare_azimuth_h5ad.py   # optional before Azimuth
  3. Rscript scripts/azimuth_annotation.R     # refresh labels if needed
  4. python scripts/extract_azimuth_markers.py

Usage:
  python scripts/extract_azimuth_markers.py
  python scripts/extract_azimuth_markers.py --level l2 --top-n 10
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import yaml
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_pipeline import PipelineLogger, load_config, setup_paths

LEVEL_COLUMNS = {
    "l1": "predicted.celltype.l1",
    "l2": "predicted.celltype.l2",
    "l3": "predicted.celltype.l3",
}


def _sanitize_label(label: str) -> str:
    key = re.sub(r"[^\w]+", "_", str(label).strip())
    key = re.sub(r"_+", "_", key).strip("_")
    return key or "unknown"


def _check_azimuth_vs_integrated(
    azimuth_path: Path, integrated_path: Path, log: PipelineLogger
) -> None:
    if not integrated_path.exists():
        return
    az_mtime = azimuth_path.stat().st_mtime
    int_mtime = integrated_path.stat().st_mtime
    if az_mtime < int_mtime:
        log.log(
            "  Note: azimuth_annotations.csv is OLDER than integrated_annotated.h5ad. "
            "If pipeline was re-run, refresh Azimuth first:"
        )
        log.log("    python scripts/prepare_azimuth_h5ad.py")
        log.log("    Rscript scripts/azimuth_annotation.R")
    else:
        log.log("  Azimuth CSV timestamp is newer or equal to integrated object — OK.")


def load_expression_for_integrated_cells(
    cfg: Dict,
    paths: Dict[str, Path],
    integrated_path: Path,
    log: PipelineLogger,
) -> sc.AnnData:
    """
    Load full-gene counts only for cells present in integrated_annotated.h5ad.
    Avoids re-loading all QC-failed cells and lowers peak RAM vs full raw merge.
    """
    if not integrated_path.exists():
        raise FileNotFoundError(f"Missing {integrated_path}. Run: python scripts/run_pipeline.py")

    log.log(f"  Reading cell barcodes from {integrated_path} ...")
    ref = sc.read_h5ad(integrated_path, backed="r")
    sample_by_cell = ref.obs["sample_id"].astype(str).to_dict()
    cells_by_sample: Dict[str, List[str]] = {}
    for cell_id, sample_id in sample_by_cell.items():
        cells_by_sample.setdefault(sample_id, []).append(cell_id)
    ref.file.close()

    adata = None
    for sample_id in cfg["samples"]:
        raw_path = paths["raw"] / f"{sample_id}.h5ad"
        if not raw_path.exists():
            raise FileNotFoundError(f"Missing {raw_path}")
        want = cells_by_sample.get(sample_id, [])
        if not want:
            continue
        log.log(f"  Loading {len(want):,} cells from {sample_id}.h5ad ...")
        raw = sc.read_h5ad(raw_path, backed="r")
        present = [c for c in want if c in raw.obs_names]
        chunk = raw[present].to_memory()
        raw.file.close()
        chunk.obs["sample_id"] = sample_id
        chunk.obs["condition"] = cfg["samples"][sample_id]["condition"]
        if sparse.issparse(chunk.X):
            chunk.layers["counts"] = sparse.csr_matrix(chunk.X, dtype=np.float32)
        else:
            chunk.layers["counts"] = np.asarray(chunk.X, dtype=np.float32)
        adata = chunk if adata is None else sc.concat([adata, chunk], join="inner")
        del chunk

    if adata is None:
        raise RuntimeError("No cells loaded from raw files")
    adata.obs_names_make_unique()
    log.log(f"  Full-gene expression: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    return adata


def _subset_hvg_for_markers(adata: sc.AnnData, n_top: int, log: PipelineLogger) -> sc.AnnData:
    """Pick HVGs from full transcriptome (includes canonical markers missing from pipeline HVG)."""
    work = adata.copy()
    counts = work.layers["counts"]
    if sparse.issparse(counts):
        work.X = counts.copy().astype(np.float32)
    else:
        work.X = np.asarray(counts, dtype=np.float32)
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    sc.pp.highly_variable_genes(work, n_top_genes=n_top, flavor="seurat", subset=True)
    work.layers["counts"] = adata[:, work.var_names].layers["counts"].copy()
    log.log(f"  HVG subset for marker search: {work.n_vars:,} genes")
    return work


def attach_azimuth_labels(
    adata: sc.AnnData,
    azimuth_path: Path,
    level: str,
    min_score: float,
    log: PipelineLogger,
) -> sc.AnnData:
    az = pd.read_csv(azimuth_path)
    label_col = LEVEL_COLUMNS[level]
    if label_col not in az.columns:
        raise ValueError(f"Column {label_col} not found in {azimuth_path}")

    az = az.set_index("cell_id")
    common = adata.obs_names.intersection(az.index)
    if len(common) < adata.n_obs * 0.9:
        log.log(
            f"  Warning: only {len(common):,}/{adata.n_obs:,} cells match Azimuth barcodes. "
            "Re-run Azimuth after the latest pipeline."
        )

    adata = adata[adata.obs_names.isin(common)].copy()
    adata.obs["azimuth_label"] = adata.obs_names.map(az[label_col])
    adata.obs["azimuth_score"] = pd.to_numeric(
        adata.obs_names.map(az["prediction.score.max"]), errors="coerce"
    )
    adata.obs["azimuth_label"] = adata.obs["azimuth_label"].fillna("unassigned")

    before = adata.n_obs
    adata = adata[adata.obs["azimuth_score"].fillna(0) >= min_score].copy()
    adata = adata[adata.obs["azimuth_label"] != "unassigned"].copy()
    log.log(
        f"  Azimuth {level.upper()}, score>={min_score}: "
        f"{adata.n_obs:,} / {before:,} cells kept"
    )
    for label, n in adata.obs["azimuth_label"].value_counts().items():
        log.log(f"    {label}: {n:,}")
    return adata


def _adata_log_normalized(adata: sc.AnnData) -> sc.AnnData:
    out = adata.copy()
    counts = out.layers["counts"]
    if sparse.issparse(counts):
        out.X = counts.copy().astype(np.float32)
    else:
        out.X = np.asarray(counts, dtype=np.float32)
    sc.pp.normalize_total(out, target_sum=1e4)
    sc.pp.log1p(out)
    return out


def _detection_fracs_for_genes(
    counts: sparse.spmatrix | np.ndarray,
    row_mask: np.ndarray,
    gene_indices: List[int],
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (frac_in_group, frac_out_group) for selected gene columns."""
    if sparse.issparse(counts):
        counts = counts.tocsr()
    in_idx = np.where(row_mask)[0]
    out_idx = np.where(~row_mask)[0]
    in_frac = np.zeros(len(gene_indices), dtype=np.float64)
    out_frac = np.zeros(len(gene_indices), dtype=np.float64)
    for i, j in enumerate(gene_indices):
        if sparse.issparse(counts):
            in_frac[i] = float((counts[in_idx, j] > 0).mean())
            out_frac[i] = float((counts[out_idx, j] > 0).mean())
        else:
            in_frac[i] = float((counts[in_idx, j] > 0).mean())
            out_frac[i] = float((counts[out_idx, j] > 0).mean())
    return in_frac, out_frac


def rank_markers_one_vs_rest(
    adata: sc.AnnData,
    min_cells: int,
    log: PipelineLogger,
    top_genes_per_type: int = 800,
) -> Tuple[pd.DataFrame, sc.AnnData]:
    valid_types = adata.obs["azimuth_label"].value_counts()
    valid_types = valid_types[valid_types >= min_cells].index.tolist()
    if len(valid_types) < 2:
        raise ValueError(f"Need >=2 Azimuth groups with >={min_cells} cells")

    adata_sub = adata[adata.obs["azimuth_label"].isin(valid_types)].copy()
    adata_sub.obs["azimuth_label"] = adata_sub.obs["azimuth_label"].astype("category")
    adata_de = _adata_log_normalized(adata_sub)

    log.log(
        f"  Wilcoxon one-vs-rest on {adata_de.n_vars:,} genes, "
        f"{len(valid_types)} Azimuth types..."
    )
    sc.tl.rank_genes_groups(
        adata_de,
        groupby="azimuth_label",
        method="wilcoxon",
        use_raw=False,
    )

    counts = adata_sub.layers["counts"]
    labels = adata_sub.obs["azimuth_label"].to_numpy()
    frames = []

    for group in valid_types:
        df = sc.get.rank_genes_groups_df(adata_de, group=group).head(top_genes_per_type)
        df["azimuth_cell_type"] = group
        df["azimuth_cell_type_key"] = _sanitize_label(group)

        genes = df["names"].astype(str).tolist()
        gene_idx = [adata_sub.var_names.get_loc(g) for g in genes if g in adata_sub.var_names]
        valid_genes = [g for g in genes if g in adata_sub.var_names]
        df = df[df["names"].isin(valid_genes)].copy()

        in_mask = labels == group
        in_frac, out_frac = _detection_fracs_for_genes(counts, in_mask, gene_idx)
        df["pct_in_group"] = in_frac
        df["pct_out_group"] = out_frac
        df["specificity"] = df["pct_in_group"] - df["pct_out_group"]
        frames.append(df)

    markers = pd.concat(frames, ignore_index=True)
    markers = markers.rename(
        columns={
            "names": "gene",
            "logfoldchanges": "logfc",
            "pvals": "pval",
            "pvals_adj": "padj",
        }
    )
    log.log(f"  Ranked {len(markers):,} gene x cell-type rows (top {top_genes_per_type} per type)")
    return markers, adata_sub


def build_marker_panels(
    markers: pd.DataFrame,
    top_n: int,
    padj_thr: float,
    min_logfc: float,
    min_pct_in_group: float,
    max_pct_out_group: float,
) -> pd.DataFrame:
    sig = markers[
        markers["padj"].notna()
        & markers["logfc"].notna()
        & (markers["padj"] < padj_thr)
        & (markers["logfc"] > min_logfc)
        & (markers["pct_in_group"] >= min_pct_in_group)
        & (markers["pct_out_group"] <= max_pct_out_group)
    ].copy()

    panels = []
    for ctype, sub in sig.groupby("azimuth_cell_type", observed=False):
        top = sub.sort_values(
            ["specificity", "logfc", "padj"],
            ascending=[False, False, True],
        ).head(top_n)
        for rank, row in enumerate(top.itertuples(index=False), start=1):
            panels.append(
                {
                    "azimuth_cell_type": ctype,
                    "azimuth_cell_type_key": row.azimuth_cell_type_key,
                    "rank": rank,
                    "gene": row.gene,
                    "logfc": row.logfc,
                    "padj": row.padj,
                    "pct_in_group": row.pct_in_group,
                    "pct_out_group": row.pct_out_group,
                    "specificity": row.specificity,
                }
            )
    return pd.DataFrame(panels)


def panels_to_yaml_dict(panels: pd.DataFrame) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for key, sub in panels.groupby("azimuth_cell_type_key", observed=False):
        out[str(key)] = sub.sort_values("rank")["gene"].astype(str).tolist()
    return out


def compare_with_config_markers(panels: pd.DataFrame, cfg: Dict) -> pd.DataFrame:
    config_markers = cfg.get("markers", {})
    rows = []
    for az_key, sub in panels.groupby("azimuth_cell_type_key", observed=False):
        az_genes = set(sub["gene"].astype(str))
        best_cfg, best_overlap = "", 0
        for cfg_key, genes in config_markers.items():
            overlap = len(az_genes & set(genes))
            if overlap > best_overlap:
                best_overlap = overlap
                best_cfg = cfg_key
        rows.append(
            {
                "azimuth_cell_type_key": az_key,
                "azimuth_cell_type": sub["azimuth_cell_type"].iloc[0],
                "n_azimuth_markers": len(az_genes),
                "best_matching_config_type": best_cfg,
                "overlap_with_config": best_overlap,
            }
        )
    return pd.DataFrame(rows)


def plot_marker_dotplot(
    adata: sc.AnnData,
    panels: pd.DataFrame,
    out_path: Path,
    max_genes_per_type: int = 5,
) -> None:
    if panels.empty:
        return

    gene_order: List[str] = []
    for _, sub in panels.groupby("azimuth_cell_type", sort=False):
        for g in sub.sort_values("rank")["gene"].head(max_genes_per_type):
            if g not in gene_order and g in adata.var_names:
                gene_order.append(g)
    if not gene_order:
        return

    counts = adata.layers["counts"]
    rows = []
    for ctype in panels["azimuth_cell_type"].drop_duplicates():
        mask = (adata.obs["azimuth_label"] == ctype).to_numpy()
        idx = np.where(mask)[0]
        for gene in gene_order:
            if gene not in adata.var_names:
                continue
            j = adata.var_names.get_loc(gene)
            if sparse.issparse(counts):
                vals = np.asarray(counts[idx, j].todense()).ravel()
            else:
                vals = counts[idx, j]
            rows.append(
                {
                    "azimuth_cell_type": ctype,
                    "gene": gene,
                    "fraction_expressing": float((vals > 0).mean()),
                    "mean_expression": float(vals.mean()),
                }
            )

    plot_df = pd.DataFrame(rows)
    plt.figure(figsize=(max(8, 0.45 * len(gene_order) + 3), max(4, 0.35 * plot_df["azimuth_cell_type"].nunique() + 2)))
    sns.scatterplot(
        data=plot_df,
        x="gene",
        y="azimuth_cell_type",
        size="fraction_expressing",
        hue="mean_expression",
        palette="viridis",
        sizes=(30, 280),
        edgecolor="0.3",
        linewidth=0.3,
        legend=False,
    )
    plt.xticks(rotation=60, ha="right")
    plt.xlabel("Azimuth-derived marker gene")
    plt.ylabel("Azimuth cell type")
    plt.title("Top Azimuth marker genes (full transcriptome)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close()


def parse_args(cfg: Dict) -> argparse.Namespace:
    ext = cfg.get("azimuth", {}).get("marker_extraction", {})
    parser = argparse.ArgumentParser(description="Extract marker genes from Azimuth labels")
    parser.add_argument("--level", choices=["l1", "l2", "l3"], default=ext.get("level", "l1"))
    parser.add_argument("--min-score", type=float, default=ext.get("min_score", 0.5))
    parser.add_argument("--min-cells", type=int, default=ext.get("min_cells", 20))
    parser.add_argument("--top-n", type=int, default=ext.get("top_n", 10))
    parser.add_argument("--padj", type=float, default=ext.get("padj", 0.05))
    parser.add_argument("--min-logfc", type=float, default=ext.get("min_logfc", 0.5))
    parser.add_argument("--min-pct", type=float, default=ext.get("min_pct_in_group", 0.25))
    parser.add_argument("--max-pct-out", type=float, default=ext.get("max_pct_out_group", 0.2))
    parser.add_argument(
        "--azimuth-csv",
        type=Path,
        default=Path("results/tables/azimuth_annotations.csv"),
    )
    parser.add_argument(
        "--integrated-h5ad",
        type=Path,
        default=Path("data/processed/integrated_annotated.h5ad"),
        help="Used only to check whether Azimuth CSV may be stale",
    )
    return parser.parse_args()


def main() -> None:
    cfg = load_config()
    args = parse_args(cfg)
    paths = setup_paths(cfg)
    log = PipelineLogger(paths["run_logs"])

    try:
        log.section("AZIMUTH MARKER EXTRACTION (full gene space)")
        if not args.azimuth_csv.exists():
            raise FileNotFoundError(
                f"Missing {args.azimuth_csv}\n"
                "Run: python scripts/prepare_azimuth_h5ad.py\n"
                "     Rscript scripts/azimuth_annotation.R"
            )

        _check_azimuth_vs_integrated(args.azimuth_csv, args.integrated_h5ad, log)

        adata = load_expression_for_integrated_cells(
            cfg, paths, args.integrated_h5ad, log
        )
        adata = attach_azimuth_labels(
            adata, args.azimuth_csv, args.level, args.min_score, log
        )
        n_hvg = int(cfg.get("azimuth", {}).get("marker_extraction", {}).get("n_hvg", 4000))
        adata = _subset_hvg_for_markers(adata, n_hvg, log)

        markers, adata_labeled = rank_markers_one_vs_rest(adata, args.min_cells, log)
        level_tag = args.level.lower()
        tables = paths["results"] / "tables"

        all_csv = tables / f"azimuth_marker_genes_{level_tag}.csv"
        markers.to_csv(all_csv, index=False)
        log.log(f"  Saved: {all_csv}")

        panels = build_marker_panels(
            markers,
            top_n=args.top_n,
            padj_thr=args.padj,
            min_logfc=args.min_logfc,
            min_pct_in_group=args.min_pct,
            max_pct_out_group=args.max_pct_out,
        )
        panel_csv = tables / f"azimuth_marker_panels_{level_tag}.csv"
        panels.to_csv(panel_csv, index=False)
        log.log(f"  Saved: {panel_csv} ({len(panels)} panel genes)")

        panel_yaml = tables / f"azimuth_marker_panels_{level_tag}.yaml"
        with open(panel_yaml, "w", encoding="utf-8") as f:
            yaml.safe_dump(
                panels_to_yaml_dict(panels),
                f,
                sort_keys=False,
                allow_unicode=True,
                default_flow_style=False,
            )
        log.log(f"  Saved: {panel_yaml}")

        cmp_csv = tables / f"azimuth_vs_config_markers_{level_tag}.csv"
        compare_with_config_markers(panels, cfg).to_csv(cmp_csv, index=False)
        log.log(f"  Saved: {cmp_csv}")

        fig_path = paths["results"] / "figures" / f"azimuth_marker_dotplot_{level_tag}.png"
        plot_marker_dotplot(adata_labeled, panels, fig_path)
        log.log(f"  Saved: {fig_path}")

        log.section("TOP MARKERS (preview)")
        for ctype, sub in panels.groupby("azimuth_cell_type", observed=False):
            genes = ", ".join(sub.sort_values("rank")["gene"].astype(str))
            log.log(f"  {ctype}: {genes}")

        log.section("DONE")
    finally:
        log.close()


if __name__ == "__main__":
    main()
