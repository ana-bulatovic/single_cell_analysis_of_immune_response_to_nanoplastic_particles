"""
Reference-based cell-type annotation (ref.Rds atlas).

Methods (professor workflow):
  - Seurat/Azimuth  (R, ref.Rds)  -> cell_type_seurat
  - CoDi            (dataset CSV) -> cell_type_codi_norm
  - Tangram         (Python)      -> cell_type_tangram
  - cell2location   (signature mapping from ref.Rds) -> cell_type_cell2location

Literature markers (config.yaml) -> cell_type_marker (validation only).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from shutil import which
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy import sparse

PRIMARY_CELL_TYPE_EXCLUDE = frozenset({"low_confidence", "unmapped", "NA"})

REF_LABEL_TO_UNIFIED = {
    "CD4 T": "CD4_T",
    "CD4+ T": "CD4_T",
    "CD4+ T cell": "CD4_T",
    "CD8 T": "CD8_T_cytotoxic",
    "CD8+ T": "CD8_T_cytotoxic",
    "Cytotoxic T cell": "CD8_T_cytotoxic",
    "T cell": "T_cell",
    "B": "B_cell",
    "B cell": "B_cell",
    "NK": "NK_cell",
    "NK cell": "NK_cell",
    "Mono": "Monocyte_CD14",
    "CD14+ Mono": "Monocyte_CD14",
    "CD14+ monocyte": "Monocyte_CD14",
    "CD16+ Mono": "Monocyte_CD16",
    "CD16+ monocyte": "Monocyte_CD16",
    "DC": "DC",
    "Dendritic cell": "DC",
    "Platelet": "Platelet",
    "other T": "other_T",
    "other": "other",
}

CODI_TO_UNIFIED = {
    "CD4+ T cell": "CD4_T",
    "Cytotoxic T cell": "CD8_T_cytotoxic",
    "B cell": "B_cell",
    "CD14+ monocyte": "Monocyte_CD14",
    "CD16+ monocyte": "Monocyte_CD16",
    "NK cell": "NK_cell",
    "Dendritic cell": "DC",
    "DC": "DC",
    "Platelet": "Platelet",
}

AZIMUTH_L1_TO_UNIFIED = REF_LABEL_TO_UNIFIED


def normalize_ref_label(label: str) -> str:
    s = str(label).strip()
    if s in REF_LABEL_TO_UNIFIED:
        return REF_LABEL_TO_UNIFIED[s]
    return s.replace(" ", "_").replace("+", "")


def ensure_ref_h5ad(paths: Dict[str, Path], log) -> Optional[Path]:
    ref_h5ad = paths["processed"] / "ref_pbmc.h5ad"
    if ref_h5ad.exists():
        log.log(f"  Reference h5ad found: {ref_h5ad}")
        return ref_h5ad

    ref_rds = paths["raw"] / "ref.Rds"
    if not ref_rds.exists():
        log.log(f"  Warning: {ref_rds} missing — run python scripts/download_data.py")
        return None

    log.log("  Exporting ref.Rds to AnnData (R export + Python build)...")
    rscript = which("Rscript")
    if rscript is None:
        log.log("  Rscript not found in PATH — skip ref export (install R or disable run_tangram/run_cell2location)")
        return None

    try:
        r_status = subprocess.run(
            [rscript, "scripts/export_ref_h5ad.R"],
            cwd=paths["project_root"],
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        log.log("  Rscript not found — skip ref export")
        return None

    if r_status.returncode != 0:
        log.log(f"  R export failed: {r_status.stderr[-500:]}")
        return None

    sys.path.insert(0, str(paths["project_root"] / "scripts"))
    from build_ref_h5ad import build_ref_h5ad

    try:
        build_ref_h5ad(
            mtx_dir=paths["processed"] / "ref_mtx",
            obs_csv=paths["processed"] / "ref_obs.csv",
            out_h5ad=ref_h5ad,
        )
    except Exception as exc:
        log.log(f"  Reference h5ad build failed: {exc}")
        return None
    return ref_h5ad if ref_h5ad.exists() else None


def load_query_for_mapping(paths: Dict[str, Path], log) -> Optional[sc.AnnData]:
    pre_hvg = paths["processed"] / "integrated_pre_hvg.h5ad"
    if pre_hvg.exists():
        log.log(f"  Query for mapping: {pre_hvg}")
        return sc.read_h5ad(pre_hvg)

    integrated = paths["processed"] / "integrated_annotated.h5ad"
    if integrated.exists():
        log.log(f"  Warning: using HVG-only object for mapping ({integrated})")
        return sc.read_h5ad(integrated)

    log.log("  No query object for reference mapping.")
    return None


def _shared_genes(query: sc.AnnData, ref: sc.AnnData) -> List[str]:
    return sorted(set(query.var_names.astype(str)).intersection(ref.var_names.astype(str)))


def _prepare_query_ref(
    query: sc.AnnData, ref: sc.AnnData, log, max_ref_cells: int = 30000
) -> Tuple[sc.AnnData, sc.AnnData, List[str]]:
    genes = _shared_genes(query, ref)
    if len(genes) < 500:
        raise ValueError(f"Too few shared genes between query and ref: {len(genes)}")

    q = query[:, genes].copy()
    r = ref[:, genes].copy()

    if "counts" in q.layers:
        q.X = q.layers["counts"].copy()
    if r.X.max() > 30:
        pass
    else:
        sc.pp.normalize_total(r, target_sum=1e4)
        sc.pp.log1p(r)

    sc.pp.normalize_total(q, target_sum=1e4)
    sc.pp.log1p(q)

    if r.n_obs > max_ref_cells:
        rng = np.random.default_rng(42)
        idx = rng.choice(r.n_obs, size=max_ref_cells, replace=False)
        r = r[idx].copy()
        log.log(f"  Reference subsampled to {max_ref_cells:,} cells for mapping")

    log.log(f"  Shared genes for mapping: {len(genes):,}")
    return q, r, genes


def run_tangram_mapping(
    query: sc.AnnData,
    ref: sc.AnnData,
    paths: Dict[str, Path],
    log,
    label_key: str = "ref_celltype",
) -> pd.Series:
    try:
        import tangram as tg
    except ImportError:
        log.log("  Tangram not installed — skip (pip install tangram-sc)")
        return pd.Series(index=query.obs_names, dtype=object)

    log.log("  [Tangram] Mapping query cells to ref.Rds atlas...")
    q, r, genes = _prepare_query_ref(query, ref, log)
    adatas = [q, r]
    tg.pp_adatas(adatas, genes=genes)
    ad_map = tg.map_cells_to_reference(adatas[0], adatas[1], mode="clusters")
    tg.project_cell_annotations(ad_map, adatas[1], annotation=label_key)
    proj = ad_map.obs["ref_celltype"].astype(str)
    proj_norm = proj.map(normalize_ref_label).fillna("unmapped")
    out = paths["results"] / "tables" / "tangram_annotations.csv"
    pd.DataFrame({"cell_id": proj_norm.index, "cell_type_tangram": proj_norm.values}).to_csv(
        out, index=False
    )
    log.log(f"  Tangram mapped {len(proj_norm):,} cells -> {out.name}")
    return proj_norm.reindex(query.obs_names).fillna("unmapped")


def run_cell2location_mapping(
    query: sc.AnnData,
    ref: sc.AnnData,
    paths: Dict[str, Path],
    log,
    label_key: str = "ref_celltype",
) -> pd.Series:
    """
    Reference signature mapping from ref.Rds (cell2location-compatible for scRNA query).

    Builds per-cell-type mean expression signatures from the reference atlas and assigns
    each query cell to the best-matching type (cosine similarity in shared gene space).
    Full cell2location spatial deconvolution requires spatial slides; this uses the same
    reference atlas for label transfer on scRNA-seq query data.
    """
    log.log("  [cell2location-style] Signature mapping from ref.Rds...")
    q, r, genes = _prepare_query_ref(query, ref, log)

    labels = r.obs[label_key].astype(str)
    unique_types = sorted(labels.unique())
    ref_mat = r.X.toarray() if hasattr(r.X, "toarray") else np.asarray(r.X)
    q_mat = q.X.toarray() if hasattr(q.X, "toarray") else np.asarray(q.X)

    signatures = []
    type_names = []
    for ct in unique_types:
        mask = (labels == ct).values
        if mask.sum() < 5:
            continue
        signatures.append(ref_mat[mask].mean(axis=0))
        type_names.append(ct)

    if len(type_names) < 2:
        log.log("  cell2location-style mapping failed: too few reference types")
        return pd.Series(index=query.obs_names, dtype=object)

    sig = np.vstack(signatures)
    sig_norm = sig / (np.linalg.norm(sig, axis=1, keepdims=True) + 1e-8)
    q_norm = q_mat / (np.linalg.norm(q_mat, axis=1, keepdims=True) + 1e-8)
    sim = q_norm @ sig_norm.T
    best_idx = sim.argmax(axis=1)
    best_score = sim.max(axis=1)
    assigned = pd.Series(
        [normalize_ref_label(type_names[i]) for i in best_idx],
        index=q.obs_names,
    )
    assigned[best_score < 0.3] = "low_confidence"

    out = paths["results"] / "tables" / "cell2location_annotations.csv"
    pd.DataFrame(
        {
            "cell_id": assigned.index,
            "cell_type_cell2location": assigned.values,
            "mapping_score": best_score,
        }
    ).to_csv(out, index=False)
    log.log(f"  cell2location-style mapped {len(assigned):,} cells -> {out.name}")
    return assigned.reindex(query.obs_names).fillna("unmapped")


def load_seurat_annotations(
    adata: sc.AnnData, paths: Dict[str, Path], cfg: Dict, log
) -> bool:
    """Load Seurat/Azimuth labels from ref.Rds run (azimuth_annotations.csv)."""
    ann_cfg = cfg.get("annotation", {})
    if not ann_cfg.get("run_seurat", True):
        return False

    az_path = paths["results"] / "tables" / "azimuth_annotations.csv"
    if not az_path.exists():
        log.log("  Seurat/Azimuth CSV not found — run Rscript scripts/azimuth_annotation.R")
        return False

    az = pd.read_csv(az_path)
    az = az.assign(cell_id_clean=az["cell_id"].astype(str).str.replace(r"-\d+$", "", regex=True))
    if "cell_id_clean" not in adata.obs.columns:
        adata.obs["cell_id_clean"] = adata.obs_names.str.replace(r"-\d+$", "", regex=True)

    az_map = az.set_index("cell_id_clean")
    adata.obs["azimuth_l1"] = adata.obs["cell_id_clean"].map(az_map["predicted.celltype.l1"]).fillna("NA")
    adata.obs["azimuth_score"] = adata.obs["cell_id_clean"].map(az_map["prediction.score.max"]).astype(float)

    min_score = float(
        cfg.get("annotation", {}).get("min_prediction_score")
        or cfg.get("azimuth", {}).get("min_prediction_score", 0.5)
    )
    mapped = adata.obs["azimuth_l1"].map(AZIMUTH_L1_TO_UNIFIED).fillna("NA")
    seurat = mapped.copy()
    seurat[adata.obs["azimuth_score"].fillna(0) < min_score] = "low_confidence"
    seurat[adata.obs["azimuth_l1"] == "NA"] = "unmapped"
    adata.obs["cell_type_seurat"] = seurat.values
    adata.obs["cell_type_ref"] = adata.obs["cell_type_seurat"].values
    adata.obs["azimuth_l1_norm"] = mapped.values

    mapped_pct = (adata.obs["cell_type_seurat"] != "unmapped").mean() * 100
    log.log(f"  Seurat/Azimuth (ref.Rds): mapped {mapped_pct:.1f}% of cells")
    return True


def load_codi_annotations(adata: sc.AnnData, paths: Dict[str, Path], log) -> None:
    log.log("  [CoDi] Loading labels from *_CoDi_KLD.csv ...")
    codi_files = sorted(paths["raw"].glob("*_CoDi_KLD.csv"))
    codi_frames = []
    for f in codi_files:
        df = pd.read_csv(f)
        if "cell_id" in df.columns and "CoDi" in df.columns:
            keep = ["cell_id", "CoDi"]
            if "CoDi_confidence" in df.columns:
                keep.append("CoDi_confidence")
            df = df[keep].copy()
            codi_frames.append(df)

    if not codi_frames:
        adata.obs["cell_type_codi"] = "NA"
        adata.obs["cell_type_codi_norm"] = "NA"
        log.log("  Warning: no CoDi CSV files found")
        return

    codi_all = pd.concat(codi_frames, ignore_index=True).drop_duplicates(subset=["cell_id"])
    codi_all = codi_all.assign(cell_id_clean=codi_all["cell_id"].str.replace(r"-\d+$", "", regex=True))
    codi_all = codi_all.set_index("cell_id_clean")
    adata.obs["cell_id_clean"] = adata.obs_names.str.replace(r"-\d+$", "", regex=True)
    adata.obs["cell_type_codi"] = adata.obs["cell_id_clean"].map(codi_all["CoDi"]).fillna("NA")
    if "CoDi_confidence" in codi_all.columns:
        adata.obs["codi_confidence"] = adata.obs["cell_id_clean"].map(codi_all["CoDi_confidence"]).astype(float)
    adata.obs["cell_type_codi_norm"] = adata.obs["cell_type_codi"].map(CODI_TO_UNIFIED).fillna("NA")
    mapped = (adata.obs["cell_type_codi"] != "NA").mean() * 100
    log.log(f"  CoDi labels mapped to {mapped:.1f}% of cells")


def run_python_reference_methods(
    adata: sc.AnnData, paths: Dict[str, Path], cfg: Dict, log
) -> None:
    ann_cfg = cfg.get("annotation", {})
    if not ann_cfg.get("run_tangram", False) and not ann_cfg.get("run_cell2location", False):
        log.log("  Tangram/cell2location disabled in config — skip Python ref mapping")
        return

    ref_path = ensure_ref_h5ad(paths, log)
    if ref_path is None:
        log.log("  Skipping Tangram/cell2location (reference h5ad unavailable)")
        return

    query = load_query_for_mapping(paths, log)
    if query is None:
        return

    ref = sc.read_h5ad(ref_path)
    q_obs = query.obs.copy()
    q_obs.index = query.obs_names.astype(str)

    if ann_cfg.get("run_tangram", False):
        tangram = run_tangram_mapping(query, ref, paths, log)
        adata.obs["cell_type_tangram"] = adata.obs_names.map(tangram.to_dict()).fillna("unmapped")

    if ann_cfg.get("run_cell2location", False):
        c2l = run_cell2location_mapping(query, ref, paths, log)
        adata.obs["cell_type_cell2location"] = adata.obs_names.map(c2l.to_dict()).fillna("unmapped")


def assign_primary_cell_type(adata: sc.AnnData, cfg: Dict, log) -> None:
    ann_cfg = cfg.get("annotation", {})
    method = str(ann_cfg.get("primary_method", "seurat")).lower()
    col_map = {
        "seurat": "cell_type_seurat",
        "ref_rds": "cell_type_seurat",
        "codi": "cell_type_codi_norm",
        "tangram": "cell_type_tangram",
        "cell2location": "cell_type_cell2location",
        "markers": "cell_type_marker",
    }
    src = col_map.get(method, "cell_type_seurat")
    if src in adata.obs.columns and adata.obs[src].astype(str).ne("NA").any():
        adata.obs["cell_type_primary"] = adata.obs[src].astype(str)
        log.log(f"  Primary annotation for DE/composition: {method} ({src})")
        return

    if "cell_type_seurat" in adata.obs.columns:
        adata.obs["cell_type_primary"] = adata.obs["cell_type_seurat"].astype(str)
        log.log("  Primary fallback: seurat (ref.Rds)")
        return

    adata.obs["cell_type_primary"] = adata.obs["cell_type_marker"].astype(str)
    log.log("  Primary fallback: literature markers")


def validate_with_literature_markers(
    adata: sc.AnnData,
    marker_dict: Dict[str, List[str]],
    paths: Dict[str, Path],
    log,
) -> None:
    """Dotplots of literature markers grouped by each annotation method."""
    fig_dir = paths["results"] / "figures" / "annotation_validation"
    fig_dir.mkdir(parents=True, exist_ok=True)

    method_cols = [
        ("literature_markers", "cell_type_marker"),
        ("seurat_ref_rds", "cell_type_seurat"),
        ("codi", "cell_type_codi_norm"),
        ("tangram", "cell_type_tangram"),
        ("cell2location", "cell_type_cell2location"),
    ]
    panels = {
        ct: [g for g in genes if g in adata.var_names]
        for ct, genes in marker_dict.items()
    }
    panels = {k: v for k, v in panels.items() if v}
    if not panels:
        log.log("  Marker validation skipped (no genes in matrix)")
        return

    log.log("  Validating annotations with literature marker dotplots...")
    for name, col in method_cols:
        if col not in adata.obs.columns:
            continue
        labels = adata.obs[col].astype(str)
        if labels.eq("NA").all() or labels.eq("unmapped").all():
            continue
        try:
            sc.pl.dotplot(
                adata,
                var_names=panels,
                groupby=col,
                standard_scale="var",
                show=False,
                dendrogram=False,
            )
            plt.gcf().suptitle(f"Literature markers vs {name} labels", y=1.02, fontsize=11)
            out = fig_dir / f"marker_validation_{name}.png"
            plt.gcf().set_size_inches(12, 5)
            plt.savefig(out, dpi=300, bbox_inches="tight")
            plt.close()
            log.log(f"  Saved: {out}")
        except Exception as exc:
            log.log(f"  Dotplot failed for {name}: {exc}")


def export_method_agreement(
    adata: sc.AnnData, paths: Dict[str, Path], log
) -> pd.DataFrame:
    marker = adata.obs.get("cell_type_marker", pd.Series(index=adata.obs_names, dtype=str)).astype(str)
    rows = []

    methods = {
        "seurat": "cell_type_seurat",
        "codi": "cell_type_codi_norm",
        "tangram": "cell_type_tangram",
        "cell2location": "cell_type_cell2location",
    }
    for name, col in methods.items():
        if col not in adata.obs.columns:
            continue
        labels = adata.obs[col].astype(str)
        valid = ~labels.isin({"NA", "unmapped", "low_confidence"})
        if valid.any():
            rows.append(
                {
                    "metric": f"{name}_marker_agreement",
                    "value": float((labels[valid] == marker[valid]).mean()),
                    "description": f"Fraction where {name} label matches literature markers",
                }
            )
            rows.append(
                {
                    "metric": f"{name}_mapped_fraction",
                    "value": float(valid.mean()),
                    "description": f"Fraction of cells with valid {name} label",
                }
            )

    primary = adata.obs.get("cell_type_primary", marker).astype(str)
    rows.insert(
        0,
        {
            "metric": "primary_marker_agreement",
            "value": float((primary == marker).mean()),
            "description": "Primary annotation vs literature marker validation",
        },
    )

    metrics = pd.DataFrame(rows)
    out = paths["results"] / "tables" / "annotation_method_agreement.csv"
    metrics.to_csv(out, index=False)
    log.log(f"  Saved: {out.name} ({len(metrics)} metrics)")

    ctab_rows = []
    for name, col in methods.items():
        if col in adata.obs.columns:
            ct = pd.crosstab(marker, adata.obs[col].astype(str), margins=True)
            p = paths["results"] / "tables" / f"annotation_crosstab_marker_{name}.csv"
            ct.to_csv(p)
            ctab_rows.append(name)
    if ctab_rows:
        log.log(f"  Crosstabs saved for: {', '.join(ctab_rows)}")

    return metrics
