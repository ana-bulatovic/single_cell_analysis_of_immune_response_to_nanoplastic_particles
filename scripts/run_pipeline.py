"""
Main scRNA-seq pipeline: immune response to nanoplastic particles.

Analysis flow:
  1. Load 4 samples (.h5ad) + QC filtering
  2. Merge, normalize, HVG, integration (Combat), UMAP, Leiden
  3. Cell-type annotation (markers) + CoDi reference labels
  4. Composition, differential expression, pathway enrichment
  5. Size-specific effects (40 nm vs 200 nm vs mix)
  6. Additional analyses (cell cycle, IFN, antigen presentation, pseudobulk)

Console and per-run log files are written in English for presentation notes.
"""

import gc
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, TextIO

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import yaml
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from annotation_methods import (
    PRIMARY_CELL_TYPE_EXCLUDE,
    assign_primary_cell_type,
    export_method_agreement,
    load_codi_annotations,
    load_seurat_annotations,
    run_python_reference_methods,
    validate_with_literature_markers,
)


class PipelineLogger:
    """Mirror messages to stdout and to a timestamped run log file."""

    def __init__(self, log_dir: Path):
        log_dir.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.log_path = log_dir / f"pipeline_run_{stamp}.txt"
        self._file: TextIO = open(self.log_path, "w", encoding="utf-8")
        self._write_header()

    def _write_header(self) -> None:
        self.log("")
        self.log("=" * 72)
        self.log("  SINGLE-CELL PIPELINE - RUN LOG")
        self.log(f"  Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self.log(f"  Log file: {self.log_path.resolve()}")
        self.log("=" * 72)

    def log(self, message: str = "") -> None:
        self._file.write(message + "\n")
        self._file.flush()
        encoding = getattr(sys.stdout, "encoding", None) or "utf-8"
        try:
            print(message)
        except UnicodeEncodeError:
            print(message.encode(encoding, errors="replace").decode(encoding))

    def section(self, title: str) -> None:
        self.log("")
        self.log("-" * 72)
        self.log(f"  {title}")
        self.log("-" * 72)

    def close(self) -> None:
        self.log("")
        self.log("=" * 72)
        self.log(f"  Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self.log("=" * 72)
        self._file.close()


CONDITION_COLORS = {
    "control": "#4C72B0",
    "PSNP_40nm": "#DD8452",
    "PSNP_200nm": "#55A868",
    "PSNP_mix_40_200": "#C44E52",
}

CONDITION_LABELS = {
    "control": "Control",
    "PSNP_40nm": "PSNP 40 nm",
    "PSNP_200nm": "PSNP 200 nm",
    "PSNP_mix_40_200": "PSNP mix",
}

EXPOSURE_COMPARISONS = [
    ("PSNP_40nm_vs_control", "40 nm"),
    ("PSNP_200nm_vs_control", "200 nm"),
    ("PSNP_mix_40_200_vs_control", "mix"),
]

CONDITION_ORDER = ["control", "PSNP_40nm", "PSNP_200nm", "PSNP_mix_40_200"]

GENE_SET_PREFIX = {
    "GO_Biological_Process_2023": "GO",
    "KEGG_2021_Human": "KEGG",
    "Reactome_2022": "Reactome",
}

FIGURE_DESCRIPTIONS = {
    "umap_condition.png": (
        "UMAP of all integrated cells colored by experimental condition "
        "(40 nm PSNP, 200 nm PSNP, mixed PSNP, or untreated control). "
        "Shows whether samples mix after batch correction and if exposure "
        "drives visible transcriptional shifts."
    ),
    "umap_clusters.png": (
        "UMAP colored by Leiden clusters (unsupervised groups). "
        "Each cluster is a group of cells with similar expression profiles; "
        "clusters are later linked to cell types via marker genes."
    ),
    "umap_celltypes_marker.png": (
        "UMAP colored by marker-based PBMC cell types (T, B, NK, monocytes, etc.). "
        "Assignment uses canonical marker gene sets from config.yaml."
    ),
    "umap_split_by_condition.png": (
        "Same UMAP embedding split into four panels (one per condition). "
        "Highlights each exposure against a gray background of all other cells — "
        "useful to compare whether any condition occupies distinct regions."
    ),
    "umap_sample_integration.png": (
        "UMAP colored by sample_id after Combat batch correction. "
        "Good integration means all four samples intermix rather than forming separate islands."
    ),
    "umap_module_scores.png": (
        "UMAP colored by cell-cycle (S, G2M) and interferon module scores. "
        "Shows spatial distribution of proliferation and innate immune activation signals."
    ),
    "umap_codi_celltypes.png": (
        "UMAP colored by external CoDi reference labels (when available). "
        "Independent validation of marker-based annotation."
    ),
    "marker_dotplot.png": (
        "Dot plot of canonical marker genes per assigned cell type. "
        "Validates that marker-based labels match expected PBMC expression patterns."
    ),
    "composition_barplot.png": (
        "Stacked bar chart of cell-type fractions per condition. "
        "Answers whether nanoplastic exposure changes immune cell proportions "
        "(e.g. more monocytes, fewer T cells) relative to control."
    ),
    "annotation_confusion_marker_codi.png": (
        "Contingency matrix (cell counts): literature marker annotation vs CoDi "
        "reference labels. Shows how each marker-based type maps to CoDi and "
        "where the methods disagree."
    ),
    "annotation_confusion_marker_codi_normalized.png": (
        "Row-normalized CoDi vs marker contingency matrix (% of each marker type). "
        "Highlights the dominant CoDi label assigned to each marker-based cell type."
    ),
    "annotation_confusion_ref_codi.png": (
        "Contingency matrix (cell counts): Seurat/ref.Rds primary annotation vs CoDi "
        "reference labels. Shows how ref.Rds-assigned types map to CoDi."
    ),
    "annotation_confusion_ref_codi_normalized.png": (
        "Row-normalized ref.Rds vs CoDi contingency matrix (% of each ref.Rds type). "
        "Primary annotation cross-validation against external CoDi labels."
    ),
}


def load_config(path: str = "config/config.yaml") -> Dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_paths(cfg: Dict) -> Dict[str, Path]:
    paths = {
        "project_root": Path(__file__).resolve().parent.parent,
        "raw": Path(cfg["data"]["raw_dir"]),
        "processed": Path(cfg["data"]["processed_dir"]),
        "results": Path(cfg["data"]["results_dir"]),
        "run_logs": Path(cfg["data"].get("run_logs_dir", "results/run_logs")),
    }
    for key in ("raw", "processed", "results", "run_logs"):
        paths[key].mkdir(parents=True, exist_ok=True)
    (paths["results"] / "figures").mkdir(parents=True, exist_ok=True)
    (paths["results"] / "tables").mkdir(parents=True, exist_ok=True)
    return paths


def _compact_adata_for_merge(adata: sc.AnnData) -> sc.AnnData:
    """Drop unused AnnData slots and use CSR float32 to lower concat peak RAM."""
    adata.raw = None
    adata.layers.clear()
    adata.obsm.clear()
    adata.varm.clear()
    adata.uns.clear()
    if sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X, dtype=np.float32)
    else:
        adata.X = np.asarray(adata.X, dtype=np.float32)
    return adata


def describe_raw_data(cfg: Dict, paths: Dict[str, Path], log: PipelineLogger) -> None:
    """Print minimal but necessary overview of input .h5ad files."""
    log.section("STEP 0 - RAW DATA OVERVIEW")
    log.log(
        "Dataset: PBMC from one donor exposed to polystyrene nanoplastic particles (PSNP) "
        "at 40 nm, 200 nm, or a 40+200 nm mix, plus an untreated control. "
        "Each file is a filtered 10x gene-cell count matrix stored as AnnData (.h5ad)."
    )
    log.log("")
    log.log("Matrix layout: rows = cells (barcodes), columns = genes; values = UMI counts.")
    log.log("")

    rows = []
    for sample_id, meta in cfg["samples"].items():
        f = paths["raw"] / f"{sample_id}.h5ad"
        if not f.exists():
            log.log(f"  [MISSING] {sample_id}.h5ad - run download_data.py first")
            continue
        ad = sc.read_h5ad(f, backed="r")
        n_cells, n_genes = ad.n_obs, ad.n_vars
        if "total_counts" in ad.obs.columns:
            total_umis = int(ad.obs["total_counts"].sum())
        else:
            total_umis = "n/a"
        if "n_genes_by_counts" in ad.obs.columns:
            median_genes = float(np.median(ad.obs["n_genes_by_counts"]))
        else:
            median_genes = float("nan")
        ad.file.close()
        rows.append(
            {
                "file": sample_id,
                "condition": meta["condition"],
                "cells": n_cells,
                "genes": n_genes,
                "total_UMIs": total_umis,
                "median_genes_per_cell": round(median_genes, 1),
            }
        )
        log.log(f"  Sample file: {sample_id}.h5ad")
        log.log(f"    Biological condition: {meta['condition']}")
        log.log(f"    Dimensions: {n_cells:,} cells x {n_genes:,} genes")
        log.log(f"    Total UMIs (mapped reads): {total_umis:,}" if isinstance(total_umis, int) else f"    Total UMIs: {total_umis}")
        if np.isnan(median_genes):
            log.log("    Median detected genes per cell: n/a")
        else:
            log.log(f"    Median detected genes per cell: {median_genes:.0f}")
        log.log("")

    if rows:
        log.log("Summary table (all samples):")
        summary = pd.DataFrame(rows)
        log.log(summary.to_string(index=False))
        log.log("")
    log.log(
        "QC filters (from config): min_genes={}, max_genes={}, min_counts={}, max_mt%={}. "
        "Cells failing these criteria are removed before integration.".format(
            cfg["qc"]["min_genes"],
            cfg["qc"]["max_genes"],
            cfg["qc"]["min_counts"],
            cfg["qc"]["max_mt_percent"],
        )
    )


def read_and_qc_sample(
    file_path: Path,
    sample_id: str,
    condition: str,
    qc_cfg: Dict,
    log: PipelineLogger,
) -> sc.AnnData:
    adata = sc.read_h5ad(file_path)
    n_before = adata.n_obs
    adata.obs["sample_id"] = sample_id
    adata.obs["condition"] = condition

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    adata = adata[
        (adata.obs["n_genes_by_counts"] >= qc_cfg["min_genes"])
        & (adata.obs["n_genes_by_counts"] <= qc_cfg["max_genes"])
        & (adata.obs["total_counts"] >= qc_cfg["min_counts"])
        & (adata.obs["pct_counts_mt"] <= qc_cfg["max_mt_percent"]),
        :
    ].copy()
    adata = _compact_adata_for_merge(adata)
    n_after = adata.n_obs
    pct_kept = 100.0 * n_after / n_before if n_before else 0
    log.log(
        f"  {sample_id} ({condition}): {n_before:,} -> {n_after:,} cells after QC "
        f"({pct_kept:.1f}% retained)"
    )
    return adata


def score_gene_modules(adata: sc.AnnData, log: PipelineLogger) -> None:
    """Score cell-cycle and immune modules on the full normalized gene space."""
    log.log("  Scoring gene modules on full normalized matrix (before HVG subset)...")

    s_genes = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "CDC6"]
    g2m_genes = ["HMGB2", "CDK1", "NUSAP1", "TOP2A", "MKI67", "BIRC5", "UBE2C", "CCNB1"]
    valid_s = [g for g in s_genes if g in adata.var_names]
    valid_g2m = [g for g in g2m_genes if g in adata.var_names]
    if valid_s and valid_g2m:
        sc.tl.score_genes_cell_cycle(adata, s_genes=valid_s, g2m_genes=valid_g2m)
        log.log(f"    Cell-cycle: {len(valid_s)} S genes, {len(valid_g2m)} G2M genes")
    else:
        log.log("  Warning: cell-cycle genes missing - scores set to NaN")
        adata.obs["S_score"] = np.nan
        adata.obs["G2M_score"] = np.nan

    ifn_genes = ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "OASL", "IFI6", "RSAD2"]
    valid_ifn = [g for g in ifn_genes if g in adata.var_names]
    if valid_ifn:
        sc.tl.score_genes(adata, gene_list=valid_ifn, score_name="IFN_score")
        log.log(f"    IFN signature: {len(valid_ifn)} genes")
    else:
        log.log("  Warning: IFN genes missing - IFN_score set to NaN")
        adata.obs["IFN_score"] = np.nan

    ag_genes = ["HLA-DRA", "HLA-DRB1", "CD74", "B2M", "TAP1", "TAP2", "HLA-DPA1", "HLA-DPB1"]
    valid_ag = [g for g in ag_genes if g in adata.var_names]
    if valid_ag:
        sc.tl.score_genes(adata, gene_list=valid_ag, score_name="antigen_presentation_score")
        log.log(f"    Antigen presentation: {len(valid_ag)} genes")
    else:
        log.log("  Warning: antigen presentation genes missing")
        adata.obs["antigen_presentation_score"] = np.nan


def merge_and_integrate(
    adata: sc.AnnData, cfg: Dict, paths: Dict[str, Path], log: PipelineLogger
) -> sc.AnnData:
    adata.obs_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    log.log(f"  Merged object: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    log.log("  Normalizing to 10,000 counts per cell and log1p transform...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    score_gene_modules(adata, log)

    pre_hvg_path = paths["processed"] / "integrated_pre_hvg.h5ad"
    log.log(f"  Saving pre-HVG object for ref.Rds mapping tools: {pre_hvg_path}")
    adata.write(pre_hvg_path)

    log.log(f"  Selecting top {cfg['preprocessing']['n_hvgs']} highly variable genes (HVG, seurat flavor)...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=cfg["preprocessing"]["n_hvgs"],
        flavor="seurat",
        batch_key="sample_id",
    )
    adata = adata[:, adata.var["highly_variable"]].copy()
    log.log(f"  Working matrix reduced to {adata.n_vars:,} HVGs")

    log.log("  Scaling expression (max_value=10)...")
    sc.pp.scale(adata, max_value=10)

    n_pcs = cfg["preprocessing"]["n_pcs"]
    log.log(f"  PCA ({n_pcs} components)...")
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")

    integration_method = cfg["preprocessing"].get("integration_method", "combat")
    if integration_method == "harmony":
        try:
            log.log("  Batch integration: Harmony on sample_id...")
            sc.external.pp.harmony_integrate(adata, key="sample_id")
            rep_key = "X_pca_harmony"
        except Exception:
            log.log("  Harmony unavailable - falling back to Combat on sample_id...")
            sc.pp.combat(adata, key="sample_id")
            rep_key = "X_pca"
    else:
        log.log("  Batch integration: Combat on sample_id (removes technical batch effects)...")
        sc.pp.combat(adata, key="sample_id")
        rep_key = "X_pca"

    k = cfg["preprocessing"]["neighbors_k"]
    log.log(f"  Neighbors graph (k={k}) and UMAP embedding...")
    sc.pp.neighbors(adata, use_rep=rep_key, n_neighbors=k)
    sc.tl.umap(adata)

    res = cfg["preprocessing"]["leiden_resolution"]
    log.log(f"  Leiden clustering (resolution={res})...")
    sc.tl.leiden(adata, resolution=res, key_added="cluster")
    n_clusters = adata.obs["cluster"].nunique()
    log.log(f"  Found {n_clusters} Leiden clusters")
    return adata


AZIMUTH_PANEL_TO_CELLTYPE = {
    "B": "B_cell",
    "CD4_T": "CD4_T",
    "CD8_T": "CD8_T_cytotoxic",
    "NK": "NK_cell",
    "Mono": "Monocyte_CD14",
    "DC": "DC",
    "other": "other",
    "other_T": "other_T",
}

AZIMUTH_L1_TO_MARKER = {
    "CD4 T": "CD4_T",
    "CD8 T": "CD8_T_cytotoxic",
    "B": "B_cell",
    "NK": "NK_cell",
    "Mono": "Monocyte_CD14",
    "DC": "DC",
    "other T": "other_T",
    "other": "other",
}

MODULE_SCORE_COLUMNS = ("S_score", "G2M_score", "IFN_score", "antigen_presentation_score")


def _curated_markers_from_config(cfg: Dict) -> Dict[str, List[str]]:
    return {
        key: genes
        for key, genes in cfg.get("markers", {}).items()
        if isinstance(genes, list)
    }


def _azimuth_marker_panel_path(cfg: Dict, paths: Dict[str, Path]) -> Path:
    az_cfg = cfg.get("azimuth", {})
    explicit = az_cfg.get("marker_panel_file")
    if explicit:
        panel_path = Path(explicit)
        if not panel_path.is_absolute():
            panel_path = Path.cwd() / panel_path
        return panel_path
    level = az_cfg.get("marker_extraction", {}).get("level", "l1")
    return paths["results"] / "tables" / f"azimuth_marker_panels_{level}.yaml"


def resolve_marker_dict(
    cfg: Dict, paths: Dict[str, Path], log: PipelineLogger
) -> Dict[str, List[str]]:
    """Load marker panels for pipeline annotation (Azimuth-derived or curated config)."""
    az_cfg = cfg.get("azimuth", {})
    curated = _curated_markers_from_config(cfg)

    if not az_cfg.get("use_panels_for_pipeline", False):
        log.log("  Marker source: curated literature panels from config.yaml (validation)")
        return curated

    panel_path = _azimuth_marker_panel_path(cfg, paths)
    if not panel_path.exists():
        log.log(
            f"  Warning: Azimuth marker panel not found ({panel_path}). "
            "Run scripts/extract_azimuth_markers.py first."
        )
        log.log("  Falling back to curated markers from config.yaml")
        return curated

    with open(panel_path, "r", encoding="utf-8") as f:
        raw_panels = yaml.safe_load(f) or {}

    marker_dict: Dict[str, List[str]] = {}
    for az_key, genes in raw_panels.items():
        if not isinstance(genes, list):
            continue
        celltype = AZIMUTH_PANEL_TO_CELLTYPE.get(az_key, az_key)
        marker_dict[celltype] = [str(g) for g in genes if g]

    for celltype in az_cfg.get("marker_fallback_types", ["Platelet"]):
        if celltype in marker_dict:
            continue
        fallback_genes = curated.get(celltype)
        if fallback_genes:
            marker_dict[celltype] = list(fallback_genes)

    log.log(f"  Marker source: Azimuth-derived panels ({panel_path.name})")
    for ct, genes in sorted(marker_dict.items()):
        log.log(f"    {ct}: {len(genes)} genes")
    return marker_dict


def marker_based_annotation(
    adata: sc.AnnData, marker_dict: Dict[str, List[str]], log: PipelineLogger
) -> None:
    log.log("  Assigning cell types by highest mean marker-gene score per cell...")
    gene_set = set(adata.var_names)
    score_matrix = {}
    for celltype, markers in marker_dict.items():
        valid = [g for g in markers if g in gene_set]
        if not valid:
            score_matrix[celltype] = np.full(adata.n_obs, -np.inf)
            continue
        score_matrix[celltype] = np.asarray(adata[:, valid].X.mean(axis=1)).ravel()
    scores_df = pd.DataFrame(score_matrix, index=adata.obs_names)
    adata.obs["cell_type_marker"] = scores_df.idxmax(axis=1).values
    counts = adata.obs["cell_type_marker"].value_counts()
    log.log("  Cell-type counts (literature markers, validation):")
    for ct, n in counts.items():
        log.log(f"    {ct}: {n:,}")


def composition_analysis(
    adata: sc.AnnData, paths: Dict[str, Path], log: PipelineLogger
) -> pd.DataFrame:
    comp = (
        adata.obs.groupby(["condition", "cell_type_primary"])
        .size()
        .reset_index(name="n_cells")
        .sort_values(["condition", "n_cells"], ascending=[True, False])
    )
    comp["fraction"] = comp["n_cells"] / comp.groupby("condition")["n_cells"].transform("sum")
    out_csv = paths["results"] / "tables" / "cell_composition_by_condition.csv"
    comp.to_csv(out_csv, index=False)
    log.log(f"  Saved table: {out_csv}")

    fig_path = paths["results"] / "figures" / "composition_barplot.png"
    plt.figure(figsize=(10, 5))
    sns.barplot(data=comp, x="condition", y="fraction", hue="cell_type_primary")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig(fig_path, dpi=300)
    plt.close()
    log.log(f"  Saved figure: {fig_path}")
    return comp


def differential_expression_by_celltype(
    adata: sc.AnnData, cfg: Dict, paths: Dict[str, Path], log: PipelineLogger
) -> pd.DataFrame:
    conditions = ["PSNP_40nm", "PSNP_200nm", "PSNP_mix_40_200"]
    all_de = []
    comparisons_run = 0
    for ctype in sorted(adata.obs["cell_type_primary"].unique()):
        if ctype in PRIMARY_CELL_TYPE_EXCLUDE:
            continue
        ad_ct = adata[adata.obs["cell_type_primary"] == ctype].copy()
        if ad_ct.n_obs < cfg["de"]["min_cells_per_group"]:
            continue

        for cond in conditions:
            n_test = (ad_ct.obs["condition"] == cond).sum()
            n_ctrl = (ad_ct.obs["condition"] == "control").sum()
            if min(n_test, n_ctrl) < cfg["de"]["min_cells_per_group"]:
                continue

            ad_ct.obs["de_group"] = np.where(ad_ct.obs["condition"] == cond, cond, "control")
            ad_sub = ad_ct[ad_ct.obs["de_group"].isin([cond, "control"])].copy()

            sc.tl.rank_genes_groups(
                ad_sub, groupby="de_group", groups=[cond], reference="control", method="wilcoxon"
            )
            de_df = sc.get.rank_genes_groups_df(ad_sub, group=cond)
            de_df["cell_type"] = ctype
            de_df["comparison"] = f"{cond}_vs_control"
            all_de.append(de_df)
            comparisons_run += 1
            log.log(f"  DE done: {ctype} - {cond} vs control (n_test={n_test}, n_ctrl={n_ctrl})")

    if not all_de:
        log.log("  No DE comparisons met minimum cell count - skipping DE table")
        return pd.DataFrame()
    de_all = pd.concat(all_de, ignore_index=True)
    out_csv = paths["results"] / "tables" / "differential_expression_all.csv"
    de_all.to_csv(out_csv, index=False)
    n_sig = (
        (de_all["pvals_adj"] < cfg["de"]["pval_adj_threshold"])
        & (de_all["logfoldchanges"].abs() > cfg["de"]["logfc_threshold"])
    ).sum()
    log.log(f"  Ran {comparisons_run} comparisons; saved {len(de_all):,} gene rows to {out_csv}")
    log.log(f"  Significant hits (padj<{cfg['de']['pval_adj_threshold']}, |logFC|>{cfg['de']['logfc_threshold']}): {n_sig:,}")
    return de_all


def _parse_overlap_fraction(overlap: str) -> float:
    try:
        hit, total = str(overlap).split("/")
        return int(hit) / max(int(total), 1)
    except (ValueError, AttributeError):
        return 0.0


def _short_pathway_label(term: str, gene_set: str, max_len: int = 52) -> str:
    prefix = GENE_SET_PREFIX.get(gene_set, "")
    label = str(term)
    if prefix and not label.upper().startswith(prefix.upper()):
        label = f"{prefix} | {label}"
    if len(label) > max_len:
        label = label[: max_len - 1] + "…"
    return label


def _significant_de_genes(
    df: pd.DataFrame, cfg: Dict, direction: str
) -> np.ndarray:
    padj = cfg["de"]["pval_adj_threshold"]
    lfc = cfg["de"]["logfc_threshold"]
    if direction == "UP":
        mask = (df["pvals_adj"] < padj) & (df["logfoldchanges"] > lfc)
    else:
        mask = (df["pvals_adj"] < padj) & (df["logfoldchanges"] < -lfc)
    return df.loc[mask, "names"].dropna().astype(str).unique()


def pathway_enrichment(
    de_all: pd.DataFrame,
    cfg: Dict,
    paths: Dict[str, Path],
    log: PipelineLogger,
    directions: Optional[tuple] = None,
) -> pd.DataFrame:
    if de_all.empty:
        log.log("  Skipped pathway enrichment (no DE results)")
        return pd.DataFrame()

    plot_cfg = cfg.get("pathway_plots", {})
    min_genes = int(plot_cfg.get("min_input_genes", 10))
    gene_sets = ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]
    run_directions = directions or ("UP", "DOWN")

    log.log(
        "  Running Enrichr (GO, KEGG, Reactome) on significant "
        f"{', '.join(run_directions)} DE genes per cell type x comparison..."
    )
    enr_frames = []
    for direction in run_directions:
        for (ctype, comp), df in de_all.groupby(["cell_type", "comparison"]):
            sig = _significant_de_genes(df, cfg, direction)
            if len(sig) < min_genes:
                continue

            for gset in gene_sets:
                try:
                    enr = gp.enrichr(
                        gene_list=list(sig), gene_sets=gset, organism="human", outdir=None
                    )
                    if enr.results is None or enr.results.empty:
                        continue
                    tmp = enr.results.copy()
                    tmp["cell_type"] = ctype
                    tmp["comparison"] = comp
                    tmp["gene_set"] = gset
                    tmp["direction"] = direction
                    enr_frames.append(tmp)
                    log.log(
                        f"    {direction}: {ctype} / {comp} / {gset} "
                        f"({len(sig)} input genes)"
                    )
                except Exception as exc:
                    log.log(f"    Enrichr failed for {ctype}/{comp}/{gset}/{direction}: {exc}")

    out_csv = paths["results"] / "tables" / "pathway_enrichment_all.csv"
    if not enr_frames:
        if out_csv.exists():
            log.log("  No new enrichment rows; keeping existing table")
            return pd.read_csv(out_csv)
        log.log("  No pathway enrichment results returned")
        return pd.DataFrame()

    new_rows = pd.concat(enr_frames, ignore_index=True)
    if out_csv.exists() and directions is not None:
        existing = pd.read_csv(out_csv)
        if "direction" not in existing.columns:
            existing["direction"] = "UP"
        combined = pd.concat([existing, new_rows], ignore_index=True)
        key_cols = ["cell_type", "comparison", "gene_set", "direction", "Term"]
        enr_all = combined.drop_duplicates(subset=key_cols, keep="last")
        log.log(f"  Merged with existing table -> {len(enr_all):,} rows total")
    else:
        enr_all = new_rows

    enr_all.to_csv(out_csv, index=False)
    log.log(f"  Saved: {out_csv} ({len(enr_all):,} rows)")
    return enr_all


def _prepare_pathway_heatmap_pivot(
    sub: pd.DataFrame,
    *,
    top_n: int,
    fdr_col: str,
    fdr_thr: float,
    column_order: List[str],
    column_labels: List[str],
) -> tuple[Optional[pd.DataFrame], List[str]]:
    """
    Build one heatmap column per exposure (aggregate GO/KEGG/Reactome per Term).
    Returns (pivot, y_labels) or (None, []) if no data.
    """
    if sub.empty:
        return None, []

    sub = sub[sub[fdr_col] < fdr_thr].copy()
    if sub.empty:
        return None, []

    sub["neg_log_fdr"] = -np.log10(sub[fdr_col].clip(lower=1e-300))

    # Best row per comparison x Term (lowest padj across gene-set databases).
    collapsed = (
        sub.sort_values(fdr_col)
        .drop_duplicates(subset=["comparison", "Term"], keep="first")
    )

    term_rank = (
        collapsed.groupby("Term")["neg_log_fdr"]
        .max()
        .sort_values(ascending=False)
        .head(top_n)
    )
    top_terms = term_rank.index.tolist()
    if not top_terms:
        return None, []

    label_rows = (
        collapsed[collapsed["Term"].isin(top_terms)]
        .sort_values(fdr_col)
        .drop_duplicates(subset=["Term"], keep="first")
    )
    label_map = {
        row["Term"]: _short_pathway_label(row["Term"], row["gene_set"])
        for _, row in label_rows.iterrows()
    }
    y_labels = [label_map.get(term, term) for term in top_terms]

    matrix = (
        collapsed[collapsed["Term"].isin(top_terms)]
        .pivot_table(
            index="Term",
            columns="comparison",
            values="neg_log_fdr",
            aggfunc="max",
        )
        .reindex(index=top_terms, columns=column_order)
        .fillna(0.0)
    )
    matrix.columns = pd.Index(column_labels, name="Exposure vs control")
    matrix.index = pd.Index(y_labels, name="Pathway")
    return matrix, y_labels


def plot_pathway_enrichment_figures(
    enr_all: pd.DataFrame, cfg: Dict, paths: Dict[str, Path], log: PipelineLogger
) -> List[Path]:
    """Heatmap: pathways x exposure (one column per condition), per cell type and direction."""
    if enr_all.empty:
        log.log("  Skipped pathway figures (no enrichment table)")
        return []

    if "direction" not in enr_all.columns:
        enr_all = enr_all.copy()
        enr_all["direction"] = "UP"

    plot_cfg = cfg.get("pathway_plots", {})
    top_n = int(plot_cfg.get("top_terms", 20))
    fdr_thr = float(plot_cfg.get("fdr_threshold", 0.05))
    fdr_col = "Adjusted P-value"

    comp_order = [c for c, _ in EXPOSURE_COMPARISONS]
    comp_labels = [label for _, label in EXPOSURE_COMPARISONS]

    fig_dir = paths["results"] / "figures" / "pathway_enrichment"
    fig_dir.mkdir(parents=True, exist_ok=True)
    saved: List[Path] = []

    for (ctype, direction), sub in enr_all.groupby(["cell_type", "direction"], observed=False):
        pivot, _ = _prepare_pathway_heatmap_pivot(
            sub,
            top_n=top_n,
            fdr_col=fdr_col,
            fdr_thr=fdr_thr,
            column_order=comp_order,
            column_labels=comp_labels,
        )
        if pivot is None or pivot.empty:
            continue

        n_terms = len(pivot)
        n_cols = pivot.shape[1]
        fig_w = max(5.0, 1.15 * n_cols + 2.5)
        fig_h = max(4.5, 0.32 * n_terms + 1.8)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        vmax = max(float(pivot.values.max()), 3.0)
        data = pivot.values.astype(float)

        im = ax.imshow(
            data,
            aspect="auto",
            cmap="YlOrRd",
            vmin=0,
            vmax=vmax,
            origin="upper",
            extent=[0, n_cols, n_terms, 0],
        )
        ax.set_xticks(np.arange(n_cols) + 0.5)
        ax.set_xticklabels(list(pivot.columns), rotation=0, ha="center")
        ax.set_yticks(np.arange(n_terms) + 0.5)
        ax.set_yticklabels(list(pivot.index), fontsize=8)
        ax.set_xlim(0, n_cols)
        ax.set_ylim(n_terms, 0)
        fig.colorbar(im, ax=ax, shrink=0.6, label="-log10(adjusted p-value)")

        ax.set_xlabel("Exposure vs control")
        ax.set_ylabel("Enriched pathway")
        ax.tick_params(axis="x", labelsize=9)

        title_ct = str(ctype).replace("_", " ")
        ax.set_title(
            f"{title_ct} — pathway enrichment ({direction}-regulated genes)",
            fontsize=11,
            pad=10,
        )
        plt.tight_layout()

        out_path = fig_dir / f"pathways_{ctype}_{direction}.png"
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        saved.append(out_path)
        log.log(f"  Saved pathway figure: {out_path} ({n_cols} exposure columns)")

    log.log(f"  Pathway enrichment figures: {len(saved)} saved under {fig_dir}")
    return saved


PATHWAY_BIOLOGY_HINTS: List[tuple] = [
    (
        ("cytokine", "interleukin", "il-10", "il-17", "tnf", "inflammatory", "nf-kappa"),
        {
            "UP": (
                "Higher expression of cytokine/inflammatory genes in exposed cells suggests "
                "amplified innate and adaptive immune signalling — a plausible systemic "
                "inflammatory response to particulate nanoplastic stress."
            ),
            "DOWN": (
                "Reduced cytokine pathway activity vs control may reflect immune suppression, "
                "cell exhaustion, or compensatory down-regulation after initial activation."
            ),
        },
    ),
    (
        ("interferon", "antiviral", "isg", "defense response to virus"),
        {
            "UP": (
                "Interferon/stress-response programmes are activated — consistent with "
                "pattern-recognition receptor sensing of foreign particulate material."
            ),
            "DOWN": (
                "Lower interferon signatures vs control indicate that IFN is not the dominant "
                "exposure response in this cell type, or that other pathways supersede it."
            ),
        },
    ),
    (
        ("apoptosis", "cell death", "p53", "necroptosis", "ferroptosis"),
        {
            "UP": (
                "Up-regulated cell-death pathways may indicate cytotoxic stress from "
                "nanoparticle uptake, oxidative damage, or activation of damage sensors."
            ),
            "DOWN": (
                "Down-regulated death pathways could reflect survival signalling or "
                "selection of stress-resistant subpopulations under PSNP exposure."
            ),
        },
    ),
    (
        ("phagocyt", "fcgr", "complement", "neutrophil degranulation", "myeloid"),
        {
            "UP": (
                "Enriched phagocytic/myeloid programmes align with particle uptake routes "
                "and monocyte/macrophage engagement — especially relevant for larger (200 nm) PSNP."
            ),
            "DOWN": (
                "Reduced phagocytic gene expression may indicate impaired clearance capacity "
                "or shifted myeloid functional states under nanoplastic exposure."
            ),
        },
    ),
    (
        ("oxidative", "reactive oxygen", "ros", "hypoxia", "mtor", "metabolism"),
        {
            "UP": (
                "Metabolic and oxidative-stress pathways suggest cellular adaptation to "
                "particle-induced ROS and energy demand during stress responses."
            ),
            "DOWN": (
                "Down-regulated metabolic pathways may reflect functional impairment or "
                "metabolic reprogramming away from homeostatic states."
            ),
        },
    ),
    (
        ("adhesion", "focal adhesion", "ecm", "integrin", "migration", "chemokine"),
        {
            "UP": (
                "Adhesion and migration programmes may promote immune-cell trafficking and "
                "tissue homing — relevant for how circulating PBMC could respond upon "
                "re-entering inflamed tissue after nanoplastic exposure."
            ),
            "DOWN": (
                "Reduced adhesion/migration signalling could limit immune surveillance "
                "or reflect altered cell–matrix interactions under stress."
            ),
        },
    ),
    (
        ("antigen", "mhc", "hla", "presentation", "t cell receptor"),
        {
            "UP": (
                "Enhanced antigen-presentation pathways support adaptive immune activation "
                "and could influence T-cell priming after particulate exposure."
            ),
            "DOWN": (
                "Lower antigen-presentation activity may reduce adaptive immune engagement "
                "and impair effective immune surveillance."
            ),
        },
    ),
]


def _pathway_biology_note(term: str, direction: str) -> str:
    """Return a short biological note for a pathway term and regulation direction."""
    text = str(term).lower()
    direction = str(direction).upper()
    for keywords, notes in PATHWAY_BIOLOGY_HINTS:
        if any(kw in text for kw in keywords):
            return notes.get(direction, notes.get("UP", ""))
    if direction == "UP":
        return (
            "Genes in this pathway are higher in the exposed condition than in control — "
            "may contribute to nanoplastic-induced transcriptional reprogramming in PBMC."
        )
    return (
        "Genes in this pathway are lower in the exposed condition than in control — "
        "may reflect suppression or redistribution of this programme under PSNP stress."
    )


def export_pathway_enrichment_summary(
    enr_all: pd.DataFrame, cfg: Dict, paths: Dict[str, Path], log: PipelineLogger
) -> pd.DataFrame:
    """Export a readable table of enriched pathways with UP/DOWN direction vs control."""
    if enr_all.empty:
        log.log("  Skipped pathway summary (no enrichment table)")
        return pd.DataFrame()

    plot_cfg = cfg.get("pathway_plots", {})
    fdr_thr = float(plot_cfg.get("fdr_threshold", 0.05))
    top_n = int(plot_cfg.get("top_terms", 20))
    fdr_col = "Adjusted P-value"
    comp_labels = {comp: label for comp, label in EXPOSURE_COMPARISONS}

    enr = enr_all.copy()
    if "direction" not in enr.columns:
        enr["direction"] = "UP"
    enr = enr[enr[fdr_col] < fdr_thr].copy()
    if enr.empty:
        log.log("  Skipped pathway summary (no terms pass FDR threshold)")
        return pd.DataFrame()

    rows: List[Dict] = []
    for (ctype, comp, direction), sub in enr.groupby(
        ["cell_type", "comparison", "direction"], observed=False
    ):
        ranked = (
            sub.sort_values(fdr_col)
            .drop_duplicates(subset=["Term"], keep="first")
            .head(top_n)
        )
        for _, r in ranked.iterrows():
            term = str(r["Term"])
            rows.append(
                {
                    "cell_type": ctype,
                    "comparison": comp,
                    "exposure": comp_labels.get(comp, comp),
                    "direction": direction,
                    "regulation_vs_control": (
                        "UP in exposure (higher than control)"
                        if direction == "UP"
                        else "DOWN in exposure (lower than control)"
                    ),
                    "gene_set": r.get("gene_set", ""),
                    "pathway": term,
                    "adjusted_p_value": r[fdr_col],
                    "overlap": r.get("Overlap", ""),
                    "genes": r.get("Genes", ""),
                    "biological_note": _pathway_biology_note(term, direction),
                }
            )

    summary = pd.DataFrame(rows)
    out = paths["results"] / "tables" / "pathway_enrichment_summary.csv"
    summary.to_csv(out, index=False)
    log.log(f"  Saved: {out.name} ({len(summary):,} rows, UP + DOWN vs control)")
    return summary


def write_pathway_enrichment_interpretation(
    summary: pd.DataFrame,
    paths: Dict[str, Path],
    log: PipelineLogger,
    *,
    top_n_per_group: int = 5,
    cfg: Optional[Dict] = None,
) -> Optional[Path]:
    """Write an English report linking enriched pathways to direction and tissue-level impact."""
    if summary is None or summary.empty:
        log.log("  Skipped pathway interpretation (empty summary)")
        return None

    deliverables = paths["project_root"] / "deliverables"
    deliverables.mkdir(parents=True, exist_ok=True)
    out = deliverables / "Pathway_Enrichment_Interpretation_EN.md"

    lines = [
        "# Pathway Enrichment — Direction and Biological Interpretation",
        "",
        "**Project:** Single-Cell Analysis of Immune Response to Nanoplastic Particles",
        "**Source table:** `results/tables/pathway_enrichment_summary.csv`",
        "**Figures:** `results/figures/pathway_enrichment/pathways_{cell_type}_{UP|DOWN}.png`",
        "",
        "---",
        "",
        "## 1. How to read these results",
        "",
        "Differential expression (Wilcoxon) compares each PSNP exposure to **untreated control** "
        "within the same cell type. Pathway enrichment (Enrichr: GO, KEGG, Reactome) is run "
        "separately on:",
        "",
        "| Direction | Input genes | Meaning vs control |",
        "|-----------|-------------|-------------------|",
        "| **UP** | Significant DE genes with log2FC > threshold | Pathway genes are **higher** "
        "in the exposed condition — activation or induction. |",
        "| **DOWN** | Significant DE genes with log2FC < −threshold | Pathway genes are **lower** "
        "in the exposed condition — suppression or loss of programme. |",
        "",
        "Each row in the summary table lists one enriched pathway for a given "
        "**cell type × exposure × direction** combination.",
        "",
        "---",
        "",
        "## 2. Linking pathways to nanoplastic tissue/organism effects",
        "",
        "PBMC are circulating immune cells; their transcriptional state reflects how the "
        "immune system senses nanoplastics **in blood**. Enriched pathways can be connected "
        "to potential organism-level consequences as follows:",
        "",
        "- **UP cytokine / TNF / IL signalling** → systemic inflammatory tone; may correlate "
        "with fever, fatigue, or tissue inflammation if activated cells traffic to organs.",
        "- **UP phagocytosis / myeloid programmes** → particle uptake and innate clearance "
        "attempts; larger PSNP (200 nm) often engage monocyte/macrophage routes.",
        "- **UP oxidative stress / metabolism** → cellular damage and energy reallocation; "
        "relevant for vascular endothelium and organ toxicity hypotheses.",
        "- **DOWN antigen presentation / HLA** → reduced capacity to activate adaptive immunity; "
        "may impair effective immune surveillance.",
        "- **DOWN interferon** → if IFN pathways are suppressed while inflammatory genes rise, "
        "the response may be skewed toward myeloid rather than antiviral programmes.",
        "",
        "---",
        "",
        "## 3. Top enriched pathways by cell type and exposure",
        "",
    ]

    top_per_group = top_n_per_group
    cfg = cfg or load_config()
    for ctype in sorted(summary["cell_type"].unique()):
        lines.append(f"### {str(ctype).replace('_', ' ')}")
        lines.append("")
        sub_ct = summary[summary["cell_type"] == ctype]
        for comp in [c for c, _ in EXPOSURE_COMPARISONS]:
            sub_comp = sub_ct[sub_ct["comparison"] == comp]
            if sub_comp.empty:
                continue
            exposure = dict(EXPOSURE_COMPARISONS).get(comp, comp)
            lines.append(f"#### {exposure} vs control")
            lines.append("")
            for direction in ("UP", "DOWN"):
                sub_dir = sub_comp[sub_comp["direction"] == direction].head(top_per_group)
                if sub_dir.empty:
                    continue
                lines.append(f"**{direction}-regulated pathways**")
                lines.append("")
                for _, row in sub_dir.iterrows():
                    padj = row["adjusted_p_value"]
                    padj_s = f"{padj:.2e}" if pd.notna(padj) else "n/a"
                    gs = GENE_SET_PREFIX.get(str(row.get("gene_set", "")), "")
                    pathway_label = str(row["pathway"])
                    if gs and not pathway_label.upper().startswith(gs.upper()):
                        pathway_label = f"{gs} | {pathway_label}"
                    lines.append(f"- **{pathway_label}** (padj = {padj_s})")
                    lines.append(f"  - {row['biological_note']}")
                lines.append("")

    lines.extend(
        [
            "---",
            "",
            "## 4. Full tabular export",
            "",
            "See `results/tables/pathway_enrichment_summary.csv` for all significant pathways "
            f"(top {int(cfg.get('pathway_plots', {}).get('top_terms', 20))} per "
            "cell type × comparison × direction) with gene overlap and biological notes.",
            "",
            "Raw Enrichr output: `results/tables/pathway_enrichment_all.csv` (includes "
            "`direction` column: UP or DOWN).",
            "",
        ]
    )

    out.write_text("\n".join(lines), encoding="utf-8")
    log.log(f"  Saved: {out}")
    return out


def _plot_annotation_crosstab_heatmap(
    ctab: pd.DataFrame,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
    out_path: Path,
    normalize_rows: bool = False,
    cbar_label: str = "Cell count",
) -> None:
    """Plot a contingency matrix as a heatmap (counts or row-normalized %)."""
    plot_df = ctab.copy()
    for label in ("All", "Total"):
        if label in plot_df.index:
            plot_df = plot_df.drop(index=label)
        if label in plot_df.columns:
            plot_df = plot_df.drop(columns=label)

    if plot_df.empty:
        return

    if normalize_rows:
        row_sums = plot_df.sum(axis=1).replace(0, np.nan)
        plot_df = plot_df.div(row_sums, axis=0) * 100.0
        plot_df = plot_df.fillna(0.0)
        fmt = ".1f"
        cbar_label = "% of marker-based type"
    else:
        fmt = ".0f"

    fig_w = max(8.0, 0.55 * plot_df.shape[1] + 3.0)
    fig_h = max(5.0, 0.45 * plot_df.shape[0] + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        plot_df,
        ax=ax,
        cmap="Blues",
        linewidths=0.3,
        linecolor="white",
        annot=True,
        fmt=fmt,
        cbar_kws={"label": cbar_label, "shrink": 0.75},
    )
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _export_codi_marker_crosstabs(
    marker: pd.Series,
    codi: pd.Series,
    paths: Dict[str, Path],
    log: PipelineLogger,
) -> pd.DataFrame:
    """Save CoDi vs marker contingency tables (counts, row %, long format)."""
    valid = codi.astype(str) != "NA"
    marker_v = marker[valid].astype(str)
    codi_v = codi[valid].astype(str)

    ctab = pd.crosstab(marker_v, codi_v, margins=True)
    ctab.to_csv(paths["results"] / "tables" / "annotation_crosstab_marker_codi.csv")
    log.log("  Saved: annotation_crosstab_marker_codi.csv")

    ctab_no_margin = pd.crosstab(marker_v, codi_v, margins=False)
    row_pct = ctab_no_margin.div(ctab_no_margin.sum(axis=1).replace(0, np.nan), axis=0) * 100.0
    row_pct = row_pct.fillna(0.0)
    row_pct.to_csv(
        paths["results"] / "tables" / "annotation_crosstab_marker_codi_row_pct.csv"
    )
    log.log("  Saved: annotation_crosstab_marker_codi_row_pct.csv")

    long_rows: List[Dict] = []
    for m_type in ctab_no_margin.index:
        row_total = int(ctab_no_margin.loc[m_type].sum())
        for c_type in ctab_no_margin.columns:
            n = int(ctab_no_margin.loc[m_type, c_type])
            if n == 0:
                continue
            long_rows.append(
                {
                    "marker_cell_type": m_type,
                    "codi_cell_type": c_type,
                    "n_cells": n,
                    "pct_of_marker_row": round(100.0 * n / row_total, 2) if row_total else 0.0,
                    "pct_of_all_mapped": round(
                        100.0 * n / len(marker_v), 2
                    ),
                }
            )
    mapping = pd.DataFrame(long_rows).sort_values(
        ["marker_cell_type", "n_cells"], ascending=[True, False]
    )
    mapping.to_csv(
        paths["results"] / "tables" / "annotation_codi_marker_mapping.csv", index=False
    )
    log.log(f"  Saved: annotation_codi_marker_mapping.csv ({len(mapping)} pairs)")
    return ctab_no_margin


def _export_codi_ref_crosstabs(
    ref: pd.Series,
    codi: pd.Series,
    paths: Dict[str, Path],
    log: PipelineLogger,
) -> pd.DataFrame:
    """Save CoDi vs ref.Rds/Seurat contingency tables (counts, row %, long format)."""
    ref_s = ref.astype(str)
    codi_s = codi.astype(str)
    valid = (codi_s != "NA") & (ref_s != "NA") & (ref_s != "unmapped")
    ref_v = ref_s[valid]
    codi_v = codi_s[valid]

    if ref_v.empty:
        log.log("  Skipped ref vs CoDi crosstabs (no mapped ref.Rds + CoDi cells)")
        return pd.DataFrame()

    ctab = pd.crosstab(ref_v, codi_v, margins=True)
    ctab.to_csv(paths["results"] / "tables" / "annotation_crosstab_ref_codi.csv")
    log.log("  Saved: annotation_crosstab_ref_codi.csv")

    ctab_no_margin = pd.crosstab(ref_v, codi_v, margins=False)
    row_pct = ctab_no_margin.div(ctab_no_margin.sum(axis=1).replace(0, np.nan), axis=0) * 100.0
    row_pct = row_pct.fillna(0.0)
    row_pct.to_csv(
        paths["results"] / "tables" / "annotation_crosstab_ref_codi_row_pct.csv"
    )
    log.log("  Saved: annotation_crosstab_ref_codi_row_pct.csv")

    long_rows: List[Dict] = []
    for r_type in ctab_no_margin.index:
        row_total = int(ctab_no_margin.loc[r_type].sum())
        for c_type in ctab_no_margin.columns:
            n = int(ctab_no_margin.loc[r_type, c_type])
            if n == 0:
                continue
            long_rows.append(
                {
                    "ref_cell_type": r_type,
                    "codi_cell_type": c_type,
                    "n_cells": n,
                    "pct_of_ref_row": round(100.0 * n / row_total, 2) if row_total else 0.0,
                    "pct_of_all_mapped": round(100.0 * n / len(ref_v), 2),
                }
            )
    mapping = pd.DataFrame(long_rows).sort_values(
        ["ref_cell_type", "n_cells"], ascending=[True, False]
    )
    mapping.to_csv(
        paths["results"] / "tables" / "annotation_codi_ref_mapping.csv", index=False
    )
    log.log(f"  Saved: annotation_codi_ref_mapping.csv ({len(mapping)} pairs)")
    return ctab_no_margin


SIZE_EFFECT_CLASSES = (
    "unique_40nm",
    "unique_200nm",
    "shared_40_200",
    "shared_all_three",
    "mix_only_emergent",
)

# Primary DE comparison used for logFC / direction when exporting genes per class.
SIZE_EFFECT_PRIMARY_COMPARISON = {
    "unique_40nm": "PSNP_40nm_vs_control",
    "unique_200nm": "PSNP_200nm_vs_control",
    "shared_40_200": "PSNP_200nm_vs_control",
    "shared_all_three": "PSNP_mix_40_200_vs_control",
    "mix_only_emergent": "PSNP_mix_40_200_vs_control",
}

SIZE_EFFECT_INTERPRETATION = {
    "unique_40nm": (
        "Response specific to 40 nm PSNP; may reflect higher surface-area-to-volume "
        "ratio and distinct uptake/signaling vs larger particles."
    ),
    "unique_200nm": (
        "Response specific to 200 nm PSNP; often aligned with stronger myeloid "
        "engagement (phagocytosis, PRR signaling)."
    ),
    "shared_40_200": (
        "Significant for both 40 nm and 200 nm solo exposures but not in the mix; "
        "mix may mask or redirect this program."
    ),
    "shared_all_three": (
        "Core PSNP response module: significant across 40 nm, 200 nm, and mix — "
        "particle-size-independent transcriptional program."
    ),
    "mix_only_emergent": (
        "Emergent mixture effect: significant only when both particle sizes are "
        "present; not a simple sum of solo exposures."
    ),
}


def _significant_de_mask(de_all: pd.DataFrame, cfg: Dict) -> pd.Series:
    return (de_all["pvals_adj"] < cfg["de"]["pval_adj_threshold"]) & (
        de_all["logfoldchanges"].abs() > cfg["de"]["logfc_threshold"]
    )


def _size_effect_sets_for_celltype(s: pd.DataFrame) -> Dict[str, set]:
    s40 = set(s[s["comparison"] == "PSNP_40nm_vs_control"]["names"])
    s200 = set(s[s["comparison"] == "PSNP_200nm_vs_control"]["names"])
    smix = set(s[s["comparison"] == "PSNP_mix_40_200_vs_control"]["names"])
    return {
        "unique_40nm": s40 - s200 - smix,
        "unique_200nm": s200 - s40 - smix,
        "shared_40_200": (s40 & s200) - smix,
        "shared_all_three": s40 & s200 & smix,
        "mix_only_emergent": smix - s40 - s200,
    }


def _run_enrichr_batches(
    gene_list: List[str],
    gene_sets: List[str],
    log: PipelineLogger,
    label: str,
) -> List[pd.DataFrame]:
    frames: List[pd.DataFrame] = []
    for gset in gene_sets:
        try:
            enr = gp.enrichr(
                gene_list=gene_list, gene_sets=gset, organism="human", outdir=None
            )
            if enr.results is None or enr.results.empty:
                continue
            tmp = enr.results.copy()
            tmp["gene_set"] = gset
            frames.append(tmp)
        except Exception as exc:
            log.log(f"    Enrichr failed for {label} / {gset}: {exc}")
    return frames


def _export_size_specific_genes(
    de_all: pd.DataFrame,
    sig: pd.DataFrame,
    cfg: Dict,
    paths: Dict[str, Path],
    log: PipelineLogger,
) -> pd.DataFrame:
    padj_thr = cfg["de"]["pval_adj_threshold"]
    lfc_thr = cfg["de"]["logfc_threshold"]
    gene_rows: List[Dict] = []

    for ctype in sorted(sig["cell_type"].unique()):
        s = sig[sig["cell_type"] == ctype]
        effect_sets = _size_effect_sets_for_celltype(s)
        de_ct = de_all[de_all["cell_type"] == ctype]

        for effect_class, genes in effect_sets.items():
            if not genes:
                continue
            primary_comp = SIZE_EFFECT_PRIMARY_COMPARISON[effect_class]
            de_primary = de_ct[de_ct["comparison"] == primary_comp].set_index("names")

            for gene in sorted(genes):
                row: Dict = {
                    "cell_type": ctype,
                    "effect_class": effect_class,
                    "gene": gene,
                    "primary_comparison": primary_comp,
                }
                if gene in de_primary.index:
                    lfc = float(de_primary.at[gene, "logfoldchanges"])
                    padj = float(de_primary.at[gene, "pvals_adj"])
                    row["logfoldchanges"] = lfc
                    row["pvals_adj"] = padj
                    if padj < padj_thr and lfc > lfc_thr:
                        row["direction"] = "UP"
                    elif padj < padj_thr and lfc < -lfc_thr:
                        row["direction"] = "DOWN"
                    else:
                        row["direction"] = "NS_in_primary"
                else:
                    row["logfoldchanges"] = np.nan
                    row["pvals_adj"] = np.nan
                    row["direction"] = "NA"
                gene_rows.append(row)

    genes_df = pd.DataFrame(gene_rows)
    out_csv = paths["results"] / "tables" / "size_specific_genes.csv"
    genes_df.to_csv(out_csv, index=False)
    log.log(f"  Saved: {out_csv} ({len(genes_df):,} gene rows)")
    return genes_df


def size_specific_pathway_enrichment(
    genes_df: pd.DataFrame, cfg: Dict, paths: Dict[str, Path], log: PipelineLogger
) -> pd.DataFrame:
    if genes_df.empty:
        log.log("  Skipped size-specific pathway enrichment (no gene lists)")
        return pd.DataFrame()

    plot_cfg = cfg.get("pathway_plots", {})
    min_genes = int(plot_cfg.get("min_input_genes", 10))
    gene_sets = ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]

    log.log(
        "  Running Enrichr on UP-regulated genes per size-effect class "
        f"(min {min_genes} genes)..."
    )
    enr_frames: List[pd.DataFrame] = []
    for (ctype, effect_class), sub in genes_df.groupby(
        ["cell_type", "effect_class"], observed=False
    ):
        up_genes = sub.loc[sub["direction"] == "UP", "gene"].dropna().astype(str).unique()
        if len(up_genes) < min_genes:
            continue

        label = f"{ctype}/{effect_class}"
        for tmp in _run_enrichr_batches(list(up_genes), gene_sets, log, label):
            tmp["cell_type"] = ctype
            tmp["effect_class"] = effect_class
            tmp["direction"] = "UP"
            tmp["n_input_genes"] = len(up_genes)
            enr_frames.append(tmp)
            log.log(f"    UP: {label} / {tmp['gene_set'].iloc[0]} ({len(up_genes)} genes)")

    out_csv = paths["results"] / "tables" / "size_specific_pathway_enrichment.csv"
    if not enr_frames:
        log.log("  No size-specific pathway enrichment results returned")
        return pd.DataFrame()

    enr_all = pd.concat(enr_frames, ignore_index=True)
    enr_all.to_csv(out_csv, index=False)
    log.log(f"  Saved: {out_csv} ({len(enr_all):,} rows)")
    return enr_all


def _build_size_specific_interpretation(
    summary: pd.DataFrame,
    genes_df: pd.DataFrame,
    enr_df: pd.DataFrame,
    cfg: Dict,
    paths: Dict[str, Path],
    log: PipelineLogger,
) -> pd.DataFrame:
    top_n_genes = int(cfg.get("size_specific", {}).get("top_genes_per_class", 8))
    top_n_paths = int(cfg.get("size_specific", {}).get("top_pathways_per_class", 5))
    fdr_thr = float(cfg.get("pathway_plots", {}).get("fdr_threshold", 0.05))
    fdr_col = "Adjusted P-value"

    rows: List[Dict] = []
    for _, srow in summary.iterrows():
        ctype = srow["cell_type"]
        effect_class = srow["effect_class"]
        n_genes = int(srow["n_genes"])

        sub_genes = genes_df[
            (genes_df["cell_type"] == ctype) & (genes_df["effect_class"] == effect_class)
        ].copy()
        sub_genes = sub_genes.dropna(subset=["logfoldchanges"])
        sub_genes["abs_lfc"] = sub_genes["logfoldchanges"].abs()
        top_genes = (
            sub_genes.sort_values("abs_lfc", ascending=False)["gene"]
            .head(top_n_genes)
            .tolist()
        )

        top_pathways: List[str] = []
        if not enr_df.empty:
            sub_enr = enr_df[
                (enr_df["cell_type"] == ctype)
                & (enr_df["effect_class"] == effect_class)
                & (enr_df[fdr_col] < fdr_thr)
            ].sort_values(fdr_col)
            for _, erow in sub_enr.head(top_n_paths).iterrows():
                top_pathways.append(
                    f"{_short_pathway_label(erow['Term'], erow['gene_set'])} "
                    f"(padj={erow[fdr_col]:.2e})"
                )

        note = SIZE_EFFECT_INTERPRETATION.get(effect_class, "")
        if n_genes == 0:
            note = "No significant genes in this class for this cell type."

        rows.append(
            {
                "cell_type": ctype,
                "effect_class": effect_class,
                "n_genes": n_genes,
                "top_genes": "; ".join(top_genes) if top_genes else "",
                "top_pathways": " | ".join(top_pathways) if top_pathways else "",
                "interpretation": note,
            }
        )

    interp_df = pd.DataFrame(rows)
    out_csv = paths["results"] / "tables" / "size_specific_interpretation.csv"
    interp_df.to_csv(out_csv, index=False)
    log.log(f"  Saved: {out_csv} ({len(interp_df)} rows)")
    return interp_df


def size_specific_effects(
    de_all: pd.DataFrame, cfg: Dict, paths: Dict[str, Path], log: PipelineLogger
) -> pd.DataFrame:
    if de_all.empty:
        log.log("  Skipped size-specific summary (no DE)")
        return pd.DataFrame()

    log.log("  Classifying DE genes by particle size (unique 40nm / 200nm / shared / mix-only)...")
    sig = de_all.loc[_significant_de_mask(de_all, cfg), ["cell_type", "comparison", "names"]].drop_duplicates()

    results = []
    for ctype in sig["cell_type"].unique():
        s = sig[sig["cell_type"] == ctype]
        rows = _size_effect_sets_for_celltype(s)
        for k, v in rows.items():
            results.append({"cell_type": ctype, "effect_class": k, "n_genes": len(v)})
        log.log(f"    {ctype}: " + ", ".join(f"{k}={len(v)}" for k, v in rows.items()))

    summary = pd.DataFrame(results)
    out_csv = paths["results"] / "tables" / "size_specific_effects_summary.csv"
    summary.to_csv(out_csv, index=False)
    log.log(f"  Saved: {out_csv}")

    genes_df = _export_size_specific_genes(de_all, sig, cfg, paths, log)
    enr_df = size_specific_pathway_enrichment(genes_df, cfg, paths, log)
    _build_size_specific_interpretation(summary, genes_df, enr_df, cfg, paths, log)
    return summary


def attach_azimuth_labels(
    adata: sc.AnnData, paths: Dict[str, Path], log: PipelineLogger, cfg: Optional[Dict] = None
) -> bool:
    """Merge Seurat/Azimuth labels (ref.Rds) for cross-validation tables."""
    if cfg is None:
        cfg = load_config()
    return load_seurat_annotations(adata, paths, cfg, log)


def _export_module_score_tables(
    adata: sc.AnnData, paths: Dict[str, Path], log: PipelineLogger
) -> Dict[str, pd.DataFrame]:
    """Aggregate module scores by condition and by condition x cell type."""
    tables: Dict[str, pd.DataFrame] = {}
    score_cols = [c for c in MODULE_SCORE_COLUMNS if c in adata.obs.columns]
    if not score_cols:
        log.log("  Warning: no module score columns in adata.obs")
        return tables

    log.log("  [1/6] Module scores by condition (global means)...")
    by_cond = adata.obs.groupby("condition", observed=False)[score_cols].mean().reset_index()
    by_cond.to_csv(paths["results"] / "tables" / "module_scores_by_condition.csv", index=False)
    tables["by_condition"] = by_cond
    log.log(f"  Saved: module_scores_by_condition.csv")

    cc = by_cond[["condition", "S_score", "G2M_score"]].copy()
    cc.to_csv(paths["results"] / "tables" / "cell_cycle_scores_by_condition.csv", index=False)
    tables["cell_cycle"] = cc

    ifn = by_cond[["condition", "IFN_score"]].copy()
    ifn.to_csv(paths["results"] / "tables" / "ifn_scores_by_condition.csv", index=False)
    tables["ifn"] = ifn

    log.log("  [2/6] Module scores by condition x cell type...")
    by_ct = (
        adata.obs.groupby(["condition", "cell_type_primary"], observed=False)[score_cols]
        .mean()
        .reset_index()
    )
    by_ct.to_csv(paths["results"] / "tables" / "module_scores_by_condition_celltype.csv", index=False)
    tables["by_condition_celltype"] = by_ct

    ag = by_ct[["condition", "cell_type_primary", "antigen_presentation_score"]].copy()
    ag = ag.rename(columns={"cell_type_primary": "cell_type_marker"})
    ag.to_csv(paths["results"] / "tables" / "antigen_presentation_scores.csv", index=False)
    tables["antigen"] = ag
    log.log("  Saved: module_scores_by_condition_celltype.csv, antigen_presentation_scores.csv")
    return tables


def _export_pseudobulk(adata: sc.AnnData, paths: Dict[str, Path], log: PipelineLogger) -> pd.DataFrame:
    log.log("  [3/6] Pseudobulk counts (sum UMIs per condition x cell type)...")
    if "counts" not in adata.layers:
        log.log("  Warning: counts layer missing — cannot build pseudobulk matrix.")
        return pd.DataFrame()

    counts_layer = adata.layers["counts"]
    matrix = counts_layer.toarray() if hasattr(counts_layer, "toarray") else counts_layer
    pseudobulk = (
        pd.DataFrame(matrix, index=adata.obs_names, columns=adata.var_names)
        .assign(
            condition=adata.obs["condition"].astype(str).values,
            cell_type=adata.obs["cell_type_primary"].astype(str).values,
        )
        .groupby(["condition", "cell_type"], observed=False)
        .sum()
    )
    p = paths["results"] / "tables" / "pseudobulk_counts_condition_celltype.csv"
    pseudobulk.to_csv(p)
    log.log(f"  Saved: {p} ({pseudobulk.shape[0]} groups)")
    return pseudobulk


def _export_annotation_agreement(
    adata: sc.AnnData, paths: Dict[str, Path], log: PipelineLogger, has_azimuth: bool
) -> pd.DataFrame:
    log.log("  [4/6] Annotation cross-validation (ref.Rds / literature markers / CoDi)...")
    marker = adata.obs["cell_type_marker"].astype(str)
    codi = adata.obs["cell_type_codi_norm"].astype(str)
    valid_codi = codi != "NA"
    ref = (
        adata.obs["cell_type_ref"].astype(str)
        if "cell_type_ref" in adata.obs.columns
        else adata.obs.get("cell_type_primary", marker).astype(str)
    )

    rows = [
        {
            "metric": "ref_marker_agreement",
            "value": float((ref == marker).mean()),
            "description": "Fraction of cells where ref.Rds primary label matches literature markers",
        },
        {
            "metric": "codi_marker_agreement",
            "value": float((codi[valid_codi] == marker[valid_codi]).mean()) if valid_codi.any() else np.nan,
            "description": "Fraction of CoDi-mapped cells with same label as literature markers",
        },
        {
            "metric": "codi_ref_agreement",
            "value": float((codi[valid_codi] == ref[valid_codi]).mean()) if valid_codi.any() else np.nan,
            "description": "Fraction of CoDi-mapped cells with same label as ref.Rds primary",
        },
        {
            "metric": "codi_mapped_fraction",
            "value": float(valid_codi.mean()),
            "description": "Fraction of cells with a CoDi reference label",
        },
    ]

    if has_azimuth and "azimuth_l1_norm" in adata.obs.columns:
        az_norm = adata.obs["azimuth_l1_norm"].astype(str)
        az_valid = az_norm != "NA"
        hi_conf = adata.obs["azimuth_score"] >= 0.5
        rows.extend(
            [
                {
                    "metric": "azimuth_marker_agreement",
                    "value": float((az_norm[az_valid] == marker[az_valid]).mean()),
                    "description": "Fraction of Azimuth-mapped cells matching marker annotation (L1 mapped)",
                },
                {
                    "metric": "azimuth_marker_agreement_score_ge_0.5",
                    "value": float((az_norm[hi_conf] == marker[hi_conf]).mean()),
                    "description": "Same as above, only cells with Azimuth score >= 0.5",
                },
                {
                    "metric": "codi_azimuth_agreement",
                    "value": float((codi[valid_codi & az_valid] == az_norm[valid_codi & az_valid]).mean())
                    if (valid_codi & az_valid).any()
                    else np.nan,
                    "description": "Fraction where CoDi normalized label matches Azimuth L1 mapped label",
                },
                {
                    "metric": "azimuth_mean_score",
                    "value": float(adata.obs["azimuth_score"].mean()),
                    "description": "Mean Azimuth prediction score across all cells",
                },
            ]
        )

        ctab = pd.crosstab(marker, adata.obs["azimuth_l1"], margins=True)
        ctab.to_csv(paths["results"] / "tables" / "annotation_crosstab_marker_azimuth.csv")
        log.log("  Saved: annotation_crosstab_marker_azimuth.csv")

        ctab_ref = pd.crosstab(ref, marker, margins=True)
        ctab_ref.to_csv(paths["results"] / "tables" / "annotation_crosstab_ref_marker.csv")
        log.log("  Saved: annotation_crosstab_ref_marker.csv")

    ctab_codi = _export_codi_marker_crosstabs(marker, codi, paths, log)

    if "cell_type_ref" in adata.obs.columns:
        ref_for_codi = adata.obs["cell_type_ref"].astype(str)
    elif "cell_type_seurat" in adata.obs.columns:
        ref_for_codi = adata.obs["cell_type_seurat"].astype(str)
    else:
        ref_for_codi = None

    if ref_for_codi is not None:
        _export_codi_ref_crosstabs(ref_for_codi, codi, paths, log)

    metrics = pd.DataFrame(rows)
    metrics.to_csv(paths["results"] / "tables" / "annotation_agreement_metrics.csv", index=False)
    log.log(f"  Saved: annotation_agreement_metrics.csv ({len(metrics)} metrics)")
    for _, r in metrics.iterrows():
        log.log(f"    {r['metric']}: {r['value']:.3f}")
    return metrics


def _interpret_module_delta(
    df: pd.DataFrame,
    score_col: str,
    *,
    group_col: str = "condition",
    control: str = "control",
    threshold: float = 0.01,
) -> List[str]:
    """Build short interpretation lines comparing each condition to control."""
    if score_col not in df.columns or control not in set(df[group_col]):
        return []

    ctrl_val = float(df.loc[df[group_col] == control, score_col].iloc[0])
    lines = []
    for cond in CONDITION_ORDER:
        if cond == control or cond not in set(df[group_col]):
            continue
        val = float(df.loc[df[group_col] == cond, score_col].iloc[0])
        delta = val - ctrl_val
        label = CONDITION_LABELS.get(cond, cond)
        if abs(delta) < threshold:
            verdict = "bez značajne promene"
        elif delta > 0:
            verdict = "povišen u odnosu na control"
        else:
            verdict = "snižen u odnosu na control"
        lines.append(f"  - {score_col} | {label}: {val:.4f} (Δ vs control {delta:+.4f}) — {verdict}")
    return lines


def write_additional_analyses_interpretation(
    module_tables: Dict[str, pd.DataFrame],
    agreement: pd.DataFrame,
    paths: Dict[str, Path],
    log: PipelineLogger,
) -> Path:
    log.log("  [5/6] Writing interpretation report (SR)...")
    out = paths["results"] / "tables" / "additional_analyses_interpretation_SR.md"
    lines = [
        "# Dodatne analize — interpretacija rezultata",
        "",
        "Ovaj fajl objašnjava šta znače tabele i figure iz `run_additional_analyses.py` / STEP 9 pipeline-a.",
        "",
        "## 1. Module score analize",
        "",
        "**Šta se radi:** Na punom normalizovanom transkriptomu (pre HVG filtriranja) računaju se",
        "skorovi genetskih programa pomoću `scanpy.tl.score_genes` / `score_genes_cell_cycle`.",
        "Svaka ćelija dobija jedan broj po programu; ovde se računaju proseci po uslovu i tipu ćelije.",
        "",
        "**Geni u panelima:**",
        "- S / G2M faza: klasični cell-cycle geni (MCM5, PCNA, MKI67, ...)",
        "- IFN: ISG15, IFIT1-3, MX1, OAS1, OASL, IFI6, RSAD2",
        "- Antigen presentation: HLA-DRA/DRB1, CD74, B2M, TAP1/2, HLA-DPA1/DPB1",
        "",
    ]

    by_cond = module_tables.get("by_condition")
    if by_cond is not None and not by_cond.empty:
        lines.append("### Globalni nalazi po uslovu")
        lines.append("")
        for col in ("S_score", "G2M_score", "IFN_score", "antigen_presentation_score"):
            if col in by_cond.columns:
                lines.append(f"**{col}**")
                lines.extend(_interpret_module_delta(by_cond, col))
                lines.append("")

    ag = module_tables.get("antigen")
    if ag is not None and not ag.empty:
        lines.extend(
            [
                "### Antigen presentation po tipu ćelije",
                "",
                "Najviši skorovi su očekivani kod B ćelija i antigen-prezentujućih tipova (DC, monociti).",
                "Poređenje uslova unutar istog tipa pokazuje da li izloženost menja MHC/HLA program.",
                "",
            ]
        )
        for ctype in ["B_cell", "Monocyte_CD14", "DC"]:
            sub = ag[ag["cell_type_marker"] == ctype]
            if sub.empty:
                continue
            lines.append(f"**{ctype.replace('_', ' ')}**")
            ctrl = sub[sub["condition"] == "control"]["antigen_presentation_score"]
            if ctrl.empty:
                continue
            ctrl_v = float(ctrl.iloc[0])
            for cond in ["PSNP_40nm", "PSNP_200nm", "PSNP_mix_40_200"]:
                row = sub[sub["condition"] == cond]
                if row.empty:
                    continue
                v = float(row["antigen_presentation_score"].iloc[0])
                d = v - ctrl_v
                lines.append(
                    f"  - {CONDITION_LABELS.get(cond, cond)}: {v:.3f} (Δ {d:+.3f} vs control)"
                )
            lines.append("")

    lines.extend(
        [
            "## 2. Pseudobulk matrica",
            "",
            "**Šta se radi:** Za svaku kombinaciju `condition × cell_type` sabiraju se sirovi UMI",
            "brojevi po genu (iz `layers['counts']`).",
            "",
        "**Šta znači:** Omogućava bulk RNA-seq stil analizu (DESeq2, edgeR) bez ponovnog",
        "učitavanja pojedinačnih ćelija. Broj grupa = broj uslova × broj tipova ćelija.",
            "",
            "## 3. Validacija anotacije",
            "",
            "**Šta se radi:** Porede se tri nezavisna izvora tipova ćelija:",
            "marker panel (pipeline), CoDi (Zenodo CSV), Azimuth PBMC referenca (opciono).",
            "",
        ]
    )

    if agreement is not None and not agreement.empty:
        lines.append("| Metrika | Vrednost | Značenje |")
        lines.append("|---------|----------|----------|")
        for _, r in agreement.iterrows():
            val = r["value"]
            val_s = f"{val:.3f}" if pd.notna(val) else "n/a"
            lines.append(f"| {r['metric']} | {val_s} | {r['description']} |")
        lines.append("")
        lines.extend(
            [
                "**Kako čitati:**",
                "- 100% slaganje je retko — NK/CD8 i granularnost tipova razlikuju metode.",
                "- ~45–65% je tipično za PBMC sa više anotacionih šema.",
                "- Azimuth mean score > 0.85 = visoko poverenje u referentno mapiranje.",
                "",
                "**Contingency matrix (CoDi vs markeri):**",
                "- `annotation_crosstab_marker_codi.csv` — broj ćelija po paru tipova",
                "- `annotation_crosstab_marker_codi_row_pct.csv` — % unutar svakog marker tipa",
                "- `annotation_codi_marker_mapping.csv` — long format (marker → CoDi mapiranje)",
                "",
                "**Contingency matrix (ref.Rds/Seurat vs CoDi):**",
                "- `annotation_crosstab_ref_codi.csv` — broj ćelija po paru tipova",
                "- `annotation_crosstab_ref_codi_row_pct.csv` — % unutar svakog ref.Rds tipa",
                "- `annotation_codi_ref_mapping.csv` — long format (ref.Rds → CoDi mapiranje)",
                "",
            ]
        )

    lines.extend(
        [
            "## 4. Figure",
            "",
            "- `results/figures/additional_analyses/module_scores_by_condition.png` — bar chart modula",
            "- `results/figures/additional_analyses/module_scores_violin.png` — raspodela po uslovu",
            "- `results/figures/additional_analyses/antigen_presentation_heatmap.png`",
            "- `results/figures/additional_analyses/annotation_agreement_bar.png`",
            "- `results/figures/additional_analyses/annotation_confusion_marker_codi.png` — CoDi vs markeri",
            "- `results/figures/additional_analyses/annotation_confusion_marker_codi_normalized.png`",
            "- `results/figures/additional_analyses/annotation_confusion_ref_codi.png` — ref.Rds vs CoDi",
            "- `results/figures/additional_analyses/annotation_confusion_ref_codi_normalized.png`",
            "- `results/figures/additional_analyses/annotation_confusion_marker_azimuth.png` (ako postoji Azimuth)",
            "",
        ]
    )

    out.write_text("\n".join(lines), encoding="utf-8")
    log.log(f"  Saved: {out}")
    return out


def plot_additional_analysis_figures(
    adata: sc.AnnData,
    module_tables: Dict[str, pd.DataFrame],
    agreement: pd.DataFrame,
    paths: Dict[str, Path],
    log: PipelineLogger,
) -> List[Path]:
    log.log("  [6/6] Plotting additional analysis figures...")
    fig_dir = paths["results"] / "figures" / "additional_analyses"
    fig_dir.mkdir(parents=True, exist_ok=True)
    saved: List[Path] = []

    by_cond = module_tables.get("by_condition")
    if by_cond is not None and not by_cond.empty:
        plot_df = by_cond.melt(
            id_vars=["condition"], var_name="module", value_name="score"
        )
        plot_df["condition"] = plot_df["condition"].astype(str)
        plot_df["condition_label"] = plot_df["condition"].map(CONDITION_LABELS)
        plot_df["condition_label"] = plot_df["condition_label"].fillna(plot_df["condition"])
        order = [CONDITION_LABELS[c] for c in CONDITION_ORDER if c in set(by_cond["condition"])]

        fig, ax = plt.subplots(figsize=(8, 4.5))
        sns.barplot(
            data=plot_df,
            x="module",
            y="score",
            hue="condition_label",
            hue_order=order,
            ax=ax,
            palette=[CONDITION_COLORS[c] for c in CONDITION_ORDER],
        )
        ax.axhline(0, color="#888888", linewidth=0.8)
        ax.set_title("Gene module scores by condition")
        ax.set_xlabel("Module")
        ax.set_ylabel("Mean score")
        ax.legend(title="Condition", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
        plt.tight_layout()
        p = fig_dir / "module_scores_by_condition.png"
        fig.savefig(p, dpi=200, bbox_inches="tight")
        plt.close(fig)
        saved.append(p)

        score_cols = [c for c in MODULE_SCORE_COLUMNS if c in adata.obs.columns]
        if score_cols:
            long_obs = adata.obs[["condition"] + list(score_cols)].copy()
            long_obs["condition"] = long_obs["condition"].astype(str)
            long_obs = long_obs.melt(
                id_vars=["condition"], var_name="module", value_name="score"
            )
            long_obs["condition_label"] = long_obs["condition"].map(CONDITION_LABELS)
            long_obs["condition_label"] = long_obs["condition_label"].fillna(long_obs["condition"])
            fig, axes = plt.subplots(1, len(score_cols), figsize=(4 * len(score_cols), 4))
            if len(score_cols) == 1:
                axes = [axes]
            for ax, mod in zip(axes, score_cols):
                sns.violinplot(
                    data=long_obs[long_obs["module"] == mod],
                    x="condition_label",
                    y="score",
                    order=order,
                    ax=ax,
                    palette=[CONDITION_COLORS[c] for c in CONDITION_ORDER],
                    cut=0,
                    inner="quartile",
                )
                ax.set_title(mod)
                ax.tick_params(axis="x", rotation=25)
                ax.set_xlabel("")
            fig.suptitle("Module score distributions per condition", y=1.02)
            plt.tight_layout()
            p = fig_dir / "module_scores_violin.png"
            fig.savefig(p, dpi=200, bbox_inches="tight")
            plt.close(fig)
            saved.append(p)

    ag = module_tables.get("antigen")
    if ag is not None and not ag.empty:
        pivot = ag.pivot_table(
            index="cell_type_marker",
            columns="condition",
            values="antigen_presentation_score",
            aggfunc="first",
        )
        pivot = pivot.reindex(columns=[c for c in CONDITION_ORDER if c in pivot.columns])
        pivot.index = [str(i).replace("_", " ") for i in pivot.index]
        pivot.columns = [CONDITION_LABELS.get(c, c) for c in pivot.columns]
        fig_h = max(4.5, 0.35 * len(pivot) + 1.5)
        fig, ax = plt.subplots(figsize=(6.5, fig_h))
        sns.heatmap(
            pivot,
            ax=ax,
            cmap="viridis",
            linewidths=0.4,
            linecolor="white",
            cbar_kws={"label": "Mean antigen presentation score", "shrink": 0.7},
        )
        ax.set_title("Antigen presentation module by cell type and condition")
        plt.tight_layout()
        p = fig_dir / "antigen_presentation_heatmap.png"
        fig.savefig(p, dpi=200, bbox_inches="tight")
        plt.close(fig)
        saved.append(p)

    if agreement is not None and not agreement.empty:
        plot_ag = agreement[agreement["metric"].str.contains("agreement|mean_score")].copy()
        fig, ax = plt.subplots(figsize=(7, 4))
        sns.barplot(data=plot_ag, x="metric", y="value", ax=ax, color="#4C72B0")
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("Fraction / mean score")
        ax.set_title("Annotation agreement metrics")
        ax.tick_params(axis="x", rotation=30, labelsize=8)
        plt.tight_layout()
        p = fig_dir / "annotation_agreement_bar.png"
        fig.savefig(p, dpi=200, bbox_inches="tight")
        plt.close(fig)
        saved.append(p)

    codi = adata.obs.get("cell_type_codi_norm", pd.Series("NA", index=adata.obs_names)).astype(str)
    valid_codi = codi != "NA"
    if valid_codi.any():
        ctab_codi = pd.crosstab(
            adata.obs.loc[valid_codi, "cell_type_marker"].astype(str),
            codi[valid_codi],
        )
        if not ctab_codi.empty:
            _plot_annotation_crosstab_heatmap(
                ctab_codi,
                title="Contingency matrix: marker genes vs CoDi annotation",
                xlabel="CoDi cell type",
                ylabel="Marker-based cell type",
                out_path=fig_dir / "annotation_confusion_marker_codi.png",
            )
            saved.append(fig_dir / "annotation_confusion_marker_codi.png")
            _plot_annotation_crosstab_heatmap(
                ctab_codi,
                title="Marker vs CoDi (% of each marker-based type)",
                xlabel="CoDi cell type",
                ylabel="Marker-based cell type",
                out_path=fig_dir / "annotation_confusion_marker_codi_normalized.png",
                normalize_rows=True,
            )
            saved.append(fig_dir / "annotation_confusion_marker_codi_normalized.png")

    ref_col: Optional[str] = None
    if "cell_type_ref" in adata.obs.columns:
        ref_col = "cell_type_ref"
    elif "cell_type_seurat" in adata.obs.columns:
        ref_col = "cell_type_seurat"

    if ref_col is not None and valid_codi.any():
        ref_s = adata.obs[ref_col].astype(str)
        plot_mask = valid_codi & (ref_s != "unmapped") & (ref_s != "NA")
        if plot_mask.any():
            ctab_ref_codi = pd.crosstab(ref_s[plot_mask], codi[plot_mask])
            if not ctab_ref_codi.empty:
                _plot_annotation_crosstab_heatmap(
                    ctab_ref_codi,
                    title="Contingency matrix: ref.Rds/Seurat vs CoDi annotation",
                    xlabel="CoDi cell type",
                    ylabel="Seurat/ref.Rds cell type",
                    out_path=fig_dir / "annotation_confusion_ref_codi.png",
                )
                saved.append(fig_dir / "annotation_confusion_ref_codi.png")
                _plot_annotation_crosstab_heatmap(
                    ctab_ref_codi,
                    title="ref.Rds/Seurat vs CoDi (% of each ref.Rds type)",
                    xlabel="CoDi cell type",
                    ylabel="Seurat/ref.Rds cell type",
                    out_path=fig_dir / "annotation_confusion_ref_codi_normalized.png",
                    normalize_rows=True,
                )
                saved.append(fig_dir / "annotation_confusion_ref_codi_normalized.png")

    if "azimuth_l1" in adata.obs.columns:
        ct = pd.crosstab(
            adata.obs["cell_type_marker"].astype(str),
            adata.obs["azimuth_l1"].astype(str),
        )
        fig, ax = plt.subplots(figsize=(9, 6))
        sns.heatmap(ct, ax=ax, cmap="Blues", linewidths=0.3, cbar_kws={"label": "Cell count"})
        ax.set_title("Marker annotation vs Azimuth L1 labels")
        ax.set_xlabel("Azimuth L1")
        ax.set_ylabel("Marker-based type")
        plt.tight_layout()
        p = fig_dir / "annotation_confusion_marker_azimuth.png"
        fig.savefig(p, dpi=200, bbox_inches="tight")
        plt.close(fig)
        saved.append(p)

    log.log(f"  Saved {len(saved)} figure(s) under {fig_dir}")
    return saved


def additional_insights(
    adata: sc.AnnData,
    paths: Dict[str, Path],
    log: PipelineLogger,
    *,
    plot_figures: bool = True,
    write_interpretation: bool = True,
) -> None:
    """
    Additional analyses (STEP 9): module scores, pseudobulk, annotation validation.

    Can also be run standalone via scripts/run_additional_analyses.py.
    Requires module scores in adata.obs (computed before HVG subset in merge_and_integrate).
    """
    cfg = load_config()
    has_azimuth = attach_azimuth_labels(adata, paths, log, cfg)
    module_tables = _export_module_score_tables(adata, paths, log)
    _export_pseudobulk(adata, paths, log)
    agreement = _export_annotation_agreement(adata, paths, log, has_azimuth)

    if write_interpretation:
        write_additional_analyses_interpretation(module_tables, agreement, paths, log)
    if plot_figures:
        plot_additional_analysis_figures(adata, module_tables, agreement, paths, log)


def load_integrated_adata(paths: Dict[str, Path], log: PipelineLogger) -> sc.AnnData:
    """Load processed object for standalone additional-analysis runs."""
    h5ad = paths["processed"] / "integrated_annotated.h5ad"
    if not h5ad.exists():
        raise FileNotFoundError(
            f"Missing {h5ad}. Run first: python scripts/run_pipeline.py"
        )
    log.log(f"  Loading {h5ad} ...")
    adata = sc.read_h5ad(h5ad)
    log.log(f"  Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    missing_scores = [c for c in MODULE_SCORE_COLUMNS if c not in adata.obs.columns]
    if missing_scores:
        log.log(
            f"  Warning: missing module scores {missing_scores}. "
            "Re-run full pipeline merge step or scores will be incomplete."
        )

    if "cell_type_codi_norm" not in adata.obs.columns:
        log.log("  CoDi labels missing from object — reloading from CSV...")
        load_codi_annotations(adata, paths, log)

    if "cell_type_primary" not in adata.obs.columns:
        log.log("  Primary annotation missing — recomputing from config + annotation CSVs...")
        cfg = load_config()
        if "cell_type_marker" not in adata.obs.columns:
            marker_dict = resolve_marker_dict(cfg, paths, log)
            marker_based_annotation(adata, marker_dict, log)
        load_seurat_annotations(adata, paths, cfg, log)
        run_python_reference_methods(adata, paths, cfg, log)
        assign_primary_cell_type(adata, cfg, log)

    return adata


def _umap_point_size(n_obs: int) -> float:
    if n_obs > 25000:
        return 6.0
    if n_obs > 15000:
        return 10.0
    return 14.0


def _save_scanpy_umap(
    adata: sc.AnnData,
    path: Path,
    *,
    color,
    title: str = "",
    palette=None,
    cmap=None,
    legend_loc: str = "right margin",
    ncols: int = 1,
    vmin=None,
    vmax=None,
) -> None:
    pt_size = _umap_point_size(adata.n_obs)
    kwargs = {
        "color": color,
        "show": False,
        "size": pt_size,
        "alpha": 0.75,
        "frameon": False,
        "legend_loc": legend_loc,
        "title": title,
        "ncols": ncols,
        "wspace": 0.35,
    }
    if palette is not None:
        kwargs["palette"] = palette
    if cmap is not None:
        kwargs["cmap"] = cmap
    if vmin is not None:
        kwargs["vmin"] = vmin
    if vmax is not None:
        kwargs["vmax"] = vmax

    width = 5.5 * ncols + 1.5
    height = 5.0 if ncols == 1 else 4.5
    sc.pl.umap(adata, **kwargs)
    fig = plt.gcf()
    fig.set_size_inches(width, height)
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _save_umap_split_by_condition(adata: sc.AnnData, path: Path) -> None:
    umap = adata.obsm["X_umap"]
    conditions = [c for c in CONDITION_COLORS if c in set(adata.obs["condition"].astype(str))]
    if not conditions:
        conditions = sorted(adata.obs["condition"].astype(str).unique())

    fig, axes = plt.subplots(2, 2, figsize=(13, 11))
    for ax, cond in zip(axes.flat, conditions):
        mask = adata.obs["condition"].astype(str) == cond
        ax.scatter(
            umap[:, 0],
            umap[:, 1],
            c="#d9d9d9",
            s=1.5,
            alpha=0.25,
            linewidths=0,
            rasterized=True,
        )
        ax.scatter(
            umap[mask, 0],
            umap[mask, 1],
            c=CONDITION_COLORS.get(cond, "#333333"),
            s=4.0,
            alpha=0.85,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(CONDITION_LABELS.get(cond, cond), fontsize=12, fontweight="bold")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle("UMAP split by experimental condition", fontsize=14, fontweight="bold", y=0.98)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_core_figures(
    adata: sc.AnnData,
    paths: Dict[str, Path],
    log: PipelineLogger,
    cfg: Optional[Dict] = None,
    marker_dict: Optional[Dict[str, List[str]]] = None,
) -> List[str]:
    saved = []
    fig_dir = paths["results"] / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    condition_palette = {
        cond: CONDITION_COLORS.get(cond, "#888888")
        for cond in sorted(adata.obs["condition"].astype(str).unique())
    }

    f1 = fig_dir / "umap_condition.png"
    _save_scanpy_umap(
        adata,
        f1,
        color=["condition"],
        title="Experimental condition (post-Combat integration)",
        palette=condition_palette,
    )
    saved.append(f1.name)
    log.log(f"  Saved: {f1}")

    f2 = fig_dir / "umap_clusters.png"
    _save_scanpy_umap(adata, f2, color=["cluster"], title="Leiden clusters")
    saved.append(f2.name)
    log.log(f"  Saved: {f2}")

    f3 = fig_dir / "umap_celltypes_marker.png"
    _save_scanpy_umap(
        adata,
        f3,
        color=["cell_type_marker"],
        title="Literaturini markeri (validacija)",
    )
    saved.append(f3.name)
    log.log(f"  Saved: {f3}")

    if "cell_type_seurat" in adata.obs.columns:
        f3b = fig_dir / "umap_celltypes_seurat.png"
        _save_scanpy_umap(
            adata,
            f3b,
            color=["cell_type_seurat"],
            title="Seurat/Azimuth anotacija (ref.Rds, primarna)",
        )
        saved.append(f3b.name)
        log.log(f"  Saved: {f3b}")

    f4 = fig_dir / "umap_split_by_condition.png"
    _save_umap_split_by_condition(adata, f4)
    saved.append(f4.name)
    log.log(f"  Saved: {f4}")

    f5 = fig_dir / "umap_sample_integration.png"
    _save_scanpy_umap(
        adata,
        f5,
        color=["sample_id"],
        title="Sample / batch integration check",
    )
    saved.append(f5.name)
    log.log(f"  Saved: {f5}")

    score_cols = [
        c
        for c in ("S_score", "G2M_score", "IFN_score")
        if c in adata.obs.columns and adata.obs[c].notna().any()
    ]
    if score_cols:
        f6 = fig_dir / "umap_module_scores.png"
        _save_scanpy_umap(
            adata,
            f6,
            color=score_cols,
            title="",
            cmap="RdYlBu_r",
            legend_loc="right margin",
            ncols=len(score_cols),
            vmin="p5",
            vmax="p95",
        )
        saved.append(f6.name)
        log.log(f"  Saved: {f6}")

    codi_mapped = (
        "cell_type_codi" in adata.obs.columns
        and (adata.obs["cell_type_codi"].astype(str) != "NA").mean() > 0.05
    )
    if codi_mapped:
        f7 = fig_dir / "umap_codi_celltypes.png"
        _save_scanpy_umap(
            adata,
            f7,
            color=["cell_type_codi"],
            title="CoDi reference cell types",
        )
        saved.append(f7.name)
        log.log(f"  Saved: {f7}")

    if marker_dict or cfg:
        panels = marker_dict or (resolve_marker_dict(cfg, paths, log) if cfg else {})
        marker_dict_filtered = {
            ct: [g for g in genes if g in adata.var_names]
            for ct, genes in panels.items()
        }
        marker_dict_filtered = {ct: genes for ct, genes in marker_dict_filtered.items() if genes}
        if marker_dict_filtered:
            f8 = fig_dir / "marker_dotplot.png"
            sc.pl.dotplot(
                adata,
                var_names=marker_dict_filtered,
                groupby="cell_type_marker",
                standard_scale="var",
                show=False,
                dendrogram=False,
            )
            plt.gcf().suptitle("Literaturini markeri po marker-anotaciji (validacija)", y=1.02)
            plt.gcf().set_size_inches(12, 5)
            plt.savefig(f8, dpi=300, bbox_inches="tight")
            plt.close()
            saved.append(f8.name)
            log.log(f"  Saved: {f8}")

    return saved


def log_output_inventory(paths: Dict[str, Path], log: PipelineLogger) -> None:
    log.section("OUTPUT FILES - TABLES")
    tables = sorted((paths["results"] / "tables").glob("*.csv"))
    if tables:
        for t in tables:
            log.log(f"  {t.name}")
    else:
        log.log("  (no CSV tables found)")

    log.section("OUTPUT FILES - FIGURES")
    figures = sorted((paths["results"] / "figures").glob("*.png"))
    for fig in figures:
        log.log(f"  {fig.name}")


def explain_figures(log: PipelineLogger) -> None:
    log.section("FIGURE GUIDE (for slides)")
    log.log("Use these captions when building your PowerPoint:")
    log.log("")
    for name, text in FIGURE_DESCRIPTIONS.items():
        log.log(f"  {name}")
        log.log(f"    {text}")
        log.log("")


def explain_tables(log: PipelineLogger) -> None:
    log.section("TABLE GUIDE (for slides)")
    guides = {
        "cell_composition_by_condition.csv": (
            "Cell counts and fractions per condition and marker-based cell type."
        ),
        "differential_expression_all.csv": (
            "Wilcoxon DE results: each exposure vs control, stratified by cell type."
        ),
        "pathway_enrichment_all.csv": (
            "Enrichr hits (GO/KEGG/Reactome) for significant UP and DOWN DE genes; "
            "includes direction column. Figures in results/figures/pathway_enrichment/."
        ),
        "pathway_enrichment_summary.csv": (
            "Top enriched pathways per cell type x exposure x direction (UP/DOWN vs control) "
            "with biological interpretation notes."
        ),
        "size_specific_effects_summary.csv": (
            "Counts of DE genes unique to 40 nm, 200 nm, shared, or mix-only effects."
        ),
        "size_specific_genes.csv": (
            "Gene lists per cell type and size-effect class with logFC from the primary comparison."
        ),
        "size_specific_pathway_enrichment.csv": (
            "Enrichr hits (GO/KEGG/Reactome) for UP genes in each size-effect class."
        ),
        "size_specific_interpretation.csv": (
            "Top genes, top pathways, and short biological notes per size-effect class."
        ),
        "cell_cycle_scores_by_condition.csv": (
            "Mean S and G2M phase scores per condition (subset of module_scores_by_condition)."
        ),
        "ifn_scores_by_condition.csv": (
            "Mean interferon module score per condition."
        ),
        "module_scores_by_condition.csv": (
            "All module scores (S, G2M, IFN, antigen presentation) averaged per condition."
        ),
        "module_scores_by_condition_celltype.csv": (
            "Module scores per condition and marker-based cell type."
        ),
        "antigen_presentation_scores.csv": (
            "Mean HLA/MHC module score by condition and cell type."
        ),
        "pseudobulk_counts_condition_celltype.csv": (
            "Summed UMI counts per condition x cell type for bulk-style follow-up."
        ),
        "annotation_agreement_metrics.csv": (
            "CoDi/marker/Azimuth agreement fractions and mean Azimuth score."
        ),
        "annotation_crosstab_marker_codi.csv": (
            "Contingency matrix: cell counts for marker-based type vs normalized CoDi label."
        ),
        "annotation_crosstab_marker_codi_row_pct.csv": (
            "Row-normalized CoDi vs marker matrix (% of each marker type)."
        ),
        "annotation_codi_marker_mapping.csv": (
            "Long-format CoDi vs marker mapping with counts and row percentages."
        ),
        "annotation_crosstab_ref_codi.csv": (
            "Contingency matrix: Seurat/ref.Rds cell type vs normalized CoDi label."
        ),
        "annotation_crosstab_ref_codi_row_pct.csv": (
            "Row-normalized ref.Rds vs CoDi matrix (% of each ref.Rds type)."
        ),
        "annotation_codi_ref_mapping.csv": (
            "Long-format ref.Rds vs CoDi mapping with counts and row percentages."
        ),
        "annotation_crosstab_marker_azimuth.csv": (
            "Cell counts: marker-based type vs Azimuth L1 label."
        ),
        "additional_analyses_interpretation_SR.md": (
            "Serbian interpretation of additional analyses tables and figures."
        ),
    }
    for name, text in guides.items():
        log.log(f"  {name}: {text}")


def main():
    sc.settings.verbosity = 1
    sc.set_figure_params(dpi=100, facecolor="white")

    cfg = load_config()
    paths = setup_paths(cfg)
    log = PipelineLogger(paths["run_logs"])

    try:
        log.section("CONFIGURATION")
        log.log(f"  Project: {cfg.get('project_name', 'scRNA-seq pipeline')}")
        log.log(f"  Integration: {cfg['preprocessing'].get('integration_method', 'combat')}")
        log.log(f"  Raw data: {paths['raw'].resolve()}")
        log.log(f"  Results: {paths['results'].resolve()}")

        describe_raw_data(cfg, paths, log)

        log.section("STEP 1 - QUALITY CONTROL (per sample)")
        adata = None
        for sample_id, meta in cfg["samples"].items():
            f = paths["raw"] / f"{sample_id}.h5ad"
            if not f.exists():
                raise FileNotFoundError(f"Missing file: {f}. Run: python scripts/download_data.py")
            sample = read_and_qc_sample(f, sample_id, meta["condition"], cfg["qc"], log)
            if adata is None:
                adata = sample
                continue
            log.log(
                f"  Merging {sample_id} into integrated object "
                f"({adata.n_obs:,} + {sample.n_obs:,} cells, inner join on genes)..."
            )
            prev = adata
            adata = sc.concat([prev, sample], join="inner", merge="same")
            adata.obs_names_make_unique()
            adata = _compact_adata_for_merge(adata)
            del prev, sample
            gc.collect()

        log.section("STEP 2 - MERGE, INTEGRATION, UMAP, CLUSTERING")
        log.log(
            f"  All samples merged: {adata.n_obs:,} cells x {adata.n_vars:,} genes "
            "(incremental concat, inner join)"
        )
        adata = merge_and_integrate(adata, cfg, paths, log)

        log.section("STEP 3 - CELL-TYPE ANNOTATION (ref.Rds reference atlas)")
        log.log("  Methods: Seurat, CoDi, Tangram, cell2location | Validation: literature markers")
        marker_dict = resolve_marker_dict(cfg, paths, log)
        marker_based_annotation(adata, marker_dict, log)
        load_codi_annotations(adata, paths, log)
        load_seurat_annotations(adata, paths, cfg, log)
        run_python_reference_methods(adata, paths, cfg, log)
        assign_primary_cell_type(adata, cfg, log)
        validate_with_literature_markers(adata, marker_dict, paths, log)
        export_method_agreement(adata, paths, log)

        log.section("STEP 4 - CORE UMAP FIGURES")
        save_core_figures(adata, paths, log, cfg, marker_dict=marker_dict)

        log.section("STEP 5 - CELL COMPOSITION")
        log.log("  Comparing immune cell proportions across PSNP conditions vs control...")
        composition_analysis(adata, paths, log)

        log.section("STEP 6 - DIFFERENTIAL EXPRESSION")
        log.log(
            f"  Wilcoxon test per cell type; min {cfg['de']['min_cells_per_group']} cells per group..."
        )
        de_all = differential_expression_by_celltype(adata, cfg, paths, log)

        log.section("STEP 7 - PATHWAY ENRICHMENT")
        enr_all = pathway_enrichment(de_all, cfg, paths, log)
        plot_pathway_enrichment_figures(enr_all, cfg, paths, log)
        pathway_summary = export_pathway_enrichment_summary(enr_all, cfg, paths, log)
        write_pathway_enrichment_interpretation(pathway_summary, paths, log, cfg=cfg)

        log.section("STEP 8 - SIZE-SPECIFIC EFFECTS")
        size_specific_effects(de_all, cfg, paths, log)

        log.section("STEP 9 - ADDITIONAL ANALYSES")
        additional_insights(adata, paths, log)

        log.section("STEP 10 - SAVE PROCESSED OBJECT")
        out_h5ad = paths["processed"] / "integrated_annotated.h5ad"
        adata.write(out_h5ad)
        log.log(f"  Saved integrated AnnData: {out_h5ad}")
        log.log(f"  Final object: {adata.n_obs:,} cells x {adata.n_vars:,} genes (HVG space)")

        log_output_inventory(paths, log)
        explain_figures(log)
        explain_tables(log)

        log.section("PIPELINE COMPLETE")
        log.log("  All steps finished successfully.")
        log.log(f"  Full run log copied to: {log.log_path.resolve()}")
        log.log("  ref.Rds annotation workflow:")
        log.log("    1) python scripts/download_data.py          # ref.Rds + idx.annoy")
        log.log("    2) python scripts/run_pipeline.py           # creates integrated_pre_hvg.h5ad")
        log.log("    3) python scripts/prepare_azimuth_h5ad.py")
        log.log("    4) Rscript scripts/azimuth_annotation.R    # Seurat + ref.Rds")
        log.log("    5) python scripts/run_reference_annotation.py  # Tangram + cell2location")
        log.log("    6) python scripts/run_pipeline.py           # full analysis with all labels")

    except Exception as exc:
        log.section("PIPELINE FAILED")
        log.log(f"  Error: {exc}")
        raise
    finally:
        log.close()


if __name__ == "__main__":
    main()
