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
from typing import Dict, List, TextIO

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import yaml
from scipy import sparse


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
    "composition_barplot.png": (
        "Stacked bar chart of cell-type fractions per condition. "
        "Answers whether nanoplastic exposure changes immune cell proportions "
        "(e.g. more monocytes, fewer T cells) relative to control."
    ),
}


def load_config(path: str = "config/config.yaml") -> Dict:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_paths(cfg: Dict) -> Dict[str, Path]:
    paths = {
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


def merge_and_integrate(adata: sc.AnnData, cfg: Dict, log: PipelineLogger) -> sc.AnnData:
    adata.obs_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    log.log(f"  Merged object: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    log.log("  Normalizing to 10,000 counts per cell and log1p transform...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

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
    log.log("  Cell-type counts (marker-based):")
    for ct, n in counts.items():
        log.log(f"    {ct}: {n:,}")


def load_codi_annotations(
    paths: Dict[str, Path], adata: sc.AnnData, log: PipelineLogger
) -> None:
    log.log("  Loading external CoDi reference labels from *_CoDi_KLD.csv files...")
    codi_files = sorted(paths["raw"].glob("*_CoDi_KLD.csv"))
    codi_frames = []
    for f in codi_files:
        df = pd.read_csv(f)
        if "cell_id" in df.columns and "CoDi" in df.columns:
            df = df[["cell_id", "CoDi"]].copy()
            df.columns = ["cell_id", "cell_type_codi"]
            codi_frames.append(df)
    if not codi_frames:
        adata.obs["cell_type_codi"] = "NA"
        log.log("  Warning: no CoDi CSV files found - cell_type_codi set to NA")
        return

    codi_all = pd.concat(codi_frames, ignore_index=True).drop_duplicates(subset=["cell_id"])
    codi_all = codi_all.assign(cell_id_clean=codi_all["cell_id"].str.replace(r"-\d+$", "", regex=True))
    codi_all = codi_all.set_index("cell_id_clean")
    adata.obs["cell_id_clean"] = adata.obs_names.str.replace(r"-\d+$", "", regex=True)
    adata.obs["cell_type_codi"] = adata.obs["cell_id_clean"].map(codi_all["cell_type_codi"]).fillna("NA")
    mapped = (adata.obs["cell_type_codi"] != "NA").mean() * 100
    log.log(f"  CoDi labels mapped to {mapped:.1f}% of cells")

    codi_to_marker = {
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
    adata.obs["cell_type_codi_norm"] = adata.obs["cell_type_codi"].map(codi_to_marker).fillna("NA")


def composition_analysis(
    adata: sc.AnnData, paths: Dict[str, Path], log: PipelineLogger
) -> pd.DataFrame:
    comp = (
        adata.obs.groupby(["condition", "cell_type_marker"])
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
    sns.barplot(data=comp, x="condition", y="fraction", hue="cell_type_marker")
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
    for ctype in sorted(adata.obs["cell_type_marker"].unique()):
        ad_ct = adata[adata.obs["cell_type_marker"] == ctype].copy()
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


def pathway_enrichment(
    de_all: pd.DataFrame, cfg: Dict, paths: Dict[str, Path], log: PipelineLogger
) -> pd.DataFrame:
    if de_all.empty:
        log.log("  Skipped pathway enrichment (no DE results)")
        return pd.DataFrame()

    log.log("  Running Enrichr (GO, KEGG, Reactome) on significant DE genes per cell type x comparison...")
    enr_frames = []
    for (ctype, comp), df in de_all.groupby(["cell_type", "comparison"]):
        sig = df[
            (df["pvals_adj"] < cfg["de"]["pval_adj_threshold"])
            & (df["logfoldchanges"] > cfg["de"]["logfc_threshold"])
        ]["names"].dropna().astype(str).unique()

        if len(sig) < 10:
            continue

        gene_sets = ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]
        for gset in gene_sets:
            try:
                enr = gp.enrichr(gene_list=list(sig), gene_sets=gset, organism="human", outdir=None)
                if enr.results is None or enr.results.empty:
                    continue
                tmp = enr.results.copy()
                tmp["cell_type"] = ctype
                tmp["comparison"] = comp
                tmp["gene_set"] = gset
                enr_frames.append(tmp)
                log.log(f"    Enriched: {ctype} / {comp} / {gset} ({len(sig)} input genes)")
            except Exception as exc:
                log.log(f"    Enrichr failed for {ctype}/{comp}/{gset}: {exc}")

    if not enr_frames:
        log.log("  No pathway enrichment results returned")
        return pd.DataFrame()
    enr_all = pd.concat(enr_frames, ignore_index=True)
    out_csv = paths["results"] / "tables" / "pathway_enrichment_all.csv"
    enr_all.to_csv(out_csv, index=False)
    log.log(f"  Saved: {out_csv} ({len(enr_all):,} rows)")
    return enr_all


def size_specific_effects(
    de_all: pd.DataFrame, cfg: Dict, paths: Dict[str, Path], log: PipelineLogger
) -> pd.DataFrame:
    if de_all.empty:
        log.log("  Skipped size-specific summary (no DE)")
        return pd.DataFrame()

    log.log("  Classifying DE genes by particle size (unique 40nm / 200nm / shared / mix-only)...")
    sig = de_all[
        (de_all["pvals_adj"] < cfg["de"]["pval_adj_threshold"])
        & (de_all["logfoldchanges"].abs() > cfg["de"]["logfc_threshold"])
    ][["cell_type", "comparison", "names"]].drop_duplicates()

    results = []
    for ctype in sig["cell_type"].unique():
        s = sig[sig["cell_type"] == ctype]
        s40 = set(s[s["comparison"] == "PSNP_40nm_vs_control"]["names"])
        s200 = set(s[s["comparison"] == "PSNP_200nm_vs_control"]["names"])
        smix = set(s[s["comparison"] == "PSNP_mix_40_200_vs_control"]["names"])

        rows = {
            "unique_40nm": len(s40 - s200 - smix),
            "unique_200nm": len(s200 - s40 - smix),
            "shared_40_200": len((s40 & s200) - smix),
            "shared_all_three": len(s40 & s200 & smix),
            "mix_only_emergent": len(smix - s40 - s200),
        }
        for k, v in rows.items():
            results.append({"cell_type": ctype, "effect_class": k, "n_genes": v})
        log.log(f"    {ctype}: " + ", ".join(f"{k}={v}" for k, v in rows.items()))

    out = pd.DataFrame(results)
    out_csv = paths["results"] / "tables" / "size_specific_effects_summary.csv"
    out.to_csv(out_csv, index=False)
    log.log(f"  Saved: {out_csv}")
    return out


def additional_insights(
    adata: sc.AnnData, paths: Dict[str, Path], log: PipelineLogger
) -> None:
    log.log("  [1/5] Cell-cycle scores (S and G2/M phases)...")
    s_genes = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4"]
    g2m_genes = ["HMGB2", "CDK1", "NUSAP1", "TOP2A", "MKI67", "BIRC5"]
    valid_s = [g for g in s_genes if g in adata.var_names]
    valid_g2m = [g for g in g2m_genes if g in adata.var_names]
    if valid_s and valid_g2m:
        sc.tl.score_genes_cell_cycle(adata, s_genes=valid_s, g2m_genes=valid_g2m)
    else:
        log.log("  Warning: cell-cycle genes missing - scores set to NaN")
        adata.obs["S_score"] = np.nan
        adata.obs["G2M_score"] = np.nan
    p = paths["results"] / "tables" / "cell_cycle_scores_by_condition.csv"
    adata.obs.groupby("condition")[["S_score", "G2M_score"]].mean().to_csv(p)
    log.log(f"  Saved: {p}")

    log.log("  [2/5] Interferon (IFN) response signature...")
    ifn_genes = ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "OASL"]
    valid_ifn = [g for g in ifn_genes if g in adata.var_names]
    if valid_ifn:
        sc.tl.score_genes(adata, gene_list=valid_ifn, score_name="IFN_score")
    else:
        log.log("  Warning: IFN genes missing - IFN_score set to NaN")
        adata.obs["IFN_score"] = np.nan
    p = paths["results"] / "tables" / "ifn_scores_by_condition.csv"
    adata.obs.groupby("condition")["IFN_score"].mean().to_csv(p)
    log.log(f"  Saved: {p}")

    log.log("  [3/5] Antigen presentation score (HLA pathway)...")
    ag_genes = ["HLA-DRA", "HLA-DRB1", "CD74", "B2M", "TAP1", "TAP2"]
    valid_ag = [g for g in ag_genes if g in adata.var_names]
    if valid_ag:
        sc.tl.score_genes(adata, gene_list=valid_ag, score_name="antigen_presentation_score")
    else:
        log.log("  Warning: antigen presentation genes missing")
        adata.obs["antigen_presentation_score"] = np.nan
    p = paths["results"] / "tables" / "antigen_presentation_scores.csv"
    adata.obs.groupby(["condition", "cell_type_marker"], observed=False)[
        "antigen_presentation_score"
    ].mean().to_csv(p)
    log.log(f"  Saved: {p}")

    log.log("  [4/5] Pseudobulk counts (sum UMIs per condition x cell type)...")
    counts_layer = adata.layers["counts"]
    matrix = counts_layer.toarray() if hasattr(counts_layer, "toarray") else counts_layer
    pseudobulk = (
        pd.DataFrame(matrix, index=adata.obs_names, columns=adata.var_names)
        .assign(
            condition=adata.obs["condition"].astype(str).values,
            cell_type=adata.obs["cell_type_marker"].astype(str).values,
        )
        .groupby(["condition", "cell_type"], observed=False)
        .sum()
    )
    p = paths["results"] / "tables" / "pseudobulk_counts_condition_celltype.csv"
    pseudobulk.to_csv(p)
    log.log(f"  Saved: {p} ({pseudobulk.shape[0]} groups)")

    log.log("  [5/5] CoDi vs marker annotation agreement...")
    agreement = (
        adata.obs["cell_type_codi_norm"].astype(str) == adata.obs["cell_type_marker"].astype(str)
    ).mean()
    p = paths["results"] / "tables" / "annotation_agreement_metrics.csv"
    pd.DataFrame({"metric": ["codi_marker_agreement"], "value": [agreement]}).to_csv(p, index=False)
    log.log(f"  Agreement fraction: {agreement:.3f} - saved {p}")


def save_core_figures(
    adata: sc.AnnData, paths: Dict[str, Path], log: PipelineLogger
) -> List[str]:
    saved = []
    fig_dir = paths["results"] / "figures"

    sc.pl.umap(adata, color=["condition"], show=False)
    f1 = fig_dir / "umap_condition.png"
    plt.savefig(f1, dpi=300, bbox_inches="tight")
    plt.close()
    saved.append(f1.name)
    log.log(f"  Saved: {f1}")

    sc.pl.umap(adata, color=["cluster"], show=False)
    f2 = fig_dir / "umap_clusters.png"
    plt.savefig(f2, dpi=300, bbox_inches="tight")
    plt.close()
    saved.append(f2.name)
    log.log(f"  Saved: {f2}")

    sc.pl.umap(adata, color=["cell_type_marker"], legend_loc="on data", show=False)
    f3 = fig_dir / "umap_celltypes_marker.png"
    plt.savefig(f3, dpi=300, bbox_inches="tight")
    plt.close()
    saved.append(f3.name)
    log.log(f"  Saved: {f3}")

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
            "Enrichr hits (GO/KEGG/Reactome) for significant upregulated genes."
        ),
        "size_specific_effects_summary.csv": (
            "Counts of DE genes unique to 40 nm, 200 nm, shared, or mix-only effects."
        ),
        "cell_cycle_scores_by_condition.csv": (
            "Mean S and G2M phase scores - proxy for proliferation/stress."
        ),
        "ifn_scores_by_condition.csv": (
            "Mean interferon signature - innate immune activation."
        ),
        "antigen_presentation_scores.csv": (
            "Mean HLA-related score by condition and cell type."
        ),
        "pseudobulk_counts_condition_celltype.csv": (
            "Summed UMI counts per conditionxcell type for bulk-style follow-up."
        ),
        "annotation_agreement_metrics.csv": (
            "Fraction of cells where CoDi reference label matches marker annotation."
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
        adata = merge_and_integrate(adata, cfg, log)

        log.section("STEP 3 - CELL-TYPE ANNOTATION")
        marker_based_annotation(adata, cfg["markers"], log)
        load_codi_annotations(paths, adata, log)

        log.section("STEP 4 - CORE UMAP FIGURES")
        save_core_figures(adata, paths, log)

        log.section("STEP 5 - CELL COMPOSITION")
        log.log("  Comparing immune cell proportions across PSNP conditions vs control...")
        composition_analysis(adata, paths, log)

        log.section("STEP 6 - DIFFERENTIAL EXPRESSION")
        log.log(
            f"  Wilcoxon test per cell type; min {cfg['de']['min_cells_per_group']} cells per group..."
        )
        de_all = differential_expression_by_celltype(adata, cfg, paths, log)

        log.section("STEP 7 - PATHWAY ENRICHMENT")
        pathway_enrichment(de_all, cfg, paths, log)

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
        log.log("  Next optional step: Rscript scripts/azimuth_annotation.R for Azimuth labels.")

    except Exception as exc:
        log.section("PIPELINE FAILED")
        log.log(f"  Error: {exc}")
        raise
    finally:
        log.close()


if __name__ == "__main__":
    main()
