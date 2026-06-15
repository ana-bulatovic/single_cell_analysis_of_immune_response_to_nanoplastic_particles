"""
Generate supplementary figures from saved pipeline CSV tables (no h5ad required).

Reads results/tables/*.csv and writes figures under results/figures/supplementary/.

Usage:
  python scripts/plot_supplementary_figures.py
  python scripts/plot_supplementary_figures.py --skip-volcano
  python scripts/plot_supplementary_figures.py --cell-types Monocyte_CD14 NK_cell
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_pipeline import (
    CONDITION_COLORS,
    CONDITION_LABELS,
    EXPOSURE_COMPARISONS,
    GENE_SET_PREFIX,
    PipelineLogger,
    SIZE_EFFECT_CLASSES,
    _short_pathway_label,
    load_config,
    setup_paths,
)

EFFECT_CLASS_LABELS = {
    "unique_40nm": "40 nm only",
    "unique_200nm": "200 nm only",
    "shared_40_200": "40+200 nm",
    "shared_all_three": "All exposures",
    "mix_only_emergent": "Mix only",
}

DEFAULT_HEATMAP_CELL_TYPES = (
    "Monocyte_CD14",
    "NK_cell",
    "B_cell",
    "DC",
    "CD4_T",
)

COMPARISON_LABELS = {comp: label for comp, label in EXPOSURE_COMPARISONS}

COMPARISON_TO_CONDITION = {
    "PSNP_40nm_vs_control": "PSNP_40nm",
    "PSNP_200nm_vs_control": "PSNP_200nm",
    "PSNP_mix_40_200_vs_control": "PSNP_mix_40_200",
}

CONDITION_ORDER = ["control", "PSNP_40nm", "PSNP_200nm", "PSNP_mix_40_200"]


def _slug(value: str) -> str:
    return re.sub(r"[^\w.-]+", "_", str(value)).strip("_")


def _sig_mask(df: pd.DataFrame, cfg: Dict) -> pd.Series:
    padj = cfg["de"]["pval_adj_threshold"]
    lfc = cfg["de"]["logfc_threshold"]
    return (df["pvals_adj"] < padj) & (df["logfoldchanges"].abs() > lfc)


def _require_table(path: Path, name: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing {name} ({path}). Run: python scripts/run_pipeline.py")
    return pd.read_csv(path)


def plot_de_summary_bar(de: pd.DataFrame, cfg: Dict, out_path: Path) -> None:
    sig = de.loc[_sig_mask(de, cfg)].copy()
    counts = (
        sig.groupby(["cell_type", "comparison"], observed=False)
        .size()
        .reset_index(name="n_sig")
    )
    if counts.empty:
        return

    counts["comparison_label"] = counts["comparison"].map(COMPARISON_LABELS).fillna(counts["comparison"])
    counts["cell_type_label"] = counts["cell_type"].str.replace("_", " ")

    order_ct = sorted(counts["cell_type"].unique())
    comp_order = [c for c, _ in EXPOSURE_COMPARISONS if c in set(counts["comparison"])]

    fig, ax = plt.subplots(figsize=(max(9, 0.55 * len(order_ct) + 4), 5.5))
    sns.barplot(
        data=counts,
        x="cell_type_label",
        y="n_sig",
        hue="comparison_label",
        order=[c.replace("_", " ") for c in order_ct],
        hue_order=[COMPARISON_LABELS[c] for c in comp_order],
        ax=ax,
        palette=[CONDITION_COLORS[COMPARISON_TO_CONDITION[c]] for c in comp_order],
    )
    ax.set_xlabel("Cell type")
    ax.set_ylabel("Significant DE genes")
    ax.set_title("Differential expression summary (padj & |logFC| thresholds)")
    ax.tick_params(axis="x", rotation=35, labelsize=9)
    ax.legend(title="Comparison", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_volcano(
    de_sub: pd.DataFrame,
    cfg: Dict,
    title: str,
    out_path: Path,
    label_top_n: int = 12,
) -> None:
    if de_sub.empty:
        return

    padj_thr = cfg["de"]["pval_adj_threshold"]
    lfc_thr = cfg["de"]["logfc_threshold"]
    df = de_sub.copy()
    df["neg_log_padj"] = -np.log10(df["pvals_adj"].clip(lower=1e-300))

    sig = _sig_mask(df, cfg)
    colors = np.where(
        df["logfoldchanges"] > lfc_thr,
        "#C44E52",
        np.where(df["logfoldchanges"] < -lfc_thr, "#4C72B0", "#BBBBBB"),
    )
    colors = np.where(sig, colors, "#DDDDDD")

    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.scatter(df["logfoldchanges"], df["neg_log_padj"], c=colors, s=8, alpha=0.75, linewidths=0)
    ax.axvline(lfc_thr, color="#666666", linestyle="--", linewidth=0.8)
    ax.axvline(-lfc_thr, color="#666666", linestyle="--", linewidth=0.8)
    ax.axhline(-np.log10(padj_thr), color="#666666", linestyle="--", linewidth=0.8)
    ax.set_xlabel("log2 fold change vs control")
    ax.set_ylabel("-log10(adjusted p-value)")
    ax.set_title(title)

    up = df[sig & (df["logfoldchanges"] > lfc_thr)].nlargest(label_top_n, "neg_log_padj")
    down = df[sig & (df["logfoldchanges"] < -lfc_thr)].nlargest(label_top_n, "neg_log_padj")
    for _, row in pd.concat([up, down]).iterrows():
        ax.annotate(
            row["names"],
            (row["logfoldchanges"], row["neg_log_padj"]),
            fontsize=7,
            alpha=0.9,
            xytext=(3, 3),
            textcoords="offset points",
        )

    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_top_de_heatmap(
    de: pd.DataFrame,
    cfg: Dict,
    cell_type: str,
    out_path: Path,
    top_n: int = 6,
) -> None:
    padj_thr = cfg["de"]["pval_adj_threshold"]
    lfc_thr = cfg["de"]["logfc_threshold"]
    sub = de[de["cell_type"] == cell_type].copy()
    if sub.empty:
        return

    selected_genes: List[str] = []
    comp_order = [c for c, _ in EXPOSURE_COMPARISONS]
    for comp in comp_order:
        comp_df = sub[sub["comparison"] == comp]
        sig_up = comp_df[
            (comp_df["pvals_adj"] < padj_thr) & (comp_df["logfoldchanges"] > lfc_thr)
        ].nlargest(top_n, "logfoldchanges")
        selected_genes.extend(sig_up["names"].astype(str).tolist())

    genes = list(dict.fromkeys(selected_genes))[: top_n * 3]
    if not genes:
        return

    pivot = (
        sub[sub["names"].isin(genes)]
        .pivot_table(index="names", columns="comparison", values="logfoldchanges", aggfunc="first")
        .reindex(index=genes, columns=comp_order)
    )
    pivot.columns = [COMPARISON_LABELS.get(c, c) for c in pivot.columns]

    fig_h = max(4.0, 0.28 * len(genes) + 1.5)
    fig, ax = plt.subplots(figsize=(5.5, fig_h))
    sns.heatmap(
        pivot,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "log2 FC vs control", "shrink": 0.7},
    )
    ax.set_xlabel("Exposure")
    ax.set_ylabel("Top upregulated genes")
    ax.set_title(f"{cell_type.replace('_', ' ')} — top DE genes")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_composition_delta_heatmap(comp: pd.DataFrame, out_path: Path) -> None:
    wide = comp.pivot_table(
        index="cell_type_marker", columns="condition", values="fraction", aggfunc="first"
    )
    wide = wide.reindex(columns=[c for c in CONDITION_ORDER if c in wide.columns])
    if "control" not in wide.columns:
        return

    delta = wide.drop(columns=["control"]).subtract(wide["control"], axis=0)
    delta.columns = [CONDITION_LABELS.get(c, c) for c in delta.columns]
    delta.index = [str(i).replace("_", " ") for i in delta.index]

    vmax = max(0.05, float(np.nanmax(np.abs(delta.values))))
    fig, ax = plt.subplots(figsize=(5.5, max(4.5, 0.35 * len(delta) + 1.5)))
    sns.heatmap(
        delta,
        ax=ax,
        cmap="PiYG",
        center=0,
        vmin=-vmax,
        vmax=vmax,
        annot=True,
        fmt=".2f",
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Δ fraction vs control", "shrink": 0.7},
    )
    ax.set_xlabel("Exposure")
    ax.set_ylabel("Cell type")
    ax.set_title("Cell composition change relative to control")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_size_specific_counts(summary: pd.DataFrame, out_path: Path) -> None:
    if summary.empty:
        return

    plot_df = summary.copy()
    plot_df["effect_label"] = plot_df["effect_class"].map(EFFECT_CLASS_LABELS).fillna(plot_df["effect_class"])
    plot_df["cell_type_label"] = plot_df["cell_type"].str.replace("_", " ")

    effect_order = [EFFECT_CLASS_LABELS[c] for c in SIZE_EFFECT_CLASSES if c in set(plot_df["effect_class"])]
    fig, ax = plt.subplots(figsize=(max(10, 0.5 * plot_df["cell_type"].nunique() + 5), 5.5))
    sns.barplot(
        data=plot_df,
        x="cell_type_label",
        y="n_genes",
        hue="effect_label",
        hue_order=effect_order,
        ax=ax,
    )
    ax.set_xlabel("Cell type")
    ax.set_ylabel("Number of significant DE genes")
    ax.set_title("Size-specific DE gene classes")
    ax.tick_params(axis="x", rotation=35, labelsize=9)
    ax.legend(title="Effect class", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_size_specific_heatmap(summary: pd.DataFrame, out_path: Path) -> None:
    if summary.empty:
        return

    pivot = summary.pivot_table(
        index="cell_type", columns="effect_class", values="n_genes", aggfunc="first"
    )
    pivot = pivot.reindex(columns=[c for c in SIZE_EFFECT_CLASSES if c in pivot.columns])
    pivot.index = [str(i).replace("_", " ") for i in pivot.index]
    pivot.columns = [EFFECT_CLASS_LABELS.get(c, c) for c in pivot.columns]

    fig, ax = plt.subplots(figsize=(7.5, max(4.5, 0.35 * len(pivot) + 1.5)))
    sns.heatmap(
        pivot,
        ax=ax,
        cmap="YlOrRd",
        annot=True,
        fmt="d",
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Gene count", "shrink": 0.7},
    )
    ax.set_xlabel("Size-effect class")
    ax.set_ylabel("Cell type")
    ax.set_title("Size-specific significant gene counts")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_module_scores(
    cc: pd.DataFrame,
    ifn: pd.DataFrame,
    out_path: Path,
) -> None:
    rows = []
    for _, r in cc.iterrows():
        cond = r["condition"]
        rows.append({"condition": cond, "score": r["S_score"], "module": "S phase"})
        rows.append({"condition": cond, "score": r["G2M_score"], "module": "G2/M phase"})
    for _, r in ifn.iterrows():
        rows.append({"condition": r["condition"], "score": r["IFN_score"], "module": "IFN signature"})

    df = pd.DataFrame(rows)
    if df.empty:
        return

    df["condition_label"] = df["condition"].map(CONDITION_LABELS).fillna(df["condition"])
    order = [CONDITION_LABELS[c] for c in CONDITION_ORDER if c in set(df["condition"])]

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    sns.barplot(
        data=df,
        x="module",
        y="score",
        hue="condition_label",
        hue_order=order,
        ax=ax,
        palette=[CONDITION_COLORS[c] for c in CONDITION_ORDER],
    )
    ax.axhline(0, color="#888888", linewidth=0.8)
    ax.set_xlabel("Gene module")
    ax.set_ylabel("Mean module score")
    ax.set_title("Cell-cycle and IFN module scores by condition")
    ax.legend(title="Condition", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_antigen_presentation_heatmap(ap: pd.DataFrame, out_path: Path) -> None:
    if ap.empty:
        return

    pivot = ap.pivot_table(
        index="cell_type_marker", columns="condition", values="antigen_presentation_score", aggfunc="first"
    )
    pivot = pivot.reindex(columns=[c for c in CONDITION_ORDER if c in pivot.columns])
    pivot.index = [str(i).replace("_", " ") for i in pivot.index]
    pivot.columns = [CONDITION_LABELS.get(c, c) for c in pivot.columns]

    fig, ax = plt.subplots(figsize=(6.5, max(4.5, 0.35 * len(pivot) + 1.5)))
    sns.heatmap(
        pivot,
        ax=ax,
        cmap="viridis",
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Antigen presentation score", "shrink": 0.7},
    )
    ax.set_xlabel("Condition")
    ax.set_ylabel("Cell type")
    ax.set_title("Antigen presentation module score")
    plt.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_size_specific_pathway_figures(
    enr: pd.DataFrame,
    cfg: Dict,
    out_dir: Path,
    cell_types: Optional[Sequence[str]] = None,
) -> List[Path]:
    if enr.empty:
        return []

    plot_cfg = cfg.get("pathway_plots", {})
    top_n = int(plot_cfg.get("top_terms", 15))
    fdr_thr = float(plot_cfg.get("fdr_threshold", 0.05))
    fdr_col = "Adjusted P-value"
    out_dir.mkdir(parents=True, exist_ok=True)
    saved: List[Path] = []

    groups = enr.groupby(["cell_type", "effect_class"], observed=False)
    for (ctype, effect_class), sub in groups:
        if cell_types and ctype not in cell_types:
            continue

        sub = sub[sub[fdr_col] < fdr_thr].copy()
        if sub.empty:
            continue

        sub["neg_log_fdr"] = -np.log10(sub[fdr_col].clip(lower=1e-300))
        best = (
            sub.sort_values(fdr_col)
            .groupby(["gene_set", "Term"], as_index=False)
            .first()
        )
        term_rank = (
            best.groupby("Term")["neg_log_fdr"]
            .max()
            .sort_values(ascending=False)
            .head(top_n)
        )
        top_terms = term_rank.index.tolist()
        if not top_terms:
            continue

        label_map = {
            row["Term"]: _short_pathway_label(row["Term"], row["gene_set"])
            for _, row in best[best["Term"].isin(top_terms)]
            .drop_duplicates("Term")
            .iterrows()
        }
        y_labels = [label_map.get(t, t) for t in top_terms]

        pivot = (
            best[best["Term"].isin(top_terms)]
            .pivot_table(index="Term", columns="gene_set", values="neg_log_fdr", aggfunc="max")
            .reindex(index=top_terms)
            .fillna(0.0)
        )
        pivot.index = y_labels
        pivot.columns = [GENE_SET_PREFIX.get(c, str(c)) for c in pivot.columns]

        fig_h = max(4.0, 0.32 * len(top_terms) + 1.8)
        fig, ax = plt.subplots(figsize=(6.5, fig_h))
        vmax = max(float(pivot.values.max()), 3.0)
        sns.heatmap(
            pivot,
            ax=ax,
            cmap="YlOrRd",
            vmin=0,
            vmax=vmax,
            linewidths=0.4,
            linecolor="white",
            cbar_kws={"label": "-log10(adjusted p-value)", "shrink": 0.6},
        )
        effect_label = EFFECT_CLASS_LABELS.get(effect_class, effect_class)
        ax.set_title(
            f"{str(ctype).replace('_', ' ')} — {effect_label} (UP genes)",
            fontsize=10,
            pad=8,
        )
        ax.set_xlabel("Database")
        ax.set_ylabel("Enriched pathway")
        ax.tick_params(axis="y", labelsize=8)
        plt.tight_layout()

        out_path = out_dir / f"pathways_{_slug(ctype)}_{effect_class}.png"
        fig.savefig(out_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        saved.append(out_path)

    return saved


def generate_supplementary_figures(
    paths: Dict[str, Path],
    cfg: Dict,
    log: PipelineLogger,
    *,
    skip_volcano: bool = False,
    cell_types: Optional[Sequence[str]] = None,
) -> List[Path]:
    tables = paths["results"] / "tables"
    fig_root = paths["results"] / "figures" / "supplementary"
    fig_root.mkdir(parents=True, exist_ok=True)
    de_dir = fig_root / "de"
    de_dir.mkdir(parents=True, exist_ok=True)
    pathway_dir = fig_root / "size_specific_pathways"
    saved: List[Path] = []

    de = _require_table(tables / "differential_expression_all.csv", "differential_expression_all.csv")
    comp = _require_table(tables / "cell_composition_by_condition.csv", "cell_composition_by_condition.csv")
    size_summary = _require_table(
        tables / "size_specific_effects_summary.csv", "size_specific_effects_summary.csv"
    )

    log.log("  DE summary bar chart...")
    p = fig_root / "de_summary_significant_genes.png"
    plot_de_summary_bar(de, cfg, p)
    saved.append(p)

    heatmap_types = list(cell_types) if cell_types else list(DEFAULT_HEATMAP_CELL_TYPES)
    for ctype in heatmap_types:
        if ctype not in set(de["cell_type"]):
            continue
        log.log(f"  Top DE heatmap: {ctype}...")
        p = de_dir / f"top_de_heatmap_{_slug(ctype)}.png"
        plot_top_de_heatmap(de, cfg, ctype, p)
        saved.append(p)

    if not skip_volcano:
        for (ctype, comp_name), sub in de.groupby(["cell_type", "comparison"], observed=False):
            if cell_types and ctype not in cell_types:
                continue
            comp_label = COMPARISON_LABELS.get(comp_name, comp_name)
            title = f"{ctype.replace('_', ' ')} — {comp_label} vs control"
            p = de_dir / f"volcano_{_slug(ctype)}_{_slug(comp_name)}.png"
            log.log(f"  Volcano: {ctype} / {comp_name}...")
            plot_volcano(sub, cfg, title, p)
            saved.append(p)

    log.log("  Composition delta heatmap...")
    p = fig_root / "composition_delta_heatmap.png"
    plot_composition_delta_heatmap(comp, p)
    saved.append(p)

    log.log("  Size-specific count plots...")
    p = fig_root / "size_specific_counts_barplot.png"
    plot_size_specific_counts(size_summary, p)
    saved.append(p)
    p = fig_root / "size_specific_counts_heatmap.png"
    plot_size_specific_heatmap(size_summary, p)
    saved.append(p)

    cc_path = tables / "cell_cycle_scores_by_condition.csv"
    ifn_path = tables / "ifn_scores_by_condition.csv"
    if cc_path.exists() and ifn_path.exists():
        log.log("  Module score bar chart...")
        p = fig_root / "module_scores_by_condition.png"
        plot_module_scores(pd.read_csv(cc_path), pd.read_csv(ifn_path), p)
        saved.append(p)

    ap_path = tables / "antigen_presentation_scores.csv"
    if ap_path.exists():
        log.log("  Antigen presentation heatmap...")
        p = fig_root / "antigen_presentation_heatmap.png"
        plot_antigen_presentation_heatmap(pd.read_csv(ap_path), p)
        saved.append(p)

    ss_path = tables / "size_specific_pathway_enrichment.csv"
    if ss_path.exists():
        log.log("  Size-specific pathway figures...")
        pathway_saved = plot_size_specific_pathway_figures(
            pd.read_csv(ss_path), cfg, pathway_dir, cell_types=cell_types
        )
        saved.extend(pathway_saved)
        log.log(f"    {len(pathway_saved)} size-specific pathway figure(s)")

    return saved


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot supplementary figures from pipeline CSV tables"
    )
    parser.add_argument(
        "--skip-volcano",
        action="store_true",
        help="Skip per-comparison volcano plots (faster)",
    )
    parser.add_argument(
        "--cell-types",
        nargs="*",
        default=None,
        help="Limit DE/pathway figures to these cell types (e.g. Monocyte_CD14 NK_cell)",
    )
    args = parser.parse_args()

    cfg = load_config()
    paths = setup_paths(cfg)
    log = PipelineLogger(paths["run_logs"])
    try:
        log.section("SUPPLEMENTARY FIGURES FROM CSV TABLES")
        saved = generate_supplementary_figures(
            paths,
            cfg,
            log,
            skip_volcano=args.skip_volcano,
            cell_types=args.cell_types,
        )
        log.log(f"  Done. {len(saved)} figure(s) written under results/figures/supplementary/")
    finally:
        log.close()


if __name__ == "__main__":
    main()
