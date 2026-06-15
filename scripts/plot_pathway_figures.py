"""
Regenerate pathway enrichment figures from saved tables.

Optionally re-runs Enrichr for missing directions (e.g. DOWN) without the full pipeline.

Usage:
  python scripts/plot_pathway_figures.py
  python scripts/plot_pathway_figures.py --enrich-missing
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_pipeline import (
    PipelineLogger,
    load_config,
    pathway_enrichment,
    plot_pathway_enrichment_figures,
    setup_paths,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot pathway enrichment figures")
    parser.add_argument(
        "--enrich-missing",
        action="store_true",
        help="Re-run Enrichr if DE table exists but enrichment is incomplete",
    )
    args = parser.parse_args()

    cfg = load_config()
    paths = setup_paths(cfg)
    tables = paths["results"] / "tables"
    enr_path = tables / "pathway_enrichment_all.csv"
    de_path = tables / "differential_expression_all.csv"

    log = PipelineLogger(paths["run_logs"])
    try:
        log.section("PATHWAY ENRICHMENT FIGURES ONLY")

        enr_all = pd.DataFrame()
        if enr_path.exists():
            enr_all = pd.read_csv(enr_path)
            if "direction" not in enr_all.columns:
                enr_all["direction"] = "UP"
                log.log("  Added direction=UP to legacy enrichment table")

        missing_directions: tuple = ()
        if enr_all.empty:
            missing_directions = ("UP", "DOWN")
        elif args.enrich_missing:
            have = set(enr_all.get("direction", pd.Series(dtype=str)).unique())
            missing_directions = tuple(d for d in ("UP", "DOWN") if d not in have)

        if missing_directions:
            if not de_path.exists():
                raise FileNotFoundError(
                    f"Missing {de_path}. Run: python scripts/run_pipeline.py"
                )
            log.log(f"  Running enrichment for: {', '.join(missing_directions)}")
            de_all = pd.read_csv(de_path)
            enr_all = pathway_enrichment(
                de_all, cfg, paths, log, directions=missing_directions
            )

        saved = plot_pathway_enrichment_figures(enr_all, cfg, paths, log)
        log.log(f"  Done. {len(saved)} figure(s) written.")
    finally:
        log.close()


if __name__ == "__main__":
    main()
