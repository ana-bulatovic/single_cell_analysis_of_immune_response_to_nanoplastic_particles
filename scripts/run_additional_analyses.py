"""
Run additional analyses standalone (STEP 9) without re-running DE/enrichment.

Requires: data/processed/integrated_annotated.h5ad
Optional: results/tables/azimuth_annotations.csv (for Azimuth cross-validation)

Usage:
  python scripts/run_additional_analyses.py
  python scripts/run_additional_analyses.py --refresh-markers   # after config.yaml marker changes
  python scripts/run_additional_analyses.py --no-plots
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_pipeline import (
    PipelineLogger,
    additional_insights,
    load_config,
    load_integrated_adata,
    marker_based_annotation,
    resolve_marker_dict,
    setup_paths,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run module scores, pseudobulk, and annotation validation only"
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Skip figure generation",
    )
    parser.add_argument(
        "--no-interpretation",
        action="store_true",
        help="Skip Serbian interpretation markdown",
    )
    parser.add_argument(
        "--refresh-markers",
        action="store_true",
        help="Re-assign cell_type_marker from config.yaml and update integrated_annotated.h5ad",
    )
    args = parser.parse_args()

    cfg = load_config()
    paths = setup_paths(cfg)
    log = PipelineLogger(paths["run_logs"])

    try:
        log.section("ADDITIONAL ANALYSES (standalone)")
        log.log("  Analyses: module scores, pseudobulk, CoDi/Azimuth validation")
        log.log("  Input: data/processed/integrated_annotated.h5ad")
        log.log("")

        adata = load_integrated_adata(paths, log)

        if args.refresh_markers:
            log.log("  Re-assigning cell_type_marker from config.yaml ...")
            marker_dict = resolve_marker_dict(cfg, paths, log)
            marker_based_annotation(adata, marker_dict, log)
            h5ad_out = paths["processed"] / "integrated_annotated.h5ad"
            adata.write(h5ad_out)
            log.log(f"  Updated: {h5ad_out}")

        additional_insights(
            adata,
            paths,
            log,
            plot_figures=not args.no_plots,
            write_interpretation=not args.no_interpretation,
        )

        log.section("OUTPUTS")
        log.log("  Tables: results/tables/module_scores_*.csv, annotation_*.csv")
        log.log("  Figures: results/figures/additional_analyses/")
        log.log("    incl. CoDi vs marker contingency heatmaps")
        log.log("  Interpretation: results/tables/additional_analyses_interpretation_SR.md")
        log.log("")
        log.log("  Done.")
    finally:
        log.close()


if __name__ == "__main__":
    main()
