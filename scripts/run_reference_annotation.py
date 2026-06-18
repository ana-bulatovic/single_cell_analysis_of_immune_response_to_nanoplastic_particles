"""
Run reference-based annotation methods (ref.Rds atlas).

Workflow (professor requirement):
  1. Export ref.Rds -> ref_pbmc.h5ad
  2. Seurat/Azimuth (R) -> azimuth_annotations.csv
  3. Tangram + cell2location-style mapping (Python)
  4. CoDi labels already in data/raw/*_CoDi_KLD.csv

Usage:
  python scripts/run_reference_annotation.py
  python scripts/run_reference_annotation.py --skip-r
  python scripts/run_reference_annotation.py --only-python
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from annotation_methods import (
    assign_primary_cell_type,
    ensure_ref_h5ad,
    export_method_agreement,
    load_codi_annotations,
    load_query_for_mapping,
    load_seurat_annotations,
    run_cell2location_mapping,
    run_python_reference_methods,
    run_tangram_mapping,
    validate_with_literature_markers,
)
from run_pipeline import (
    PipelineLogger,
    load_config,
    marker_based_annotation,
    resolve_marker_dict,
    setup_paths,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Reference annotation with ref.Rds atlas")
    parser.add_argument("--skip-r", action="store_true", help="Skip Seurat/Azimuth R step")
    parser.add_argument("--only-python", action="store_true", help="Only Tangram + cell2location")
    parser.add_argument("--no-plots", action="store_true", help="Skip marker validation dotplots")
    args = parser.parse_args()

    cfg = load_config()
    paths = setup_paths(cfg)
    log = PipelineLogger(paths["run_logs"])

    try:
        log.section("REFERENCE ANNOTATION (ref.Rds)")
        log.log("  Methods: Seurat, CoDi, Tangram, cell2location")
        log.log("  Validation: literature markers from config.yaml")
        log.log("")

        ref_path = ensure_ref_h5ad(paths, log)
        if ref_path is None:
            log.log("  ERROR: Could not build reference from ref.Rds")
            return

        if not args.only_python:
            pre_hvg = paths["processed"] / "integrated_pre_hvg.h5ad"
            if not pre_hvg.exists():
                log.log(f"  Missing {pre_hvg} — run: python scripts/run_pipeline.py first")
                return

            if not args.skip_r:
                log.log("  Preparing Azimuth inputs...")
                subprocess.run([sys.executable, "scripts/prepare_azimuth_h5ad.py"], check=False)
                log.log("  Running Seurat/Azimuth (ref.Rds)...")
                r = subprocess.run(["Rscript", "scripts/azimuth_annotation.R"], cwd=paths["project_root"])
                if r.returncode != 0:
                    log.log("  Warning: Azimuth R step failed — continue with available methods")

        query = load_query_for_mapping(paths, log)
        if query is None:
            return

        import scanpy as sc

        adata = sc.read_h5ad(paths["processed"] / "integrated_annotated.h5ad")
        marker_dict = resolve_marker_dict(cfg, paths, log)
        marker_based_annotation(adata, marker_dict, log)
        load_codi_annotations(adata, paths, log)
        load_seurat_annotations(adata, paths, cfg, log)

        ref = sc.read_h5ad(ref_path)
        if cfg.get("annotation", {}).get("run_tangram", True):
            tangram = run_tangram_mapping(query, ref, paths, log)
            adata.obs["cell_type_tangram"] = adata.obs_names.map(tangram.to_dict()).fillna("unmapped")
        if cfg.get("annotation", {}).get("run_cell2location", True):
            c2l = run_cell2location_mapping(query, ref, paths, log)
            adata.obs["cell_type_cell2location"] = adata.obs_names.map(c2l.to_dict()).fillna("unmapped")

        assign_primary_cell_type(adata, cfg, log)
        if not args.no_plots:
            validate_with_literature_markers(adata, marker_dict, paths, log)
        export_method_agreement(adata, paths, log)

        log.section("DONE")
        log.log("  Tables: results/tables/*_annotations.csv, annotation_method_agreement.csv")
        log.log("  Figures: results/figures/annotation_validation/")
        log.log("  Re-run full pipeline: python scripts/run_pipeline.py")
    finally:
        log.close()


if __name__ == "__main__":
    main()
