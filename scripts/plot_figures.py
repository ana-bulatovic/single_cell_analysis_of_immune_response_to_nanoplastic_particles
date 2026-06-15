"""
Regenerate UMAP figures and dotplot from a saved integrated AnnData object.

Use this after tweaking plot styling without re-running the full pipeline
(DE, enrichment, etc. can take 30–90 minutes).

Requires: data/processed/integrated_annotated.h5ad
"""

from pathlib import Path
import sys

import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_pipeline import PipelineLogger, load_config, save_core_figures, setup_paths


def main() -> None:
    sc.settings.verbosity = 1
    sc.set_figure_params(dpi=100, facecolor="white")

    cfg = load_config()
    paths = setup_paths(cfg)
    h5ad_path = paths["processed"] / "integrated_annotated.h5ad"
    if not h5ad_path.exists():
        raise FileNotFoundError(
            f"Missing {h5ad_path}. Run first: python scripts/run_pipeline.py"
        )

    log = PipelineLogger(paths["run_logs"])
    try:
        log.section("REGENERATE FIGURES ONLY")
        log.log(f"  Loading: {h5ad_path}")
        adata = sc.read_h5ad(h5ad_path)
        log.log(f"  Object: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
        saved = save_core_figures(adata, paths, log, cfg)
        log.log(f"  Wrote {len(saved)} figure(s) to {paths['results'] / 'figures'}")
    finally:
        log.close()


if __name__ == "__main__":
    main()
