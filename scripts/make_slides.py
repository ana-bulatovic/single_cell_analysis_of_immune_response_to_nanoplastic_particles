"""
Automatsko generisanje PowerPoint prezentacije iz rezultata analize.

Čita figure iz results/figures/ i tabele iz results/tables/,
pa pravi deliverables/nanoplastic_scRNA_results.pptx.
"""

from pathlib import Path

import pandas as pd
from pptx import Presentation
from pptx.util import Inches, Pt


def add_title_slide(prs, title, subtitle):
    """Dodaje naslovni slajd (naziv projekta + kratak opis)."""
    slide = prs.slides.add_slide(prs.slide_layouts[0])
    slide.shapes.title.text = title
    slide.placeholders[1].text = subtitle


def add_bullet_slide(prs, title, bullets):
    """Dodaje slajd sa naslovom i listom bullet tačaka."""
    slide = prs.slides.add_slide(prs.slide_layouts[1])
    slide.shapes.title.text = title
    tf = slide.shapes.placeholders[1].text_frame
    tf.clear()
    for i, b in enumerate(bullets):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = b
        p.level = 0
        p.font.size = Pt(20)


def add_image_slide(prs, title, image_path):
    """
    Dodaje slajd sa slikom (UMAP, barplot, itd.).

    Ako slika ne postoji (pipeline još nije pokrenut), prikazuje poruku umesto slike.
    """
    slide = prs.slides.add_slide(prs.slide_layouts[5])
    slide.shapes.title.text = title
    if Path(image_path).exists():
        slide.shapes.add_picture(str(image_path), Inches(0.8), Inches(1.4), width=Inches(11.5))
    else:
        tx = slide.shapes.add_textbox(Inches(1), Inches(2), Inches(10), Inches(1))
        tx.text_frame.text = f"Image not found: {image_path}"


def main():
    """
    Sastavlja kompletnu prezentaciju slajd po slajd.

    Struktura:
      1. Naslov
      2. Dizajn studije
      3–5. UMAP i kompozicija (slike iz pipeline-a)
      6. Size-specific nalazi (auto iz CSV tabele)
      7. Dodatne analize
      8. Zaključak
    """
    results = Path("results")
    figs = results / "figures"
    tables = results / "tables"
    out = Path("deliverables")
    out.mkdir(parents=True, exist_ok=True)

    prs = Presentation()

    add_title_slide(
        prs,
        "Single-Cell Analysis of Immune Response to Nanoplastic Particles",
        "Donor PBMC exposed to 40 nm, 200 nm, mixture, and control",
    )

    add_bullet_slide(
        prs,
        "Study Design",
        [
            "4 conditions: 40 nm, 200 nm, 40+200 nm mix, and unexposed control",
            "Data modality: scRNA-seq (AnnData .h5ad)",
            "Goal: identify size-dependent and shared immune responses",
        ],
    )

    add_image_slide(prs, "UMAP by Condition", figs / "umap_condition.png")
    add_image_slide(prs, "UMAP by Clusters", figs / "umap_clusters.png")
    add_image_slide(prs, "Cell Type Composition", figs / "composition_barplot.png")

    # Automatski popunjava slajd sa top size-specific nalazima iz CSV-a
    summary_bullets = []
    size_file = tables / "size_specific_effects_summary.csv"
    if size_file.exists():
        df = pd.read_csv(size_file)
        top = df.sort_values("n_genes", ascending=False).head(5)
        summary_bullets = [f"{r.cell_type}: {r.effect_class} = {int(r.n_genes)} genes" for _, r in top.iterrows()]
    if not summary_bullets:
        summary_bullets = [
            "Run scripts/run_pipeline.py to generate size-specific summary table.",
            "Then rerun this script to auto-populate key findings.",
        ]
    add_bullet_slide(prs, "Size-Specific Findings (Auto-summary)", summary_bullets)

    add_bullet_slide(
        prs,
        "Additional Insights (3–5 Analyses)",
        [
            "Cell-cycle scoring across conditions",
            "Interferon response signature score by condition",
            "Antigen presentation signature by cell type and condition",
            "Pseudobulk expression matrix for downstream modeling",
            "CoDi label agreement with marker-based annotation",
        ],
    )

    add_bullet_slide(
        prs,
        "Conclusions",
        [
            "Nanoplastic size changes both cell composition and gene programs.",
            "Mixture can induce emergent transcriptional responses not seen in single-size exposure.",
            "Pathway-level results support immune activation and stress-related biology.",
        ],
    )

    output_file = out / "nanoplastic_scRNA_results.pptx"
    prs.save(output_file)
    print(f"Saved slides to: {output_file}")


if __name__ == "__main__":
    main()
