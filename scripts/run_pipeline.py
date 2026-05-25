"""
Glavni scRNA-seq pipeline za analizu imunog odgovora na nanoplastične čestice.

Tok analize (redom):
  1. Učitavanje 4 uzorka (.h5ad) + QC filtriranje ćelija
  2. Spajanje uzoraka, normalizacija, HVG, integracija (Combat), UMAP, Leiden
  3. Anotacija tipova ćelija (marker geni) + učitavanje CoDi referenci
  4. Analiza kompozicije, diferencijalna ekspresija, pathway enrichment
  5. Size-specific efekti (40 nm vs 200 nm vs mix)
  6. Dodatne analize (cell cycle, IFN, antigen presentation, pseudobulk)
"""

from pathlib import Path
from typing import Dict, List

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import yaml


def load_config(path: str = "config/config.yaml") -> Dict:
    """
    Učitava YAML konfiguraciju (putanje, QC pragove, markere, DE parametre).

    Svi pragovi i liste gena su centralizovani u config/config.yaml
    da bi se analiza mogla ponoviti bez menjanja koda.
    """
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def setup_paths(cfg: Dict) -> Dict[str, Path]:
    """
    Kreira folder strukturu za raw/processed/results podatke.

    Vraća rečnik putanja:
      - raw: originalni .h5ad i CSV fajlovi sa Zenodo-a
      - processed: integrisani AnnData objekat
      - results: figure/ i tables/ podfolderi za izlaze
    """
    paths = {
        "raw": Path(cfg["data"]["raw_dir"]),
        "processed": Path(cfg["data"]["processed_dir"]),
        "results": Path(cfg["data"]["results_dir"]),
    }
    for p in paths.values():
        p.mkdir(parents=True, exist_ok=True)
    (paths["results"] / "figures").mkdir(parents=True, exist_ok=True)
    (paths["results"] / "tables").mkdir(parents=True, exist_ok=True)
    return paths


def read_and_qc_sample(file_path: Path, sample_id: str, condition: str, qc_cfg: Dict) -> sc.AnnData:
    """
    Učitava jedan uzorak i primenjuje kontrolu kvaliteta (QC).

    Koraci:
      1. Čita AnnData (.h5ad) fajl
      2. Dodaje metadata: sample_id i condition (40nm / 200nm / mix / control)
      3. Računa QC metrike: broj gena, ukupni readovi, % mitohondrijskih gena (MT-)
      4. Filtrira ćelije koje ne prolaze pragove iz config.yaml

    Zašto ovi filteri:
      - min_genes: uklanja prazne/low-quality ćelije
      - max_genes: uklanja potencijalne duplete
      - min_counts: uklanja ćelije sa premalo RNA signala
      - max_mt_percent: uklanja oštećene/stresirane ćelije (visok MT signal)
    """
    adata = sc.read_h5ad(file_path)
    adata.obs["sample_id"] = sample_id
    adata.obs["condition"] = condition

    # Mitohondrijski geni počinju sa "MT-" (ljudski genom)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    adata = adata[
        (adata.obs["n_genes_by_counts"] >= qc_cfg["min_genes"])
        & (adata.obs["n_genes_by_counts"] <= qc_cfg["max_genes"])
        & (adata.obs["total_counts"] >= qc_cfg["min_counts"])
        & (adata.obs["pct_counts_mt"] <= qc_cfg["max_mt_percent"]),
        :
    ].copy()
    return adata


def merge_and_integrate(adatas: List[sc.AnnData], cfg: Dict) -> sc.AnnData:
    """
    Spaja 4 uzorka u jedan objekat i radi integraciju + klasterizaciju.

    Koraci:
      1. concat — spaja sve uzorke u jedan AnnData (outer join po genima)
      2. Čuva sirove brojeve u layer["counts"] pre normalizacije
      3. normalize_total + log1p — standardna normalizacija za scRNA-seq
      4. highly_variable_genes — bira ~3000 najinformativnijih gena (HVG)
      5. scale — skalira ekspresiju (ograničeno na max_value=10)
      6. PCA — smanjuje dimenzionalnost
      7. Combat — uklanja batch efekat između 4 uzorka (sample_id)
         (Harmony je opcija, ali na Windows-u često ne radi)
      8. neighbors + UMAP — 2D vizualizacija ćelija
      9. Leiden — nelinearna klasterizacija ćelija u grupe
    """
    adata = sc.concat(
        adatas,
        join="outer",
        label="batch",
        keys=[a.obs["sample_id"].iloc[0] for a in adatas],
    )
    # Nakon spajanja, barcode-ovi moraju biti jedinstveni
    adata.obs_names_make_unique()
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=cfg["preprocessing"]["n_hvgs"],
        flavor="seurat",
        batch_key="sample_id",
    )
    # Dalje radimo samo sa HVG genima (manje šuma, brža analiza)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=cfg["preprocessing"]["n_pcs"], svd_solver="arpack")

    integration_method = cfg["preprocessing"].get("integration_method", "combat")
    if integration_method == "harmony":
        try:
            sc.external.pp.harmony_integrate(adata, key="sample_id")
            rep_key = "X_pca_harmony"
        except Exception:
            print("Harmony nije dostupan, prelazimo na Combat.")
            sc.pp.combat(adata, key="sample_id")
            rep_key = "X_pca"
    else:
        # Combat koriguje sistematske razlike između uzoraka u PCA prostoru
        sc.pp.combat(adata, key="sample_id")
        rep_key = "X_pca"

    sc.pp.neighbors(adata, use_rep=rep_key, n_neighbors=cfg["preprocessing"]["neighbors_k"])
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=cfg["preprocessing"]["leiden_resolution"], key_added="cluster")
    return adata


def marker_based_annotation(adata: sc.AnnData, marker_dict: Dict[str, List[str]]) -> None:
    """
    Dodeljuje tip ćelije na osnovu poznatih marker gena (PBMC).

    Za svaki tip (T, B, NK, monocyte, ...) računa prosečnu ekspresiju
    njegovih markera po ćeliji. Ćelija dobija tip sa najvišim skorom.

    Primer: CD3D/CD3E → T ćelija; MS4A1/CD79A → B ćelija; LYZ/S100A8 → monocyte.

    Rezultat se čuva u adata.obs["cell_type_marker"].
    """
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


def load_codi_annotations(paths: Dict[str, Path], adata: sc.AnnData) -> None:
    """
    Učitava CoDi anotacije iz referentnih CSV fajlova (Zenodo dataset).

    CoDi je spoljašnji model za predikciju tipa ćelije. Koristimo ga da:
      - uporedimo naše marker-based rezultate sa referencom
      - izračunamo slaganje (agreement) između dve metode

    Mapiranje: barcode ćelije → CoDi label u adata.obs["cell_type_codi"].
    """
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
        return

    codi_all = pd.concat(codi_frames, ignore_index=True).drop_duplicates(subset=["cell_id"])
    codi_all = codi_all.assign(cell_id_clean=codi_all["cell_id"].str.replace(r"-\d+$", "", regex=True))
    codi_all = codi_all.set_index("cell_id_clean")
    adata.obs["cell_id_clean"] = adata.obs_names.str.replace(r"-\d+$", "", regex=True)
    adata.obs["cell_type_codi"] = adata.obs["cell_id_clean"].map(codi_all["cell_type_codi"]).fillna("NA")

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


def composition_analysis(adata: sc.AnnData, paths: Dict[str, Path]) -> pd.DataFrame:
    """
    Poredi proporcije tipova ćelija između 4 uslova (40nm, 200nm, mix, control).

    Za svaki condition računa:
      - n_cells: broj ćelija po tipu
      - fraction: udeo tipa u odnosu na sve ćelije tog uslova

    Izlaz:
      - CSV tabela: results/tables/cell_composition_by_condition.csv
      - Barplot: results/figures/composition_barplot.png

    Biološki smisao: nanoplastici mogu menjati broj monocita, T ćelija, itd.
    """
    comp = (
        adata.obs.groupby(["condition", "cell_type_marker"])
        .size()
        .reset_index(name="n_cells")
        .sort_values(["condition", "n_cells"], ascending=[True, False])
    )
    comp["fraction"] = comp["n_cells"] / comp.groupby("condition")["n_cells"].transform("sum")
    comp.to_csv(paths["results"] / "tables" / "cell_composition_by_condition.csv", index=False)

    plt.figure(figsize=(10, 5))
    sns.barplot(data=comp, x="condition", y="fraction", hue="cell_type_marker")
    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig(paths["results"] / "figures" / "composition_barplot.png", dpi=300)
    plt.close()
    return comp


def differential_expression_by_celltype(adata: sc.AnnData, cfg: Dict, paths: Dict[str, Path]) -> pd.DataFrame:
    """
    Diferencijalna ekspresija (DE): exposure vs control, odvojeno po tipu ćelije.

    Za svaki major cell type i svaki uslov (40nm, 200nm, mix):
      - uzima samo ćelije tog tipa
      - poredi exposure grupu sa control grupom (Wilcoxon test)
      - vraća gene sa logFC, p-value, adjusted p-value

    Zašto po cell type: mešanje tipova bi dalo pogrešan signal
    (npr. promena u monocytima bi izgledala kao promena u T ćelijama).

    Izlaz: results/tables/differential_expression_all.csv
    """
    conditions = ["PSNP_40nm", "PSNP_200nm", "PSNP_mix_40_200"]
    all_de = []
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

            sc.tl.rank_genes_groups(ad_sub, groupby="de_group", groups=[cond], reference="control", method="wilcoxon")
            de_df = sc.get.rank_genes_groups_df(ad_sub, group=cond)
            de_df["cell_type"] = ctype
            de_df["comparison"] = f"{cond}_vs_control"
            all_de.append(de_df)

    if not all_de:
        return pd.DataFrame()
    de_all = pd.concat(all_de, ignore_index=True)
    de_all.to_csv(paths["results"] / "tables" / "differential_expression_all.csv", index=False)
    return de_all


def pathway_enrichment(de_all: pd.DataFrame, cfg: Dict, paths: Dict[str, Path]) -> pd.DataFrame:
    """
    Pathway enrichment na značajno promenjenim genima iz DE analize.

    Za svaku kombinaciju (cell_type, comparison):
      1. Uzima gene sa padj < 0.05 i logFC > 0.25
      2. Šalje listu gena na Enrichr (preko gseapy)
      3. Traži obogaćenje u GO, KEGG i Reactome bazama

    Biološki smisao: ne gledamo samo pojedinačne gene, već funkcionalne
    puteve (npr. inflamatorni odgovor, interferon signaling, fagocitoza).
    """
    if de_all.empty:
        return pd.DataFrame()

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
                enr = gp.enrichr(gene_list=list(sig), gene_sets=gset, organism="Human", outdir=None)
                if enr.results is None or enr.results.empty:
                    continue
                tmp = enr.results.copy()
                tmp["cell_type"] = ctype
                tmp["comparison"] = comp
                tmp["gene_set"] = gset
                enr_frames.append(tmp)
            except Exception:
                continue

    if not enr_frames:
        return pd.DataFrame()
    enr_all = pd.concat(enr_frames, ignore_index=True)
    enr_all.to_csv(paths["results"] / "tables" / "pathway_enrichment_all.csv", index=False)
    return enr_all


def size_specific_effects(de_all: pd.DataFrame, cfg: Dict, paths: Dict[str, Path]) -> pd.DataFrame:
    """
    Klasifikuje DE gene po veličini čestica (size-specific efekti).

    Za svaki cell type pravi skupove gena:
      - unique_40nm: promenjeni samo kod 40 nm (ne i 200 nm / mix)
      - unique_200nm: promenjeni samo kod 200 nm
      - shared_40_200: zajednički za 40 i 200 nm (ali ne mix-only)
      - shared_all_three: promenjeni u sva tri exposure uslova
      - mix_only_emergent: promenjeni samo u mešavini (emergentni efekat)

    Ovo direktno odgovara na zadatak o biološkoj interpretaciji veličine čestica.
    """
    if de_all.empty:
        return pd.DataFrame()

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

    out = pd.DataFrame(results)
    out.to_csv(paths["results"] / "tables" / "size_specific_effects_summary.csv", index=False)
    return out


def additional_insights(adata: sc.AnnData, paths: Dict[str, Path]) -> None:
    """
    5 dodatnih analiza (za dodatne bodove na projektu).

    1) Cell-cycle score — da li nanoplastici menjaju proliferaciju ćelija
    2) IFN (interferon) signature — indikator antiviralnog/imunog aktiviranja
    3) Antigen presentation score — sposobnost prezentacije antigena (HLA geni)
    4) Pseudobulk — sabira readove po (condition, cell_type) za bulk-like analizu
    5) CoDi vs marker agreement — koliko se naša anotacija slaže sa referencom
    """
    # 1) Cell-cycle scoring (S faza i G2/M faza)
    s_genes = ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4"]
    g2m_genes = ["HMGB2", "CDK1", "NUSAP1", "TOP2A", "MKI67", "BIRC5"]
    valid_s = [g for g in s_genes if g in adata.var_names]
    valid_g2m = [g for g in g2m_genes if g in adata.var_names]
    if valid_s and valid_g2m:
        sc.tl.score_genes_cell_cycle(adata, s_genes=valid_s, g2m_genes=valid_g2m)
    else:
        print(
            "Warning: skipped cell cycle scoring because no valid cell cycle genes were found in the dataset."
        )
        adata.obs["S_score"] = np.nan
        adata.obs["G2M_score"] = np.nan
    adata.obs.groupby("condition")[["S_score", "G2M_score"]].mean().to_csv(
        paths["results"] / "tables" / "cell_cycle_scores_by_condition.csv"
    )

    # 2) Interferon signature (imuni odgovor na stres/infekciju)
    ifn_genes = ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "OASL"]
    valid_ifn = [g for g in ifn_genes if g in adata.var_names]
    if valid_ifn:
        sc.tl.score_genes(adata, gene_list=valid_ifn, score_name="IFN_score")
    else:
        print("Warning: skipped IFN scoring because no valid IFN genes were found in the dataset.")
        adata.obs["IFN_score"] = np.nan
    adata.obs.groupby("condition")["IFN_score"].mean().to_csv(
        paths["results"] / "tables" / "ifn_scores_by_condition.csv"
    )

    # 3) Antigen presentation (važno za APC ćelije: monocytes, DC, B)
    ag_genes = ["HLA-DRA", "HLA-DRB1", "CD74", "B2M", "TAP1", "TAP2"]
    valid_ag = [g for g in ag_genes if g in adata.var_names]
    if valid_ag:
        sc.tl.score_genes(adata, gene_list=valid_ag, score_name="antigen_presentation_score")
    else:
        print(
            "Warning: skipped antigen presentation scoring because no valid antigen presentation genes were found."
        )
        adata.obs["antigen_presentation_score"] = np.nan
    adata.obs.groupby(["condition", "cell_type_marker"], observed=False)["antigen_presentation_score"].mean().to_csv(
        paths["results"] / "tables" / "antigen_presentation_scores.csv"
    )

    # 4) Pseudobulk: sabira sirove count-e po uslovu i tipu ćelije
    pseudobulk = (
        pd.DataFrame(
            adata.layers["counts"].toarray() if hasattr(adata.layers["counts"], "toarray") else adata.layers["counts"],
            index=adata.obs_names,
            columns=adata.var_names,
        )
        .assign(condition=adata.obs["condition"].astype(str).values, cell_type=adata.obs["cell_type_marker"].astype(str).values)
        .groupby(["condition", "cell_type"], observed=False)
        .sum()
    )
    pseudobulk.to_csv(paths["results"] / "tables" / "pseudobulk_counts_condition_celltype.csv")

    # 5) Slaganje CoDi i marker anotacije (0–1, više = bolje slaganje)
    agreement = (
        adata.obs["cell_type_codi_norm"].astype(str) == adata.obs["cell_type_marker"].astype(str)
    ).mean()
    pd.DataFrame({"metric": ["codi_marker_agreement"], "value": [agreement]}).to_csv(
        paths["results"] / "tables" / "annotation_agreement_metrics.csv", index=False
    )


def save_core_figures(adata: sc.AnnData, paths: Dict[str, Path]) -> None:
    """
    Čuva ključne UMAP figure za prezentaciju i izveštaj.

    Tri grafikona:
      - umap_condition.png: bojenje po uslovu (40nm / 200nm / mix / control)
      - umap_clusters.png: Leiden klasteri (nelabeled grupe ćelija)
      - umap_celltypes_marker.png: tipovi ćelija po marker anotaciji
    """
    sc.pl.umap(adata, color=["condition"], show=False)
    plt.savefig(paths["results"] / "figures" / "umap_condition.png", dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.umap(adata, color=["cluster"], show=False)
    plt.savefig(paths["results"] / "figures" / "umap_clusters.png", dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.umap(adata, color=["cell_type_marker"], legend_loc="on data", show=False)
    plt.savefig(paths["results"] / "figures" / "umap_celltypes_marker.png", dpi=300, bbox_inches="tight")
    plt.close()


def main():
    """
    Glavna funkcija — pokreće ceo pipeline od početka do kraja.

    Redosled izvršavanja:
      1. Učitaj config i napravi foldere
      2. QC za svaki od 4 uzorka
      3. Integracija + klasterizacija
      4. Anotacija + CoDi
      5. Figure (UMAP)
      6. Kompozicija, DE, pathway, size-specific
      7. Dodatne analize
      8. Sačuvaj finalni .h5ad objekat
    """
    sc.settings.verbosity = 2
    sc.set_figure_params(dpi=100, facecolor="white")

    cfg = load_config()
    paths = setup_paths(cfg)

    # Korak 1–2: učitavanje i QC po uzorku
    adatas = []
    for sample_id, meta in cfg["samples"].items():
        f = paths["raw"] / f"{sample_id}.h5ad"
        if not f.exists():
            raise FileNotFoundError(f"Missing file: {f}")
        adatas.append(read_and_qc_sample(f, sample_id, meta["condition"], cfg["qc"]))

    # Korak 3: spajanje, integracija, UMAP, klasteri
    adata = merge_and_integrate(adatas, cfg)
    marker_based_annotation(adata, cfg["markers"])
    load_codi_annotations(paths, adata)
    save_core_figures(adata, paths)

    # Korak 4–7: statističke analize i dodatni uvidi
    composition_analysis(adata, paths)
    de_all = differential_expression_by_celltype(adata, cfg, paths)
    pathway_enrichment(de_all, cfg, paths)
    size_specific_effects(de_all, cfg, paths)
    additional_insights(adata, paths)

    # Korak 8: finalni integrisani objekat za kasniju upotrebu (npr. Azimuth u R)
    adata.write(paths["processed"] / "integrated_annotated.h5ad")
    print("Pipeline finished successfully.")


if __name__ == "__main__":
    main()
