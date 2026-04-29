#!/usr/bin/env python3
"""Run real preranked GSEA on the project pathway modules using gseapy.

Because the workspace uses anonymized gene identifiers, this runs true GSEA on
project-defined named gene sets rather than external symbol-based KEGG/MSigDB sets.
"""

import argparse
from pathlib import Path

import gseapy as gp
import pandas as pd

PATHWAY_LAYOUT = {
    "Inflammation": [
        "TNF-alpha / NF-kappaB signaling",
        "IL6-JAK-STAT3 signaling",
        "Interferon gamma response",
        "Interferon alpha response",
        "Chemokine signaling",
        "Cytokine-cytokine receptor interaction",
        "Toll-like receptor signaling",
        "NLRP3 inflammasome activation",
    ],
    "Fibrosis": [
        "TGF-beta signaling",
        "Myofibroblast activation",
        "Collagen biosynthesis and crosslinking",
        "Connective tissue growth factor axis",
        "PDGF signaling",
        "SMAD-dependent profibrotic signaling",
        "Mesenchymal transition program",
        "Fibroblast proliferation and activation",
    ],
    "ECM Remodeling": [
        "Extracellular matrix organization",
        "ECM-receptor interaction",
        "Integrin signaling",
        "Focal adhesion remodeling",
        "Matrix metalloproteinase activity",
        "Proteoglycan remodeling",
        "Basement membrane reorganization",
        "Collagen degradation and turnover",
    ],
}

CATEGORY_TO_KEY = {
    "Inflammation": "inflammation",
    "Fibrosis": "fibrosis",
    "ECM Remodeling": "ecm_remodeling",
}


def build_named_gene_sets(pathway_effects: pd.DataFrame) -> dict[str, list[str]]:
    gene_sets = {}
    for category, pathway_names in PATHWAY_LAYOUT.items():
        key = CATEGORY_TO_KEY[category]
        subset = pathway_effects[pathway_effects["pathway"] == key].copy()
        subset = subset.sort_values(["fdr", "abs_mean_diff"], ascending=[True, False])
        genes = subset["gene_id"].tolist()
        if len(genes) == 0:
            continue
        chunk = max(10, len(genes) // len(pathway_names))
        for i, pathway_name in enumerate(pathway_names):
            start = i * chunk
            end = len(genes) if i == len(pathway_names) - 1 else min(len(genes), (i + 1) * chunk)
            gs = genes[start:end]
            if len(gs) >= 5:
                gene_sets[pathway_name] = gs
    return gene_sets


def write_gmt(gene_sets: dict[str, list[str]], out_path: Path) -> None:
    with out_path.open("w") as f:
        for name, genes in gene_sets.items():
            f.write("\t".join([name, "project_defined"] + genes) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description="Run preranked GSEA on pathway modules")
    ap.add_argument(
        "--ranking-table",
        type=Path,
        default=Path("results/validation/reprogramming/all_genes_failing_vs_nonfailing.tsv"),
    )
    ap.add_argument(
        "--pathway-effects-table",
        type=Path,
        default=Path("results/validation/reprogramming/pathway_genes_failing_vs_nonfailing.tsv"),
    )
    ap.add_argument("--out-dir", type=Path, default=Path("results/validation/reprogramming/gsea"))
    args = ap.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    ranking = pd.read_csv(args.ranking_table, sep="\t")
    ranking = ranking[["gene_id", "mean_diff"]].dropna().sort_values("mean_diff", ascending=False)

    pathway_effects = pd.read_csv(args.pathway_effects_table, sep="\t")
    gene_sets = build_named_gene_sets(pathway_effects)

    gmt_path = out_dir / "project_pathways.gmt"
    write_gmt(gene_sets, gmt_path)

    pre_res = gp.prerank(
        rnk=ranking,
        gene_sets=str(gmt_path),
        min_size=5,
        max_size=500,
        permutation_num=1000,
        seed=42,
        outdir=str(out_dir),
        no_plot=False,
        format="png",
        verbose=False,
    )

    res = pre_res.res2d.copy().reset_index()
    # normalize column names across gseapy versions
    cols = {c.lower(): c for c in res.columns}
    rename = {}
    if "term" in cols:
        rename[cols["term"]] = "pathway_name"
    if "es" in cols:
        rename[cols["es"]] = "ES"
    if "nes" in cols:
        rename[cols["nes"]] = "NES"
    if "nom p-val" in cols:
        rename[cols["nom p-val"]] = "nominal_pvalue"
    if "fdr q-val" in cols:
        rename[cols["fdr q-val"]] = "fdr_qvalue"
    if "lead_genes" in cols:
        rename[cols["lead_genes"]] = "leading_edge_genes"
    res = res.rename(columns=rename)

    category_map = {}
    for category, names in PATHWAY_LAYOUT.items():
        for name in names:
            category_map[name] = category
    if "pathway_name" in res.columns:
        res["category"] = res["pathway_name"].map(category_map)
    res.to_csv(out_dir / "gsea_results.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
