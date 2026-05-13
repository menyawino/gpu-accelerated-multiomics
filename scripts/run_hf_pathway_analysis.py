#!/usr/bin/env python3
"""Run HF vs non-HF pathway analyses across ORA, FCS, and topology views.

This script uses the real processed GSE116250 expression matrix and the existing
HF vs non-HF differential table to produce:
- ORA on HF-up and HF-down genes using Enrichr libraries
- GSEA-style functional class scoring on the ranked full gene list
- Topology-aware KEGG overlays via pathview for enriched KEGG pathways
"""

import argparse
import re
import subprocess
from io import StringIO
from pathlib import Path

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"
ENRICHR_LIBRARIES = [
    "KEGG_2021_Human",
    "Reactome_2022",
    "MSigDB_Hallmark_2020",
    "GO_Biological_Process_2023",
]
FCS_LIBRARIES = [
    "KEGG_2021_Human",
    "Reactome_2022",
    "MSigDB_Hallmark_2020",
]


def load_inputs(de_path: Path, expression_path: Path, metadata_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    de = pd.read_csv(de_path, sep="\t")
    expr = pd.read_csv(expression_path, sep="\t", index_col=0)
    meta = pd.read_csv(metadata_path, sep="\t")
    return de, expr, meta


def split_hf_gene_lists(de: pd.DataFrame, fdr_cutoff: float) -> tuple[list[str], list[str], list[str]]:
    de = de.dropna(subset=["gene_id", "mean_diff", "fdr"]).copy()
    de = de.sort_values(["fdr", "abs_mean_diff"], ascending=[True, False])
    universe = de["gene_id"].astype(str).drop_duplicates().tolist()
    up = de.loc[(de["fdr"] <= fdr_cutoff) & (de["mean_diff"] > 0), "gene_id"].astype(str).drop_duplicates().tolist()
    down = de.loc[(de["fdr"] <= fdr_cutoff) & (de["mean_diff"] < 0), "gene_id"].astype(str).drop_duplicates().tolist()
    return up, down, universe


def download_gmt(library_name: str, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{library_name}.gmt"
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path

    response = requests.get(
        f"{ENRICHR_BASE}/geneSetLibrary",
        params={"mode": "text", "libraryName": library_name},
        timeout=120,
    )
    response.raise_for_status()
    out_path.write_text(response.text)
    return out_path


def parse_gmt(path: Path) -> dict[str, list[str]]:
    gene_sets: dict[str, list[str]] = {}
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = re.split(r"\t+", line)
            if len(parts) < 2:
                parts = re.split(r"\s{2,}", line)
            term = parts[0].strip()
            genes = [part.strip() for part in parts[1:] if part.strip()]
            if len(genes) >= 5:
                gene_sets[term] = genes
    return gene_sets


def run_ora(
    gene_list: list[str],
    background: list[str],
    libraries: list[str],
    direction: str,
    gmt_dir: Path,
) -> pd.DataFrame:
    gene_set = set(gene_list)
    background_set = set(background)
    rows = []

    for library in libraries:
        gmt_path = download_gmt(library, gmt_dir)
        term_map = parse_gmt(gmt_path)

        for term, members in term_map.items():
            term_genes = background_set.intersection(members)
            if len(term_genes) < 5:
                continue

            overlap_genes = sorted(gene_set.intersection(term_genes))
            overlap = len(overlap_genes)
            if overlap == 0:
                continue

            universe_size = len(background_set)
            term_size = len(term_genes)
            input_size = len(gene_set)

            pvalue = hypergeom.sf(overlap - 1, universe_size, term_size, input_size)
            contingency_a = overlap
            contingency_b = max(input_size - overlap, 0)
            contingency_c = max(term_size - overlap, 0)
            contingency_d = max(universe_size - term_size - contingency_b, 0)
            odds_ratio = ((contingency_a + 0.5) * (contingency_d + 0.5)) / ((contingency_b + 0.5) * (contingency_c + 0.5))

            rows.append(
                {
                    "Gene_set": library,
                    "Term": term,
                    "Overlap": f"{overlap}/{term_size}",
                    "P-value": float(pvalue),
                    "Odds Ratio": float(odds_ratio),
                    "Combined Score": float(-np.log10(max(pvalue, 1e-300)) * np.log2(odds_ratio + 1.0)),
                    "Genes": ";".join(overlap_genes),
                    "direction": direction,
                    "input_gene_count": input_size,
                    "background_gene_count": universe_size,
                }
            )

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["Adjusted P-value"] = multipletests(out["P-value"].values, method="fdr_bh")[1]
    return out.sort_values(["Gene_set", "Adjusted P-value", "P-value", "Odds Ratio"], ascending=[True, True, True, False])


def prepare_ranking(de: pd.DataFrame) -> pd.DataFrame:
    rank = de[["gene_id", "mean_diff"]].dropna().copy()
    rank["gene_id"] = rank["gene_id"].astype(str)
    rank = rank.sort_values("mean_diff", ascending=False).drop_duplicates("gene_id")
    rank["mean_diff"] = rank["mean_diff"] + rank.groupby("mean_diff").cumcount() * 1e-12
    return rank


def run_fcs(rank: pd.DataFrame, gmt_dir: Path, libraries: list[str], permutations: int) -> pd.DataFrame:
    results = []
    for library in libraries:
        gmt_path = download_gmt(library, gmt_dir)
        gene_sets = parse_gmt(gmt_path)
        prerank = gp.prerank(
            rnk=rank,
            gene_sets=gene_sets,
            min_size=5,
            max_size=500,
            permutation_num=permutations,
            seed=42,
            outdir=None,
            no_plot=True,
            verbose=False,
        )
        result = prerank.res2d.copy().reset_index()
        if "Term" in result.columns and "index" in result.columns:
            result = result.rename(columns={"index": "rank"})
        elif "index" in result.columns:
            result = result.rename(columns={"index": "Term"})
        result["Gene_set"] = library
        results.append(result)

    out = pd.concat(results, ignore_index=True)
    rename_map = {
        "Name": "comparison",
        "Term": "Term",
        "ES": "ES",
        "NES": "NES",
        "NOM p-val": "NOM p-val",
        "FDR q-val": "FDR q-val",
        "FWER p-val": "FWER p-val",
        "Lead_genes": "Lead_genes",
        "Tag %": "Tag %",
        "Gene %": "Gene %",
    }
    out = out.rename(columns={k: v for k, v in rename_map.items() if k in out.columns})
    return out.sort_values(["Gene_set", "FDR q-val", "NOM p-val", "NES"], ascending=[True, True, True, False])


def phenotype_groups(meta: pd.DataFrame, expr: pd.DataFrame) -> tuple[list[str], list[str]]:
    sample_col = meta.columns[0]
    pheno_col = "phenotype" if "phenotype" in meta.columns else meta.columns[1]
    labels = dict(zip(meta[sample_col].astype(str), meta[pheno_col].astype(str)))
    non_hf = [sample for sample, label in labels.items() if label == "NF" and sample in expr.columns]
    hf = [sample for sample, label in labels.items() if label in {"DCM", "ICM"} and sample in expr.columns]
    return hf, non_hf


def score_gene_set(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    valid = [gene for gene in genes if gene in expr.index]
    if not valid:
        return pd.Series(dtype=float)
    subset = expr.loc[valid]
    z = subset.sub(subset.mean(axis=1), axis=0).div(subset.std(axis=1).replace(0, np.nan), axis=0)
    return z.fillna(0.0).mean(axis=0)


def summarize_top_fcs_pathways(fcs: pd.DataFrame, expr: pd.DataFrame, meta: pd.DataFrame, top_n: int) -> pd.DataFrame:
    hf_samples, non_hf_samples = phenotype_groups(meta, expr)
    top = fcs.loc[fcs["FDR q-val"] <= 0.25].copy()
    if top.empty:
        top = fcs.copy()
    top = top.sort_values(["NES", "FDR q-val"], ascending=[False, True]).head(top_n)

    rows = []
    for _, row in top.iterrows():
        genes = [gene.strip() for gene in str(row.get("Lead_genes", "")).split(";") if gene.strip()]
        scores = score_gene_set(expr, genes)
        if scores.empty:
            continue
        rows.append(
            {
                "Gene_set": row["Gene_set"],
                "Term": row["Term"],
                "NES": row["NES"],
                "FDR q-val": row["FDR q-val"],
                "hf_mean_score": float(scores[hf_samples].mean()),
                "non_hf_mean_score": float(scores[non_hf_samples].mean()),
                "score_diff": float(scores[hf_samples].mean() - scores[non_hf_samples].mean()),
                "n_leading_edge_genes": len(genes),
            }
        )
    return pd.DataFrame(rows)


def plot_top_terms(ora: pd.DataFrame, fcs: pd.DataFrame, out_dir: Path, top_n: int) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    ora_plot = ora.sort_values("Adjusted P-value").head(top_n).copy()
    ora_plot = ora_plot.iloc[::-1]
    ora_plot["Adjusted P-value"] = pd.to_numeric(ora_plot["Adjusted P-value"], errors="coerce")
    axes[0].barh(ora_plot["Term"], -np.log10(ora_plot["Adjusted P-value"].clip(lower=1e-300)), color="#b2182b")
    axes[0].set_xlabel(r"$-\log_{10}$(adjusted p)")
    axes[0].set_title("ORA top pathways")

    fcs_plot = fcs.sort_values("FDR q-val").head(top_n).copy()
    fcs_plot = fcs_plot.iloc[::-1]
    fcs_plot["NES"] = pd.to_numeric(fcs_plot["NES"], errors="coerce")
    colors = ["#2166ac" if nes < 0 else "#b2182b" for nes in fcs_plot["NES"]]
    axes[1].barh(fcs_plot["Term"], fcs_plot["NES"], color=colors)
    axes[1].set_xlabel("Normalized enrichment score")
    axes[1].set_title("FCS / preranked GSEA top pathways")

    fig.tight_layout()
    fig.savefig(str(out_dir / "hf_vs_nonhf_pathway_summary.png"), dpi=400, bbox_inches="tight")
    fig.savefig(str(out_dir / "hf_vs_nonhf_pathway_summary.pdf"), bbox_inches="tight")
    plt.close(fig)


def fetch_kegg_name_map() -> dict[str, str]:
    response = requests.get("https://rest.kegg.jp/list/pathway/hsa", timeout=120)
    response.raise_for_status()
    mapping = {}
    for line in response.text.splitlines():
        pathway_id, name = line.split("\t", 1)
        clean = re.sub(r"\s+- Homo sapiens.*$", "", name).strip().lower()
        mapping[clean] = pathway_id
    return mapping


def pick_top_kegg_pathways(ora: pd.DataFrame, fcs: pd.DataFrame, top_n: int) -> pd.DataFrame:
    kegg_map = fetch_kegg_name_map()

    ora_kegg = ora.loc[ora["Gene_set"] == "KEGG_2021_Human", ["direction", "Term", "Adjusted P-value"]].copy()
    ora_kegg["source"] = "ORA"
    ora_kegg["score"] = -np.log10(ora_kegg["Adjusted P-value"].clip(lower=1e-300))

    fcs_kegg = fcs.loc[fcs["Gene_set"] == "KEGG_2021_Human", ["Term", "NES", "FDR q-val"]].copy()
    fcs_kegg["direction"] = np.where(fcs_kegg["NES"] >= 0, "hf_up", "hf_down")
    fcs_kegg["source"] = "FCS"
    fcs_kegg["score"] = fcs_kegg["NES"].abs()

    combined = pd.concat([
        ora_kegg[["source", "direction", "Term", "score"]],
        fcs_kegg[["source", "direction", "Term", "score"]],
    ], ignore_index=True)
    combined["pathway_id"] = combined["Term"].str.lower().map(kegg_map)
    combined = combined.dropna(subset=["pathway_id"])
    combined = combined.sort_values("score", ascending=False).drop_duplicates("pathway_id")
    return combined.head(top_n)


def run_pathview_subset(de_path: Path, topology_df: pd.DataFrame, out_dir: Path) -> None:
    if topology_df.empty:
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = topology_df.copy()
    manifest["category"] = manifest["direction"] + "_" + manifest["source"].str.lower()
    manifest["pathway_name"] = manifest["Term"]
    manifest[["category", "pathway_id", "pathway_name"]].to_csv(
        out_dir / "selected_pathways.tsv", sep="\t", index=False
    )

    cmd = [
        "Rscript",
        "scripts/run_pathview_kegg_selected.R",
        "--de-table",
        str(de_path),
        "--pathway-table",
        str(out_dir / "selected_pathways.tsv"),
        "--out-dir",
        str(out_dir),
    ]
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run HF vs non-HF pathway analyses")
    parser.add_argument(
        "--de-table",
        type=Path,
        default=Path("results/real_processed/reprogramming/all_genes_failing_vs_nonfailing.tsv"),
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=Path("results/real_processed/inputs/gse116250_expression.tsv"),
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=Path("results/real_processed/inputs/gse116250_metadata.tsv"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results/real_processed/hf_vs_nonhf_pathways"),
    )
    parser.add_argument("--fdr-cutoff", type=float, default=0.05)
    parser.add_argument("--gsea-permutations", type=int, default=1000)
    parser.add_argument("--top-n", type=int, default=12)
    parser.add_argument("--topology-top-n", type=int, default=6)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    gmt_dir = args.out_dir / "gene_sets"
    topology_dir = args.out_dir / "topology"

    de, expr, meta = load_inputs(args.de_table, args.expression, args.metadata)
    hf_up, hf_down, universe = split_hf_gene_lists(de, args.fdr_cutoff)

    pd.DataFrame(
        {
            "direction": ["hf_up", "hf_down", "background"],
            "gene_count": [len(hf_up), len(hf_down), len(universe)],
        }
    ).to_csv(args.out_dir / "gene_list_sizes.tsv", sep="\t", index=False)

    ora = pd.concat(
        [
            run_ora(hf_up, universe, ENRICHR_LIBRARIES, "hf_up", gmt_dir),
            run_ora(hf_down, universe, ENRICHR_LIBRARIES, "hf_down", gmt_dir),
        ],
        ignore_index=True,
    )
    ora.to_csv(args.out_dir / "ora_results.tsv", sep="\t", index=False)

    rank = prepare_ranking(de)
    fcs = run_fcs(rank, gmt_dir, FCS_LIBRARIES, args.gsea_permutations)
    fcs.to_csv(args.out_dir / "fcs_prerank_results.tsv", sep="\t", index=False)

    fcs_scores = summarize_top_fcs_pathways(fcs, expr, meta, args.top_n)
    fcs_scores.to_csv(args.out_dir / "fcs_top_pathway_scores.tsv", sep="\t", index=False)

    top_ora = ora.sort_values("Adjusted P-value").groupby(["direction", "Gene_set"], as_index=False).head(max(1, args.top_n // 2))
    top_fcs = fcs.sort_values("FDR q-val").head(args.top_n)
    plot_top_terms(top_ora, top_fcs, args.out_dir, args.top_n)

    topology_df = pick_top_kegg_pathways(ora, fcs, args.topology_top_n)
    topology_df.to_csv(args.out_dir / "topology_selected_pathways.tsv", sep="\t", index=False)
    run_pathview_subset(args.de_table, topology_df, topology_dir)

    print(f"Wrote HF vs non-HF pathway analyses to: {args.out_dir}")


if __name__ == "__main__":
    main()