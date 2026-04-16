#!/usr/bin/env python3
import argparse

import networkx as nx
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from statsmodels.stats.multitest import multipletests


def align_by_pairs(df: pd.DataFrame, ordered_runs: list[str]) -> pd.DataFrame:
    cols = [c for c in ordered_runs if c in df.columns]
    return df[cols] if cols else pd.DataFrame(index=df.index)


def correlate_rows(meth: pd.DataFrame, expr: pd.DataFrame, method: str, min_pairs: int):
    shared_features = sorted(set(meth.index).intersection(expr.index))
    results = []

    for feat in shared_features:
        x = meth.loc[feat]
        y = expr.loc[feat]
        common = pd.concat([x, y], axis=1, join="inner").dropna()
        if common.shape[0] < min_pairs:
            continue

        if method == "pearson":
            r, p = pearsonr(common.iloc[:, 0], common.iloc[:, 1])
        else:
            r, p = spearmanr(common.iloc[:, 0], common.iloc[:, 1])

        results.append((feat, float(r), float(p), int(common.shape[0])))

    out = pd.DataFrame(results, columns=["feature_id", "correlation", "pvalue", "n_pairs"])
    if out.empty:
        return out

    out["fdr"] = multipletests(out["pvalue"].values, method="fdr_bh")[1]
    return out.sort_values(["fdr", "correlation"], ascending=[True, False])


def write_graph(edges_df: pd.DataFrame, out_graphml: str, layer: str):
    g = nx.Graph()
    for _, row in edges_df.iterrows():
        feat = row["feature_id"]
        methyl_node = f"methyl::{layer}::{feat}"
        expr_node = f"expr::{layer}::{feat}"
        g.add_node(methyl_node, omic="methylation", level=layer, feature=feat)
        g.add_node(expr_node, omic="expression", level=layer, feature=feat)
        g.add_edge(
            methyl_node,
            expr_node,
            weight=float(row["correlation"]),
            fdr=float(row["fdr"]),
            n_pairs=int(row["n_pairs"]),
        )
    nx.write_graphml(g, out_graphml)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-expression", required=True)
    ap.add_argument("--iso-expression", required=True)
    ap.add_argument("--gene-methylation", required=True)
    ap.add_argument("--iso-methylation", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--method", choices=["pearson", "spearman"], default="spearman")
    ap.add_argument("--min-abs-correlation", type=float, default=0.35)
    ap.add_argument("--max-fdr", type=float, default=0.05)
    ap.add_argument("--min-pairs", type=int, default=6)
    ap.add_argument("--out-gene-edges", required=True)
    ap.add_argument("--out-iso-edges", required=True)
    ap.add_argument("--out-gene-graphml", required=True)
    ap.add_argument("--out-iso-graphml", required=True)
    args = ap.parse_args()

    pairs = pd.read_csv(args.pairs, sep="\t") if args.pairs else pd.DataFrame()
    if not pairs.empty and {"rnaseq_run", "wgbs_run"}.issubset(set(pairs.columns)):
        # Order by pair rows and relabel WGBS to RNA IDs for direct column alignment.
        ordered_rna = pairs["rnaseq_run"].astype(str).tolist()
        wgbs_to_rna = dict(zip(pairs["wgbs_run"].astype(str), ordered_rna))
    else:
        ordered_rna = []
        wgbs_to_rna = {}

    gene_expr = pd.read_csv(args.gene_expression, sep="\t", index_col=0)
    iso_expr = pd.read_csv(args.iso_expression, sep="\t", index_col=0)
    gene_meth = pd.read_csv(args.gene_methylation, sep="\t", index_col=0)
    iso_meth = pd.read_csv(args.iso_methylation, sep="\t", index_col=0)

    if wgbs_to_rna:
        gene_meth = gene_meth.rename(columns=wgbs_to_rna)
        iso_meth = iso_meth.rename(columns=wgbs_to_rna)

    if ordered_rna:
        gene_expr = align_by_pairs(gene_expr, ordered_rna)
        iso_expr = align_by_pairs(iso_expr, ordered_rna)
        gene_meth = align_by_pairs(gene_meth, ordered_rna)
        iso_meth = align_by_pairs(iso_meth, ordered_rna)

    gene_edges = correlate_rows(gene_meth, gene_expr, args.method, args.min_pairs)
    iso_edges = correlate_rows(iso_meth, iso_expr, args.method, args.min_pairs)

    if not gene_edges.empty:
        gene_edges = gene_edges[
            (gene_edges["fdr"] <= args.max_fdr)
            & (gene_edges["correlation"].abs() >= args.min_abs_correlation)
        ]
    if not iso_edges.empty:
        iso_edges = iso_edges[
            (iso_edges["fdr"] <= args.max_fdr)
            & (iso_edges["correlation"].abs() >= args.min_abs_correlation)
        ]

    gene_edges.to_csv(args.out_gene_edges, sep="\t", index=False)
    iso_edges.to_csv(args.out_iso_edges, sep="\t", index=False)

    write_graph(gene_edges, args.out_gene_graphml, layer="gene")
    write_graph(iso_edges, args.out_iso_graphml, layer="isoform")


if __name__ == "__main__":
    main()
