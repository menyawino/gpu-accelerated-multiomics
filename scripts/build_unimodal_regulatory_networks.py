#!/usr/bin/env python3
"""Build a unimodal regulatory network from a gene x sample matrix.

The network is inferred from pairwise gene-gene correlations on the top-variance
features, with FDR correction across tested edges.
"""

import argparse
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import t as t_dist
from statsmodels.stats.multitest import multipletests


def select_top_variable_features(df: pd.DataFrame, max_features: int) -> pd.DataFrame:
    if max_features <= 0 or df.shape[0] <= max_features:
        return df
    variances = df.var(axis=1, skipna=True)
    keep = variances.nlargest(max_features).index
    return df.loc[keep]


def build_edges(df: pd.DataFrame, min_abs_corr: float, max_fdr: float, min_samples: int) -> pd.DataFrame:
    # Drop features with too many missing values for robust pairwise estimation.
    keep = (df.notna().sum(axis=1) >= min_samples)
    df = df.loc[keep]
    if df.shape[0] < 2:
        return pd.DataFrame(columns=["source", "target", "correlation", "pvalue", "fdr", "n_samples"])

    # Fill remaining missing values per feature with feature mean.
    centered = df.sub(df.mean(axis=1), axis=0)
    centered = centered.T.fillna(centered.T.mean()).T

    n_samples = centered.shape[1]
    corr = np.corrcoef(centered.values)

    iu = np.triu_indices(corr.shape[0], k=1)
    r = corr[iu]

    # Guard against exact +/-1 values causing division errors.
    r = np.clip(r, -0.999999, 0.999999)
    dof = max(n_samples - 2, 1)
    t_stat = r * np.sqrt(dof / np.maximum(1e-12, 1.0 - r * r))
    pvals = 2.0 * t_dist.sf(np.abs(t_stat), dof)
    fdr = multipletests(pvals, method="fdr_bh")[1]

    genes = centered.index.to_numpy()
    edges = pd.DataFrame(
        {
            "source": genes[iu[0]],
            "target": genes[iu[1]],
            "correlation": r,
            "pvalue": pvals,
            "fdr": fdr,
            "n_samples": n_samples,
        }
    )

    edges = edges[(edges["correlation"].abs() >= min_abs_corr) & (edges["fdr"] <= max_fdr)].copy()
    return edges.sort_values(["fdr", "correlation"], ascending=[True, False])


def write_graph(edges: pd.DataFrame, graphml_path: Path, omic_label: str) -> None:
    g = nx.Graph()
    for _, row in edges.iterrows():
        s = row["source"]
        t = row["target"]
        g.add_node(s, omic=omic_label)
        g.add_node(t, omic=omic_label)
        g.add_edge(
            s,
            t,
            weight=float(row["correlation"]),
            correlation=float(row["correlation"]),
            fdr=float(row["fdr"]),
            pvalue=float(row["pvalue"]),
            n_samples=int(row["n_samples"]),
        )
    nx.write_graphml(g, graphml_path)


def main() -> None:
    ap = argparse.ArgumentParser(description="Build unimodal (RNA or methylation) regulatory network")
    ap.add_argument("--matrix", required=True, type=Path, help="Gene x sample matrix TSV")
    ap.add_argument("--omic-label", required=True, choices=["rna", "methylation"])
    ap.add_argument("--max-features", type=int, default=1000, help="Top variable genes to use")
    ap.add_argument("--min-abs-correlation", type=float, default=0.6)
    ap.add_argument("--max-fdr", type=float, default=0.05)
    ap.add_argument("--min-samples", type=int, default=8)
    ap.add_argument("--out-edges", required=True, type=Path)
    ap.add_argument("--out-graphml", required=True, type=Path)
    ap.add_argument("--out-summary", required=True, type=Path)
    args = ap.parse_args()

    mat = pd.read_csv(args.matrix, sep="\t", index_col=0)
    mat = select_top_variable_features(mat, args.max_features)

    edges = build_edges(mat, args.min_abs_correlation, args.max_fdr, args.min_samples)

    args.out_edges.parent.mkdir(parents=True, exist_ok=True)
    args.out_graphml.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)

    edges.to_csv(args.out_edges, sep="\t", index=False)
    write_graph(edges, args.out_graphml, args.omic_label)

    summary = pd.DataFrame(
        [
            {
                "omic": args.omic_label,
                "features_input": pd.read_csv(args.matrix, sep="\t", index_col=0).shape[0],
                "features_used": mat.shape[0],
                "samples_used": mat.shape[1],
                "edges_retained": edges.shape[0],
                "min_abs_correlation": args.min_abs_correlation,
                "max_fdr": args.max_fdr,
            }
        ]
    )
    summary.to_csv(args.out_summary, sep="\t", index=False)


if __name__ == "__main__":
    main()
