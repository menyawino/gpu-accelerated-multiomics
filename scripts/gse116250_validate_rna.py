#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from processed_analysis_common import differential_expression, read_matrix


def running_enrichment_score(ranked_genes: list[str], gene_set: set[str]) -> float:
    hits = np.array([g in gene_set for g in ranked_genes], dtype=float)
    if hits.sum() == 0:
        return 0.0
    miss = 1.0 - hits
    hit_step = hits / hits.sum()
    miss_step = miss / miss.sum() if miss.sum() else miss
    running = np.cumsum(hit_step - miss_step)
    return float(running[np.argmax(np.abs(running))])


def gsea_like(contrast: str, de: pd.DataFrame, sets: dict[str, set[str]], permutations: int, seed: int) -> pd.DataFrame:
    de = de.sort_values("log2fc", ascending=False)
    ranked = de["gene_id"].astype(str).tolist()
    rng = np.random.default_rng(seed)
    rows = []
    for set_name, genes in sets.items():
        genes = genes.intersection(set(ranked))
        if not genes:
            continue
        observed = running_enrichment_score(ranked, genes)
        null = []
        for _ in range(permutations):
            perm = ranked.copy()
            rng.shuffle(perm)
            null.append(running_enrichment_score(perm, genes))
        pvalue = (sum(abs(x) >= abs(observed) for x in null) + 1) / (len(null) + 1)
        rows.append((contrast, set_name, len(genes), observed, pvalue))
    out = pd.DataFrame(rows, columns=["contrast", "gene_set", "n_genes", "enrichment_score", "pvalue"])
    if out.empty:
        out["fdr"] = []
        return out
    from processed_analysis_common import adjust_pvalues

    out["fdr"] = adjust_pvalues(out["pvalue"])
    return out.sort_values(["contrast", "fdr", "gene_set"])


def sample_signature_scores(expr: pd.DataFrame, sets: dict[str, set[str]]) -> pd.DataFrame:
    zexpr = expr.apply(lambda col: (col - col.mean()) / (col.std(ddof=0) + 1e-9), axis=0)
    rows = []
    genes_available = set(zexpr.index.astype(str))
    for set_name, genes in sets.items():
        valid = [g for g in genes if g in genes_available]
        if not valid:
            continue
        vals = zexpr.loc[valid].mean(axis=0)
        for sample_id, score in vals.items():
            rows.append((sample_id, set_name, float(score)))
    return pd.DataFrame(rows, columns=["sample_id", "signature", "score"])


def concordance_rate(discovery: pd.DataFrame, de: pd.DataFrame, class_name: str, expected_sign: int) -> float:
    genes = discovery.loc[discovery["integrated_class"] == class_name, "gene_id"].astype(str)
    sub = de[de["gene_id"].astype(str).isin(genes)]
    if sub.empty:
        return np.nan
    return float((np.sign(sub["log2fc"]) == expected_sign).mean())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna-matrix", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--discovery-integrated", required=True)
    ap.add_argument("--control-label", required=True)
    ap.add_argument("--dcm-label", required=True)
    ap.add_argument("--icm-label", required=True)
    ap.add_argument("--fdr-threshold", type=float, required=True)
    ap.add_argument("--log2fc-threshold", type=float, required=True)
    ap.add_argument("--permutations", type=int, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--min-group-size", type=int, required=True)
    ap.add_argument("--out-dcm-deg", required=True)
    ap.add_argument("--out-icm-deg", required=True)
    ap.add_argument("--out-gsea", required=True)
    ap.add_argument("--out-validation-stats", required=True)
    ap.add_argument("--out-signature-scores", required=True)
    ap.add_argument("--out-plot", required=True)
    args = ap.parse_args()

    expr = read_matrix(args.rna_matrix, index_col="gene_id")
    metadata = pd.read_csv(args.metadata, sep="\t")
    discovery = pd.read_csv(args.discovery_integrated, sep="\t")

    dcm_deg = differential_expression(expr, metadata, args.dcm_label, args.control_label, args.min_group_size)
    icm_deg = differential_expression(expr, metadata, args.icm_label, args.control_label, args.min_group_size)

    for de in (dcm_deg, icm_deg):
        de["significant"] = (de["fdr"] <= args.fdr_threshold) & (de["log2fc"].abs() >= args.log2fc_threshold)

    sets = {
        "Hyper-down": set(discovery.loc[discovery["integrated_class"] == "Hyper-down", "gene_id"].astype(str)),
        "Hypo-up": set(discovery.loc[discovery["integrated_class"] == "Hypo-up", "gene_id"].astype(str)),
    }

    gsea = pd.concat(
        [
            gsea_like("DCM_vs_NF", dcm_deg, sets, args.permutations, args.seed),
            gsea_like("ICM_vs_NF", icm_deg, sets, args.permutations, args.seed + 1),
        ],
        ignore_index=True,
    )

    stats = pd.DataFrame(
        {
            "contrast": ["DCM_vs_NF", "ICM_vs_NF"],
            "n_deg_significant": [int(dcm_deg["significant"].sum()), int(icm_deg["significant"].sum())],
            "hyper_down_concordance": [
                concordance_rate(discovery, dcm_deg, "Hyper-down", expected_sign=-1),
                concordance_rate(discovery, icm_deg, "Hyper-down", expected_sign=-1),
            ],
            "hypo_up_concordance": [
                concordance_rate(discovery, dcm_deg, "Hypo-up", expected_sign=1),
                concordance_rate(discovery, icm_deg, "Hypo-up", expected_sign=1),
            ],
        }
    )

    scores = sample_signature_scores(np.log2(expr + 1.0), sets)
    scores = scores.merge(metadata[["sample_id", "group"]], on="sample_id", how="left")

    score_plot = (
        scores.groupby(["group", "signature"], as_index=False)["score"].mean()
        if not scores.empty
        else pd.DataFrame({"group": [], "signature": [], "score": []})
    )
    plt.figure(figsize=(7, 4))
    if not score_plot.empty:
        x = np.arange(score_plot["group"].nunique())
        groups = sorted(score_plot["group"].dropna().unique())
        width = 0.35
        for i, signature in enumerate(sorted(score_plot["signature"].unique())):
            vals = [
                score_plot.loc[
                    (score_plot["group"] == g) & (score_plot["signature"] == signature),
                    "score",
                ].mean()
                for g in groups
            ]
            plt.bar(x + (i - 0.5) * width, vals, width=width, label=signature)
        plt.xticks(x, groups)
        plt.ylabel("Mean z-score")
        plt.legend(frameon=False)
    plt.title("GSE116250 signature score validation")
    plt.tight_layout()
    Path(args.out_plot).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out_plot, dpi=150)
    plt.close()

    dcm_deg.to_csv(args.out_dcm_deg, sep="\t", index=False)
    icm_deg.to_csv(args.out_icm_deg, sep="\t", index=False)
    gsea.to_csv(args.out_gsea, sep="\t", index=False)
    stats.to_csv(args.out_validation_stats, sep="\t", index=False)
    scores.to_csv(args.out_signature_scores, sep="\t", index=False)


if __name__ == "__main__":
    main()
