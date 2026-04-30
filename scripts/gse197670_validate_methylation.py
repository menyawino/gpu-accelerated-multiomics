#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, ttest_ind

from processed_analysis_common import adjust_pvalues, read_matrix


def pick_contrast(metadata: pd.DataFrame, control_label: str, disease_labels: list[str], min_group_size: int) -> tuple[str, list[str], list[str]]:
    controls = metadata.loc[metadata["group"] == control_label, "sample_id"].astype(str).tolist()
    for disease in disease_labels:
        cases = metadata.loc[metadata["group"] == disease, "sample_id"].astype(str).tolist()
        if len(cases) >= min_group_size and len(controls) >= min_group_size:
            return disease, cases, controls
    raise ValueError("No disease/control contrast with enough samples for methylation validation.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta-matrix", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--probe-annotation", required=True)
    ap.add_argument("--discovery-integrated", required=True)
    ap.add_argument("--control-label", required=True)
    ap.add_argument("--disease-labels", required=True)
    ap.add_argument("--fdr-threshold", type=float, required=True)
    ap.add_argument("--permutations", type=int, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--min-group-size", type=int, required=True)
    ap.add_argument("--out-dmp", required=True)
    ap.add_argument("--out-gene-promoter", required=True)
    ap.add_argument("--out-concordance", required=True)
    ap.add_argument("--out-plot", required=True)
    ap.add_argument("--out-qc", required=True)
    args = ap.parse_args()

    metadata = pd.read_csv(args.metadata, sep="\t")
    beta = read_matrix(args.beta_matrix, index_col="probe_id")
    annot = pd.read_csv(args.probe_annotation, sep="\t")
    annot["is_promoter"] = (
        annot["is_promoter"]
        .astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1", "yes", "y"})
    )

    disease_label, cases, controls = pick_contrast(
        metadata,
        args.control_label,
        [x.strip() for x in args.disease_labels.split(",") if x.strip()],
        args.min_group_size,
    )

    shared = [c for c in cases + controls if c in beta.columns]
    beta = beta[shared]
    cases = [c for c in cases if c in beta.columns]
    controls = [c for c in controls if c in beta.columns]

    case_vals = beta[cases].astype(float)
    control_vals = beta[controls].astype(float)
    stat, pvals = ttest_ind(case_vals.values, control_vals.values, axis=1, equal_var=False, nan_policy="omit")
    dmp = pd.DataFrame(
        {
            "probe_id": beta.index.astype(str),
            "contrast": f"{disease_label}_vs_{args.control_label}",
            "delta_beta": case_vals.mean(axis=1).values - control_vals.mean(axis=1).values,
            "statistic": stat,
            "pvalue": pvals,
        }
    )
    dmp["fdr"] = adjust_pvalues(dmp["pvalue"])
    dmp["direction"] = np.where(dmp["delta_beta"] >= 0, "Hyper", "Hypo")
    dmp["significant"] = dmp["fdr"] <= args.fdr_threshold

    merged = dmp.merge(annot[["probe_id", "gene_id", "is_promoter"]], on="probe_id", how="left")
    promoter = merged[merged["is_promoter"]].copy()
    gene_summary = (
        promoter.groupby("gene_id", as_index=False)
        .agg(
            promoter_delta_beta=("delta_beta", "mean"),
            promoter_min_fdr=("fdr", "min"),
            n_promoter_probes=("probe_id", "nunique"),
        )
        .sort_values(["promoter_min_fdr", "gene_id"])
    )
    gene_summary["promoter_direction"] = np.where(gene_summary["promoter_delta_beta"] >= 0, "Hyper", "Hypo")

    discovery = pd.read_csv(args.discovery_integrated, sep="\t")
    discovery_sets = discovery[discovery["integrated_class"].isin(["Hyper-down", "Hypo-up"])][
        ["gene_id", "integrated_class"]
    ]
    concordance = discovery_sets.merge(gene_summary, on="gene_id", how="left")
    concordance["expected_direction"] = np.where(
        concordance["integrated_class"] == "Hyper-down", "Hyper", "Hypo"
    )
    concordance = concordance.dropna(subset=["promoter_direction"])
    concordance["is_concordant"] = concordance["expected_direction"] == concordance["promoter_direction"]

    rng = np.random.default_rng(args.seed)
    observed = concordance["is_concordant"].mean() if not concordance.empty else np.nan
    null_vals = []
    if not concordance.empty:
        dirs = concordance["promoter_direction"].to_numpy()
        expected = concordance["expected_direction"].to_numpy()
        for _ in range(args.permutations):
            shuffled = rng.permutation(dirs)
            null_vals.append(float(np.mean(shuffled == expected)))
    perm_p = (
        (sum(v >= observed for v in null_vals) + 1) / (len(null_vals) + 1)
        if null_vals and np.isfinite(observed)
        else np.nan
    )
    sign_p = (
        binomtest(int(concordance["is_concordant"].sum()), n=int(concordance.shape[0]), p=0.5, alternative="greater").pvalue
        if not concordance.empty
        else np.nan
    )

    out_concordance = pd.DataFrame(
        {
            "contrast": [f"{disease_label}_vs_{args.control_label}"],
            "n_genes_tested": [int(concordance.shape[0])],
            "n_concordant": [int(concordance["is_concordant"].sum()) if not concordance.empty else 0],
            "concordance_rate": [float(observed) if np.isfinite(observed) else np.nan],
            "sign_test_pvalue": [sign_p],
            "permutation_pvalue": [perm_p],
        }
    )

    qc = pd.DataFrame(
        {
            "metric": ["n_probes", "n_promoter_probes", "n_genes_with_promoter_summary"],
            "value": [int(dmp.shape[0]), int(promoter["probe_id"].nunique()), int(gene_summary.shape[0])],
        }
    )

    plt.figure(figsize=(5, 4))
    class_rates = (
        concordance.groupby("integrated_class", as_index=False)["is_concordant"].mean()
        if not concordance.empty
        else pd.DataFrame({"integrated_class": ["Hyper-down", "Hypo-up"], "is_concordant": [np.nan, np.nan]})
    )
    plt.bar(class_rates["integrated_class"], class_rates["is_concordant"], color=["#2166ac", "#b2182b"])
    plt.ylim(0, 1)
    plt.ylabel("Concordance rate")
    plt.title("GSE197670 promoter methylation concordance")
    plt.tight_layout()
    Path(args.out_plot).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out_plot, dpi=150)
    plt.close()

    dmp.sort_values(["fdr", "probe_id"]).to_csv(args.out_dmp, sep="\t", index=False)
    gene_summary.to_csv(args.out_gene_promoter, sep="\t", index=False)
    out_concordance.to_csv(args.out_concordance, sep="\t", index=False)
    qc.to_csv(args.out_qc, sep="\t", index=False)


if __name__ == "__main__":
    main()
