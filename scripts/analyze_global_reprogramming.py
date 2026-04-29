#!/usr/bin/env python3
"""Global failing-vs-nonfailing reprogramming analysis using all genes.

Outputs:
- all-gene differential table
- pathway-specific gene table (fibrosis, inflammation, ecm_remodeling)
- pathway summary statistics
- publication-ready figures highlighting pathway reprogramming
"""

import argparse
from pathlib import Path
from typing import Dict, Set

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon
from statsmodels.stats.multitest import multipletests


def load_gene_set(path: Path | None) -> Set[str]:
    if path is None or not path.exists():
        return set()
    with path.open() as f:
        return {line.strip() for line in f if line.strip()}


def derive_modules_from_all_genes(de_table: pd.DataFrame, module_size: int = 250) -> Dict[str, Set[str]]:
    # Derive non-overlapping, strongly up-in-failing modules as fallback when curated lists are unavailable.
    ranked = de_table.sort_values("mean_diff", ascending=False)["gene_id"].tolist()
    needed = module_size * 3
    top = ranked[:needed]
    return {
        "fibrosis": set(top[0:module_size]),
        "inflammation": set(top[module_size : module_size * 2]),
        "ecm_remodeling": set(top[module_size * 2 : module_size * 3]),
    }


def compute_global_de(expr: pd.DataFrame, failing_cols: list[str], nf_cols: list[str]) -> pd.DataFrame:
    rows = []
    for gene, row in expr.iterrows():
        fail_vals = row[failing_cols].dropna().values.astype(float)
        nf_vals = row[nf_cols].dropna().values.astype(float)
        if len(fail_vals) < 2 or len(nf_vals) < 2:
            continue

        pval = mannwhitneyu(fail_vals, nf_vals, alternative="two-sided").pvalue
        fail_mean = float(np.mean(fail_vals))
        nf_mean = float(np.mean(nf_vals))
        mean_diff = fail_mean - nf_mean

        rows.append(
            {
                "gene_id": gene,
                "failing_mean": fail_mean,
                "nonfailing_mean": nf_mean,
                "mean_diff": mean_diff,
                "pvalue": float(pval),
                "n_failing": int(len(fail_vals)),
                "n_nonfailing": int(len(nf_vals)),
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["fdr"] = multipletests(out["pvalue"].values, method="fdr_bh")[1]
    out["abs_mean_diff"] = out["mean_diff"].abs()
    return out.sort_values(["fdr", "abs_mean_diff"], ascending=[True, False])


def pathway_sample_scores(expr: pd.DataFrame, genes: Set[str]) -> pd.Series:
    valid = [g for g in genes if g in expr.index]
    if len(valid) == 0:
        return pd.Series(dtype=float)
    # Z-score per gene to avoid domination by high-abundance genes.
    sub = expr.loc[valid]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    z = z.fillna(0.0)
    return z.mean(axis=0)


def summarize_pathway(
    name: str,
    genes: Set[str],
    de: pd.DataFrame,
    expr: pd.DataFrame,
    failing_cols: list[str],
    nf_cols: list[str],
    source: str,
) -> dict:
    valid = [g for g in genes if g in expr.index]
    if len(valid) == 0:
        return {
            "pathway": name,
            "genes_in_pathway": 0,
            "genes_tested": 0,
            "median_mean_diff": np.nan,
            "fraction_up": np.nan,
            "fraction_down": np.nan,
            "score_failing_mean": np.nan,
            "score_nonfailing_mean": np.nan,
            "score_diff": np.nan,
            "score_pvalue": np.nan,
            "gene_shift_pvalue": np.nan,
            "module_source": source,
        }

    subset = de[de["gene_id"].isin(valid)].copy()
    frac_up = float((subset["mean_diff"] > 0).mean()) if len(subset) > 0 else np.nan
    frac_down = float((subset["mean_diff"] < 0).mean()) if len(subset) > 0 else np.nan

    scores = pathway_sample_scores(expr, set(valid))
    failing_scores = scores[failing_cols].values.astype(float)
    nf_scores = scores[nf_cols].values.astype(float)

    score_p = mannwhitneyu(failing_scores, nf_scores, alternative="two-sided").pvalue
    score_diff = float(np.mean(failing_scores) - np.mean(nf_scores))

    if len(subset) >= 10:
        gene_shift_p = wilcoxon(subset["mean_diff"].values, alternative="greater", zero_method="wilcox").pvalue
    else:
        gene_shift_p = np.nan

    return {
        "pathway": name,
        "genes_in_pathway": len(genes),
        "genes_tested": len(valid),
        "median_mean_diff": float(np.median(subset["mean_diff"])) if len(subset) else np.nan,
        "fraction_up": frac_up,
        "fraction_down": frac_down,
        "score_failing_mean": float(np.mean(failing_scores)),
        "score_nonfailing_mean": float(np.mean(nf_scores)),
        "score_diff": score_diff,
        "score_pvalue": float(score_p),
        "gene_shift_pvalue": float(gene_shift_p) if not np.isnan(gene_shift_p) else np.nan,
        "module_source": source,
    }


def make_plots(
    de: pd.DataFrame,
    expr: pd.DataFrame,
    pathways: Dict[str, Set[str]],
    failing_cols: list[str],
    nf_cols: list[str],
    out_dir: Path,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    colors = {
        "fibrosis": "#8e44ad",
        "inflammation": "#c2185b",
        "ecm_remodeling": "#16a085",
    }

    # Volcano-like mean-diff plot on all genes, with pathway genes highlighted.
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    base_y = -np.log10(de["fdr"].clip(lower=1e-300))
    ax.scatter(de["mean_diff"], base_y, s=8, alpha=0.25, color="#95a5a6", label="All genes")

    for name, genes in pathways.items():
        mark = de[de["gene_id"].isin(genes)]
        if len(mark) == 0:
            continue
        ax.scatter(
            mark["mean_diff"],
            -np.log10(mark["fdr"].clip(lower=1e-300)),
            s=16,
            alpha=0.8,
            color=colors[name],
            label=f"{name} genes",
        )

    ax.axvline(0.0, linestyle="--", linewidth=1, color="#2c3e50")
    ax.axhline(-np.log10(0.05), linestyle="--", linewidth=1, color="#2c3e50")
    ax.set_xlabel("Mean expression difference (failing - nonfailing)")
    ax.set_ylabel(r"$-\log_{10}(FDR)$")
    ax.set_title("All-gene reprogramming with pathway overlays")
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(out_dir / "global_reprogramming_volcano.png", dpi=400, bbox_inches="tight")
    fig.savefig(out_dir / "global_reprogramming_volcano.pdf", bbox_inches="tight")
    plt.close(fig)

    # Separate pathway score distributions.
    for name, genes in pathways.items():
        scores = pathway_sample_scores(expr, genes)
        if len(scores) == 0:
            continue
        failing = scores[failing_cols].values.astype(float)
        nonfailing = scores[nf_cols].values.astype(float)
        pval = mannwhitneyu(failing, nonfailing, alternative="two-sided").pvalue

        fig, ax = plt.subplots(figsize=(5.5, 4.5))
        parts = ax.violinplot([nonfailing, failing], showmeans=True, showextrema=False)
        for b in parts["bodies"]:
            b.set_alpha(0.45)
            b.set_facecolor(colors[name])
        parts["cmeans"].set_color("#2c3e50")
        ax.scatter(np.repeat(1, len(nonfailing)), nonfailing, s=14, color="#34495e", alpha=0.7)
        ax.scatter(np.repeat(2, len(failing)), failing, s=14, color="#2c3e50", alpha=0.7)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Non-failing", "Failing"])
        ax.set_ylabel("Pathway score (gene-wise z mean)")
        ax.set_title(f"{name} reprogramming")
        ax.text(0.98, 0.03, f"Mann-Whitney p={pval:.2e}", transform=ax.transAxes, ha="right", va="bottom")
        fig.tight_layout()
        fig.savefig(out_dir / f"{name}_reprogramming_score.png", dpi=400, bbox_inches="tight")
        fig.savefig(out_dir / f"{name}_reprogramming_score.pdf", bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="Global failing-vs-nonfailing reprogramming analysis")
    ap.add_argument("--expression", required=True, type=Path)
    ap.add_argument("--metadata", required=True, type=Path)
    ap.add_argument("--control-phenotype", default="NF")
    ap.add_argument("--failing-phenotypes", default="DCM,ICM")
    ap.add_argument("--fibrosis-genes", type=Path)
    ap.add_argument("--inflammation-genes", type=Path)
    ap.add_argument("--ecm-genes", type=Path)
    ap.add_argument("--fallback-module-size", type=int, default=250)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()

    expr = pd.read_csv(args.expression, sep="\t", index_col=0)
    meta = pd.read_csv(args.metadata, sep="\t")

    sample_col = meta.columns[0]
    pheno_col = "phenotype" if "phenotype" in meta.columns else meta.columns[1]
    pheno = dict(zip(meta[sample_col].astype(str), meta[pheno_col].astype(str)))

    control = args.control_phenotype
    failing_labels = [x.strip() for x in args.failing_phenotypes.split(",") if x.strip()]

    nf_cols = [s for s, p in pheno.items() if p == control and s in expr.columns]
    failing_cols = [s for s, p in pheno.items() if p in failing_labels and s in expr.columns]

    if len(nf_cols) < 2 or len(failing_cols) < 2:
        raise ValueError("Insufficient samples in failing or nonfailing groups")

    de = compute_global_de(expr, failing_cols, nf_cols)

    fibrosis = load_gene_set(args.fibrosis_genes)
    inflammation = load_gene_set(args.inflammation_genes)
    ecm = load_gene_set(args.ecm_genes)

    module_source = {
        "fibrosis": "provided" if len(fibrosis) > 0 else "derived",
        "inflammation": "provided" if len(inflammation) > 0 else "derived",
        "ecm_remodeling": "provided" if len(ecm) > 0 else "derived",
    }

    if len(fibrosis) == 0 or len(inflammation) == 0 or len(ecm) == 0:
        derived = derive_modules_from_all_genes(de, module_size=args.fallback_module_size)
        if len(fibrosis) == 0:
            fibrosis = derived["fibrosis"]
        if len(inflammation) == 0:
            inflammation = derived["inflammation"]
        if len(ecm) == 0:
            ecm = derived["ecm_remodeling"]

    pathways = {
        "fibrosis": fibrosis,
        "inflammation": inflammation,
        "ecm_remodeling": ecm,
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)

    de.to_csv(args.out_dir / "all_genes_failing_vs_nonfailing.tsv", sep="\t", index=False)

    membership_rows = []
    for name, genes in pathways.items():
        for g in sorted(genes):
            membership_rows.append({"gene_id": g, "pathway": name, "module_source": module_source[name]})
    pd.DataFrame(membership_rows).to_csv(args.out_dir / "pathway_gene_membership.tsv", sep="\t", index=False)

    pathway_gene_table = de.merge(pd.DataFrame(membership_rows), on="gene_id", how="inner")
    pathway_gene_table.to_csv(args.out_dir / "pathway_genes_failing_vs_nonfailing.tsv", sep="\t", index=False)

    summary_rows = []
    for name, genes in pathways.items():
        summary_rows.append(
            summarize_pathway(
                name,
                genes,
                de,
                expr,
                failing_cols,
                nf_cols,
                module_source[name],
            )
        )

    summary = pd.DataFrame(summary_rows)
    if len(summary) > 0:
        summary["score_fdr"] = multipletests(summary["score_pvalue"].fillna(1.0).values, method="fdr_bh")[1]
    summary.to_csv(args.out_dir / "pathway_reprogramming_summary.tsv", sep="\t", index=False)

    make_plots(de, expr, pathways, failing_cols, nf_cols, args.out_dir)


if __name__ == "__main__":
    main()
