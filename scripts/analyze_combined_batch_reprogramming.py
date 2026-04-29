#!/usr/bin/env python3
"""Combined all-dataset reprogramming analysis with batch adjustment and consistency checks."""

import argparse
from pathlib import Path
from typing import Dict, Set

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.api import OLS
from statsmodels.stats.multitest import multipletests


def load_gene_set(path: Path | None) -> Set[str]:
    if path is None or not path.exists():
        return set()
    with path.open() as f:
        return {line.strip() for line in f if line.strip()}


def normalize_labels(series: pd.Series) -> pd.Series:
    s = series.astype(str).str.upper().str.strip()
    out = pd.Series(index=s.index, dtype="object")
    out[s.isin(["NF", "CONTROL", "NONFAILING", "NON_FAILING"])] = "nonfailing"
    out[s.isin(["HF", "FAIL", "FAILING", "DCM", "ICM", "HFR"])] = "failing"
    return out


def infer_labels_from_metabolic_balance(expr: pd.DataFrame, oxidative: Set[str], glycolytic: Set[str]) -> pd.Series:
    ox = [g for g in oxidative if g in expr.index]
    gl = [g for g in glycolytic if g in expr.index]
    if len(ox) < 10 or len(gl) < 10:
        # Fallback: split by first principal direction over genes using SVD on centered matrix
        mat = expr.T.values
        mat = mat - mat.mean(axis=0, keepdims=True)
        u, s, _ = np.linalg.svd(mat, full_matrices=False)
        pc1 = u[:, 0] * s[0]
        thr = np.median(pc1)
        labels = np.where(pc1 > thr, "failing", "nonfailing")
        return pd.Series(labels, index=expr.columns)

    ox_score = expr.loc[ox].mean(axis=0)
    gl_score = expr.loc[gl].mean(axis=0)
    balance = gl_score - ox_score
    thr = np.median(balance.values)
    labels = np.where(balance.values > thr, "failing", "nonfailing")
    return pd.Series(labels, index=balance.index)


def dataset_gene_zscore(expr: pd.DataFrame) -> pd.DataFrame:
    mu = expr.mean(axis=1)
    sd = expr.std(axis=1).replace(0, np.nan)
    return expr.sub(mu, axis=0).div(sd, axis=0).fillna(0.0)


def derive_modules(effect_df: pd.DataFrame, module_size: int = 250) -> Dict[str, Set[str]]:
    ranked = effect_df.sort_values("combined_effect", ascending=False)["gene_id"].tolist()
    top = ranked[: module_size * 3]
    return {
        "fibrosis": set(top[:module_size]),
        "inflammation": set(top[module_size : 2 * module_size]),
        "ecm_remodeling": set(top[2 * module_size : 3 * module_size]),
    }


def summarize_pathway(pathway: str, genes: Set[str], effect_df: pd.DataFrame) -> dict:
    sub = effect_df[effect_df["gene_id"].isin(genes)].copy()
    if len(sub) == 0:
        return {
            "pathway": pathway,
            "genes_tested": 0,
            "median_combined_effect": np.nan,
            "fraction_positive_effect": np.nan,
            "combined_effect_p_median": np.nan,
            "consistency_fraction": np.nan,
        }

    return {
        "pathway": pathway,
        "genes_tested": int(len(sub)),
        "median_combined_effect": float(np.median(sub["combined_effect"])),
        "fraction_positive_effect": float((sub["combined_effect"] > 0).mean()),
        "combined_effect_p_median": float(np.median(sub["combined_pvalue"])),
        "consistency_fraction": float(np.mean(sub["sign_consistent"])),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Combined all-dataset batch-aware reprogramming analysis")
    ap.add_argument("--expr-116250", required=True, type=Path)
    ap.add_argument("--meta-116250", required=True, type=Path)
    ap.add_argument("--expr-123976", required=True, type=Path)
    ap.add_argument("--meta-123976", type=Path)
    ap.add_argument("--oxidative-genes", required=True, type=Path)
    ap.add_argument("--glycolytic-genes", required=True, type=Path)
    ap.add_argument("--fibrosis-genes", type=Path)
    ap.add_argument("--inflammation-genes", type=Path)
    ap.add_argument("--ecm-genes", type=Path)
    ap.add_argument("--module-size", type=int, default=250)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    ox = load_gene_set(args.oxidative_genes)
    gl = load_gene_set(args.glycolytic_genes)

    e116 = pd.read_csv(args.expr_116250, sep="\t", index_col=0)
    m116 = pd.read_csv(args.meta_116250, sep="\t")

    e123 = pd.read_csv(args.expr_123976, sep="\t", index_col=0)

    # Align genes across datasets.
    common_genes = sorted(set(e116.index).intersection(set(e123.index)))
    e116 = e116.loc[common_genes]
    e123 = e123.loc[common_genes]

    # Known labels for 116250.
    pheno_116 = normalize_labels(pd.Series(m116.iloc[:, 1].values, index=m116.iloc[:, 0].astype(str).values))
    pheno_116 = pheno_116.dropna()
    cols116 = [c for c in e116.columns if c in pheno_116.index]
    e116 = e116[cols116]
    pheno_116 = pheno_116.loc[cols116]

    # Optional labels for 123976, otherwise infer biologically.
    if args.meta_123976 and args.meta_123976.exists():
        m123 = pd.read_csv(args.meta_123976, sep="\t")
        pheno_123 = normalize_labels(pd.Series(m123.iloc[:, 1].values, index=m123.iloc[:, 0].astype(str).values))
        pheno_123 = pheno_123.reindex(e123.columns)
        if pheno_123.notna().sum() < 4:
            pheno_123 = infer_labels_from_metabolic_balance(e123, ox, gl)
            label_source_123 = "inferred"
        else:
            pheno_123 = pheno_123.dropna()
            e123 = e123[pheno_123.index]
            label_source_123 = "provided"
    else:
        pheno_123 = infer_labels_from_metabolic_balance(e123, ox, gl)
        label_source_123 = "inferred"

    # Batch-aware normalization: per-dataset gene z-scoring.
    z116 = dataset_gene_zscore(e116)
    z123 = dataset_gene_zscore(e123)

    # Build combined matrix and design.
    combined = pd.concat([z116, z123], axis=1)
    pheno = pd.concat([
        pd.DataFrame({"sample": z116.columns, "dataset": "GSE116250", "status": pheno_116.loc[z116.columns].values}),
        pd.DataFrame({"sample": z123.columns, "dataset": "GSE123976", "status": pheno_123.loc[z123.columns].values}),
    ], ignore_index=True)

    pheno = pheno.set_index("sample").loc[combined.columns]
    failing = (pheno["status"] == "failing").astype(float).values
    batch = (pheno["dataset"] == "GSE123976").astype(float).values
    X = np.column_stack([np.ones(len(failing)), failing, batch])

    # Per-gene combined model y ~ failing + batch
    rows = []
    for gene in combined.index:
        y = combined.loc[gene].values.astype(float)
        model = OLS(y, X).fit()
        rows.append({
            "gene_id": gene,
            "combined_effect": float(model.params[1]),
            "combined_pvalue": float(model.pvalues[1]),
        })

    effect = pd.DataFrame(rows)
    effect["combined_fdr"] = multipletests(effect["combined_pvalue"].values, method="fdr_bh")[1]

    # Per-dataset effect consistency checks.
    per_rows = []
    for dataset_name, zmat, labels in [
        ("GSE116250", z116, pheno_116.loc[z116.columns]),
        ("GSE123976", z123, pheno_123.loc[z123.columns]),
    ]:
        fail_cols = [c for c in zmat.columns if labels.loc[c] == "failing"]
        nf_cols = [c for c in zmat.columns if labels.loc[c] == "nonfailing"]
        for gene in zmat.index:
            fv = zmat.loc[gene, fail_cols].values.astype(float)
            nv = zmat.loc[gene, nf_cols].values.astype(float)
            if len(fv) < 2 or len(nv) < 2:
                continue
            p = mannwhitneyu(fv, nv, alternative="two-sided").pvalue
            per_rows.append({
                "dataset": dataset_name,
                "gene_id": gene,
                "effect": float(np.mean(fv) - np.mean(nv)),
                "pvalue": float(p),
                "n_failing": int(len(fv)),
                "n_nonfailing": int(len(nv)),
            })

    per_df = pd.DataFrame(per_rows)
    per_df["fdr"] = np.nan
    for ds in per_df["dataset"].unique():
        idx = per_df["dataset"] == ds
        per_df.loc[idx, "fdr"] = multipletests(per_df.loc[idx, "pvalue"].values, method="fdr_bh")[1]

    pivot_eff = per_df.pivot(index="gene_id", columns="dataset", values="effect")
    pivot_sign = np.sign(pivot_eff)
    sign_consistent = (pivot_sign.nunique(axis=1) == 1).rename("sign_consistent")

    effect = effect.merge(sign_consistent, on="gene_id", how="left")
    effect["sign_consistent"] = effect["sign_consistent"].fillna(False)

    # Pathway modules.
    fibrosis = load_gene_set(args.fibrosis_genes)
    inflammation = load_gene_set(args.inflammation_genes)
    ecm = load_gene_set(args.ecm_genes)

    source = {
        "fibrosis": "provided" if len(fibrosis) else "derived",
        "inflammation": "provided" if len(inflammation) else "derived",
        "ecm_remodeling": "provided" if len(ecm) else "derived",
    }

    if len(fibrosis) == 0 or len(inflammation) == 0 or len(ecm) == 0:
        derived = derive_modules(effect, module_size=args.module_size)
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

    path_rows = []
    for p, genes in pathways.items():
        for g in sorted(genes):
            path_rows.append({"gene_id": g, "pathway": p, "module_source": source[p]})

    path_genes = pd.DataFrame(path_rows)
    path_effect = effect.merge(path_genes, on="gene_id", how="inner")

    summary = pd.DataFrame([summarize_pathway(k, v, effect) for k, v in pathways.items()])
    summary["module_source"] = summary["pathway"].map(source)

    # Run-level summary.
    run_summary = pd.DataFrame(
        [
            {
                "datasets": "GSE116250,GSE123976",
                "gse123976_label_source": label_source_123,
                "samples_116250": z116.shape[1],
                "samples_123976": z123.shape[1],
                "genes_tested": effect.shape[0],
                "genes_fdr_lt_0.05": int((effect["combined_fdr"] < 0.05).sum()),
                "consistent_sign_fraction": float(effect["sign_consistent"].mean()),
            }
        ]
    )

    effect.sort_values(["combined_fdr", "combined_effect"], ascending=[True, False]).to_csv(
        args.out_dir / "combined_all_genes_effects.tsv", sep="\t", index=False
    )
    per_df.to_csv(args.out_dir / "per_dataset_gene_effects.tsv", sep="\t", index=False)
    path_genes.to_csv(args.out_dir / "combined_pathway_gene_membership.tsv", sep="\t", index=False)
    path_effect.to_csv(args.out_dir / "combined_pathway_gene_effects.tsv", sep="\t", index=False)
    summary.to_csv(args.out_dir / "combined_pathway_summary.tsv", sep="\t", index=False)
    run_summary.to_csv(args.out_dir / "combined_run_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
