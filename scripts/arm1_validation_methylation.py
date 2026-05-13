#!/usr/bin/env python3
"""
Arm 1 Validation: Test methylation signatures in GSE197670

Tests whether:
1. Same genes show concordant promoter methylation direction (HF vs NF)
2. Pathway-level methylation shifts replicate
3. Metabolic gene stratification replicates (oxidative hyper, glycolytic hypo)

Platform note: GSE197670 includes both GPL13534 (Illumina 450K) and GPL21145
(Illumina EPIC 850K) arrays, which differ in CpG coverage.  When the input
methylation matrix has a leading '__platform__' index row (added by
prepare_real_processed_geo_inputs.py), this script extracts it and reports
concordance stratified by platform as well as combined, so cross-platform
consistency of signatures can be assessed.
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, Optional, Set

import pandas as pd
import numpy as np
from scipy.stats import binomtest, mannwhitneyu
from statsmodels.stats.multitest import multipletests

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def load_signatures(hyper_path: Path, hypo_path: Path) -> tuple:
    """Load discovery signatures."""
    hyper = pd.read_csv(hyper_path, sep="\t")
    hypo = pd.read_csv(hypo_path, sep="\t")
    return hyper, hypo


def _strip_platform_row(meth: pd.DataFrame) -> tuple[pd.DataFrame, Optional[pd.Series]]:
    """Remove the '__platform__' metadata row if present; return (matrix, platform_series)."""
    if "__platform__" in meth.index:
        platform_info = meth.loc["__platform__"].copy()
        meth = meth.drop(index="__platform__")
        meth = meth.apply(pd.to_numeric, errors="coerce")
        return meth, platform_info
    return meth, None


def test_concordance(
    hyper_discovery: pd.DataFrame,
    hypo_discovery: pd.DataFrame,
    validation_meth: pd.DataFrame,
    platform_info: Optional[pd.Series] = None,
) -> pd.DataFrame:
    """
    Test 1: Concordance of methylation direction with FDR correction.

    For each signature gene present in the validation cohort, perform a one-sample
    Wilcoxon signed-rank test vs. 0.5 to determine whether methylation deviates in
    the expected direction.  The 0.5 threshold is the neutral midpoint for beta
    values; hypermethylated genes should be > 0.5, hypomethylated < 0.5.

    If *platform_info* is provided (a Series mapping sample IDs to platform IDs),
    concordance is computed per platform as well as combined.
    """
    results = []

    for sig_type, sig_df, expected_above_half in [
        ("hypermethylated", hyper_discovery, True),
        ("hypomethylated", hypo_discovery, False),
    ]:
        for gene in sig_df["gene_id"]:
            if gene not in validation_meth.index:
                continue

            val_beta = validation_meth.loc[gene].values.astype(float)
            val_beta = val_beta[~np.isnan(val_beta)]
            if len(val_beta) == 0:
                continue

            val_mean = float(np.nanmean(val_beta))

            # Wilcoxon signed-rank vs. 0.5 as the neutral methylation reference.
            from scipy.stats import wilcoxon as _wilcoxon
            if len(val_beta) >= 2 and not np.all(val_beta == val_beta[0]):
                alt = "greater" if expected_above_half else "less"
                try:
                    _, direction_pval = _wilcoxon(val_beta - 0.5, alternative=alt)
                except ValueError:
                    direction_pval = np.nan
            else:
                direction_pval = np.nan

            concordant = (val_mean > 0.5) if expected_above_half else (val_mean < 0.5)

            row: dict = {
                "gene_id": gene,
                "signature_direction": sig_type,
                "validation_mean_beta": val_mean,
                "validation_median_beta": float(np.nanmedian(val_beta)),
                "validation_n_samples": int((~np.isnan(validation_meth.loc[gene].values)).sum()),
                "direction_pvalue": direction_pval,
                "concordant": concordant,
                "platform": "combined",
            }
            results.append(row)

            # Per-platform breakdown when platform metadata is available.
            if platform_info is not None:
                for plat in platform_info.unique():
                    plat_cols = [c for c in validation_meth.columns if platform_info.get(c) == plat]
                    if not plat_cols:
                        continue
                    plat_betas = validation_meth.loc[gene, plat_cols].values.astype(float)
                    plat_betas = plat_betas[~np.isnan(plat_betas)]
                    if len(plat_betas) == 0:
                        continue
                    plat_mean = float(np.nanmean(plat_betas))
                    results.append({
                        **row,
                        "validation_mean_beta": plat_mean,
                        "validation_median_beta": float(np.nanmedian(plat_betas)),
                        "validation_n_samples": len(plat_betas),
                        "concordant": (plat_mean > 0.5) if expected_above_half else (plat_mean < 0.5),
                        "platform": str(plat),
                        "direction_pvalue": np.nan,  # Underpowered at per-platform level; report combined.
                    })

    concordance_df = pd.DataFrame(results)

    if len(concordance_df) > 0:
        # FDR correct within the combined-platform rows only.
        combined_mask = concordance_df["platform"] == "combined"
        combined_pvals = concordance_df.loc[combined_mask, "direction_pvalue"].fillna(1.0).values
        if len(combined_pvals) > 0:
            concordance_df.loc[combined_mask, "direction_fdr"] = multipletests(
                combined_pvals, method="fdr_bh"
            )[1]
        concordance_rate = concordance_df.loc[combined_mask, "concordant"].mean()
        logger.info(
            f"Combined concordance rate: {concordance_rate:.2%} "
            f"({concordance_df.loc[combined_mask, 'concordant'].sum()}/{combined_mask.sum()})"
        )

    return concordance_df


def test_pathway_replication(hyper_discovery: pd.DataFrame,
                             hypo_discovery: pd.DataFrame,
                             validation_meth: pd.DataFrame) -> pd.DataFrame:
    """
    Test 2: Pathway-level methylation shift replication
    
    For genes annotated with pathways, check if the pathway-level signature
    (mean methylation trend) is preserved.
    """
    results = []
    
    for sig_type, genes_df in [("hypermethylated", hyper_discovery), 
                               ("hypomethylated", hypo_discovery)]:
        
        if "pathways" not in genes_df.columns:
            continue
        
        # Group by pathway
        pathway_genes = {}
        for _, row in genes_df.iterrows():
            pathways = str(row.get("pathways", "")).split(",")
            for pw in pathways:
                pw = pw.strip()
                if pw and pw != "unclassified":
                    if pw not in pathway_genes:
                        pathway_genes[pw] = []
                    pathway_genes[pw].append(row["gene_id"])
        
        # Test each pathway
        for pathway, genes in pathway_genes.items():
            valid_genes = [g for g in genes if g in validation_meth.index]
            if len(valid_genes) < 3:  # Need at least 3 genes
                continue
            
            pathway_betas = validation_meth.loc[valid_genes].values.flatten()
            pathway_betas = pathway_betas[~np.isnan(pathway_betas)]
            
            if len(pathway_betas) == 0:
                continue
            
            result = {
                "pathway": pathway,
                "signature_direction": sig_type,
                "n_genes_in_pathway": len(genes),
                "n_genes_validated": len(valid_genes),
                "mean_beta_validation": np.mean(pathway_betas),
                "median_beta_validation": np.median(pathway_betas),
                "expected_direction": "high" if sig_type == "hypermethylated" else "low",
                "direction_match": (np.median(pathway_betas) > 0.5) if sig_type == "hypermethylated" else (np.median(pathway_betas) < 0.5)
            }
            results.append(result)
    
    pathway_df = pd.DataFrame(results)
    
    if len(pathway_df) > 0:
        match_rate = pathway_df["direction_match"].sum() / len(pathway_df)
        logger.info(f"Pathway direction match rate: {match_rate:.2%}")
    
    return pathway_df


def test_metabolic_stratification(validation_meth: pd.DataFrame,
                                  oxidative_genes: Set[str],
                                  glycolytic_genes: Set[str]) -> pd.DataFrame:
    """
    Test 3: Metabolic gene stratification
    
    Oxidative metabolism genes should skew toward hypermethylation.
    Glycolytic/stress response genes should skew toward hypomethylation.
    """
    results = []
    
    # Test oxidative genes
    oxidative_present = [g for g in oxidative_genes if g in validation_meth.index]
    if len(oxidative_present) >= 3:
        ox_betas = validation_meth.loc[oxidative_present].values.flatten()
        ox_betas = ox_betas[~np.isnan(ox_betas)]
        
        ox_hyper_count = (ox_betas > 0.6).sum()
        ox_hypo_count = (ox_betas < 0.4).sum()
        ox_neutral_count = len(ox_betas) - ox_hyper_count - ox_hypo_count
        
        # Test if oxidative genes are skewed toward hyper
        p_value_ox = binomtest(ox_hyper_count, len(ox_betas), 0.5, alternative='greater').pvalue
        
        results.append({
            "metabolic_program": "oxidative_phosphorylation",
            "n_genes": len(oxidative_present),
            "n_cpgs": len(ox_betas),
            "hypermethylated_count": ox_hyper_count,
            "hypomethylated_count": ox_hypo_count,
            "neutral_count": ox_neutral_count,
            "mean_beta": np.mean(ox_betas),
            "median_beta": np.median(ox_betas),
            "expected_skew": "hypermethylated",
            "p_value": p_value_ox,
            "skew_confirmed": p_value_ox < 0.05 and np.mean(ox_betas) > 0.5
        })
    
    # Test glycolytic genes
    glycolytic_present = [g for g in glycolytic_genes if g in validation_meth.index]
    if len(glycolytic_present) >= 3:
        glc_betas = validation_meth.loc[glycolytic_present].values.flatten()
        glc_betas = glc_betas[~np.isnan(glc_betas)]
        
        glc_hyper_count = (glc_betas > 0.6).sum()
        glc_hypo_count = (glc_betas < 0.4).sum()
        glc_neutral_count = len(glc_betas) - glc_hyper_count - glc_hypo_count
        
        # Test if glycolytic genes are skewed toward hypo
        p_value_glc = binomtest(glc_hypo_count, len(glc_betas), 0.5, alternative='greater').pvalue
        
        results.append({
            "metabolic_program": "glycolysis",
            "n_genes": len(glycolytic_present),
            "n_cpgs": len(glc_betas),
            "hypermethylated_count": glc_hyper_count,
            "hypomethylated_count": glc_hypo_count,
            "neutral_count": glc_neutral_count,
            "mean_beta": np.mean(glc_betas),
            "median_beta": np.median(glc_betas),
            "expected_skew": "hypomethylated",
            "p_value": p_value_glc,
            "skew_confirmed": p_value_glc < 0.05 and np.mean(glc_betas) < 0.5
        })
    
    metabolic_df = pd.DataFrame(results)
    
    if len(metabolic_df) > 0:
        confirmed = metabolic_df["skew_confirmed"].sum()
        logger.info(f"Metabolic stratification confirmed: {confirmed}/{len(metabolic_df)}")
    
    return metabolic_df


def main():
    ap = argparse.ArgumentParser(
        description="Validate methylation signatures in GSE197670"
    )
    ap.add_argument("--hyper-signature", required=True, type=Path,
                   help="Hypermethylated genes from Arm1 discovery")
    ap.add_argument("--hypo-signature", required=True, type=Path,
                   help="Hypomethylated genes from Arm1 discovery")
    ap.add_argument("--validation-methylation", required=True, type=Path,
                   help="Gene methylation matrix from GSE197670")
    ap.add_argument("--oxidative-genes", type=Path,
                   help="One gene per line: oxidative phosphorylation genes")
    ap.add_argument("--glycolytic-genes", type=Path,
                   help="One gene per line: glycolysis/stress genes")
    ap.add_argument("--out-concordance", required=True, type=Path)
    ap.add_argument("--out-pathway-replication", required=True, type=Path)
    ap.add_argument("--out-metabolic-stratification", required=True, type=Path)
    ap.add_argument("--out-summary", required=True, type=Path)
    
    args = ap.parse_args()
    
    logger.info("Loading signatures and validation data...")
    hyper, hypo = load_signatures(args.hyper_signature, args.hypo_signature)
    validation_meth_raw = pd.read_csv(args.validation_methylation, sep="\t", index_col=0)

    # Strip platform metadata row if present (added by prepare_real_processed_geo_inputs.py).
    validation_meth, platform_info = _strip_platform_row(validation_meth_raw)

    logger.info(f"Hyper genes: {len(hyper)}, Hypo genes: {len(hypo)}")
    logger.info(f"Validation matrix: {validation_meth.shape[0]} genes x {validation_meth.shape[1]} samples")
    if platform_info is not None:
        plat_counts = platform_info.value_counts().to_dict()
        logger.info(f"Platforms in validation data: {plat_counts}")

    # Test 1: Concordance (with FDR correction and per-platform breakdown)
    logger.info("\nTest 1: Testing concordance of methylation direction...")
    concordance = test_concordance(hyper, hypo, validation_meth, platform_info=platform_info)
    concordance.to_csv(args.out_concordance, sep="\t", index=False)
    logger.info(f"Saved to {args.out_concordance}")
    
    # Test 2: Pathway replication
    logger.info("\nTest 2: Testing pathway-level replication...")
    pathway_rep = test_pathway_replication(hyper, hypo, validation_meth)
    pathway_rep.to_csv(args.out_pathway_replication, sep="\t", index=False)
    logger.info(f"Saved to {args.out_pathway_replication}")
    
    # Test 3: Metabolic stratification
    logger.info("\nTest 3: Testing metabolic gene stratification...")
    oxidative_genes = set()
    glycolytic_genes = set()
    
    if args.oxidative_genes and args.oxidative_genes.exists():
        with args.oxidative_genes.open() as f:
            oxidative_genes = set(line.strip() for line in f if line.strip())
        logger.info(f"Loaded {len(oxidative_genes)} oxidative genes")
    
    if args.glycolytic_genes and args.glycolytic_genes.exists():
        with args.glycolytic_genes.open() as f:
            glycolytic_genes = set(line.strip() for line in f if line.strip())
        logger.info(f"Loaded {len(glycolytic_genes)} glycolytic genes")
    
    metabolic_strat = test_metabolic_stratification(validation_meth, oxidative_genes, glycolytic_genes)
    metabolic_strat.to_csv(args.out_metabolic_stratification, sep="\t", index=False)
    logger.info(f"Saved to {args.out_metabolic_stratification}")
    
    # Summary
    summary = {
        "test": ["concordance", "pathway_replication", "metabolic_stratification"],
        "genes_tested": [len(concordance), len(pathway_rep), len(metabolic_strat)],
        "test_passed": [
            concordance["concordant"].sum() if len(concordance) > 0 else 0,
            pathway_rep["direction_match"].sum() if len(pathway_rep) > 0 else 0,
            metabolic_strat["skew_confirmed"].sum() if len(metabolic_strat) > 0 else 0
        ]
    }
    
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(args.out_summary, sep="\t", index=False)
    logger.info(f"\nValidation summary:\n{summary_df.to_string()}")


if __name__ == "__main__":
    main()
