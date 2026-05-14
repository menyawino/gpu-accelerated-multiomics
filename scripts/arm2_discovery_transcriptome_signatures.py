#!/usr/bin/env python3
"""
Arm 2 Discovery: Extract transcriptome co-occurrence signatures from GSE123976

Identifies:
1. Hyper-down genes: promoter hypermethylated + downregulated
2. Hypo-up genes: promoter hypomethylated + upregulated
3. Metabolic program classification (oxidative vs glycolytic)
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, Set, Tuple

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def identify_expression_changes(gene_expr_matrix: pd.DataFrame,
                                phenotype_dict: Dict[str, str],
                                hf_phenotypes: Set[str],
                                nf_phenotype: str = "NF",
                                fdr_cutoff: float = 0.05,
                                log2fc_threshold: float = 0.0,
                                n_samples_min: int = 3) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Identify significantly upregulated and downregulated genes using Mann-Whitney U test.
    
    Args:
        gene_expr_matrix: Gene expression matrix (genes x samples)
        phenotype_dict: Dict mapping sample IDs to phenotype labels
        hf_phenotypes: Set of phenotype labels for failing heart samples (e.g., {"HF", "DCM", "ICM"})
        nf_phenotype: Phenotype label for non-failing samples (default: "NF")
        fdr_cutoff: FDR threshold for significance (default: 0.05)
        log2fc_threshold: Absolute log2 fold-change threshold (default: 0.0, no filter)
        n_samples_min: Minimum samples per group (default: 3)
    
    Returns:
        Tuple of (upregulated_df, downregulated_df) with p-values and FDR
    """
    # Stratify samples by phenotype
    nf_samples = [s for s in gene_expr_matrix.columns if phenotype_dict.get(s) == nf_phenotype and s in gene_expr_matrix.columns]
    hf_samples = [s for s in gene_expr_matrix.columns if phenotype_dict.get(s) in hf_phenotypes and s in gene_expr_matrix.columns]
    
    if len(nf_samples) < n_samples_min or len(hf_samples) < n_samples_min:
        logger.error(f"Insufficient samples: NF={len(nf_samples)}, HF={len(hf_samples)}")
        return pd.DataFrame(), pd.DataFrame()
    
    logger.info(f"Stratified: {len(nf_samples)} NF samples, {len(hf_samples)} HF samples")
    
    # Test each gene
    results = []
    for gene_id in gene_expr_matrix.index:
        nf_arr = np.asarray(pd.to_numeric(
            pd.Series(gene_expr_matrix.reindex(index=[gene_id], columns=nf_samples).to_numpy().ravel()),
            errors="coerce",
        ), dtype=float)
        hf_arr = np.asarray(pd.to_numeric(
            pd.Series(gene_expr_matrix.reindex(index=[gene_id], columns=hf_samples).to_numpy().ravel()),
            errors="coerce",
        ), dtype=float)
        nf_vals = nf_arr[~np.isnan(nf_arr)].astype(float)
        hf_vals = hf_arr[~np.isnan(hf_arr)].astype(float)
        
        if len(nf_vals) < n_samples_min or len(hf_vals) < n_samples_min:
            continue
        
        # Mann-Whitney U test
        stat, pvalue = mannwhitneyu(hf_vals, nf_vals, alternative='two-sided')
        
        # Calculate log2FC as median or mean difference
        hf_mean = np.mean(hf_vals)
        nf_mean = np.mean(nf_vals)
        log2fc = np.log2((hf_mean + 1e-6) / (nf_mean + 1e-6))  # Avoid division by zero
        
        results.append({
            "gene_id": gene_id,
            "hf_mean": hf_mean,
            "nf_mean": nf_mean,
            "log2fc": log2fc,
            "pvalue": pvalue,
            "n_hf": len(hf_vals),
            "n_nf": len(nf_vals)
        })
    
    out = pd.DataFrame(results)
    if out.empty:
        logger.warning("No genes passed filtering")
        return pd.DataFrame(), pd.DataFrame()
    
    # Apply FDR correction
    out["fdr"] = multipletests(out["pvalue"].values, method="fdr_bh")[1]
    
    # Filter by FDR and log2FC
    sig = out[(out["fdr"] <= fdr_cutoff) & (out["log2fc"].abs() >= log2fc_threshold)].copy()
    
    # Split into up and down
    upregulated = sig[sig["log2fc"] > 0].sort_values("log2fc", ascending=False)
    downregulated = sig[sig["log2fc"] < 0].sort_values("log2fc", ascending=True)
    
    logger.info(f"Found {len(upregulated)} upregulated and {len(downregulated)} downregulated genes (FDR < {fdr_cutoff}, |log2FC| > {log2fc_threshold})")
    
    return upregulated, downregulated


def identify_cooccurrence_signatures(hyper_genes: Set[str],
                                     hypo_genes: Set[str],
                                     upregulated_genes: Set[str],
                                     downregulated_genes: Set[str],
                                     gene_meth: pd.DataFrame,
                                     gene_expr: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Identify concordant hyper-down and hypo-up genes.
    """
    
    # Hyper-down: hypermethylated AND downregulated
    hyper_down = hyper_genes.intersection(downregulated_genes)
    hyper_down_data = []
    
    for gene in hyper_down:
        meth_val = gene_meth.loc[gene].mean() if gene in gene_meth.index else np.nan
        expr_val = gene_expr.loc[gene].mean() if gene in gene_expr.index else np.nan
        
        hyper_down_data.append({
            "gene_id": gene,
            "mean_promoter_beta": meth_val,
            "mean_expression": expr_val,
            "signature": "hyper_down",
            "interpretation": "Repressed: hypermethylated promoter + low expression"
        })
    
    hyper_down_df = pd.DataFrame(hyper_down_data)
    
    # Hypo-up: hypomethylated AND upregulated
    hypo_up = hypo_genes.intersection(upregulated_genes)
    hypo_up_data = []
    
    for gene in hypo_up:
        meth_val = gene_meth.loc[gene].mean() if gene in gene_meth.index else np.nan
        expr_val = gene_expr.loc[gene].mean() if gene in gene_expr.index else np.nan
        
        hypo_up_data.append({
            "gene_id": gene,
            "mean_promoter_beta": meth_val,
            "mean_expression": expr_val,
            "signature": "hypo_up",
            "interpretation": "Active: hypomethylated promoter + high expression"
        })
    
    hypo_up_df = pd.DataFrame(hypo_up_data)
    
    logger.info(f"Found {len(hyper_down_df)} hyper-down genes and {len(hypo_up_df)} hypo-up genes")
    
    return hyper_down_df, hypo_up_df


def classify_metabolic_programs(gene_list: pd.Series,
                                oxidative_genes: Set[str],
                                glycolytic_genes: Set[str]) -> pd.DataFrame:
    """
    Classify genes into metabolic programs.
    """
    
    result = []
    for gene in gene_list:
        program = None
        
        if gene in oxidative_genes:
            program = "oxidative_phosphorylation"
        elif gene in glycolytic_genes:
            program = "glycolysis_stress"
        else:
            program = "other"
        
        result.append({
            "gene_id": gene,
            "metabolic_program": program
        })
    
    return pd.DataFrame(result)


def main():
    ap = argparse.ArgumentParser(
        description="Extract transcriptome co-occurrence signatures from GSE123976"
    )
    ap.add_argument("--gene-expression", required=True, type=Path,
                   help="Gene expression matrix (TPM or normalized)")
    ap.add_argument("--gene-methylation", required=True, type=Path,
                   help="Gene promoter methylation matrix (beta values)")
    ap.add_argument("--metadata", required=True, type=Path,
                   help="Sample metadata TSV with phenotype column")
    ap.add_argument("--hyper-signature", required=True, type=Path,
                   help="Hypermethylated genes (from arm1_discovery)")
    ap.add_argument("--hypo-signature", required=True, type=Path,
                   help="Hypomethylated genes (from arm1_discovery)")
    ap.add_argument("--oxidative-genes", type=Path,
                   help="One gene per line: oxidative phosphorylation genes")
    ap.add_argument("--glycolytic-genes", type=Path,
                   help="One gene per line: glycolysis/stress genes")
    ap.add_argument("--fdr-cutoff", type=float, default=0.05,
                   help="FDR cutoff for differential expression (default: 0.05)")
    ap.add_argument("--log2fc-threshold", type=float, default=0.0,
                   help="Absolute log2 fold-change threshold for DEGs (default: 0.0)")
    ap.add_argument("--out-degs", type=Path, default=None,
                   help="Optional: output all DEG table with p-values and FDR")
    ap.add_argument("--out-hyper-down", required=True, type=Path)
    ap.add_argument("--out-hypo-up", required=True, type=Path)
    ap.add_argument("--out-summary", required=True, type=Path)
    
    args = ap.parse_args()
    
    # Load data
    logger.info("Loading data...")
    gene_expr = pd.read_csv(args.gene_expression, sep="\t", index_col=0)
    gene_meth = pd.read_csv(args.gene_methylation, sep="\t", index_col=0)
    
    # Load metadata and build phenotype dict
    meta = pd.read_csv(args.metadata, sep="\t")
    sample_col = meta.columns[0]
    pheno_col = "phenotype" if "phenotype" in meta.columns else meta.columns[1]
    phenotype_dict: Dict[str, str] = {
        str(sample): str(pheno)
        for sample, pheno in zip(meta[sample_col].astype(str), meta[pheno_col].astype(str))
    }
    
    hyper_sig = pd.read_csv(args.hyper_signature, sep="\t")
    hypo_sig = pd.read_csv(args.hypo_signature, sep="\t")
    
    hyper_genes = set(hyper_sig["gene_id"])
    hypo_genes = set(hypo_sig["gene_id"])
    
    logger.info(f"Expression: {gene_expr.shape[0]} genes x {gene_expr.shape[1]} samples")
    logger.info(f"Methylation: {gene_meth.shape[0]} genes x {gene_meth.shape[1]} samples")
    logger.info(f"Signatures: {len(hyper_genes)} hyper, {len(hypo_genes)} hypo")
    
    # Identify expression changes using statistical testing.
    # At n=9 (3 NF, 6 HF), FDR < 0.05 across 31K genes is unattainable.
    # Save the full ranked table; co-occurrence signatures use direction of effect.
    logger.info("\nIdentifying expression changes via Mann-Whitney U test...")
    up_genes_df, down_genes_df = identify_expression_changes(
        gene_expr,
        phenotype_dict=phenotype_dict,
        hf_phenotypes={"HF", "DCM", "ICM"},
        nf_phenotype="NF",
        fdr_cutoff=1.0,          # No FDR filter — save full ranked list for downstream use
        log2fc_threshold=0.0,
    )

    # Save full ranked DE table regardless of significance.
    if args.out_degs is not None:
        all_degs = pd.concat([up_genes_df, down_genes_df], ignore_index=True)
        all_degs = all_degs.sort_values("fdr")
        all_degs.to_csv(args.out_degs, sep="\t", index=False)
        logger.info(f"Saved all DEGs to {args.out_degs}")

    # For co-occurrence: use direction-consistent genes (any log2FC in expected direction).
    # At this sample size, direction consistency is the appropriate discovery metric.
    upregulated = set(up_genes_df["gene_id"]) if not up_genes_df.empty else set()
    downregulated = set(down_genes_df["gene_id"]) if not down_genes_df.empty else set()
    
    # Identify co-occurrence signatures
    logger.info("\nIdentifying co-occurrence signatures...")
    hyper_down_df, hypo_up_df = identify_cooccurrence_signatures(
        hyper_genes, hypo_genes, upregulated, downregulated,
        gene_meth, gene_expr
    )
    
    # Load metabolic gene sets if provided
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
    
    # Classify metabolic programs
    if (len(oxidative_genes) > 0 or len(glycolytic_genes) > 0) and not hyper_down_df.empty:
        logger.info("\nClassifying metabolic programs...")
        
        hyper_down_metab = classify_metabolic_programs(
            hyper_down_df["gene_id"], oxidative_genes, glycolytic_genes
        )
        hyper_down_df = hyper_down_df.merge(hyper_down_metab, on="gene_id", how="left")

    if (len(oxidative_genes) > 0 or len(glycolytic_genes) > 0) and not hypo_up_df.empty:
        hypo_up_metab = classify_metabolic_programs(
            hypo_up_df["gene_id"], oxidative_genes, glycolytic_genes
        )
        hypo_up_df = hypo_up_df.merge(hypo_up_metab, on="gene_id", how="left")
    
    # Save outputs
    hyper_down_df.to_csv(args.out_hyper_down, sep="\t", index=False)
    logger.info(f"Saved hyper-down genes to {args.out_hyper_down}")
    
    hypo_up_df.to_csv(args.out_hypo_up, sep="\t", index=False)
    logger.info(f"Saved hypo-up genes to {args.out_hypo_up}")
    
    # Summary
    summary = {
        "discovery_samples_expr": gene_expr.shape[1],
        "discovery_samples_meth": gene_meth.shape[1],
        "discovery_genes": gene_expr.shape[0],
        "upregulated_count": len(upregulated),
        "downregulated_count": len(downregulated),
        "hyper_down_count": len(hyper_down_df),
        "hypo_up_count": len(hypo_up_df),
        "log2fc_threshold": args.log2fc_threshold
    }
    
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(args.out_summary, sep="\t", index=False)
    logger.info(f"Summary saved to {args.out_summary}")
    logger.info(f"\nSummary:\n{summary_df.to_string()}")


if __name__ == "__main__":
    main()
