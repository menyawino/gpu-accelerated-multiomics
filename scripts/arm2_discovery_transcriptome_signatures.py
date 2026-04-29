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
from scipy.stats import zscore, pearsonr
from statsmodels.stats.multitest import multipletests

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def identify_expression_changes(gene_expr_matrix: pd.DataFrame, 
                                log2fc_threshold: float = 1.0,
                                n_samples_min: int = 3) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Identify significantly upregulated and downregulated genes.
    
    Uses mean expression and z-score to define DEGs.
    """
    
    # Calculate mean expression per gene
    mean_expr = gene_expr_matrix.mean(axis=1)
    
    # Calculate z-scores across genes and keep gene_id indexing
    z_scores = pd.Series(zscore(mean_expr, nan_policy='omit'), index=mean_expr.index)
    global_mean = float(mean_expr.mean())
    
    # Identify up and downregulated
    upregulated = []
    downregulated = []
    
    for gene_id, idx_row in gene_expr_matrix.iterrows():
        n_valid = (~idx_row.isna()).sum()
        if n_valid < n_samples_min:
            continue
        
        if mean_expr[gene_id] >= (global_mean + log2fc_threshold):
            upregulated.append({
                "gene_id": gene_id,
                "mean_expr": mean_expr[gene_id],
                "z_score": z_scores[gene_id],
                "n_samples": n_valid,
                "direction": "upregulated"
            })
        elif mean_expr[gene_id] <= (global_mean - log2fc_threshold):
            downregulated.append({
                "gene_id": gene_id,
                "mean_expr": mean_expr[gene_id],
                "z_score": z_scores[gene_id],
                "n_samples": n_valid,
                "direction": "downregulated"
            })
    
    up_df = pd.DataFrame(upregulated).sort_values("mean_expr", ascending=False)
    down_df = pd.DataFrame(downregulated).sort_values("mean_expr", ascending=True)
    
    logger.info(f"Found {len(up_df)} upregulated and {len(down_df)} downregulated genes")
    
    return up_df, down_df


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
    ap.add_argument("--hyper-signature", required=True, type=Path,
                   help="Hypermethylated genes (from arm1_discovery)")
    ap.add_argument("--hypo-signature", required=True, type=Path,
                   help="Hypomethylated genes (from arm1_discovery)")
    ap.add_argument("--oxidative-genes", type=Path,
                   help="One gene per line: oxidative phosphorylation genes")
    ap.add_argument("--glycolytic-genes", type=Path,
                   help="One gene per line: glycolysis/stress genes")
    ap.add_argument("--log2fc-threshold", type=float, default=1.0,
                   help="Log2 fold-change threshold for DEGs")
    ap.add_argument("--out-hyper-down", required=True, type=Path)
    ap.add_argument("--out-hypo-up", required=True, type=Path)
    ap.add_argument("--out-summary", required=True, type=Path)
    
    args = ap.parse_args()
    
    # Load data
    logger.info("Loading data...")
    gene_expr = pd.read_csv(args.gene_expression, sep="\t", index_col=0)
    gene_meth = pd.read_csv(args.gene_methylation, sep="\t", index_col=0)
    
    hyper_sig = pd.read_csv(args.hyper_signature, sep="\t")
    hypo_sig = pd.read_csv(args.hypo_signature, sep="\t")
    
    hyper_genes = set(hyper_sig["gene_id"])
    hypo_genes = set(hypo_sig["gene_id"])
    
    logger.info(f"Expression: {gene_expr.shape[0]} genes x {gene_expr.shape[1]} samples")
    logger.info(f"Methylation: {gene_meth.shape[0]} genes x {gene_meth.shape[1]} samples")
    logger.info(f"Signatures: {len(hyper_genes)} hyper, {len(hypo_genes)} hypo")
    
    # Identify expression changes
    logger.info("\nIdentifying expression changes...")
    up_genes_df, down_genes_df = identify_expression_changes(
        gene_expr,
        log2fc_threshold=args.log2fc_threshold
    )
    
    upregulated = set(up_genes_df["gene_id"])
    downregulated = set(down_genes_df["gene_id"])
    
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
    if len(oxidative_genes) > 0 or len(glycolytic_genes) > 0:
        logger.info("\nClassifying metabolic programs...")
        
        hyper_down_metab = classify_metabolic_programs(
            hyper_down_df["gene_id"], oxidative_genes, glycolytic_genes
        )
        hyper_down_df = hyper_down_df.merge(hyper_down_metab, on="gene_id", how="left")
        
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
