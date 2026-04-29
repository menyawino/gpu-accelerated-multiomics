#!/usr/bin/env python3
"""
Arm 2 Validation: Test transcriptome signatures in GSE116250 (DCM/ICM vs NF)

Tests whether:
1. Hyper-down genes are downregulated in DCM vs NF and/or ICM vs NF
2. Hypo-up genes are upregulated in DCM vs NF and/or ICM vs NF
3. Program-level stratification replicates (metabolic + fibrosis/inflammation/ECM)
"""

import argparse
import logging
from pathlib import Path
from typing import Dict, Set, List, Optional

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def parse_sample_phenotypes(metadata_path: Path) -> Dict[str, str]:
    """
    Parse sample metadata to extract phenotype (NF, DCM, ICM, etc).
    
    Assumes TSV with columns: sample_id, phenotype (or similar)
    """
    if not metadata_path.exists():
        logger.warning(f"Metadata file not found: {metadata_path}")
        return {}
    
    try:
        metadata = pd.read_csv(metadata_path, sep="\t")
        
        # Try to find phenotype column
        pheno_col = None
        for col in ["phenotype", "disease_state", "condition", "group"]:
            if col in metadata.columns:
                pheno_col = col
                break
        
        if pheno_col is None:
            logger.warning(f"Could not find phenotype column in {metadata_path}")
            return {}
        
        sample_col = metadata.columns[0]  # Assume first column is sample ID
        result = dict(zip(metadata[sample_col].astype(str), metadata[pheno_col].astype(str)))
        
        logger.info(f"Loaded phenotypes for {len(result)} samples")
        return result
    
    except Exception as e:
        logger.warning(f"Error parsing metadata: {e}")
        return {}


def stratify_samples(expression_matrix: pd.DataFrame,
                     sample_phenotypes: Dict[str, str]) -> Dict[str, pd.DataFrame]:
    """
    Stratify expression matrix by phenotype.
    
    Returns dict: phenotype -> expression matrix (genes x samples of that phenotype)
    """
    stratified = {}
    
    for phenotype in set(sample_phenotypes.values()):
        samples = [s for s, p in sample_phenotypes.items() if p == phenotype]
        valid_samples = [s for s in samples if s in expression_matrix.columns]
        
        if len(valid_samples) > 0:
            stratified[phenotype] = expression_matrix[valid_samples]
            logger.info(f"Phenotype {phenotype}: {len(valid_samples)} samples")
    
    return stratified


def test_gene_direction(signature_genes: List[str],
                       expected_direction: str,
                       control_expr: pd.DataFrame,
                       case_expr: pd.DataFrame) -> pd.DataFrame:
    """
    Test whether signature genes show expected direction in case vs control.
    
    Args:
        signature_genes: Gene IDs to test
        expected_direction: "down" or "up"
        control_expr: Expression matrix for control samples (genes x samples)
        case_expr: Expression matrix for case samples (genes x samples)
    
    Returns:
        DataFrame with test results per gene
    """
    results = []
    
    for gene in signature_genes:
        if gene not in control_expr.index or gene not in case_expr.index:
            continue
        
        ctrl_vals = control_expr.loc[gene].values
        ctrl_vals = ctrl_vals[~np.isnan(ctrl_vals)]
        
        case_vals = case_expr.loc[gene].values
        case_vals = case_vals[~np.isnan(case_vals)]
        
        if len(ctrl_vals) < 2 or len(case_vals) < 2:
            continue
        
        # Mann-Whitney U test
        stat, pval = mannwhitneyu(case_vals, ctrl_vals, alternative='two-sided')
        
        case_mean = np.mean(case_vals)
        ctrl_mean = np.mean(ctrl_vals)
        case_median = np.median(case_vals)
        ctrl_median = np.median(ctrl_vals)
        
        # Determine if direction matches
        if expected_direction == "down":
            direction_match = case_mean < ctrl_mean
        else:  # up
            direction_match = case_mean > ctrl_mean
        
        results.append({
            "gene_id": gene,
            "expected_direction": expected_direction,
            "case_mean": case_mean,
            "control_mean": ctrl_mean,
            "case_median": case_median,
            "control_median": ctrl_median,
            "fold_change_mean": case_mean / ctrl_mean if ctrl_mean != 0 else np.nan,
            "pvalue": pval,
            "direction_match": direction_match,
            "n_case_samples": len(case_vals),
            "n_control_samples": len(ctrl_vals)
        })
    
    return pd.DataFrame(results)


def test_program_stratification(stratified_expr: Dict[str, pd.DataFrame],
                               program_gene_sets: Dict[str, Set[str]]) -> pd.DataFrame:
    """
    Test 3: Program stratification across phenotypes
    
    For each phenotype, compute program scores and test enrichment.
    """
    results = []
    
    for program_name, program_genes in program_gene_sets.items():
        
        if len(program_genes) == 0:
            continue
        
        program_data = []
        
        for phenotype, expr_matrix in stratified_expr.items():
            # Get expression for program genes present in this dataset
            valid_genes = [g for g in program_genes if g in expr_matrix.index]
            
            if len(valid_genes) < 3:
                continue
            
            program_expr = expr_matrix.loc[valid_genes]
            
            # Compute sample-level program scores (mean of log-expr)
            program_scores = program_expr.mean(axis=0)
            
            program_data.append({
                "phenotype": phenotype,
                "program": program_name,
                "n_genes": len(valid_genes),
                "n_samples": len(program_scores),
                "mean_program_score": np.mean(program_scores),
                "median_program_score": np.median(program_scores),
                "std_program_score": np.std(program_scores)
            })
        
        results.extend(program_data)
    
    result_df = pd.DataFrame(results)
    
    if len(result_df) > 0:
        logger.info(f"\nProgram scores:\n{result_df.to_string()}")
    
    return result_df


def load_gene_set(path: Optional[Path]) -> Set[str]:
    """Load one-gene-per-line gene-set file."""
    genes = set()
    if path and path.exists():
        with path.open() as f:
            genes = set(line.strip() for line in f if line.strip())
    return genes


def derive_case_up_modules(validation_expr: pd.DataFrame,
                           stratified_expr: Dict[str, pd.DataFrame],
                           control_pheno: str,
                           case_phenos: List[str],
                           module_size: int = 60) -> Dict[str, Set[str]]:
    """
    Derive fallback program modules from strongest case-up genes if files are not provided.
    """
    if control_pheno not in stratified_expr:
        return {"fibrosis": set(), "inflammation": set(), "ecm_remodeling": set()}

    case_frames = [stratified_expr[p] for p in case_phenos if p in stratified_expr]
    if len(case_frames) == 0:
        return {"fibrosis": set(), "inflammation": set(), "ecm_remodeling": set()}

    control_mean = stratified_expr[control_pheno].mean(axis=1)
    case_concat = pd.concat(case_frames, axis=1)
    case_mean = case_concat.mean(axis=1)
    delta = (case_mean - control_mean).sort_values(ascending=False)

    # Use non-overlapping windows for reproducible fallback modules.
    top = delta.index.tolist()[: module_size * 3]
    return {
        "fibrosis": set(top[0:module_size]),
        "inflammation": set(top[module_size: module_size * 2]),
        "ecm_remodeling": set(top[module_size * 2: module_size * 3]),
    }


def main():
    ap = argparse.ArgumentParser(
        description="Validate transcriptome signatures in GSE116250"
    )
    ap.add_argument("--validation-expression", required=True, type=Path,
                   help="Gene expression matrix from GSE116250 (genes x samples)")
    ap.add_argument("--sample-metadata", required=True, type=Path,
                   help="Sample metadata TSV with phenotype information")
    ap.add_argument("--hyper-down-signature", required=True, type=Path,
                   help="Hyper-down genes from Arm2 discovery")
    ap.add_argument("--hypo-up-signature", required=True, type=Path,
                   help="Hypo-up genes from Arm2 discovery")
    ap.add_argument("--oxidative-genes", type=Path,
                   help="One gene per line: oxidative phosphorylation genes")
    ap.add_argument("--glycolytic-genes", type=Path,
                   help="One gene per line: glycolysis/stress genes")
    ap.add_argument("--fibrosis-genes", type=Path,
                   help="One gene per line: fibrosis-associated genes (optional)")
    ap.add_argument("--inflammation-genes", type=Path,
                   help="One gene per line: inflammation-associated genes (optional)")
    ap.add_argument("--ecm-remodeling-genes", type=Path,
                   help="One gene per line: ECM remodeling-associated genes (optional)")
    ap.add_argument("--control-phenotype", default="NF",
                   help="Control phenotype label (default: NF)")
    ap.add_argument("--case-phenotypes", default="DCM,ICM",
                   help="Case phenotype labels (comma-separated, default: DCM,ICM)")
    ap.add_argument("--out-hyper-down-validation", required=True, type=Path)
    ap.add_argument("--out-hypo-up-validation", required=True, type=Path)
    ap.add_argument("--out-metabolic-programs", required=True, type=Path)
    ap.add_argument("--out-summary", required=True, type=Path)
    
    args = ap.parse_args()
    
    # Load data
    logger.info("Loading validation data...")
    validation_expr = pd.read_csv(args.validation_expression, sep="\t", index_col=0)
    
    hyper_down_sig = pd.read_csv(args.hyper_down_signature, sep="\t")
    hypo_up_sig = pd.read_csv(args.hypo_up_signature, sep="\t")
    
    hyper_down_genes = hyper_down_sig["gene_id"].tolist()
    hypo_up_genes = hypo_up_sig["gene_id"].tolist()
    
    logger.info(f"Validation expression: {validation_expr.shape[0]} genes x {validation_expr.shape[1]} samples")
    logger.info(f"Signatures: {len(hyper_down_genes)} hyper-down, {len(hypo_up_genes)} hypo-up")
    
    # Parse sample phenotypes
    sample_phenotypes = parse_sample_phenotypes(args.sample_metadata)
    
    if not sample_phenotypes:
        logger.warning("Could not parse sample phenotypes; skipping phenotype-stratified analysis")
        summary = {
            "status": "error",
            "message": "No valid sample phenotypes parsed"
        }
        pd.DataFrame([summary]).to_csv(args.out_summary, sep="\t", index=False)
        return
    
    # Stratify samples
    stratified_expr = stratify_samples(validation_expr, sample_phenotypes)
    
    control_pheno = args.control_phenotype
    case_phenos = args.case_phenotypes.split(",")
    
    # Test 1 & 2: Gene-level direction
    logger.info("\n" + "="*70)
    logger.info("Test 1 & 2: Gene-level expression direction")
    logger.info("="*70)
    
    hyper_down_results = []
    hypo_up_results = []
    
    if control_pheno in stratified_expr:
        ctrl_expr = stratified_expr[control_pheno]
        
        for case_pheno in case_phenos:
            if case_pheno not in stratified_expr:
                logger.warning(f"Case phenotype {case_pheno} not found in data")
                continue
            
            case_expr = stratified_expr[case_pheno]
            
            logger.info(f"\nTesting {case_pheno} vs {control_pheno}...")
            
            # Test hyper-down genes (should be DOWN in case)
            hyper_down_test = test_gene_direction(hyper_down_genes, "down", ctrl_expr, case_expr)
            if len(hyper_down_test) > 0:
                hyper_down_test["comparison"] = f"{case_pheno}_vs_{control_pheno}"
                hyper_down_results.append(hyper_down_test)
            
            # Test hypo-up genes (should be UP in case)
            hypo_up_test = test_gene_direction(hypo_up_genes, "up", ctrl_expr, case_expr)
            if len(hypo_up_test) > 0:
                hypo_up_test["comparison"] = f"{case_pheno}_vs_{control_pheno}"
                hypo_up_results.append(hypo_up_test)
    
    if hyper_down_results:
        hyper_down_df = pd.concat(hyper_down_results, ignore_index=True)
        hyper_down_df.to_csv(args.out_hyper_down_validation, sep="\t", index=False)
        logger.info(f"\nHyper-down validation results:")
        match_rate = hyper_down_df["direction_match"].sum() / len(hyper_down_df)
        logger.info(f"  Direction match rate: {match_rate:.2%} ({hyper_down_df['direction_match'].sum()}/{len(hyper_down_df)})")
    else:
        hyper_down_df = pd.DataFrame()
        pd.DataFrame().to_csv(args.out_hyper_down_validation, sep="\t")
    
    if hypo_up_results:
        hypo_up_df = pd.concat(hypo_up_results, ignore_index=True)
        hypo_up_df.to_csv(args.out_hypo_up_validation, sep="\t", index=False)
        logger.info(f"\nHypo-up validation results:")
        match_rate = hypo_up_df["direction_match"].sum() / len(hypo_up_df)
        logger.info(f"  Direction match rate: {match_rate:.2%} ({hypo_up_df['direction_match'].sum()}/{len(hypo_up_df)})")
    else:
        hypo_up_df = pd.DataFrame()
        pd.DataFrame().to_csv(args.out_hypo_up_validation, sep="\t")
    
    # Test 3: Metabolic program stratification
    logger.info("\n" + "="*70)
    logger.info("Test 3: Program stratification")
    logger.info("="*70)

    program_gene_sets = {
        "oxidative": load_gene_set(args.oxidative_genes),
        "glycolytic": load_gene_set(args.glycolytic_genes),
        "fibrosis": load_gene_set(args.fibrosis_genes),
        "inflammation": load_gene_set(args.inflammation_genes),
        "ecm_remodeling": load_gene_set(args.ecm_remodeling_genes),
    }

    fallback = derive_case_up_modules(validation_expr, stratified_expr, control_pheno, case_phenos)
    fallback_used = []
    for key in ["fibrosis", "inflammation", "ecm_remodeling"]:
        if len(program_gene_sets[key]) == 0:
            program_gene_sets[key] = fallback[key]
            if len(program_gene_sets[key]) > 0:
                fallback_used.append(key)

    for name, genes in program_gene_sets.items():
        if len(genes) > 0:
            logger.info(f"Loaded {len(genes)} genes for program: {name}")

    program_scores = test_program_stratification(stratified_expr, program_gene_sets)
    if len(program_scores) > 0:
        program_scores["module_source"] = "provided"
        if len(fallback_used) > 0:
            program_scores.loc[program_scores["program"].isin(fallback_used), "module_source"] = "derived"

    program_scores.to_csv(args.out_metabolic_programs, sep="\t", index=False)
    logger.info(f"Saved to {args.out_metabolic_programs}")
    
    # Summary
    summary = {
        "validation_samples": validation_expr.shape[1],
        "validation_genes": validation_expr.shape[0],
        "phenotypes_found": len(stratified_expr),
        "hyper_down_genes_tested": len(hyper_down_df) if len(hyper_down_df) > 0 else 0,
        "hyper_down_direction_match": hyper_down_df["direction_match"].sum() if len(hyper_down_df) > 0 else 0,
        "hypo_up_genes_tested": len(hypo_up_df) if len(hypo_up_df) > 0 else 0,
        "hypo_up_direction_match": hypo_up_df["direction_match"].sum() if len(hypo_up_df) > 0 else 0,
        "programs_tested": len(program_scores),
        "fibrosis_module_source": "derived" if "fibrosis" in fallback_used else "provided_or_missing",
        "inflammation_module_source": "derived" if "inflammation" in fallback_used else "provided_or_missing",
        "ecm_module_source": "derived" if "ecm_remodeling" in fallback_used else "provided_or_missing"
    }
    
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(args.out_summary, sep="\t", index=False)
    logger.info(f"\n" + "="*70)
    logger.info("VALIDATION SUMMARY")
    logger.info("="*70)
    logger.info(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
