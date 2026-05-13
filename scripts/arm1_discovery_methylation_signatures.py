#!/usr/bin/env python3
"""
Arm 1 Discovery: Extract methylation signatures from GSE123976

Derives three signature types:
1. Promoter hypermethylated genes
2. Promoter hypomethylated genes
3. Pathway-linked differentially methylated regions (DMRs) with gene associations

These are used for validation in GSE197670.
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
import numpy as np
from scipy.stats import zscore, mannwhitneyu
from statsmodels.stats.multitest import multipletests


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def parse_gtf_promoters(gtf_path: Path, upstream=2000, downstream=500):
    """
    Extract promoter regions (-upstream to +downstream from TSS) for all genes.
    
    Args:
        gtf_path: Path to GTF file
        upstream: bp upstream of TSS
        downstream: bp downstream of TSS
        
    Returns:
        DataFrame with columns: gene_id, chrom, tss_pos, strand, 
                                prom_start, prom_end
    """
    promoters = []
    
    with gtf_path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
                
            chrom = fields[0]
            start = int(fields[3]) - 1  # 0-based
            end = int(fields[4])
            strand = fields[6]
            
            # Extract gene_id from attributes
            attrs = fields[8]
            gene_id = None
            for attr in attrs.split(";"):
                if "gene_id" in attr:
                    gene_id = attr.split('"')[1]
                    break
            
            if not gene_id:
                continue
            
            # Determine TSS and promoter boundaries
            if strand == "+":
                tss = start
                prom_start = max(0, tss - upstream)
                prom_end = tss + downstream
            else:
                tss = end
                prom_start = max(0, tss - downstream)
                prom_end = tss + upstream
            
            promoters.append({
                "gene_id": gene_id,
                "chrom": chrom,
                "tss_pos": tss,
                "strand": strand,
                "prom_start": prom_start,
                "prom_end": prom_end
            })
    
    return pd.DataFrame(promoters)


def overlap_regions(promoter_row, meth_df, within=True):
    """
    Find CpGs overlapping a promoter region.
    
    Args:
        promoter_row: Row with chrom, prom_start, prom_end
        meth_df: CpG methylation data (chrom, start, end, beta columns)
        within: If True, CpG must be completely within; else any overlap
        
    Returns:
        List of beta values for overlapping CpGs
    """
    if meth_df.empty:
        return []
    
    subset = meth_df[meth_df["chrom"] == promoter_row["chrom"]]
    if within:
        overlap = subset[
            (subset["start"] >= promoter_row["prom_start"]) &
            (subset["end"] <= promoter_row["prom_end"])
        ]
    else:
        overlap = subset[
            (subset["start"] < promoter_row["prom_end"]) &
            (subset["end"] > promoter_row["prom_start"])
        ]
    
    return overlap["beta"].values.tolist()


def aggregate_promoter_methylation(gene_meth_matrix, promoters, meth_cpg_df=None):
    """
    Aggregate promoter CpG methylation to gene level.
    
    If gene_meth_matrix exists, use directly.
    Otherwise, use CpG-level data and promoter regions.
    
    Args:
        gene_meth_matrix: DataFrame with genes x samples (may be pre-aggregated)
        promoters: DataFrame with promoter regions
        meth_cpg_df: Optional CpG-level data (chrom, start, end, beta)
        
    Returns:
        DataFrame: promoter_methyl (genes x samples)
    """
    if not gene_meth_matrix.empty:
        # Already aggregated
        promoters_subset = promoters[promoters["gene_id"].isin(gene_meth_matrix.index)]
        return gene_meth_matrix.loc[promoters_subset["gene_id"]]
    
    if meth_cpg_df is None or meth_cpg_df.empty:
        logger.warning("No methylation data provided")
        return pd.DataFrame()
    
    # Aggregate CpGs to promoters per sample
    gene_meth = {}
    for _, prom_row in promoters.iterrows():
        gene_id = prom_row["gene_id"]
        betas = overlap_regions(prom_row, meth_cpg_df)
        if betas:
            gene_meth[gene_id] = np.mean(betas)
    
    return pd.Series(gene_meth)


def identify_hypermethylated_genes(gene_meth_matrix, delta_beta_threshold=0.2, n_samples_min=3):
    """
    Identify promoter hypermethylated genes.
    
    Genes with mean promoter beta > (global median + delta_beta_threshold)
    
    Returns:
        DataFrame: genes, mean_beta, z_score, classification
    """
    if gene_meth_matrix.empty:
        return pd.DataFrame()
    
    mean_beta = gene_meth_matrix.mean(axis=1)
    global_median = gene_meth_matrix.values.flatten()
    global_median = np.nanmedian(global_median)
    
    z_scores = zscore(mean_beta, nan_policy='omit')
    
    # Hypermethylated: mean_beta > median + threshold
    hyper = mean_beta > (global_median + delta_beta_threshold)
    
    result = pd.DataFrame({
        "gene_id": gene_meth_matrix.index,
        "mean_beta": mean_beta.values,
        "z_score": z_scores,
        "n_samples_nonmissing": (~gene_meth_matrix.isna()).sum(axis=1).values,
        "classification": "hypermethylated"
    }).set_index("gene_id")

    sample_mask = result["n_samples_nonmissing"] >= n_samples_min
    result = result[hyper & sample_mask]
    return result.reset_index().sort_values("mean_beta", ascending=False)


def identify_hypomethylated_genes(gene_meth_matrix, delta_beta_threshold=0.2, n_samples_min=3):
    """
    Identify promoter hypomethylated genes.
    
    Genes with mean promoter beta < (global median - delta_beta_threshold)
    """
    if gene_meth_matrix.empty:
        return pd.DataFrame()
    
    mean_beta = gene_meth_matrix.mean(axis=1)
    global_median = gene_meth_matrix.values.flatten()
    global_median = np.nanmedian(global_median)
    
    z_scores = zscore(mean_beta, nan_policy='omit')
    
    # Hypomethylated: mean_beta < median - threshold
    hypo = mean_beta < (global_median - delta_beta_threshold)
    
    result = pd.DataFrame({
        "gene_id": gene_meth_matrix.index,
        "mean_beta": mean_beta.values,
        "z_score": z_scores,
        "n_samples_nonmissing": (~gene_meth_matrix.isna()).sum(axis=1).values,
        "classification": "hypomethylated"
    }).set_index("gene_id")

    sample_mask = result["n_samples_nonmissing"] >= n_samples_min
    result = result[hypo & sample_mask]
    return result.reset_index().sort_values("mean_beta", ascending=True)


def identify_differential_methylation(
    gene_meth_matrix: pd.DataFrame,
    phenotype_dict: dict[str, str],
    hf_phenotypes: set[str],
    nf_phenotype: str = "NF",
    delta_beta_threshold: float = 0.2,
    fdr_cutoff: float = 0.05,
    n_samples_min: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Identify hyper/hypo methylated genes using HF-vs-NF differential testing."""
    nf_samples = [s for s in gene_meth_matrix.columns if phenotype_dict.get(s) == nf_phenotype]
    hf_samples = [s for s in gene_meth_matrix.columns if phenotype_dict.get(s) in hf_phenotypes]

    if len(nf_samples) < n_samples_min or len(hf_samples) < n_samples_min:
        raise ValueError(
            f"Insufficient samples for differential methylation: NF={len(nf_samples)}, HF={len(hf_samples)}"
        )

    rows = []
    for gene_id in gene_meth_matrix.index:
        nf_vals = pd.to_numeric(gene_meth_matrix.reindex(index=[gene_id], columns=nf_samples).iloc[0], errors="coerce")
        hf_vals = pd.to_numeric(gene_meth_matrix.reindex(index=[gene_id], columns=hf_samples).iloc[0], errors="coerce")

        nf = nf_vals.dropna().values.astype(float)
        hf = hf_vals.dropna().values.astype(float)

        if len(nf) < n_samples_min or len(hf) < n_samples_min:
            continue

        stat, pvalue = mannwhitneyu(hf, nf, alternative="two-sided")
        delta_beta = float(np.mean(hf) - np.mean(nf))

        rows.append(
            {
                "gene_id": gene_id,
                "hf_mean_beta": float(np.mean(hf)),
                "nf_mean_beta": float(np.mean(nf)),
                "delta_beta": delta_beta,
                "pvalue": float(pvalue),
                "n_hf": int(len(hf)),
                "n_nf": int(len(nf)),
            }
        )

    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(), pd.DataFrame()

    out["fdr"] = multipletests(out["pvalue"].values, method="fdr_bh")[1]
    sig = out[(out["fdr"] <= fdr_cutoff) & (out["delta_beta"].abs() >= delta_beta_threshold)].copy()

    hyper = sig[sig["delta_beta"] > 0].copy()
    hyper["classification"] = "hypermethylated"
    hyper = hyper.sort_values("delta_beta", ascending=False)

    hypo = sig[sig["delta_beta"] < 0].copy()
    hypo["classification"] = "hypomethylated"
    hypo = hypo.sort_values("delta_beta", ascending=True)

    return hyper, hypo


def assign_pathway_context(genes, pathway_genesets_dict):
    """
    Annotate genes with pathway membership.
    
    Args:
        genes: List or Series of gene IDs
        pathway_genesets_dict: Dict[pathway_name] = set(gene_ids)
        
    Returns:
        DataFrame with gene_id, pathways (comma-separated)
    """
    result = []
    for gene in genes:
        pathways = [pw for pw, genes_set in pathway_genesets_dict.items() 
                   if gene in genes_set]
        result.append({
            "gene_id": gene,
            "pathways": ",".join(pathways) if pathways else "unclassified"
        })
    return pd.DataFrame(result)


def main():
    ap = argparse.ArgumentParser(
        description="Extract promoter methylation signatures from GSE123976"
    )
    ap.add_argument("--gene-methylation", required=True, type=Path,
                   help="Gene-level methylation matrix (genes x samples)")
    ap.add_argument("--metadata", required=True, type=Path,
                   help="Sample metadata TSV with phenotype column")
    ap.add_argument("--gtf", required=True, type=Path,
                   help="GTF annotation file")
    ap.add_argument("--pathway-genesets", type=Path,
                   help="TSV: pathway_name<TAB>gene1,gene2,... (optional)")
    ap.add_argument("--delta-beta-threshold", type=float, default=0.2,
                   help="Methylation change threshold for classification")
    ap.add_argument("--fdr-cutoff", type=float, default=0.05,
                   help="FDR cutoff for differential methylation (default: 0.05)")
    ap.add_argument("--min-samples", type=int, default=3,
                   help="Minimum non-missing samples per gene")
    ap.add_argument("--out-hypermethylated", required=True, type=Path,
                   help="Output: promoter hypermethylated genes")
    ap.add_argument("--out-hypomethylated", required=True, type=Path,
                   help="Output: promoter hypomethylated genes")
    ap.add_argument("--out-summary", required=True, type=Path,
                   help="Output: summary statistics")
    
    args = ap.parse_args()
    
    # Load methylation matrix
    logger.info(f"Loading gene methylation matrix from {args.gene_methylation}")
    gene_meth = pd.read_csv(args.gene_methylation, sep="\t", index_col=0)
    logger.info(f"Loaded {gene_meth.shape[0]} genes x {gene_meth.shape[1]} samples")

    meta = pd.read_csv(args.metadata, sep="\t")
    sample_col = meta.columns[0]
    pheno_col = "phenotype" if "phenotype" in meta.columns else meta.columns[1]
    phenotype_dict = {
        str(sample): str(pheno)
        for sample, pheno in zip(meta[sample_col].astype(str), meta[pheno_col].astype(str))
    }
    
    # Load pathway genesets if provided
    pathway_genesets = {}
    if args.pathway_genesets and args.pathway_genesets.exists():
        logger.info(f"Loading pathway genesets from {args.pathway_genesets}")
        with args.pathway_genesets.open() as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    pathway = parts[0]
                    genes = set(parts[1].split(","))
                    pathway_genesets[pathway] = genes
        logger.info(f"Loaded {len(pathway_genesets)} pathways")
    
    # Identify hyper/hypomethylated genes via differential testing
    logger.info("Identifying differential methylation signatures (HF vs NF)...")
    hyper, hypo = identify_differential_methylation(
        gene_meth,
        phenotype_dict=phenotype_dict,
        hf_phenotypes={"HF", "DCM", "ICM"},
        nf_phenotype="NF",
        delta_beta_threshold=args.delta_beta_threshold,
        fdr_cutoff=args.fdr_cutoff,
        n_samples_min=args.min_samples,
    )

    logger.info(f"Found {len(hyper)} hypermethylated genes")
    
    if pathway_genesets and not hyper.empty:
        hyper_pathways = assign_pathway_context(hyper["gene_id"], pathway_genesets)
        hyper = hyper.merge(hyper_pathways, on="gene_id", how="left")
    
    hyper.to_csv(args.out_hypermethylated, sep="\t", index=False)
    logger.info(f"Saved to {args.out_hypermethylated}")
    
    logger.info(f"Found {len(hypo)} hypomethylated genes")
    
    if pathway_genesets and not hypo.empty:
        hypo_pathways = assign_pathway_context(hypo["gene_id"], pathway_genesets)
        hypo = hypo.merge(hypo_pathways, on="gene_id", how="left")
    
    hypo.to_csv(args.out_hypomethylated, sep="\t", index=False)
    logger.info(f"Saved to {args.out_hypomethylated}")
    
    # Summary statistics
    summary = {
        "discovery_samples": gene_meth.shape[1],
        "discovery_genes": gene_meth.shape[0],
        "hypermethylated_count": len(hyper),
        "hypomethylated_count": len(hypo),
        "delta_beta_threshold": args.delta_beta_threshold,
        "fdr_cutoff": args.fdr_cutoff,
        "min_samples": args.min_samples,
        "pathways_included": len(pathway_genesets)
    }
    
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(args.out_summary, sep="\t", index=False)
    logger.info(f"Summary saved to {args.out_summary}")
    logger.info(f"Hypermethylated: {summary['hypermethylated_count']}, " \
                f"Hypomethylated: {summary['hypomethylated_count']}")


if __name__ == "__main__":
    main()
