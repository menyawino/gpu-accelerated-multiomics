#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from processed_analysis_common import (
    differential_expression,
    enrichment_by_hypergeom,
    load_gmt,
    read_matrix,
)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rna-matrix", required=True)
    ap.add_argument("--rna-metadata", required=True)
    ap.add_argument("--methylation-table", required=True)
    ap.add_argument("--methylation-gene-column", required=True)
    ap.add_argument("--methylation-effect-column", required=True)
    ap.add_argument("--methylation-pvalue-column", required=True)
    ap.add_argument("--hf-label", required=True)
    ap.add_argument("--control-label", required=True)
    ap.add_argument("--pathways-gmt", required=True)
    ap.add_argument("--fdr-threshold", type=float, required=True)
    ap.add_argument("--log2fc-threshold", type=float, required=True)
    ap.add_argument("--min-group-size", type=int, required=True)
    ap.add_argument("--out-deg", required=True)
    ap.add_argument("--out-integrated", required=True)
    ap.add_argument("--out-enrichment", required=True)
    ap.add_argument("--out-qc-plot", required=True)
    ap.add_argument("--out-qc-summary", required=True)
    args = ap.parse_args()

    expr = read_matrix(args.rna_matrix, index_col="gene_id")
    metadata = pd.read_csv(args.rna_metadata, sep="\t")

    deg = differential_expression(
        expr=expr,
        metadata=metadata,
        case_label=args.hf_label,
        control_label=args.control_label,
        min_group_size=args.min_group_size,
    )
    deg["de_direction"] = np.where(deg["log2fc"] > 0, "Up", "Down")
    deg["de_significant"] = (
        (deg["fdr"] <= args.fdr_threshold) & (deg["log2fc"].abs() >= args.log2fc_threshold)
    )

    meth = pd.read_csv(args.methylation_table, sep="\t")
    meth = meth.rename(
        columns={
            args.methylation_gene_column: "gene_id",
            args.methylation_effect_column: "delta_methylation",
            args.methylation_pvalue_column: "methylation_pvalue",
        }
    )[["gene_id", "delta_methylation", "methylation_pvalue"]]
    meth["methylation_direction"] = np.where(meth["delta_methylation"] >= 0, "Hyper", "Hypo")
    meth["methylation_significant"] = meth["methylation_pvalue"] <= args.fdr_threshold

    integrated = deg.merge(meth, on="gene_id", how="inner")
    integrated["integrated_class"] = "Other"
    integrated.loc[
        integrated["methylation_direction"].eq("Hyper")
        & integrated["methylation_significant"]
        & integrated["de_significant"]
        & (integrated["log2fc"] < 0),
        "integrated_class",
    ] = "Hyper-down"
    integrated.loc[
        integrated["methylation_direction"].eq("Hypo")
        & integrated["methylation_significant"]
        & integrated["de_significant"]
        & (integrated["log2fc"] > 0),
        "integrated_class",
    ] = "Hypo-up"
    integrated.loc[
        integrated["methylation_direction"].eq("Hyper")
        & integrated["methylation_significant"]
        & integrated["de_significant"]
        & (integrated["log2fc"] > 0),
        "integrated_class",
    ] = "Hyper-up"
    integrated.loc[
        integrated["methylation_direction"].eq("Hypo")
        & integrated["methylation_significant"]
        & integrated["de_significant"]
        & (integrated["log2fc"] < 0),
        "integrated_class",
    ] = "Hypo-down"

    pathways = load_gmt(args.pathways_gmt)
    universe = set(integrated["gene_id"].astype(str))
    enrich_rows = []
    for class_name in ["Hyper-down", "Hypo-up", "Hyper-up", "Hypo-down"]:
        genes = set(integrated.loc[integrated["integrated_class"] == class_name, "gene_id"].astype(str))
        tbl = enrichment_by_hypergeom(genes, pathways, universe)
        if tbl.empty:
            continue
        tbl.insert(0, "gene_class", class_name)
        enrich_rows.append(tbl)
    enrichment = pd.concat(enrich_rows, ignore_index=True) if enrich_rows else pd.DataFrame(
        columns=["gene_class", "pathway", "overlap", "pathway_size", "query_size", "pvalue", "fdr"]
    )

    plt.figure(figsize=(7, 5))
    y = -np.log10(deg["fdr"].clip(lower=1e-300))
    plt.scatter(deg["log2fc"], y, s=8, alpha=0.4, color="#999999")
    sig = deg["de_significant"]
    plt.scatter(deg.loc[sig, "log2fc"], y[sig], s=10, alpha=0.8, color="#b2182b")
    plt.axvline(args.log2fc_threshold, ls="--", lw=0.8, color="black")
    plt.axvline(-args.log2fc_threshold, ls="--", lw=0.8, color="black")
    plt.axhline(-np.log10(args.fdr_threshold), ls="--", lw=0.8, color="black")
    plt.xlabel("log2 fold-change (HF vs NF)")
    plt.ylabel("-log10(FDR)")
    plt.title("GSE123976 differential expression")
    plt.tight_layout()
    Path(args.out_qc_plot).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(args.out_qc_plot, dpi=150)
    plt.close()

    qc_summary = pd.DataFrame(
        {
            "metric": [
                "n_genes_tested",
                "n_de_significant",
                "n_hyper_down",
                "n_hypo_up",
            ],
            "value": [
                int(deg.shape[0]),
                int(deg["de_significant"].sum()),
                int((integrated["integrated_class"] == "Hyper-down").sum()),
                int((integrated["integrated_class"] == "Hypo-up").sum()),
            ],
        }
    )

    deg.sort_values(["fdr", "gene_id"]).to_csv(args.out_deg, sep="\t", index=False)
    integrated.sort_values(["integrated_class", "fdr", "gene_id"]).to_csv(args.out_integrated, sep="\t", index=False)
    enrichment.to_csv(args.out_enrichment, sep="\t", index=False)
    qc_summary.to_csv(args.out_qc_summary, sep="\t", index=False)


if __name__ == "__main__":
    main()
