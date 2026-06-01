#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def markdown_table(df: pd.DataFrame, n: int = 10) -> str:
    if df.empty:
        return "No rows available."
    return df.head(n).to_markdown(index=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--discovery", required=True)
    ap.add_argument("--methylation-stats", required=True)
    ap.add_argument("--rna-stats", required=True)
    ap.add_argument("--integrated-summary", required=True)
    ap.add_argument("--out-report", required=True)
    args = ap.parse_args()

    discovery = pd.read_csv(args.discovery, sep="\t")
    meth_stats = pd.read_csv(args.methylation_stats, sep="\t")
    rna_stats = pd.read_csv(args.rna_stats, sep="\t")
    integrated = pd.read_csv(args.integrated_summary, sep="\t")

    class_counts = discovery["integrated_class"].value_counts().rename_axis("class").reset_index(name="n_genes")

    report = f"""---
title: "Reproduce and validate GSE123976 methylation–expression program"
format: gfm
---

## Discovery reproduction (GSE123976)

Integrated class counts:

{markdown_table(class_counts, n=10)}

Top integrated genes:

{markdown_table(integrated[[c for c in ["gene_id", "integrated_class", "log2fc", "delta_methylation", "dcm_log2fc", "icm_log2fc"] if c in integrated.columns]], n=20)}

## External methylation validation (GSE197670/GSE197672 processed arm)

{markdown_table(meth_stats, n=10)}

## External transcriptome validation (GSE116250)

{markdown_table(rna_stats, n=10)}

## Figures

- Discovery QC volcano: `results/gse123976/qc_volcano.png`
- Methylation concordance: `results/gse197670/concordance_plot.png`
- Transcriptome signature scores: `results/gse116250/signature_score_plot.png`

## Notes

- Validation is performed at gene/program/pathway level and is not designed for exact CpG-for-CpG replication between WGBS and array platforms.
- GSE197671 is intentionally excluded from myocardial validation.
"""

    out = Path(args.out_report)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(report)


if __name__ == "__main__":
    main()
