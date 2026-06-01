#!/usr/bin/env python3
from __future__ import annotations

import argparse

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--discovery", required=True)
    ap.add_argument("--methylation-validation", required=True)
    ap.add_argument("--dcm-validation", required=True)
    ap.add_argument("--icm-validation", required=True)
    ap.add_argument("--out-summary", required=True)
    args = ap.parse_args()

    discovery = pd.read_csv(args.discovery, sep="\t")
    meth = pd.read_csv(args.methylation_validation, sep="\t")
    dcm = pd.read_csv(args.dcm_validation, sep="\t")[["gene_id", "log2fc", "fdr"]].rename(
        columns={"log2fc": "dcm_log2fc", "fdr": "dcm_fdr"}
    )
    icm = pd.read_csv(args.icm_validation, sep="\t")[["gene_id", "log2fc", "fdr"]].rename(
        columns={"log2fc": "icm_log2fc", "fdr": "icm_fdr"}
    )

    out = (
        discovery.merge(meth, on="gene_id", how="left")
        .merge(dcm, on="gene_id", how="left")
        .merge(icm, on="gene_id", how="left")
        .sort_values(["integrated_class", "gene_id"])
    )
    out.to_csv(args.out_summary, sep="\t", index=False)


if __name__ == "__main__":
    main()
