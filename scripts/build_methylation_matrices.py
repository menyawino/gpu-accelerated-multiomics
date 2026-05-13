#!/usr/bin/env python3
import argparse
import gzip
import re
from pathlib import Path

import pandas as pd
from pybedtools import BedTool

# Promoter window relative to TSS (strand-aware):
#   upstream (bp before TSS) and downstream (bp after TSS).
PROMOTER_UPSTREAM_BP = 2000
PROMOTER_DOWNSTREAM_BP = 500


def parse_gtf_promoters(gtf_path: Path):
    """Extract strand-aware promoter windows for each gene from a GTF file.

    Promoter is defined as [TSS - PROMOTER_UPSTREAM_BP, TSS + PROMOTER_DOWNSTREAM_BP].
    For '+' strand genes TSS = start; for '-' strand genes TSS = end.
    One promoter region is emitted per unique (gene_id, strand) pair using the most
    extreme TSS across all transcripts of that gene.
    """
    pat_gid = re.compile(r'gene_id "([^"]+)"')

    # Collect the min-start and max-end per gene to find canonical TSS.
    gene_info: dict[str, tuple[str, int, int, str]] = {}  # gene_id -> (chrom, min_start, max_end, strand)

    with gtf_path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = fields
            if feature != "gene":
                continue
            m_gid = pat_gid.search(attrs)
            if not m_gid:
                continue
            gid = m_gid.group(1)
            s, e = int(start) - 1, int(end)  # GTF is 1-based; convert to 0-based half-open
            if gid not in gene_info:
                gene_info[gid] = (chrom, s, e, strand)
            else:
                prev = gene_info[gid]
                gene_info[gid] = (chrom, min(prev[1], s), max(prev[2], e), strand)

    rows = []
    for gid, (chrom, gstart, gend, strand) in gene_info.items():
        if strand == "+":
            tss = gstart
        else:
            tss = gend  # 0-based exclusive end = last base of gene body = TSS for '-' strand
        prom_start = max(0, tss - PROMOTER_UPSTREAM_BP)
        prom_end = tss + PROMOTER_DOWNSTREAM_BP
        rows.append((chrom, prom_start, prom_end, gid, strand))

    promoters = pd.DataFrame(rows, columns=["chrom", "start", "end", "gene_id", "strand"])
    return promoters


def methyl_cov_to_bed(cov_path: Path) -> pd.DataFrame:
    opener = gzip.open if cov_path.suffix == ".gz" else open
    rows = []
    with opener(cov_path, "rt") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            chrom, start, end, _, meth, unmeth = fields[:6]
            meth = float(meth)
            unmeth = float(unmeth)
            total = meth + unmeth
            if total <= 0:
                continue
            beta = meth / total
            rows.append((chrom, int(start) - 1, int(end), beta))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "beta"])


def summarize_beta(cpg_df: pd.DataFrame, region_df: pd.DataFrame, id_col: str) -> pd.Series:
    if cpg_df.empty or region_df.empty:
        return pd.Series(dtype=float)

    cpg_bt = BedTool.from_dataframe(cpg_df[["chrom", "start", "end", "beta"]])
    reg_bt = BedTool.from_dataframe(region_df[["chrom", "start", "end", id_col]])

    inter = reg_bt.intersect(cpg_bt, wa=True, wb=True)
    vals = {}
    for iv in inter:
        rid = iv[3]
        beta = float(iv[7])
        vals.setdefault(rid, []).append(beta)

    out = {k: sum(v) / len(v) for k, v in vals.items() if v}
    return pd.Series(out)


def load_pairs(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["subject_id", "rnaseq_run", "wgbs_run"])
    df = pd.read_csv(path, sep="\t")
    if "wgbs_run" not in df.columns:
        return pd.DataFrame(columns=["subject_id", "rnaseq_run", "wgbs_run"])
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wgbs-runs", required=True, type=Path)
    ap.add_argument("--matched-pairs", required=True, type=Path)
    ap.add_argument("--gtf", required=True, type=Path)
    ap.add_argument("--methyl-dir", required=True, type=Path)
    ap.add_argument("--out-gene", required=True, type=Path)
    ap.add_argument(
        "--out-isoform",
        required=False,
        type=Path,
        default=None,
        help=(
            "Deprecated. Isoform-level methylation is not computed (promoter windows "
            "are gene-level). Providing this argument writes an empty placeholder."
        ),
    )
    args = ap.parse_args()

    wgbs = pd.read_csv(args.wgbs_runs, sep="\t") if args.wgbs_runs.exists() else pd.DataFrame()
    if wgbs.empty or "wgbs_run" not in wgbs.columns:
        pd.DataFrame().to_csv(args.out_gene, sep="\t")
        if args.out_isoform:
            pd.DataFrame().to_csv(args.out_isoform, sep="\t")
        return

    pairs = load_pairs(args.matched_pairs)
    keep_runs = set(pairs["wgbs_run"].astype(str)) if not pairs.empty else set(wgbs["wgbs_run"].astype(str))

    # Use promoter-restricted regions only (TSS-2000 to TSS+500, strand-aware).
    promoters = parse_gtf_promoters(args.gtf)

    gene_frames = []

    for run in sorted(keep_runs):
        cov = args.methyl_dir / run / f"{run}.bismark.cov.gz"
        if not cov.exists():
            continue
        cpg = methyl_cov_to_bed(cov)

        gser = summarize_beta(cpg, promoters, "gene_id").rename(run)
        gene_frames.append(gser)

    gene_mat = pd.concat(gene_frames, axis=1).fillna(0.0) if gene_frames else pd.DataFrame()
    gene_mat.to_csv(args.out_gene, sep="\t")

    # Write empty isoform placeholder so downstream rules that depend on the file don't fail.
    if args.out_isoform:
        pd.DataFrame(
            columns=gene_mat.columns,
            dtype=float,
        ).to_csv(args.out_isoform, sep="\t")


if __name__ == "__main__":
    main()
