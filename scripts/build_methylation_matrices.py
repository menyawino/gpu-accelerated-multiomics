#!/usr/bin/env python3
import argparse
import gzip
import re
from pathlib import Path

import pandas as pd
from pybedtools import BedTool


def parse_gtf_regions(gtf_path: Path):
    gene_rows = []
    tx_rows = []

    pat_gid = re.compile(r'gene_id "([^"]+)"')
    pat_tid = re.compile(r'transcript_id "([^"]+)"')

    with gtf_path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = fields
            m_gid = pat_gid.search(attrs)
            m_tid = pat_tid.search(attrs)
            if feature == "gene" and m_gid:
                gid = m_gid.group(1)
                gene_rows.append((chrom, int(start) - 1, int(end), gid, strand))
            elif feature == "transcript" and m_gid and m_tid:
                gid = m_gid.group(1)
                tid = m_tid.group(1)
                tx_rows.append((chrom, int(start) - 1, int(end), tid, gid, strand))

    genes = pd.DataFrame(gene_rows, columns=["chrom", "start", "end", "gene_id", "strand"])
    tx = pd.DataFrame(tx_rows, columns=["chrom", "start", "end", "isoform_id", "gene_id", "strand"])

    return genes, tx


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
    ap.add_argument("--out-isoform", required=True, type=Path)
    args = ap.parse_args()

    wgbs = pd.read_csv(args.wgbs_runs, sep="\t") if args.wgbs_runs.exists() else pd.DataFrame()
    if wgbs.empty or "wgbs_run" not in wgbs.columns:
        pd.DataFrame().to_csv(args.out_gene, sep="\t")
        pd.DataFrame().to_csv(args.out_isoform, sep="\t")
        return

    pairs = load_pairs(args.matched_pairs)
    keep_runs = set(pairs["wgbs_run"].astype(str)) if not pairs.empty else set(wgbs["wgbs_run"].astype(str))

    genes, isoforms = parse_gtf_regions(args.gtf)

    gene_frames = []
    iso_frames = []

    for run in sorted(keep_runs):
        cov = args.methyl_dir / run / f"{run}.bismark.cov.gz"
        if not cov.exists():
            continue
        cpg = methyl_cov_to_bed(cov)

        gser = summarize_beta(cpg, genes, "gene_id").rename(run)
        iser = summarize_beta(cpg, isoforms, "isoform_id").rename(run)

        gene_frames.append(gser)
        iso_frames.append(iser)

    gene_mat = pd.concat(gene_frames, axis=1).fillna(0.0) if gene_frames else pd.DataFrame()
    iso_mat = pd.concat(iso_frames, axis=1).fillna(0.0) if iso_frames else pd.DataFrame()

    gene_mat.to_csv(args.out_gene, sep="\t")
    iso_mat.to_csv(args.out_isoform, sep="\t")


if __name__ == "__main__":
    main()
