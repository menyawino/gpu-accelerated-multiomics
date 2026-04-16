#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import pandas as pd


def parse_gtf_transcript_to_gene(gtf_path: Path) -> dict:
    mapping = {}
    pat_tid = re.compile(r'transcript_id "([^"]+)"')
    pat_gid = re.compile(r'gene_id "([^"]+)"')

    with gtf_path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature not in {"transcript", "exon"}:
                continue
            attrs = fields[8]
            m_tid = pat_tid.search(attrs)
            m_gid = pat_gid.search(attrs)
            if m_tid and m_gid:
                mapping[m_tid.group(1)] = m_gid.group(1)
    return mapping


def load_pairs(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["subject_id", "rnaseq_run", "wgbs_run"])
    df = pd.read_csv(path, sep="\t")
    if "rnaseq_run" not in df.columns:
        return pd.DataFrame(columns=["subject_id", "rnaseq_run", "wgbs_run"])
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rnaseq-runs", required=True, type=Path)
    ap.add_argument("--matched-pairs", required=True, type=Path)
    ap.add_argument("--gtf", required=True, type=Path)
    ap.add_argument("--salmon-dir", required=True, type=Path)
    ap.add_argument("--out-gene", required=True, type=Path)
    ap.add_argument("--out-isoform", required=True, type=Path)
    args = ap.parse_args()

    rnaseq = pd.read_csv(args.rnaseq_runs, sep="\t") if args.rnaseq_runs.exists() else pd.DataFrame()
    if rnaseq.empty or "rnaseq_run" not in rnaseq.columns:
        pd.DataFrame().to_csv(args.out_gene, sep="\t")
        pd.DataFrame().to_csv(args.out_isoform, sep="\t")
        return

    pairs = load_pairs(args.matched_pairs)
    keep_runs = set(pairs["rnaseq_run"].astype(str)) if not pairs.empty else set(rnaseq["rnaseq_run"].astype(str))

    t2g = parse_gtf_transcript_to_gene(args.gtf)

    iso_frames = []
    gene_frames = []

    for run in sorted(keep_runs):
        qsf = args.salmon_dir / run / "quant.sf"
        if not qsf.exists():
            continue
        q = pd.read_csv(qsf, sep="\t")
        if not {"Name", "TPM"}.issubset(q.columns):
            continue

        iso = q[["Name", "TPM"]].rename(columns={"Name": "isoform_id", "TPM": run})
        iso_frames.append(iso.set_index("isoform_id"))

        q["gene_id"] = q["Name"].map(t2g)
        g = q.dropna(subset=["gene_id"]).groupby("gene_id", as_index=False)["TPM"].sum()
        gene_frames.append(g.rename(columns={"TPM": run}).set_index("gene_id"))

    if iso_frames:
        iso_mat = pd.concat(iso_frames, axis=1).fillna(0.0)
    else:
        iso_mat = pd.DataFrame()

    if gene_frames:
        gene_mat = pd.concat(gene_frames, axis=1).fillna(0.0)
    else:
        gene_mat = pd.DataFrame()

    iso_mat.to_csv(args.out_isoform, sep="\t")
    gene_mat.to_csv(args.out_gene, sep="\t")


if __name__ == "__main__":
    main()
