#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import pandas as pd


def _normalize_subject(text: str) -> str:
    if pd.isna(text):
        return ""
    s = str(text).upper()
    # Heuristic extraction of subject identifiers from free text.
    patterns = [
        r"(PATIENT[_\- ]?\d+)",
        r"(SUBJECT[_\- ]?\d+)",
        r"(SAMPLE[_\- ]?\d+)",
        r"(HF[_\- ]?\d+)",
        r"(CTRL[_\- ]?\d+)",
    ]
    for pat in patterns:
        m = re.search(pat, s)
        if m:
            return re.sub(r"[^A-Z0-9]", "", m.group(1))
    # Fallback: use first alphanumeric token of at least length 4.
    for token in re.findall(r"[A-Z0-9]{4,}", s):
        if token not in {"RNA", "WGBS", "SEQ", "SAMPLE", "TREAT", "CONTROL"}:
            return token
    return ""


def discover(gse: str, metadata_dir: Path, manual_pairs_path: Path):
    metadata_dir.mkdir(parents=True, exist_ok=True)

    gse_map = pd.read_csv(
        Path("/dev/stdin"),
        sep="\t",
        engine="python",
    ) if False else None

    # Use pysradb directly via shell-safe import to avoid subprocess complexity.
    import subprocess

    gse_to_srp_path = metadata_dir / "gse_to_srp.tsv"
    with gse_to_srp_path.open("w") as f:
        subprocess.run(["pysradb", "gse-to-srp", gse], check=True, stdout=f)

    gse_map = pd.read_csv(gse_to_srp_path, sep="\t")
    if gse_map.empty:
        raise RuntimeError(f"No SRP found for {gse}")

    srp = gse_map.iloc[0].dropna().astype(str).tolist()[-1]

    detailed_path = metadata_dir / "srp_detailed.tsv"
    with detailed_path.open("w") as f:
        subprocess.run(["pysradb", "metadata", srp, "--detailed", "--expand"], check=True, stdout=f)

    meta = pd.read_csv(detailed_path, sep="\t")
    if "run_accession" not in meta.columns:
        raise RuntimeError("pysradb metadata output missing run_accession")

    strategy_col = "library_strategy" if "library_strategy" in meta.columns else None
    if strategy_col is None:
        raise RuntimeError("Cannot determine library_strategy for RNA-seq vs WGBS split")

    rna = meta[meta[strategy_col].str.contains("RNA", case=False, na=False)].copy()
    wgbs = meta[
        meta[strategy_col].str.contains("BISULFITE|WGBS", case=False, na=False)
    ].copy()

    if wgbs.empty:
        # Fallback for datasets that annotate WGBS as WGS with bisulfite tags in source names.
        text_cols = [c for c in ["sample_attribute", "experiment_title", "study_title"] if c in meta.columns]
        if text_cols:
            joined = meta[text_cols].fillna("").agg(" ".join, axis=1)
            wgbs = meta[joined.str.contains("BISULFITE|WGBS|METHYL", case=False, na=False)].copy()

    for df in (rna, wgbs):
        for c in ["sample_title", "experiment_title", "sample_alias", "run_alias"]:
            if c not in df.columns:
                df[c] = ""
        df["subject_id"] = (
            df["sample_title"].fillna("")
            + " "
            + df["experiment_title"].fillna("")
            + " "
            + df["sample_alias"].fillna("")
            + " "
            + df["run_alias"].fillna("")
        ).map(_normalize_subject)

    keep_cols = [
        c
        for c in [
            "run_accession",
            "experiment_accession",
            "sample_accession",
            "study_accession",
            "library_strategy",
            "library_layout",
            "sample_title",
            "experiment_title",
            "subject_id",
        ]
        if c in meta.columns or c == "subject_id"
    ]

    rna = rna[keep_cols].drop_duplicates().rename(columns={"run_accession": "rnaseq_run"})
    wgbs = wgbs[keep_cols].drop_duplicates().rename(columns={"run_accession": "wgbs_run"})

    rnaseq_runs_tsv = metadata_dir / "rnaseq_runs.tsv"
    wgbs_runs_tsv = metadata_dir / "wgbs_runs.tsv"
    rna.to_csv(rnaseq_runs_tsv, sep="\t", index=False)
    wgbs.to_csv(wgbs_runs_tsv, sep="\t", index=False)

    (metadata_dir / "rnaseq_runs.txt").write_text(
        "\n".join(sorted(rna["rnaseq_run"].dropna().astype(str).unique())) + "\n"
    )
    (metadata_dir / "wgbs_runs.txt").write_text(
        "\n".join(sorted(wgbs["wgbs_run"].dropna().astype(str).unique())) + "\n"
    )

    pairs = []
    if manual_pairs_path.exists():
        manual = pd.read_csv(manual_pairs_path, sep="\t")
        required = {"subject_id", "rnaseq_run", "wgbs_run"}
        if not required.issubset(set(manual.columns)):
            raise RuntimeError(
                f"Manual pairing file {manual_pairs_path} must have columns: {sorted(required)}"
            )
        pairs = manual[["subject_id", "rnaseq_run", "wgbs_run"]].drop_duplicates()
    else:
        if not rna.empty and not wgbs.empty:
            pairs = (
                rna[["subject_id", "rnaseq_run"]]
                .merge(wgbs[["subject_id", "wgbs_run"]], on="subject_id", how="inner")
                .dropna()
                .drop_duplicates()
            )
        else:
            pairs = pd.DataFrame(columns=["subject_id", "rnaseq_run", "wgbs_run"])

    pairs_path = metadata_dir / "matched_pairs.tsv"
    if isinstance(pairs, list):
        pd.DataFrame(columns=["subject_id", "rnaseq_run", "wgbs_run"]).to_csv(
            pairs_path, sep="\t", index=False
        )
    else:
        pairs.to_csv(pairs_path, sep="\t", index=False)

    if len(pairs) == 0:
        print(
            "WARNING: No matched RNA/WGBS pairs inferred. Provide config/manual_pairs.tsv for deterministic pairing."
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gse", required=True)
    parser.add_argument("--metadata-dir", required=True, type=Path)
    parser.add_argument("--manual-pairs", required=True, type=Path)
    args = parser.parse_args()

    discover(args.gse, args.metadata_dir, args.manual_pairs)


if __name__ == "__main__":
    main()
