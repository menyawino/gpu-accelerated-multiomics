#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml


def must_exist(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing required input file: {path}")


def require_columns(path: Path, required: set[str]) -> None:
    df = pd.read_csv(path, sep="\t", nrows=5)
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--touch", required=True)
    args = ap.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text())
    inputs = cfg["inputs"]

    g123 = inputs["gse123976"]
    must_exist(Path(g123["rna_matrix"]))
    must_exist(Path(g123["rna_metadata"]))
    must_exist(Path(g123["methylation_table"]))
    require_columns(Path(g123["rna_metadata"]), {"sample_id", "group"})
    require_columns(
        Path(g123["methylation_table"]),
        {g123["methylation_gene_column"], g123["methylation_effect_column"], g123["methylation_pvalue_column"]},
    )

    g197670 = inputs["gse197670"]
    must_exist(Path(g197670["beta_matrix"]))
    must_exist(Path(g197670["metadata"]))
    must_exist(Path(g197670["probe_annotation"]))
    require_columns(Path(g197670["metadata"]), {"sample_id", "group"})
    require_columns(Path(g197670["probe_annotation"]), {"probe_id", "gene_id", "is_promoter"})

    g116250 = inputs["gse116250"]
    must_exist(Path(g116250["rna_matrix"]))
    must_exist(Path(g116250["metadata"]))
    require_columns(Path(g116250["metadata"]), {"sample_id", "group"})

    touch = Path(args.touch)
    touch.parent.mkdir(parents=True, exist_ok=True)
    touch.write_text("ok\n")


if __name__ == "__main__":
    main()
