#!/usr/bin/env python3
"""Run preranked GSEA on external curated pathway gene sets using gseapy.

Gene sets MUST be provided as external, pre-defined files (GMT format or one-per-line
plain text files organised into categories). Gene sets derived from the same effect table
being ranked would be circular and are explicitly forbidden here.
"""

import argparse
from pathlib import Path
from typing import Dict, List

import gseapy as gp
import pandas as pd

# Human-readable display names for known curated gene-set categories.
CATEGORY_LABELS: Dict[str, str] = {
    "fibrosis": "Fibrosis",
    "inflammation": "Inflammation",
    "ecm_remodeling": "ECM Remodeling",
}


def load_gmt(gmt_path: Path) -> Dict[str, List[str]]:
    """Parse a standard GMT file (name TAB description TAB gene ...) into a dict."""
    gene_sets: Dict[str, List[str]] = {}
    with gmt_path.open() as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0].strip()
            genes = [g.strip() for g in parts[2:] if g.strip()]
            if genes:
                gene_sets[name] = genes
    return gene_sets


def load_category_files(
    fibrosis_path: Path | None,
    inflammation_path: Path | None,
    ecm_path: Path | None,
) -> Dict[str, List[str]]:
    """Build gene-set dict from per-category plain-text files (one gene per line)."""
    gene_sets: Dict[str, List[str]] = {}
    for label, path in [
        ("Fibrosis", fibrosis_path),
        ("Inflammation", inflammation_path),
        ("ECM Remodeling", ecm_path),
    ]:
        if path is None or not path.exists():
            continue
        genes = [line.strip() for line in path.read_text().splitlines() if line.strip()]
        if len(genes) >= 5:
            gene_sets[label] = genes
    return gene_sets


def write_gmt(gene_sets: Dict[str, List[str]], out_path: Path) -> None:
    with out_path.open("w") as f:
        for name, genes in gene_sets.items():
            f.write("\t".join([name, "curated_external"] + genes) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Run preranked GSEA using external, curated gene sets. "
            "Providing gene sets derived from the same effect table as the ranking "
            "is circular and not supported."
        )
    )
    ap.add_argument(
        "--ranking-table",
        type=Path,
        default=Path("results/validation/reprogramming/all_genes_failing_vs_nonfailing.tsv"),
        help="TSV with columns gene_id and mean_diff (or combined_effect) used as the ranking metric.",
    )
    ap.add_argument(
        "--gene-sets-gmt",
        type=Path,
        default=None,
        help="External GMT file with curated gene sets (e.g., from MSigDB or manual curation). "
             "Mutually exclusive with --fibrosis-genes / --inflammation-genes / --ecm-genes.",
    )
    ap.add_argument(
        "--fibrosis-genes",
        type=Path,
        default=Path("resources/genesets/fibrosis.txt"),
        help="One gene per line: curated fibrosis gene set.",
    )
    ap.add_argument(
        "--inflammation-genes",
        type=Path,
        default=Path("resources/genesets/inflammation.txt"),
        help="One gene per line: curated inflammation gene set.",
    )
    ap.add_argument(
        "--ecm-genes",
        type=Path,
        default=Path("resources/genesets/ecm_remodeling.txt"),
        help="One gene per line: curated ECM-remodeling gene set.",
    )
    ap.add_argument("--rank-col", default="mean_diff",
                    help="Column to use as ranking metric (default: mean_diff).")
    ap.add_argument("--out-dir", type=Path, default=Path("results/validation/reprogramming/gsea"))
    args = ap.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load ranking table.
    ranking = pd.read_csv(args.ranking_table, sep="\t")
    rank_col = args.rank_col if args.rank_col in ranking.columns else "mean_diff"
    if "gene_id" not in ranking.columns or rank_col not in ranking.columns:
        raise ValueError(
            f"--ranking-table must contain 'gene_id' and '{rank_col}' columns; "
            f"found: {list(ranking.columns)}"
        )
    ranking = ranking[["gene_id", rank_col]].dropna().sort_values(rank_col, ascending=False)

    # Load external gene sets — never derive from the ranking/effect table.
    if args.gene_sets_gmt is not None:
        if not args.gene_sets_gmt.exists():
            raise FileNotFoundError(f"GMT file not found: {args.gene_sets_gmt}")
        gene_sets = load_gmt(args.gene_sets_gmt)
    else:
        gene_sets = load_category_files(args.fibrosis_genes, args.inflammation_genes, args.ecm_genes)

    if not gene_sets:
        raise RuntimeError(
            "No gene sets loaded. Provide --gene-sets-gmt or at least one of "
            "--fibrosis-genes / --inflammation-genes / --ecm-genes with ≥5 genes each.\n"
            "Gene sets must be external and independent of the ranking table — "
            "using effect-derived sets is circular and scientifically invalid."
        )

    gmt_path = out_dir / "external_gene_sets.gmt"
    write_gmt(gene_sets, gmt_path)

    pre_res = gp.prerank(
        rnk=ranking,
        gene_sets=str(gmt_path),
        min_size=5,
        max_size=500,
        permutation_num=1000,
        seed=42,
        outdir=str(out_dir),
        no_plot=False,
        format="png",
        verbose=False,
    )

    res = pre_res.res2d.copy().reset_index()
    # Normalise column names across gseapy versions.
    cols = {c.lower(): c for c in res.columns}
    rename: Dict[str, str] = {}
    if "term" in cols:
        rename[cols["term"]] = "pathway_name"
    if "es" in cols:
        rename[cols["es"]] = "ES"
    if "nes" in cols:
        rename[cols["nes"]] = "NES"
    if "nom p-val" in cols:
        rename[cols["nom p-val"]] = "nominal_pvalue"
    if "fdr q-val" in cols:
        rename[cols["fdr q-val"]] = "fdr_qvalue"
    if "lead_genes" in cols:
        rename[cols["lead_genes"]] = "leading_edge_genes"
    res = res.rename(columns=rename)

    # Attach category label based on the gene-set names that were loaded.
    known_labels = set(gene_sets.keys())
    if "pathway_name" in res.columns:
        res["gene_set_source"] = "curated_external"
        res["gene_set_loaded"] = res["pathway_name"].apply(
            lambda n: "yes" if n in known_labels else "no"
        )
    res.to_csv(out_dir / "gsea_results.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
