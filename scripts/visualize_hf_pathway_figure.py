#!/usr/bin/env python3
"""Create a composite HF vs non-HF pathway figure modeled on the reference layout."""

import argparse
from pathlib import Path
from textwrap import fill

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.stats import mannwhitneyu


GROUP_ORDER = ["NF", "DCM", "ICM"]
REPRESENTATIVE_PATHWAYS = [
    "Interferon Alpha Response",
    "Extracellular Matrix Organization R-HSA-1474244",
    "Citric Acid (TCA) Cycle And Respiratory Electron Transport R-HSA-1428517",
]


def configure_style() -> dict[str, str]:
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 12.5,
            "axes.labelsize": 9,
            "axes.titleweight": "semibold",
            "figure.titlesize": 16,
            "figure.titleweight": "semibold",
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "savefig.dpi": 400,
        }
    )
    return {
        "NF": "#4ca6a8",
        "DCM": "#e6a23c",
        "ICM": "#b33b46",
        "low": "#364b9a",
        "mid": "#f7f7f7",
        "high": "#d73027",
    }


def load_inputs(ora_path: Path, expression_path: Path, metadata_path: Path, de_path: Path):
    ora = pd.read_csv(ora_path, sep="\t")
    expr = pd.read_csv(expression_path, sep="\t", index_col=0)
    meta = pd.read_csv(metadata_path, sep="\t")
    de = pd.read_csv(de_path, sep="\t")
    return ora, expr, meta, de


def clean_term_label(term: str) -> str:
    replacements = {
        "Extracellular Matrix Organization R-HSA-1474244": "Extracellular Matrix Organization",
        "Citric Acid (TCA) Cycle And Respiratory Electron Transport R-HSA-1428517": "TCA Cycle And Respiratory Electron Transport",
        "Respiratory Electron Transport, ATP Synthesis By Chemiosmotic Coupling, Heat Production By Uncoupling Proteins R-HSA-163200": "Respiratory Electron Transport",
        "Interferon Alpha Response": "Interferon Alpha Response",
    }
    return replacements.get(term, term)


def phenotype_sample_map(meta: pd.DataFrame, expr: pd.DataFrame) -> dict[str, list[str]]:
    sample_col = meta.columns[0]
    phenotype_col = "phenotype" if "phenotype" in meta.columns else meta.columns[1]
    mapping = {group: [] for group in GROUP_ORDER}
    for _, row in meta.iterrows():
        sample = str(row[sample_col])
        phenotype = str(row[phenotype_col])
        if phenotype in mapping and sample in expr.columns:
            mapping[phenotype].append(sample)
    return mapping


def score_gene_set(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    valid = [gene for gene in genes if gene in expr.index]
    if not valid:
        return pd.Series(dtype=float)
    subset = expr.loc[valid]
    centered = subset.sub(subset.mean(axis=1), axis=0)
    scaled = centered.div(subset.std(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    return scaled.mean(axis=0)


def pathway_pvalue(scores: pd.Series, sample_groups: dict[str, list[str]]) -> float:
    hf = scores[sample_groups["DCM"] + sample_groups["ICM"]].dropna().values.astype(float)
    nf = scores[sample_groups["NF"]].dropna().values.astype(float)
    if len(hf) < 2 or len(nf) < 2:
        return np.nan
    return float(mannwhitneyu(hf, nf, alternative="two-sided").pvalue)


def select_hallmark_pathways(ora: pd.DataFrame) -> pd.DataFrame:
    hallmark = ora[ora["Gene_set"] == "MSigDB_Hallmark_2020"].copy()
    selected = []
    for direction, count in [("hf_up", 6), ("hf_down", 6)]:
        top = hallmark[hallmark["direction"] == direction].sort_values("Adjusted P-value").head(count)
        selected.append(top)
    out = pd.concat(selected, ignore_index=True)
    out["display_term"] = out["Term"].str.replace(" Signaling", " signaling", regex=False)
    return out


def build_score_matrix(selected: pd.DataFrame, expr: pd.DataFrame, sample_groups: dict[str, list[str]]) -> tuple[pd.DataFrame, dict[str, pd.Series], dict[str, float]]:
    matrix_rows = []
    scores_by_term = {}
    pvalues = {}
    for _, row in selected.iterrows():
        genes = [gene for gene in str(row["Genes"]).split(";") if gene]
        scores = score_gene_set(expr, genes)
        scores_by_term[row["Term"]] = scores
        pvalues[row["Term"]] = pathway_pvalue(scores, sample_groups)
        group_means = {group: float(scores[sample_groups[group]].mean()) for group in GROUP_ORDER}
        matrix_rows.append({"Term": row["Term"], **group_means})
    matrix = pd.DataFrame(matrix_rows).set_index("Term")
    matrix = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    return matrix, scores_by_term, pvalues


def select_rows_by_terms(ora: pd.DataFrame, terms: list[str]) -> pd.DataFrame:
    selected = []
    for term in terms:
        matches = ora[ora["Term"] == term].sort_values("Adjusted P-value")
        if matches.empty:
            raise ValueError(f"Required pathway term not found in ORA results: {term}")
        selected.append(matches.iloc[0])
    return pd.DataFrame(selected).reset_index(drop=True)


def top_genes_for_term(term: str, selected: pd.DataFrame, de: pd.DataFrame, n_genes: int = 18) -> list[str]:
    gene_str = selected.loc[selected["Term"] == term, "Genes"].iloc[0]
    genes = [gene for gene in str(gene_str).split(";") if gene]
    de_sub = de[de["gene_id"].isin(genes)].copy()
    if de_sub.empty:
        return genes[:n_genes]
    de_sub["abs_mean_diff"] = de_sub["mean_diff"].abs()
    return de_sub.sort_values(["abs_mean_diff", "fdr"], ascending=[False, True])["gene_id"].head(n_genes).tolist()


def gene_group_heatmap(expr: pd.DataFrame, genes: list[str], sample_groups: dict[str, list[str]]) -> pd.DataFrame:
    valid = [gene for gene in genes if gene in expr.index]
    group_means = pd.DataFrame(
        {
            group: expr.loc[valid, sample_groups[group]].mean(axis=1)
            for group in GROUP_ORDER
        }
    )
    scaled = group_means.sub(group_means.mean(axis=1), axis=0).div(group_means.std(axis=1).replace(0, np.nan), axis=0).fillna(0.0)
    return scaled


def significance_marker(pvalue: float) -> str:
    if np.isnan(pvalue):
        return ""
    if pvalue < 0.001:
        return "***"
    if pvalue < 0.01:
        return "**"
    if pvalue < 0.05:
        return "*"
    return ""


def add_panel_label(ax: plt.Axes, label: str, x: float = -0.06, y: float = 1.05) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=15, fontweight="bold", va="bottom")


def wrap_term(term: str, width: int = 28) -> str:
    return fill(clean_term_label(term), width=width)


def plot_composite(selected: pd.DataFrame, matrix: pd.DataFrame, scores_by_term: dict[str, pd.Series], pvalues: dict[str, float], representative: pd.DataFrame, representative_scores: dict[str, pd.Series], representative_pvalues: dict[str, float], expr: pd.DataFrame, sample_groups: dict[str, list[str]], de: pd.DataFrame, out_dir: Path, palette: dict[str, str]) -> None:
    row_linkage = linkage(matrix.values, method="average", metric="euclidean")
    row_order = dendrogram(row_linkage, no_plot=True)["leaves"]
    ordered_terms = matrix.index[row_order].tolist()
    ordered_matrix = matrix.loc[ordered_terms]

    figure = plt.figure(figsize=(15.5, 10.4))
    grid = figure.add_gridspec(3, 3, width_ratios=[1.55, 0.92, 0.9], height_ratios=[1, 1, 1], wspace=0.42, hspace=0.46)

    left_grid = grid[:, 0].subgridspec(1, 2, width_ratios=[0.16, 0.84], wspace=0.035)
    dend_ax = figure.add_subplot(left_grid[0, 0])
    heat_ax = figure.add_subplot(left_grid[0, 1])

    dendrogram(row_linkage, orientation="left", no_labels=True, color_threshold=0, above_threshold_color="#2c3e50", ax=dend_ax)
    dend_ax.invert_yaxis()
    dend_ax.axis("off")

    heatmap = heat_ax.imshow(ordered_matrix.values, aspect="auto", cmap="coolwarm", vmin=-1.8, vmax=1.8)
    heat_ax.set_xticks(np.arange(len(GROUP_ORDER)))
    heat_ax.set_xticklabels(GROUP_ORDER)
    heat_ax.set_yticks(np.arange(len(ordered_terms)))
    heat_ax.set_yticklabels([wrap_term(term) for term in ordered_terms])
    heat_ax.tick_params(axis="y", length=0, labelsize=8.5, pad=2)
    heat_ax.set_title("Hallmark Pathway Score Shifts", pad=10)
    for i, term in enumerate(ordered_terms):
        marker = significance_marker(pvalues.get(term, np.nan))
        if marker:
            heat_ax.text(2.45, i, marker, va="center", ha="left", fontsize=9, color="#111111")
    add_panel_label(heat_ax, "A", x=-0.08, y=1.06)

    cbar = figure.colorbar(heatmap, ax=heat_ax, fraction=0.046, pad=0.03)
    cbar.set_label("Group Mean Pathway Z-score")
    cbar.ax.tick_params(length=0)

    panel_labels = [("B", "C"), ("D", "E"), ("F", "G")]
    for row_index, term in enumerate(REPRESENTATIVE_PATHWAYS):
        violin_ax = figure.add_subplot(grid[row_index, 1])
        gene_ax = figure.add_subplot(grid[row_index, 2])

        pathway_scores = representative_scores[term]
        violin_data = [pathway_scores[sample_groups[group]].dropna().values.astype(float) for group in GROUP_ORDER]
        parts = violin_ax.violinplot(violin_data, showmeans=True, showextrema=False, widths=0.82)
        for body, group in zip(parts["bodies"], GROUP_ORDER):
            body.set_facecolor(palette[group])
            body.set_alpha(0.75)
        parts["cmeans"].set_color("#1f2933")
        for pos, group in enumerate(GROUP_ORDER, start=1):
            jitter = np.linspace(-0.06, 0.06, len(violin_data[pos - 1])) if len(violin_data[pos - 1]) > 1 else np.array([0.0])
            violin_ax.scatter(np.full(len(violin_data[pos - 1]), pos) + jitter, violin_data[pos - 1], s=12, color="#1f2933", alpha=0.45)
        violin_ax.set_xticks([1, 2, 3])
        violin_ax.set_xticklabels(GROUP_ORDER)
        violin_ax.set_ylabel("Pathway score")
        violin_ax.set_title(fill(clean_term_label(term), width=32), pad=14)
        marker = significance_marker(representative_pvalues.get(term, np.nan))
        if marker:
            violin_ax.text(0.98, 0.98, marker, transform=violin_ax.transAxes, ha="right", va="top", fontsize=11)
        add_panel_label(violin_ax, panel_labels[row_index][0], x=-0.08, y=1.08)
        violin_ax.spines["top"].set_visible(False)
        violin_ax.spines["right"].set_visible(False)

        genes = top_genes_for_term(term, representative, de)
        gene_matrix = gene_group_heatmap(expr, genes, sample_groups)
        gene_ax.imshow(gene_matrix.values, aspect="auto", cmap="coolwarm", vmin=-2.0, vmax=2.0)
        gene_ax.set_xticks(np.arange(len(GROUP_ORDER)))
        gene_ax.set_xticklabels(GROUP_ORDER)
        gene_ax.set_yticks(np.arange(gene_matrix.shape[0]))
        gene_ax.set_yticklabels(gene_matrix.index.tolist())
        gene_ax.tick_params(axis="y", labelsize=8.5, length=0)
        gene_ax.set_title("Leading Genes", pad=14)
        add_panel_label(gene_ax, panel_labels[row_index][1], x=-0.08, y=1.08)
        gene_ax.spines["top"].set_visible(False)
        gene_ax.spines["right"].set_visible(False)

    figure.suptitle("HF Pathway Remodeling Across Non-failing, DCM, and ICM Samples", y=0.985)
    out_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(str(out_dir / "hf_vs_nonhf_pathway_multiplot.png"), bbox_inches="tight")
    figure.savefig(str(out_dir / "hf_vs_nonhf_pathway_multiplot.pdf"), bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a composite HF pathway figure")
    parser.add_argument("--ora-results", type=Path, default=Path("results/real_processed/hf_vs_nonhf_pathways/ora_results.tsv"))
    parser.add_argument("--expression", type=Path, default=Path("results/real_processed/inputs/gse116250_expression.tsv"))
    parser.add_argument("--metadata", type=Path, default=Path("results/real_processed/inputs/gse116250_metadata.tsv"))
    parser.add_argument("--de-table", type=Path, default=Path("results/real_processed/reprogramming/all_genes_failing_vs_nonfailing.tsv"))
    parser.add_argument("--output-dir", type=Path, default=Path("results/real_processed/hf_vs_nonhf_pathways"))
    args = parser.parse_args()

    palette = configure_style()
    ora, expr, meta, de = load_inputs(args.ora_results, args.expression, args.metadata, args.de_table)
    sample_groups = phenotype_sample_map(meta, expr)
    selected = select_hallmark_pathways(ora)
    matrix, scores_by_term, pvalues = build_score_matrix(selected, expr, sample_groups)
    representative = select_rows_by_terms(ora, REPRESENTATIVE_PATHWAYS)
    representative_matrix, representative_scores, representative_pvalues = build_score_matrix(representative, expr, sample_groups)
    plot_composite(selected, matrix, scores_by_term, pvalues, representative, representative_scores, representative_pvalues, expr, sample_groups, de, args.output_dir, palette)

    print(f"Wrote composite pathway figure to: {args.output_dir}")


if __name__ == "__main__":
    main()