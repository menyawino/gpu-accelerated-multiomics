#!/usr/bin/env python3
"""Render original pathway-style gene maps for fibrosis, inflammation, and ECM remodeling.

This is an original schematic figure inspired by pathway-map layouts, not a reproduction
of proprietary KEGG artwork.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

PATHWAY_LAYOUT = {
    "Inflammation": [
        "TNF-alpha / NF-kappaB signaling",
        "IL6-JAK-STAT3 signaling",
        "Interferon gamma response",
        "Interferon alpha response",
        "Chemokine signaling",
        "Cytokine-cytokine receptor interaction",
        "Toll-like receptor signaling",
        "NLRP3 inflammasome activation",
    ],
    "Fibrosis": [
        "TGF-beta signaling",
        "Myofibroblast activation",
        "Collagen biosynthesis and crosslinking",
        "Connective tissue growth factor axis",
        "PDGF signaling",
        "SMAD-dependent profibrotic signaling",
        "Mesenchymal transition program",
        "Fibroblast proliferation and activation",
    ],
    "ECM Remodeling": [
        "Extracellular matrix organization",
        "ECM-receptor interaction",
        "Integrin signaling",
        "Focal adhesion remodeling",
        "Matrix metalloproteinase activity",
        "Proteoglycan remodeling",
        "Basement membrane reorganization",
        "Collagen degradation and turnover",
    ],
}

COLOR_MAP = {
    "Inflammation": "#c2185b",
    "Fibrosis": "#8e44ad",
    "ECM Remodeling": "#16a085",
}

PATHWAY_KEY = {
    "Inflammation": "inflammation",
    "Fibrosis": "fibrosis",
    "ECM Remodeling": "ecm_remodeling",
}


def assign_genes_to_modules(
    pathway_gene_effects: pd.DataFrame,
    category_key: str,
    pathway_names: list[str],
    genes_per_module: int = 3,
    gsea_results: pd.DataFrame | None = None,
) -> list[tuple[str, list[str]]]:
    modules = []

    if gsea_results is not None and not gsea_results.empty:
        gsea_subset = gsea_results[gsea_results["category"].str.lower().str.replace(" ", "_") == category_key].copy()
        gsea_subset = gsea_subset.set_index("pathway_name")
        for pathway_name in pathway_names:
            if pathway_name in gsea_subset.index:
                lead = str(gsea_subset.loc[pathway_name, "leading_edge_genes"])
                module_genes = [g for g in lead.split(";") if g][:genes_per_module]
            else:
                module_genes = []
            modules.append((pathway_name, module_genes))
        return modules

    subset = pathway_gene_effects[pathway_gene_effects["pathway"] == category_key].copy()
    subset = subset.sort_values(["fdr", "abs_mean_diff"], ascending=[True, False])
    genes = subset["gene_id"].tolist()[: len(pathway_names) * genes_per_module]
    idx = 0
    for pathway_name in pathway_names:
        module_genes = genes[idx : idx + genes_per_module]
        modules.append((pathway_name, module_genes))
        idx += genes_per_module
    return modules


def draw_category_panel(ax, title: str, modules: list[tuple[str, list[str]]], color: str) -> None:
    ax.set_axis_off()
    ax.text(0.5, 0.965, title, ha="center", va="top", fontsize=14, fontweight="bold", color=color, transform=ax.transAxes)

    positions = [
        (0.06, 0.72), (0.53, 0.72),
        (0.06, 0.50), (0.53, 0.50),
        (0.06, 0.28), (0.53, 0.28),
        (0.06, 0.06), (0.53, 0.06),
    ]
    box_w = 0.40
    box_h = 0.16

    for (x, y), (pathway_name, genes) in zip(positions, modules):
        patch = FancyBboxPatch(
            (x, y), box_w, box_h,
            boxstyle="round,pad=0.02,rounding_size=0.02",
            linewidth=1.2,
            edgecolor=color,
            facecolor="#fbfbfb",
            transform=ax.transAxes,
        )
        ax.add_patch(patch)
        ax.text(x + 0.02, y + box_h - 0.03, pathway_name, ha="left", va="top", fontsize=9.5, color=color, fontweight="bold", transform=ax.transAxes, wrap=True)
        gene_text = "\n".join(genes) if genes else "No genes assigned"
        ax.text(x + 0.02, y + box_h - 0.085, gene_text, ha="left", va="top", fontsize=8.5, color="#2c3e50", transform=ax.transAxes)

    # Add simple directional flow arrows to give pathway-map feel.
    arrow_pairs = [
        ((0.26, 0.70), (0.26, 0.64)),
        ((0.73, 0.70), (0.73, 0.64)),
        ((0.26, 0.48), (0.26, 0.42)),
        ((0.73, 0.48), (0.73, 0.42)),
        ((0.26, 0.26), (0.26, 0.20)),
        ((0.73, 0.26), (0.73, 0.20)),
    ]
    for start, end in arrow_pairs:
        ax.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=10,
                linewidth=1.0,
                color="#7f8c8d",
                alpha=0.7,
                transform=ax.transAxes,
            )
        )


def render() -> None:
    out_dir = Path("results/validation/reprogramming")
    pathway_gene_effects = pd.read_csv(out_dir / "pathway_genes_failing_vs_nonfailing.tsv", sep="\t")
    gsea_path = out_dir / "gsea" / "gsea_results.tsv"
    gsea_results = pd.read_csv(gsea_path, sep="\t") if gsea_path.exists() else None

    fig, axes = plt.subplots(1, 3, figsize=(18, 9))
    fig.suptitle("Pathway-Style Gene Maps for Cardiac Remodeling Programs", fontsize=17, y=0.995)

    gene_map_rows = []
    for ax, title in zip(axes, ["Inflammation", "Fibrosis", "ECM Remodeling"]):
        key = PATHWAY_KEY[title]
        modules = assign_genes_to_modules(
            pathway_gene_effects,
            key,
            PATHWAY_LAYOUT[title],
            genes_per_module=3,
            gsea_results=gsea_results,
        )
        for pathway_name, genes in modules:
            for gene in genes:
                gene_map_rows.append({"category": title, "pathway_name": pathway_name, "gene_id": gene})
        draw_category_panel(ax, title, modules, COLOR_MAP[title])

    fig.tight_layout(rect=[0, 0, 1, 0.975])
    fig.savefig(out_dir / "pathway_gene_map_figure.png", dpi=400, bbox_inches="tight")
    fig.savefig(out_dir / "pathway_gene_map_figure.pdf", bbox_inches="tight")
    pd.DataFrame(gene_map_rows).to_csv(out_dir / "pathway_gene_map_assignments.tsv", sep="\t", index=False)


if __name__ == "__main__":
    render()
