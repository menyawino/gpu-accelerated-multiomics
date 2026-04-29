#!/usr/bin/env python3
"""Render a labeled figure of canonical pathway names for fibrosis, inflammation, and ECM remodeling."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyBboxPatch


PATHWAY_MAP = {
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
        "Epithelial/Endothelial-to-mesenchymal transition",
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

COLORS = {
    "Inflammation": "#c2185b",
    "Fibrosis": "#8e44ad",
    "ECM Remodeling": "#16a085",
}


def write_catalog(out_dir: Path) -> Path:
    rows = []
    for category, pathways in PATHWAY_MAP.items():
        for pathway in pathways:
            rows.append({"category": category, "pathway_name": pathway})
    out_path = out_dir / "pathway_name_catalog.tsv"
    pd.DataFrame(rows).to_csv(out_path, sep="\t", index=False)
    return out_path


def render(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    write_catalog(out_dir)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "figure.titlesize": 16,
            "axes.titlesize": 12,
        }
    )

    fig, axes = plt.subplots(1, 3, figsize=(15, 7))
    fig.suptitle("Canonical Pathway Names by Remodeling Program", y=0.98)

    for ax, (category, pathways) in zip(axes, PATHWAY_MAP.items()):
        ax.set_axis_off()
        patch = FancyBboxPatch(
            (0.02, 0.02),
            0.96,
            0.96,
            boxstyle="round,pad=0.02,rounding_size=0.03",
            linewidth=1.5,
            edgecolor=COLORS[category],
            facecolor="#fafafa",
            transform=ax.transAxes,
        )
        ax.add_patch(patch)
        ax.text(
            0.5,
            0.93,
            category,
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold",
            color=COLORS[category],
            transform=ax.transAxes,
        )

        y = 0.84
        for pathway in pathways:
            ax.text(
                0.07,
                y,
                f"• {pathway}",
                ha="left",
                va="top",
                fontsize=10,
                color="#2c3e50",
                transform=ax.transAxes,
                wrap=True,
            )
            y -= 0.095

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_dir / "pathway_name_figure.png", dpi=400, bbox_inches="tight")
    fig.savefig(out_dir / "pathway_name_figure.pdf", bbox_inches="tight")


if __name__ == "__main__":
    render(Path("results/validation/reprogramming"))
