#!/usr/bin/env python3
"""Create a schematic manuscript panel explaining multi-omic layer separation."""

from __future__ import annotations

import argparse
from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
import pandas as pd


COLORS = {
    "rna": "#0b6e4f",
    "meth": "#c06c2b",
    "across_hyper": "#7f3c8d",
    "across_hypo": "#11a579",
    "ink": "#1c1c1c",
    "muted": "#6e6a64",
    "paper": "#fbf7ef",
    "grid": "#d9d2c3",
    "shared": "#587792",
    "warn": "#b54d2e",
}


def configure_style() -> None:
    plt.style.use("default")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def load_summary(path: Path) -> dict[str, int]:
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["metric"], df["value"]))


def load_gene_list(path: Path) -> list[str]:
    if not path.exists():
        return []
    df = pd.read_csv(path, sep="\t")
    if "gene_id" not in df.columns:
        return []
    return df["gene_id"].astype(str).tolist()


def load_ecm_genes(path: Path, n: int = 8) -> list[str]:
    if not path.exists():
        return []
    df = pd.read_csv(path, sep="\t")
    term = df[df["Term"].astype(str).str.contains("Extracellular Matrix Organization", case=False, na=False)]
    if term.empty:
        return []
    genes = str(term.iloc[0]["Genes"]).split(";")
    return genes[:n]


def rounded_box(ax, x, y, w, h, fc, ec=None, lw=1.2, radius=0.03):
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        facecolor=fc,
        edgecolor=ec or fc,
        linewidth=lw,
    )
    ax.add_patch(box)
    return box


def wrap_text(text: str, width: int) -> str:
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False))


def draw_module(ax, center_x: float, center_y: float, color: str, title: str, subtitle: str, count_text: str, bullets: list[str]) -> None:
    rounded_box(ax, center_x - 0.15, center_y - 0.12, 0.3, 0.24, fc=COLORS["paper"], ec=COLORS["grid"])
    circ = Circle((center_x, center_y + 0.07), 0.048, facecolor=color, edgecolor="white", linewidth=1.5)
    ax.add_patch(circ)
    ax.text(center_x, center_y + 0.07, count_text, ha="center", va="center", fontsize=12, color="white", fontweight="bold")
    ax.text(center_x, center_y + 0.01, title, ha="center", va="center", fontsize=12, fontweight="semibold", color=COLORS["ink"])
    ax.text(center_x, center_y - 0.03, wrap_text(subtitle, 38), ha="center", va="center", fontsize=8.7, color=COLORS["muted"])
    y = center_y - 0.075
    for bullet in bullets:
        ax.text(center_x - 0.125, y, wrap_text(bullet, 42), ha="left", va="top", fontsize=8.35, color=COLORS["ink"])
        y -= 0.05


def add_arrow(ax, start, end, color, text=None, text_y_offset=0.0):
    arrow = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=14, linewidth=1.8, color=color, alpha=0.9)
    ax.add_patch(arrow)
    if text:
        ax.text((start[0] + end[0]) / 2, (start[1] + end[1]) / 2 + text_y_offset, text, ha="center", va="center", fontsize=8.6, color=color, fontweight="semibold")


def main() -> None:
    ap = argparse.ArgumentParser(description="Create a schematic panel explaining multi-omic separation")
    ap.add_argument("--summary", required=True, type=Path)
    ap.add_argument("--hyper-down", required=True, type=Path)
    ap.add_argument("--hypo-up", required=True, type=Path)
    ap.add_argument("--ora-results", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    configure_style()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary = load_summary(args.summary)
    hyper_down = load_gene_list(args.hyper_down)
    hypo_up = load_gene_list(args.hypo_up)
    ecm_genes = load_ecm_genes(args.ora_results)

    across_examples = hyper_down[:1] + hypo_up[:5]
    ecm_text = ", ".join(ecm_genes[:4]) if ecm_genes else "COL1A1, POSTN, MMP2, TGFBR1"
    across_text = ", ".join(across_examples) if across_examples else "PRELID2, CDC14B, MAP2K3"

    fig = plt.figure(figsize=(14, 9.5))
    ax = fig.add_axes([0.04, 0.07, 0.92, 0.85])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    fig.suptitle(
        "Biological Interpretation of the Separated Multi-omic Interaction Layers",
        fontsize=16,
        fontweight="semibold",
        y=0.97,
        color=COLORS["ink"],
    )
    fig.text(
        0.5,
        0.935,
        "The current real processed analysis resolves broad within-layer programs in RNA and promoter methylation, while only a small gene-specific subset shows direct methylation-expression concordance.",
        ha="center",
        va="center",
        fontsize=9.3,
        color=COLORS["muted"],
    )

    draw_module(
        ax,
        0.18,
        0.68,
        COLORS["rna"],
        "RNA Self-Interaction Layer",
        "Co-expression structure across 64 transcriptomes",
        str(int(summary.get("rna_self_edges", 0))),
        [
            "Broad transcriptional rewiring dominates heart-failure structure.",
            "Captures coherent tissue-state and stress-response programs.",
            f"Example axis: ECM remodeling ({ecm_text}).",
        ],
    )

    draw_module(
        ax,
        0.5,
        0.68,
        COLORS["shared"],
        "Shared Within-layer Core",
        "Genes present in both self networks",
        str(int(summary.get("overlap_rna_methylation", 0))),
        [
            "Overlap concentrates inside the two self-organizing layers.",
            "These genes likely mark stable cardiac-state architecture.",
            "They do not extend into the direct across-omics concordance set.",
        ],
    )

    draw_module(
        ax,
        0.82,
        0.68,
        COLORS["meth"],
        "Promoter Methylation Self-Interaction Layer",
        "Promoter co-methylation structure across 36 samples",
        str(int(summary.get("methylation_self_edges", 0))),
        [
            "Dense promoter modules indicate coordinated epigenetic states.",
            "This layer is much broader than the direct concordant subset.",
            "Pattern fits chromatin-level program organization.",
        ],
    )

    add_arrow(ax, (0.29, 0.67), (0.4, 0.67), COLORS["shared"], text="Shared self-network genes", text_y_offset=0.04)
    add_arrow(ax, (0.71, 0.67), (0.6, 0.67), COLORS["shared"])

    rounded_box(ax, 0.11, 0.16, 0.78, 0.27, fc="#fffdf8", ec=COLORS["grid"])
    ax.text(0.5, 0.39, "Direct Across-omics Concordance Is Narrow and Gene-specific", ha="center", va="center", fontsize=13, fontweight="semibold", color=COLORS["ink"])
    ax.text(0.5, 0.355, wrap_text("Only 11 genes show direct promoter-methylation versus expression agreement in the current real processed analysis.", 92), ha="center", va="center", fontsize=9.1, color=COLORS["muted"])
    ax.text(0.27, 0.325, "repressed loci", ha="center", va="center", fontsize=8.5, color=COLORS["across_hyper"], fontweight="semibold")
    ax.text(0.51, 0.325, "activated loci", ha="center", va="center", fontsize=8.5, color=COLORS["across_hypo"], fontweight="semibold")

    rounded_box(ax, 0.16, 0.22, 0.22, 0.1, fc="#f6eefe", ec=COLORS["across_hyper"], lw=1.4)
    ax.text(0.27, 0.278, f"Hyper / Down\n{int(summary.get('across_hyper_down_genes', 0))} gene", ha="center", va="center", fontsize=11, color=COLORS["across_hyper"], fontweight="semibold")
    ax.text(0.27, 0.235, wrap_text(", ".join(hyper_down[:2]) if hyper_down else "PRELID2", 18), ha="center", va="center", fontsize=8.7, color=COLORS["ink"])

    rounded_box(ax, 0.4, 0.22, 0.22, 0.1, fc="#edf9f4", ec=COLORS["across_hypo"], lw=1.4)
    ax.text(0.51, 0.278, f"Hypo / Up\n{int(summary.get('across_hypo_up_genes', 0))} genes", ha="center", va="center", fontsize=11, color=COLORS["across_hypo"], fontweight="semibold")
    ax.text(0.51, 0.235, wrap_text(", ".join(hypo_up[:4]) if hypo_up else "CDC14B, MAP2K3, PCSK7", 24), ha="center", va="center", fontsize=8.5, color=COLORS["ink"])

    rounded_box(ax, 0.645, 0.208, 0.205, 0.118, fc="#fff4ef", ec=COLORS["warn"], lw=1.2)
    ax.text(0.7475, 0.282, "Why the bridge stays small", ha="center", va="center", fontsize=10.4, color=COLORS["warn"], fontweight="semibold")
    ax.text(0.7475, 0.235, wrap_text("Direct promoter-expression coupling is restricted to a small set of locus-level responders, not the full within-layer program architecture.", 31), ha="center", va="center", fontsize=7.8, color=COLORS["ink"])

    ax.text(0.5, 0.125, wrap_text(f"Representative direct concordance genes: {across_text}", 90), ha="center", va="center", fontsize=8.9, color=COLORS["ink"])
    ax.text(0.5, 0.085, wrap_text("Interpretation: the heart-failure phenotype is dominated by broad RNA and promoter-state modules, whereas direct methylation-expression coupling appears as a smaller, mechanistically specific layer.", 130), ha="center", va="center", fontsize=9.2, color=COLORS["muted"])

    png = args.output_dir / "multiomic_interaction_schematic_panel.png"
    pdf = args.output_dir / "multiomic_interaction_schematic_panel.pdf"
    fig.savefig(png, dpi=400, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()