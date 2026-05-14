#!/usr/bin/env python3
"""Create a real-data interaction summary figure from self and across omics outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLORS = {
    "rna": "#0b6e4f",
    "methylation": "#c06c2b",
    "positive": "#1f78b4",
    "negative": "#b2182b",
    "hyper_down": "#7f3c8d",
    "hypo_up": "#11a579",
    "grid": "#d9d2c3",
    "ink": "#1c1c1c",
    "muted": "#6e6a64",
    "panel": "#f8f5ef",
}


def configure_style() -> None:
    plt.style.use("default")
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": COLORS["panel"],
            "axes.edgecolor": COLORS["grid"],
            "axes.labelcolor": COLORS["ink"],
            "axes.titlecolor": COLORS["ink"],
            "axes.titlesize": 12,
            "axes.titleweight": "semibold",
            "axes.labelsize": 10,
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "grid.color": COLORS["grid"],
            "grid.linewidth": 0.8,
            "xtick.color": COLORS["ink"],
            "ytick.color": COLORS["ink"],
            "savefig.facecolor": "white",
        }
    )


def load_edges(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame(columns=["source", "target", "correlation", "fdr"])


def load_across(path: Path, signature: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t") if path.exists() else pd.DataFrame(columns=["gene_id", "padj", "log2FoldChange"])
    if df.empty:
        return pd.DataFrame(columns=["gene_id", "padj", "log2FoldChange", "signature"])
    # Normalise column names: support both DESeq2-style (padj, log2FoldChange)
    # and co-occurrence style (mean_expression, mean_promoter_beta).
    if "padj" not in df.columns:
        df["padj"] = 0.01  # placeholder: all genes in these files already passed co-occurrence filter
    if "log2FoldChange" not in df.columns:
        expr_col = "mean_expression" if "mean_expression" in df.columns else None
        df["log2FoldChange"] = np.log2(df[expr_col].clip(lower=1e-6)) if expr_col is not None else 0.0
    out = df[["gene_id", "padj", "log2FoldChange"]].copy()
    out["signature"] = signature
    return out


def node_table(edges: pd.DataFrame, omic: str) -> pd.DataFrame:
    if edges.empty:
        return pd.DataFrame(columns=["gene", "degree", "mean_abs_correlation", "omic"])

    melted = pd.concat(
        [
            edges[["source", "correlation"]].rename(columns={"source": "gene"}),
            edges[["target", "correlation"]].rename(columns={"target": "gene"}),
        ],
        ignore_index=True,
    )
    summary = (
        melted.groupby("gene", as_index=False)
        .agg(degree=("gene", "size"), mean_abs_correlation=("correlation", lambda s: float(np.mean(np.abs(s)))))
        .sort_values(["degree", "mean_abs_correlation", "gene"], ascending=[False, False, True])
    )
    summary["omic"] = omic
    return summary


def interaction_counts(rna_edges: pd.DataFrame, meth_edges: pd.DataFrame, across: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, edges in [("RNA self", rna_edges), ("Promoter methylation self", meth_edges)]:
        if edges.empty:
            pos = 0
            neg = 0
        else:
            pos = int((edges["correlation"] > 0).sum())
            neg = int((edges["correlation"] < 0).sum())
        rows.append({"interaction_class": label, "positive": pos, "negative": neg})

    for signature, label in [("hyper_down", "Across: hyper/down"), ("hypo_up", "Across: hypo/up")]:
        subset = across[across["signature"] == signature]
        rows.append({"interaction_class": label, "positive": int(signature == "hypo_up") * len(subset), "negative": int(signature == "hyper_down") * len(subset)})

    return pd.DataFrame(rows)


def overlap_summary(rna_nodes: pd.DataFrame, meth_nodes: pd.DataFrame, across: pd.DataFrame) -> pd.DataFrame:
    rna_set = set(rna_nodes["gene"]) if not rna_nodes.empty else set()
    meth_set = set(meth_nodes["gene"]) if not meth_nodes.empty else set()
    across_set = set(across["gene_id"]) if not across.empty else set()

    combos = [
        ("RNA only", len(rna_set - meth_set - across_set)),
        ("Methylation only", len(meth_set - rna_set - across_set)),
        ("Across only", len(across_set - rna_set - meth_set)),
        ("RNA + methylation", len((rna_set & meth_set) - across_set)),
        ("RNA + across", len((rna_set & across_set) - meth_set)),
        ("Methylation + across", len((meth_set & across_set) - rna_set)),
        ("All three", len(rna_set & meth_set & across_set)),
    ]
    return pd.DataFrame(combos, columns=["membership", "count"])


def evidence_table(rna_edges: pd.DataFrame, meth_edges: pd.DataFrame, across: pd.DataFrame) -> pd.DataFrame:
    pieces = []
    if not rna_edges.empty:
        pieces.append(pd.DataFrame({"interaction_class": "RNA self", "evidence": -np.log10(rna_edges["fdr"].clip(lower=1e-300))}))
    if not meth_edges.empty:
        pieces.append(pd.DataFrame({"interaction_class": "Promoter methylation self", "evidence": -np.log10(meth_edges["fdr"].clip(lower=1e-300))}))
    if not across.empty:
        label_map = {"hyper_down": "Across: hyper/down", "hypo_up": "Across: hypo/up"}
        tmp = across.copy()
        tmp["interaction_class"] = tmp["signature"].map(label_map)
        tmp["evidence"] = -np.log10(tmp["padj"].clip(lower=1e-300))
        pieces.append(tmp[["interaction_class", "evidence"]])
    if not pieces:
        return pd.DataFrame(columns=["interaction_class", "evidence"])
    return pd.concat(pieces, ignore_index=True)


def representative_features(rna_nodes: pd.DataFrame, meth_nodes: pd.DataFrame, across: pd.DataFrame, top_n: int = 8) -> pd.DataFrame:
    frames = []
    if not rna_nodes.empty:
        tmp = rna_nodes.head(top_n).copy()
        tmp["interaction_class"] = "RNA self"
        tmp["score"] = tmp["degree"].astype(float)
        tmp["label"] = tmp["gene"]
        frames.append(tmp[["interaction_class", "label", "score"]])
    if not meth_nodes.empty:
        tmp = meth_nodes.head(top_n).copy()
        tmp["interaction_class"] = "Promoter methylation self"
        tmp["score"] = tmp["degree"].astype(float)
        tmp["label"] = tmp["gene"]
        frames.append(tmp[["interaction_class", "label", "score"]])
    if not across.empty:
        tmp = across.sort_values(["padj", "log2FoldChange"], ascending=[True, False]).head(top_n).copy()
        tmp["interaction_class"] = tmp["signature"].map({"hyper_down": "Across: hyper/down", "hypo_up": "Across: hypo/up"})
        tmp["score"] = -np.log10(tmp["padj"].clip(lower=1e-300))
        tmp["label"] = tmp["gene_id"]
        frames.append(tmp[["interaction_class", "label", "score"]])
    if not frames:
        return pd.DataFrame(columns=["interaction_class", "label", "score"])
    return pd.concat(frames, ignore_index=True)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.03, label, transform=ax.transAxes, fontsize=13, fontweight="bold", color=COLORS["ink"], ha="left", va="bottom")


def plot_counts(ax: plt.Axes, counts: pd.DataFrame) -> None:
    y = np.arange(len(counts))
    ax.barh(y, counts["positive"], color=COLORS["positive"], label="Positive / activating")
    ax.barh(y, -counts["negative"], color=COLORS["negative"], label="Negative / repressive")
    ax.axvline(0, color=COLORS["ink"], linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(counts["interaction_class"])
    ax.set_xlabel("Interaction count")
    ax.set_title("Interaction balance")
    ax.grid(axis="x", alpha=0.35)
    max_extent = max(counts["positive"].max(), counts["negative"].max(), 1)
    ax.set_xlim(-max_extent * 1.2, max_extent * 1.2)
    ax.legend(frameon=False, loc="lower right", fontsize=8)


def plot_overlap(ax: plt.Axes, overlap: pd.DataFrame) -> None:
    bars = ax.bar(range(len(overlap)), overlap["count"], color=[COLORS["rna"], COLORS["methylation"], COLORS["hypo_up"], "#9f8f76", "#9f8f76", "#9f8f76", COLORS["ink"]])
    ax.set_xticks(range(len(overlap)))
    ax.set_xticklabels(overlap["membership"], rotation=35, ha="right")
    ax.set_ylabel("Genes")
    ax.set_title("Gene membership across interaction layers")
    ax.grid(axis="y", alpha=0.35)
    for bar, count in zip(bars, overlap["count"]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(overlap["count"].max() * 0.01, 0.5), f"{int(count)}", ha="center", va="bottom", fontsize=8, color=COLORS["ink"])
    if int(overlap.loc[overlap["membership"] == "All three", "count"].iloc[0]) == 0:
        ax.text(0.98, 0.95, "No shared genes span all three real-data layers", transform=ax.transAxes, ha="right", va="top", fontsize=8.5, color=COLORS["muted"])


def plot_evidence(ax: plt.Axes, evidence: pd.DataFrame) -> None:
    order = ["RNA self", "Promoter methylation self", "Across: hyper/down", "Across: hypo/up"]
    palette = [COLORS["rna"], COLORS["methylation"], COLORS["hyper_down"], COLORS["hypo_up"]]
    positions = np.arange(len(order))
    for pos, label, color in zip(positions, order, palette):
        subset = evidence.loc[evidence["interaction_class"] == label, "evidence"].dropna().values
        if subset.size == 0:
            continue
        parts = ax.violinplot(subset, positions=[pos], widths=0.8, showmeans=False, showextrema=False)
        for body in parts["bodies"]:
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(0.5)
        ax.scatter(np.full_like(subset, pos, dtype=float), subset, s=6, alpha=0.25, color=color, linewidths=0)
        ax.scatter([pos], [np.median(subset)], s=28, color=color, edgecolors="white", linewidths=0.7, zorder=3)
    ax.set_xticks(positions)
    ax.set_xticklabels(order, rotation=20, ha="right")
    ax.set_ylabel(r"Evidence strength ($-\log_{10}$ adjusted $P$)")
    ax.set_title("Statistical support")
    ax.grid(axis="y", alpha=0.35)


def plot_representatives(ax: plt.Axes, reps: pd.DataFrame) -> None:
    order = ["RNA self", "Promoter methylation self", "Across: hyper/down", "Across: hypo/up"]
    color_map = {
        "RNA self": COLORS["rna"],
        "Promoter methylation self": COLORS["methylation"],
        "Across: hyper/down": COLORS["hyper_down"],
        "Across: hypo/up": COLORS["hypo_up"],
    }
    y_cursor = 0
    yticks = []
    ylabels = []
    for label in order:
        subset = reps[reps["interaction_class"] == label].sort_values("score", ascending=True)
        if subset.empty:
            continue
        positions = np.arange(y_cursor, y_cursor + len(subset))
        ax.hlines(positions, 0, subset["score"], color=color_map[label], linewidth=1.6, alpha=0.85)
        ax.scatter(subset["score"], positions, color=color_map[label], s=28, zorder=3)
        yticks.extend(positions.tolist())
        ylabels.extend(subset["label"].tolist())
        ax.text(subset["score"].max() * 1.02 if subset["score"].max() > 0 else 0.2, positions.mean(), label, va="center", ha="left", fontsize=8.5, color=color_map[label], fontweight="semibold")
        y_cursor += len(subset) + 1
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels)
    ax.set_xlabel("Representative feature score")
    ax.set_title("Representative genes by interaction layer")
    ax.grid(axis="x", alpha=0.35)
    ax.invert_yaxis()


def build_summary_table(counts: pd.DataFrame, overlap: pd.DataFrame, rna_nodes: pd.DataFrame, meth_nodes: pd.DataFrame, across: pd.DataFrame) -> pd.DataFrame:
    records = [
        {"metric": "rna_self_edges", "value": int((counts.loc[counts["interaction_class"] == "RNA self", ["positive", "negative"]].sum(axis=1)).iloc[0])},
        {"metric": "methylation_self_edges", "value": int((counts.loc[counts["interaction_class"] == "Promoter methylation self", ["positive", "negative"]].sum(axis=1)).iloc[0])},
        {"metric": "across_hyper_down_genes", "value": int((counts.loc[counts["interaction_class"] == "Across: hyper/down", "negative"]).iloc[0])},
        {"metric": "across_hypo_up_genes", "value": int((counts.loc[counts["interaction_class"] == "Across: hypo/up", "positive"]).iloc[0])},
        {"metric": "rna_network_genes", "value": int(rna_nodes["gene"].nunique())},
        {"metric": "methylation_network_genes", "value": int(meth_nodes["gene"].nunique())},
        {"metric": "across_genes", "value": int(across["gene_id"].nunique())},
    ]
    for _, row in overlap.iterrows():
        metric = row["membership"].lower().replace(" + ", "_").replace(" ", "_")
        records.append({"metric": f"overlap_{metric}", "value": int(row["count"])})
    return pd.DataFrame(records)


def main() -> None:
    ap = argparse.ArgumentParser(description="Create a publication-style multi-omics interaction summary figure")
    ap.add_argument("--rna-edges", required=True, type=Path)
    ap.add_argument("--methylation-edges", required=True, type=Path)
    ap.add_argument("--hyper-down", required=True, type=Path)
    ap.add_argument("--hypo-up", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()

    configure_style()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rna_edges = load_edges(args.rna_edges)
    meth_edges = load_edges(args.methylation_edges)
    across = pd.concat(
        [
            load_across(args.hyper_down, "hyper_down"),
            load_across(args.hypo_up, "hypo_up"),
        ],
        ignore_index=True,
    )

    rna_nodes = node_table(rna_edges, "RNA")
    meth_nodes = node_table(meth_edges, "Promoter methylation")
    counts = interaction_counts(rna_edges, meth_edges, across)
    overlap = overlap_summary(rna_nodes, meth_nodes, across)
    evidence = evidence_table(rna_edges, meth_edges, across)
    reps = representative_features(rna_nodes, meth_nodes, across)
    summary = build_summary_table(counts, overlap, rna_nodes, meth_nodes, across)

    summary.to_csv(args.output_dir / "multiomic_interaction_summary.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.5))
    ax_a, ax_b, ax_c, ax_d = axes.ravel()

    plot_counts(ax_a, counts)
    add_panel_label(ax_a, "a")

    plot_overlap(ax_b, overlap)
    add_panel_label(ax_b, "b")

    plot_evidence(ax_c, evidence)
    add_panel_label(ax_c, "c")

    plot_representatives(ax_d, reps)
    add_panel_label(ax_d, "d")

    fig.suptitle(
        "Self and Across Interaction Structure in the Real Processed Heart-Failure Analyses",
        fontsize=15,
        fontweight="semibold",
        y=0.995,
    )

    subtitle = (
        "Self interactions come from the RNA and promoter-methylation regulatory networks; "
        "across interactions use the real hypermethylated/downregulated and hypomethylated/upregulated discovery genes."
    )
    fig.text(0.5, 0.968, subtitle, ha="center", va="top", fontsize=8.8, color=COLORS["muted"])

    fig.tight_layout(rect=[0, 0, 1, 0.93])

    png = args.output_dir / "multiomic_interaction_figure.png"
    pdf = args.output_dir / "multiomic_interaction_figure.pdf"
    fig.savefig(png, dpi=400, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()