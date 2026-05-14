#!/usr/bin/env python3
"""
Create publication-ready figures for each validation pipeline step.

Generates four figure panels:
1) Arm 1 discovery
2) Arm 1 validation
3) Arm 2 discovery
4) Arm 2 validation
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def configure_style() -> dict:
    """Set plotting defaults and return a consistent color palette."""
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "figure.titlesize": 14,
            "legend.fontsize": 9,
            "savefig.dpi": 400,
        }
    )
    return {
        "hyper": "#c0392b",
        "hypo": "#1f78b4",
        "neutral": "#7f8c8d",
        "oxidative": "#0b6e4f",
        "glycolytic": "#d35400",
        "fibrosis": "#8e44ad",
        "inflammation": "#c2185b",
        "ecm_remodeling": "#16a085",
        "other": "#95a5a6",
    }


def save_figure(fig: plt.Figure, out_dir: Path, stem: str) -> None:
    """Save a figure as both PNG and PDF."""
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{stem}.png", bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def load_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def plot_arm1_discovery(data_dir: Path, inputs_dir: Path, out_dir: Path, palette: dict) -> None:
    hyper = load_tsv(data_dir / "arm1_hypermethylated.tsv")
    hypo = load_tsv(data_dir / "arm1_hypomethylated.tsv")
    meth = pd.read_csv(inputs_dir / "gse123976_methylation.tsv", sep="\t", index_col=0)

    mean_beta_all = meth.mean(axis=1)
    non_sig = mean_beta_all[~mean_beta_all.index.isin(set(hyper["gene_id"]).union(set(hypo["gene_id"])))].values

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Compute per-gene mean beta across HF+NF for violin; actual columns are hf_mean_beta / nf_mean_beta
    hyper["mean_beta"] = (hyper["hf_mean_beta"] + hyper["nf_mean_beta"]) / 2
    hypo["mean_beta"] = (hypo["hf_mean_beta"] + hypo["nf_mean_beta"]) / 2
    violin_data = [hyper["mean_beta"].values, hypo["mean_beta"].values, non_sig]
    vp = axes[0].violinplot(violin_data, showmeans=True, showextrema=False, widths=0.8)
    for body, c in zip(vp["bodies"], [palette["hyper"], palette["hypo"], palette["neutral"]]):
        body.set_facecolor(c)
        body.set_alpha(0.5)
    vp["cmeans"].set_color("#2c3e50")
    axes[0].set_xticks([1, 2, 3])
    axes[0].set_xticklabels(["Hypermethylated", "Hypomethylated", "Background"])
    axes[0].set_ylabel("Mean promoter beta")
    axes[0].set_title("Discovery methylation distributions")
    axes[0].axhline(0.5, linestyle="--", color="#2c3e50", linewidth=1)

    top_n = 20
    top_hyper = hyper.nlargest(top_n, "mean_beta")
    top_hypo = hypo.nsmallest(top_n, "mean_beta")
    axes[1].scatter(
        top_hyper["delta_beta"],
        top_hyper["mean_beta"],
        s=40,
        color=palette["hyper"],
        label=f"Top {top_n} hyper",
        alpha=0.9,
    )
    axes[1].scatter(
        top_hypo["delta_beta"],
        top_hypo["mean_beta"],
        s=40,
        color=palette["hypo"],
        label=f"Top {top_n} hypo",
        alpha=0.9,
    )
    axes[1].set_xlabel("Delta beta (HF - NF)")
    axes[1].set_ylabel("Mean promoter beta")
    axes[1].set_title("Most extreme discovery signatures")
    axes[1].legend(frameon=True)

    fig.suptitle("Arm 1 Discovery: GSE123976 Methylation Signatures", y=1.03)
    fig.tight_layout()
    save_figure(fig, out_dir, "figure1_arm1_discovery")


def plot_arm1_validation(data_dir: Path, out_dir: Path, palette: dict) -> None:
    concordance_all = load_tsv(data_dir / "arm1_validation_concordance.tsv")
    # Use only combined-platform rows for summary bar charts
    concordance = concordance_all[concordance_all["platform"] == "combined"].copy()
    metabolic = load_tsv(data_dir / "arm1_validation_metabolic.tsv")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    dir_rates = (
        concordance.groupby("signature_direction")["concordant"]
        .agg(["sum", "count"])
        .reset_index()
    )
    dir_rates["rate"] = 100.0 * dir_rates["sum"] / dir_rates["count"]

    colors = [palette["hyper"], palette["hypo"]]
    axes[0].bar(dir_rates["signature_direction"], dir_rates["rate"], color=colors, alpha=0.85)
    axes[0].set_ylim(0, 105)
    axes[0].set_ylabel("Concordance (%)")
    axes[0].set_title("Direction concordance")
    for i, r in dir_rates.iterrows():
        axes[0].text(i, r["rate"] + 1.5, f"{int(r['sum'])}/{int(r['count'])}", ha="center", fontsize=9)

    beta_color_map = {
        "hypermethylated_signature": palette["hyper"],
        "hypomethylated_signature": palette["hypo"],
        "oxidative_phosphorylation": palette["oxidative"],
        "glycolysis": palette["glycolytic"],
    }
    beta_colors = [beta_color_map.get(p, palette["other"]) for p in metabolic["metabolic_program"]]
    axes[1].bar(metabolic["metabolic_program"], metabolic["mean_beta"], color=beta_colors, alpha=0.9)
    axes[1].axhline(0.5, linestyle="--", color="#2c3e50", linewidth=1)
    axes[1].set_ylabel("Mean beta")
    axes[1].set_title("Metabolic program methylation")
    axes[1].tick_params(axis="x", rotation=15)

    comp = metabolic[["metabolic_program", "hypermethylated_count", "hypomethylated_count", "neutral_count"]].copy()
    x = np.arange(len(comp))
    axes[2].bar(x, comp["hypermethylated_count"], color=palette["hyper"], label="Hyper")
    axes[2].bar(x, comp["hypomethylated_count"], bottom=comp["hypermethylated_count"], color=palette["hypo"], label="Hypo")
    axes[2].bar(
        x,
        comp["neutral_count"],
        bottom=comp["hypermethylated_count"] + comp["hypomethylated_count"],
        color=palette["neutral"],
        label="Neutral",
    )
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(comp["metabolic_program"], rotation=15)
    axes[2].set_ylabel("CpG count")
    axes[2].set_title("Program-level methylation composition")
    axes[2].legend(frameon=True)

    fig.suptitle("Arm 1 Validation: GSE197670 Replication", y=1.03)
    fig.tight_layout()
    save_figure(fig, out_dir, "figure2_arm1_validation")


def plot_arm1_concordance_venn(data_dir: Path, inputs_dir: Path, out_dir: Path, palette: dict) -> None:
    hyper = load_tsv(data_dir / "arm1_hypermethylated.tsv")
    hypo = load_tsv(data_dir / "arm1_hypomethylated.tsv")
    concordance_all = load_tsv(data_dir / "arm1_validation_concordance.tsv")
    concordance = concordance_all[concordance_all["platform"] == "combined"].copy()
    validation_meth = pd.read_csv(inputs_dir / "gse197670_methylation.tsv", sep="\t", index_col=0)
    validation_genes = set(validation_meth.index.astype(str))

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

    specs = [
        ("hypermethylated", hyper, palette["hyper"], "Hypermethylated signature"),
        ("hypomethylated", hypo, palette["hypo"], "Hypomethylated signature"),
    ]

    for ax, (direction, sig_df, color, title) in zip(axes, specs):
        discovered = set(sig_df["gene_id"].astype(str))
        tested_df = concordance[concordance["signature_direction"] == direction].copy()
        tested = set(tested_df["gene_id"].astype(str))
        concordant = set(tested_df.loc[tested_df["concordant"].astype(bool), "gene_id"].astype(str))
        not_tested = len(discovered - validation_genes)
        tested_count = len(tested)
        concordant_count = len(concordant)
        discordant_count = tested_count - concordant_count

        ax.set_xlim(0, 10)
        ax.set_ylim(0, 8)
        ax.set_aspect("equal")
        ax.axis("off")

        left = Circle((4.0, 4.0), 2.55, facecolor=color, edgecolor=color, alpha=0.32, linewidth=2)
        right = Circle((6.0, 4.0), 2.55, facecolor="#95a5a6", edgecolor="#7f8c8d", alpha=0.24, linewidth=2)
        ax.add_patch(left)
        ax.add_patch(right)

        ax.text(3.0, 6.85, title, ha="center", va="bottom", fontsize=11, color=color, fontweight="bold")
        ax.text(7.0, 6.85, "Validation matrix", ha="center", va="bottom", fontsize=11, color="#2c3e50", fontweight="bold")

        ax.text(2.4, 4.1, f"Not tested\n{not_tested}", ha="center", va="center", fontsize=13, color="#2c3e50")
        ax.text(5.0, 4.55, f"Tested\n{tested_count}", ha="center", va="center", fontsize=13, color="#2c3e50", fontweight="bold")
        ax.text(5.0, 3.15, f"Concordant\n{concordant_count}", ha="center", va="center", fontsize=12, color="#2c3e50")
        ax.text(5.0, 1.25, f"Discordant among tested: {discordant_count}", ha="center", va="center", fontsize=10, color="#2c3e50")
        ax.text(5.0, 0.55, f"Discovery total: {len(discovered)}", ha="center", va="center", fontsize=10, color="#2c3e50")

    fig.suptitle("Arm 1 Concordance Venn-Style Summary", y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    save_figure(fig, out_dir, "figure3_arm1_concordance_venn")


def plot_arm2_discovery(data_dir: Path, out_dir: Path, palette: dict) -> None:
    hyper_down = load_tsv(data_dir / "arm2_hyper_down.tsv")
    hypo_up = load_tsv(data_dir / "arm2_hypo_up.tsv")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    axes[0].scatter(
        hyper_down["mean_promoter_beta"],
        hyper_down["mean_expression"],
        s=24,
        alpha=0.7,
        color=palette["hyper"],
        label="Hyper-down",
    )
    axes[0].scatter(
        hypo_up["mean_promoter_beta"],
        hypo_up["mean_expression"],
        s=24,
        alpha=0.7,
        color=palette["hypo"],
        label="Hypo-up",
    )
    axes[0].axvline(0.5, linestyle="--", color="#2c3e50", linewidth=1)
    axes[0].set_xlabel("Mean promoter beta")
    axes[0].set_ylabel("Mean expression")
    axes[0].set_title("Joint methylation-expression space")
    axes[0].legend(frameon=True)

    bins = np.linspace(0, 1, 25)
    axes[1].hist(hyper_down["mean_promoter_beta"], bins=bins, alpha=0.75, color=palette["hyper"], label="Hyper-down")
    axes[1].hist(hypo_up["mean_promoter_beta"], bins=bins, alpha=0.75, color=palette["hypo"], label="Hypo-up")
    axes[1].set_xlabel("Mean promoter beta")
    axes[1].set_ylabel("Gene count")
    axes[1].set_title("Promoter beta distributions")
    axes[1].legend(frameon=True)

    # Show signature sizes and mean expression per signature
    sig_labels = ["Hyper-down", "Hypo-up"]
    sig_counts = [len(hyper_down), len(hypo_up)]
    sig_colors = [palette["hyper"], palette["hypo"]]
    bars = axes[2].bar(sig_labels, sig_counts, color=sig_colors, alpha=0.85)
    for bar, cnt in zip(bars, sig_counts):
        axes[2].text(bar.get_x() + bar.get_width() / 2, cnt + 0.3, str(cnt), ha="center", va="bottom", fontsize=11)
    axes[2].set_ylabel("Gene count")
    axes[2].set_title("Co-occurrence signature sizes")

    fig.suptitle("Arm 2 Discovery: Co-occurrence Transcriptome Signatures", y=1.03)
    fig.tight_layout()
    save_figure(fig, out_dir, "figure3_arm2_discovery")


def plot_arm2_validation(data_dir: Path, out_dir: Path, palette: dict) -> None:
    hd_val = load_tsv(data_dir / "arm2_validation_hyper_down.tsv")
    hu_val = load_tsv(data_dir / "arm2_validation_hypo_up.tsv")
    programs = load_tsv(data_dir / "arm2_validation_metabolic.tsv")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    hd_rate = hd_val.groupby("comparison")["direction_match"].mean() * 100
    hu_rate = hu_val.groupby("comparison")["direction_match"].mean() * 100
    comparisons = sorted(set(hd_rate.index).union(set(hu_rate.index)))
    x = np.arange(len(comparisons))
    width = 0.38

    axes[0].bar(x - width / 2, [hd_rate.get(c, np.nan) for c in comparisons], width, label="Hyper-down", color=palette["hyper"], alpha=0.85)
    axes[0].bar(x + width / 2, [hu_rate.get(c, np.nan) for c in comparisons], width, label="Hypo-up", color=palette["hypo"], alpha=0.85)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(comparisons, rotation=15)
    axes[0].set_ylim(0, 105)
    axes[0].set_ylabel("Direction match (%)")
    axes[0].set_title("Gene-direction replication")
    axes[0].legend(frameon=True)

    hd_plot = hd_val.assign(signature="Hyper-down")
    hu_plot = hu_val.assign(signature="Hypo-up")
    combo = pd.concat([hd_plot, hu_plot], ignore_index=True)
    combo["neglog10_p"] = -np.log10(combo["pvalue"].clip(lower=1e-300))

    axes[1].scatter(
        combo.loc[combo["signature"] == "Hyper-down", "fold_change_mean"],
        combo.loc[combo["signature"] == "Hyper-down", "neglog10_p"],
        s=14,
        alpha=0.45,
        color=palette["hyper"],
        label="Hyper-down",
    )
    axes[1].scatter(
        combo.loc[combo["signature"] == "Hypo-up", "fold_change_mean"],
        combo.loc[combo["signature"] == "Hypo-up", "neglog10_p"],
        s=14,
        alpha=0.45,
        color=palette["hypo"],
        label="Hypo-up",
    )
    axes[1].axvline(1.0, linestyle="--", color="#2c3e50", linewidth=1)
    axes[1].set_xlabel("Case / control mean expression")
    axes[1].set_ylabel(r"$-\log_{10}(p)$")
    axes[1].set_title("Gene-level differential signal")
    axes[1].legend(frameon=True)

    phenotype_order = ["NF", "DCM", "ICM"]
    program_order = [
        p for p in ["oxidative", "glycolytic", "fibrosis", "inflammation", "ecm_remodeling"]
        if p in set(programs["program"])
    ]

    score_pivot = programs.pivot_table(
        index="program",
        columns="phenotype",
        values="mean_program_score",
        aggfunc="mean",
    )
    score_pivot = score_pivot.reindex(index=program_order, columns=phenotype_order)

    im = axes[2].imshow(score_pivot.values, aspect="auto", cmap="YlOrRd")
    axes[2].set_xticks(np.arange(len(phenotype_order)))
    axes[2].set_xticklabels(phenotype_order)
    axes[2].set_yticks(np.arange(len(program_order)))
    axes[2].set_yticklabels(program_order)
    axes[2].set_title("Program activity heatmap")
    plt.colorbar(im, ax=axes[2], fraction=0.046, pad=0.04, label="Mean program score")

    fig.suptitle("Arm 2 Validation: GSE116250 Replication", y=1.03)
    fig.tight_layout()
    save_figure(fig, out_dir, "figure4_arm2_validation")


def plot_step2_pathway_separate(data_dir: Path, out_dir: Path, palette: dict) -> None:
    """Create separate Step 2 pathway figures for fibrosis, inflammation, and ECM remodeling."""
    programs = load_tsv(data_dir / "arm2_validation_metabolic.tsv")
    phenotype_order = ["NF", "DCM", "ICM"]

    targets = [
        ("fibrosis", "Step 2: Fibrosis Program"),
        ("inflammation", "Step 2: Inflammation Program"),
        ("ecm_remodeling", "Step 2: ECM Remodeling Program"),
    ]

    color_map = {
        "fibrosis": palette["fibrosis"],
        "inflammation": palette["inflammation"],
        "ecm_remodeling": palette["ecm_remodeling"],
    }

    for program_name, title in targets:
        subset = programs[programs["program"] == program_name].copy()
        if len(subset) == 0:
            continue

        values = [subset.loc[subset["phenotype"] == ph, "mean_program_score"].mean() for ph in phenotype_order]

        fig, ax = plt.subplots(1, 1, figsize=(6.5, 4.5))
        bars = ax.bar(phenotype_order, values, color=color_map[program_name], alpha=0.9, width=0.62)

        for bar, val in zip(bars, values):
            if np.isfinite(val):
                ax.text(bar.get_x() + bar.get_width() / 2.0, val + 0.03, f"{val:.2f}", ha="center", va="bottom", fontsize=9)

        ax.set_ylabel("Mean program score")
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.25)

        module_source = subset["module_source"].dropna().iloc[0] if "module_source" in subset.columns and subset["module_source"].notna().any() else "unknown"
        ax.text(0.99, 0.02, f"module_source: {module_source}", transform=ax.transAxes, ha="right", va="bottom", fontsize=8)

        fig.tight_layout()
        save_figure(fig, out_dir, f"step2_{program_name}_separate")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate publication-ready validation figures.")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("results/validation/outputs"),
        help="Directory containing validation output TSV files.",
    )
    parser.add_argument(
        "--inputs-dir",
        type=Path,
        default=Path("results/real_processed/inputs"),
        help="Directory containing prepared real input matrices.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/validation/figures"),
        help="Directory where figures will be written.",
    )
    parser.add_argument(
        "--arms",
        choices=["all", "arm1"],
        default="all",
        help="Select whether to render all figures or only Arm 1 figures.",
    )

    args = parser.parse_args()
    palette = configure_style()

    plot_arm1_discovery(args.input_dir, args.inputs_dir, args.output_dir, palette)
    plot_arm1_validation(args.input_dir, args.output_dir, palette)
    plot_arm1_concordance_venn(args.input_dir, args.inputs_dir, args.output_dir, palette)
    if args.arms == "all":
        plot_arm2_discovery(args.input_dir, args.output_dir, palette)
        plot_arm2_validation(args.input_dir, args.output_dir, palette)
        plot_step2_pathway_separate(args.input_dir, args.output_dir, palette)

    print(f"Saved publication-ready figures to: {args.output_dir}")


if __name__ == "__main__":
    main()
