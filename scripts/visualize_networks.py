#!/usr/bin/env python3
"""Visualize RNA and methylation regulatory networks from edge tables."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


def load_edges(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["source", "target", "correlation", "pvalue", "fdr", "n_samples"])
    return pd.read_csv(path, sep="\t")


def build_graph(edges: pd.DataFrame) -> nx.Graph:
    g = nx.Graph()
    for _, r in edges.iterrows():
        s = str(r["source"])
        t = str(r["target"])
        w = float(r["correlation"])
        g.add_edge(s, t, weight=w)
    return g


def save_hubs_table(g: nx.Graph, out_path: Path, top_n: int = 30) -> None:
    if g.number_of_nodes() == 0:
        pd.DataFrame(columns=["gene", "degree", "weighted_degree", "betweenness"]).to_csv(out_path, sep="\t", index=False)
        return

    deg = dict(g.degree())
    wdeg = dict(g.degree(weight="weight"))
    btw = nx.betweenness_centrality(g, weight=None)
    rows = []
    for node in g.nodes():
        rows.append(
            {
                "gene": node,
                "degree": int(deg[node]),
                "weighted_degree": float(wdeg[node]),
                "betweenness": float(btw[node]),
            }
        )
    hubs = pd.DataFrame(rows).sort_values(["degree", "betweenness"], ascending=[False, False]).head(top_n)
    hubs.to_csv(out_path, sep="\t", index=False)


def plot_single_network(edges: pd.DataFrame, title: str, out_png: Path, out_pdf: Path, max_edges_draw: int = 250) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    if edges.empty:
        axes[0].text(0.5, 0.5, "No edges", ha="center", va="center")
        axes[1].text(0.5, 0.5, "No edges", ha="center", va="center")
        for ax in axes:
            ax.set_axis_off()
        fig.suptitle(title)
        fig.tight_layout()
        fig.savefig(out_png, dpi=400, bbox_inches="tight")
        fig.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)
        return

    # Panel 1: effect distributions.
    axes[0].hist(edges["correlation"], bins=25, color="#1f78b4", alpha=0.8)
    axes[0].axvline(0, linestyle="--", linewidth=1, color="#2c3e50")
    axes[0].set_xlabel("Correlation")
    axes[0].set_ylabel("Edge count")
    axes[0].set_title("Edge correlation distribution")

    # Panel 2: spring-layout graph for strongest edges.
    draw = edges.reindex(edges["correlation"].abs().sort_values(ascending=False).index).head(max_edges_draw)
    g = build_graph(draw)
    pos = nx.spring_layout(g, seed=42, k=0.35)

    weights = np.array([g[u][v]["weight"] for u, v in g.edges()])
    edge_colors = ["#c0392b" if w > 0 else "#1f78b4" for w in weights]
    edge_widths = 0.6 + 2.2 * (np.abs(weights) - np.min(np.abs(weights))) / (np.ptp(np.abs(weights)) + 1e-9)

    degrees = dict(g.degree())
    node_sizes = [20 + 35 * degrees[n] for n in g.nodes()]

    nx.draw_networkx_nodes(g, pos, node_size=node_sizes, node_color="#95a5a6", alpha=0.9, ax=axes[1], linewidths=0)
    nx.draw_networkx_edges(g, pos, width=edge_widths, edge_color=edge_colors, alpha=0.75, ax=axes[1])

    # Label top hubs only.
    top_nodes = sorted(degrees, key=degrees.get, reverse=True)[:12]
    labels = {n: n for n in top_nodes}
    nx.draw_networkx_labels(g, pos, labels=labels, font_size=7, ax=axes[1])

    axes[1].set_title(f"Top {len(draw)} edges network view")
    axes[1].set_axis_off()

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=400, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description="Visualize regulatory networks")
    ap.add_argument("--rna-edges", required=True, type=Path)
    ap.add_argument("--meth-edges", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    rna = load_edges(args.rna_edges)
    meth = load_edges(args.meth_edges)

    plot_single_network(
        rna,
        title="RNA Regulatory Network",
        out_png=args.out_dir / "rna_network_overview.png",
        out_pdf=args.out_dir / "rna_network_overview.pdf",
        max_edges_draw=300,
    )
    plot_single_network(
        meth,
        title="Methylation Regulatory Network",
        out_png=args.out_dir / "methylation_network_overview.png",
        out_pdf=args.out_dir / "methylation_network_overview.pdf",
        max_edges_draw=120,
    )

    save_hubs_table(build_graph(rna), args.out_dir / "rna_top_hubs.tsv", top_n=30)
    save_hubs_table(build_graph(meth), args.out_dir / "methylation_top_hubs.tsv", top_n=30)


if __name__ == "__main__":
    main()
