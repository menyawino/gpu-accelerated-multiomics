#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


def read_matrix(path: str | Path, index_col: str = "gene_id") -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if index_col in df.columns:
        df = df.set_index(index_col)
    elif df.columns[0] != index_col:
        df = df.set_index(df.columns[0])
    return df


def adjust_pvalues(pvalues: pd.Series) -> pd.Series:
    if pvalues.empty:
        return pvalues
    vals = pvalues.fillna(1.0).clip(lower=0.0, upper=1.0)
    return pd.Series(multipletests(vals, method="fdr_bh")[1], index=pvalues.index)


def differential_expression(
    expr: pd.DataFrame,
    metadata: pd.DataFrame,
    case_label: str,
    control_label: str,
    min_group_size: int,
) -> pd.DataFrame:
    md = metadata.copy()
    md["sample_id"] = md["sample_id"].astype(str)
    md = md[md["sample_id"].isin(expr.columns)]
    case_samples = md.loc[md["group"] == case_label, "sample_id"].tolist()
    control_samples = md.loc[md["group"] == control_label, "sample_id"].tolist()
    if len(case_samples) < min_group_size or len(control_samples) < min_group_size:
        raise ValueError(
            f"Not enough samples for contrast {case_label} vs {control_label}. "
            f"Need >= {min_group_size} per group."
        )

    log_expr = np.log2(expr[case_samples + control_samples].astype(float) + 1.0)
    case = log_expr[case_samples]
    control = log_expr[control_samples]

    stat, pvals = ttest_ind(case.values, control.values, axis=1, equal_var=False, nan_policy="omit")
    res = pd.DataFrame(index=log_expr.index)
    res["mean_case"] = case.mean(axis=1)
    res["mean_control"] = control.mean(axis=1)
    res["log2fc"] = res["mean_case"] - res["mean_control"]
    res["statistic"] = stat
    res["pvalue"] = pvals
    res["fdr"] = adjust_pvalues(res["pvalue"])
    res.index.name = "gene_id"
    return res.reset_index().sort_values(["fdr", "pvalue", "gene_id"])


def load_gmt(path: str | Path) -> dict[str, set[str]]:
    gmt = {}
    p = Path(path)
    if not p.exists():
        return gmt
    for line in p.read_text().splitlines():
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        gmt[parts[0]] = {g for g in parts[2:] if g}
    return gmt


def enrichment_by_hypergeom(gene_set: set[str], pathways: dict[str, set[str]], universe: set[str]) -> pd.DataFrame:
    from scipy.stats import hypergeom

    rows = []
    if not gene_set or not pathways or not universe:
        return pd.DataFrame(columns=["pathway", "overlap", "pathway_size", "query_size", "pvalue", "fdr"])

    q = len(gene_set)
    n = len(universe)
    for pathway, members in pathways.items():
        members = members.intersection(universe)
        if not members:
            continue
        overlap = len(gene_set.intersection(members))
        if overlap == 0:
            continue
        m = len(members)
        p = hypergeom.sf(overlap - 1, n, m, q)
        rows.append((pathway, overlap, m, q, p))

    out = pd.DataFrame(rows, columns=["pathway", "overlap", "pathway_size", "query_size", "pvalue"])
    if out.empty:
        out["fdr"] = []
        return out
    out["fdr"] = adjust_pvalues(out["pvalue"])
    return out.sort_values(["fdr", "pvalue", "pathway"]).reset_index(drop=True)
