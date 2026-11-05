"""Packaged for revise/; source: project/review/R2.5/scripts/01_R2_5_edge_FDR_CLR.py
Figures: Methods
Note: edge FDR/CLR sensitivity
"""
from pathlib import Path
import os, sys
PKG_ROOT = Path(__file__).resolve().parents[1]
SHARED = PKG_ROOT / 'shared_input'
OUTDIR = PKG_ROOT / 'output' / Path(__file__).stem
OUTDIR.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(Path(__file__).resolve().parent))
try:
    from paths_config import PATHS, SHARED as _SHARED, PKG_ROOT as _PKG
    SHARED, PKG_ROOT = _SHARED, _PKG
except Exception:
    pass
# --- original script below; prefer SHARED/PATHS over hard-coded host paths ---

#!/usr/bin/env python3
"""R2.5: edge-level BH-FDR + CLR-Spearman sensitivity for MAG co-occurrence networks.

Canonical analysis edges: |r|>=0.5, P<0.05 (no FDR) — matches manuscript/Gephi.
Also report: (i) BH-FDR on pairwise Spearman P; (ii) CLR + Spearman same thresholds.
seed = 20261028
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

SEED = 20261028
np.random.seed(SEED)

ROOT = PKG_ROOT
OUT = Path(__file__).resolve().parents[1] / "tables"
OUT.mkdir(parents=True, exist_ok=True)

ABUND = ROOT / "01this_study_dataset/input/bin_abundance_table.tsv"
GROUP = ROOT / "01this_study_dataset/sgb/species_network/this.study.group.tsv"
GEPHI = {
    "Normal": ROOT
    / "01this_study_dataset/sgb/species_network/result_bins/Normal/Normal_Gephi_edge.csv",
    "PoorPlaque": ROOT
    / "01this_study_dataset/sgb/species_network/result_bins/PoorPlaque/PoorPlaque_Gephi_edge.csv",
}


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        val = ranked[i] * n / (i + 1)
        prev = min(prev, val)
        adj[i] = prev
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adj, 0, 1)
    return out


def clr_transform(X: np.ndarray, pseudo: float = 1e-6) -> np.ndarray:
    """Multiplicative replacement for zeros then CLR (samples x features)."""
    X = np.asarray(X, dtype=float)
    X = X.copy()
    X[X <= 0] = pseudo
    X = X / X.sum(axis=1, keepdims=True)
    logX = np.log(X)
    return logX - logX.mean(axis=1, keepdims=True)


def pairwise_spearman(mat: np.ndarray, labels: list[str]) -> pd.DataFrame:
    """mat: samples x features. Return upper-triangle edges."""
    n_feat = mat.shape[1]
    rows = []
    for i in range(n_feat):
        xi = mat[:, i]
        for j in range(i + 1, n_feat):
            xj = mat[:, j]
            # need variation
            if np.nanstd(xi) == 0 or np.nanstd(xj) == 0:
                continue
            r, p = stats.spearmanr(xi, xj)
            if np.isnan(r) or np.isnan(p):
                continue
            rows.append(
                {
                    "source": labels[i],
                    "target": labels[j],
                    "r": float(r),
                    "p": float(p),
                }
            )
    df = pd.DataFrame(rows)
    if len(df) == 0:
        return df
    df["FDR"] = bh_fdr(df["p"].values)
    return df


def summarize_edges(df: pd.DataFrame, group: str, pipeline: str) -> dict:
    n_nodes = len(set(df["source"]).union(df["target"])) if len(df) else 0
    n_edges = len(df)
    dens = (
        2 * n_edges / (n_nodes * (n_nodes - 1)) if n_nodes > 1 else np.nan
    )
    mean_deg = 2 * n_edges / n_nodes if n_nodes else np.nan
    return {
        "group": group,
        "pipeline": pipeline,
        "n_edges": n_edges,
        "n_pos": int((df["r"] > 0).sum()) if n_edges else 0,
        "n_neg": int((df["r"] < 0).sum()) if n_edges else 0,
        "n_nodes": n_nodes,
        "density": dens,
        "mean_degree": mean_deg,
    }


def main():
    abund = pd.read_csv(ABUND, sep="\t", index_col=0)
    # columns = samples; index = MAG
    group = pd.read_csv(GROUP, sep="\t", header=None, names=["SampleID", "Group"])
    gmap = dict(zip(group["SampleID"], group["Group"]))

    summary_rows = []
    display_rows = []

    for grp in ["Normal", "PoorPlaque"]:
        samples = [s for s in abund.columns if gmap.get(s) == grp]
        # features x samples -> transpose to samples x features
        sub = abund[samples].T
        # drop features with all zeros in group
        keep = (sub > 0).sum(axis=0) >= 1
        sub = sub.loc[:, keep]
        labels = list(sub.columns)

        # P0: relative abundance Spearman
        edges0 = pairwise_spearman(sub.values, labels)
        e_p05 = edges0[(edges0["r"].abs() >= 0.5) & (edges0["p"] < 0.05)].copy()
        e_fdr = edges0[(edges0["r"].abs() >= 0.5) & (edges0["FDR"] < 0.05)].copy()
        summary_rows.append(summarize_edges(e_p05, grp, "P0_rel_Spearman_r0.5_p0.05"))
        summary_rows.append(summarize_edges(e_fdr, grp, "P0_rel_Spearman_r0.5_FDR0.05"))

        e_p05.to_csv(OUT / f"Table_R2.5_edges_P0_{grp}.tsv", sep="\t", index=False)
        e_p05.to_csv(OUT / f"Table_R2.5_edges_P0_{grp}.csv", index=False)
        e_fdr.to_csv(OUT / f"Table_R2.5_edges_P0_FDR_{grp}.tsv", sep="\t", index=False)
        e_fdr.to_csv(OUT / f"Table_R2.5_edges_P0_FDR_{grp}.csv", index=False)

        # P1: CLR Spearman
        clr = clr_transform(sub.values)
        edges1 = pairwise_spearman(clr, labels)
        e1_p05 = edges1[(edges1["r"].abs() >= 0.5) & (edges1["p"] < 0.05)].copy()
        e1_fdr = edges1[(edges1["r"].abs() >= 0.5) & (edges1["FDR"] < 0.05)].copy()
        summary_rows.append(summarize_edges(e1_p05, grp, "P1_CLR_Spearman_r0.5_p0.05"))
        summary_rows.append(summarize_edges(e1_fdr, grp, "P1_CLR_Spearman_r0.5_FDR0.05"))
        e1_p05.to_csv(OUT / f"Table_R2.5_edges_P1_CLR_{grp}.tsv", sep="\t", index=False)
        e1_p05.to_csv(OUT / f"Table_R2.5_edges_P1_CLR_{grp}.csv", index=False)

        # Gephi + display median filter reconciliation
        ge = pd.read_csv(GEPHI[grp])
        v = ge["correlation"].astype(float).values
        thr = float(np.quantile(np.abs(v), 0.5))
        n_disp = int((np.abs(v) >= thr).sum())
        display_rows.append(
            {
                "group": grp,
                "source": "Gephi_analysis_edge_list",
                "n_edges": len(ge),
                "n_pos": int((v > 0).sum()),
                "n_neg": int((v < 0).sum()),
                "n_nodes_in_edges": len(set(ge["source"]).union(ge["target"])),
                "note": "canonical manuscript network",
            }
        )
        display_rows.append(
            {
                "group": grp,
                "source": "Fig5A_display_median_abs_r_filter",
                "n_edges": n_disp,
                "n_pos": int(((v > 0) & (np.abs(v) >= thr)).sum()),
                "n_neg": int(((v < 0) & (np.abs(v) >= thr)).sum()),
                "n_nodes_in_edges": np.nan,
                "note": f"|r| >= median(|r|) = {thr:.6f}; visualization only",
            }
        )
        display_rows.append(
            {
                "group": grp,
                "source": "manuscript_text",
                "n_edges": 771 if grp == "Normal" else 1426,
                "n_pos": 571 if grp == "Normal" else 850,
                "n_neg": 200 if grp == "Normal" else 576,
                "n_nodes_in_edges": 123 if grp == "Normal" else 125,
                "note": "main-text reported counts",
            }
        )

    summ = pd.DataFrame(summary_rows)
    summ.to_csv(OUT / "Table_R2.5_pipeline_FDR_CLR_summary.tsv", sep="\t", index=False)
    summ.to_csv(OUT / "Table_R2.5_pipeline_FDR_CLR_summary.csv", index=False)

    disp = pd.DataFrame(display_rows)
    disp.to_csv(OUT / "Table_R2.5_edge_count_reconciliation.tsv", sep="\t", index=False)
    disp.to_csv(OUT / "Table_R2.5_edge_count_reconciliation.csv", index=False)

    # copy pointer note for R1.4 permutation/bootstrap reuse
    note = OUT / "README_reuse_R1.4.txt"
    note.write_text(
        "Participant-level bootstrap/LOO and label-permutation results are in "
        "../R1.4/tables/Table_B3_* and Table_B4_* (seed=20261021).\n"
        "This folder adds FDR/CLR sensitivity and explicit Fig5A median-|r| reconciliation.\n"
    )
    print(summ.to_string(index=False))
    print("\nDisplay reconciliation:")
    print(disp.to_string(index=False))
    print("Wrote", OUT)


if __name__ == "__main__":
    main()
