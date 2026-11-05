"""Packaged for revise/; source: project/review/R1.4/Fig5_MAG_network_stability/scripts/03_Fig5_stability_prev50_absR.py
Figures: Fig5A;FigS4
Note: MAG network stability prev50
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
"""Rebuild Fig5 |r|-binned bootstrap stability on the SAME input as
Fig_Fig5_rho080_prev50_AB_CNS:
  - MAG abundance: bin_abundance_table.tsv
  - within-group prevalence ≥ 50% (ceil(0.5*n) nonzero)
  - Spearman |r| ≥ 0.5, P < 0.05 primary screen
  - then stratify observed edges by |r| bin (Fig.2D-style dual-axis plot)

Bootstrap B=500 + LOO; robust = boot≥0.8 AND LOO≥0.8.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = PKG_ROOT
OUT = ROOT / "review/R1.4/Fig5_MAG_network_stability"
TAB = OUT / "tables"
FIG = OUT / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

MAG_AB = ROOT / "01this_study_dataset/input/bin_abundance_table.tsv"
GROUP = ROOT / "01this_study_dataset/sgb/species_network/this.study.group.tsv"
REF_N = ROOT / "review/R1.4/Fig5_MAG_network_rho080_prev50/tables/Table_Fig5_rho080_prev50_edges_Normal.tsv"
REF_P = ROOT / "review/R1.4/Fig5_MAG_network_rho080_prev50/tables/Table_Fig5_rho080_prev50_edges_PoorPlaque.tsv"

SEED = 20261021
N_BOOT = 500
R_MIN = 0.5
P_MAX = 0.05
PREV_MIN = 0.5
STEM = "Fig_Fig5_absR_vs_bootstrap_binned"
BIN_ORDER = ["0.50–0.60", "0.60–0.70", "0.70–0.80", "0.80–0.90", "≥0.90"]


def dual_write(df: pd.DataFrame, stem: str) -> None:
    df.to_csv(TAB / f"{stem}.tsv", sep="\t", index=False)
    df.to_csv(TAB / f"{stem}.csv", index=False)
    print(f"  wrote {stem} ({len(df)} rows)")


def corr_edges(X: np.ndarray, names: list[str], r_min=R_MIN, p_max=P_MAX):
    n_f, n_s = X.shape
    if n_s < 4 or n_f < 2:
        return []
    ranks = np.apply_along_axis(stats.rankdata, 1, X)
    sd = ranks.std(axis=1)
    valid = sd > 0
    ranks2 = ranks.copy()
    ranks2[~valid] = np.nan
    idx = np.where(valid)[0]
    if len(idx) < 2:
        return []
    R = np.corrcoef(ranks2[idx])
    dof = n_s - 2
    with np.errstate(divide="ignore", invalid="ignore"):
        t = R * np.sqrt(dof / np.clip(1 - R**2, 1e-12, None))
        Pmat = 2 * stats.t.sf(np.abs(t), dof)
    edges = []
    m = len(idx)
    for i in range(m):
        for j in range(i + 1, m):
            r, p = R[i, j], Pmat[i, j]
            if np.isfinite(r) and np.isfinite(p) and abs(r) >= r_min and p < p_max:
                edges.append((names[idx[i]], names[idx[j]], float(r), float(p)))
    return edges


def abs_r_bin(a: float) -> str:
    if a < 0.6:
        return "0.50–0.60"
    if a < 0.7:
        return "0.60–0.70"
    if a < 0.8:
        return "0.70–0.80"
    if a < 0.9:
        return "0.80–0.90"
    return "≥0.90"


def edge_key(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted([a, b]))


def style_one_panel(ax, summ: pd.DataFrame, title: str) -> object:
    x = np.arange(len(summ))
    ax.bar(
        x,
        summ["frac_robust_edge"],
        width=0.55,
        color="#E07A3D",
        alpha=0.85,
        label="Proportion robust",
        zorder=2,
    )
    ax.set_ylabel("Proportion of robust edges", fontsize=11)
    ax.set_ylim(0, 1.28)
    ax.set_xticks(x)
    ax.set_xticklabels(
        summ["abs_r_bin"],
        fontsize=9,
        rotation=45,
        ha="right",
        rotation_mode="anchor",
    )
    ax.tick_params(axis="y", labelsize=9)
    ax.set_xlabel(r"$|r|$ bin", fontsize=11)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax2 = ax.twinx()
    ax2.fill_between(
        x,
        summ["q25_bootstrap_support"],
        summ["q75_bootstrap_support"],
        color="#2F4B7C",
        alpha=0.15,
        zorder=1,
    )
    ax2.plot(
        x,
        summ["median_bootstrap_support"],
        "o-",
        color="#2F4B7C",
        lw=1.5,
        ms=5.5,
        label="Median bootstrap support",
        zorder=3,
    )
    ax2.set_ylabel("Bootstrap support", fontsize=11)
    ax2.set_ylim(0, 1.28)
    ax2.tick_params(axis="y", labelsize=9)
    ax2.axhline(0.8, color="#C0392B", ls="--", lw=0.9, zorder=1)
    ax2.spines["top"].set_visible(False)

    for i, n in enumerate(summ["n_edges"].astype(int)):
        ax.text(i, 1.15, f"n={n}", ha="center", va="center", fontsize=8.5, zorder=10, clip_on=False)
    ax.set_title(title, fontsize=11, pad=6)
    return ax2


def main() -> None:
    print("=== Fig5 absR stability on rho080_prev50 input ===")
    grp = pd.read_csv(GROUP, sep="\t", header=None, names=["Sample", "Group"])
    mag = pd.read_csv(MAG_AB, sep="\t", index_col=0)
    samples = [s for s in grp["Sample"].tolist() if s in mag.columns]
    meta = grp[grp["Sample"].isin(samples)].copy()
    mag = mag[samples]

    kept: dict[str, list[str]] = {}
    obs_r: dict[str, dict] = {}
    obs_edges: dict[str, list] = {}
    recon_rows = []

    for label in ["Normal", "PoorPlaque"]:
        smp = meta.loc[meta["Group"] == label, "Sample"].tolist()
        Xg = mag[smp]
        n_s = len(smp)
        n_nz = (Xg > 0).sum(axis=1)
        keep = n_nz[n_nz >= int(np.ceil(PREV_MIN * n_s))].index.tolist()
        kept[label] = keep
        edges = corr_edges(Xg.loc[keep].values.astype(float), keep)
        obs_edges[label] = edges
        obs_r[label] = {edge_key(a, c): (r, p) for a, c, r, p in edges}
        n_ge08 = sum(1 for a, c, r, p in edges if abs(r) >= 0.8)
        recon_rows.append(
            {
                "group": label,
                "n_MAG_pass_prev50": len(keep),
                "n_edges_r0.5": len(edges),
                "n_edges_r0.8": n_ge08,
            }
        )
        print(f"  {label}: MAGs pass prev50={len(keep)}; |r|≥0.5={len(edges)}; |r|≥0.8={n_ge08}")

    # reconcile with Fig5A edge tables
    for label, ref_path in [("Normal", REF_N), ("PoorPlaque", REF_P)]:
        ref = pd.read_csv(ref_path, sep="\t")
        ref_keys = {edge_key(a, b) for a, b in zip(ref["source"], ref["target"])}
        our_ge08 = {k for k, (r, p) in obs_r[label].items() if abs(r) >= 0.8}
        inter = len(ref_keys & our_ge08)
        print(
            f"  reconcile {label}: Fig5A |r|≥0.8={len(ref_keys)}; "
            f"rebuilt={len(our_ge08)}; intersect={inter}"
        )
        if inter != len(ref_keys) or len(our_ge08) != len(ref_keys):
            raise SystemExit(f"ERROR: |r|≥0.8 edge set mismatch vs {ref_path.name}")

    dual_write(pd.DataFrame(recon_rows), "Table_Fig5_prev50_edgecount_reconcile")
    for label in ["Normal", "PoorPlaque"]:
        dual_write(
            pd.DataFrame(obs_edges[label], columns=["source", "target", "r", "p"]),
            f"Table_Fig5_rebuilt_edges_prev50_{label}",
        )

    # Bootstrap + LOO on prevalence-filtered MAG matrix (fixed MAG set)
    print(f"  bootstrap B={N_BOOT} ...")
    rng = np.random.default_rng(SEED)
    edge_support = {g: {} for g in ["Normal", "PoorPlaque"]}
    loo_support = {g: {} for g in ["Normal", "PoorPlaque"]}
    loo_n = {}

    for label in ["Normal", "PoorPlaque"]:
        smp = meta.loc[meta["Group"] == label, "Sample"].tolist()
        keep = kept[label]
        Xg = mag.loc[keep, smp].values.astype(float)
        n = Xg.shape[1]
        loo_n[label] = n
        for b in range(N_BOOT):
            idx = rng.choice(n, size=n, replace=True)
            for a, c, r, p in corr_edges(Xg[:, idx], keep):
                k = edge_key(a, c)
                edge_support[label][k] = edge_support[label].get(k, 0) + 1
            if (b + 1) % 100 == 0:
                print(f"    {label} bootstrap {b+1}/{N_BOOT}")
        for i in range(n):
            keep_i = [j for j in range(n) if j != i]
            for a, c, r, p in corr_edges(Xg[:, keep_i], keep):
                k = edge_key(a, c)
                loo_support[label][k] = loo_support[label].get(k, 0) + 1

    summary_rows = []
    absr_scatter_rows = []
    binned_all = []
    for label in ["Normal", "PoorPlaque"]:
        rows = []
        for (a, c), (r_obs, p_obs) in obs_r[label].items():
            bs = edge_support[label].get((a, c), 0) / N_BOOT
            lo = loo_support[label].get((a, c), 0) / float(loo_n[label])
            rows.append(
                {
                    "group": label,
                    "source": a,
                    "target": c,
                    "in_observed_network": 1,
                    "r_obs": r_obs,
                    "abs_r_obs": abs(r_obs),
                    "p_obs": p_obs,
                    "bootstrap_support": bs,
                    "loo_support": lo,
                    "robust_edge": int(bs >= 0.8 and lo >= 0.8),
                    "prevalence_filter": PREV_MIN,
                }
            )
        edf = pd.DataFrame(rows)
        dual_write(edf, f"Table_Fig5_edge_support_prev50_{label}")
        n_rob = int(edf["robust_edge"].sum())
        summary_rows.append(
            {
                "group": label,
                "prevalence_filter": PREV_MIN,
                "n_MAG_pass_prev50": len(kept[label]),
                "n_observed_edges": len(edf),
                "n_robust_edges_boot0.8_loo0.8": n_rob,
                "robust_fraction": n_rob / len(edf) if len(edf) else np.nan,
                "median_bootstrap_support": float(edf["bootstrap_support"].median()),
                "median_loo_support": float(edf["loo_support"].median()),
                "frac_boot_ge_0.8": float((edf["bootstrap_support"] >= 0.8).mean()),
            }
        )
        print(f"  {label} robust: {n_rob}/{len(edf)} ({n_rob/len(edf):.1%})")
        for _, r in edf.iterrows():
            absr_scatter_rows.append(
                {
                    "group": label,
                    "source": r["source"],
                    "target": r["target"],
                    "r_obs": r["r_obs"],
                    "abs_r_obs": r["abs_r_obs"],
                    "bootstrap_support": r["bootstrap_support"],
                    "loo_support": r["loo_support"],
                    "robust_edge": r["robust_edge"],
                    "prevalence_filter": PREV_MIN,
                }
            )

        edf = edf.copy()
        edf["abs_r_bin"] = edf["abs_r_obs"].map(abs_r_bin)
        binned = (
            edf.groupby("abs_r_bin", observed=False)
            .agg(
                n_edges=("robust_edge", "size"),
                n_robust=("robust_edge", "sum"),
                median_bootstrap_support=("bootstrap_support", "median"),
                q25_bootstrap_support=("bootstrap_support", lambda s: float(np.percentile(s, 25))),
                q75_bootstrap_support=("bootstrap_support", lambda s: float(np.percentile(s, 75))),
                frac_robust_edge=("robust_edge", "mean"),
                frac_boot_ge_0_8=("bootstrap_support", lambda s: float((s >= 0.8).mean())),
            )
            .reindex(BIN_ORDER)
            .dropna(how="all")
            .reset_index()
        )
        binned.insert(0, "group", label)
        binned["prevalence_filter"] = PREV_MIN
        dual_write(binned, f"Table_Fig5_stability_by_absR_prev50_{label}")
        # also overwrite the primary by_absR tables used by replot
        dual_write(binned.drop(columns=["prevalence_filter"]), f"Table_Fig5_stability_by_absR_{label}")
        binned_all.append(binned)

    comb = pd.concat(binned_all, ignore_index=True)
    dual_write(comb, "Table_Fig5_stability_by_absR")
    dual_write(comb, "Table_Fig5_stability_by_absR_prev50")
    dual_write(pd.DataFrame(summary_rows), "Table_Fig5_stability_summary_prev50")
    dual_write(pd.DataFrame(absr_scatter_rows), "Table_Fig5_absR_vs_bootstrap_scatter_prev50")

    # ---- styled figure (match 02_replot) ----
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif", "Times"],
            "font.weight": "normal",
            "axes.labelweight": "normal",
            "axes.titleweight": "normal",
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "savefig.dpi": 300,
            "axes.spines.top": False,
        }
    )

    def prep(df, group):
        g = df[df["group"] == group].copy()
        return g.set_index("abs_r_bin").reindex(BIN_ORDER).dropna(subset=["n_edges"]).reset_index()

    summN = prep(comb, "Normal")
    summP = prep(comb, "PoorPlaque")

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.6))
    ax2N = style_one_panel(axes[0], summN, "Normal")
    style_one_panel(axes[1], summP, "PoorPlaque")
    h1, l1 = axes[0].get_legend_handles_labels()
    h2, l2 = ax2N.get_legend_handles_labels()
    fig.legend(
        h1 + h2,
        l1 + l2,
        frameon=False,
        fontsize=8,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.995),
        ncol=2,
    )
    fig.subplots_adjust(left=0.09, right=0.91, bottom=0.18, top=0.86, wspace=0.42)

    plot_df = pd.concat(
        [summN.assign(group="Normal"), summP.assign(group="PoorPlaque")],
        ignore_index=True,
    )
    dual_write(plot_df, f"{STEM}_plotdata")

    for ext in ("pdf", "svg", "png", "jpg"):
        fig.savefig(FIG / f"{STEM}.{ext}", dpi=300, bbox_inches=None, pad_inches=0.05)
    plt.close(fig)

    for group, summ in [("Normal", summN), ("PoorPlaque", summP)]:
        stem_g = f"{STEM}_{group}"
        fig, ax = plt.subplots(figsize=(4.2, 5.0))
        ax2 = style_one_panel(ax, summ, group)
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8, loc="upper left")
        fig.subplots_adjust(left=0.16, right=0.84, bottom=0.16, top=0.92)
        dual_write(summ.assign(group=group), f"{stem_g}_plotdata")
        for ext in ("pdf", "svg", "png", "jpg"):
            fig.savefig(FIG / f"{stem_g}.{ext}", dpi=300, bbox_inches=None, pad_inches=0.05)
        plt.close(fig)

    caption = (
        "Figure. Bootstrap stability of MAG co-occurrence edges increases with absolute Spearman |r| "
        "(same input as Fig. 5A: within-group prevalence ≥ 50%).\n\n"
        f"Group-specific MAG–MAG edges after within-group prevalence ≥ 50% "
        f"(Normal n = {int(summN['n_edges'].sum())}; PoorPlaque n = {int(summP['n_edges'].sum())}; "
        "|r| ≥ 0.5, P < 0.05; n = 10 participants per group) were stratified by |r| bins. "
        "Bars: proportion of robust edges (bootstrap support ≥ 0.8 and leave-one-out support ≥ 0.8). "
        "Blue line/points: median bootstrap support with interquartile range (shaded). "
        "Numbers above bars indicate edge counts per bin. Dashed line: support = 0.8. "
        "Bootstrap: participant-level resampling (B = 500). "
        "The |r| ≥ 0.8 subset matches Fig. 5A edge counts (Normal 199; PoorPlaque 428). "
        "Co-occurrence denotes statistical association, not direct microbial interaction.\n"
    )
    (TAB / f"{STEM}_CAPTION.txt").write_text(caption, encoding="utf-8")
    (FIG / f"{STEM}_CAPTION.txt").write_text(caption, encoding="utf-8")
    print(f"saved {STEM}")


if __name__ == "__main__":
    main()
