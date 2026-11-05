"""Packaged for revise/; source: project/review/R1.4/scripts/12_prevalence_filter_rebuild_Fig2D.py
Figures: Fig2E
Note: prevalence-filtered rebuild
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
"""Rebuild Fig.2D no-PI KO–metabolite network AFTER prevalence filtering.

Filter: keep features with nonzero abundance in ≥ MIN_NONZERO samples (default 10/20 = 50%).
Then: Spearman |ρ|≥0.80, P<0.05, FDR≤0.05; bootstrap B=500 + LOO; export tables for
fr_tight network + absRho-vs-bootstrap figure.
"""
from __future__ import annotations

from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = PKG_ROOT
OUT = ROOT / "review/R1.4/Fig2D_KOmetab_noPI"
TAB = OUT / "tables"
FIG = OUT / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SEED = 20261021
N_BOOT = 500
RHO_MIN = 0.80
P_MAX = 0.05
FDR_MAX = 0.05
MIN_NONZERO = 10  # ≥50% of n=20
RHO_DISPLAY = 0.85

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
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def dual_write(df: pd.DataFrame, stem: str) -> None:
    df.to_csv(TAB / f"{stem}.tsv", sep="\t", index=False)
    df.to_csv(TAB / f"{stem}.csv", index=False)
    print(f"  wrote {stem} ({len(df)} rows)")


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = p.size
    out = np.full(n, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    q = ranked * n / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[order] = np.clip(q, 0, 1)
    return out


def load_matrices():
    group = pd.read_csv(
        ROOT / "01this_study_dataset/sgb/metab_corr/this.study.group.tsv",
        sep="\t",
        header=None,
        names=["sample", "group"],
    )
    samples = group["sample"].tolist()

    diff_ko = pd.read_csv(
        ROOT / "01this_study_dataset/M2-humann/v2/Figure3/tables/Table_differential_KO_full.tsv",
        sep="\t",
    )
    diff_ko = diff_ko[diff_ko["Robust_DA"].astype(str).isin(["True", "TRUE", "1", "true"])].copy()
    if "function_name" in diff_ko.columns:
        fn = diff_ko["function_name"].astype(str)
        diff_ko = diff_ko[~fn.str.contains(r"unclassified|NO_NAME", case=False, na=False)]
    ko_ids = diff_ko["ko_id"].astype(str).tolist()
    ko_lab = dict(zip(diff_ko["ko_id"].astype(str), diff_ko.get("function_name", diff_ko["ko_id"]).astype(str)))
    ko_dir = {}
    if "Robust_direction" in diff_ko.columns:
        ko_dir = dict(zip(diff_ko["ko_id"].astype(str), diff_ko["Robust_direction"].astype(str)))

    metab_rob = pd.read_csv(
        ROOT
        / "01this_study_dataset/sgb/metab_corr/paper_result/05_RobustDA_metab_DA/tables/metabolite_RobustDA_only_with_HMDB_Class.tsv",
        sep="\t",
    )
    metab_ids = metab_rob["metab_id"].astype(str).tolist()
    metab_lab = {}
    for _, r in metab_rob.iterrows():
        mid = str(r["metab_id"])
        name = r["Metabolite"] if "Metabolite" in metab_rob.columns else mid
        if pd.isna(name) or str(name) in ("-", ""):
            name = mid
        metab_lab[mid] = str(name)
    metab_dir = {}
    if "Robust_direction" in metab_rob.columns:
        metab_dir = dict(zip(metab_rob["metab_id"].astype(str), metab_rob["Robust_direction"].astype(str)))

    ko_path = ROOT / "01this_study_dataset/M2-humann/comb_genefamilies-cpm.uniref90_ko.normed.named.tsv"
    header = pd.read_csv(ko_path, sep="\t", nrows=0).columns.tolist()
    feat_col = header[0]
    colmap = {}
    for s in samples:
        hits = [h for h in header if h.startswith(f"{s}_") or h == s]
        if hits:
            colmap[s] = hits[0]
    usecols = [feat_col] + list(colmap.values())
    chunks = []
    want = set(ko_ids)
    for chunk in pd.read_csv(ko_path, sep="\t", usecols=usecols, chunksize=8000):
        feat = chunk[feat_col].astype(str)
        chunk = chunk[~feat.str.contains(r"\|", regex=True)].copy()
        chunk["ko_id"] = chunk[feat_col].astype(str).str.extract(r"(K\d{5})", expand=False)
        chunk = chunk[chunk["ko_id"].isin(want)]
        if len(chunk):
            chunks.append(chunk)
    ko = pd.concat(chunks, ignore_index=True)
    ko = ko.rename(columns={v: k for k, v in colmap.items()})
    ko = ko.groupby("ko_id")[samples].sum()
    ko = ko.loc[[k for k in ko_ids if k in ko.index]]

    metab = pd.read_csv(
        ROOT / "01this_study_dataset/sgb/metab_corr/output/tables/metab_abund_AB_select.tsv",
        sep="\t",
        index_col=0,
    )
    metab = metab[[c for c in samples if c in metab.columns]]
    metab = metab.loc[[m for m in metab_ids if m in metab.index]]

    return samples, ko, metab, ko_lab, metab_lab, ko_dir, metab_dir


def prevalence_filter(df: pd.DataFrame, samples, min_nonzero: int, kind: str):
    x = df[samples].astype(float)
    n_nz = (x > 0).sum(axis=1)
    frac = n_nz / len(samples)
    keep = n_nz >= min_nonzero
    report = pd.DataFrame(
        {
            "id": df.index.astype(str),
            "type": kind,
            "n_nonzero": n_nz.values,
            "frac_nonzero": frac.values,
            "pass_prevalence": keep.values.astype(int),
            "min_nonzero_threshold": min_nonzero,
            "n_samples": len(samples),
        }
    )
    return df.loc[keep].copy(), report


def cross_spearman_edges(ko_mat, metab_mat, ko_names, metab_names):
    n_s = ko_mat.shape[1]
    rk = np.apply_along_axis(stats.rankdata, 1, ko_mat)
    rm = np.apply_along_axis(stats.rankdata, 1, metab_mat)
    rk = (rk - rk.mean(axis=1, keepdims=True)) / (rk.std(axis=1, keepdims=True) + 1e-12)
    rm = (rm - rm.mean(axis=1, keepdims=True)) / (rm.std(axis=1, keepdims=True) + 1e-12)
    R = (rk @ rm.T) / (n_s - 1)
    dof = n_s - 2
    with np.errstate(divide="ignore", invalid="ignore"):
        t = R * np.sqrt(dof / np.clip(1 - R**2, 1e-12, None))
        P = 2 * stats.t.sf(np.abs(t), dof)
    ii, jj = np.where(np.isfinite(R))
    rho = R[ii, jj]
    p = P[ii, jj]
    fdr = bh_fdr(p)
    keep = (np.abs(rho) >= RHO_MIN) & (p < P_MAX) & (fdr <= FDR_MAX)
    edges = []
    for i, j, r, pv, q, k in zip(ii, jj, rho, p, fdr, keep):
        if k:
            edges.append((ko_names[i], metab_names[j], float(r), float(pv), float(q)))
    return edges


def edge_key(a, b):
    return (a, b)


def main():
    print("=== Prevalence-filtered Fig.2D rebuild ===")
    print(f"  MIN_NONZERO={MIN_NONZERO}/{20} (≥{MIN_NONZERO/20:.0%})")
    samples, ko0, metab0, ko_lab, metab_lab, ko_dir, metab_dir = load_matrices()
    print(f"  before prevalence: KO {ko0.shape}, metab {metab0.shape}")

    ko, ko_rep = prevalence_filter(ko0, samples, MIN_NONZERO, "KO")
    metab, metab_rep = prevalence_filter(metab0, samples, MIN_NONZERO, "Metabolite")
    # also drop residual zero-variance (should be rare after prevalence)
    ko = ko.loc[ko[samples].std(axis=1) > 0]
    metab = metab.loc[metab[samples].std(axis=1) > 0]
    prev_rep = pd.concat([ko_rep, metab_rep], ignore_index=True)
    dual_write(prev_rep, "Table_Fig2D_noPI_prevalence_filter_report")
    summ_prev = pd.DataFrame(
        [
            {
                "type": "KO",
                "n_before": int(len(ko_rep)),
                "n_after": int(ko_rep["pass_prevalence"].sum()),
                "n_removed": int((ko_rep["pass_prevalence"] == 0).sum()),
                "min_nonzero": MIN_NONZERO,
            },
            {
                "type": "Metabolite",
                "n_before": int(len(metab_rep)),
                "n_after": int(metab_rep["pass_prevalence"].sum()),
                "n_removed": int((metab_rep["pass_prevalence"] == 0).sum()),
                "min_nonzero": MIN_NONZERO,
            },
        ]
    )
    dual_write(summ_prev, "Table_Fig2D_noPI_prevalence_filter_summary")
    print(f"  after prevalence: KO {ko.shape}, metab {metab.shape}")
    print(summ_prev.to_string(index=False))

    ko_names = ko.index.tolist()
    metab_names = metab.index.tolist()
    Xk = ko[samples].values.astype(float)
    Xm = metab[samples].values.astype(float)

    obs_edges = cross_spearman_edges(Xk, Xm, ko_names, metab_names)
    print(f"  observed edges |ρ|≥{RHO_MIN}: {len(obs_edges)}")
    obs_df = pd.DataFrame(obs_edges, columns=["KO", "metabolite", "rho", "p", "FDR"])
    dual_write(obs_df, "Table_Fig2D_noPI_observed_edges_rebuilt")

    # Enriched edge table (replace Module13-filtered noPI table used downstream)
    rows_e = []
    for a, b, r, p, q in obs_edges:
        rows_e.append(
            {
                "from": a,
                "to": b,
                "rho": r,
                "p": p,
                "type_from": "KO",
                "type_to": "Metabolite",
                "edge_type": "KO–Metabolite",
                "FDR": q,
                "edge_sign": "positive" if r > 0 else "negative",
                "involves_PI": False,
                "enrich_from": ko_dir.get(a, ""),
                "label_from": a,
                "type_from2": "KO",
                "enrich_to": metab_dir.get(b, ""),
                "label_to": metab_lab.get(b, b),
                "type_to2": "Metabolite",
                "same_direction": ko_dir.get(a, "") == metab_dir.get(b, "") if ko_dir.get(a) and metab_dir.get(b) else "",
                "prevalence_filtered": True,
                "min_nonzero": MIN_NONZERO,
            }
        )
    edges_tab = pd.DataFrame(rows_e)
    dual_write(edges_tab, "Table_Fig2D_KOmetab_edges_noPI")

    # Bootstrap + LOO
    print(f"  bootstrap B={N_BOOT} ...")
    rng = np.random.default_rng(SEED)
    edge_support = Counter()
    n = len(samples)
    boot_n_edges = []
    for b in range(N_BOOT):
        idx = rng.choice(n, size=n, replace=True)
        edges = cross_spearman_edges(Xk[:, idx], Xm[:, idx], ko_names, metab_names)
        boot_n_edges.append(len(edges))
        for a, c, *_ in edges:
            edge_support[edge_key(a, c)] += 1
        if (b + 1) % 50 == 0:
            print(f"    bootstrap {b+1}/{N_BOOT}")

    print("  leave-one-out ...")
    loo_support = Counter()
    for i in range(n):
        keep = [j for j in range(n) if j != i]
        edges = cross_spearman_edges(Xk[:, keep], Xm[:, keep], ko_names, metab_names)
        for a, c, *_ in edges:
            loo_support[edge_key(a, c)] += 1

    support_rows = []
    for a, b, r, p, q in obs_edges:
        k = edge_key(a, b)
        bs = edge_support[k] / N_BOOT
        lo = loo_support[k] / float(n)
        support_rows.append(
            {
                "KO": a,
                "metabolite": b,
                "rho_obs": r,
                "FDR_obs": q,
                "bootstrap_support": bs,
                "loo_support": lo,
                "robust_edge": int(bs >= 0.8 and lo >= 0.8),
            }
        )
    edf = pd.DataFrame(support_rows).sort_values(
        ["robust_edge", "bootstrap_support"], ascending=[False, False]
    )
    dual_write(edf, "Table_Fig2D_noPI_edge_support")
    dual_write(edf, "Table_Fig2D_noPI_edge_support_module13edges")

    n_rob = int(edf["robust_edge"].sum())
    print(f"  robust edges: {n_rob}/{len(edf)} ({n_rob/len(edf):.1%})")

    dual_write(
        pd.DataFrame(
            [
                {
                    "network": "Fig2D_KOmetab_noPI_prevalenceFiltered",
                    "min_nonzero": MIN_NONZERO,
                    "n_KO": len(ko_names),
                    "n_metabolite": len(metab_names),
                    "n_observed_edges": len(edf),
                    "n_robust_edges_boot0.8_loo0.8": n_rob,
                    "robust_fraction": n_rob / len(edf) if len(edf) else 0,
                    "median_bootstrap_support": float(edf["bootstrap_support"].median()),
                    "median_loo_support": float(edf["loo_support"].median()),
                    "boot_n_edges_median": float(np.median(boot_n_edges)),
                    "boot_n_edges_ci2.5": float(np.percentile(boot_n_edges, 2.5)),
                    "boot_n_edges_ci97.5": float(np.percentile(boot_n_edges, 97.5)),
                }
            ]
        ),
        "Table_Fig2D_noPI_stability_summary",
    )

    # Display edges |ρ|≥0.85 with support merge
    display = edges_tab.copy()
    display["key"] = list(zip(display["from"], display["to"]))
    edf2 = edf.copy()
    edf2["key"] = list(zip(edf2["KO"], edf2["metabolite"]))
    display = display.merge(
        edf2[["key", "bootstrap_support", "loo_support", "robust_edge"]], on="key", how="left"
    )
    display = display[display["rho"].abs() >= RHO_DISPLAY].copy()
    dual_write(display, "Fig2D_KOmetab_network_display_edges")
    print(f"  display |ρ|≥{RHO_DISPLAY}: {len(display)} edges")

    # absRho vs bootstrap tables + figure
    edf = edf.copy()
    edf["abs_rho"] = edf["rho_obs"].abs()
    bins = [0.80, 0.85, 0.90, 0.95, 1.01]
    labels = ["0.80–0.85", "0.85–0.90", "0.90–0.95", "≥0.95"]
    edf["abs_rho_bin"] = pd.cut(edf["abs_rho"], bins=bins, labels=labels, right=False, include_lowest=True)
    # last bin include 1.0
    edf.loc[edf["abs_rho"] >= 0.95, "abs_rho_bin"] = "≥0.95"

    scatter = edf[["KO", "metabolite", "rho_obs", "abs_rho", "bootstrap_support", "loo_support", "robust_edge"]].copy()
    dual_write(scatter, "FigS_Fig2D_noPI_absRho_vs_bootstrap_scatter_plotdata")
    dual_write(scatter, "Table_S_Fig2D_noPI_network_stability_absRho_vs_bootstrap")

    summ = (
        edf.groupby("abs_rho_bin", observed=True)
        .agg(
            bin_mid=("abs_rho", "mean"),
            n_edges=("abs_rho", "size"),
            median_bootstrap_support=("bootstrap_support", "median"),
            q25_bootstrap_support=("bootstrap_support", lambda s: float(np.percentile(s, 25))),
            q75_bootstrap_support=("bootstrap_support", lambda s: float(np.percentile(s, 75))),
            frac_support_ge_0_8=("bootstrap_support", lambda s: float((s >= 0.8).mean())),
            frac_robust_edge=("robust_edge", "mean"),
        )
        .reset_index()
        .rename(columns={"frac_support_ge_0_8": "frac_support_ge_0.8"})
    )
    # ensure bin order
    summ["abs_rho_bin"] = pd.Categorical(summ["abs_rho_bin"], categories=labels, ordered=True)
    summ = summ.sort_values("abs_rho_bin")
    dual_write(summ, "FigS_Fig2D_noPI_absRho_vs_bootstrap_binned_plotdata")
    dual_write(summ, "Table_S_Fig2D_noPI_stability_by_absRho")

    stem = "FigS_Fig2D_noPI_absRho_vs_bootstrap_support"
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    x = np.arange(len(summ))
    ax.bar(x, summ["frac_robust_edge"], width=0.55, color="#E07A3D", alpha=0.85, label="Proportion robust")
    ax.set_ylabel("Proportion of robust edges")
    ax.set_ylim(0, 1.08)
    ax.set_xticks(x)
    ax.set_xticklabels(summ["abs_rho_bin"], fontsize=8)
    ax.set_xlabel(r"$|\rho|$ bin")
    ax.spines["right"].set_visible(False)

    ax2 = ax.twinx()
    ax2.plot(x, summ["median_bootstrap_support"], "o-", color="#2F4B7C", lw=1.4, ms=5, label="Median bootstrap support")
    ax2.fill_between(
        x,
        summ["q25_bootstrap_support"],
        summ["q75_bootstrap_support"],
        color="#2F4B7C",
        alpha=0.15,
    )
    ax2.set_ylabel("Bootstrap support")
    ax2.set_ylim(0, 1.08)
    ax2.axhline(0.8, color="#C0392B", ls="--", lw=0.9)
    ax2.spines["top"].set_visible(False)

    for i, n_e in enumerate(summ["n_edges"]):
        ax.text(i, min(0.98, float(summ["frac_robust_edge"].iloc[i]) + 0.04), f"n={n_e}", ha="center", fontsize=7)

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=7, loc="upper left")
    ax.set_title("Stability by |ρ| bin\n(prevalence-filtered KO–metabolite network, no PI)", fontsize=10)
    fig.suptitle(
        f"Stability of Fig. 2D KO–metabolite edges increases with |ρ| (n = {len(edf)}; prevalence ≥{MIN_NONZERO}/{n})",
        fontsize=10.5,
        y=1.02,
    )
    fig.tight_layout()
    for ext in ("pdf", "svg", "png", "jpg"):
        fig.savefig(FIG / f"{stem}.{ext}")
    plt.close(fig)

    (TAB / f"{stem}_CAPTION.txt").write_text(
        f"""Supplementary Figure. Bootstrap stability of Plaque Index–removed Fig. 2D KO–metabolite edges after prevalence filtering.

Features were restricted to those detected (abundance > 0) in ≥{MIN_NONZERO} of {n} samples before Spearman correlation. Edges (n = {len(edf)}; |ρ|≥0.80, P<0.05, FDR≤0.05) are stratified by |ρ| bins. Bars: proportion of robust edges (bootstrap support ≥0.8 and leave-one-out support ≥0.8). Blue line/points: median bootstrap support with interquartile range (shaded). Numbers above bars indicate edge counts per bin. Dashed line: support = 0.8. Bootstrap: participant-level resampling (B={N_BOOT}).
""",
        encoding="utf-8",
    )
    print(f"  fig {stem}")

    # note
    (OUT / "Fig2D_noPI_prevalence_filter_NOTE.md").write_text(
        f"""# Prevalence filter applied before Fig.2D network

- Rule: nonzero abundance in ≥ **{MIN_NONZERO}/{n}** samples (≥{MIN_NONZERO/n:.0%})
- Applied to Robust_DA KO and Robust_DA metabolites (then SD>0)
- KO: {int(summ_prev.loc[summ_prev.type=='KO','n_before'].iloc[0])} → {len(ko_names)}
- Metabolite: {int(summ_prev.loc[summ_prev.type=='Metabolite','n_before'].iloc[0])} → {len(metab_names)}
- Edges |ρ|≥0.80: {len(edf)}; robust: {n_rob} ({n_rob/len(edf):.1%})
- Display |ρ|≥0.85: {len(display)}
""",
        encoding="utf-8",
    )
    print("DONE prevalence rebuild")


if __name__ == "__main__":
    main()
