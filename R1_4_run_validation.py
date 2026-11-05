"""Packaged for revise/; source: project/review/R1.4/scripts/01_run_R1_4_validation.py
Figures: shared
Note: R1.4 validation runner
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
"""R1.4 validation: PI circularity + MAG network robustness + cross-omics partial corr."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module

P = import_module("00_paths")

# ---- plotting ----
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif", "Times"],
        "font.weight": "normal",
        "axes.labelweight": "normal",
        "axes.titleweight": "normal",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def dual_write(df: pd.DataFrame, stem: str) -> None:
    tsv = P.TABLES / f"{stem}.tsv"
    csv = P.TABLES / f"{stem}.csv"
    df.to_csv(tsv, sep="\t", index=False)
    df.to_csv(csv, index=False)
    print(f"  wrote {stem} ({len(df)} rows)")


def save_fig(fig: plt.Figure, stem: str, plot_df: pd.DataFrame | None = None) -> None:
    if plot_df is not None:
        dual_write(plot_df, f"{stem}_plotdata")
    for ext in ("pdf", "svg", "png", "jpg"):
        fig.savefig(P.FIGS / f"{stem}.{ext}")
    plt.close(fig)
    print(f"  fig {stem}")


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    out = np.full(n, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    q = ranked * n / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out[order] = q
    return out


def spearman_pair(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return np.nan, np.nan
    if np.nanstd(x[m]) == 0 or np.nanstd(y[m]) == 0:
        return np.nan, np.nan
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p)


def partial_spearman_adj_group(x: np.ndarray, y: np.ndarray, g: np.ndarray) -> tuple[float, float]:
    """Residualize ranks of x,y on binary/group indicator, then Spearman."""
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(g.astype(float))
    if m.sum() < 6:
        return np.nan, np.nan
    xr = stats.rankdata(x[m]).astype(float)
    yr = stats.rankdata(y[m]).astype(float)
    g1 = g[m].astype(float)
    # OLS residual vs group
    gdm = np.column_stack([np.ones(len(g1)), g1])
    bx, *_ = np.linalg.lstsq(gdm, xr, rcond=None)
    by, *_ = np.linalg.lstsq(gdm, yr, rcond=None)
    rx = xr - gdm @ bx
    ry = yr - gdm @ by
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan, np.nan
    r, p = stats.spearmanr(rx, ry)
    return float(r), float(p)


def load_meta() -> pd.DataFrame:
    env = pd.read_csv(P.ENV, sep="\t", index_col=0)
    meta = env.T.copy()
    meta.index.name = "Sample"
    meta = meta.reset_index()
    grp = pd.read_csv(P.GROUP, sep="\t", header=None, names=["Sample", "Group"])
    meta = meta.merge(grp, on="Sample", how="inner")
    for c in ["Plaque_Index", "Daily_Brushing_Frequency", "Enamel_Demineralization", "Decayed_Teeth"]:
        meta[c] = pd.to_numeric(meta[c], errors="coerce")
    meta["group_bin"] = (meta["Group"] == "PoorPlaque").astype(int)
    return meta


def load_ko_unstratified(samples: list[str]) -> pd.DataFrame:
    """Load unstratified KO rows, rename columns to sample IDs."""
    # map sample -> column name
    header = pd.read_csv(P.KO_MAT, sep="\t", nrows=0).columns.tolist()
    feat_col = header[0]
    colmap = {}
    for s in samples:
        hits = [h for h in header if h.startswith(f"{s}_") or h == s]
        if hits:
            colmap[s] = hits[0]
    usecols = [feat_col] + list(colmap.values())
    # stream filter unstratified
    chunks = []
    for chunk in pd.read_csv(P.KO_MAT, sep="\t", usecols=usecols, chunksize=5000):
        chunk = chunk[~chunk[feat_col].astype(str).str.contains(r"\|", regex=True)]
        chunk = chunk[chunk[feat_col].astype(str).str.startswith("K")]
        if len(chunk):
            chunks.append(chunk)
    ko = pd.concat(chunks, ignore_index=True)
    ko = ko.rename(columns={feat_col: "KO_ID"})
    ko = ko.rename(columns={v: k for k, v in colmap.items()})
    ko = ko.set_index("KO_ID")
    # collapse duplicate KO if any
    ko = ko.groupby(level=0).sum()
    return ko[samples]


def extract_ko_id(s: str) -> str:
    return str(s).split(":")[0].strip().split(" ")[0]


# ===================== Module A =====================
def run_module_A(meta: pd.DataFrame) -> None:
    print("\n=== Module A: PI associations ===")
    samples = meta["Sample"].tolist()
    pi = meta.set_index("Sample")["Plaque_Index"]
    gbin = meta.set_index("Sample")["group_bin"]
    group = meta.set_index("Sample")["Group"]

    # A0 flagellar ready matrix
    ready = pd.read_csv(P.PI_READY, sep="\t")
    rows = []
    for feat in ["K02416", "K02404", "FAD"]:
        x = ready.set_index("sample")[feat].reindex(samples).values.astype(float)
        y = pi.reindex(samples).values.astype(float)
        gb = gbin.reindex(samples).values
        r_full, p_full = spearman_pair(x, y)
        r_part, p_part = partial_spearman_adj_group(x, y, gb)
        for grp_name in ["Normal", "PoorPlaque"]:
            m = group.reindex(samples).values == grp_name
            r_w, p_w = spearman_pair(x[m], y[m])
            rows.append(
                {
                    "feature": feat,
                    "scope": f"within_{grp_name}",
                    "n": int(m.sum()),
                    "rho": r_w,
                    "pval": p_w,
                }
            )
        rows.append({"feature": feat, "scope": "full", "n": len(samples), "rho": r_full, "pval": p_full})
        rows.append({"feature": feat, "scope": "partial_adj_group", "n": len(samples), "rho": r_part, "pval": p_part})
    flag_df = pd.DataFrame(rows)
    dual_write(flag_df, "Table_A0_flagellar_FAD_PI_scopes")

    # Load KO + metab DA
    print("  loading KO matrix...")
    ko = load_ko_unstratified(samples)
    da_ko = pd.read_csv(P.DA_KO)
    da_ko["KO"] = da_ko["KO_ID"].map(extract_ko_id)
    # match by KO prefix in index
    ko.index = [extract_ko_id(i) for i in ko.index]
    ko = ko.groupby(level=0).sum()
    da_kos = [k for k in da_ko["KO"].unique() if k in ko.index]
    print(f"  Robust_DA KO matched: {len(da_kos)}/{da_ko['KO'].nunique()}")

    metab = pd.read_csv(P.METAB_AB, sep="\t", index_col=0)
    metab = metab[[c for c in samples if c in metab.columns]]
    da_m = pd.read_csv(P.DA_METAB, sep="\t")
    da_metabs = [m for m in da_m["metab_id"] if m in metab.index]
    print(f"  Robust_DA metab matched: {len(da_metabs)}/{len(da_m)}")

    def feature_scopes(mat: pd.DataFrame, feats: list[str], ftype: str) -> pd.DataFrame:
        out = []
        for f in feats:
            x = mat.loc[f, samples].values.astype(float)
            y = pi.reindex(samples).values.astype(float)
            gb = gbin.reindex(samples).values
            r_f, p_f = spearman_pair(x, y)
            r_p, p_p = partial_spearman_adj_group(x, y, gb)
            out.append({"feature": f, "type": ftype, "scope": "full", "n": 20, "rho": r_f, "pval": p_f})
            out.append({"feature": f, "type": ftype, "scope": "partial_adj_group", "n": 20, "rho": r_p, "pval": p_p})
            for grp_name in ["Normal", "PoorPlaque"]:
                m = group.reindex(samples).values == grp_name
                r_w, p_w = spearman_pair(x[m], y[m])
                out.append({"feature": f, "type": ftype, "scope": f"within_{grp_name}", "n": int(m.sum()), "rho": r_w, "pval": p_w})
        df = pd.DataFrame(out)
        # FDR within scope+type for full and partial
        for sc in ["full", "partial_adj_group"]:
            idx = df["scope"] == sc
            df.loc[idx, "FDR"] = bh_fdr(df.loc[idx, "pval"].values)
        return df

    a_ko = feature_scopes(ko, da_kos, "KO")
    a_mb = feature_scopes(metab, da_metabs, "metabolite")
    a_all = pd.concat([a_ko, a_mb], ignore_index=True)
    dual_write(a_all, "Table_A1_DA_feature_PI_all_scopes")

    # Summary: |rho|>=0.8 & P<0.05 counts by scope
    summ_rows = []
    for ftype in ["KO", "metabolite"]:
        sub = a_all[a_all["type"] == ftype]
        nfeat = sub["feature"].nunique()
        for sc in ["full", "partial_adj_group", "within_Normal", "within_PoorPlaque"]:
            s = sub[sub["scope"] == sc]
            hit = s[(s["rho"].abs() >= P.RHO_PI) & (s["pval"] < 0.05)]
            hit_fdr = s[(s["rho"].abs() >= P.RHO_PI) & (s.get("FDR", pd.Series(np.nan, index=s.index)) <= 0.05)] if "FDR" in s else hit.iloc[0:0]
            summ_rows.append(
                {
                    "type": ftype,
                    "n_features": nfeat,
                    "scope": sc,
                    "n_abs_rho_ge_0.8_p_lt_0.05": int(len(hit)),
                    "n_abs_rho_ge_0.8_FDR_le_0.05": int(len(hit_fdr)) if sc in ("full", "partial_adj_group") else np.nan,
                    "median_abs_rho": float(np.nanmedian(s["rho"].abs())),
                    "median_rho": float(np.nanmedian(s["rho"])),
                }
            )
    dual_write(pd.DataFrame(summ_rows), "Table_A2_PI_association_scope_summary")

    # Figure: full vs within |rho| for KO
    plot_rows = []
    for ftype, label in [("KO", "Robust_DA KO"), ("metabolite", "Robust_DA metabolite")]:
        wide = a_all[a_all["type"] == ftype].pivot_table(index="feature", columns="scope", values="rho", aggfunc="first")
        if "full" not in wide.columns:
            continue
        for sc in ["within_Normal", "within_PoorPlaque", "partial_adj_group"]:
            if sc not in wide.columns:
                continue
            for feat, row in wide.iterrows():
                plot_rows.append(
                    {
                        "type": ftype,
                        "feature": feat,
                        "rho_full": row["full"],
                        "rho_alt": row[sc],
                        "alt_scope": sc,
                    }
                )
    plot_df = pd.DataFrame(plot_rows)
    dual_write(plot_df, "Fig_A2_rho_full_vs_alt_plotdata")

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
    for ax, ftype, title in zip(axes, ["KO", "metabolite"], ["Robust_DA KOs", "Robust_DA metabolites"]):
        d = plot_df[(plot_df["type"] == ftype) & (plot_df["alt_scope"] == "partial_adj_group")]
        ax.scatter(d["rho_full"], d["rho_alt"], s=10, alpha=0.5, c="#4C72B0", edgecolors="none")
        ax.axhline(0, color="0.5", lw=0.6)
        ax.axvline(0, color="0.5", lw=0.6)
        ax.plot([-1, 1], [-1, 1], ls="--", color="0.6", lw=0.8)
        ax.set_xlim(-1.05, 1.05)
        ax.set_ylim(-1.05, 1.05)
        ax.set_xlabel(r"Spearman $\rho$ (full cohort)")
        ax.set_ylabel(r"Partial $\rho$ (adj. group)")
        ax.set_title(title)
        # annotate collapse
        n_full = int(((d["rho_full"].abs() >= 0.8)).sum())
        n_part = int(((d["rho_alt"].abs() >= 0.8)).sum())
        ax.text(0.02, 0.98, f"|ρ|≥0.8: full={n_full}\npartial={n_part}", transform=ax.transAxes, va="top", fontsize=8)
    fig.tight_layout()
    save_fig(fig, "Fig_A2_rho_full_vs_partial_group")

    # within-group collapse figure
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
    for ax, sc, title in zip(
        axes,
        ["within_Normal", "within_PoorPlaque"],
        ["Within Normal (n=10)", "Within PoorPlaque (n=10)"],
    ):
        d = plot_df[(plot_df["type"] == "KO") & (plot_df["alt_scope"] == sc)]
        ax.scatter(d["rho_full"], d["rho_alt"], s=10, alpha=0.45, c="#55A868", edgecolors="none")
        ax.axhline(0, color="0.5", lw=0.6)
        ax.axvline(0, color="0.5", lw=0.6)
        ax.set_xlim(-1.05, 1.05)
        ax.set_ylim(-1.05, 1.05)
        ax.set_xlabel(r"$\rho$ full cohort")
        ax.set_ylabel(r"$\rho$ within group")
        ax.set_title(title)
        n_full = int((d["rho_full"].abs() >= 0.8).sum())
        n_w = int((d["rho_alt"].abs() >= 0.8).sum())
        ax.text(0.02, 0.98, f"KO |ρ|≥0.8:\nfull={n_full}, within={n_w}", transform=ax.transAxes, va="top", fontsize=8)
    fig.tight_layout()
    save_fig(fig, "Fig_A2_rho_full_vs_within_KO")


# ===================== Module B =====================
def corr_network_from_matrix(X: np.ndarray, names: list[str], r_min=P.R_MIN, p_max=P.P_MAX):
    """X shape (n_features, n_samples). Return edge list and density stats."""
    n_f, n_s = X.shape
    if n_s < 4:
        return [], {"n_nodes": n_f, "n_edges": 0, "density": 0.0, "mean_degree": 0.0}
    # spearman via rank + pearson
    ranks = np.apply_along_axis(stats.rankdata, 1, X)
    # constant rows
    sd = ranks.std(axis=1)
    valid = sd > 0
    ranks2 = ranks.copy()
    ranks2[~valid] = np.nan
    # pairwise
    # use numpy corrcoef on valid subset
    idx = np.where(valid)[0]
    if len(idx) < 2:
        return [], {"n_nodes": int(valid.sum()), "n_edges": 0, "density": 0.0, "mean_degree": 0.0}
    R = np.corrcoef(ranks2[idx])
    # p-values approx via t distribution
    dof = n_s - 2
    with np.errstate(divide="ignore", invalid="ignore"):
        t = R * np.sqrt(dof / np.clip(1 - R**2, 1e-12, None))
        Pmat = 2 * stats.t.sf(np.abs(t), dof)
    edges = []
    m = len(idx)
    for i in range(m):
        for j in range(i + 1, m):
            r = R[i, j]
            p = Pmat[i, j]
            if np.isfinite(r) and np.isfinite(p) and abs(r) >= r_min and p < p_max:
                edges.append((names[idx[i]], names[idx[j]], float(r), float(p)))
    n_nodes = int(valid.sum())
    n_edges = len(edges)
    max_e = n_nodes * (n_nodes - 1) / 2
    density = n_edges / max_e if max_e else 0.0
    mean_deg = 2 * n_edges / n_nodes if n_nodes else 0.0
    return edges, {"n_nodes": n_nodes, "n_edges": n_edges, "density": density, "mean_degree": mean_deg}


def network_from_edges(edges, all_nodes=None):
    G = nx.Graph()
    if all_nodes is not None:
        G.add_nodes_from(all_nodes)
    for a, b, r, p in edges:
        G.add_edge(a, b, weight=r, pval=p)
    return G


def top_hubs(G: nx.Graph, k=20) -> list[str]:
    if G.number_of_nodes() == 0:
        return []
    deg = dict(G.degree())
    return [n for n, _ in sorted(deg.items(), key=lambda x: (-x[1], x[0]))[:k]]


def run_module_B(meta: pd.DataFrame) -> None:
    print("\n=== Module B: MAG network robustness ===")
    mag = pd.read_csv(P.MAG_AB, sep="\t", index_col=0)
    samples = meta["Sample"].tolist()
    mag = mag[[c for c in samples if c in mag.columns]]
    names = mag.index.tolist()
    X_all = mag.values.astype(float)

    # B0 reproduce manuscript edges
    eN = pd.read_csv(P.EDGE_N)
    eP = pd.read_csv(P.EDGE_P)
    b0 = pd.DataFrame(
        [
            {"group": "Normal", "source": "Gephi_edge", "n_edges": len(eN), "n_nodes_manuscript": 123},
            {"group": "PoorPlaque", "source": "Gephi_edge", "n_edges": len(eP), "n_nodes_manuscript": 125},
            {"group": "Normal", "source": "manuscript_text", "n_edges": 771, "n_nodes_manuscript": 123},
            {"group": "PoorPlaque", "source": "manuscript_text", "n_edges": 1426, "n_nodes_manuscript": 125},
            {"group": "Normal", "source": "Fig5A_reported", "n_edges": 387, "n_nodes_manuscript": np.nan},
            {"group": "PoorPlaque", "source": "Fig5A_reported", "n_edges": 785, "n_nodes_manuscript": np.nan},
        ]
    )
    # rebuild from abundance
    for grp, label in [("Normal", "Normal"), ("PoorPlaque", "PoorPlaque")]:
        smp = meta.loc[meta["Group"] == grp, "Sample"].tolist()
        edges, st = corr_network_from_matrix(mag[smp].values.astype(float), names)
        b0 = pd.concat(
            [
                b0,
                pd.DataFrame(
                    [
                        {
                            "group": label,
                            "source": "rebuilt_spearman_r0.5_p0.05",
                            "n_edges": st["n_edges"],
                            "n_nodes_manuscript": st["n_nodes"],
                            "density": st["density"],
                            "mean_degree": st["mean_degree"],
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )
        # write rebuilt edge list
        edf = pd.DataFrame(edges, columns=["source", "target", "r", "p"])
        dual_write(edf, f"Table_B0_rebuilt_edges_{label}")
    dual_write(b0, "Table_B0_network_summary_reproduce")

    # Observed topology
    obs = {}
    obs_edges = {}
    obs_hubs = {}
    for grp in ["Normal", "PoorPlaque"]:
        smp = meta.loc[meta["Group"] == grp, "Sample"].tolist()
        edges, st = corr_network_from_matrix(mag[smp].values.astype(float), names)
        G = network_from_edges(edges, all_nodes=names)
        obs[grp] = st
        obs_edges[grp] = edges
        obs_hubs[grp] = top_hubs(G, 20)
        print(f"  observed {grp}: edges={st['n_edges']} density={st['density']:.4f} mean_deg={st['mean_degree']:.2f}")

    delta_obs = {
        "delta_density": obs["PoorPlaque"]["density"] - obs["Normal"]["density"],
        "delta_n_edges": obs["PoorPlaque"]["n_edges"] - obs["Normal"]["n_edges"],
        "delta_mean_degree": obs["PoorPlaque"]["mean_degree"] - obs["Normal"]["mean_degree"],
    }

    # B3 Bootstrap edge support + hub stability
    print(f"  bootstrap B={P.N_BOOT} ...")
    rng = np.random.default_rng(P.SEED)
    edge_support = {grp: {} for grp in ["Normal", "PoorPlaque"]}
    hub_support = {grp: {h: 0 for h in names} for grp in ["Normal", "PoorPlaque"]}
    boot_stats = {grp: [] for grp in ["Normal", "PoorPlaque"]}

    for grp in ["Normal", "PoorPlaque"]:
        smp = meta.loc[meta["Group"] == grp, "Sample"].tolist()
        Xg = mag[smp].values.astype(float)
        n = Xg.shape[1]
        for b in range(P.N_BOOT):
            idx = rng.choice(n, size=n, replace=True)
            edges, st = corr_network_from_matrix(Xg[:, idx], names)
            boot_stats[grp].append(st)
            for a, c, r, p in edges:
                key = tuple(sorted([a, c]))
                edge_support[grp][key] = edge_support[grp].get(key, 0) + 1
            G = network_from_edges(edges, all_nodes=names)
            for h in top_hubs(G, 20):
                hub_support[grp][h] += 1
            if (b + 1) % 100 == 0:
                print(f"    {grp} bootstrap {b+1}/{P.N_BOOT}")

    # LOO
    print("  leave-one-out ...")
    loo_support = {grp: {} for grp in ["Normal", "PoorPlaque"]}
    for grp in ["Normal", "PoorPlaque"]:
        smp = meta.loc[meta["Group"] == grp, "Sample"].tolist()
        Xg = mag[smp].values.astype(float)
        n = Xg.shape[1]
        for i in range(n):
            keep = [j for j in range(n) if j != i]
            edges, st = corr_network_from_matrix(Xg[:, keep], names)
            for a, c, r, p in edges:
                key = tuple(sorted([a, c]))
                loo_support[grp][key] = loo_support[grp].get(key, 0) + 1

    # Edge support tables
    for grp in ["Normal", "PoorPlaque"]:
        rows = []
        # union of observed + bootstrap keys
        keys = set(edge_support[grp]) | {tuple(sorted([a, c])) for a, c, _, _ in obs_edges[grp]}
        for a, c in keys:
            bs = edge_support[grp].get((a, c), 0) / P.N_BOOT
            lo = loo_support[grp].get((a, c), 0) / 10.0
            in_obs = (a, c) in {tuple(sorted([x, y])) for x, y, _, _ in obs_edges[grp]} or (
                c,
                a,
            ) in {
                tuple(sorted([x, y])) for x, y, _, _ in obs_edges[grp]
            }
            # fix in_obs
            in_obs = tuple(sorted([a, c])) in {tuple(sorted([x, y])) for x, y, _, _ in obs_edges[grp]}
            rows.append(
                {
                    "group": grp,
                    "source": a,
                    "target": c,
                    "in_observed_network": int(in_obs),
                    "bootstrap_support": bs,
                    "loo_support": lo,
                    "robust_edge": int(bs >= 0.8 and lo >= 0.8),
                }
            )
        edf = pd.DataFrame(rows).sort_values(["robust_edge", "bootstrap_support"], ascending=[False, False])
        dual_write(edf, f"Table_B3_edge_support_{grp}")
        # observed edge support summary
        obs_sup = edf[edf["in_observed_network"] == 1]
        print(
            f"  {grp} observed edges with boot>=0.8 & loo>=0.8: {int(obs_sup['robust_edge'].sum())}/{len(obs_sup)}"
        )

    # Hub stability
    hub_rows = []
    for grp in ["Normal", "PoorPlaque"]:
        for rank, h in enumerate(obs_hubs[grp], 1):
            hub_rows.append(
                {
                    "group": grp,
                    "MAG": h,
                    "obs_degree_rank": rank,
                    "bootstrap_top20_freq": hub_support[grp][h] / P.N_BOOT,
                }
            )
    hub_df = pd.DataFrame(hub_rows)
    dual_write(hub_df, "Table_B3_hub_stability")

    # Bootstrap CI for topology
    ci_rows = []
    for grp in ["Normal", "PoorPlaque"]:
        for metric in ["n_edges", "density", "mean_degree"]:
            vals = np.array([s[metric] for s in boot_stats[grp]], dtype=float)
            ci_rows.append(
                {
                    "group": grp,
                    "metric": metric,
                    "observed": obs[grp][metric],
                    "boot_median": float(np.median(vals)),
                    "boot_ci2.5": float(np.percentile(vals, 2.5)),
                    "boot_ci97.5": float(np.percentile(vals, 97.5)),
                }
            )
    dual_write(pd.DataFrame(ci_rows), "Table_B4_bootstrap_topology_CI")

    # B4 label permutation
    print(f"  label permutation B={P.N_PERM} ...")
    rng2 = np.random.default_rng(P.SEED + 1)
    labels = meta.set_index("Sample").loc[samples, "Group"].values.copy()
    X = mag[samples].values.astype(float)
    null_delta = {k: [] for k in delta_obs}
    for b in range(P.N_PERM):
        lab = labels.copy()
        rng2.shuffle(lab)
        smpN = [s for s, g in zip(samples, lab) if g == "Normal"]
        smpP = [s for s, g in zip(samples, lab) if g == "PoorPlaque"]
        # map to column indices
        idxN = [samples.index(s) for s in smpN]
        idxP = [samples.index(s) for s in smpP]
        _, stN = corr_network_from_matrix(X[:, idxN], names)
        _, stP = corr_network_from_matrix(X[:, idxP], names)
        null_delta["delta_density"].append(stP["density"] - stN["density"])
        null_delta["delta_n_edges"].append(stP["n_edges"] - stN["n_edges"])
        null_delta["delta_mean_degree"].append(stP["mean_degree"] - stN["mean_degree"])
        if (b + 1) % 100 == 0:
            print(f"    perm {b+1}/{P.N_PERM}")

    perm_rows = []
    for k, obs_v in delta_obs.items():
        null = np.array(null_delta[k], dtype=float)
        p_emp = float((np.abs(null) >= abs(obs_v)).mean())
        perm_rows.append(
            {
                "metric": k,
                "observed_Poor_minus_Normal": obs_v,
                "null_mean": float(null.mean()),
                "null_sd": float(null.std()),
                "P_empirical": p_emp,
                "null_ci2.5": float(np.percentile(null, 2.5)),
                "null_ci97.5": float(np.percentile(null, 97.5)),
            }
        )
    perm_df = pd.DataFrame(perm_rows)
    dual_write(perm_df, "Table_B4_permutation_topology")

    # Figures
    # edge support hist for observed edges
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))
    plot_hist = []
    for ax, grp in zip(axes, ["Normal", "PoorPlaque"]):
        edf = pd.read_csv(P.TABLES / f"Table_B3_edge_support_{grp}.tsv", sep="\t")
        obs_e = edf[edf["in_observed_network"] == 1]
        ax.hist(obs_e["bootstrap_support"], bins=20, range=(0, 1), color="#4C72B0", edgecolor="white", lw=0.4)
        ax.axvline(0.8, color="#C44E52", ls="--", lw=1)
        ax.set_xlabel("Bootstrap edge support")
        ax.set_ylabel("Count (observed edges)")
        ax.set_title(grp)
        frac = (obs_e["bootstrap_support"] >= 0.8).mean() if len(obs_e) else 0
        ax.text(0.02, 0.98, f"support≥0.8: {frac:.1%}", transform=ax.transAxes, va="top", fontsize=8)
        for _, r in obs_e.iterrows():
            plot_hist.append({"group": grp, "bootstrap_support": r["bootstrap_support"], "robust_edge": r["robust_edge"]})
    fig.tight_layout()
    save_fig(fig, "Fig_B3_edge_support_hist", pd.DataFrame(plot_hist))

    # hub stability
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    hub_plot = hub_df.copy()
    x = np.arange(20)
    w = 0.38
    nvals = hub_plot[hub_plot.group == "Normal"].sort_values("obs_degree_rank")["bootstrap_top20_freq"].values
    pvals = hub_plot[hub_plot.group == "PoorPlaque"].sort_values("obs_degree_rank")["bootstrap_top20_freq"].values
    ax.bar(x - w / 2, nvals, w, label="Normal hubs", color="#4C72B0")
    ax.bar(x + w / 2, pvals, w, label="PoorPlaque hubs", color="#C44E52")
    ax.axhline(0.5, color="0.4", ls="--", lw=0.8)
    ax.set_xlabel("Observed degree rank (top 20)")
    ax.set_ylabel("Bootstrap frequency in top-20")
    ax.set_xticks(x)
    ax.set_xticklabels([str(i) for i in range(1, 21)], fontsize=7)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    save_fig(fig, "Fig_B3_hub_stability", hub_plot)

    # permutation null
    fig, axes = plt.subplots(1, 3, figsize=(8.5, 3.0))
    metrics = ["delta_density", "delta_n_edges", "delta_mean_degree"]
    titles = [r"$\Delta$ density", r"$\Delta$ n edges", r"$\Delta$ mean degree"]
    null_plot = []
    for ax, k, title in zip(axes, metrics, titles):
        null = np.array(null_delta[k])
        ax.hist(null, bins=30, color="0.75", edgecolor="white", lw=0.3)
        ax.axvline(delta_obs[k], color="#C44E52", lw=1.2, label="observed")
        p = float((np.abs(null) >= abs(delta_obs[k])).mean())
        ax.set_title(f"{title}\nP={p:.3f}")
        ax.set_xlabel("PoorPlaque − Normal")
        for v in null:
            null_plot.append({"metric": k, "null_value": v, "observed": delta_obs[k], "P_empirical": p})
    axes[0].set_ylabel("Permutation count")
    fig.tight_layout()
    save_fig(fig, "Fig_B4_permutation_null", pd.DataFrame(null_plot))

    # decision summary bits for B
    return perm_df, hub_df


# ===================== Module C =====================
def run_module_C(meta: pd.DataFrame) -> None:
    print("\n=== Module C: cross-omics partial correlations ===")
    samples = meta["Sample"].tolist()
    gbin = meta.set_index("Sample")["group_bin"]
    group = meta.set_index("Sample")["Group"]

    # MAG-metab sig edges from manuscript
    mag_metab = pd.read_csv(
        P.ROOT
        / "01this_study_dataset/sgb/metab_corr/paper_result/06_RobustDA_cooccurrence_network/tables/network_MAG_metab_sig_edges.tsv",
        sep="\t",
    )
    print("  MAG-metab edge cols:", mag_metab.columns.tolist()[:12])
    # try to detect columns
    cols = {c.lower(): c for c in mag_metab.columns}
    # common patterns
    def pick(*cands):
        for c in cands:
            if c in mag_metab.columns:
                return c
            if c.lower() in cols:
                return cols[c.lower()]
        return None

    # load abundances
    mag = pd.read_csv(P.MAG_AB, sep="\t", index_col=0)[samples]
    metab = pd.read_csv(P.METAB_AB, sep="\t", index_col=0)
    metab = metab[[c for c in samples if c in metab.columns]]

    # Focus on MAG3 / MAG95 opposite associations if present
    focus = []
    for mag_id in ["MAG3", "MAG95"]:
        if mag_id not in mag.index:
            continue
        # metabolites connected in sig table — flexible
        sub = mag_metab[
            mag_metab.apply(lambda r: mag_id in r.astype(str).values, axis=1)
        ]
        # find metabolite ids in row
        for _, r in sub.iterrows():
            vals = r.astype(str).tolist()
            metab_hits = [v for v in vals if v.startswith("metab_")]
            if not metab_hits:
                continue
            mid = metab_hits[0]
            if mid not in metab.index:
                continue
            x = mag.loc[mag_id, samples].values.astype(float)
            y = metab.loc[mid, samples].values.astype(float)
            r_f, p_f = spearman_pair(x, y)
            r_p, p_p = partial_spearman_adj_group(x, y, gbin.reindex(samples).values)
            focus.append(
                {
                    "MAG": mag_id,
                    "metabolite": mid,
                    "rho_full": r_f,
                    "pval_full": p_f,
                    "rho_partial_adj_group": r_p,
                    "pval_partial": p_p,
                    "sign_full": np.sign(r_f) if np.isfinite(r_f) else np.nan,
                    "sign_partial": np.sign(r_p) if np.isfinite(r_p) else np.nan,
                }
            )
    if focus:
        fdf = pd.DataFrame(focus).drop_duplicates(subset=["MAG", "metabolite"])
        fdf["sign_consistent"] = (fdf["sign_full"] == fdf["sign_partial"]).astype(int)
        dual_write(fdf, "Table_C1_MAG3_MAG95_metab_partial")
        print(f"  MAG3/95 edges evaluated: {len(fdf)}; sign-consistent after group adj: {fdf['sign_consistent'].sum()}")
    else:
        print("  WARNING: no MAG3/MAG95 metab edges resolved; writing empty note")
        dual_write(pd.DataFrame([{"note": "no edges resolved"}]), "Table_C1_MAG3_MAG95_metab_partial")

    # KO-metab-PI: use flagellar already in A; also PI ego edges re-test
    ego = pd.read_csv(
        P.ROOT
        / "01this_study_dataset/sgb/metab_corr/paper_result/13_RobustDA_KO_metab_PI_network/data/Fig_KOmetabPI_PI_ego_edges.tsv",
        sep="\t",
    )
    print("  PI ego edge cols:", ego.columns.tolist())
    # Re-evaluate features linked to Plaque_Index
    ready = pd.read_csv(P.PI_READY, sep="\t").set_index("sample")
    # load KO for ego KOs
    ko = load_ko_unstratified(samples)
    ko.index = [extract_ko_id(i) for i in ko.index]
    ko = ko.groupby(level=0).sum()
    pi = meta.set_index("Sample")["Plaque_Index"]

    # parse ego nodes connected to PI
    node_cols = [c for c in ego.columns if c.lower() in ("source", "target", "from", "to", "node1", "node2", "a", "b")]
    rows = []
    if len(ego):
        # get unique non-PI nodes from first two identifier-like columns
        id_like = []
        for c in ego.columns:
            if ego[c].dtype == object:
                id_like.append(c)
        if len(id_like) >= 2:
            nodes = set(ego[id_like[0]].astype(str)) | set(ego[id_like[1]].astype(str))
            nodes = [n for n in nodes if n not in ("Plaque_Index", "Plaque Index", "PI") and not n.lower().startswith("plaque")]
            for n in nodes:
                # try KO
                kid = extract_ko_id(n)
                if kid in ko.index:
                    x = ko.loc[kid, samples].values.astype(float)
                    ftype = "KO"
                    feat = kid
                elif n in metab.index:
                    x = metab.loc[n, samples].values.astype(float)
                    ftype = "metabolite"
                    feat = n
                elif n in ready.columns:
                    x = ready.reindex(samples)[n].values.astype(float)
                    ftype = "ready"
                    feat = n
                else:
                    continue
                y = pi.reindex(samples).values.astype(float)
                r_f, p_f = spearman_pair(x, y)
                r_p, p_p = partial_spearman_adj_group(x, y, gbin.reindex(samples).values)
                rows.append(
                    {
                        "feature": feat,
                        "type": ftype,
                        "rho_full": r_f,
                        "pval_full": p_f,
                        "rho_partial_adj_group": r_p,
                        "pval_partial": p_p,
                        "abs_rho_full_ge_0.8": int(abs(r_f) >= 0.8) if np.isfinite(r_f) else 0,
                        "abs_rho_partial_ge_0.8": int(abs(r_p) >= 0.8) if np.isfinite(r_p) else 0,
                    }
                )
    if rows:
        cdf = pd.DataFrame(rows).drop_duplicates(subset=["feature", "type"])
        dual_write(cdf, "Table_C2_PI_ego_partial_corr")
        print(
            f"  PI ego features: {len(cdf)}; |ρ|≥0.8 full={cdf['abs_rho_full_ge_0.8'].sum()} partial={cdf['abs_rho_partial_ge_0.8'].sum()}"
        )
    else:
        dual_write(pd.DataFrame([{"note": "could not parse PI ego nodes"}]), "Table_C2_PI_ego_partial_corr")


def write_decision_summary() -> None:
    print("\n=== Decision summary ===")
    rows = []
    # A
    a2 = pd.read_csv(P.TABLES / "Table_A2_PI_association_scope_summary.tsv", sep="\t")
    for _, r in a2.iterrows():
        rows.append(
            {
                "module": "A",
                "analysis": "PI_association",
                "metric": f"{r['type']}|{r['scope']}|n_ge0.8",
                "manuscript_claim": "PI correlations confirm group differences",
                "validation_result": f"n={int(r['n_abs_rho_ge_0.8_p_lt_0.05'])}, median|ρ|={r['median_abs_rho']:.3f}",
                "action_on_manuscript": "SOFTEN: exploratory only; not independent confirmation",
            }
        )
    # B
    b0 = pd.read_csv(P.TABLES / "Table_B0_network_summary_reproduce.tsv", sep="\t")
    rows.append(
        {
            "module": "B0",
            "analysis": "edge_count",
            "metric": "Gephi_vs_text_vs_Fig5A",
            "manuscript_claim": "771/1426 edges",
            "validation_result": b0.to_csv(index=False).replace("\n", " | "),
            "action_on_manuscript": "KEEP text edge counts; clarify Fig5A as display subset if rebuilt matches Gephi",
        }
    )
    perm = pd.read_csv(P.TABLES / "Table_B4_permutation_topology.tsv", sep="\t")
    for _, r in perm.iterrows():
        action = "KEEP exploratory density contrast" if r["P_empirical"] < 0.05 else "DOWNGRADE density difference claim"
        rows.append(
            {
                "module": "B4",
                "analysis": "label_permutation",
                "metric": r["metric"],
                "manuscript_claim": "PoorPlaque network denser (node Wilcoxon)",
                "validation_result": f"obs={r['observed_Poor_minus_Normal']:.4g}, P_emp={r['P_empirical']:.4g}",
                "action_on_manuscript": action + "; replace node Wilcoxon with sample permutation",
            }
        )
    for grp in ["Normal", "PoorPlaque"]:
        edf = pd.read_csv(P.TABLES / f"Table_B3_edge_support_{grp}.tsv", sep="\t")
        obs = edf[edf["in_observed_network"] == 1]
        frac = float(obs["robust_edge"].mean()) if len(obs) else 0
        rows.append(
            {
                "module": "B3",
                "analysis": "edge_bootstrap_loo",
                "metric": f"{grp}_robust_edge_fraction",
                "manuscript_claim": "stable co-occurrence network",
                "validation_result": f"robust={int(obs['robust_edge'].sum())}/{len(obs)} ({frac:.1%})",
                "action_on_manuscript": "Present as exploratory; report robust subset",
            }
        )
    # C
    c1p = P.TABLES / "Table_C1_MAG3_MAG95_metab_partial.tsv"
    if c1p.exists():
        c1 = pd.read_csv(c1p, sep="\t")
        if "sign_consistent" in c1.columns:
            rows.append(
                {
                    "module": "C",
                    "analysis": "MAG3_MAG95_partial",
                    "metric": "sign_consistency",
                    "manuscript_claim": "antagonistic hubs",
                    "validation_result": f"{int(c1['sign_consistent'].sum())}/{len(c1)} sign-consistent after group adj",
                    "action_on_manuscript": "RENAME to opposite-signed associations; soften mechanism language",
                }
            )
    dual_write(pd.DataFrame(rows), "Table_R1.4_decision_summary")


def main():
    print("R1.4 validation start")
    print("OUT:", P.OUT)
    meta = load_meta()
    print("samples:", len(meta), meta.Group.value_counts().to_dict())
    run_module_A(meta)
    run_module_B(meta)
    run_module_C(meta)
    write_decision_summary()
    print("\nDONE")


if __name__ == "__main__":
    main()
