"""Packaged for revise/; source: project/review/R1.4/scripts/03_Fig2D_KOmetab_noPI_stability.py
Figures: Fig2D
Note: edge bootstrap/LOO stability
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
"""Stability of Fig.2D KO–metabolite network (no PI): bootstrap + LOO + hub support."""
from __future__ import annotations

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
    df.to_csv(TAB / f"{stem}.tsv", sep="\t", index=False)
    df.to_csv(TAB / f"{stem}.csv", index=False)
    print(f"  wrote {stem} ({len(df)} rows)")


def save_fig(fig: plt.Figure, stem: str, plot_df: pd.DataFrame | None = None) -> None:
    if plot_df is not None:
        dual_write(plot_df, f"{stem}_plotdata")
    for ext in ("pdf", "svg", "png", "jpg"):
        fig.savefig(FIG / f"{stem}.{ext}")
    plt.close(fig)
    print(f"  fig {stem}")


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

    metab_rob = pd.read_csv(
        ROOT
        / "01this_study_dataset/sgb/metab_corr/paper_result/05_RobustDA_metab_DA/tables/metabolite_RobustDA_only_with_HMDB_Class.tsv",
        sep="\t",
    )
    metab_ids = metab_rob["metab_id"].astype(str).tolist()

    # KO matrix
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

    metab = pd.read_csv(
        ROOT / "01this_study_dataset/sgb/metab_corr/output/tables/metab_abund_AB_select.tsv",
        sep="\t",
        index_col=0,
    )
    metab = metab[[c for c in samples if c in metab.columns]]
    metab = metab.loc[[m for m in metab_ids if m in metab.index]]

    # align KO list to available
    ko = ko.loc[[k for k in ko_ids if k in ko.index]]
    # drop zero variance
    ko = ko.loc[ko.std(axis=1) > 0]
    metab = metab.loc[metab.std(axis=1) > 0]
    print(f"  KO matrix {ko.shape}, metab matrix {metab.shape}")
    return samples, ko, metab


def cross_spearman_edges(ko_mat: np.ndarray, metab_mat: np.ndarray, ko_names, metab_names):
    """ko_mat/metab_mat: features × samples. Return list of (ko, metab, rho, p, fdr) passing thresholds."""
    n_s = ko_mat.shape[1]
    # rank
    rk = np.apply_along_axis(stats.rankdata, 1, ko_mat)
    rm = np.apply_along_axis(stats.rankdata, 1, metab_mat)
    # standardize rows
    rk = (rk - rk.mean(axis=1, keepdims=True)) / (rk.std(axis=1, keepdims=True) + 1e-12)
    rm = (rm - rm.mean(axis=1, keepdims=True)) / (rm.std(axis=1, keepdims=True) + 1e-12)
    # correlation matrix KO × metab
    R = (rk @ rm.T) / (n_s - 0)  # since already z-scored on ranks; for pearson of ranks use / (n-1)?
    # With z-scored ranks (population std), corr = (z1·z2)/n ; scipy uses sample. Use n-1 for closer match:
    R = (rk @ rm.T) / (n_s - 1)
    dof = n_s - 2
    with np.errstate(divide="ignore", invalid="ignore"):
        t = R * np.sqrt(dof / np.clip(1 - R**2, 1e-12, None))
        P = 2 * stats.t.sf(np.abs(t), dof)
    # flatten
    ii, jj = np.where(np.isfinite(R))
    rho = R[ii, jj]
    p = P[ii, jj]
    fdr = bh_fdr(p)
    keep = (np.abs(rho) >= RHO_MIN) & (p < P_MAX) & (fdr <= FDR_MAX)
    edges = []
    for i, j, r, pv, q, k in zip(ii, jj, rho, p, fdr, keep):
        if k:
            edges.append((ko_names[i], metab_names[j], float(r), float(pv), float(q)))
    return edges, {"n_edges": len(edges), "n_nodes": len({a for a, b, *_ in edges} | {b for a, b, *_ in edges})}


def edge_key(a, b):
    return (a, b)  # directed type pair KO-metab; keep ordered ko,metab


def main():
    print("=== Fig.2D KO–metab (no PI) stability ===")
    samples, ko_df, metab_df = load_matrices()
    ko_names = ko_df.index.tolist()
    metab_names = metab_df.index.tolist()
    Xk = ko_df[samples].values.astype(float)
    Xm = metab_df[samples].values.astype(float)

    # Observed network
    obs_edges, obs_st = cross_spearman_edges(Xk, Xm, ko_names, metab_names)
    print(f"  observed edges: {obs_st['n_edges']}, nodes~{obs_st['n_nodes']}")
    obs_set = {edge_key(a, b) for a, b, *_ in obs_edges}
    obs_df = pd.DataFrame(obs_edges, columns=["KO", "metabolite", "rho", "p", "FDR"])
    dual_write(obs_df, "Table_Fig2D_noPI_observed_edges_rebuilt")

    # Compare to filtered Module-13 edges
    prev = pd.read_csv(TAB / "Table_Fig2D_KOmetab_edges_noPI.tsv", sep="\t")
    prev_set = {(r["from"], r["to"]) for _, r in prev.iterrows()}
    # Module13 from/to may be either order; normalize KO first
    def norm_pair(a, b):
        if str(a).startswith("K"):
            return (a, b)
        return (b, a)

    prev_set = {norm_pair(a, b) for a, b in prev_set}
    overlap = len(obs_set & prev_set)
    print(f"  overlap with prior filtered table: {overlap} / rebuilt {len(obs_set)} / prior {len(prev_set)}")

    # Degrees / hubs observed
    from collections import Counter

    deg = Counter()
    for a, b, *_ in obs_edges:
        deg[a] += 1
        deg[b] += 1
    hubs = [n for n, _ in deg.most_common(20)]

    # Bootstrap
    print(f"  bootstrap B={N_BOOT} ...")
    rng = np.random.default_rng(SEED)
    edge_support = Counter()
    hub_support = Counter({h: 0 for h in deg})  # all nodes init 0; update hubs of interest
    boot_n_edges = []
    n = len(samples)
    for b in range(N_BOOT):
        idx = rng.choice(n, size=n, replace=True)
        edges, st = cross_spearman_edges(Xk[:, idx], Xm[:, idx], ko_names, metab_names)
        boot_n_edges.append(st["n_edges"])
        for a, c, *_ in edges:
            edge_support[edge_key(a, c)] += 1
        # hub top20
        d2 = Counter()
        for a, c, *_ in edges:
            d2[a] += 1
            d2[c] += 1
        top = {n for n, _ in d2.most_common(20)}
        for h in top:
            hub_support[h] += 1
        if (b + 1) % 50 == 0:
            print(f"    bootstrap {b+1}/{N_BOOT}")

    # LOO
    print("  leave-one-out ...")
    loo_support = Counter()
    for i in range(n):
        keep = [j for j in range(n) if j != i]
        edges, st = cross_spearman_edges(Xk[:, keep], Xm[:, keep], ko_names, metab_names)
        for a, c, *_ in edges:
            loo_support[edge_key(a, c)] += 1
        print(f"    LOO left-out sample {samples[i]}: edges={st['n_edges']}")

    # Edge support table for observed edges
    rows = []
    for a, b, r, p, q in obs_edges:
        k = edge_key(a, b)
        bs = edge_support[k] / N_BOOT
        lo = loo_support[k] / float(n)
        rows.append(
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
    edf = pd.DataFrame(rows).sort_values(["robust_edge", "bootstrap_support"], ascending=[False, False])
    dual_write(edf, "Table_Fig2D_noPI_edge_support")

    n_rob = int(edf["robust_edge"].sum())
    frac = n_rob / len(edf) if len(edf) else 0
    summ = pd.DataFrame(
        [
            {
                "network": "Fig2D_KOmetab_noPI",
                "n_observed_edges": len(edf),
                "n_robust_edges_boot0.8_loo0.8": n_rob,
                "robust_fraction": frac,
                "median_bootstrap_support": float(edf["bootstrap_support"].median()),
                "median_loo_support": float(edf["loo_support"].median()),
                "boot_n_edges_median": float(np.median(boot_n_edges)),
                "boot_n_edges_ci2.5": float(np.percentile(boot_n_edges, 2.5)),
                "boot_n_edges_ci97.5": float(np.percentile(boot_n_edges, 97.5)),
                "overlap_with_module13_filtered": overlap,
                "n_rebuilt": len(obs_set),
                "n_module13_filtered": len(prev_set),
            }
        ]
    )
    dual_write(summ, "Table_Fig2D_noPI_stability_summary")

    # Hub stability for observed top20
    hub_rows = []
    for rank, h in enumerate(hubs, 1):
        hub_rows.append(
            {
                "id": h,
                "obs_degree": deg[h],
                "obs_degree_rank": rank,
                "bootstrap_top20_freq": hub_support[h] / N_BOOT,
                "type": "KO" if str(h).startswith("K") else "Metabolite",
            }
        )
    hub_df = pd.DataFrame(hub_rows)
    dual_write(hub_df, "Table_Fig2D_noPI_hub_stability")

    # Focus nodes
    focus = ["K02416", "K02404", "metab_27961"]
    focus_rows = []
    for h in focus:
        sub = edf[(edf["KO"] == h) | (edf["metabolite"] == h)]
        focus_rows.append(
            {
                "id": h,
                "obs_degree": int(deg.get(h, 0)),
                "n_obs_edges_in_table": len(sub),
                "n_robust_edges": int(sub["robust_edge"].sum()) if len(sub) else 0,
                "median_boot_support_of_its_edges": float(sub["bootstrap_support"].median()) if len(sub) else np.nan,
                "bootstrap_top20_freq": hub_support[h] / N_BOOT if h in hub_support else 0.0,
            }
        )
    dual_write(pd.DataFrame(focus_rows), "Table_Fig2D_noPI_focus_node_stability")

    # Figures
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.hist(edf["bootstrap_support"], bins=20, range=(0, 1), color="#3D5A80", edgecolor="white", lw=0.4)
    ax.axvline(0.8, color="#C44E52", ls="--", lw=1)
    ax.set_xlabel("Bootstrap edge support")
    ax.set_ylabel("Count (observed KO–metabolite edges)")
    ax.set_title(f"Fig.2D no-PI edge stability (robust={n_rob}/{len(edf)} = {frac:.1%})")
    fig.tight_layout()
    save_fig(fig, "Fig2D_noPI_edge_support_hist", edf[["KO", "metabolite", "bootstrap_support", "loo_support", "robust_edge"]])

    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    x = np.arange(len(hub_df))
    ax.bar(x, hub_df["bootstrap_top20_freq"], color=np.where(hub_df["type"] == "KO", "#3D5A80", "#E09F3E"))
    ax.axhline(0.5, color="0.4", ls="--", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([str(i) for i in hub_df["obs_degree_rank"]], fontsize=7)
    ax.set_xlabel("Observed degree rank (top 20)")
    ax.set_ylabel("Bootstrap frequency in top-20")
    ax.set_title("Fig.2D no-PI hub stability")
    fig.tight_layout()
    save_fig(fig, "Fig2D_noPI_hub_stability", hub_df)

    # decision note
    note = OUT / "Fig2D_noPI_stability_结论.md"
    note.write_text(
        f"""# 去 PI 后 Fig.2D（KO–代谢物）网络稳定性结论

## 方法
- 重建规则与 Module 13 一致：Robust_DA KO × Robust_DA 代谢物，Spearman，`|ρ|≥0.80`, `P<0.05`, `FDR≤0.05`；**不含 PI**
- Bootstrap：有放回重采样 n=20，B={N_BOOT}（seed={SEED}）
- LOO：逐一去掉 1 个样本
- 稳健边：bootstrap support ≥0.8 且 LOO support ≥0.8

## 结果
| 指标 | 数值 |
|------|------|
| 观察边数（本脚本重建） | {len(edf)} |
| 与先前过滤表重叠 | {overlap} / prior {len(prev_set)} |
| 稳健边 | **{n_rob} / {len(edf)} ({frac:.1%})** |
| 中位 bootstrap support | {float(edf['bootstrap_support'].median()):.3f} |
| 中位 LOO support | {float(edf['loo_support'].median()):.3f} |
| Bootstrap 边数中位 [95% CI] | {float(np.median(boot_n_edges)):.0f} [{float(np.percentile(boot_n_edges,2.5)):.0f}, {float(np.percentile(boot_n_edges,97.5)):.0f}] |

## 解读（对修回）
- 去 PI 后跨组学关联**主体仍在**，但边/hub 在 n=20 下同样**不够稳**（稳健边比例见上）。
- 因此 Fig.2D 即使去掉 PI，仍应标为 **exploratory / hypothesis-generating**，与 MAG 网络（Fig.5）口径一致。
- FliM/FlhF/FAD 相关边的稳健性见 `Table_Fig2D_noPI_focus_node_stability`；主文 Fig.2E 可保留 FAD 散点，但勿过度声称稳定网络枢纽。

## 文件
- `tables/Table_Fig2D_noPI_edge_support.tsv`
- `tables/Table_Fig2D_noPI_hub_stability.tsv`
- `tables/Table_Fig2D_noPI_stability_summary.tsv`
- `figures/Fig2D_noPI_edge_support_hist.*`
- `figures/Fig2D_noPI_hub_stability.*`
""",
        encoding="utf-8",
    )
    print("wrote", note)
    print("DONE")


if __name__ == "__main__":
    main()
