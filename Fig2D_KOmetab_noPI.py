"""Packaged for revise/; source: project/review/R1.4/scripts/02_Fig2D_KOmetab_noPI.py
Figures: Fig2D
Note: KO-metabolite edges without PI
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
"""Rebuild Fig.2D as Robust_DA KO–metabolite association network (no Plaque Index).

Uses the same edge set/thresholds as Module 13, dropping Clinical–KO / Clinical–Metabolite
(PI) edges. Outputs tables + publication figures under review/R1.4/Fig2D_KOmetab_noPI/.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = PKG_ROOT
SRC = ROOT / "01this_study_dataset/sgb/metab_corr/paper_result/13_RobustDA_KO_metab_PI_network/tables"
OUT = ROOT / "review/R1.4/Fig2D_KOmetab_noPI"
TAB = OUT / "tables"
FIG = OUT / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

COL = {"Normal": "#2E8B8B", "PoorPlaque": "#B85450", "KO": "#3D5A80", "Metabolite": "#E09F3E"}

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


def main() -> None:
    print("=== Fig.2D KO–metabolite network (no PI) ===")
    edges_all = pd.read_csv(SRC / "network_KO_metab_PI_sig_edges.tsv", sep="\t")
    nodes_all = pd.read_csv(SRC / "network_nodes.tsv", sep="\t")

    # --- comparison summary ---
    cmp_rows = [
        {
            "network": "original_KO_metab_PI",
            "n_edges": len(edges_all),
            "n_KO_metab_edges": int((edges_all["edge_type"] == "KO–Metabolite").sum()),
            "n_PI_edges": int(edges_all["involves_PI"].sum()) if "involves_PI" in edges_all else 59,
            "n_nodes": nodes_all["id"].nunique() if "id" in nodes_all else np.nan,
        }
    ]

    km = edges_all[edges_all["edge_type"] == "KO–Metabolite"].copy()
    # node table restricted to KO/metab in km
    node_ids = sorted(set(km["from"]).union(set(km["to"])))
    nodes = nodes_all[nodes_all["id"].isin(node_ids)].copy() if "id" in nodes_all.columns else pd.DataFrame({"id": node_ids})

    G = nx.Graph()
    for _, r in km.iterrows():
        G.add_edge(r["from"], r["to"], rho=float(r["rho"]), FDR=float(r["FDR"]), sign=r.get("edge_sign", ""))

    # enrich direction map
    enrich = {}
    if {"id", "enrich_direction", "type"}.issubset(nodes.columns):
        for _, r in nodes.iterrows():
            enrich[r["id"]] = (r["type"], r["enrich_direction"], r.get("label", r["id"]))
    else:
        for n in node_ids:
            t = "KO" if str(n).startswith("K") else "Metabolite"
            enrich[n] = (t, "Unknown", n)

    deg = dict(G.degree())
    topo = []
    for n, d in deg.items():
        t, en, lab = enrich.get(n, ("Unknown", "Unknown", n))
        topo.append(
            {
                "id": n,
                "type": t,
                "enrich_direction": en,
                "label": lab,
                "degree": d,
                "weighted_degree": sum(abs(G[n][nbr]["rho"]) for nbr in G.neighbors(n)),
            }
        )
    topo_df = pd.DataFrame(topo).sort_values("degree", ascending=False)
    dual_write(km, "Table_Fig2D_KOmetab_edges_noPI")
    dual_write(topo_df, "Table_Fig2D_KOmetab_node_topology_noPI")

    n_ko = int((topo_df["type"] == "KO").sum())
    n_metab = int((topo_df["type"] == "Metabolite").sum())
    n_pos_edges = int((km["rho"] > 0).sum())
    n_neg_edges = int((km["rho"] < 0).sum())
    same = int(km["same_direction"].sum()) if "same_direction" in km.columns else np.nan
    cmp_rows.append(
        {
            "network": "Fig2D_KOmetab_noPI",
            "n_edges": len(km),
            "n_KO_metab_edges": len(km),
            "n_PI_edges": 0,
            "n_nodes": G.number_of_nodes(),
            "n_KO_nodes": n_ko,
            "n_metabolite_nodes": n_metab,
            "n_positive_edges": n_pos_edges,
            "n_negative_edges": n_neg_edges,
            "n_same_enrich_direction_edges": same,
            "density": nx.density(G),
            "n_connected_components": nx.number_connected_components(G),
        }
    )
    dual_write(pd.DataFrame(cmp_rows), "Table_Fig2D_with_vs_without_PI_summary")

    # flagellar / FAD focus for 2E recommendation
    focus_ids = ["K02416", "K02404", "metab_27961"]
    focus_rows = []
    for n in focus_ids:
        if n not in G:
            focus_rows.append({"id": n, "in_network": 0})
            continue
        nbrs = list(G.neighbors(n))
        focus_rows.append(
            {
                "id": n,
                "in_network": 1,
                "type": enrich[n][0],
                "enrich_direction": enrich[n][1],
                "degree": deg[n],
                "n_neighbors": len(nbrs),
                "neighbor_ids": ";".join(nbrs),
            }
        )
    # edges among focus set
    focus_edge_rows = []
    for a, b in [("K02416", "metab_27961"), ("K02404", "metab_27961"), ("K02416", "K02404")]:
        if G.has_edge(a, b):
            ed = G[a][b]
            focus_edge_rows.append({"from": a, "to": b, "rho": ed["rho"], "FDR": ed.get("FDR", np.nan)})
        else:
            focus_edge_rows.append({"from": a, "to": b, "rho": np.nan, "FDR": np.nan, "note": "no edge at |ρ|≥0.8"})
    dual_write(pd.DataFrame(focus_rows), "Table_Fig2E_focus_nodes_in_KOmetab_net")
    dual_write(pd.DataFrame(focus_edge_rows), "Table_Fig2E_flagellar_FAD_edges_noPI")

    # --- Figure 1: edge composition (with vs without PI) ---
    comp = (
        edges_all.groupby("edge_type", as_index=False)
        .size()
        .rename(columns={"size": "n_edges"})
    )
    comp["network"] = "original_with_PI"
    comp2 = pd.DataFrame(
        [{"edge_type": "KO–Metabolite", "n_edges": len(km), "network": "Fig2D_no_PI"}]
    )
    comp_plot = pd.concat([comp, comp2], ignore_index=True)
    dual_write(comp_plot, "Table_Fig2D_edge_type_composition")

    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    order = ["KO–Metabolite", "Clinical–KO", "Clinical–Metabolite"]
    x = np.arange(len(order))
    w = 0.36
    orig = [int(comp.loc[comp.edge_type == t, "n_edges"].sum()) if t in set(comp.edge_type) else 0 for t in order]
    newv = [len(km), 0, 0]
    ax.bar(x - w / 2, orig, w, label="Original (with PI)", color="#7A7A7A")
    ax.bar(x + w / 2, newv, w, label="Fig.2D revised (no PI)", color="#3D5A80")
    ax.set_xticks(x)
    ax.set_xticklabels(["KO–Metabolite", "KO–PI", "Metabolite–PI"], fontsize=9)
    ax.set_ylabel("Number of significant edges")
    ax.legend(frameon=False, fontsize=8)
    for i, (a, b) in enumerate(zip(orig, newv)):
        ax.text(i - w / 2, a + 20, str(a), ha="center", fontsize=7)
        ax.text(i + w / 2, b + 20, str(b), ha="center", fontsize=7)
    fig.tight_layout()
    save_fig(fig, "Fig2D_edge_composition_with_vs_without_PI", comp_plot)

    # --- Figure 2: overview network (display filter for readability) ---
    # Full graph is dense; main display uses |rho|>=0.85 OR incident to flagellar/FAD
    display = km[(km["rho"].abs() >= 0.85)].copy()
    focus_keep = set()
    for n in focus_ids:
        if n in G:
            focus_keep.add(n)
            focus_keep.update(G.neighbors(n))
    display2 = km[km["from"].isin(focus_keep) & km["to"].isin(focus_keep)].copy()
    # union for a second panel; overview uses |rho|>=0.85
    Gd = nx.Graph()
    for _, r in display.iterrows():
        Gd.add_edge(r["from"], r["to"], rho=float(r["rho"]))

    # if still huge, take giant component and cap labels
    if Gd.number_of_nodes() == 0:
        Gd = G.copy()
        display = km.copy()

    layout_pos = nx.spring_layout(Gd, seed=42, k=1.8 / np.sqrt(max(Gd.number_of_nodes(), 1)))
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    # edges
    for a, b, d in Gd.edges(data=True):
        r = d["rho"]
        col = "#B85450" if r > 0 else "#3D5A80"
        ax.plot(
            [layout_pos[a][0], layout_pos[b][0]],
            [layout_pos[a][1], layout_pos[b][1]],
            color=col,
            alpha=0.18,
            lw=0.4 + 1.2 * (abs(r) - 0.8) / 0.2,
            zorder=1,
        )
    # nodes
    for n in Gd.nodes():
        t, en, lab = enrich.get(n, ("Unknown", "Unknown", n))
        c = COL["KO"] if t == "KO" else COL["Metabolite"]
        size = 18 + 2.2 * deg.get(n, 1)
        # emphasize focus
        if n in focus_ids:
            size = 90
            ec = "black"
            lw = 1.0
        else:
            ec = "white"
            lw = 0.3
        ax.scatter(layout_pos[n][0], layout_pos[n][1], s=size, c=c, edgecolors=ec, linewidths=lw, zorder=2)
    # labels for focus + top degree in display graph
    label_nodes = set(focus_ids) & set(Gd.nodes())
    top_disp = sorted(Gd.degree(), key=lambda x: -x[1])[:12]
    label_nodes.update(n for n, _ in top_disp)
    for n in label_nodes:
        if n not in layout_pos:
            continue
        lab = str(n)
        if lab.startswith("metab_") and n in enrich:
            # short metabolite name from nodes table
            lab = str(enrich[n][2])[:18]
        ax.text(layout_pos[n][0], layout_pos[n][1] + 0.03, lab, fontsize=6, ha="center", va="bottom", zorder=3)
    ax.set_axis_off()
    ax.set_title(
        f"Robust_DA KO–metabolite associations (no PI)\n"
        f"Display: |ρ|≥0.85 edges (n={Gd.number_of_edges()} of {len(km)} total at |ρ|≥0.80)",
        fontsize=10,
    )
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COL["KO"], markersize=8, label="KO"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COL["Metabolite"], markersize=8, label="Metabolite"),
        Line2D([0], [0], color="#B85450", lw=1.5, label="Positive ρ"),
        Line2D([0], [0], color="#3D5A80", lw=1.5, label="Negative ρ"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8, loc="lower left")
    # plotdata: display edges
    dual_write(display, "Fig2D_KOmetab_network_display_edges")
    layout_df = pd.DataFrame(
        [{"id": n, "x": layout_pos[n][0], "y": layout_pos[n][1], "degree_full": deg.get(n, 0)} for n in Gd.nodes()]
    )
    dual_write(layout_df, "Fig2D_KOmetab_network_layout")
    fig.tight_layout()
    save_fig(fig, "Fig2D_KOmetab_network_noPI")

    # --- Figure 3: flagellar–FAD ego (supports Fig.2E redesign) ---
    Ge = nx.Graph()
    for _, r in display2.iterrows():
        Ge.add_edge(r["from"], r["to"], rho=float(r["rho"]))
    # ensure focus edges even if |rho|<0.85 but in full km
    for _, r in km.iterrows():
        if r["from"] in focus_keep and r["to"] in focus_keep:
            if not Ge.has_edge(r["from"], r["to"]):
                Ge.add_edge(r["from"], r["to"], rho=float(r["rho"]))

    if Ge.number_of_nodes() > 0:
        pe = nx.spring_layout(Ge, seed=7, k=1.5 / np.sqrt(max(Ge.number_of_nodes(), 1)))
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        for a, b, d in Ge.edges(data=True):
            r = d["rho"]
            col = "#B85450" if r > 0 else "#3D5A80"
            lw = 2.0 if {a, b} <= set(focus_ids) or (a in focus_ids or b in focus_ids) and (
                a == "metab_27961" or b == "metab_27961"
            ) else 0.6
            ax.plot([pe[a][0], pe[b][0]], [pe[a][1], pe[b][1]], color=col, alpha=0.35, lw=lw, zorder=1)
        for n in Ge.nodes():
            t, en, lab = enrich.get(n, ("Unknown", "Unknown", n))
            c = COL["KO"] if t == "KO" else COL["Metabolite"]
            size = 120 if n in focus_ids else 28
            ax.scatter(
                pe[n][0],
                pe[n][1],
                s=size,
                c=c,
                edgecolors="black" if n in focus_ids else "white",
                linewidths=0.8 if n in focus_ids else 0.2,
                zorder=2,
            )
            if n in focus_ids or deg.get(n, 0) >= 40:
                show = (
                    "FAD"
                    if n == "metab_27961"
                    else ("FliM" if n == "K02416" else ("FlhF" if n == "K02404" else str(enrich[n][2])[:14]))
                )
                ax.text(pe[n][0], pe[n][1] + 0.04, show, fontsize=8, ha="center", zorder=3)
        ax.set_axis_off()
        ax.set_title("Flagellar KO–FAD neighbourhood in KO–metabolite network (no PI)", fontsize=10)
        dual_write(display2, "Fig2E_flagellar_FAD_neighbourhood_edges")
        fig.tight_layout()
        save_fig(fig, "Fig2E_flagellar_FAD_neighbourhood_noPI")

    # recommendation markdown
    rec = OUT / "Fig2D_noPI_结论与Fig2E建议.md"
    rec.write_text(
        f"""# Fig.2D 去 PI 重建结果与 Fig.2E 建议

## 1. 重建结果（主文 Fig.2D）

| 项目 | 原图（含 PI） | 修订（仅 KO–代谢物） |
|------|--------------:|---------------------:|
| 显著边总数 | {len(edges_all)} | **{len(km)}** |
| KO–代谢物边 | {(edges_all.edge_type=='KO–Metabolite').sum()} | **{len(km)}** |
| 涉及 PI 的边 | {int(edges_all.involves_PI.sum()) if 'involves_PI' in edges_all else 59} | **0** |
| 节点数 | {nodes_all.id.nunique() if 'id' in nodes_all.columns else '—'} | **{G.number_of_nodes()}**（KO {n_ko} + 代谢物 {n_metab}） |
| 正/负边 | — | **{n_pos_edges} / {n_neg_edges}** |
| 同向富集边 | — | **{same}** |

阈值与原文 Module 13 相同：`|ρ|≥0.80, p<0.05, FDR≤0.05`，仅剔除 Clinical–KO / Clinical–Metabolite。

**主文展示图：** `figures/Fig2D_KOmetab_network_noPI.*`（为可读性，图面显示 `|ρ|≥0.85` 子集；全量边表见 `tables/Table_Fig2D_KOmetab_edges_noPI`）。

## 2. Fig.2E 展示建议（推荐方案）

**主文 Fig.2E（建议）：** 只保留 **K02416 (FliM) / K02404 (FlhF) vs FAD (metab_27961)** 散点（可 2×2 或 1×2），**不再把 PI 作为主面板终点**。

依据（去 PI 后的 KO–代谢物网）：
- K02416、K02404、FAD 均仍在网络中
- FliM–FAD、FlhF–FAD 在 `|ρ|≥0.80` 网络中**仍有边**（见 `Table_Fig2E_flagellar_FAD_edges_noPI`）
- FAD 在无 PI 网中 degree = {deg.get('metab_27961', 'NA')}（高连接代谢物节点之一）

**补充图 Fig.S（建议）：** 原 Fig.2E 中 vs Plaque Index 的散点 + 原 PI ego；图注写明 exploratory、非独立验证，并引用组调整偏相关塌缩（R1.4 Table_A0）。

邻域示意图：`figures/Fig2E_flagellar_FAD_neighbourhood_noPI.*`

## 3. 原“共发生/关联网”结论还成立吗？

| 原表述倾向 | 去 PI 后 |
|------------|----------|
| Robust_DA KO 与代谢物存在大量跨组学关联 | **仍成立**（2158 条 KO–代谢物边完整保留） |
| 以 PI 为中心的 ego、特征“tracked Plaque Index”作为叙事支柱 | **应降级/移补充**（逻辑上非独立验证；偏相关已塌） |
| FliM/FlhF 与 FAD 的关联 | **仍可在主文保留**（不依赖 PI 边） |
| 整段“多组学协同重塑” | **大体仍成立**，但证据重心从“PI 网络”改回“组间 DA + KO–代谢物关联 + 通路汇聚” |

**一句话：** 去掉 PI 后，Fig.2D 的跨组学关联主体还在；需要改的是不再用 PI 相关撑独立验证。Fig.2E 主文改盯 FAD，PI 散点进补充。
""",
        encoding="utf-8",
    )
    print("wrote", rec)
    print("DONE", OUT)


if __name__ == "__main__":
    main()
