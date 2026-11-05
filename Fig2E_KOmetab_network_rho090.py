"""Packaged for revise/; source: project/review/R1.4/scripts/04_replot_Fig2D_network_rho09.py
Figures: Fig2E
Note: high-|rho| network
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
"""Replot Fig.2D KO–metabolite network (no PI) at |ρ|≥0.9 with clearer layout/labels."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.lines import Line2D

ROOT = PKG_ROOT
OUT = ROOT / "review/R1.4/Fig2D_KOmetab_noPI"
TAB = OUT / "tables"
FIG = OUT / "figures"
RHO_THR = 0.90

COL = {"KO": "#2F4B7C", "Metabolite": "#E07A3D", "pos": "#C0392B", "neg": "#2874A6"}

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
    }
)


def dual_write(df: pd.DataFrame, stem: str) -> None:
    df.to_csv(TAB / f"{stem}.tsv", sep="\t", index=False)
    df.to_csv(TAB / f"{stem}.csv", index=False)
    print(f"  wrote {stem} ({len(df)} rows)")


def short_label(n: str, enrich: dict, max_len: int = 20) -> str:
    if n == "metab_27961":
        return "FAD"
    if n == "K02416":
        return "FliM (K02416)"
    if n == "K02404":
        return "FlhF (K02404)"
    _t, _en, lab = enrich.get(n, ("Unknown", "Unknown", n))
    s = str(lab) if lab else str(n)
    if str(n).startswith("K"):
        # keep KO id primarily; append short name if compact
        kid = str(n)
        name = s.split(":", 1)[-1].strip() if ":" in s else ""
        if name and len(name) <= 12:
            s = f"{kid}"
        else:
            s = kid
    else:
        s = str(lab)
        s = s.replace("Flavin adenine dinucleotide", "FAD")
        s = s.replace("N-Acetyl-", "NAc-").replace("N-acetyl-", "NAc-")
    if len(s) > max_len:
        s = s[: max_len - 1] + "…"
    return s


def layout_by_components(Gd: nx.Graph) -> dict:
    """Place each connected component with size-adaptive spacing."""
    comps = sorted(nx.connected_components(Gd), key=len, reverse=True)
    layout_pos: dict = {}

    # Place components left-to-right, wrapping by cumulative width
    x_cursor, y_cursor = 0.0, 0.0
    row_height = 0.0
    max_row_width = 14.0

    for comp in comps:
        sub = Gd.subgraph(comp).copy()
        m = sub.number_of_nodes()
        # adaptive footprint
        if m <= 2:
            target = 1.8
        elif m <= 4:
            target = 2.6
        elif m <= 7:
            target = 4.2
        else:
            target = 5.5

        if m == 1:
            pos = {next(iter(comp)): np.array([0.0, 0.0])}
        elif m == 2:
            nodes = list(comp)
            pos = {nodes[0]: np.array([-0.7, 0.0]), nodes[1]: np.array([0.7, 0.0])}
        elif m <= 5:
            pos = nx.circular_layout(sub, scale=1.0)
        else:
            try:
                pos = nx.kamada_kawai_layout(sub, weight=None, scale=1.0)
            except Exception:
                pos = nx.spring_layout(sub, seed=42, k=2.4, iterations=500, weight=None)
            pos = nx.spring_layout(sub, pos=pos, seed=42, k=2.8, iterations=250, weight=None)

        arr = np.vstack([pos[n] for n in pos])
        arr = arr - arr.mean(axis=0)
        span = np.ptp(arr, axis=0)
        scale = target / np.maximum(span.max(), 1e-6)
        arr = arr * scale

        width = float(np.ptp(arr[:, 0])) + 1.8
        height = float(np.ptp(arr[:, 1])) + 1.6
        if x_cursor + width > max_row_width and x_cursor > 0:
            x_cursor = 0.0
            y_cursor -= row_height + 0.8
            row_height = 0.0

        cx = x_cursor + width / 2.0
        cy = y_cursor - height / 2.0
        for j, n in enumerate(pos):
            layout_pos[n] = arr[j] + np.array([cx, cy])

        x_cursor += width + 0.6
        row_height = max(row_height, height)

    all_xy = np.vstack(list(layout_pos.values()))
    center = all_xy.mean(axis=0)
    for n in layout_pos:
        layout_pos[n] = layout_pos[n] - center
    return layout_pos


def main() -> None:
    print(f"=== Replot Fig.2D at |ρ|≥{RHO_THR} ===")
    km = pd.read_csv(TAB / "Table_Fig2D_KOmetab_edges_noPI.tsv", sep="\t")
    nodes = pd.read_csv(
        ROOT
        / "01this_study_dataset/sgb/metab_corr/paper_result/13_RobustDA_KO_metab_PI_network/tables/network_nodes.tsv",
        sep="\t",
    )
    enrich = {r["id"]: (r["type"], r.get("enrich_direction", ""), r.get("label", r["id"])) for _, r in nodes.iterrows()}

    Gfull = nx.Graph()
    for _, r in km.iterrows():
        Gfull.add_edge(r["from"], r["to"], rho=float(r["rho"]))
    deg_full = dict(Gfull.degree())

    display = km[km["rho"].abs() >= RHO_THR].copy().sort_values("rho", key=lambda s: s.abs(), ascending=False)
    Gd = nx.Graph()
    for _, r in display.iterrows():
        Gd.add_edge(r["from"], r["to"], rho=float(r["rho"]))
    print(f"  display edges={Gd.number_of_edges()} nodes={Gd.number_of_nodes()} components={nx.number_connected_components(Gd)}")

    layout_pos = layout_by_components(Gd)
    focus = {"K02416", "K02404", "metab_27961"} & set(Gd.nodes())

    fig, ax = plt.subplots(figsize=(12.0, 9.0))
    ax.set_facecolor("white")

    # edges under nodes
    for a, b, d in sorted(Gd.edges(data=True), key=lambda x: abs(x[2]["rho"])):
        r = float(d["rho"])
        col = COL["pos"] if r > 0 else COL["neg"]
        lw = 1.6 + 3.0 * (abs(r) - RHO_THR) / max(1e-6, 1.0 - RHO_THR)
        ax.plot(
            [layout_pos[a][0], layout_pos[b][0]],
            [layout_pos[a][1], layout_pos[b][1]],
            color=col,
            alpha=0.72,
            lw=lw,
            solid_capstyle="round",
            zorder=1,
        )

    # nodes
    for n in Gd.nodes():
        is_ko = str(n).startswith("K")
        c = COL["KO"] if is_ko else COL["Metabolite"]
        d0 = deg_full.get(n, Gd.degree(n))
        size = 140 + 6.0 * min(d0, 40)
        if n in focus:
            size = max(size, 280)
            ec, lw = "#111111", 1.35
        else:
            ec, lw = "#222222", 0.6
        ax.scatter(
            layout_pos[n][0],
            layout_pos[n][1],
            s=size,
            c=c,
            edgecolors=ec,
            linewidths=lw,
            zorder=3,
            alpha=1.0,
        )

    # radial seed offsets by local angle, then adjust_text
    texts = []
    label_rows = []
    nodes_list = list(Gd.nodes())
    for i, n in enumerate(nodes_list):
        lab = short_label(n, enrich, max_len=18)
        is_ko = str(n).startswith("K")
        fs = 9.5 if n in focus else 8.2
        x0, y0 = layout_pos[n]
        ang = 2 * np.pi * i / max(len(nodes_list), 1)
        # push labels outward from component centroid
        rad = 0.55 if Gd.degree(n) >= 3 else 0.42
        tx = x0 + rad * np.cos(ang)
        ty = y0 + rad * np.sin(ang)
        txt = ax.text(
            tx,
            ty,
            lab,
            fontsize=fs,
            ha="center",
            va="center",
            color="#111111",
            zorder=5,
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#CCCCCC", lw=0.45, alpha=0.95),
        )
        texts.append(txt)
        label_rows.append(
            {
                "id": n,
                "label": lab,
                "type": "KO" if is_ko else "Metabolite",
                "x": float(x0),
                "y": float(y0),
                "degree_full_network": int(deg_full.get(n, 0)),
                "degree_display": int(Gd.degree(n)),
            }
        )

    x = [layout_pos[n][0] for n in Gd.nodes()]
    y = [layout_pos[n][1] for n in Gd.nodes()]
    adjust_text(
        texts,
        x=x,
        y=y,
        ax=ax,
        expand_points=(1.8, 1.9),
        expand_text=(1.35, 1.45),
        force_points=(0.9, 1.0),
        force_text=(1.3, 1.5),
        lim=800,
        arrowprops=dict(arrowstyle="-", color="#777777", lw=0.5, shrinkA=2, shrinkB=4),
    )

    ax.set_axis_off()
    pad = 0.85
    ax.set_xlim(min(x) - pad, max(x) + pad)
    ax.set_ylim(min(y) - pad, max(y) + pad)

    ax.set_title(
        "Robust_DA KO–metabolite associations (no Plaque Index)\n"
        f"Display threshold |ρ|≥{RHO_THR:.2f} "
        f"(n_edges={Gd.number_of_edges()}, n_nodes={Gd.number_of_nodes()}; "
        f"full network at |ρ|≥0.80 has {len(km)} edges)",
        fontsize=12,
        pad=12,
    )
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COL["KO"], markeredgecolor="#222", markersize=10, label="KO"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=COL["Metabolite"], markeredgecolor="#222", markersize=10, label="Metabolite"),
        Line2D([0], [0], color=COL["pos"], lw=2.2, label="Positive ρ"),
        Line2D([0], [0], color=COL["neg"], lw=2.2, label="Negative ρ"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=9.5, loc="lower left", bbox_to_anchor=(0.0, -0.03), ncol=2)

    dual_write(display, "Fig2D_KOmetab_network_display_edges")
    dual_write(pd.DataFrame(label_rows), "Fig2D_KOmetab_network_layout")

    fig.tight_layout()
    stem = "Fig2D_KOmetab_network_noPI"
    for ext in ("pdf", "svg", "png", "jpg"):
        fig.savefig(FIG / f"{stem}.{ext}")
    plt.close(fig)
    print(f"  updated {stem}.*  (|ρ|≥{RHO_THR}, {Gd.number_of_edges()} edges / {Gd.number_of_nodes()} nodes)")
    print("DONE")


if __name__ == "__main__":
    main()
