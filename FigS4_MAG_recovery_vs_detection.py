"""Packaged for revise/; source: project/review/R1.6/scripts/R1.6A_fig_redesign.py
Figures: FigS4
Note: recovery vs abundance
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
"""Fig_R1.6A — compact single-panel: larger fonts, tight limits, less whitespace."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(str(PKG_ROOT) + "/review/R1.6")
PLOTDATA = ROOT / "tables/Fig_R1.6A_recovery_vs_abundance_plotdata.tsv"
OUT_F = ROOT / "figures" / "Fig_R1.6A_recovery_vs_detection"
OUT_T = ROOT / "tables"

C_NORMAL = "#2F6B8A"
C_POOR = "#C45C26"
C_DIAG = "#6E747C"
C_BAND = "#E4EFE8"
C_BAND2 = "#F5F2EE"

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.weight": "normal",
        "axes.labelweight": "normal",
        "axes.titleweight": "normal",
        "axes.labelsize": 11.5,
        "xtick.labelsize": 10.5,
        "ytick.labelsize": 10.5,
        "legend.fontsize": 9.5,
        "axes.titlesize": 12.5,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.85,
        "xtick.major.width": 0.75,
        "ytick.major.width": 0.75,
        "xtick.major.size": 3.8,
        "ytick.major.size": 3.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
    }
)


def write_tsv_csv(df: pd.DataFrame, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(str(stem) + ".tsv", sep="\t", index=False)
    df.to_csv(str(stem) + ".csv", index=False)


def save_all(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "svg", "png", "jpg"):
        out = Path(str(stem) + f".{ext}")
        kw = dict(bbox_inches="tight", pad_inches=0.015, facecolor="white")
        if ext in ("png", "jpg"):
            fig.savefig(out, dpi=600, **kw)
        else:
            fig.savefig(out, **kw)
    print("Saved", stem)


def main():
    df = pd.read_csv(PLOTDATA, sep="\t")
    if "mean_abund_other_group" not in df.columns:
        df["mean_abund_other_group"] = df["mean_abund_opposite_group"]
    df["x"] = np.log10(df["mean_abund_other_group"] + 1e-3)
    df["y"] = np.log10(df["mean_abund_recovery_group"] + 1e-3)
    df["higher_in_recovery"] = df["mean_abund_recovery_group"] > df["mean_abund_other_group"]
    df["is_one_only"] = df["concordance_class"] != "recovered_one_group_but_detected_both"

    n_both = int((~df["is_one_only"]).sum())
    n_one = int(df["is_one_only"].sum())
    n_tot = len(df)
    n_normal = int((df["recovery_group"] == "Normal").sum())
    n_poor = int((df["recovery_group"] == "PoorPlaque").sum())
    pct_higher = 100.0 * df["higher_in_recovery"].mean()
    med_l2 = float(df["log2_mean_recovery_over_opposite"].median())

    write_tsv_csv(
        pd.DataFrame(
            [
                {"metric": "n_MAG", "value": n_tot},
                {"metric": "n_detected_both", "value": n_both},
                {"metric": "n_detected_one_only", "value": n_one},
                {"metric": "pct_higher_in_recovery", "value": round(pct_higher, 1)},
                {"metric": "median_log2_ratio", "value": round(med_l2, 3)},
            ]
        ),
        OUT_T / "Fig_R1.6A_redesign_panel_summary",
    )
    write_tsv_csv(df, OUT_T / "Fig_R1.6A_recovery_vs_abundance_plotdata")

    both = df[~df["is_one_only"]]
    one = df[df["is_one_only"]]

    # Tight limits from the main cloud (both-group MAGs)
    x_lo, x_hi = both["x"].quantile([0.0, 1.0])
    y_lo, y_hi = both["y"].quantile([0.0, 1.0])
    span = max(x_hi - x_lo, y_hi - y_lo)
    pad = 0.12 * span
    lim_lo = float(min(x_lo, y_lo) - pad)
    lim_hi = float(max(x_hi, y_hi) + pad)

    fig, ax = plt.subplots(figsize=(4.0, 4.0), facecolor="white")
    fig.subplots_adjust(left=0.155, right=0.98, top=0.92, bottom=0.135)

    xx = np.linspace(lim_lo, lim_hi, 400)
    ax.fill_between(xx, xx, lim_hi, color=C_BAND, alpha=0.95, zorder=0, lw=0)
    ax.fill_between(xx, lim_lo, xx, color=C_BAND2, alpha=0.95, zorder=0, lw=0)
    ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], ls=(0, (3.5, 2)), lw=1.05, color=C_DIAG, zorder=1)

    for g, color in [("Normal", C_NORMAL), ("PoorPlaque", C_POOR)]:
        sub = both[both["recovery_group"] == g]
        ax.scatter(
            sub["x"],
            sub["y"],
            s=42,
            c=color,
            marker="o",
            alpha=0.88,
            edgecolors="white",
            linewidths=0.4,
            zorder=3,
        )

    # Place one-group-only MAG at left margin of zoomed view (true other-group abund = 0)
    if len(one):
        for _, r in one.iterrows():
            y_true = float(r["y"])
            x_plot = lim_lo + 0.04 * (lim_hi - lim_lo)
            y_plot = float(np.clip(y_true, lim_lo + 0.05 * (lim_hi - lim_lo), lim_hi - 0.05 * (lim_hi - lim_lo)))
            ax.scatter(
                [x_plot],
                [y_plot],
                s=85,
                c=C_POOR if r["recovery_group"] == "PoorPlaque" else C_NORMAL,
                marker="D",
                alpha=1.0,
                edgecolors="#1a1a1a",
                linewidths=0.7,
                zorder=6,
            )

    ax.set_xlim(lim_lo, lim_hi)
    ax.set_ylim(lim_lo, lim_hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("log$_{10}$ mean abundance in the other group", labelpad=4)
    ax.set_ylabel("log$_{10}$ mean abundance in the recovery group", labelpad=4)

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=C_NORMAL,
            markersize=8.5,
            label=f"Recovered from Normal (n={n_normal})",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=C_POOR,
            markersize=8.5,
            label=f"Recovered from PoorPlaque (n={n_poor})",
        ),
        Line2D(
            [0],
            [0],
            marker="D",
            color="w",
            markerfacecolor=C_POOR,
            markeredgecolor="#1a1a1a",
            markersize=8,
            label=f"Detected in one group only ({n_one}/{n_tot})",
        ),
    ]
    leg = ax.legend(
        handles=handles,
        loc="lower right",
        frameon=True,
        fancybox=False,
        edgecolor="#D0D0D0",
        framealpha=0.94,
        handletextpad=0.4,
        borderpad=0.5,
        labelspacing=0.4,
        borderaxespad=0.4,
        fontsize=9.5,
    )
    leg.get_frame().set_linewidth(0.55)

    ax.text(
        0.03,
        0.97,
        f"{pct_higher:.0f}% higher in recovery group\nmedian log$_2$ ratio = {med_l2:.2f}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9.5,
        color="#1f1f1f",
        linespacing=1.3,
        bbox=dict(
            boxstyle="round,pad=0.32",
            facecolor="white",
            edgecolor="#D8D8D8",
            linewidth=0.55,
            alpha=0.94,
        ),
    )

    ax.set_title(
        "MAG recovery origin versus MAG abundance",
        fontsize=12.5,
        pad=6,
        loc="left",
        color="#1a1a1a",
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", which="major", pad=2.5)

    save_all(fig, OUT_F)
    plt.close(fig)


if __name__ == "__main__":
    main()
