"""Packaged for revise/; source: project/review/R1.4/scripts/07_replot_FigS_absRho_vs_bootstrap_panelB_only.py
Figures: Fig2D
Note: stability vs |rho| plot
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
"""Replot FigS absRho vs bootstrap: panel B only (binned stability)."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = PKG_ROOT
OUT = ROOT / "review/R1.4/Fig2D_KOmetab_noPI"
TAB = OUT / "tables"
FIG = OUT / "figures"
STEM = "FigS_Fig2D_noPI_absRho_vs_bootstrap_support"


def main() -> None:
    summ = pd.read_csv(TAB / "FigS_Fig2D_noPI_absRho_vs_bootstrap_binned_plotdata.tsv", sep="\t")

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

    # Tall/narrow canvas for vertical panel compositing (avoid bbox_inches='tight'
    # which expands width around long titles and flattens the intended aspect).
    fig, ax = plt.subplots(figsize=(3.4, 5.4))
    x = np.arange(len(summ))
    ax.bar(
        x,
        summ["frac_robust_edge"],
        width=0.55,
        color="#E07A3D",
        alpha=0.85,
        label="Proportion robust",
    )
    ax.set_ylabel("Proportion of robust edges")
    ax.set_ylim(0, 1.08)
    ax.set_xticks(x)
    ax.set_xticklabels(
        summ["abs_rho_bin"],
        fontsize=8,
        rotation=45,
        ha="right",
        rotation_mode="anchor",
    )
    ax.set_xlabel(r"$|\rho|$ bin")
    ax.spines["right"].set_visible(False)

    ax2 = ax.twinx()
    ax2.plot(
        x,
        summ["median_bootstrap_support"],
        "o-",
        color="#2F4B7C",
        lw=1.4,
        ms=5,
        label="Median bootstrap support",
    )
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

    for i, n in enumerate(summ["n_edges"]):
        ax.text(
            i,
            min(0.98, summ["frac_robust_edge"].iloc[i] + 0.04),
            f"n={n}",
            ha="center",
            fontsize=7,
        )

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=6.5, loc="upper left")
    ax.set_title("Stability by |ρ| bin\n(full KO–metabolite network, no PI)", fontsize=9)
    fig.suptitle(
        "Stability of Fig. 2D KO–metabolite edges\nincreases with |ρ| (n = 2,158)",
        fontsize=10,
        y=0.98,
    )
    fig.subplots_adjust(left=0.18, right=0.82, bottom=0.18, top=0.86)

    for ext in ("pdf", "svg", "png", "jpg"):
        fig.savefig(FIG / f"{STEM}.{ext}", dpi=300, bbox_inches=None, pad_inches=0.05)
    plt.close(fig)

    (TAB / f"{STEM}_CAPTION.txt").write_text(
        """Supplementary Figure. Bootstrap stability of Plaque Index–removed Fig. 2D KO–metabolite edges increases with absolute Spearman ρ.

Edges (n = 2,158 significant KO–metabolite associations; |ρ|≥0.80, P<0.05, FDR≤0.05) are stratified by |ρ| bins. Bars: proportion of robust edges (bootstrap support ≥0.8 and leave-one-out support ≥0.8). Blue line/points: median bootstrap support with interquartile range (shaded). Numbers above bars indicate edge counts per bin. Dashed line: support = 0.8. Bootstrap support was estimated by participant-level resampling (B=500, n=20).

Interpretation: stronger associations are substantially more reproducible under resampling, indicating a stable high-|ρ| core despite the small cohort size. This does not remove small-sample limitations or establish causality, but it argues against the network being uniformly noise-driven.
""",
        encoding="utf-8",
    )
    print(f"saved {STEM}.{{pdf,svg,png,jpg}}")


if __name__ == "__main__":
    main()
