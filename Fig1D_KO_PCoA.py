"""Packaged for revise/; source: project/01this_study_dataset/M2-humann/v1/Figure2A/M2_PanelA.py
Figures: Fig1D
Note: KO Bray-Curtis PCoA
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
"""
M2_PanelA.py — Figure 2 Panel A: Functional Alpha Diversity & Functional Redundancy
Project:  P2026-OralMB
Author:   Bioinformatics Analyst (A)

Usage:
    python3 M2_PanelA.py \\
        --humann   this_study_humann_input.tsv \\
        --group    this_study_group.tsv \\
        --species  this_study_s_relative_pop.xls \\
        --output   ./Fig2_PanelA_out \\
        --seed     42 \\
        --dpi      300 \\
        --n_boot   5000

Output files:
    Fig2_PanelA.pdf / .svg / .png
    Table_S5_normalization_stats.csv
    Table_PanelA_alpha_diversity.csv
    analysis_log_PanelA.txt
"""

import argparse
import logging
import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from scipy.stats import shapiro, mannwhitneyu, ttest_ind

warnings.filterwarnings("ignore")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)

NORMAL_COLOR = "#4DBBD5"
POOR_COLOR   = "#E64B35"


# ══════════════════════════════════════════════════════════════════════════════
# 1. DATA LOADING & CPM NORMALISATION
# ══════════════════════════════════════════════════════════════════════════════

def load_and_normalize(humann_path: str, group_path: str):
    """
    Load HUMAnN3 RPK table and convert to CPM.

    CPM_ij = RPK_ij / total_RPK_j × 1,000,000
    total_RPK_j = UNMAPPED_j + UNGROUPED_j + Σ(all unstratified KO RPK_j)

    Returns
    -------
    cpm_full : DataFrame, 3836 KO × 20 samples, CPM values
    unstrat  : DataFrame, same rows but retains ko_id / ko_ann columns
    normal_samples, poor_samples, all_samps : lists
    mapping_stats : DataFrame, per-sample normalisation statistics
    """
    logging.info("[1/5] Loading HUMAnN3 data ...")
    raw   = pd.read_csv(humann_path, sep="\t")
    grp   = pd.read_csv(group_path, sep="\t", header=None,
                        names=["Sample", "Group"])
    grp.set_index("Sample", inplace=True)

    normal_samples = grp[grp["Group"] == "Normal"].index.tolist()
    poor_samples   = grp[grp["Group"] == "PoorPlaque"].index.tolist()
    all_samps      = normal_samples + poor_samples

    # Verify all expected samples are present
    missing = [s for s in all_samps if s not in raw.columns]
    if missing:
        logging.error(f"Missing samples in HUMAnN3 table: {missing}")
        sys.exit(1)

    feat = raw["Features"]

    # Separate layers
    unmapped_row  = raw[feat == "UNMAPPED"].iloc[0]
    ungrouped_row = raw[feat == "UNGROUPED"].iloc[0]
    unstrat = raw[
        ~feat.str.contains("|", regex=False) &
        ~feat.isin(["UNMAPPED", "UNGROUPED"])
    ].copy()

    logging.info(f"  Unstratified KOs detected: {len(unstrat)}")

    # CPM normalisation
    total_rpk = (
        unmapped_row[all_samps].values.astype(float) +
        ungrouped_row[all_samps].values.astype(float) +
        unstrat[all_samps].values.astype(float).sum(axis=0)
    )

    cpm_vals = unstrat[all_samps].values.astype(float) / total_rpk * 1e6
    cpm_full = pd.DataFrame(cpm_vals, columns=all_samps)

    # Parse KO ID and annotation
    def _parse(row):
        parts = row.split(":", 1)
        return parts[0].strip(), parts[1].strip() if len(parts) > 1 else ""

    ko_ids, ko_anns = zip(*unstrat["Features"].apply(_parse))
    unstrat = unstrat.copy()
    unstrat["ko_id"]  = ko_ids
    unstrat["ko_ann"] = ko_anns
    cpm_full.index    = unstrat["ko_id"].values

    # Normalisation statistics
    mapped_rpk   = unstrat[all_samps].values.astype(float).sum(axis=0)
    mapping_stats = pd.DataFrame({
        "Sample":           all_samps,
        "Group":            [grp.loc[s, "Group"] for s in all_samps],
        "UNMAPPED_RPK":     unmapped_row[all_samps].values.astype(float),
        "UNGROUPED_RPK":    ungrouped_row[all_samps].values.astype(float),
        "MappedKO_RPK":     mapped_rpk,
        "Total_RPK":        total_rpk,
        "Mapped_fraction":  mapped_rpk / total_rpk,
    })

    logging.info(
        f"  Mapped fraction range: "
        f"{mapping_stats['Mapped_fraction'].min():.4f} – "
        f"{mapping_stats['Mapped_fraction'].max():.4f}"
    )
    logging.info(
        f"  Total RPK range: "
        f"{total_rpk.min():.0f} – {total_rpk.max():.0f}"
    )

    return cpm_full, unstrat, normal_samples, poor_samples, all_samps, mapping_stats


# ══════════════════════════════════════════════════════════════════════════════
# 2. DIVERSITY METRIC COMPUTATION
# ══════════════════════════════════════════════════════════════════════════════

def compute_shannon(col_vals: np.ndarray) -> float:
    """H' = -Σ p_i ln p_i  (zero-abundance KOs excluded)."""
    v = col_vals[col_vals > 0].astype(float)
    if v.sum() == 0:
        return np.nan
    p = v / v.sum()
    return float(-np.sum(p * np.log(p)))


def compute_diversity_metrics(cpm_full: pd.DataFrame,
                              all_samps: list,
                              sp_richness: dict) -> pd.DataFrame:
    """
    Compute three functional alpha diversity metrics per sample.

    Parameters
    ----------
    cpm_full    : CPM matrix (KO × samples), full 3836-KO set
    all_samps   : ordered list of sample names
    sp_richness : dict {sample: observed_species_count} from Figure 1

    Returns
    -------
    alpha_df : DataFrame with one row per sample
    """
    logging.info("[2/5] Computing diversity metrics ...")

    rows = []
    for s in all_samps:
        vals        = cpm_full[s].values.astype(float)
        ko_rich     = int((vals > 0).sum())
        shannon_fd  = compute_shannon(vals)
        sp_rich     = sp_richness.get(s, np.nan)
        ratio       = sp_rich / ko_rich if ko_rich > 0 else np.nan
        rows.append({
            "Sample":       s,
            "KO_Richness":  ko_rich,
            "Shannon_FD":   shannon_fd,
            "Sp_KO_Ratio":  ratio,
            "Sp_Richness":  sp_rich,
        })

    alpha_df = pd.DataFrame(rows)
    logging.info(f"  KO Richness: {alpha_df['KO_Richness'].describe().to_dict()}")
    return alpha_df


# ══════════════════════════════════════════════════════════════════════════════
# 3. STATISTICAL TESTING
# ══════════════════════════════════════════════════════════════════════════════

def cohens_d_bootstrap(a: np.ndarray, b: np.ndarray,
                       n_boot: int = 5000, seed: int = 42) -> tuple:
    """
    Cohen's d (pooled SD, denominator n1+n2-2) with 95% bootstrap CI.

    Returns
    -------
    d, ci_lo, ci_hi : float
    """
    n1, n2 = len(a), len(b)
    pool   = np.sqrt(
        ((n1 - 1) * np.var(a, ddof=1) + (n2 - 1) * np.var(b, ddof=1)) /
        (n1 + n2 - 2)
    )
    d = abs(np.mean(a) - np.mean(b)) / pool if pool > 0 else 0.0

    rng = np.random.default_rng(seed)
    ds  = []
    for _ in range(n_boot):
        ra = rng.choice(a, n1, replace=True)
        rb = rng.choice(b, n2, replace=True)
        pl = np.sqrt(
            ((n1 - 1) * np.var(ra, ddof=1) + (n2 - 1) * np.var(rb, ddof=1)) /
            (n1 + n2 - 2)
        )
        ds.append(abs(np.mean(ra) - np.mean(rb)) / pl if pl > 0 else 0.0)

    ci = np.percentile(ds, [2.5, 97.5])
    return d, float(ci[0]), float(ci[1])


def run_stats(alpha_df: pd.DataFrame,
              normal_samples: list,
              poor_samples: list,
              n_boot: int = 5000,
              seed: int = 42) -> dict:
    """
    For each of the three metrics, run:
      - Shapiro-Wilk normality test (both groups)
      - Welch t-test (if both normal) or Mann-Whitney U (otherwise)
      - Cohen's d + 95% bootstrap CI

    Returns
    -------
    stats : dict[metric_key -> result_dict]
    """
    logging.info("[3/5] Running statistical tests ...")

    metrics = {
        "KO_Richness": "KO Richness",
        "Shannon_FD":  "Shannon Functional Diversity (H')",
        "Sp_KO_Ratio": "Species/KO Richness Ratio",
    }
    stats = {}

    for col, label in metrics.items():
        nv = alpha_df.loc[alpha_df["Sample"].isin(normal_samples), col].values.astype(float)
        pv = alpha_df.loc[alpha_df["Sample"].isin(poor_samples),   col].values.astype(float)

        sw_n = shapiro(nv).pvalue
        sw_p = shapiro(pv).pvalue

        if sw_n > 0.05 and sw_p > 0.05:
            _, pval = ttest_ind(nv, pv, equal_var=False)
            test_name = "Welch's t-test"
        else:
            _, pval = mannwhitneyu(pv, nv, alternative="two-sided")
            test_name = "Mann-Whitney U"

        d, ci_lo, ci_hi = cohens_d_bootstrap(nv, pv, n_boot=n_boot, seed=seed)

        stats[col] = dict(
            label      = label,
            nv         = nv,
            pv         = pv,
            n_mean     = float(np.mean(nv)),
            n_sd       = float(np.std(nv, ddof=1)),
            p_mean     = float(np.mean(pv)),
            p_sd       = float(np.std(pv, ddof=1)),
            sw_n       = sw_n,
            sw_p       = sw_p,
            test       = test_name,
            pval       = float(pval),
            d          = d,
            ci_lo      = ci_lo,
            ci_hi      = ci_hi,
        )

        logging.info(
            f"  {col}: Normal {np.mean(nv):.3f}±{np.std(nv,ddof=1):.3f} | "
            f"Poor {np.mean(pv):.3f}±{np.std(pv,ddof=1):.3f} | "
            f"{test_name} p={pval:.3e} | d={d:.3f} [{ci_lo:.3f}–{ci_hi:.3f}]"
        )

    return stats


# ══════════════════════════════════════════════════════════════════════════════
# 4. FIGURE
# ══════════════════════════════════════════════════════════════════════════════

def _pval_stars(p: float) -> str:
    if   p < 0.001: return "***"
    elif p < 0.01:  return "**"
    elif p < 0.05:  return "*"
    return "ns"


def draw_panel_A(alpha_df: pd.DataFrame,
                 stats: dict,
                 normal_samples: list,
                 poor_samples: list,
                 seed: int = 42,
                 dpi: int = 300,
                 out_base: str = "Fig2_PanelA") -> None:
    """
    Draw Panel A: three side-by-side box + jitter plots.
    Saves .pdf, .svg, and .png.
    """
    logging.info("[4/5] Drawing figure ...")

    plt.rcParams.update({
        "font.family":       "DejaVu Sans",
        "axes.labelsize":    10,
        "xtick.labelsize":   9,
        "ytick.labelsize":   9,
        "axes.linewidth":    0.9,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "pdf.fonttype":      42,
        "svg.fonttype":      "none",
    })
    sns.set_style("ticks")

    fig, axes = plt.subplots(1, 3, figsize=(7.09, 3.70))
    fig.subplots_adjust(
        left=0.09, right=0.97,
        bottom=0.18, top=0.88,
        wspace=0.52,
    )

    metric_keys = ["KO_Richness", "Shannon_FD", "Sp_KO_Ratio"]
    ylabels = [
        "KO Richness",
        "Shannon Functional Diversity (H')",
        "Species / KO Richness Ratio",
    ]
    titles = [
        "KO Richness",
        "Shannon\nFunctional H'",
        "Species/KO\nRichness Ratio",
    ]

    for i, (col, ylabel, title) in enumerate(zip(metric_keys, ylabels, titles)):
        ax = axes[i]
        st = stats[col]
        nv, pv = st["nv"], st["pv"]

        rng = np.random.default_rng(seed + i)

        # ── Box plots ──────────────────────────────────────────────────────
        for grp_vals, pos, col_hex in [
            (nv, 1, NORMAL_COLOR),
            (pv, 2, POOR_COLOR),
        ]:
            bp = ax.boxplot(
                [grp_vals],
                positions=[pos],
                widths=0.42,
                patch_artist=True,
                medianprops=dict(color="white", linewidth=1.5),
                whiskerprops=dict(linewidth=0.9, color="#444444"),
                capprops=dict(linewidth=0.9, color="#444444"),
                flierprops=dict(marker=""),
                boxprops=dict(linewidth=0.9),
                zorder=3,
            )
            bp["boxes"][0].set_facecolor(col_hex)
            bp["boxes"][0].set_alpha(0.75)

            # Jitter
            jx = rng.normal(pos, 0.07, size=len(grp_vals))
            ax.scatter(
                jx, grp_vals,
                color=col_hex, alpha=0.82, s=30,
                zorder=5, edgecolors="white", linewidths=0.4,
            )

        # ── Significance annotation ────────────────────────────────────────
        y_max  = max(float(nv.max()), float(pv.max()))
        y_min  = min(float(nv.min()), float(pv.min()))
        y_rng  = max(y_max - y_min, abs(y_max) * 0.05, 1e-9)

        y_bar  = y_max + y_rng * 0.18
        y_star = y_max + y_rng * 0.25
        y_d    = y_max + y_rng * 0.40

        # Bracket
        ax.plot([1, 2], [y_bar, y_bar], "k-", lw=0.9, zorder=6)
        ax.plot([1, 1], [y_max + y_rng * 0.02, y_bar], "k-", lw=0.9, zorder=6)
        ax.plot([2, 2], [float(pv.max()) + y_rng * 0.02, y_bar], "k-", lw=0.9, zorder=6)

        stars = _pval_stars(st["pval"])
        p_fmt = (
            f"p={st['pval']:.2e}" if st["pval"] < 0.001
            else f"p={st['pval']:.4f}"
        )
        ax.text(1.5, y_star, f"{stars}  {p_fmt}",
                ha="center", va="bottom", fontsize=7.5, fontweight="bold")
        ax.text(
            1.5, y_d,
            f"d={st['d']:.2f} [{st['ci_lo']:.2f}–{st['ci_hi']:.2f}]",
            ha="center", va="bottom", fontsize=7.0, color="#444444",
        )

        # Mean ± SD below x-axis
        fmt = ".1f" if col == "KO_Richness" else ".3f"
        ax.text(
            1, y_min - y_rng * 0.13,
            f"{st['n_mean']:{fmt}}±{st['n_sd']:{fmt}}",
            ha="center", fontsize=6.5, color=NORMAL_COLOR, style="italic",
        )
        ax.text(
            2, y_min - y_rng * 0.13,
            f"{st['p_mean']:{fmt}}±{st['p_sd']:{fmt}}",
            ha="center", fontsize=6.5, color=POOR_COLOR, style="italic",
        )

        # Test label (lower right)
        ax.text(
            0.97, 0.02, st["test"],
            transform=ax.transAxes,
            fontsize=6.0, ha="right", va="bottom", color="#666666",
        )

        # Axes formatting
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Normal", "PoorPlaque"], fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9.5, labelpad=4)
        ax.set_title(title, fontsize=9.0, pad=3)
        ax.set_xlim(0.40, 2.60)
        ax.set_ylim(
            bottom=(max(y_min - y_rng * 0.18, 0) if col != "Sp_KO_Ratio"
                    else y_min - y_rng * 0.18),
            top=y_d + y_rng * 0.25,
        )
        sns.despine(ax=ax)

    # Panel label
    axes[0].text(
        -0.24, 1.08, "A",
        transform=axes[0].transAxes,
        fontsize=14, fontweight="bold", va="top",
    )

    # Footer note
    fig.text(
        0.50, 0.01,
        ("Species/KO Richness Ratio: functional redundancy proxy (Louca et al. 2018). "
         "Observed Species from Fig. 1 (MetaPhlAn v4). "
         "Bootstrap 95% CI for Cohen's d, 5,000 iter., seed=42."),
        ha="center", va="bottom", fontsize=5.8, color="#555555", style="italic",
    )

    # ── Save ──────────────────────────────────────────────────────────────
    for ext, kw in [
        (".pdf", dict(format="pdf")),
        (".svg", dict(format="svg")),
        (".png", dict(format="png", dpi=dpi)),
    ]:
        path = out_base + ext
        fig.savefig(path, bbox_inches="tight", **kw)
        sz = os.path.getsize(path) // 1024
        logging.info(f"  Saved {os.path.basename(path)} ({sz} KB)")

    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# 5. OUTPUT TABLES & LOG
# ══════════════════════════════════════════════════════════════════════════════

def save_tables(alpha_df: pd.DataFrame,
                stats: dict,
                mapping_stats: pd.DataFrame,
                normal_samples: list,
                poor_samples: list,
                tables_dir: str) -> None:
    """Write Table_S5 and Table_PanelA."""
    os.makedirs(tables_dir, exist_ok=True)

    # Table S5 – normalisation stats
    mapping_stats.to_csv(
        os.path.join(tables_dir, "Table_S5_normalization_stats.csv"),
        index=False,
    )

    # Table PanelA – per-sample diversity values + group summaries
    # Add Group column
    grp_col = (
        ["Normal"] * len(normal_samples) +
        ["PoorPlaque"] * len(poor_samples)
    )
    alpha_out = alpha_df.copy()
    alpha_out.insert(1, "Group", grp_col)
    alpha_out.to_csv(
        os.path.join(tables_dir, "Table_PanelA_alpha_diversity.csv"),
        index=False,
    )

    logging.info(f"  Tables written to {tables_dir}/")


def write_log(stats: dict,
              mapping_stats: pd.DataFrame,
              n_total_ko: int,
              args,
              out_dir: str) -> None:
    """Write plain-text analysis log."""
    lines = [
        "# Figure 2 Panel A — Analysis Log",
        f"Script:      M2_PanelA.py",
        f"HUMAnN3:     {args.humann}",
        f"Group:       {args.group}",
        f"Species:     {args.species}",
        f"Seed:        {args.seed}",
        f"DPI:         {args.dpi}",
        f"n_boot:      {args.n_boot}",
        "",
        "## Normalisation",
        "Method:      CPM  (total_RPK = UNMAPPED + UNGROUPED + Σ KO RPK)",
        f"Total KOs:   {n_total_ko}",
        f"Mapped frac: {mapping_stats['Mapped_fraction'].min():.4f} – "
        f"{mapping_stats['Mapped_fraction'].max():.4f}",
        "",
        "## Note on pre-filtering",
        "No prevalence/mean-CPM pre-filter applied for diversity metrics.",
        "Full 3,836-KO matrix used (filtering would underestimate richness).",
        "",
        "## Statistical Results",
    ]
    for col, st in stats.items():
        lines += [
            f"\n[{col}]",
            f"  Normal:     {st['n_mean']:.4f} ± {st['n_sd']:.4f}",
            f"  PoorPlaque: {st['p_mean']:.4f} ± {st['p_sd']:.4f}",
            f"  Shapiro-Wilk (Normal/Poor): p={st['sw_n']:.4f} / p={st['sw_p']:.4f}",
            f"  Test:       {st['test']}",
            f"  p-value:    {st['pval']:.4e}",
            f"  Cohen's d:  {st['d']:.4f}  95%CI [{st['ci_lo']:.4f}, {st['ci_hi']:.4f}]",
        ]
    log_path = os.path.join(out_dir, "analysis_log_PanelA.txt")
    with open(log_path, "w") as f:
        f.write("\n".join(lines))
    logging.info(f"  Log written: {log_path}")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(
        description="Figure 2 Panel A: Functional Alpha Diversity — P2026-OralMB",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--humann", required=True,
        help="HUMAnN3 KO abundance table (TSV, RPK units)",
    )
    p.add_argument(
        "--group", required=True,
        help="Sample group metadata (TSV, no header; col1=sample, col2=group)",
    )
    p.add_argument(
        "--species", required=True,
        help="MetaPhlAn species relative abundance table (TSV) from Figure 1",
    )
    p.add_argument(
        "--output", default="./Fig2_PanelA_out",
        help="Output directory",
    )
    p.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for bootstrap CI",
    )
    p.add_argument(
        "--dpi", type=int, default=300,
        help="PNG output resolution (DPI)",
    )
    p.add_argument(
        "--n_boot", type=int, default=5000,
        help="Bootstrap iterations for Cohen's d CI",
    )
    return p.parse_args()


def main():
    args = parse_args()

    os.makedirs(args.output, exist_ok=True)
    tables_dir = os.path.join(args.output, "tables")
    panels_dir = os.path.join(args.output, "panels")
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(panels_dir, exist_ok=True)

    # ── 1. Load & normalise ───────────────────────────────────────────────
    (cpm_full, unstrat,
     normal_samples, poor_samples, all_samps,
     mapping_stats) = load_and_normalize(args.humann, args.group)

    # ── 2. Species richness from Figure 1 ─────────────────────────────────
    logging.info("[2/5] Loading Figure 1 species data ...")
    abund = pd.read_csv(args.species, sep="\t", index_col=0)
    missing_sp = [s for s in all_samps if s not in abund.columns]
    if missing_sp:
        logging.warning(f"  Samples not found in species table: {missing_sp}")
    sp_richness = {
        s: int((abund[s] > 0).sum())
        for s in all_samps if s in abund.columns
    }
    logging.info(
        f"  Species richness range: "
        f"{min(sp_richness.values())} – {max(sp_richness.values())}"
    )

    # ── 3. Diversity metrics ───────────────────────────────────────────────
    alpha_df = compute_diversity_metrics(cpm_full, all_samps, sp_richness)

    # ── 4. Statistics ──────────────────────────────────────────────────────
    stats = run_stats(
        alpha_df, normal_samples, poor_samples,
        n_boot=args.n_boot, seed=args.seed,
    )

    # ── 5. Figure ──────────────────────────────────────────────────────────
    out_base = os.path.join(panels_dir, "Fig2_PanelA")
    draw_panel_A(
        alpha_df, stats,
        normal_samples, poor_samples,
        seed=args.seed, dpi=args.dpi,
        out_base=out_base,
    )

    # ── 6. Tables & log ────────────────────────────────────────────────────
    save_tables(alpha_df, stats, mapping_stats,
                normal_samples, poor_samples, tables_dir)
    write_log(stats, mapping_stats, len(cpm_full), args, args.output)

    # ── 7. Terminal summary ────────────────────────────────────────────────
    logging.info("[5/5] Done. Summary:")
    for col, st in stats.items():
        logging.info(
            f"  {col:20s}  N={st['n_mean']:.3f}±{st['n_sd']:.3f}  "
            f"P={st['p_mean']:.3f}±{st['p_sd']:.3f}  "
            f"p={st['pval']:.2e}  d={st['d']:.3f}"
        )
    logging.info(f"  Output: {args.output}/")


if __name__ == "__main__":
    main()
