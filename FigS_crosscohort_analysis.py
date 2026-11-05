"""Packaged for revise/; source: project/02dataset/analyse/plaque_cross_cohort_analysis.py
Figures: FigS7;FigS8;FigS9
Note: cross-cohort MetaPhlAn-MAG
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
# -*- coding: utf-8 -*-
"""
================================================================================
Adolescent Dental Plaque Metagenome — Cross-cohort Core Species & DA Analysis
================================================================================

Integrated pipeline for:
  M1. Core species identification across 3 independent cohorts
  M2. PRJNA766357 differential abundance analysis (Normal vs Caries)
  M3. PRJNA766357 alpha & beta diversity (skbio)
  M4. All-sample beta diversity & batch effect assessment
  M5. Cross-cohort DA concordance (strict + relaxed thresholds)
  M6. Publication-ready visualisation (9 figures, PDF/SVG/JPG)

Usage:
  python3 plaque_cross_cohort_analysis.py \\
      --abund   s_relative_sort.xls \\
      --group   group.tsv \\
      --da-ts   Table_S_differential_species.csv \\
      --output  ./results \\
      --seed    42 --dpi 300

Dependencies:
  numpy pandas scipy statsmodels scikit-bio seaborn matplotlib
"""

import argparse, os, sys, math, warnings, logging
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist, squareform
from statsmodels.stats.multitest import multipletests

from skbio.stats.ordination import pcoa as sk_pcoa
from skbio.stats.distance import DistanceMatrix, permanova, anosim

import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = 'serif'
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
import seaborn as sns

matplotlib.rcParams.update({
    'font.size': 9, 'axes.titlesize': 9, 'axes.labelsize': 10,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
    'figure.dpi': 300, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'axes.linewidth': 0.8, 'lines.linewidth': 1.0,
})

# ══════════════════════════════════════════════════════════════
#  GLOBAL PALETTE  (consistent with v2.py / Figure 1)
# ══════════════════════════════════════════════════════════════
COL = {
    'Normal': '#2E8B8B', 'PoorPlaque': '#B85450',
    'Caries': '#D4A843', 'Unknown': '#888888',
}
DS_COL = {
    'ThisStudy': '#6A5ACD', 'PRJNA597465': '#888888', 'PRJNA766357': '#D4A843',
}
VOL_ROBUST_UP  = '#B85450'
VOL_ROBUST_DN  = '#2E8B8B'
VOL_WIL_UP     = '#F4A582'
VOL_WIL_DN     = '#92C5C5'
VOL_NS         = '#DDDDDD'


# ══════════════════════════════════════════════════════════════
#  UTILITY FUNCTIONS
# ══════════════════════════════════════════════════════════════
def save_fig(fig, base, dpi=300):
    """Save figure as PDF + SVG + JPG."""
    fig.savefig(f"{base}.pdf", format='pdf', bbox_inches='tight')
    fig.savefig(f"{base}.jpg", format='jpg', dpi=dpi, bbox_inches='tight')
    fig.savefig(f"{base}.svg", format='svg', bbox_inches='tight')
    logging.info(f"  [SAVED] {base} (.pdf/.svg/.jpg)")


def fmt_pval(p):
    return f"p = {p:.2e}" if p < 0.001 else f"p = {p:.4f}"


def calc_shannon(row):
    p = row[row > 0]; p = p / p.sum()
    return -(p * np.log(p)).sum()


def calc_simpson(row):
    p = row[row > 0]; p = p / p.sum()
    return 1 - (p ** 2).sum()


def cov_ellipse(x, y, ax, n_std=1.96, facecolor='none', **kwargs):
    """Draw 95% confidence ellipse."""
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    evals, evecs = np.linalg.eigh(cov)
    order = evals.argsort()[::-1]
    evals, evecs = evals[order], evecs[:, order]
    angle = np.degrees(np.arctan2(evecs[1, 0], evecs[0, 0]))
    width  = 2 * n_std * np.sqrt(max(evals[0], 0))
    height = 2 * n_std * np.sqrt(max(evals[1], 0))
    ax.add_patch(Ellipse(xy=(np.mean(x), np.mean(y)),
                         width=width, height=height, angle=angle,
                         facecolor=facecolor, **kwargs))


def clr_mc(abund_f, sn, sp, n_boot=128, seed=42):
    """CLR-Monte-Carlo Wilcoxon for compositional robustness."""
    rng = np.random.default_rng(seed)
    mat = abund_f.T.values
    n_sp = mat.shape[1]
    scaled = mat * 10 + 0.01
    cols = list(abund_f.columns)
    ni = [cols.index(s) for s in sn]
    pi = [cols.index(s) for s in sp]
    p_mat = np.zeros((n_boot, n_sp))
    for b in range(n_boot):
        sim = np.array([rng.dirichlet(scaled[i] + 1e-6) * scaled[i].sum()
                        for i in range(scaled.shape[0])])
        sim[sim == 0] = np.finfo(float).tiny
        lg  = np.log(sim)
        clr = lg - lg.mean(axis=1, keepdims=True)
        for j in range(n_sp):
            _, p = stats.mannwhitneyu(clr[ni, j], clr[pi, j],
                                      alternative='two-sided')
            p_mat[b, j] = p
    return np.median(p_mat, axis=0)


def boxplot_two_groups(ax, vals_a, vals_b, labels, colors, ylabel, p_val, seed,
                       title=None, title_color='#333'):
    """Draw a styled 2-group boxplot on *ax*."""
    for i, (lab, vals, color) in enumerate(zip(labels, [vals_a, vals_b], colors)):
        ax.boxplot(vals, positions=[i], widths=0.45,
                   patch_artist=True,
                   medianprops=dict(color='white', linewidth=1.5),
                   whiskerprops=dict(color=color, linewidth=0.8),
                   capprops=dict(color=color, linewidth=0.8),
                   flierprops=dict(marker='', markersize=0),
                   boxprops=dict(facecolor=color, alpha=0.75,
                                 linewidth=0.8, edgecolor=color))
        rng = np.random.default_rng(seed)
        jx = rng.uniform(-0.12, 0.12, len(vals))
        ax.scatter(np.full(len(vals), i) + jx, vals,
                   color=color, alpha=0.75, s=22, zorder=5, linewidths=0)
    all_v = np.concatenate([vals_a, vals_b])
    y_max, y_min = all_v.max(), all_v.min()
    y_range = max(y_max - y_min, abs(y_max) * 0.1, 0.001)
    y_line = y_max + y_range * 0.10
    y_text = y_line + y_range * 0.04
    ax.plot([0, 0, 1, 1],
            [y_line - y_range*0.015, y_line, y_line, y_line - y_range*0.015],
            lw=0.8, color='#333')
    ax.text(0.5, y_text, fmt_pval(p_val), ha='center', va='bottom', fontsize=6.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels, fontsize=7, rotation=45, ha='right')
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_xlim(-0.65, 1.65)
    ax.set_ylim(bottom=y_min - y_range*0.12, top=y_text + y_range*0.14)
    sns.despine(ax=ax)
    if title:
        ax.set_title(title, fontsize=6.5, style='italic', pad=3, color=title_color)


# ══════════════════════════════════════════════════════════════
#  MAIN PIPELINE
# ══════════════════════════════════════════════════════════════
def run(args):
    OUT = args.output
    SEED = args.seed
    DPI  = args.dpi
    for sub in ['figures', 'tables']:
        os.makedirs(os.path.join(OUT, sub), exist_ok=True)

    fig_dir = os.path.join(OUT, 'figures')
    tab_dir = os.path.join(OUT, 'tables')

    # ────────────────────────────────────────────────────────
    #  LOAD DATA
    # ────────────────────────────────────────────────────────
    logging.info("=" * 70)
    logging.info("Loading data ...")
    abund = pd.read_csv(args.abund, sep='\t', index_col=0)
    if 'TAX_INFO' in abund.columns:
        abund = abund.drop(columns=['TAX_INFO'])

    group_df = pd.read_csv(args.group, sep='\t')
    group_map   = dict(zip(group_df['Sample'], group_df['Group']))
    dataset_map = dict(zip(group_df['Sample'], group_df['Dataset']))

    da_thisstudy = pd.read_csv(args.da_ts)

    ds_samples = {}
    for ds in ['ThisStudy', 'PRJNA597465', 'PRJNA766357']:
        ds_samples[ds] = [s for s in group_df[group_df['Dataset'] == ds]['Sample']
                          if s in abund.columns]

    ts_normal  = [s for s in ds_samples['ThisStudy']    if group_map[s] == 'Normal']
    ts_poor    = [s for s in ds_samples['ThisStudy']    if group_map[s] == 'PoorPlaque']
    p766_normal = [s for s in ds_samples['PRJNA766357'] if group_map[s] == 'Normal']
    p766_caries = [s for s in ds_samples['PRJNA766357'] if group_map[s] == 'Caries']

    logging.info(f"  Species × Samples: {abund.shape}")
    for ds, ss in ds_samples.items():
        logging.info(f"  {ds}: n={len(ss)}")
    logging.info(f"  ThisStudy   Normal={len(ts_normal)}, PoorPlaque={len(ts_poor)}")
    logging.info(f"  PRJNA766357 Normal={len(p766_normal)}, Caries={len(p766_caries)}")

    # ────────────────────────────────────────────────────────
    # M1. CORE SPECIES
    # ────────────────────────────────────────────────────────
    logging.info("\n" + "=" * 70)
    logging.info("[M1] Core species analysis")

    prev_data = {}
    for ds, ss in ds_samples.items():
        prev_data[ds] = (abund[ss] > 0).sum(axis=1) / len(ss)
    prev_df = pd.DataFrame(prev_data)
    prev_df.columns = ['ThisStudy', 'PRJNA597465', 'PRJNA766357']

    thresholds = [0.3, 0.5, 0.7]
    core_results = {}
    for thr in thresholds:
        n_ds = (prev_df >= thr).sum(axis=1)
        core_results[thr] = {
            'in_all_3': prev_df.index[n_ds == 3].tolist(),
            'in_2':     prev_df.index[n_ds >= 2].tolist(),
            'in_1':     prev_df.index[n_ds >= 1].tolist(),
        }
        logging.info(f"  Prevalence >= {thr*100:.0f}%:  3/3={len(core_results[thr]['in_all_3'])}  "
                     f">=2/3={len(core_results[thr]['in_2'])}  >=1/3={len(core_results[thr]['in_1'])}")

    core70 = core_results[0.7]['in_all_3']
    core50 = core_results[0.5]['in_all_3']
    core30 = core_results[0.3]['in_all_3']

    all_samp = [s for ss in ds_samples.values() for s in ss]
    grand_mean = abund[all_samp].mean(axis=1)

    # Disease impact on core
    core_da_ts = da_thisstudy[da_thisstudy['species'].isin(core50) &
                               (da_thisstudy['Robust_DA'] == True)]
    logging.info(f"  Core50 species that are Robust DA in ThisStudy: {len(core_da_ts)}")

    # ────────────────────────────────────────────────────────
    # M2. PRJNA766357 DIFFERENTIAL ABUNDANCE
    # ────────────────────────────────────────────────────────
    logging.info("\n" + "=" * 70)
    logging.info("[M2] PRJNA766357 differential abundance (Normal vs Caries)")

    abund_766 = abund[p766_normal + p766_caries]
    n_samp766 = len(p766_normal) + len(p766_caries)
    prev_766  = (abund_766 > 0).sum(axis=1) / n_samp766
    mn_766    = abund_766[p766_normal].mean(axis=1)
    mc_766    = abund_766[p766_caries].mean(axis=1)
    keep      = (prev_766 >= 0.30) & ((mn_766 > 0.01) | (mc_766 > 0.01))
    abund_766f = abund_766.loc[keep]
    logging.info(f"  Pre-filter: {abund_766.shape[0]} -> {abund_766f.shape[0]} species")

    species_766 = abund_766f.index.tolist()
    p_wil_766 = np.array([
        stats.mannwhitneyu(abund_766f.loc[s, p766_normal].values,
                           abund_766f.loc[s, p766_caries].values,
                           alternative='two-sided')[1]
        for s in species_766])
    _, fdr_766, _, _ = multipletests(p_wil_766, method='fdr_bh')

    mean_n766 = abund_766f[p766_normal].mean(axis=1).values
    mean_c766 = abund_766f[p766_caries].mean(axis=1).values
    sd_n766   = abund_766f[p766_normal].std(axis=1).values
    sd_c766   = abund_766f[p766_caries].std(axis=1).values
    l2fc_766  = np.log2((mean_c766 + 0.01) / (mean_n766 + 0.01))

    logging.info("  Running CLR-MC-Wilcoxon (128 draws) ...")
    p_clrmc_766 = clr_mc(abund_766f, p766_normal, p766_caries, n_boot=128, seed=SEED)
    _, fdr_clrmc_766, _, _ = multipletests(p_clrmc_766, method='fdr_bh')

    robust_766    = (fdr_766 < 0.05) & (fdr_clrmc_766 < 0.05)
    robust_up_766 = robust_766 & (l2fc_766 > 0)
    robust_dn_766 = robust_766 & (l2fc_766 < 0)
    logging.info(f"  Robust DA: total={robust_766.sum()} up={robust_up_766.sum()} down={robust_dn_766.sum()}")

    da_766 = pd.DataFrame({
        'species': species_766,
        'mean_Normal': mean_n766, 'sd_Normal': sd_n766,
        'mean_Caries': mean_c766, 'sd_Caries': sd_c766,
        'log2FC': l2fc_766,
        'Wilcoxon_p': p_wil_766, 'Wilcoxon_FDR': fdr_766,
        'CLRMC_p': p_clrmc_766, 'CLRMC_FDR': fdr_clrmc_766,
        'Robust_DA': robust_766,
        'Robust_direction': np.where(robust_up_766, 'Caries',
                             np.where(robust_dn_766, 'Normal', 'NS')),
    })
    da_766.to_csv(os.path.join(tab_dir, 'Table_PRJNA766357_differential_species.csv'), index=False)

    # ────────────────────────────────────────────────────────
    # M3. PRJNA766357 ALPHA & BETA DIVERSITY
    # ────────────────────────────────────────────────────────
    logging.info("\n" + "=" * 70)
    logging.info("[M3] PRJNA766357 diversity analysis")

    samples_766 = p766_normal + p766_caries
    mat_766     = np.clip(abund[samples_766].T.values, 0, None)
    groups_766  = [group_map[s] for s in samples_766]

    obs_sp  = (mat_766 > 0).sum(axis=1)
    shannon = np.array([calc_shannon(mat_766[i]) for i in range(mat_766.shape[0])])
    simpson = np.array([calc_simpson(mat_766[i]) for i in range(mat_766.shape[0])])

    alpha_df = pd.DataFrame({
        'sample': samples_766, 'group': groups_766,
        'Observed_Species': obs_sp, 'Shannon_H': shannon, 'Simpson_D': simpson,
    })
    alpha_df.to_csv(os.path.join(tab_dir, 'Table_PRJNA766357_alpha_diversity.csv'), index=False)

    for metric in ['Observed_Species', 'Shannon_H', 'Simpson_D']:
        vn = alpha_df.loc[alpha_df.group == 'Normal', metric].values
        vc = alpha_df.loc[alpha_df.group == 'Caries', metric].values
        _, p = stats.mannwhitneyu(vn, vc, alternative='two-sided')
        logging.info(f"  {metric}: Normal median={np.median(vn):.3f}, "
                     f"Caries median={np.median(vc):.3f}, p={p:.4g}")

    # Beta diversity
    bc_766 = squareform(pdist(mat_766, metric='braycurtis'))
    dm_766 = DistanceMatrix(bc_766, ids=samples_766)
    res_766 = sk_pcoa(dm_766)
    pc_766  = res_766.samples[['PC1', 'PC2']].values
    var1_766 = res_766.proportion_explained['PC1'] * 100
    var2_766 = res_766.proportion_explained['PC2'] * 100

    perm_766 = permanova(dm_766, groups_766, permutations=9999)
    pF766, pP766 = perm_766['test statistic'], perm_766['p-value']
    pR2_766 = pF766 / (pF766 + len(samples_766) - 2)

    anos_766 = anosim(dm_766, groups_766, permutations=9999)
    aR766, aP766 = anos_766['test statistic'], anos_766['p-value']

    ni766 = [i for i, s in enumerate(samples_766) if group_map[s] == 'Normal']
    ci766 = [i for i, s in enumerate(samples_766) if group_map[s] == 'Caries']
    dn766 = np.linalg.norm(pc_766[ni766] - pc_766[ni766].mean(0), axis=1)
    dc766 = np.linalg.norm(pc_766[ci766] - pc_766[ci766].mean(0), axis=1)
    lF766, lP766 = stats.levene(dn766, dc766)

    logging.info(f"  PERMANOVA R2={pR2_766:.4f} p={pP766:.4g}")
    logging.info(f"  ANOSIM    R={aR766:.4f}  p={aP766:.4g}")
    logging.info(f"  PERMDISP  F={lF766:.4f}  p={lP766:.4g}")

    beta_stats_766 = pd.DataFrame([
        {'test': 'PERMANOVA', 'stat': 'F', 'value': round(pF766, 4), 'R2': round(pR2_766, 4), 'p': pP766},
        {'test': 'ANOSIM',    'stat': 'R', 'value': round(aR766, 4), 'R2': None,               'p': aP766},
        {'test': 'PERMDISP',  'stat': 'F', 'value': round(lF766, 4), 'R2': None,               'p': lP766},
    ])
    beta_stats_766.to_csv(os.path.join(tab_dir, 'Table_PRJNA766357_beta_stats.csv'), index=False)

    # ────────────────────────────────────────────────────────
    # M4. ALL-SAMPLE BETA DIVERSITY & BATCH EFFECT
    # ────────────────────────────────────────────────────────
    logging.info("\n" + "=" * 70)
    logging.info("[M4] All-sample beta diversity & batch effect")

    mat_all      = np.clip(abund[all_samp].T.values, 0, None)
    groups_all   = [group_map[s] for s in all_samp]
    datasets_all = [dataset_map[s] for s in all_samp]

    bc_all  = squareform(pdist(mat_all, metric='braycurtis'))
    dm_all  = DistanceMatrix(bc_all, ids=all_samp)
    res_all = sk_pcoa(dm_all)
    pc_all  = res_all.samples[['PC1', 'PC2']].values
    var1_all = res_all.proportion_explained['PC1'] * 100
    var2_all = res_all.proportion_explained['PC2'] * 100

    perm_ds  = permanova(dm_all, datasets_all, permutations=9999)
    pF_ds, pP_ds = perm_ds['test statistic'], perm_ds['p-value']
    pR2_ds = pF_ds / (pF_ds + len(all_samp) - 2)

    perm_grp = permanova(dm_all, groups_all, permutations=9999)
    pF_grp, pP_grp = perm_grp['test statistic'], perm_grp['p-value']
    pR2_grp = pF_grp / (pF_grp + len(all_samp) - 2)

    logging.info(f"  PERMANOVA by Dataset (batch): R2={pR2_ds:.4f}, p={pP_ds:.4g}")
    logging.info(f"  PERMANOVA by Group (disease): R2={pR2_grp:.4f}, p={pP_grp:.4g}")
    logging.info(f"  R2_batch / R2_disease = {pR2_ds / max(pR2_grp, 1e-10):.2f}")

    pd.DataFrame([
        {'factor': 'Dataset (Batch)', 'F': round(pF_ds, 4),  'R2': round(pR2_ds, 4),  'p': pP_ds},
        {'factor': 'Group (Disease)', 'F': round(pF_grp, 4), 'R2': round(pR2_grp, 4), 'p': pP_grp},
    ]).to_csv(os.path.join(tab_dir, 'Table_AllSample_PERMANOVA_batch_vs_disease.csv'), index=False)

    # ────────────────────────────────────────────────────────
    # M5. CROSS-COHORT DA CONCORDANCE  (strict + relaxed)
    # ────────────────────────────────────────────────────────
    logging.info("\n" + "=" * 70)
    logging.info("[M5] Cross-cohort DA concordance")

    ts_merge = da_thisstudy[['species','log2FC','Wilcoxon_p','Wilcoxon_FDR',
                              'CLRMC_FDR','Robust_DA','Robust_direction']].copy()
    ts_merge.columns = ['species','log2FC_TS','p_TS','FDR_TS','CLRMC_FDR_TS','Robust_TS','Dir_TS']

    p766_merge = da_766[['species','log2FC','Wilcoxon_p','Wilcoxon_FDR',
                          'CLRMC_FDR','Robust_DA','Robust_direction']].copy()
    p766_merge.columns = ['species','log2FC_766','p_766','FDR_766','CLRMC_FDR_766','Robust_766','Dir_766']

    merged = pd.merge(ts_merge, p766_merge, on='species', how='inner')
    logging.info(f"  Shared species in both tests: {len(merged)}")

    # Strict: Robust in both
    strict = merged[(merged['Robust_TS'] == True) & (merged['Robust_766'] == True)]
    logging.info(f"  [Strict] Robust in both: {len(strict)}")

    # Relaxed levels
    for label, col_ts, col_766, thr in [
        ('FDR<0.05', 'FDR_TS', 'FDR_766', 0.05),
        ('FDR<0.10', 'FDR_TS', 'FDR_766', 0.10),
        ('p<0.05',   'p_TS',   'p_766',   0.05),
    ]:
        both = merged[(merged[col_ts] < thr) & (merged[col_766] < thr)]
        same = both[np.sign(both['log2FC_TS']) == np.sign(both['log2FC_766'])] if len(both) else pd.DataFrame()
        diff = both[np.sign(both['log2FC_TS']) != np.sign(both['log2FC_766'])] if len(both) else pd.DataFrame()
        logging.info(f"  [{label}] Both sig: {len(both)}, concordant: {len(same)}, discordant: {len(diff)}")

    # TS-Robust + 766 trend (p<0.10, same direction)
    ts_robust_766_trend = merged[
        (merged['Robust_TS'] == True) &
        (merged['p_766'] < 0.10) &
        (np.sign(merged['log2FC_TS']) == np.sign(merged['log2FC_766']))
    ]
    logging.info(f"  [Cross] TS-Robust + 766-trend(p<0.10, same dir): {len(ts_robust_766_trend)}")

    # Compile all concordant
    all_concordant = pd.DataFrame()
    for label, col_ts, col_766, thr in [
        ('FDR<0.05', 'FDR_TS', 'FDR_766', 0.05),
        ('FDR<0.10', 'FDR_TS', 'FDR_766', 0.10),
        ('p<0.05',   'p_TS',   'p_766',   0.05),
    ]:
        both = merged[(merged[col_ts] < thr) & (merged[col_766] < thr)]
        same = both[np.sign(both['log2FC_TS']) == np.sign(both['log2FC_766'])]
        if len(same):
            temp = same.copy(); temp['threshold'] = label
            all_concordant = pd.concat([all_concordant, temp])
    # Add cross-validation
    if len(ts_robust_766_trend):
        temp = ts_robust_766_trend.copy(); temp['threshold'] = 'TS-Robust+766-trend'
        all_concordant = pd.concat([all_concordant, temp])

    if len(all_concordant):
        all_concordant = all_concordant.drop_duplicates(subset='species')
        all_concordant['direction'] = np.where(all_concordant['log2FC_TS'] > 0,
                                               'Disease-enriched', 'Health-enriched')
        all_concordant = all_concordant.sort_values('log2FC_TS')
        all_concordant.to_csv(os.path.join(tab_dir, 'Table_relaxed_concordant_species.csv'), index=False)
        logging.info(f"  Total unique concordant (any relaxed level): {len(all_concordant)}")
        for _, r in all_concordant.iterrows():
            d = "Disease↑" if r['log2FC_TS'] > 0 else "Health↑"
            logging.info(f"    {r['species']}: {d} | TS FC={r['log2FC_TS']:.2f} 766 FC={r['log2FC_766']:.2f}")

    merged.to_csv(os.path.join(tab_dir, 'Table_cross_cohort_DA_comparison.csv'), index=False)

    # Save core species table
    core_table = []
    for sp in core30:
        core_table.append({
            'species': sp,
            'Prev_ThisStudy':   prev_df.loc[sp, 'ThisStudy'],
            'Prev_PRJNA597465': prev_df.loc[sp, 'PRJNA597465'],
            'Prev_PRJNA766357': prev_df.loc[sp, 'PRJNA766357'],
            'Core_70': sp in core70, 'Core_50': sp in core50, 'Core_30': sp in core30,
            'MeanAbund': grand_mean[sp],
            'DA_ThisStudy':   (da_thisstudy[da_thisstudy['species']==sp]['Robust_direction'].values[0]
                               if sp in da_thisstudy['species'].values else 'Not tested'),
            'DA_PRJNA766357': (da_766[da_766['species']==sp]['Robust_direction'].values[0]
                               if sp in da_766['species'].values else 'Not tested'),
        })
    pd.DataFrame(core_table).sort_values('MeanAbund', ascending=False).to_csv(
        os.path.join(tab_dir, 'Table_core_species_summary.csv'), index=False)

    # ════════════════════════════════════════════════════════
    # M6. VISUALISATIONS  (9 figures)
    # ════════════════════════════════════════════════════════
    logging.info("\n" + "=" * 70)
    logging.info("[M6] Generating figures ...")
    sns.set_style('ticks')

    # ── Fig 1: Core Species Prevalence Heatmap ──────────────
    logging.info("  Fig 1: Core species prevalence heatmap")
    core_sp = sorted(core30, key=lambda x: -grand_mean[x])[:40]

    fig, axes = plt.subplots(1, 3, figsize=(10, 10),
                              gridspec_kw={'width_ratios': [1, 1, 0.4]})
    prev_sub = prev_df.loc[core_sp]
    im = axes[0].imshow(prev_sub.values, aspect='auto', cmap='YlOrRd', vmin=0, vmax=1)
    axes[0].set_yticks(range(len(core_sp)))
    axes[0].set_yticklabels([s.replace('_', ' ') for s in core_sp], fontsize=6.5, style='italic')
    axes[0].set_xticks(range(3))
    axes[0].set_xticklabels([f'ThisStudy\n(n={len(all_samp)-len(ds_samples["PRJNA597465"])-len(ds_samples["PRJNA766357"])})',
                             f'PRJNA597465\n(n={len(ds_samples["PRJNA597465"])})',
                             f'PRJNA766357\n(n={len(ds_samples["PRJNA766357"])})'], fontsize=7)
    axes[0].set_title('Prevalence', fontsize=10, fontweight='bold')
    for i in range(len(core_sp)):
        for j in range(3):
            v = prev_sub.values[i, j]
            axes[0].text(j, i, f'{v:.0%}', ha='center', va='center', fontsize=5,
                         color='white' if v > 0.6 else 'black')

    mean_ab = pd.DataFrame({ds: abund.loc[core_sp, ss].mean(axis=1) for ds, ss in ds_samples.items()})
    im2 = axes[1].imshow(np.log10(mean_ab.values + 0.001), aspect='auto', cmap='viridis')
    axes[1].set_yticks([]); axes[1].set_xticks(range(3))
    axes[1].set_xticklabels(['ThisStudy', 'PRJNA597465', 'PRJNA766357'], fontsize=7)
    axes[1].set_title('Mean Abundance (log₁₀)', fontsize=10, fontweight='bold')

    n_ds50 = (prev_df.loc[core_sp] >= 0.5).sum(axis=1)
    tier_c = {3: '#2E8B8B', 2: '#D4A843', 1: '#B85450'}
    axes[2].barh(range(len(core_sp)), n_ds50.values,
                 color=[tier_c.get(v, '#ccc') for v in n_ds50.values],
                 height=0.8, edgecolor='white', linewidth=0.3)
    axes[2].set_yticks([]); axes[2].set_xlim(0, 3.5); axes[2].set_xticks([1,2,3])
    axes[2].set_xlabel('# Cohorts (≥50%)', fontsize=8)
    axes[2].set_title('Core Tier', fontsize=10, fontweight='bold')
    axes[2].invert_yaxis()
    axes[2].legend(handles=[mpatches.Patch(fc=tier_c[k], label=f'{k}/3') for k in [3,2,1]],
                   fontsize=6, loc='lower right')

    cbar_ax = fig.add_axes([0.05, 0.02, 0.25, 0.012])
    plt.colorbar(im, cax=cbar_ax, orientation='horizontal', label='Prevalence')
    fig.suptitle('Core Species of Adolescent Dental Plaque Microbiome\nAcross Three Independent Cohorts',
                 fontsize=11, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0.04, 1, 0.95])
    save_fig(fig, os.path.join(fig_dir, 'Fig1_CoreSpecies_Prevalence'), DPI)
    plt.close(fig)

    # ── Fig 2: Dual Volcano ─────────────────────────────────
    logging.info("  Fig 2: Dual volcano plot")
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    for ax_i, (ax, da_df, fc_col, fdr_col, rob_col, title, xlabel) in enumerate([
        (axes[0], da_thisstudy, 'log2FC', 'Wilcoxon_FDR', 'Robust_DA',
         'A. ThisStudy (Normal vs PoorPlaque)', r'$\log_2$FC (PoorPlaque vs Normal)'),
        (axes[1], da_766, 'log2FC', 'Wilcoxon_FDR', 'Robust_DA',
         'B. PRJNA766357 (Normal vs Caries)', r'$\log_2$FC (Caries vs Normal)'),
    ]):
        fc  = da_df[fc_col].values
        fdr = da_df[fdr_col].values
        rob = da_df[rob_col].values
        yv  = -np.log10(fdr + 1e-300)
        r_up = rob & (fc > 0); r_dn = rob & (fc < 0)
        w_up = (~rob) & (fdr < 0.05) & (fc > 0); w_dn = (~rob) & (fdr < 0.05) & (fc < 0)
        ns = ~(r_up | r_dn | w_up | w_dn)
        for mask, c, s, a, z in [
            (ns, VOL_NS, 4, 0.55, 2), (w_up, VOL_WIL_UP, 9, 0.7, 3),
            (w_dn, VOL_WIL_DN, 9, 0.7, 3), (r_up, VOL_ROBUST_UP, 20, 0.9, 5),
            (r_dn, VOL_ROBUST_DN, 20, 0.9, 5),
        ]:
            ax.scatter(fc[mask], yv[mask], c=c, s=s, alpha=a, linewidths=0, zorder=z)
        ax.axhline(-np.log10(0.05), color='#666', lw=0.8, ls='--')
        ax.axvline(1, color='#666', lw=0.8, ls='--')
        ax.axvline(-1, color='#666', lw=0.8, ls='--')
        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel(r'$-\log_{10}$(Wilcoxon FDR)', fontsize=10)
        ax.set_title(title, fontsize=10, fontweight='bold')
        ax.text(0.03, 0.97, f'Robust UP: {r_up.sum()}\nRobust DOWN: {r_dn.sum()}',
                transform=ax.transAxes, fontsize=7, va='top',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, ec='#aaa'))
        sns.despine(ax=ax)

    leg = [Line2D([0],[0], marker='o', color='w', markerfacecolor=c, markersize=ms, label=l)
           for c, ms, l in [
               (VOL_ROBUST_UP, 6, 'Robust UP (Disease)'), (VOL_WIL_UP, 5, 'Wilcoxon-only UP'),
               (VOL_ROBUST_DN, 6, 'Robust DOWN (Normal)'), (VOL_WIL_DN, 5, 'Wilcoxon-only DOWN'),
               (VOL_NS, 4, 'Not Significant')]]
    fig.legend(handles=leg, loc='lower center', ncol=5, fontsize=7,
               frameon=True, edgecolor='#ccc', bbox_to_anchor=(0.5, -0.02))
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    save_fig(fig, os.path.join(fig_dir, 'Fig2_DualVolcano'), DPI)
    plt.close(fig)

    # ── Fig 3: Cross-cohort concordance scatter ─────────────
    logging.info("  Fig 3: Cross-cohort concordance scatter")
    fig, ax = plt.subplots(figsize=(6, 6))
    rho, rpv = spearmanr(merged['log2FC_TS'], merged['log2FC_766'])

    def get_cat(row):
        ts_sig = row['FDR_TS'] < 0.05; p766_sig = row['FDR_766'] < 0.05
        ts_nom = row['p_TS'] < 0.05;   p766_nom = row['p_766'] < 0.05
        same = np.sign(row['log2FC_TS']) == np.sign(row['log2FC_766'])
        if ts_sig and p766_sig and same:   return 'FDR<0.05 both (concordant)'
        if ts_nom and p766_nom and same:   return 'p<0.05 both (concordant)'
        if ts_nom and p766_nom and not same: return 'p<0.05 both (discordant)'
        if ts_sig or p766_sig:             return 'FDR<0.05 in one'
        return 'NS'

    merged['category'] = merged.apply(get_cat, axis=1)
    cat_style = {
        'NS': ('#DDDDDD', 10, 0.4, 'o', 1),
        'FDR<0.05 in one': ('#AAAAAA', 18, 0.6, 'o', 2),
        'p<0.05 both (discordant)': ('#E74C3C', 35, 0.8, 'X', 4),
        'p<0.05 both (concordant)': ('#3498DB', 35, 0.8, 'D', 4),
        'FDR<0.05 both (concordant)': ('#2C3E50', 60, 0.95, '*', 5),
    }
    for cat, (color, size, alpha, marker, zorder) in cat_style.items():
        m = merged['category'] == cat
        if m.sum() == 0: continue
        ax.scatter(merged.loc[m, 'log2FC_TS'], merged.loc[m, 'log2FC_766'],
                   c=color, s=size, alpha=alpha, marker=marker, zorder=zorder,
                   linewidths=0.3, edgecolors='white', label=f'{cat} (n={m.sum()})')
    for cat_name in ['FDR<0.05 both (concordant)', 'p<0.05 both (concordant)']:
        for _, r in merged[merged['category'] == cat_name].iterrows():
            ax.annotate(r['species'].replace('_', ' '), (r['log2FC_TS'], r['log2FC_766']),
                        fontsize=5, style='italic', xytext=(4, 4), textcoords='offset points',
                        bbox=dict(boxstyle='round,pad=0.15', facecolor='lightyellow', alpha=0.8, ec='#ddd', lw=0.3))
    for _, r in ts_robust_766_trend.iterrows():
        if r['species'] not in merged[merged['category'].str.contains('concordant')]['species'].values:
            ax.annotate(r['species'].replace('_', ' '), (r['log2FC_TS'], r['log2FC_766']),
                        fontsize=4.5, style='italic', xytext=(3,3), textcoords='offset points', color='#666',
                        bbox=dict(boxstyle='round,pad=0.1', facecolor='white', alpha=0.7, edgecolor='#ccc', lw=0.2))

    lims = [min(merged['log2FC_TS'].min(), merged['log2FC_766'].min()) - 0.5,
            max(merged['log2FC_TS'].max(), merged['log2FC_766'].max()) + 0.5]
    ax.plot(lims, lims, 'k--', lw=0.6, alpha=0.3, zorder=0)
    ax.axhline(0, color='#ccc', lw=0.5, zorder=0); ax.axvline(0, color='#ccc', lw=0.5, zorder=0)
    ax.fill_between([0, lims[1]], 0, lims[1], alpha=0.03, color=COL['PoorPlaque'], zorder=0)
    ax.fill_between([lims[0], 0], lims[0], 0, alpha=0.03, color=COL['Normal'], zorder=0)
    ax.text(0.03, 0.97, f'Spearman ρ = {rho:.3f}\np = {rpv:.2e}\nn = {len(merged)}',
            transform=ax.transAxes, fontsize=7, va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85, ec='#aaa'))
    ax.set_xlabel(r'$\log_2$FC ThisStudy', fontsize=10)
    ax.set_ylabel(r'$\log_2$FC PRJNA766357', fontsize=10)
    ax.set_title('Cross-cohort Effect Size Concordance (Relaxed)', fontsize=11, fontweight='bold')
    ax.legend(fontsize=5.5, loc='lower right', frameon=True, edgecolor='#ccc')
    sns.despine(ax=ax)
    plt.tight_layout()
    save_fig(fig, os.path.join(fig_dir, 'Fig3_CrossCohort_Concordance'), DPI)
    plt.close(fig)

    # ── Fig 4: PRJNA766357 Alpha Diversity ──────────────────
    logging.info("  Fig 4: PRJNA766357 alpha diversity")
    fig, axes = plt.subplots(1, 3, figsize=(9, 3.5))
    for idx, (metric, ylabel) in enumerate([
        ('Observed_Species', 'Observed Species'),
        ('Shannon_H', "Shannon H'"),
        ('Simpson_D', "Simpson's D"),
    ]):
        vn = alpha_df.loc[alpha_df.group == 'Normal', metric].values
        vc = alpha_df.loc[alpha_df.group == 'Caries', metric].values
        _, pv = stats.mannwhitneyu(vn, vc, alternative='two-sided')
        boxplot_two_groups(axes[idx], vn, vc, ['Normal', 'Caries'],
                           [COL['Normal'], COL['Caries']], ylabel, pv, SEED)
        axes[idx].text(-0.25, 1.10, chr(65+idx), transform=axes[idx].transAxes,
                       fontsize=13, fontweight='bold', va='top')
    fig.suptitle('PRJNA766357: Alpha Diversity (Normal vs Caries)', fontsize=11,
                 fontweight='bold', y=1.02)
    plt.tight_layout()
    save_fig(fig, os.path.join(fig_dir, 'Fig4_PRJNA766357_AlphaDiversity'), DPI)
    plt.close(fig)

    # ── Fig 5: PRJNA766357 PCoA ─────────────────────────────
    logging.info("  Fig 5: PRJNA766357 beta diversity PCoA")
    fig, ax = plt.subplots(figsize=(5, 4.5))
    ax.axhline(0, color='#CCC', lw=0.5, zorder=0); ax.axvline(0, color='#CCC', lw=0.5, zorder=0)
    for grp, idx_list in [('Normal', ni766), ('Caries', ci766)]:
        pts = pc_766[idx_list]
        ax.scatter(pts[:,0], pts[:,1], c=COL[grp], s=35, alpha=0.85, label=grp, zorder=5, linewidths=0)
        cov_ellipse(pts[:,0], pts[:,1], ax, edgecolor=COL[grp], linewidth=1.0, facecolor=COL[grp], alpha=0.12)
    ax.set_xlabel(f'PCo1 ({var1_766:.1f}%)'); ax.set_ylabel(f'PCo2 ({var2_766:.1f}%)')
    ax.text(0.03, 0.03, f'PERMANOVA\n$R^2$={pR2_766:.4f}\np={pP766:.4f}',
            transform=ax.transAxes, fontsize=6.5, va='bottom',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', ec='#AAA', alpha=0.85))
    ax.legend(fontsize=8, frameon=False, loc='upper right')
    ax.set_title('PRJNA766357: PCoA (Bray-Curtis)', fontsize=11, fontweight='bold')
    sns.despine(ax=ax); plt.tight_layout()
    save_fig(fig, os.path.join(fig_dir, 'Fig5_PRJNA766357_BetaDiversity'), DPI)
    plt.close(fig)

    # ── Fig 6: All-sample PCoA (batch effect) ──────────────
    logging.info("  Fig 6: All-sample PCoA batch effect")
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax in axes:
        ax.axhline(0, color='#CCC', lw=0.5, zorder=0); ax.axvline(0, color='#CCC', lw=0.5, zorder=0)

    # Panel A: by dataset
    for ds in ['ThisStudy', 'PRJNA597465', 'PRJNA766357']:
        idx = [i for i, s in enumerate(all_samp) if dataset_map[s] == ds]
        pts = pc_all[idx]
        mk = {'ThisStudy':'o', 'PRJNA597465':'s', 'PRJNA766357':'^'}[ds]
        axes[0].scatter(pts[:,0], pts[:,1], c=DS_COL[ds], s=30, alpha=0.8, label=ds, marker=mk,
                        zorder=5, linewidths=0.3, edgecolors='white')
        cov_ellipse(pts[:,0], pts[:,1], axes[0], edgecolor=DS_COL[ds], linewidth=1.2,
                    facecolor=DS_COL[ds], alpha=0.08)
    axes[0].set_xlabel(f'PCo1 ({var1_all:.1f}%)'); axes[0].set_ylabel(f'PCo2 ({var2_all:.1f}%)')
    axes[0].text(0.03, 0.03, f'PERMANOVA (Dataset)\n$R^2$={pR2_ds:.4f}\np={pP_ds:.4f}',
                 transform=axes[0].transAxes, fontsize=7, va='bottom',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='white', ec='#AAA', alpha=0.85))
    axes[0].legend(fontsize=7, frameon=True, edgecolor='#ccc', loc='upper right')
    axes[0].set_title('A. Colored by Dataset (Batch)', fontsize=10, fontweight='bold')
    sns.despine(ax=axes[0])

    # Panel B: by group
    gm_map = {
        ('ThisStudy','Normal'):('o',COL['Normal']), ('ThisStudy','PoorPlaque'):('o',COL['PoorPlaque']),
        ('PRJNA597465','Unknown'):('s',COL['Unknown']),
        ('PRJNA766357','Normal'):('^',COL['Normal']), ('PRJNA766357','Caries'):('^',COL['Caries']),
    }
    lh = []
    for (ds, grp), (mk, color) in gm_map.items():
        idx = [i for i, s in enumerate(all_samp) if dataset_map[s]==ds and group_map[s]==grp]
        if not idx: continue
        axes[1].scatter(pc_all[idx,0], pc_all[idx,1], c=color, s=30, alpha=0.8, marker=mk,
                        zorder=5, linewidths=0.3, edgecolors='white')
        lh.append(Line2D([0],[0], marker=mk, color='w', markerfacecolor=color, markersize=7,
                         label=f'{ds} – {grp}'))
    axes[1].set_xlabel(f'PCo1 ({var1_all:.1f}%)'); axes[1].set_ylabel(f'PCo2 ({var2_all:.1f}%)')
    axes[1].text(0.03, 0.03, f'PERMANOVA (Group)\n$R^2$={pR2_grp:.4f}\np={pP_grp:.4f}',
                 transform=axes[1].transAxes, fontsize=7, va='bottom',
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='white', ec='#AAA', alpha=0.85))
    axes[1].legend(handles=lh, fontsize=5.5, frameon=True, edgecolor='#ccc', loc='upper right')
    axes[1].set_title('B. Colored by Disease Group', fontsize=10, fontweight='bold')
    sns.despine(ax=axes[1])
    fig.suptitle('All-sample PCoA: Batch Effect vs Disease Effect', fontsize=11, fontweight='bold', y=1.02)
    plt.tight_layout()
    save_fig(fig, os.path.join(fig_dir, 'Fig6_AllSample_PCoA_BatchEffect'), DPI)
    plt.close(fig)

    # ── Fig 7: Relaxed DA concordance bar + scatter ─────────
    logging.info("  Fig 7: Relaxed DA concordance")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    # reuse merged with category from Fig3
    ax = axes[0]
    for cat, (color, size, alpha, marker, zorder) in cat_style.items():
        m = merged['category'] == cat
        if m.sum() == 0: continue
        ax.scatter(merged.loc[m, 'log2FC_TS'], merged.loc[m, 'log2FC_766'],
                   c=color, s=size, alpha=alpha, marker=marker, zorder=zorder,
                   linewidths=0.3, edgecolors='white', label=f'{cat} (n={m.sum()})')
    ax.plot(lims, lims, 'k--', lw=0.6, alpha=0.3, zorder=0)
    ax.axhline(0, color='#ccc', lw=0.5); ax.axvline(0, color='#ccc', lw=0.5)
    ax.text(0.03, 0.97, f'Spearman ρ = {rho:.3f}', transform=ax.transAxes, fontsize=7, va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.85, ec='#aaa'))
    ax.set_xlabel(r'$\log_2$FC ThisStudy'); ax.set_ylabel(r'$\log_2$FC PRJNA766357')
    ax.set_title('A. Cross-cohort Effect Size', fontsize=10, fontweight='bold')
    ax.legend(fontsize=5.5, loc='lower right', frameon=True, edgecolor='#ccc')
    sns.despine(ax=ax)

    # Bar chart
    ax = axes[1]
    # Recompute for bar
    levels_data = []
    for lab, ct, c7, thr in [('Robust\nboth','Robust_TS','Robust_766',None),
                              ('FDR<0.05\nboth','FDR_TS','FDR_766',0.05),
                              ('FDR<0.10\nboth','FDR_TS','FDR_766',0.10),
                              ('p<0.05\nboth','p_TS','p_766',0.05)]:
        if thr is None:
            both = merged[(merged['Robust_TS']==True)&(merged['Robust_766']==True)]
        else:
            both = merged[(merged[ct]<thr)&(merged[c7]<thr)]
        nc = len(both[np.sign(both['log2FC_TS'])==np.sign(both['log2FC_766'])]) if len(both) else 0
        nd = len(both[np.sign(both['log2FC_TS'])!=np.sign(both['log2FC_766'])]) if len(both) else 0
        levels_data.append((lab, nc, nd))
    levels_data.append(('TS-Robust\n+766 trend', len(ts_robust_766_trend), 0))

    x = np.arange(len(levels_data)); w = 0.35
    nc_vals = [d[1] for d in levels_data]; nd_vals = [d[2] for d in levels_data]
    b1 = ax.bar(x-w/2, nc_vals, w, color=COL['Normal'], alpha=0.85, label='Concordant', ec='white')
    b2 = ax.bar(x+w/2, nd_vals, w, color=COL['PoorPlaque'], alpha=0.85, label='Discordant', ec='white')
    for bar, val in zip(b1, nc_vals):
        if val: ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.2, str(val), ha='center', fontsize=8, fontweight='bold')
    for bar, val in zip(b2, nd_vals):
        if val: ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+0.2, str(val), ha='center', fontsize=8, fontweight='bold')
    ax.set_xticks(x); ax.set_xticklabels([d[0] for d in levels_data], fontsize=7)
    ax.set_ylabel('Number of Species')
    ax.set_title('B. Concordance at Relaxed Thresholds', fontsize=10, fontweight='bold')
    ax.legend(fontsize=7, frameon=True, edgecolor='#ccc')
    ax.annotate('', xy=(4.3, max(nc_vals)*0.7), xytext=(-0.3, max(nc_vals)*0.7),
                arrowprops=dict(arrowstyle='->', color='#666', lw=1.5))
    ax.text(2, max(nc_vals)*0.78, 'Decreasing stringency →', ha='center', fontsize=7, color='#666', style='italic')
    sns.despine(ax=ax)
    plt.tight_layout()
    save_fig(fig, os.path.join(fig_dir, 'Fig7_RelaxedDA_Concordance'), DPI)
    plt.close(fig)

    # ── Fig 8: Concordant species boxplots ──────────────────
    if len(all_concordant) > 0:
        logging.info("  Fig 8: Concordant species boxplots")
        n_sp = len(all_concordant)
        ncols = min(4, n_sp); nrows = max(1, math.ceil(n_sp / ncols))
        fig, axg = plt.subplots(nrows, ncols, figsize=(ncols*2.8, nrows*3.5))
        axf = np.array(axg).flatten() if n_sp > 1 else [axg]
        for idx, (_, row) in enumerate(all_concordant.iterrows()):
            ax = axf[idx]; sp = row['species']
            if sp not in abund.index: ax.set_visible(False); continue
            v_tn = abund.loc[sp, ts_normal].values.astype(float)
            v_tp = abund.loc[sp, ts_poor].values.astype(float)
            v_7n = abund.loc[sp, p766_normal].values.astype(float)
            v_7c = abund.loc[sp, p766_caries].values.astype(float)
            pos = [0, 1, 2.5, 3.5]; dg = [v_tn, v_tp, v_7n, v_7c]
            cb = [COL['Normal'], COL['PoorPlaque'], COL['Normal'], COL['Caries']]
            rng = np.random.default_rng(SEED)
            for i, (p, vals, c) in enumerate(zip(pos, dg, cb)):
                h = '//' if i >= 2 else ''
                ax.boxplot(vals, positions=[p], widths=0.45, patch_artist=True,
                           medianprops=dict(color='white', lw=1.5),
                           whiskerprops=dict(color=c, lw=0.6), capprops=dict(color=c, lw=0.6),
                           flierprops=dict(marker=''), boxprops=dict(facecolor=c, alpha=0.6, lw=0.6, ec=c, hatch=h))
                jx = rng.uniform(-0.10, 0.10, len(vals))
                ax.scatter(np.full(len(vals), p)+jx, vals, color=c, alpha=0.6, s=15, zorder=5, linewidths=0)
            _, p_ts = stats.mannwhitneyu(v_tn, v_tp, alternative='two-sided')
            _, p_766 = stats.mannwhitneyu(v_7n, v_7c, alternative='two-sided')
            av = np.concatenate(dg); ym = av.max(); yr = max(ym-av.min(), abs(ym)*0.1, 0.001)
            yb1 = ym + yr*0.08
            ax.plot([0,0,1,1],[yb1-yr*0.01,yb1,yb1,yb1-yr*0.01], lw=0.6, color='#333')
            ax.text(0.5, yb1+yr*0.02, f'p={p_ts:.3f}', ha='center', fontsize=5)
            yb2 = ym + yr*0.25
            ax.plot([2.5,2.5,3.5,3.5],[yb2-yr*0.01,yb2,yb2,yb2-yr*0.01], lw=0.6, color='#333')
            ax.text(3.0, yb2+yr*0.02, f'p={p_766:.3f}', ha='center', fontsize=5)
            ax.set_xticks(pos); ax.set_xticklabels(['TS\nNorm','TS\nPoor','766\nNorm','766\nCar'], fontsize=5.5)
            ax.set_ylabel('Rel. abund. (%)', fontsize=6)
            dc = COL['PoorPlaque'] if row['log2FC_TS'] > 0 else COL['Normal']
            ax.set_title(f"{sp.replace('_',' ')}\n(TS FC={row['log2FC_TS']:.2f}, 766 FC={row['log2FC_766']:.2f})",
                         fontsize=6, style='italic', color=dc, pad=3)
            ax.set_xlim(-0.7, 4.2); sns.despine(ax=ax)
        for i in range(n_sp, len(axf)): axf[i].set_visible(False)
        fig.suptitle('Concordant DA Species Across Two Cohorts', fontsize=10, fontweight='bold', y=1.01)
        plt.tight_layout(pad=0.8, h_pad=2.0, w_pad=0.8)
        save_fig(fig, os.path.join(fig_dir, 'Fig8_Concordant_Species_Boxplots'), DPI)
        plt.close(fig)

    # ── Fig 9: Core species threshold summary ──────────────
    logging.info("  Fig 9: Core species threshold summary")
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    thrs_plot = [0.3, 0.5, 0.7]
    cc = {'3/3 cohorts': [len(core_results[t]['in_all_3']) for t in thrs_plot],
          '≥2/3 cohorts': [len(core_results[t]['in_2']) for t in thrs_plot],
          '≥1/3 cohorts': [len(core_results[t]['in_1']) for t in thrs_plot]}
    x = np.arange(len(thrs_plot)); w = 0.25
    clrs = ['#2E8B8B', '#D4A843', '#F4A582']
    for i, (lab, vals) in enumerate(cc.items()):
        axes[0].bar(x + i*w, vals, w, label=lab, color=clrs[i], edgecolor='white')
        for xi, c in zip(x + i*w, vals):
            axes[0].text(xi, c+1, str(c), ha='center', va='bottom', fontsize=7)
    axes[0].set_xticks(x+w); axes[0].set_xticklabels([f'≥{int(t*100)}%' for t in thrs_plot])
    axes[0].set_xlabel('Prevalence Threshold'); axes[0].set_ylabel('Number of Species')
    axes[0].set_title('A. Core Species Count by Threshold', fontweight='bold')
    axes[0].legend(fontsize=7, frameon=True); sns.despine(ax=axes[0])

    # Venn
    ax = axes[1]; ax.axis('off')
    robust_ts_sp  = set(da_thisstudy[da_thisstudy['Robust_DA']==True]['species'])
    robust_766_sp = set(da_766[da_766['Robust_DA']==True]['species'])
    c1 = plt.Circle((0.35,0.55),0.25,fill=True,fc=COL['PoorPlaque'],alpha=0.2,ec=COL['PoorPlaque'],lw=2)
    c2 = plt.Circle((0.65,0.55),0.25,fill=True,fc=COL['Caries'],alpha=0.2,ec=COL['Caries'],lw=2)
    ax.add_patch(c1); ax.add_patch(c2)
    ts_o = robust_ts_sp - robust_766_sp; p7_o = robust_766_sp - robust_ts_sp; bo = robust_ts_sp & robust_766_sp
    ax.text(0.22,0.55,f'{len(ts_o)}',fontsize=16,fontweight='bold',ha='center',va='center',color=COL['PoorPlaque'])
    ax.text(0.78,0.55,f'{len(p7_o)}',fontsize=16,fontweight='bold',ha='center',va='center',color=COL['Caries'])
    ax.text(0.50,0.55,f'{len(bo)}',fontsize=16,fontweight='bold',ha='center',va='center',color='#333')
    ax.text(0.22,0.80,'ThisStudy\nonly',fontsize=8,ha='center',color=COL['PoorPlaque'],fontweight='bold')
    ax.text(0.78,0.80,'PRJNA766357\nonly',fontsize=8,ha='center',color=COL['Caries'],fontweight='bold')
    ax.text(0.50,0.80,'Shared',fontsize=8,ha='center',color='#333',fontweight='bold')
    ax.text(0.50,0.18,f'Robust DA: TS={len(robust_ts_sp)} | 766={len(robust_766_sp)} | Shared={len(bo)}',
            fontsize=8,ha='center',bbox=dict(boxstyle='round,pad=0.5',facecolor='#f9f9f9',edgecolor='#ccc'))
    ax.set_xlim(0,1); ax.set_ylim(0,1)
    ax.set_title('B. Cross-cohort Robust DA Overlap', fontweight='bold', fontsize=10)
    plt.tight_layout()
    save_fig(fig, os.path.join(fig_dir, 'Fig9_CoreSpecies_Summary'), DPI)
    plt.close(fig)

    # ────────────────────────────────────────────────────────
    logging.info("\n" + "=" * 70)
    logging.info(f"[DONE] All outputs → {OUT}/")
    logging.info("=" * 70)


# ══════════════════════════════════════════════════════════════
#  CLI
# ══════════════════════════════════════════════════════════════
def parse_args():
    p = argparse.ArgumentParser(
        description='Adolescent dental plaque metagenome: cross-cohort analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    p.add_argument('--abund',  '-a', required=True,
                   help='Species-level relative abundance table (MetaPhlAn, TSV, species × samples)')
    p.add_argument('--group',  '-g', required=True,
                   help='Sample metadata (TSV: Sample, Group, Dataset)')
    p.add_argument('--da-ts',  '-d', required=True,
                   help='ThisStudy differential species table (CSV from v2.py)')
    p.add_argument('--output', '-o', default='./results',
                   help='Output directory (default: ./results)')
    p.add_argument('--seed',   type=int, default=42,
                   help='Random seed (default: 42)')
    p.add_argument('--dpi',    type=int, default=300,
                   help='Figure DPI (default: 300)')
    p.add_argument('--verbose','-v', action='store_true')
    return p.parse_args()


if __name__ == '__main__':
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)
    log_path = os.path.join(args.output, 'analysis_log.txt')
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.StreamHandler(sys.stdout),
                  logging.FileHandler(log_path, mode='w')])
    np.random.seed(args.seed)
    run(args)
