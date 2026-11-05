"""Packaged for revise/; source: project/01this_study_dataset/M1-metaphalan/v2/v2.py
Figures: Fig1A;Fig1B;Fig1C;FigS1A;FigS1B
Note: MetaPhlAn alpha/beta/DA composite
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
Module 1: Figure 1 - Oral Microbiome Dysbiosis in Dental Plaque
Version 3.0  --  100% vector SVG output (no PNG-embed hack)

All panel drawing functions accept axes parameters.
Composite figures call drawing functions directly -> no raster embedding.

Usage:
    python3 M1_Figure1_v3.py \
        --abund   this_study_s_relative_pop.xls \
        --group   this_study_group.tsv \
        --output  ./Figure1 \
        --seed    42 \
        --dpi     300
"""

import argparse, os, logging, sys, math, warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.spatial.distance import pdist, squareform
from statsmodels.stats.multitest import multipletests
import seaborn as sns

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
from adjustText import adjust_text

matplotlib.rcParams.update({
    'font.family':     'DejaVu Sans',
    'font.size':        9,
    'axes.titlesize':   9,
    'axes.labelsize':  10,
    'xtick.labelsize':  8,
    'ytick.labelsize':  8,
    'legend.fontsize':  8,
    'figure.dpi':      300,
    'savefig.dpi':     300,
    'savefig.bbox':   'tight',
    'axes.linewidth':   0.8,
    'lines.linewidth':  1.0,
    'pdf.fonttype':    42,
    'ps.fonttype':     42,
    'svg.fonttype':   'none',   # text as <text> elements, editable in Illustrator
})

# ── global palette ────────────────────────────────────────────
COL  = {'Normal': '#2E8B8B', 'PoorPlaque': '#B85450'}
SEED = 42

# ── literature-selected DA species ───────────────────────────
MAIN_UP = [
    'Leptotrichia_wadei',
    'Tannerella_serpentiformis',
    'Prevotella_melaninogenica',
    'Capnocytophaga_leadbetteri',
    'Actinomyces_dentalis',
]
MAIN_DOWN = [
    'Haemophilus_parahaemolyticus',
    'Streptococcus_mitis',
    'Neisseria_subflava',
    'Actinomyces_oris',
    'Streptococcus_oralis',
]
MAIN_SPECIES = MAIN_UP + MAIN_DOWN


# ════════════════════════════════════════════════════════════════
# I/O
# ════════════════════════════════════════════════════════════════
def load_data(abund_path, group_path):
    abund    = pd.read_csv(abund_path, sep='\t', index_col=0)
    group_df = pd.read_csv(group_path, sep='\t', header=None,
                           names=['sample', 'group'])
    group    = dict(zip(group_df['sample'], group_df['group']))
    sn = [s for s in abund.columns if group.get(s) == 'Normal']
    sp = [s for s in abund.columns if group.get(s) == 'PoorPlaque']
    abund  = abund[sn + sp]
    groups = pd.Series({s: group[s] for s in sn + sp}, name='group')
    return abund, groups, sn, sp


def save_fig(fig, base_path, dpi=300):
    """Save as PDF + SVG + PNG.  SVG uses svg.fonttype='none' -> pure vector."""
    fig.savefig(f"{base_path}.pdf", format='pdf', bbox_inches='tight')
    fig.savefig(f"{base_path}.svg", format='svg',  bbox_inches='tight')
    fig.savefig(f"{base_path}.png", format='png',  dpi=dpi, bbox_inches='tight')
    logging.info(f"  [SAVED] {base_path} (.pdf/.svg/.png)")


def fmt_pval(p):
    return f"p = {p:.2e}" if p < 0.001 else f"p = {p:.4f}"


# ════════════════════════════════════════════════════════════════
# Shared drawing primitives
# ════════════════════════════════════════════════════════════════
def _boxplot_on_ax(ax, vals_n, vals_p, ylabel, p_val, seed,
                   title=None, title_color='#333333',
                   panel_label=None, fontsize_title=6.5):
    """
    Draw a single two-group boxplot on *ax*.
    - No numeric annotations
    - x-tick labels rotated 45 deg
    - significance bracket + p-value on top
    This is the single source of truth for boxplot style.
    """
    for i, (grp, vals) in enumerate([('Normal', vals_n),
                                       ('PoorPlaque', vals_p)]):
        color = COL[grp]
        ax.boxplot(vals, positions=[i], widths=0.45,
                   patch_artist=True,
                   medianprops=dict(color='white', linewidth=1.5),
                   whiskerprops=dict(color=color, linewidth=0.8),
                   capprops=dict(color=color, linewidth=0.8),
                   flierprops=dict(marker='', markersize=0),
                   boxprops=dict(facecolor=color, alpha=0.75,
                                 linewidth=0.8, edgecolor=color))
        rng = np.random.default_rng(seed)
        jx  = rng.uniform(-0.12, 0.12, len(vals))
        ax.scatter(np.full(len(vals), i) + jx, vals,
                   color=color, alpha=0.75, s=22, zorder=5,
                   linewidths=0, rasterized=False)

    all_v   = np.concatenate([vals_n, vals_p])
    y_max   = all_v.max()
    y_min   = all_v.min()
    y_range = max(y_max - y_min, abs(y_max) * 0.1, 0.001)
    y_line  = y_max + y_range * 0.10
    y_text  = y_line + y_range * 0.04

    ax.plot([0, 0, 1, 1],
            [y_line - y_range * 0.015, y_line,
             y_line, y_line - y_range * 0.015],
            lw=0.8, color='#333333')
    ax.text(0.5, y_text, fmt_pval(p_val),
            ha='center', va='bottom', fontsize=6.5)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Normal', 'PoorPlaque'], fontsize=7,
                       rotation=45, ha='right')
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_xlim(-0.65, 1.65)
    ax.set_ylim(bottom=y_min  - y_range * 0.12,
                top=y_text    + y_range * 0.14)
    sns.despine(ax=ax)

    if title is not None:
        ax.set_title(title, fontsize=fontsize_title, style='italic',
                     pad=3, color=title_color)
    if panel_label is not None:
        ax.text(-0.30, 1.18, panel_label,
                transform=ax.transAxes,
                fontsize=13, fontweight='bold', va='top')


def _cov_ellipse(x, y, ax, n_std=1.96, facecolor='none', **kwargs):
    """95% confidence ellipse via covariance eigendecomposition."""
    if len(x) < 2:
        return
    cov   = np.cov(x, y)
    evals, evecs = np.linalg.eigh(cov)
    order  = evals.argsort()[::-1]
    evals, evecs = evals[order], evecs[:, order]
    angle  = np.degrees(np.arctan2(evecs[1, 0], evecs[0, 0]))
    width  = 2 * n_std * np.sqrt(max(evals[0], 0))
    height = 2 * n_std * np.sqrt(max(evals[1], 0))
    ax.add_patch(Ellipse(xy=(np.mean(x), np.mean(y)),
                         width=width, height=height, angle=angle,
                         facecolor=facecolor, **kwargs))


# ════════════════════════════════════════════════════════════════
# Panel A: draw function (accepts axes list)
# ════════════════════════════════════════════════════════════════
def _compute_alpha(abund, groups, sn, sp):
    p    = abund.div(abund.sum(axis=0), axis=1).replace(0, np.nan)
    obs  = (abund > 0).sum(axis=0)
    sh   = -(p * np.log(p)).sum(axis=0, skipna=True)
    alpha_df = pd.DataFrame({
        'sample': abund.columns, 'group': groups.values,
        'Observed_Species': obs.values, 'Shannon_H': sh.values})

    n_obs = alpha_df.loc[alpha_df.group == 'Normal',     'Observed_Species'].values
    p_obs = alpha_df.loc[alpha_df.group == 'PoorPlaque', 'Observed_Species'].values
    n_sh  = alpha_df.loc[alpha_df.group == 'Normal',     'Shannon_H'].values
    p_sh  = alpha_df.loc[alpha_df.group == 'PoorPlaque', 'Shannon_H'].values

    _, pv_obs = stats.mannwhitneyu(n_obs, p_obs, alternative='two-sided')
    _, pv_sh  = stats.mannwhitneyu(n_sh,  p_sh,  alternative='two-sided')

    return alpha_df, (n_obs, p_obs, pv_obs), (n_sh, p_sh, pv_sh)


def draw_panel_a(ax_obs, ax_sh, alpha_data, seed, panel_label='A'):
    """Draw Panel A onto two provided axes."""
    (n_obs, p_obs, pv_obs), (n_sh, p_sh, pv_sh) = alpha_data
    _boxplot_on_ax(ax_obs, n_obs, p_obs, 'Observed Species', pv_obs, seed,
                   panel_label=panel_label)
    _boxplot_on_ax(ax_sh,  n_sh,  p_sh,  "Shannon Diversity (H')", pv_sh, seed)


# ════════════════════════════════════════════════════════════════
# Panel B: draw function (accepts single ax)
# ════════════════════════════════════════════════════════════════
def _compute_beta(abund, groups, sn, sp):
    mat   = abund.T.values
    names = abund.columns.tolist()
    gvec  = groups.values

    from skbio.stats.ordination import pcoa as sk_pcoa
    from skbio.stats.distance   import DistanceMatrix, permanova, anosim

    bc    = squareform(pdist(mat, metric='braycurtis'))
    dm    = DistanceMatrix(bc, ids=names)
    res   = sk_pcoa(dm)
    pc    = res.samples[['PC1', 'PC2']].copy()
    var1  = res.proportion_explained['PC1'] * 100
    var2  = res.proportion_explained['PC2'] * 100

    perm  = permanova(dm, gvec, permutations=9999)
    pF    = perm['test statistic']
    pP    = perm['p-value']
    pR2   = pF / (pF + len(names) - 2)

    anos  = anosim(dm, gvec, permutations=9999)
    aR, aP = anos['test statistic'], anos['p-value']

    coords = pc.values
    ni = [i for i, s in enumerate(names) if groups[s] == 'Normal']
    pi = [i for i, s in enumerate(names) if groups[s] == 'PoorPlaque']
    dn = np.linalg.norm(coords[ni] - coords[ni].mean(0), axis=1)
    dp = np.linalg.norm(coords[pi] - coords[pi].mean(0), axis=1)
    lF, lP = stats.levene(dn, dp)

    logging.info(f"  PC1={var1:.1f}%  PC2={var2:.1f}%")
    logging.info(f"  PERMANOVA R2={pR2:.4f} p={pP:.4g} | "
                 f"ANOSIM R={aR:.4f} p={aP:.4g} | "
                 f"PERMDISP F={lF:.4f} p={lP:.4g}")

    return (coords, ni, pi, var1, var2,
            pR2, pP, aR, aP, lF, lP,
            pd.DataFrame([
                {'test':'PERMANOVA','stat':'F',       'value':round(pF, 4), 'R2':round(pR2,4),'p':pP},
                {'test':'ANOSIM',   'stat':'R',       'value':round(aR, 4), 'R2':None,        'p':aP},
                {'test':'PERMDISP', 'stat':'F_Levene','value':round(lF, 4), 'R2':None,        'p':lP},
                {'test':'PCoA_var', 'stat':'PC1_pct', 'value':round(var1,3),'R2':None,        'p':None},
                {'test':'PCoA_var', 'stat':'PC2_pct', 'value':round(var2,3),'R2':None,        'p':None},
            ]))


def draw_panel_b(ax, beta_data, panel_label='B'):
    """Draw Panel B (PCoA) onto a provided ax."""
    coords, ni, pi, var1, var2, pR2, pP, aR, aP, lF, lP, _ = beta_data

    ax.axhline(0, color='#CCCCCC', lw=0.5, zorder=0)
    ax.axvline(0, color='#CCCCCC', lw=0.5, zorder=0)

    for grp, idx in [('Normal', ni), ('PoorPlaque', pi)]:
        pts = coords[idx]
        ax.scatter(pts[:, 0], pts[:, 1], c=COL[grp], s=28, alpha=0.85,
                   label=grp, zorder=5, linewidths=0, rasterized=False)
        _cov_ellipse(pts[:, 0], pts[:, 1], ax, n_std=1.96,
                     edgecolor=COL[grp], linewidth=1.0,
                     facecolor=COL[grp], alpha=0.12)

    ax.set_xlabel(f'PCo1 ({var1:.1f}%)', fontsize=10)
    ax.set_ylabel(f'PCo2 ({var2:.1f}%)', fontsize=10)

    ax.text(0.03, 0.03,
            f'PERMANOVA\n$R^2$={pR2:.3f}\np={pP:.4f}',
            transform=ax.transAxes, ha='left', va='bottom', fontsize=6.5,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='#AAAAAA', linewidth=0.6, alpha=0.85))

    handles = [mpatches.Patch(facecolor=COL[g], label=g)
               for g in ['Normal', 'PoorPlaque']]
    ax.legend(handles=handles, fontsize=7, frameon=False,
              loc='center left', bbox_to_anchor=(1.02, 0.5),
              handlelength=1.0, handleheight=0.8)

    if panel_label:
        ax.text(-0.22, 1.08, panel_label,
                transform=ax.transAxes,
                fontsize=13, fontweight='bold', va='top')
    sns.despine(ax=ax)


# ════════════════════════════════════════════════════════════════
# Panel C: volcano draw function (accepts single ax)
# ════════════════════════════════════════════════════════════════
def _compute_da(abund, groups, sn, sp, seed):
    n_total = abund.shape[0]
    n_samp  = abund.shape[1]
    prev    = (abund > 0).sum(axis=1) / n_samp
    mn      = abund[sn].mean(axis=1)
    mp      = abund[sp].mean(axis=1)
    keep    = (prev >= 0.30) & ((mn > 0.01) | (mp > 0.01))
    abund_f = abund.loc[keep]
    n_filt  = abund_f.shape[0]
    logging.info(f"  Pre-filter: {n_total} -> {n_filt} species")

    species = abund_f.index.tolist()
    p_wil   = np.array([
        stats.mannwhitneyu(abund_f.loc[s, sn].values,
                           abund_f.loc[s, sp].values,
                           alternative='two-sided')[1]
        for s in species])
    _, fdr_wil, _, _ = multipletests(p_wil, method='fdr_bh')

    mean_n = abund_f[sn].mean(axis=1).values
    mean_p = abund_f[sp].mean(axis=1).values
    sd_n   = abund_f[sn].std(axis=1).values
    sd_p   = abund_f[sp].std(axis=1).values
    l2fc   = np.log2((mean_p + 0.01) / (mean_n + 0.01))

    logging.info("  Running CLR-MC-Wilcoxon (128 draws) ...")
    p_clrmc = _clr_mc(abund_f, sn, sp, n_boot=128, seed=seed)
    _, fdr_clrmc, _, _ = multipletests(p_clrmc, method='fdr_bh')

    robust    = (fdr_wil < 0.05) & (fdr_clrmc < 0.05)
    robust_up = robust & (l2fc > 0)
    robust_dn = robust & (l2fc < 0)
    logging.info(f"  Robust DA: total={robust.sum()} up={robust_up.sum()} down={robust_dn.sum()}")

    result_df = pd.DataFrame({
        'species':          species,
        'mean_Normal':      mean_n,  'sd_Normal':      sd_n,
        'mean_PoorPlaque':  mean_p,  'sd_PoorPlaque':  sd_p,
        'log2FC':           l2fc,
        'Wilcoxon_p':       p_wil,   'Wilcoxon_FDR':   fdr_wil,
        'CLRMC_p':          p_clrmc, 'CLRMC_FDR':      fdr_clrmc,
        'Robust_DA':        robust,
        'Robust_direction': np.where(robust_up, 'PoorPlaque',
                             np.where(robust_dn, 'Normal', 'NS')),
    })
    return result_df, l2fc, fdr_wil, robust, robust_up, robust_dn, species


def _clr_mc(abund_f, sn, sp, n_boot=128, seed=42):
    rng    = np.random.default_rng(seed)
    mat    = abund_f.T.values
    n_sp   = mat.shape[1]
    scaled = mat * 10 + 0.01
    cols   = list(abund_f.columns)
    ni     = [cols.index(s) for s in sn]
    pi     = [cols.index(s) for s in sp]
    p_mat  = np.zeros((n_boot, n_sp))
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


def draw_panel_c(ax, da_data, panel_label='C'):
    """Draw volcano plot with 4-color scheme onto a provided ax."""
    from matplotlib.lines import Line2D
    _, l2fc, fdr_wil, robust, robust_up, robust_dn, _ = da_data

    robust_up  = robust & (l2fc > 0)
    robust_dn  = robust & (l2fc < 0)
    wilonly_up = (~robust) & (fdr_wil < 0.05) & (l2fc > 0)
    wilonly_dn = (~robust) & (fdr_wil < 0.05) & (l2fc < 0)

    n_rup = int(robust_up.sum());  n_rdn = int(robust_dn.sum())
    n_wup = int(wilonly_up.sum()); n_wdn = int(wilonly_dn.sum())

    y_vals = -np.log10(fdr_wil + 1e-300)

    for mask, c, s, a, z in [
        (~robust & ~wilonly_up & ~wilonly_dn, '#DDDDDD',  4, 0.55, 2),
        (wilonly_up,  '#F4A582',  9, 0.70, 3),
        (wilonly_dn,  '#92C5C5',  9, 0.70, 3),
        (robust_up,   '#B85450', 20, 0.90, 5),
        (robust_dn,   '#2E8B8B', 20, 0.90, 5),
    ]:
        ax.scatter(l2fc[mask], y_vals[mask],
                   c=c, s=s, alpha=a, linewidths=0,
                   zorder=z, rasterized=False)

    ax.axhline(-np.log10(0.05), color='#666666', lw=0.8, ls='--', zorder=1)
    ax.axvline( 1, color='#666666', lw=0.8, ls='--', zorder=1)
    ax.axvline(-1, color='#666666', lw=0.8, ls='--', zorder=1)

    leg_elements = [
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#B85450',
               markersize=5, label=f'Robust UP (n={n_rup})'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#F4A582',
               markersize=4, label=f'Wilcoxon UP (n={n_wup})'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#2E8B8B',
               markersize=5, label=f'Robust DOWN (n={n_rdn})'),
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#92C5C5',
               markersize=4, label=f'Wilcoxon DOWN (n={n_wdn})'),
    ]
    ax.legend(handles=leg_elements, fontsize=5.5, frameon=True,
              loc='lower right', framealpha=0.85, edgecolor='#CCCCCC',
              handletextpad=0.3, borderpad=0.5)

    ax.set_xlabel(r'$\log_2$FC (PoorPlaque vs Normal)', fontsize=10)
    ax.set_ylabel(r'$-\log_{10}$(Wilcoxon FDR)', fontsize=10)
    if panel_label:
        ax.text(-0.22, 1.08, panel_label,
                transform=ax.transAxes,
                fontsize=13, fontweight='bold', va='top')
    sns.despine(ax=ax)


def draw_panel_c_qval(ax, da_data, panel_label=None):
    """q-value version of volcano."""
    _, l2fc, fdr_wil, robust, robust_up, robust_dn, _ = da_data
    n_up = int(robust_up.sum())
    n_dn = int(robust_dn.sum())

    y_vals = -np.log10(fdr_wil + 1e-300)
    colors = np.where(robust_up, '#B85450',
             np.where(robust_dn, '#2E8B8B',
             np.where((fdr_wil < 0.05) & ~robust, '#AAAAAA', '#DDDDDD')))
    sizes  = np.where(robust, 16, 6)

    ax.scatter(l2fc, y_vals, c=colors, s=sizes, alpha=0.75,
               linewidths=0, zorder=3, rasterized=False)
    ax.axhline(-np.log10(0.05), color='#666666', lw=0.8, ls='--', zorder=1)
    ax.axvline( 1, color='#666666', lw=0.8, ls='--', zorder=1)
    ax.axvline(-1, color='#666666', lw=0.8, ls='--', zorder=1)

    ax.text(0.97, 0.97, f'[UP] {n_up} enriched\n(PoorPlaque)',
            transform=ax.transAxes, ha='right', va='top', fontsize=6.5,
            color='#B85450')
    ax.text(0.03, 0.97, f'[DOWN] {n_dn} enriched\n(Normal)',
            transform=ax.transAxes, ha='left', va='top', fontsize=6.5,
            color='#2E8B8B')

    ax.set_xlabel(r'$\log_2$FC (PoorPlaque vs Normal)', fontsize=10)
    ax.set_ylabel(r'$-\log_{10}$(BH $q$-value)', fontsize=10)
    if panel_label:
        ax.text(-0.22, 1.08, panel_label,
                transform=ax.transAxes, fontsize=13, fontweight='bold', va='top')
    sns.despine(ax=ax)


# ════════════════════════════════════════════════════════════════
# Panel D: 10 selected DA species boxplots (draw into grid of axes)
# ════════════════════════════════════════════════════════════════
def draw_panel_d_into_gridspec(fig, gs_slot, abund, sn, sp,
                               da_result, seed,
                               row_label_d='D', row_label_e='E'):
    """
    Draw 2 rows x 5 cols of DA species boxplots into gs_slot.
    gs_slot should be a SubplotSpec that spans the desired area.
    Returns the inner GridSpec object.
    """
    up_sel = [s for s in MAIN_UP   if s in abund.index]
    dn_sel = [s for s in MAIN_DOWN if s in abund.index]

    gs_inner = gridspec.GridSpecFromSubplotSpec(
        2, 5, subplot_spec=gs_slot,
        hspace=0.90, wspace=0.55)

    for row, (sp_list, direction, rl) in enumerate([
        (up_sel,  'PoorPlaque', row_label_d),
        (dn_sel,  'Normal',     row_label_e),
    ]):
        t_color = '#B85450' if direction == 'PoorPlaque' else '#2E8B8B'
        for col, species in enumerate(sp_list):
            ax  = fig.add_subplot(gs_inner[row, col])
            vn  = abund.loc[species, sn].values.astype(float)
            vp  = abund.loc[species, sp].values.astype(float)
            row_info = da_result[da_result['species'] == species]
            l2fc = row_info['log2FC'].values[0] if len(row_info) else np.nan
            _, pv = stats.mannwhitneyu(vn, vp, alternative='two-sided')
            sp_label  = species.replace('_', ' ')
            full_title = f"{sp_label}\n(log\u2082FC = {l2fc:+.2f})"

            _boxplot_on_ax(ax, vn, vp,
                           ylabel='Rel. abund. (%)',
                           p_val=pv, seed=seed,
                           title=full_title,
                           title_color=t_color,
                           fontsize_title=6.0,
                           panel_label=rl if col == 0 else None)

    return gs_inner


# ════════════════════════════════════════════════════════════════
# Standalone panel wrappers
# ════════════════════════════════════════════════════════════════
def save_panel_a(alpha_data, out_dir, dpi, seed):
    sns.set_style('ticks')
    fig, axes = plt.subplots(1, 2, figsize=(3.5, 2.56))
    draw_panel_a(axes[0], axes[1], alpha_data, seed, panel_label='A')
    plt.tight_layout(pad=0.8)
    save_fig(fig, os.path.join(out_dir, 'panels', 'Fig1_PanelA'), dpi)
    plt.close(fig)


def save_panel_b(beta_data, out_dir, dpi):
    sns.set_style('ticks')
    fig, ax = plt.subplots(figsize=(3.0, 2.36))
    draw_panel_b(ax, beta_data, panel_label='B')
    plt.tight_layout(pad=0.8)
    save_fig(fig, os.path.join(out_dir, 'panels', 'Fig1_PanelB'), dpi)
    plt.close(fig)


def save_panel_c(da_data, out_dir, dpi):
    sns.set_style('ticks')
    fig, ax = plt.subplots(figsize=(2.36, 2.56))
    draw_panel_c(ax, da_data, panel_label='C')
    plt.tight_layout(pad=0.8)
    save_fig(fig, os.path.join(out_dir, 'panels', 'Fig1_PanelC'), dpi)
    plt.close(fig)

    # q-value version
    fig2, ax2 = plt.subplots(figsize=(2.36, 2.56))
    draw_panel_c_qval(ax2, da_data, panel_label='C')
    plt.tight_layout(pad=0.8)
    save_fig(fig2, os.path.join(out_dir, 'panels', 'Fig1_PanelC_qval'), dpi)
    plt.close(fig2)


def save_panel_d(abund, sn, sp, da_result, out_dir, dpi, seed):
    sns.set_style('ticks')
    fig = plt.figure(figsize=(7.09, 4.10))
    gs  = gridspec.GridSpec(1, 1, figure=fig,
                            left=0.07, right=0.98,
                            top=0.90, bottom=0.12)
    draw_panel_d_into_gridspec(fig, gs[0, 0], abund, sn, sp,
                               da_result, seed)
    # row labels
    fig.text(0.01, 0.78, 'Enriched in\nPoorPlaque',
             fontsize=7, fontweight='bold', color='#B85450',
             va='center', ha='left', rotation=90)
    fig.text(0.01, 0.30, 'Enriched in\nNormal',
             fontsize=7, fontweight='bold', color='#2E8B8B',
             va='center', ha='left', rotation=90)
    save_fig(fig, os.path.join(out_dir, 'panels', 'Fig1_PanelD'), dpi)
    plt.close(fig)


# ════════════════════════════════════════════════════════════════
# Supplementary DA boxplots
# ════════════════════════════════════════════════════════════════
def save_supp_da(abund, sn, sp, da_result, out_dir, dpi, seed):
    robust_df = da_result[da_result['Robust_DA']].copy()
    supp_df   = robust_df[~robust_df['species'].isin(MAIN_SPECIES)].copy()
    supp_up   = supp_df[supp_df['Robust_direction'] == 'PoorPlaque'].sort_values('log2FC', ascending=False)
    supp_dn   = supp_df[supp_df['Robust_direction'] == 'Normal'].sort_values('log2FC')
    all_supp  = pd.concat([supp_up, supp_dn]).reset_index(drop=True)
    n_sp      = len(all_supp)
    logging.info(f"  Supplementary DA: {n_sp} species ({len(supp_up)} up + {len(supp_dn)} down)")

    ncols = 4
    nrows = math.ceil(n_sp / ncols)

    sns.set_style('ticks')
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 1.55, nrows * 2.10))
    axes = np.array(axes).flatten()

    for idx, (_, row) in enumerate(all_supp.iterrows()):
        ax        = axes[idx]
        species   = row['species']
        direction = row['Robust_direction']
        t_color   = '#B85450' if direction == 'PoorPlaque' else '#2E8B8B'
        vn  = abund.loc[species, sn].values.astype(float)
        vp  = abund.loc[species, sp].values.astype(float)
        _, pv = stats.mannwhitneyu(vn, vp, alternative='two-sided')
        l2fc = row['log2FC']
        full_title = f"{species.replace('_', ' ')}\n(log\u2082FC = {l2fc:+.2f})"
        _boxplot_on_ax(ax, vn, vp,
                       ylabel='Rel. abund. (%)', p_val=pv, seed=seed,
                       title=full_title, title_color=t_color)

    for idx in range(n_sp, len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle('Figure S1: Remaining Differentially Abundant Species',
                 fontsize=9, fontweight='bold', y=1.01)
    plt.tight_layout(pad=0.6, h_pad=1.8, w_pad=0.7)
    save_fig(fig, os.path.join(out_dir, 'FigS1_DA_supplementary'), dpi)
    plt.close(fig)


# ════════════════════════════════════════════════════════════════
# Composite Figure 1  —  100% vector, no PNG embedding
# ════════════════════════════════════════════════════════════════
def save_composite(alpha_data, beta_data, da_data, da_result,
                   abund, sn, sp, out_dir, dpi, seed):
    """
    Full Figure 1 drawn directly — no imread/imshow.
    Layout:
      Row 0 (40%): [Panel A: 2 sub-axes] | [Panel B] | [Panel C]
      Row 1 (60%): [Panel D: 2 rows x 5 cols, full width]
    """
    logging.info("[Composite] Assembling Figure 1 (fully vector) ...")
    sns.set_style('ticks')

    fig = plt.figure(figsize=(7.09, 5.71))  # 180mm x 145mm

    # outer: 2 rows
    gs_outer = gridspec.GridSpec(
        2, 1, figure=fig,
        height_ratios=[1, 1.55],
        hspace=0.10,
        left=0.07, right=0.97,
        top=0.97, bottom=0.05)

    # ── top row: A (2 subaxes) + B + C ──────────────────────────
    gs_top = gridspec.GridSpecFromSubplotSpec(
        1, 3, subplot_spec=gs_outer[0],
        width_ratios=[2, 1.6, 1.4],
        wspace=0.55)

    # Panel A: 2 sub-axes inside gs_top[0]
    gs_a = gridspec.GridSpecFromSubplotSpec(
        1, 2, subplot_spec=gs_top[0, 0], wspace=0.55)
    ax_obs = fig.add_subplot(gs_a[0, 0])
    ax_sh  = fig.add_subplot(gs_a[0, 1])
    draw_panel_a(ax_obs, ax_sh, alpha_data, seed, panel_label='A')

    # Panel B
    ax_b = fig.add_subplot(gs_top[0, 1])
    draw_panel_b(ax_b, beta_data, panel_label='B')

    # Panel C
    ax_c = fig.add_subplot(gs_top[0, 2])
    draw_panel_c(ax_c, da_data, panel_label='C')

    # ── bottom row: Panel D full width ──────────────────────────
    draw_panel_d_into_gridspec(fig, gs_outer[1],
                               abund, sn, sp, da_result, seed,
                               row_label_d='D', row_label_e='E')

    save_fig(fig, os.path.join(out_dir, 'Figure1_combined'), dpi)
    plt.close(fig)
    logging.info("[Composite] Done.")


# ════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════
def parse_args():
    p = argparse.ArgumentParser(
        description='M1 Figure 1 v3 — 100% vector SVG',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    p.add_argument('--abund',  '-a', required=True)
    p.add_argument('--group',  '-g', required=True)
    p.add_argument('--output', '-o', required=True)
    p.add_argument('--seed',   type=int, default=42)
    p.add_argument('--dpi',    type=int, default=300)
    p.add_argument('--verbose','-v', action='store_true')
    return p.parse_args()


def main():
    args = parse_args()
    for sub in ['panels', 'tables']:
        os.makedirs(os.path.join(args.output, sub), exist_ok=True)

    log_path = os.path.join(args.output, 'analysis_log_Fig1.txt')
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.StreamHandler(sys.stdout),
                  logging.FileHandler(log_path, mode='w')])

    np.random.seed(args.seed)
    logging.info(f"[START] P2026-OralMB Figure 1 v3 | seed={args.seed} | dpi={args.dpi}")

    abund, groups, sn, sp = load_data(args.abund, args.group)
    logging.info(f"[DATA] {abund.shape[0]} x {abund.shape[1]} | "
                 f"Normal n={len(sn)}, PoorPlaque n={len(sp)}")

    # ── compute ──────────────────────────────────────────────────
    logging.info("[Panel A] Computing alpha diversity ...")
    alpha_df, obs_data, sh_data = _compute_alpha(abund, groups, sn, sp)
    alpha_data = obs_data, sh_data
    logging.info(f"  Observed p={obs_data[2]:.4g}  Shannon p={sh_data[2]:.4g}")

    logging.info("[Panel B] Computing beta diversity ...")
    beta_data = _compute_beta(abund, groups, sn, sp)

    logging.info("[Panel C] Differential abundance ...")
    da_result, l2fc, fdr_wil, robust, robust_up, robust_dn, species = \
        _compute_da(abund, groups, sn, sp, args.seed)
    da_data = (da_result, l2fc, fdr_wil, robust, robust_up, robust_dn, species)

    # ── save tables ──────────────────────────────────────────────
    alpha_df.to_csv(os.path.join(args.output, 'tables', 'Table_PanelA_alpha.csv'), index=False)
    beta_data[-1].to_csv(os.path.join(args.output, 'tables', 'Table_PanelB_beta_stats.csv'), index=False)
    da_result.to_csv(os.path.join(args.output, 'tables', 'Table_S_differential_species.csv'), index=False)
    logging.info("  [SAVED] All tables")

    # ── individual panels ────────────────────────────────────────
    save_panel_a(alpha_data, args.output, args.dpi, args.seed)
    save_panel_b(beta_data,  args.output, args.dpi)
    save_panel_c(da_data,    args.output, args.dpi)
    save_panel_d(abund, sn, sp, da_result, args.output, args.dpi, args.seed)

    # ── supplementary ────────────────────────────────────────────
    save_supp_da(abund, sn, sp, da_result, args.output, args.dpi, args.seed)

    # ── composite (fully vector) ─────────────────────────────────
    save_composite(alpha_data, beta_data, da_data, da_result,
                   abund, sn, sp, args.output, args.dpi, args.seed)

    logging.info(f"[DONE] All outputs in: {args.output}")


if __name__ == '__main__':
    main()
