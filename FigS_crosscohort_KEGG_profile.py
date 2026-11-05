"""Packaged for revise/; source: project/02dataset/kegg_lunkuo/plot_kegg_profile.py
Figures: FigS8C
Note: KEGG Level-2 profiles
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
KEGG Functional Profile Plotter
================================
跨队列牙菌斑MAG KEGG功能基因profile分析与可视化

Usage
-----
python3 plot_kegg_profile.py \\
    --kegg   kegg_output.tsv \\
    --meta   group.tsv \\
    --outdir ./results

Arguments
---------
--kegg    KEGG功能频数矩阵 (行=功能类别, 列=MAG, 第一列为功能名称)
--meta    样本元数据 TSV (列: Sample, Group, Dataset)
--outdir  输出目录 (不存在则自动创建)
--font    字体名称 (默认: Liberation Serif)
--dpi     图片分辨率 (默认: 300)

Output
------
<outdir>/
  Fig_A_by_dataset.pdf/.jpg/.svg       -- Panel A: 按Dataset分组
  Fig_B_by_phenotype.pdf/.jpg/.svg     -- Panel B: 按表型分组
  STable_A_by_dataset_stats.tsv        -- Panel A统计数据
  STable_B_by_phenotype_stats.tsv      -- Panel B统计数据
"""

import argparse
import os
import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ─────────────────────────────────────────────────────────────────────
# 常量定义
# ─────────────────────────────────────────────────────────────────────
KEEP_FUNCTIONS = [
    'Carbohydrate metabolism',
    'Energy metabolism',
    'Lipid metabolism',
    'Nucleotide metabolism',
    'Amino acid metabolism',
    'Metabolism of other amino acids',
    'Glycan biosynthesis and metabolism',
    'Metabolism of cofactors and vitamins',
    'Metabolism of terpenoids and polyketides',
    'Biosynthesis of other secondary metabolites',
    'Xenobiotics biodegradation and metabolism',
    'Transcription',
    'Translation',
    'Folding, sorting and degradation',
    'Replication and repair',
    'Membrane transport',
    'Signal transduction',
    'Cell motility',
    'Cellular community - prokaryotes',
    'Environmental adaptation',
    'Drug resistance: antimicrobial',
    'Infectious disease: bacterial',
]

RENAME_MAP = {
    'Carbohydrate metabolism':                      'Carbohydrate metabolism',
    'Energy metabolism':                            'Energy metabolism',
    'Lipid metabolism':                             'Lipid metabolism',
    'Nucleotide metabolism':                        'Nucleotide metabolism',
    'Amino acid metabolism':                        'Amino acid metabolism',
    'Metabolism of other amino acids':              'Other amino acid metabolism',
    'Glycan biosynthesis and metabolism':           'Glycan biosynthesis',
    'Metabolism of cofactors and vitamins':         'Cofactors & vitamins',
    'Metabolism of terpenoids and polyketides':     'Terpenoids & polyketides',
    'Biosynthesis of other secondary metabolites':  'Secondary metabolites',
    'Xenobiotics biodegradation and metabolism':    'Xenobiotics degradation',
    'Transcription':                                'Transcription',
    'Translation':                                  'Translation',
    'Folding, sorting and degradation':             'Folding & degradation',
    'Replication and repair':                       'Replication & repair',
    'Membrane transport':                           'Membrane transport',
    'Signal transduction':                          'Signal transduction',
    'Cell motility':                                'Cell motility',
    'Cellular community - prokaryotes':             'Cellular community',
    'Environmental adaptation':                     'Environmental adaptation',
    'Drug resistance: antimicrobial':               'Antimicrobial resistance',
    'Infectious disease: bacterial':                'Bacterial infectious disease',
}

CATEGORY_MAP = {
    'Carbohydrate metabolism':      'Metabolism',
    'Energy metabolism':            'Metabolism',
    'Lipid metabolism':             'Metabolism',
    'Nucleotide metabolism':        'Metabolism',
    'Amino acid metabolism':        'Metabolism',
    'Other amino acid metabolism':  'Metabolism',
    'Glycan biosynthesis':          'Metabolism',
    'Cofactors & vitamins':         'Metabolism',
    'Terpenoids & polyketides':     'Metabolism',
    'Secondary metabolites':        'Metabolism',
    'Xenobiotics degradation':      'Metabolism',
    'Transcription':                'Genetic Info Processing',
    'Translation':                  'Genetic Info Processing',
    'Folding & degradation':        'Genetic Info Processing',
    'Replication & repair':         'Genetic Info Processing',
    'Membrane transport':           'Environmental Info Processing',
    'Signal transduction':          'Environmental Info Processing',
    'Cell motility':                'Cellular Processes',
    'Cellular community':           'Cellular Processes',
    'Environmental adaptation':     'Cellular Processes',
    'Antimicrobial resistance':     'Human Diseases',
    'Bacterial infectious disease': 'Human Diseases',
}

CAT_ORDER = [
    'Metabolism',
    'Genetic Info Processing',
    'Environmental Info Processing',
    'Cellular Processes',
    'Human Diseases',
]

CAT_COLORS = {
    'Metabolism':                    '#2D7D9A',
    'Genetic Info Processing':       '#E07B39',
    'Environmental Info Processing': '#5BAD6F',
    'Cellular Processes':            '#9B59B6',
    'Human Diseases':                '#C0392B',
}


# ─────────────────────────────────────────────────────────────────────
# 工具函数
# ─────────────────────────────────────────────────────────────────────
def get_sample_id(col: str) -> str:
    """从MAG列名中提取样本ID（去掉 .bin.N 后缀）"""
    parts = col.split('.')
    for i, p in enumerate(parts):
        if p.lower().startswith('bin'):
            return '.'.join(parts[:i])
    return col


def sig_label(p: float) -> str:
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    return ''


def kruskal_posthoc(df_filt, sorted_funcs, groups):
    """
    对每个功能类别做 Kruskal-Wallis 检验，
    并对所有组对做 Mann-Whitney U 事后检验（Bonferroni校正）。
    返回 DataFrame。
    """
    from itertools import combinations
    group_names = list(groups.keys())
    pairs = list(combinations(group_names, 2))
    n_tests = len(sorted_funcs) * len(pairs)

    rows = []
    for f in sorted_funcs:
        vals_all = [df_filt.loc[f, cols].values for cols in groups.values()]
        kw_stat, kw_p = stats.kruskal(*vals_all)

        row = {
            'Function':     f,
            'KEGG_category': CATEGORY_MAP.get(f, ''),
            'KW_statistic': round(kw_stat, 4),
            'KW_pvalue':    kw_p,
            'KW_sig':       sig_label(kw_p),
        }
        # 各组均值 / SD / n
        for gname, gcols in groups.items():
            v = df_filt.loc[f, gcols].values
            row[f'{gname}_n']    = len(gcols)
            row[f'{gname}_mean'] = round(float(np.mean(v)), 3)
            row[f'{gname}_sd']   = round(float(np.std(v, ddof=1)), 3)
            row[f'{gname}_median'] = round(float(np.median(v)), 3)

        # 事后两两检验（原始p，Bonferroni校正）
        for g1, g2 in pairs:
            v1 = df_filt.loc[f, groups[g1]].values
            v2 = df_filt.loc[f, groups[g2]].values
            _, mwu_p = stats.mannwhitneyu(v1, v2, alternative='two-sided')
            bonf_p = min(mwu_p * n_tests, 1.0)
            key = f'{g1}_vs_{g2}'
            row[f'{key}_MWU_p']      = round(mwu_p, 6)
            row[f'{key}_Bonf_p']     = round(bonf_p, 6)
            row[f'{key}_Bonf_sig']   = sig_label(bonf_p)

        rows.append(row)

    return pd.DataFrame(rows)


def draw_panel(ax, df_filt, sorted_funcs, groups, colors, labels,
               title, font_size, show_ylabel=True):
    """
    绘制单个横向分组条形图 Panel。
    返回 cat_mid dict（大类中点位置，用于辅助轴标签）。
    """
    means = {g: df_filt[s].mean(axis=1) for g, s in groups.items()}
    sems  = {g: df_filt[s].sem(axis=1)  for g, s in groups.items()}

    # Kruskal-Wallis 显著性
    sigs = []
    for f in sorted_funcs:
        vals = [df_filt.loc[f, s].values for s in groups.values()]
        _, p = stats.kruskal(*vals)
        sigs.append(sig_label(p))

    n      = len(sorted_funcs)
    y      = np.arange(n)
    bar_h  = 0.24
    n_grp  = len(groups)
    offsets = np.linspace(-(n_grp - 1) * bar_h / 2,
                           (n_grp - 1) * bar_h / 2, n_grp)

    # ── 背景色带 + 大类分界线 ──────────────────────────────────────
    band_col   = ['#F2F7FB', '#FFFFFF']
    band_idx   = 0
    prev_cat   = None
    cat_mid    = {}
    cat_start  = 0
    boundaries = []

    for i, f in enumerate(sorted_funcs):
        cat = CATEGORY_MAP.get(f, '')
        if cat != prev_cat:
            if prev_cat is not None:
                cat_mid[prev_cat] = (cat_start + i - 1) / 2
                boundaries.append(i - 0.5)
            band_idx  = 1 - band_idx
            cat_start = i
            prev_cat  = cat
        # 用 axhspan 替代超界 barh，避免 AI 导入时元素飞出画板
        ax.axhspan(i - 0.5, i + 0.5,
                   facecolor=band_col[band_idx], alpha=1.0, zorder=0)
    cat_mid[prev_cat] = (cat_start + n - 1) / 2

    for b in boundaries:
        ax.axhline(b, color='#CCCCCC', lw=0.8, ls='--', zorder=1)

    # ── 左侧大类色条 ───────────────────────────────────────────────
    bar_offset = -17 * (font_size / 19)
    for i, f in enumerate(sorted_funcs):
        ax.barh(i, 2.0, height=0.92,
                left=bar_offset,
                color=CAT_COLORS[CATEGORY_MAP.get(f, '')],
                alpha=1.0, zorder=3)

    # ── 条形图 + 误差线 ────────────────────────────────────────────
    for gname, offset in zip(groups.keys(), offsets):
        vals = [means[gname][f] for f in sorted_funcs]
        errs = [sems[gname][f]  for f in sorted_funcs]
        ax.barh(y + offset, vals, height=bar_h * 0.88,
                color=colors[gname], alpha=0.82,
                label=labels[gname], zorder=4)
        ax.errorbar(vals, y + offset, xerr=errs,
                    fmt='none', color='#444444',
                    capsize=2.0, lw=0.8, zorder=5)

    # ── 显著性高亮 + 星号 ──────────────────────────────────────────
    max_vals = [
        max(means[g][f] + sems[g][f] for g in groups)
        for f in sorted_funcs
    ]
    for i, (f, sig) in enumerate(zip(sorted_funcs, sigs)):
        if sig:
            # axhspan 替代超界 barh，AI 可正常导入
            ax.axhspan(i - 0.5, i + 0.5,
                       facecolor='#FFF3CD', alpha=0.50, zorder=2)
            ax.text(max_vals[i] + 1.2, i, sig,
                    va='center', ha='left',
                    fontsize=font_size * 1.1,
                    color='#C0392B', fontweight='bold')

    # ── 轴格式 ────────────────────────────────────────────────────
    ax.xaxis.grid(True, alpha=0.3, lw=0.5, color='#BBBBBB', zorder=1)
    ax.axvline(0, color='#888888', lw=0.5, zorder=2)
    ax.set_xlim(bar_offset - 1, max(max_vals) * 1.15)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel('Mean no. of KEGG functional genes per MAG',
                  fontsize=font_size)
    ax.set_title(title, fontsize=font_size * 1.15,
                 fontweight='bold', pad=12, loc='left')

    if show_ylabel:
        ax.set_yticks(y)
        ax.set_yticklabels(sorted_funcs, fontsize=font_size)
        ax.tick_params(axis='y', length=0, pad=26)
    else:
        ax.set_yticks(y)
        ax.set_yticklabels([], fontsize=0)

    ax.tick_params(axis='x', labelsize=font_size * 0.9)

    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_color('#AAAAAA')
        ax.spines[sp].set_linewidth(0.8)

    # ── 图例 ──────────────────────────────────────────────────────
    handles = [
        mpatches.Patch(color=colors[g], alpha=0.82, label=labels[g])
        for g in groups
    ]
    ax.legend(handles=handles, loc='lower right',
              fontsize=font_size * 0.85,
              framealpha=0.95, edgecolor='#CCCCCC',
              title='Group', title_fontsize=font_size * 0.9)

    return cat_mid


def add_category_axis(ax, cat_mid, font_size):
    """在右侧添加大类标签辅助轴"""
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ticks  = [cat_mid[c] for c in CAT_ORDER if c in cat_mid]
    tlabels = [c for c in CAT_ORDER if c in cat_mid]
    ax2.set_yticks(ticks)
    ax2.set_yticklabels(tlabels, fontsize=font_size * 0.85, fontweight='bold')
    for tick, cat in zip(ax2.get_yticklabels(), tlabels):
        tick.set_color(CAT_COLORS[cat])
    ax2.tick_params(axis='y', length=0, pad=5)
    for sp in ax2.spines.values():
        sp.set_visible(False)
    return ax2


def save_figure(fig, base_path, dpi):
    """同时保存 PDF / JPG / SVG 三种格式"""
    fig.savefig(base_path + '.pdf', dpi=dpi, bbox_inches='tight')
    fig.savefig(base_path + '.jpg', dpi=dpi, bbox_inches='tight')
    fig.savefig(base_path + '.svg', bbox_inches='tight')
    print(f'  Saved: {base_path}.pdf / .jpg / .svg')


# ─────────────────────────────────────────────────────────────────────
# 主流程
# ─────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description='KEGG functional profile plotter for dental plaque MAGs'
    )
    parser.add_argument('--kegg',   required=True, help='KEGG频数矩阵 TSV')
    parser.add_argument('--meta',   required=True, help='样本元数据 TSV (Sample/Group/Dataset)')
    parser.add_argument('--outdir', default='./kegg_profile_output', help='输出目录')
    parser.add_argument('--font',   default='Liberation Serif', help='字体名称')
    parser.add_argument('--dpi',    type=int, default=300, help='图片分辨率')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ── 设置字体 ──────────────────────────────────────────────────
    matplotlib.rcParams['font.family'] = args.font
    BASE_FS = 19   # 基础字号（原9.5的2倍）

    # ── 加载数据 ──────────────────────────────────────────────────
    print('[1/5] Loading data...')
    df   = pd.read_csv(args.kegg, sep='\t', index_col=0)
    meta = pd.read_csv(args.meta, sep='\t')

    meta_dict = {
        str(r['Sample']): (str(r['Group']), str(r['Dataset']))
        for _, r in meta.iterrows()
    }
    col_meta = {
        col: meta_dict.get(get_sample_id(col), (None, None))
        for col in df.columns
    }

    n_matched = sum(1 for v in col_meta.values() if v[0] is not None)
    print(f'  MAGs matched to metadata: {n_matched}/{len(df.columns)}')
    if n_matched == 0:
        sys.exit('ERROR: No MAGs matched metadata. Check sample ID format.')

    # ── 按 Dataset 分组 ───────────────────────────────────────────
    datasets = sorted(set(d for _, (g, d) in col_meta.items() if d))
    groups_dataset = {
        d: [c for c, (g, dd) in col_meta.items() if dd == d]
        for d in datasets
    }
    print('[2/5] Dataset groups:')
    for k, v in groups_dataset.items():
        print(f'  {k}: {len(v)} MAGs')

    # ── 按 Phenotype 分组（排除Unknown） ──────────────────────────
    phenotypes = sorted(set(
        g for _, (g, d) in col_meta.items()
        if g and g.lower() not in ('unknown', 'none', 'nan')
    ))
    groups_pheno = {
        p: [c for c, (g, d) in col_meta.items() if g == p]
        for p in phenotypes
    }
    print('[2/5] Phenotype groups (Unknown excluded):')
    for k, v in groups_pheno.items():
        print(f'  {k}: {len(v)} MAGs')

    # ── 过滤 & 重命名功能类别 ─────────────────────────────────────
    available = [f for f in KEEP_FUNCTIONS if f in df.index]
    missing   = [f for f in KEEP_FUNCTIONS if f not in df.index]
    if missing:
        print(f'  Warning: {len(missing)} functions not found in KEGG file: {missing}')
    df_filt = df.loc[available].rename(index=RENAME_MAP)
    funcs   = list(df_filt.index)
    print(f'  Functional categories used: {len(funcs)}')

    # ── 排序（大类顺序 + ThisStudy均值降序） ──────────────────────
    ref_group  = list(groups_dataset.keys())[0]
    means_ref  = df_filt[groups_dataset[ref_group]].mean(axis=1)
    sorted_funcs = []
    for cat in CAT_ORDER:
        cf = sorted(
            [f for f in funcs if CATEGORY_MAP.get(f) == cat],
            key=lambda f: means_ref[f], reverse=True
        )
        sorted_funcs.extend(cf)

    # ── 颜色方案 ──────────────────────────────────────────────────
    palette_d = ['#2D7D9A', '#C0392B', '#27AE60', '#8E44AD', '#E67E22']
    palette_p = {'Normal': '#2D7D9A', 'PoorPlaque': '#C0392B',
                 'Caries': '#E8A838', 'Healthy': '#27AE60'}

    colors_d = {d: palette_d[i % len(palette_d)]
                for i, d in enumerate(datasets)}
    labels_d = {d: f'{d} (n={len(groups_dataset[d])})' for d in datasets}

    colors_p = {p: palette_p.get(p, palette_d[i % len(palette_d)])
                for i, p in enumerate(phenotypes)}
    labels_p = {p: f'{p} (n={len(groups_pheno[p])})' for p in phenotypes}

    # ── 统计检验 & 输出 STable ───────────────────────────────────
    print('[3/5] Running statistical tests...')

    stats_d = kruskal_posthoc(df_filt, sorted_funcs, groups_dataset)
    out_sta = os.path.join(args.outdir, 'STable_A_by_dataset_stats.tsv')
    stats_d.to_csv(out_sta, sep='\t', index=False)
    print(f'  Saved: {out_sta}')

    stats_p = kruskal_posthoc(df_filt, sorted_funcs, groups_pheno)
    out_stb = os.path.join(args.outdir, 'STable_B_by_phenotype_stats.tsv')
    stats_p.to_csv(out_stb, sep='\t', index=False)
    print(f'  Saved: {out_stb}')

    # ── Panel A：按 Dataset ───────────────────────────────────────
    print('[4/5] Plotting Panel A (by dataset)...')
    fig_w = max(16, len(sorted_funcs) * 0.4 + 8)
    fig_h = max(10, len(sorted_funcs) * 0.42)
    figA, axA = plt.subplots(figsize=(fig_w, fig_h))
    figA.patch.set_facecolor('white')

    cat_mid = draw_panel(
        axA, df_filt, sorted_funcs,
        groups_dataset, colors_d, labels_d,
        title='KEGG functional gene profile by dataset',
        font_size=BASE_FS,
        show_ylabel=True
    )
    add_category_axis(axA, cat_mid, BASE_FS)

    figA.text(
        0.5, -0.02,
        'Highlighted rows: significant difference (Kruskal-Wallis)  '
        '* p<0.05   ** p<0.01   *** p<0.001',
        ha='center', fontsize=BASE_FS * 0.75, color='#555555', style='italic'
    )
    figA.tight_layout()
    save_figure(figA, os.path.join(args.outdir, 'Fig_A_by_dataset'), args.dpi)
    plt.close(figA)

    # ── Panel B：按 Phenotype ─────────────────────────────────────
    print('[5/5] Plotting Panel B (by phenotype)...')
    figB, axB = plt.subplots(figsize=(fig_w, fig_h))
    figB.patch.set_facecolor('white')

    cat_mid_p = draw_panel(
        axB, df_filt, sorted_funcs,
        groups_pheno, colors_p, labels_p,
        title='KEGG functional gene profile by phenotype (cross-cohort)',
        font_size=BASE_FS,
        show_ylabel=True
    )
    add_category_axis(axB, cat_mid_p, BASE_FS)

    figB.text(
        0.5, -0.02,
        'Highlighted rows: significant difference (Kruskal-Wallis)  '
        '* p<0.05   ** p<0.01   *** p<0.001\n'
        'Note: PoorPlaque = dental plaque disease (This Study); '
        'Caries = dental caries (PRJNA766357); '
        'Unknown phenotype samples excluded',
        ha='center', fontsize=BASE_FS * 0.75, color='#555555', style='italic'
    )
    figB.tight_layout()
    save_figure(figB, os.path.join(args.outdir, 'Fig_B_by_phenotype'), args.dpi)
    plt.close(figB)

    print('\n✓ All done. Outputs in:', args.outdir)


if __name__ == '__main__':
    main()
