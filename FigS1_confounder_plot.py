"""Packaged for revise/; source: project/review/R1.2/scripts/02_plot_FigS_PERMANOVA_confounder_R2.py
Figures: FigS1
Note: confounder R2 plot
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
"""Plot FigS: PERMANOVA R2 after removing covariates, then comparing groups."""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
IN = os.path.join(ROOT, 'tables', 'Table_PERMANOVA_confounder_sensitivity.tsv')
OUTF = os.path.join(ROOT, 'figures')
OUTT = os.path.join(ROOT, 'tables')

candidates = ['Times New Roman', 'Liberation Serif', 'Nimbus Roman', 'STIX', 'DejaVu Serif']
available = {f.name for f in fm.fontManager.ttflist}
font_name = next((c for c in candidates if c in available), 'DejaVu Serif')

mpl.rcParams.update({
    'font.family': font_name,
    'font.weight': 'normal',
    'axes.labelweight': 'normal',
    'axes.titleweight': 'normal',
    'figure.titleweight': 'normal',
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'font.size': 10,
    'axes.labelsize': 11,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 7.5,
    'axes.spines.right': False,
    'axes.spines.top': False,
    'axes.linewidth': 0.8,
    'legend.frameon': False,
    'savefig.dpi': 600,
})


def sig_label(p):
    if p <= 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return 'ns'


def main():
    df = pd.read_csv(IN, sep='\t')
    keep_models = [
        'Group_only',
        'Group_resid_Brushing',
        'Group_resid_Caries',
        'Group_resid_Demin',
        'Group_resid_Brushing_Caries_Demin',
    ]
    model_labels_data = {
        'Group_only': 'Group only (no covariate removed)',
        'Group_resid_Brushing': 'Remove Daily Brushing Frequency, then compare groups',
        'Group_resid_Caries': 'Remove Decayed Teeth, then compare groups',
        'Group_resid_Demin': 'Remove Enamel Demineralization, then compare groups',
        'Group_resid_Brushing_Caries_Demin': (
            'Remove Daily Brushing Frequency + Decayed Teeth + Enamel Demineralization, then compare groups'
        ),
    }
    model_labels_legend = {
        'Group_only': 'Group only (no covariate removed)',
        'Group_resid_Brushing': 'Remove Daily Brushing Frequency,\nthen compare groups',
        'Group_resid_Caries': 'Remove Decayed Teeth,\nthen compare groups',
        'Group_resid_Demin': 'Remove Enamel Demineralization,\nthen compare groups',
        'Group_resid_Brushing_Caries_Demin': (
            'Remove Daily Brushing Frequency\n+ Decayed Teeth\n+ Enamel Demineralization,\nthen compare groups'
        ),
    }
    layer_labels = {'species': 'Species', 'KO': 'KO', 'metabolite': 'Metabolite'}
    layer_order = ['species', 'KO', 'metabolite']

    plot = df[df['model'].isin(keep_models)].copy()
    plot['layer_label'] = plot['layer'].map(layer_labels)
    plot['model_label'] = plot['model'].map(model_labels_data)
    plot['R2'] = plot['R2_from_F']
    plot['sig'] = plot['P'].map(sig_label)
    plot['layer_ord'] = plot['layer'].map({k: i for i, k in enumerate(layer_order)})
    plot['model_ord'] = plot['model'].map({k: i for i, k in enumerate(keep_models)})
    plot = plot.sort_values(['layer_ord', 'model_ord'])
    export = plot[['layer', 'layer_label', 'model', 'model_label', 'metric', 'n',
                   'pseudo_F', 'R2', 'P', 'sig', 'note']]
    export.to_csv(os.path.join(OUTT, 'FigS_PERMANOVA_confounder_R2_plotdata.tsv'),
                  sep='\t', index=False)
    export.to_csv(os.path.join(OUTT, 'FigS_PERMANOVA_confounder_R2_plotdata.csv'),
                  index=False)

    colors = {
        'Group_only': '#1f4e79',
        'Group_resid_Brushing': '#5b8fbc',
        'Group_resid_Caries': '#7ba3c9',
        'Group_resid_Demin': '#9bb8d4',
        'Group_resid_Brushing_Caries_Demin': '#b5674a',
    }

    fig, ax = plt.subplots(figsize=(6.8, 3.4))
    x = np.arange(len(layer_order))
    n_models = len(keep_models)
    width = 0.15
    offsets = (np.arange(n_models) - (n_models - 1) / 2) * width

    for i, model in enumerate(keep_models):
        ys, sigs = [], []
        for layer in layer_order:
            row = plot[(plot['layer'] == layer) & (plot['model'] == model)].iloc[0]
            ys.append(row['R2'])
            sigs.append(row['sig'])
        ax.bar(x + offsets[i], ys, width=width * 0.92, color=colors[model],
               edgecolor='none', label=model_labels_legend[model], zorder=3)
        for xi, y, sig in zip(x + offsets[i], ys, sigs):
            ax.text(xi, y + 0.008, sig, ha='center', va='bottom', fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels([layer_labels[l] for l in layer_order])
    ax.set_ylabel(r'PERMANOVA $R^{2}$')
    ax.set_ylim(0, 0.40)
    ax.set_xlim(-0.45, 2.45)
    ax.legend(loc='upper right', ncol=1, handlelength=1.0, borderaxespad=0.2,
              labelspacing=0.45, borderpad=0.25, fontsize=7.5)
    fig.tight_layout(pad=0.3)

    stem = os.path.join(OUTF, 'FigS_PERMANOVA_confounder_R2')
    for ext in ['pdf', 'svg', 'png', 'jpg']:
        kw = dict(bbox_inches='tight', pad_inches=0.05)
        if ext in ('png', 'jpg'):
            kw['dpi'] = 600
        fig.savefig(f'{stem}.{ext}', **kw)
    plt.close(fig)
    print('Saved', stem, 'font=', font_name)


if __name__ == '__main__':
    main()
