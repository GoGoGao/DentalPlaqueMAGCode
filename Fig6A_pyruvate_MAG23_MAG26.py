"""Fig6A: Opposite-signed MAG23/MAG26 associations with pyruvate.
Packaged for revise/manuscript.md Figure 6A.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PKG_ROOT = Path(__file__).resolve().parents[1]
SHARED = PKG_ROOT / 'shared_input'
OUTDIR = PKG_ROOT / 'output' / 'Fig6A_pyruvate_MAG23_MAG26'
OUTDIR.mkdir(parents=True, exist_ok=True)

COL_G = {'Normal': '#2E8B8B', 'PoorPlaque': '#B85450'}
COL_POS, COL_NEG = '#C44E52', '#4C72B0'
STEM = 'Fig_Pyruvate_MAG23_pos_MAG26_neg_scatter'
PAIRS = [
    dict(MAG='MAG23', metab_id='metab_26768'),
    dict(MAG='MAG26', metab_id='metab_26768'),
]

group = pd.read_csv(SHARED/'sample_group.tsv', sep='\t', header=None, names=['sample','group'])
samples = group['sample'].tolist()
gmap = dict(zip(group['sample'], group['group']))
mag_set = pd.read_csv(SHARED/'robust_DA_MAGs.tsv', sep='\t')
mag_set = mag_set[mag_set['Robust_DA'].astype(str).isin(['True','TRUE','1','true'])]
sp_map = {}
for _, r in mag_set.iterrows():
    sp = r.get('species_plain')
    if pd.isna(sp) or str(sp).strip() == '':
        sp = str(r.get('species', 'Unknown')).replace('_', ' ')
    sp_map[r['MAG']] = sp
mag_ab = pd.read_csv(SHARED/'mag_bin_abundance.tsv', sep='\t')
mag_ab = mag_ab.set_index(mag_ab.columns[0])[samples].apply(pd.to_numeric, errors='coerce')
metab_path = SHARED/'metab_abundance_AB_select.tsv'
if not metab_path.exists():
    metab_path = SHARED/'metab_abundance_AB.tsv'
metab = pd.read_csv(metab_path, sep='\t')
if 'metab_id' in metab.columns:
    metab = metab.set_index('metab_id')
else:
    metab = metab.set_index(metab.columns[0])
metab = metab[samples].apply(pd.to_numeric, errors='coerce')

def log10_halfmin(v):
    v = np.asarray(v, float)
    pos = v[v > 0]
    floor = (pos.min() / 2) if len(pos) else 1e-6
    return np.log10(np.where(v > 0, v, floor))

def ols_ci(x, y, nxx=120):
    x = np.asarray(x, float); y = np.asarray(y, float); n = len(x)
    slope, intercept, r, p, se = stats.linregress(x, y)
    xx = np.linspace(x.min(), x.max(), nxx)
    yy = intercept + slope * xx
    resid = y - (intercept + slope * x)
    dof = max(n - 2, 1)
    s2 = np.sum(resid**2) / dof
    xbar = x.mean(); Sxx = np.sum((x - xbar)**2)
    se_fit = np.sqrt(s2 * (1.0/n + (xx - xbar)**2 / max(Sxx, 1e-12)))
    t = stats.t.ppf(0.975, dof)
    return slope, intercept, r, p, xx, yy, yy - t*se_fit, yy + t*se_fit

def short_sp(sp, mag):
    toks = [t for t in str(sp).replace('_', ' ').split() if t]
    if not toks or toks[0] == 'Unknown':
        return f'Unknown ({mag})'
    return f"{toks[0]} {toks[1] if len(toks)>1 else ''} ({mag})".replace('  ', ' ')

plot_rows, fit_rows = [], []
fig, axes = plt.subplots(2, 1, figsize=(3.4, 4.8), sharex=True)
for ax, pair in zip(axes, PAIRS):
    mag, mid = pair['MAG'], pair['metab_id']
    x = metab.loc[mid, samples].astype(float).values
    y = mag_ab.loc[mag, samples].astype(float).values
    xl = log10_halfmin(x); yl = log10_halfmin(y)
    rho, p_sp = stats.spearmanr(x, y)
    slope, intercept, rpe, ppe, xx, yy, lo, hi = ols_ci(xl, yl)
    col_line = COL_POS if rho >= 0 else COL_NEG
    for s, xv, yv in zip(samples, xl, yl):
        ax.scatter(xv, yv, c=COL_G[gmap[s]], s=28, edgecolors='white', linewidths=0.4, zorder=3)
        plot_rows.append(dict(sample=s, group=gmap[s], MAG=mag, metab_id=mid,
                              pyruvate_log10=float(xv), mag_log10=float(yv),
                              pyruvate_raw=float(metab.loc[mid, s]), mag_raw=float(mag_ab.loc[mag, s])))
    ax.fill_between(xx, lo, hi, color=col_line, alpha=0.18, linewidth=0)
    ax.plot(xx, yy, color=col_line, lw=1.6)
    sp = sp_map.get(mag, mag)
    ax.set_ylabel(short_sp(sp, mag), fontsize=8)
    ax.text(0.02, 0.95, f'ρ={rho:.3f}\nP={p_sp:.1e}', transform=ax.transAxes,
            va='top', ha='left', fontsize=7)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    fit_rows.append(dict(MAG=mag, species=sp, metab_id=mid, metabolite='Pyruvate',
                         spearman_rho=float(rho), spearman_p=float(p_sp),
                         ols_slope_log10=float(slope), ols_intercept_log10=float(intercept),
                         pearson_r_log10=float(rpe), pearson_p_log10=float(ppe)))
axes[-1].set_xlabel('Pyruvate (log10 abundance)', fontsize=8)
fig.tight_layout()
for ext in ('pdf', 'svg', 'png', 'jpg'):
    fig.savefig(OUTDIR/f'{STEM}.{ext}', dpi=600, bbox_inches='tight', pad_inches=0.04)
plt.close(fig)
pd.DataFrame(plot_rows).to_csv(OUTDIR/f'{STEM}_plotdata.tsv', sep='\t', index=False)
pd.DataFrame(plot_rows).to_csv(OUTDIR/f'{STEM}_plotdata.csv', index=False)
pd.DataFrame(fit_rows).to_csv(OUTDIR/f'{STEM}_fit_stats.tsv', sep='\t', index=False)
pd.DataFrame(fit_rows).to_csv(OUTDIR/f'{STEM}_fit_stats.csv', index=False)
print('Wrote', OUTDIR)
