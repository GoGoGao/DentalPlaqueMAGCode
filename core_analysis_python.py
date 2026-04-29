#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Core Analysis Functions for Oral Microbiome Metagenomic Study
=============================================================

Reusable statistical and visualization functions for:
  - Alpha / Beta diversity analysis
  - Differential abundance analysis (Wilcoxon + CLR-MC dual validation)
  - MAG-level differential KO comparison (Fisher exact test)

Dependencies:
  numpy, pandas, scipy, statsmodels, scikit-bio, matplotlib, seaborn
"""

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.spatial.distance import pdist, squareform
from statsmodels.stats.multitest import multipletests

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = 'Times New Roman'


# ══════════════════════════════════════════════════════════════════════════════
#  1. Alpha Diversity
# ══════════════════════════════════════════════════════════════════════════════

def compute_alpha_diversity(abund_df):
    """
    Compute observed species richness and Shannon diversity index.

    Parameters
    ----------
    abund_df : pd.DataFrame
        Species × samples abundance table (rows = species, cols = samples).

    Returns
    -------
    pd.DataFrame
        Columns: sample, Observed_Species, Shannon_H.
    """
    prop = abund_df.div(abund_df.sum(axis=0), axis=1).replace(0, np.nan)
    observed = (abund_df > 0).sum(axis=0)
    shannon = -(prop * np.log(prop)).sum(axis=0, skipna=True)

    return pd.DataFrame({
        'sample': abund_df.columns,
        'Observed_Species': observed.values,
        'Shannon_H': shannon.values,
    })


def compare_alpha(alpha_df, group_map, group_a='Normal', group_b='PoorPlaque'):
    """
    Compare alpha diversity between two groups using Mann-Whitney U test.

    Parameters
    ----------
    alpha_df : pd.DataFrame
        Output from compute_alpha_diversity().
    group_map : dict
        {sample_name: group_label}.
    group_a, group_b : str
        Group labels to compare.

    Returns
    -------
    dict
        Keys: metric -> (U_stat, p_value).
    """
    alpha_df = alpha_df.copy()
    alpha_df['group'] = alpha_df['sample'].map(group_map)

    results = {}
    for metric in ['Observed_Species', 'Shannon_H']:
        vals_a = alpha_df.loc[alpha_df.group == group_a, metric].values
        vals_b = alpha_df.loc[alpha_df.group == group_b, metric].values
        u, p = stats.mannwhitneyu(vals_a, vals_b, alternative='two-sided')
        results[metric] = {'U_statistic': u, 'p_value': p}

    return results


# ══════════════════════════════════════════════════════════════════════════════
#  2. Beta Diversity (PCoA + PERMANOVA)
# ══════════════════════════════════════════════════════════════════════════════

def compute_beta_diversity(abund_df, group_map, metric='braycurtis',
                           n_permutations=9999):
    """
    Compute Bray-Curtis PCoA, PERMANOVA, ANOSIM, and PERMDISP.

    Parameters
    ----------
    abund_df : pd.DataFrame
        Species × samples abundance table.
    group_map : dict
        {sample_name: group_label}.
    metric : str
        Distance metric for scipy.spatial.distance.pdist.
    n_permutations : int
        Number of permutations for PERMANOVA/ANOSIM.

    Returns
    -------
    dict with keys:
        'coords'    : np.ndarray (n_samples, 2) — PCoA coordinates
        'var_explained' : (float, float) — % variance for PC1, PC2
        'permanova' : dict(F, R2, p_value)
        'anosim'    : dict(R, p_value)
        'permdisp'  : dict(F, p_value)
        'sample_ids': list
        'groups'    : np.ndarray
    """
    from skbio.stats.ordination import pcoa as sk_pcoa
    from skbio.stats.distance import DistanceMatrix, permanova, anosim

    sample_ids = [s for s in abund_df.columns if s in group_map]
    mat = abund_df[sample_ids].T.values
    groups = np.array([group_map[s] for s in sample_ids])

    bc = squareform(pdist(mat, metric=metric))
    dm = DistanceMatrix(bc, ids=sample_ids)

    # PCoA
    res = sk_pcoa(dm)
    coords = res.samples[['PC1', 'PC2']].values
    var1 = res.proportion_explained['PC1'] * 100
    var2 = res.proportion_explained['PC2'] * 100

    # PERMANOVA
    perm = permanova(dm, groups, permutations=n_permutations)
    pF = perm['test statistic']
    pP = perm['p-value']
    pR2 = pF / (pF + len(sample_ids) - 2)

    # ANOSIM
    anos = anosim(dm, groups, permutations=n_permutations)

    # PERMDISP (Levene's test on distances to centroids)
    unique_groups = np.unique(groups)
    centroid_dists = {}
    for g in unique_groups:
        idx = np.where(groups == g)[0]
        centroid = coords[idx].mean(axis=0)
        centroid_dists[g] = np.linalg.norm(coords[idx] - centroid, axis=1)
    dist_arrays = [centroid_dists[g] for g in unique_groups]
    levene_F, levene_p = stats.levene(*dist_arrays)

    return {
        'coords': coords,
        'var_explained': (var1, var2),
        'permanova': {'F': pF, 'R2': pR2, 'p_value': pP},
        'anosim': {'R': anos['test statistic'], 'p_value': anos['p-value']},
        'permdisp': {'F': levene_F, 'p_value': levene_p},
        'sample_ids': sample_ids,
        'groups': groups,
    }


# ══════════════════════════════════════════════════════════════════════════════
#  3. Differential Abundance Analysis (Wilcoxon + CLR-MC)
# ══════════════════════════════════════════════════════════════════════════════

def filter_features(abund_df, samples_a, samples_b,
                    min_prevalence=0.30, min_mean_abund=0.01):
    """
    Pre-filter features by prevalence and mean abundance.

    Parameters
    ----------
    abund_df : pd.DataFrame
        Features × samples.
    samples_a, samples_b : list
        Sample IDs for each group.
    min_prevalence : float
        Minimum detection rate across all samples.
    min_mean_abund : float
        Minimum mean abundance in at least one group.

    Returns
    -------
    pd.DataFrame
        Filtered subset.
    """
    all_samples = samples_a + samples_b
    sub = abund_df[all_samples]
    n_total = len(all_samples)

    prev = (sub > 0).sum(axis=1) / n_total
    mean_a = sub[samples_a].mean(axis=1)
    mean_b = sub[samples_b].mean(axis=1)
    keep = (prev >= min_prevalence) & ((mean_a > min_mean_abund) | (mean_b > min_mean_abund))

    return abund_df.loc[keep]


def clr_monte_carlo(abund_df, samples_a, samples_b,
                    n_boot=128, seed=42):
    """
    CLR Monte Carlo test for compositional bias correction.

    For each bootstrap iteration:
      1. Draw Dirichlet-distributed samples to simulate count uncertainty
      2. Apply centered log-ratio (CLR) transformation
      3. Perform Mann-Whitney U test on CLR-transformed values
    Return the median p-value across iterations per feature.

    Parameters
    ----------
    abund_df : pd.DataFrame
        Features × samples (pre-filtered).
    samples_a, samples_b : list
        Sample IDs for each group.
    n_boot : int
        Number of Monte Carlo iterations.
    seed : int
        Random seed.

    Returns
    -------
    np.ndarray
        Median p-value per feature (same order as abund_df.index).
    """
    rng = np.random.default_rng(seed)
    all_samples = samples_a + samples_b
    mat = abund_df[all_samples].T.values  # (n_samples, n_features)
    n_features = mat.shape[1]

    scaled = mat * 10 + 0.01
    idx_a = [all_samples.index(s) for s in samples_a]
    idx_b = [all_samples.index(s) for s in samples_b]

    p_mat = np.zeros((n_boot, n_features))
    for b in range(n_boot):
        sim = np.array([
            rng.dirichlet(scaled[i] + 1e-6) * scaled[i].sum()
            for i in range(scaled.shape[0])
        ])
        sim[sim == 0] = np.finfo(float).tiny
        lg = np.log(sim)
        clr = lg - lg.mean(axis=1, keepdims=True)
        for j in range(n_features):
            _, p = stats.mannwhitneyu(
                clr[idx_a, j], clr[idx_b, j], alternative='two-sided'
            )
            p_mat[b, j] = p

    return np.median(p_mat, axis=0)


def differential_abundance(abund_df, samples_a, samples_b,
                           group_a_name='Normal', group_b_name='PoorPlaque',
                           min_prevalence=0.30, min_mean_abund=0.01,
                           fdr_threshold=0.05, n_boot=128, seed=42,
                           pseudocount=0.01):
    """
    Full differential abundance pipeline: filter → Wilcoxon → CLR-MC → classify.

    Parameters
    ----------
    abund_df : pd.DataFrame
        Features × samples abundance table.
    samples_a, samples_b : list
        Sample IDs for each group.
    group_a_name, group_b_name : str
        Group labels for output.
    fdr_threshold : float
        FDR cutoff for significance.
    n_boot : int
        CLR-MC iterations.
    seed : int
        Random seed.
    pseudocount : float
        Added to means before log2FC computation.

    Returns
    -------
    pd.DataFrame
        One row per tested feature with columns:
        feature, mean_A, mean_B, log2FC, Wilcoxon_p, Wilcoxon_FDR,
        CLRMC_p, CLRMC_FDR, Robust_DA, Robust_direction.
    """
    filtered = filter_features(
        abund_df, samples_a, samples_b,
        min_prevalence=min_prevalence, min_mean_abund=min_mean_abund
    )

    features = filtered.index.tolist()
    mean_a = filtered[samples_a].mean(axis=1).values
    mean_b = filtered[samples_b].mean(axis=1).values
    log2fc = np.log2((mean_b + pseudocount) / (mean_a + pseudocount))

    # Wilcoxon rank-sum test
    p_wil = np.array([
        stats.mannwhitneyu(
            filtered.loc[f, samples_a].values,
            filtered.loc[f, samples_b].values,
            alternative='two-sided'
        )[1] for f in features
    ])
    _, fdr_wil, _, _ = multipletests(p_wil, method='fdr_bh')

    # CLR Monte Carlo
    p_clrmc = clr_monte_carlo(filtered, samples_a, samples_b,
                              n_boot=n_boot, seed=seed)
    _, fdr_clrmc, _, _ = multipletests(p_clrmc, method='fdr_bh')

    # Classify
    robust = (fdr_wil < fdr_threshold) & (fdr_clrmc < fdr_threshold)
    robust_up = robust & (log2fc > 0)
    robust_dn = robust & (log2fc < 0)

    result = pd.DataFrame({
        'feature': features,
        f'mean_{group_a_name}': mean_a,
        f'mean_{group_b_name}': mean_b,
        'log2FC': log2fc,
        'Wilcoxon_p': p_wil,
        'Wilcoxon_FDR': fdr_wil,
        'CLRMC_p': p_clrmc,
        'CLRMC_FDR': fdr_clrmc,
        'Robust_DA': robust,
        'Robust_direction': np.where(
            robust_up, group_b_name,
            np.where(robust_dn, group_a_name, 'NS')
        ),
    })

    return result


# ══════════════════════════════════════════════════════════════════════════════
#  4. MAG-Level Differential KO (Fisher Exact Test)
# ══════════════════════════════════════════════════════════════════════════════

def fisher_differential_features(binary_matrix, mags_a, mags_b,
                                 group_a_name='GroupA',
                                 group_b_name='GroupB',
                                 p_threshold=0.05):
    """
    Fisher exact test for differential feature prevalence between two MAG groups.

    Parameters
    ----------
    binary_matrix : pd.DataFrame
        Features × MAGs binary (0/1) matrix.
    mags_a, mags_b : list
        MAG IDs for each group.
    group_a_name, group_b_name : str
        Group labels.
    p_threshold : float
        Nominal p-value cutoff for reporting.

    Returns
    -------
    pd.DataFrame
        Columns: feature, prevalence_A, prevalence_B, mean_count_A,
        mean_count_B, odds_ratio, p_value, enriched_in.
    """
    avail_a = [m for m in mags_a if m in binary_matrix.columns]
    avail_b = [m for m in mags_b if m in binary_matrix.columns]
    n_a, n_b = len(avail_a), len(avail_b)

    results = []
    for feat in binary_matrix.index:
        a_present = int((binary_matrix.loc[feat, avail_a] > 0).sum())
        b_present = int((binary_matrix.loc[feat, avail_b] > 0).sum())
        a_absent = n_a - a_present
        b_absent = n_b - b_present

        table = np.array([[a_present, a_absent],
                          [b_present, b_absent]])
        odds, pval = stats.fisher_exact(table, alternative='two-sided')

        prev_a = a_present / max(n_a, 1)
        prev_b = b_present / max(n_b, 1)

        enriched = group_a_name if prev_a > prev_b else group_b_name

        results.append({
            'feature': feat,
            f'{group_a_name}_prevalence': prev_a,
            f'{group_b_name}_prevalence': prev_b,
            'odds_ratio': odds,
            'p_value': pval,
            'enriched_in': enriched,
        })

    result_df = pd.DataFrame(results)
    result_df['padj'] = multipletests(result_df['p_value'], method='fdr_bh')[1]

    return result_df


# ══════════════════════════════════════════════════════════════════════════════
#  5. Utility Functions
# ══════════════════════════════════════════════════════════════════════════════

def save_figure(fig, prefix, dpi=300):
    """Save figure in PDF, SVG, and JPG formats simultaneously."""
    fig.savefig(f"{prefix}.pdf", format='pdf', bbox_inches='tight')
    fig.savefig(f"{prefix}.svg", format='svg', bbox_inches='tight')
    fig.savefig(f"{prefix}.jpg", format='jpg', dpi=dpi, bbox_inches='tight')


def format_pvalue(p):
    """Format p-value for display."""
    if p < 0.001:
        return f"p = {p:.2e}"
    else:
        return f"p = {p:.4f}"
