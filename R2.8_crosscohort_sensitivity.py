"""Packaged for revise/; source: project/review/R2.8/scripts/01_R2_8_sensitivity.py
Figures: FigS8;FigS9
Note: cross-cohort sensitivity
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
"""R2.8 sensitivity: equal-n MAG subsample, HQ/near-complete subsets,
key-KO direction consistency, cross-cohort fine-feature concordance,
depth/recovery proxies. seed=20261028
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr

SEED = 20261028
B_SUB = 1000
OUT = Path(str(PKG_ROOT) + "/review/R2.8")
TAB = OUT / "tables"
FIG = OUT / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

ROOT = PKG_ROOT
KO_PATH = ROOT / "01this_study_dataset/sgb/compare_genomics/total.KO.xls"
RD_PATH = ROOT / "01this_study_dataset/sgb/metab_corr/paper_result/00_RobustDA_MAG_set/tables/RobustDA_MAGs.tsv"
QC_PATH = ROOT / "01this_study_dataset/sgb/tree/plot.input2.tsv"
ABUND_PATH = ROOT / "01this_study_dataset/input/bin_abundance_table.tsv"
ENV_PATH = ROOT / "01this_study_dataset/input/env.tsv"
DET_PATH = ROOT / "review/R1.6/tables/Table_R1.6_summary_counts.tsv"
CROSS_PATH = ROOT / "02dataset/analyse/results/tables/Table_cross_cohort_DA_comparison.csv"
FISHER_PATH = (
    ROOT
    / "01this_study_dataset/sgb/metab_corr/paper_result/03_RobustDA_diffKO_heatmap/tables/fisher_KO_significant.tsv"
)

# Highlighted KOs from Fig.4E / Fig.6D (manuscript)
KEY_KOS = {
    "K00174": "korA/oorA",
    "K00175": "korB/oorB",
    "K23351": "mmdC",
    "K00244": "frdA",
    "K00245": "frdB",
    "K00246": "frdC",
    "K00247": "frdD",
    "K11258": "ilvM",
    "K01652": "ilvB/G/I",
    "K01653": "ilvH/N",
    "K00169": "porA",  # por family if present
    "K00170": "porB",
    "K00625": "pta",
    "K03386": "ahpC",
    "K02238": "comEB",
    "K01938": "purK",
    "K00548": "metX",
    "K01798": "dapD",
    "K00161": "aceE",
}


def dual_write(df: pd.DataFrame, stem: str) -> None:
    df.to_csv(TAB / f"{stem}.tsv", sep="\t", index=False)
    df.to_csv(TAB / f"{stem}.csv", index=False)


def is_hq(c, cont):
    return (c >= 90) & (cont < 5)


def is_near(c, cont):
    # stricter near-complete draft used for sensitivity
    return (c >= 95) & (cont < 5)


def ko_union_stats(mat: pd.DataFrame, mags_a: list, mags_b: list) -> dict:
    a = [m for m in mags_a if m in mat.columns]
    b = [m for m in mags_b if m in mat.columns]
    set_a = set(mat.index[mat[a].gt(0).any(axis=1)]) if a else set()
    set_b = set(mat.index[mat[b].gt(0).any(axis=1)]) if b else set()
    shared = set_a & set_b
    only_a = set_a - set_b
    only_b = set_b - set_a
    union = shared | only_a | only_b
    n_u = len(union) or 1
    return {
        "n_A": len(a),
        "n_B": len(b),
        "n_KO_A": len(set_a),
        "n_KO_B": len(set_b),
        "n_shared": len(shared),
        "n_A_only": len(only_a),
        "n_B_only": len(only_b),
        "n_union": len(union),
        "pct_shared": 100 * len(shared) / n_u,
        "pct_group_specific": 100 * (len(only_a) + len(only_b)) / n_u,
        "pct_A_only": 100 * len(only_a) / n_u,
        "pct_B_only": 100 * len(only_b) / n_u,
    }


def prop_delta(mat: pd.DataFrame, mags_n: list, mags_p: list, ko: str) -> float | None:
    if ko not in mat.index:
        return None
    n = [m for m in mags_n if m in mat.columns]
    p = [m for m in mags_p if m in mat.columns]
    if not n or not p:
        return None
    pn = (mat.loc[ko, n] > 0).mean()
    pp = (mat.loc[ko, p] > 0).mean()
    return float(pp - pn)


def main() -> None:
    rng = np.random.default_rng(SEED)

    ko_raw = pd.read_csv(KO_PATH, sep="\t", index_col=0)
    # Matrix rownames are "Kxxxxx:gene"; collapse to KO id (max copy if duplicates)
    ko_raw = ko_raw.copy()
    ko_raw.index = ko_raw.index.astype(str).str.split(":").str[0]
    ko = ko_raw.groupby(level=0).max()

    rd = pd.read_csv(RD_PATH, sep="\t")
    qc = pd.read_csv(QC_PATH, sep="\t")
    m = rd.merge(qc, on="genome", how="left", suffixes=("", "_qc"))
    # enrichment column
    if "enrichment" not in m.columns:
        m["enrichment"] = m["Robust_direction"]

    normal = m.loc[m["enrichment"] == "Normal", "genome"].tolist()
    poor = m.loc[m["enrichment"] == "PoorPlaque", "genome"].tolist()

    # ---- 1. quality definitions + counts ----
    qc_all = qc.copy()
    qc_all["HQ_MIMAG_CheckM"] = is_hq(qc_all["completeness"], qc_all["contamination"])
    qc_all["near_complete_ge95"] = is_near(qc_all["completeness"], qc_all["contamination"])
    qc_all["medium_or_better"] = (qc_all["completeness"] >= 50) & (qc_all["contamination"] < 10)
    dual_write(
        qc_all.rename(columns={"genome": "MAG"}),
        "Table_R2.8_MAG_quality_definitions_all125",
    )

    m["HQ"] = is_hq(m["completeness"], m["contamination"])
    m["near95"] = is_near(m["completeness"], m["contamination"])
    dual_write(
        m[
            [
                "genome",
                "enrichment",
                "species",
                "completeness",
                "contamination",
                "HQ",
                "near95",
            ]
        ].rename(columns={"genome": "MAG"}),
        "Table_R2.8_RobustDA_MAG_quality",
    )

    def_counts = pd.DataFrame(
        [
            {
                "definition": "medium_or_better (completeness>=50%, contamination<10%)",
                "scope": "all_125_MAGs",
                "n": int(qc_all["medium_or_better"].sum()),
            },
            {
                "definition": "high_quality (completeness>=90%, contamination<5%; manuscript/MIMAG CheckM)",
                "scope": "all_125_MAGs",
                "n": int(qc_all["HQ_MIMAG_CheckM"].sum()),
            },
            {
                "definition": "near_complete_sensitivity (completeness>=95%, contamination<5%)",
                "scope": "all_125_MAGs",
                "n": int(qc_all["near_complete_ge95"].sum()),
            },
            {
                "definition": "high_quality",
                "scope": "Robust_DA_Normal",
                "n": int(m.loc[m["enrichment"] == "Normal", "HQ"].sum()),
            },
            {
                "definition": "high_quality",
                "scope": "Robust_DA_PoorPlaque",
                "n": int(m.loc[m["enrichment"] == "PoorPlaque", "HQ"].sum()),
            },
            {
                "definition": "near_complete_ge95",
                "scope": "Robust_DA_Normal",
                "n": int(m.loc[m["enrichment"] == "Normal", "near95"].sum()),
            },
            {
                "definition": "near_complete_ge95",
                "scope": "Robust_DA_PoorPlaque",
                "n": int(m.loc[m["enrichment"] == "PoorPlaque", "near95"].sum()),
            },
        ]
    )
    dual_write(def_counts, "Table_R2.8_quality_definition_counts")

    # ---- 2. baseline + sensitivity Venn scenarios ----
    scenarios = []
    base = ko_union_stats(ko, normal, poor)
    base.update(
        {
            "scenario": "baseline_RobustDA_10v35",
            "note": "Primary Fig.4C comparison",
        }
    )
    scenarios.append(base)

    # HQ-only (unequal n)
    n_hq = m.loc[(m["enrichment"] == "Normal") & m["HQ"], "genome"].tolist()
    p_hq = m.loc[(m["enrichment"] == "PoorPlaque") & m["HQ"], "genome"].tolist()
    s = ko_union_stats(ko, n_hq, p_hq)
    s.update({"scenario": "HQ_only_unequal", "note": f"n={len(n_hq)} vs {len(p_hq)}"})
    scenarios.append(s)

    # near95-only
    n_nc = m.loc[(m["enrichment"] == "Normal") & m["near95"], "genome"].tolist()
    p_nc = m.loc[(m["enrichment"] == "PoorPlaque") & m["near95"], "genome"].tolist()
    s = ko_union_stats(ko, n_nc, p_nc)
    s.update(
        {
            "scenario": "near95_only_unequal",
            "note": f"n={len(n_nc)} vs {len(p_nc)}; small Normal n",
        }
    )
    scenarios.append(s)

    # equal-n subsample: draw 10 PP, keep all 10 Normal
    sub_rows = []
    for i in range(B_SUB):
        draw = list(rng.choice(poor, size=len(normal), replace=False))
        st = ko_union_stats(ko, normal, draw)
        st["iter"] = i
        sub_rows.append(st)
    sub_df = pd.DataFrame(sub_rows)
    dual_write(sub_df, "Table_R2.8_equalN_subsample_iterations")

    summ = {
        "scenario": "equalN_subsample_10v10",
        "note": f"B={B_SUB}; subsample PoorPlaque to n=10; seed={SEED}",
        "n_A": 10,
        "n_B": 10,
        "n_KO_A_median": float(sub_df["n_KO_A"].median()),
        "n_KO_B_median": float(sub_df["n_KO_B"].median()),
        "n_shared_median": float(sub_df["n_shared"].median()),
        "n_A_only_median": float(sub_df["n_A_only"].median()),
        "n_B_only_median": float(sub_df["n_B_only"].median()),
        "n_union_median": float(sub_df["n_union"].median()),
        "pct_shared_median": float(sub_df["pct_shared"].median()),
        "pct_shared_q025": float(sub_df["pct_shared"].quantile(0.025)),
        "pct_shared_q975": float(sub_df["pct_shared"].quantile(0.975)),
        "pct_group_specific_median": float(sub_df["pct_group_specific"].median()),
        "pct_group_specific_q025": float(sub_df["pct_group_specific"].quantile(0.025)),
        "pct_group_specific_q975": float(sub_df["pct_group_specific"].quantile(0.975)),
        "pct_B_only_median": float(sub_df["pct_B_only"].median()),
        "pct_B_only_q025": float(sub_df["pct_B_only"].quantile(0.025)),
        "pct_B_only_q975": float(sub_df["pct_B_only"].quantile(0.975)),
    }
    # store medians also under common keys for scenario table
    scenarios.append(
        {
            "scenario": summ["scenario"],
            "note": summ["note"],
            "n_A": 10,
            "n_B": 10,
            "n_KO_A": summ["n_KO_A_median"],
            "n_KO_B": summ["n_KO_B_median"],
            "n_shared": summ["n_shared_median"],
            "n_A_only": summ["n_A_only_median"],
            "n_B_only": summ["n_B_only_median"],
            "n_union": summ["n_union_median"],
            "pct_shared": summ["pct_shared_median"],
            "pct_group_specific": summ["pct_group_specific_median"],
            "pct_A_only": float(sub_df["pct_A_only"].median()),
            "pct_B_only": summ["pct_B_only_median"],
        }
    )

    # equal HQ subsample: min(n_hq, p_hq)
    n_eq = min(len(n_hq), len(p_hq))
    hq_rows = []
    if n_eq >= 2:
        for i in range(B_SUB):
            dn = list(rng.choice(n_hq, size=n_eq, replace=False))
            dp = list(rng.choice(p_hq, size=n_eq, replace=False))
            st = ko_union_stats(ko, dn, dp)
            st["iter"] = i
            hq_rows.append(st)
        hq_df = pd.DataFrame(hq_rows)
        dual_write(hq_df, "Table_R2.8_equalHQ_subsample_iterations")
        scenarios.append(
            {
                "scenario": f"equalHQ_subsample_{n_eq}v{n_eq}",
                "note": f"B={B_SUB}; HQ only; seed={SEED}",
                "n_A": n_eq,
                "n_B": n_eq,
                "n_KO_A": float(hq_df["n_KO_A"].median()),
                "n_KO_B": float(hq_df["n_KO_B"].median()),
                "n_shared": float(hq_df["n_shared"].median()),
                "n_A_only": float(hq_df["n_A_only"].median()),
                "n_B_only": float(hq_df["n_B_only"].median()),
                "n_union": float(hq_df["n_union"].median()),
                "pct_shared": float(hq_df["pct_shared"].median()),
                "pct_group_specific": float(hq_df["pct_group_specific"].median()),
                "pct_A_only": float(hq_df["pct_A_only"].median()),
                "pct_B_only": float(hq_df["pct_B_only"].median()),
            }
        )
        summ_hq = {
            "scenario": f"equalHQ_subsample_{n_eq}v{n_eq}",
            "pct_shared_median": float(hq_df["pct_shared"].median()),
            "pct_shared_q025": float(hq_df["pct_shared"].quantile(0.025)),
            "pct_shared_q975": float(hq_df["pct_shared"].quantile(0.975)),
            "pct_group_specific_median": float(hq_df["pct_group_specific"].median()),
            "pct_group_specific_q025": float(hq_df["pct_group_specific"].quantile(0.025)),
            "pct_group_specific_q975": float(hq_df["pct_group_specific"].quantile(0.975)),
            "pct_B_only_median": float(hq_df["pct_B_only"].median()),
            "pct_B_only_q025": float(hq_df["pct_B_only"].quantile(0.025)),
            "pct_B_only_q975": float(hq_df["pct_B_only"].quantile(0.975)),
        }
    else:
        summ_hq = {}

    scen_df = pd.DataFrame(scenarios)
    dual_write(scen_df, "Table_R2.8_KO_union_sensitivity_summary")
    dual_write(pd.DataFrame([summ]), "Table_R2.8_equalN_subsample_summary")
    if summ_hq:
        dual_write(pd.DataFrame([summ_hq]), "Table_R2.8_equalHQ_subsample_summary")

    # ---- 3. key KO Δprop under equal-n subsample ----
    # map matrix index: KO ids without suffix
    key_rows = []
    present_keys = [k for k in KEY_KOS if k in ko.index]
    # also pull top fisher padj KOs
    fisher = pd.read_csv(FISHER_PATH, sep="\t")
    fisher["KO_id"] = fisher["KO"].astype(str).str.split(":").str[0]
    top_fisher = fisher.nsmallest(30, "padj")["KO_id"].tolist()
    focus = list(dict.fromkeys(present_keys + [k for k in top_fisher if k in ko.index]))

    base_deltas = {k: prop_delta(ko, normal, poor, k) for k in focus}
    for k in focus:
        d0 = base_deltas[k]
        if d0 is None:
            continue
        deltas = []
        for i in range(B_SUB):
            draw = list(rng.choice(poor, size=len(normal), replace=False))
            d = prop_delta(ko, normal, draw, k)
            if d is not None:
                deltas.append(d)
        arr = np.array(deltas)
        same_sign = np.mean(np.sign(arr) == np.sign(d0)) if d0 != 0 else np.nan
        key_rows.append(
            {
                "KO": k,
                "gene": KEY_KOS.get(k, ""),
                "delta_prop_baseline_10v35": d0,
                "delta_prop_equalN_median": float(np.median(arr)),
                "delta_prop_equalN_q025": float(np.quantile(arr, 0.025)),
                "delta_prop_equalN_q975": float(np.quantile(arr, 0.975)),
                "sign_consistency_rate": float(same_sign),
                "in_manuscript_highlight": k in KEY_KOS,
            }
        )
    key_df = pd.DataFrame(key_rows)
    if len(key_df) == 0:
        key_df = pd.DataFrame(
            columns=[
                "KO",
                "gene",
                "delta_prop_baseline_10v35",
                "delta_prop_equalN_median",
                "delta_prop_equalN_q025",
                "delta_prop_equalN_q975",
                "sign_consistency_rate",
                "in_manuscript_highlight",
            ]
        )
    else:
        key_df = key_df.sort_values("sign_consistency_rate", ascending=False)
    dual_write(key_df, "Table_R2.8_keyKO_delta_prop_equalN")

    # HQ equal for key KOs (if possible)
    if n_eq >= 2 and len(n_hq) and len(p_hq):
        hq_key = []
        for k in focus:
            d0 = prop_delta(ko, n_hq, p_hq, k)
            if d0 is None:
                continue
            deltas = []
            for i in range(min(B_SUB, 1000)):
                dn = list(rng.choice(n_hq, size=n_eq, replace=False))
                dp = list(rng.choice(p_hq, size=n_eq, replace=False))
                d = prop_delta(ko, dn, dp, k)
                if d is not None:
                    deltas.append(d)
            arr = np.array(deltas)
            hq_key.append(
                {
                    "KO": k,
                    "gene": KEY_KOS.get(k, ""),
                    "delta_prop_HQ_unequal": d0,
                    "delta_prop_equalHQ_median": float(np.median(arr)),
                    "sign_consistency_vs_HQunequal": float(
                        np.mean(np.sign(arr) == np.sign(d0))
                    )
                    if d0 != 0
                    else np.nan,
                }
            )
        dual_write(pd.DataFrame(hq_key), "Table_R2.8_keyKO_delta_prop_equalHQ")

    # ---- 4. recovery vs detection (cite R1.6) + depth/recovery proxy ----
    abund = pd.read_csv(ABUND_PATH, sep="\t", index_col=0)
    # columns are samples
    samples = list(abund.columns)
    # map A* Normal, B* PoorPlaque from env / naming
    env = pd.read_csv(ENV_PATH, sep="\t", index_col=0)
    # Sample row may exist; use column names: A=Normal, B=PoorPlaque convention in this project
    sample_group = {
        s: ("Normal" if s.startswith("A") else "PoorPlaque" if s.startswith("B") else "NA")
        for s in samples
    }
    # library-size proxy: sum of MAG abundances
    lib = abund.sum(axis=0)
    n_detected = (abund > 0).sum(axis=0)
    # recovery origin counts from QC group GP/PP
    recovery = (
        qc.groupby("group")
        .size()
        .rename({"GP": "Normal_assembly_origin", "PP": "PoorPlaque_assembly_origin"})
    )
    depth_df = pd.DataFrame(
        {
            "Sample": samples,
            "Group": [sample_group[s] for s in samples],
            "MAG_abundance_sum_proxy": lib.values,
            "n_MAG_detected_gt0": n_detected.values,
        }
    )
    dual_write(depth_df, "Table_R2.8_depth_recovery_proxy_by_sample")

    # correlate proxy depth with n detected within cohort
    corr_rows = []
    for g, sub in depth_df.groupby("Group"):
        if len(sub) < 3:
            continue
        r, p = spearmanr(sub["MAG_abundance_sum_proxy"], sub["n_MAG_detected_gt0"])
        corr_rows.append(
            {
                "Group": g,
                "spearman_rho_libsize_vs_nDetected": r,
                "p": p,
                "n": len(sub),
                "median_libsize_proxy": float(sub["MAG_abundance_sum_proxy"].median()),
                "median_n_detected": float(sub["n_MAG_detected_gt0"].median()),
            }
        )
    # between-group recovery origin
    corr_rows.append(
        {
            "Group": "assembly_origin_counts",
            "spearman_rho_libsize_vs_nDetected": np.nan,
            "p": np.nan,
            "n": np.nan,
            "median_libsize_proxy": np.nan,
            "median_n_detected": np.nan,
            "n_MAGs_recovered_from_Normal_samples_GP": int(
                (qc["group"] == "GP").sum()
            ),
            "n_MAGs_recovered_from_PoorPlaque_samples_PP": int(
                (qc["group"] == "PP").sum()
            ),
        }
    )
    # fix columns for first rows
    corr_df = pd.DataFrame(corr_rows)
    dual_write(corr_df, "Table_R2.8_depth_vs_detection_summary")

    # R1.6 detection summary copy-forward
    if DET_PATH.exists():
        det = pd.read_csv(DET_PATH, sep="\t")
        dual_write(det, "Table_R2.8_from_R1.6_detection_summary")

    # ---- 5. cross-cohort fine features: Robust_DA ThisStudy directional consistency ----
    cross = pd.read_csv(CROSS_PATH)
    rob = cross[cross["Robust_TS"] == True].copy()  # noqa: E712
    rob["same_direction"] = rob["Dir_TS"] == rob["Dir_766"]
    # Dir_766 may be NS
    rob["concordant_nonNS"] = (rob["Dir_766"] != "NS") & (
        rob["Dir_TS"] == rob["Dir_766"]
    )
    rob["discordant_nonNS"] = (rob["Dir_766"] != "NS") & (
        rob["Dir_TS"] != rob["Dir_766"]
    )
    summary_cross = pd.DataFrame(
        [
            {
                "n_Robust_DA_ThisStudy_in_table": len(rob),
                "n_caries_also_Robust": int((rob["Robust_766"] == True).sum()),  # noqa: E712
                "n_caries_Dir_NS": int((rob["Dir_766"] == "NS").sum()),
                "n_same_Dir_including_when_766_NS": int(rob["same_direction"].sum()),
                "n_concordant_both_directional": int(rob["concordant_nonNS"].sum()),
                "n_discordant_both_directional": int(rob["discordant_nonNS"].sum()),
                "frac_ThisStudy_Robust_NS_in_caries": float(
                    (rob["Dir_766"] == "NS").mean()
                ),
            }
        ]
    )
    dual_write(summary_cross, "Table_R2.8_cross_cohort_direction_summary")
    dual_write(
        rob[
            [
                "species",
                "log2FC_TS",
                "Dir_TS",
                "Robust_TS",
                "log2FC_766",
                "Dir_766",
                "Robust_766",
                "same_direction",
                "concordant_nonNS",
                "discordant_nonNS",
            ]
        ],
        "Table_R2.8_cross_cohort_RobustTS_direction",
    )

    # ---- 6. decision one-pager ----
    highlight = key_df[key_df["in_manuscript_highlight"]].copy()
    if len(highlight):
        med_cons = float(highlight["sign_consistency_rate"].median())
        min_cons = float(highlight["sign_consistency_rate"].min())
    else:
        med_cons = min_cons = np.nan

    decision = pd.DataFrame(
        [
            {
                "item": "baseline_group_specific_KO_pct",
                "value": round(base["pct_group_specific"], 2),
                "interpretation": "Fig.4C exclusive fraction under unequal 10 vs 35",
            },
            {
                "item": "equalN_group_specific_KO_pct_median_95CI",
                "value": f"{summ['pct_group_specific_median']:.1f} ({summ['pct_group_specific_q025']:.1f}-{summ['pct_group_specific_q975']:.1f})",
                "interpretation": "Exclusive fraction shrinks when n matched → inflated by MAG count",
            },
            {
                "item": "equalN_PoorPlaque_only_pct_median_95CI",
                "value": f"{summ['pct_B_only_median']:.1f} ({summ['pct_B_only_q025']:.1f}-{summ['pct_B_only_q975']:.1f})",
                "interpretation": "PP-only union KOs drop sharply under equal n",
            },
            {
                "item": "HQ_RobustDA_counts",
                "value": f"Normal {len(n_hq)} / PoorPlaque {len(p_hq)}",
                "interpretation": "Use HQ subset + equal HQ subsample as sensitivity",
            },
            {
                "item": "key_highlight_KO_sign_consistency_median_min",
                "value": f"{med_cons:.2f} / {min_cons:.2f}",
                "interpretation": "Fraction of equal-n draws keeping baseline Δprop sign",
            },
            {
                "item": "cross_cohort_RobustTS_NS_in_caries_frac",
                "value": round(
                    float(summary_cross["frac_ThisStudy_Robust_NS_in_caries"].iloc[0]), 3
                ),
                "interpretation": "Most plaque-control Robust species are NS in caries contrast → not validation",
            },
            {
                "item": "shared_MAG_criterion",
                "value": "dRep secondary_cluster ANI>=95% present in all 3 collections (15/208=7.2%)",
                "interpretation": "Catalog overlap, not phenotype replication",
            },
        ]
    )
    dual_write(decision, "Table_R2.8_decision_summary")

    # compact printout
    print("=== R2.8 sensitivity done ===")
    print(scen_df[["scenario", "n_A", "n_B", "pct_shared", "pct_group_specific", "pct_B_only"]].to_string(index=False))
    print("equalN summary:", summ)
    print("highlight KO consistency median/min:", med_cons, min_cons)
    print("cross:", summary_cross.to_dict("records")[0])
    print("wrote tables to", TAB)


if __name__ == "__main__":
    main()
