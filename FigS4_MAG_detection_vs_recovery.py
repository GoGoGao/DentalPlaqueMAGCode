"""Packaged for revise/; source: project/review/R1.6/scripts/R1.6_mag_detection_vs_recovery.py
Figures: FigS4
Note: detection vs recovery tables
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
"""R1.6: MAG assembly-recovery vs read-recruitment detection across groups."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = PKG_ROOT
ABUND = ROOT / "01this_study_dataset/input/bin_abundance_table.tsv"
GROUP = ROOT / "01this_study_dataset/sgb/species_network/this.study.group.tsv"
GTDB = ROOT / "01this_study_dataset/input/gtdbtk.bac120.summary.tsv"
OUTDIR = ROOT / "review/R1.6/tables"

# Manuscript-emphasized taxa (GTDB species suffix / genus match)
HIGHLIGHT_SPECIES = [
    "Porphyromonas gingivalis",
    "Tannerella forsythia",
    "Prevotella intermedia",
    "Parvimonas micra",
    "Filifactor alocis",
    "Peptostreptococcus stomatis",
    "Porphyromonas endodontalis",
    "Phocaeicola abscessus",
    "Lachnoanaerobaculum sp938041725",
    "Lachnoanaerobaculum gingivalis",
    "Rothia dentocariosa",
]

THRESHOLDS = [0.0, 0.01, 0.1]


def write_tsv_csv(df: pd.DataFrame, stem: str) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTDIR / f"{stem}.tsv", sep="\t", index=False)
    df.to_csv(OUTDIR / f"{stem}.csv", index=False)
    print(f"Wrote {stem} ({df.shape[0]} rows)")


def species_from_gtdb(classification: str) -> str:
    parts = str(classification).split(";")
    for p in reversed(parts):
        p = p.strip()
        if p.startswith("s__"):
            return p.replace("s__", "") or "unassigned"
        if p.startswith("g__") and (not parts[-1].startswith("s__") or parts[-1] in ("s__", "s__")):
            return p.replace("g__", "") + " (genus)"
    return "unassigned"


def load():
    abund = pd.read_csv(ABUND, sep="\t")
    abund = abund.rename(columns={abund.columns[0]: "MAG"})
    group = pd.read_csv(GROUP, sep="\t", header=None, names=["Sample", "Group"])
    gm = dict(zip(group["Sample"], group["Group"]))
    sample_cols = [c for c in abund.columns if c != "MAG"]
    missing = [c for c in sample_cols if c not in gm]
    if missing:
        raise SystemExit(f"Missing group labels: {missing}")

    gtdb = pd.read_csv(GTDB, sep="\t")
    gtdb["MAG"] = gtdb["user_genome"].astype(str).str.extract(r"(MAG\d+)", expand=False)
    gtdb["gtdb_species"] = gtdb["classification"].map(species_from_gtdb)
    gtdb["gtdb_classification"] = gtdb["classification"]
    tax = gtdb[["MAG", "gtdb_species", "gtdb_classification"]].drop_duplicates("MAG")
    return abund, gm, sample_cols, tax


def summarize(abund, gm, sample_cols, thr: float) -> pd.DataFrame:
    normal = [s for s in sample_cols if gm[s] == "Normal"]
    poor = [s for s in sample_cols if gm[s] == "PoorPlaque"]
    rows = []
    for _, r in abund.iterrows():
        n_vals = r[normal].astype(float)
        p_vals = r[poor].astype(float)
        n_det = int((n_vals > thr).sum())
        p_det = int((p_vals > thr).sum())
        if n_det > 0 and p_det > 0:
            label = "detected_both_groups"
        elif n_det > 0:
            label = "detected_Normal_only"
        elif p_det > 0:
            label = "detected_PoorPlaque_only"
        else:
            label = "undetected"
        rows.append(
            {
                "MAG": r["MAG"],
                "threshold": thr,
                "n_detected_Normal": n_det,
                "n_detected_PoorPlaque": p_det,
                "mean_abund_Normal": float(n_vals.mean()),
                "mean_abund_PoorPlaque": float(p_vals.mean()),
                "max_abund_Normal": float(n_vals.max()),
                "max_abund_PoorPlaque": float(p_vals.max()),
                "detection_label": label,
            }
        )
    return pd.DataFrame(rows)


def main():
    abund, gm, sample_cols, tax = load()
    frames = [summarize(abund, gm, sample_cols, thr) for thr in THRESHOLDS]
    full = pd.concat(frames, ignore_index=True).merge(tax, on="MAG", how="left")
    write_tsv_csv(full, "Table_R1.6_MAG_detection_by_group")

    primary = full[full["threshold"] == 0.0].copy()
    # highlight manuscript taxa
    mask = primary["gtdb_species"].fillna("").isin(HIGHLIGHT_SPECIES)
    # also genus-level loose match for any missed
    for name in HIGHLIGHT_SPECIES:
        mask = mask | primary["gtdb_species"].fillna("").str.contains(name.split()[-1], case=False, na=False)
    # keep unique by preferred list first
    hi = primary[primary["gtdb_species"].fillna("").isin(HIGHLIGHT_SPECIES)].copy()
    if hi.empty:
        hi = primary.loc[mask].copy()
    hi.insert(0, "manuscript_role", "named_in_Results_group_association_sentence")
    write_tsv_csv(hi.sort_values(["gtdb_species", "MAG"]), "Table_R1.6_highlight_pathogens_and_named")

    # exclusive at each threshold
    excl = full[full["detection_label"] != "detected_both_groups"].copy()
    write_tsv_csv(excl, "Table_R1.6_abundance_one_group_only")

    summary_rows = []
    for thr, sub in full.groupby("threshold"):
        vc = sub["detection_label"].value_counts()
        both = int(vc.get("detected_both_groups", 0))
        summary_rows.append(
            {
                "threshold": thr,
                "n_MAG": len(sub),
                "detected_both_groups": both,
                "detected_Normal_only": int(vc.get("detected_Normal_only", 0)),
                "detected_PoorPlaque_only": int(vc.get("detected_PoorPlaque_only", 0)),
                "undetected": int(vc.get("undetected", 0)),
                "pct_both": round(100.0 * both / len(sub), 1),
            }
        )
    write_tsv_csv(pd.DataFrame(summary_rows), "Table_R1.6_summary_counts")

    s0 = summary_rows[0]
    print("\n=== Takeaway (threshold=0) ===")
    print(
        f"Of {s0['n_MAG']} MAGs: both={s0['detected_both_groups']} ({s0['pct_both']}%), "
        f"Normal-only={s0['detected_Normal_only']}, PoorPlaque-only={s0['detected_PoorPlaque_only']}"
    )
    if not hi.empty:
        print("\nHighlight taxa (thr=0):")
        for _, r in hi.iterrows():
            print(
                f"  {r['MAG']} {r['gtdb_species']}: "
                f"N={r['n_detected_Normal']}/10 P={r['n_detected_PoorPlaque']}/10 "
                f"meanN={r['mean_abund_Normal']:.3f} meanP={r['mean_abund_PoorPlaque']:.3f} "
                f"-> {r['detection_label']}"
            )


if __name__ == "__main__":
    main()
