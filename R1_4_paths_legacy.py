"""Packaged for revise/; source: project/review/R1.4/scripts/00_paths.py
Figures: shared
Note: legacy absolute paths (reference)
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

"""Shared paths and constants for R1.4 validation."""
from pathlib import Path

ROOT = PKG_ROOT
OUT = ROOT / "review" / "R1.4"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
TABLES.mkdir(parents=True, exist_ok=True)
FIGS.mkdir(parents=True, exist_ok=True)

SEED = 20261021
N_BOOT = 500
N_PERM = 999

ENV = ROOT / "01this_study_dataset/input/env.tsv"
GROUP = ROOT / "01this_study_dataset/sgb/species_network/this.study.group.tsv"
MAG_AB = ROOT / "01this_study_dataset/input/bin_abundance_table.tsv"
METAB_AB = ROOT / "01this_study_dataset/sgb/metab_corr/input/metab_abund.AB.txt"
KO_MAT = ROOT / "01this_study_dataset/M2-humann/comb_genefamilies-cpm.uniref90_ko.normed.named.tsv"
DA_KO = ROOT / "01this_study_dataset/M2-humann/v1/Figure2C/Fig2_PanelC.results/robust_DA_KO_list.csv"
DA_METAB = ROOT / "01this_study_dataset/sgb/metab_corr/paper_result/05_RobustDA_metab_DA/tables/metabolite_RobustDA_only_with_HMDB_Class.tsv"
DA_MAG = ROOT / "01this_study_dataset/sgb/metab_corr/paper_result/00_RobustDA_MAG_set/tables/RobustDA_MAGs.tsv"
EDGE_N = ROOT / "01this_study_dataset/sgb/species_network/result_bins/Normal/Normal_Gephi_edge.csv"
EDGE_P = ROOT / "01this_study_dataset/sgb/species_network/result_bins/PoorPlaque/PoorPlaque_Gephi_edge.csv"
PI_READY = ROOT / "01this_study_dataset/sgb/metab_corr/paper_result/13_RobustDA_KO_metab_PI_network/tables/flagellar_KO_FAD_PI_sample_matrix.tsv"
NET_SUMMARY = ROOT / "01this_study_dataset/sgb/species_network/result_bins/summary/network_property_summary.tsv"

# Manuscript network thresholds
R_MIN = 0.5
P_MAX = 0.05
# Cross-omics PI network (manuscript)
RHO_PI = 0.80
