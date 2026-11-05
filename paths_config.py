from pathlib import Path
import os
PKG_ROOT = Path(__file__).resolve().parents[1]
SHARED = PKG_ROOT / "shared_input"
RSCRIPT = os.environ.get("RSCRIPT", "Rscript")

PATHS = {
  "sample_group": SHARED / "sample_group.tsv",
  "clinical_env": SHARED / "clinical_env.tsv",
  "metaphlan_species_abundance": SHARED / "metaphlan_species_abundance.xls",
  "humann_ko_cpm": SHARED / "humann_ko_cpm.tsv",
  "humann_ko_cpm_named_all": SHARED / "humann_ko_cpm_named_all.tsv",
  "robust_DA_KO_list": SHARED / "robust_DA_KO_list.csv",
  "metab_abundance_AB": SHARED / "metab_abundance_AB.tsv",
  "metab_abundance_AB_select": SHARED / "metab_abundance_AB_select.tsv",
  "metab_desc": SHARED / "metab_desc.tsv",
  "metab_select": SHARED / "metab_select.tsv",
  "robust_DA_metabolites": SHARED / "robust_DA_metabolites.tsv",
  "mag_bin_abundance": SHARED / "mag_bin_abundance.tsv",
  "mag_gtdb_summary": SHARED / "mag_gtdb_summary.tsv",
  "mag_plot_input": SHARED / "mag_plot_input.tsv",
  "mag_ko_presence": SHARED / "mag_ko_presence.xls",
  "mag_vfdb": SHARED / "mag_vfdb.xls",
  "mag_card": SHARED / "mag_card.xls",
  "mag_phylogeny_unrooted": SHARED / "mag_phylogeny_unrooted.nwk",
  "robust_DA_MAGs": SHARED / "robust_DA_MAGs.tsv",
  "crosscohort_group": SHARED / "crosscohort_group.tsv",
  "crosscohort_gtdb_summary": SHARED / "crosscohort_gtdb_summary.tsv",
  "crosscohort_total_KO": SHARED / "crosscohort_total_KO.xls",
}
