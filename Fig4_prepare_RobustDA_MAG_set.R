## Packaged for revise/; source: project/01this_study_dataset/sgb/metab_corr/paper_result/00_RobustDA_MAG_set/scripts/00_prepare_RobustDA_MAG_set.R
## Figures: Fig4B;Fig4C;Fig4D
## Note: prepare Robust_DA MAG set
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig4_prepare_RobustDA_MAG_set.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# Prepare Robust_DA MAG set for KO comparative genomics (paper_result)
# Source: sig_mag/.../Table_SGB_differential_bins.tsv (Robust_DA == TRUE)
# Taxonomy: metab_corr/gtdbtk.bac120.summary.tsv

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
})

base_dir   <- file.path(.pkg_root, "01this_study_dataset/sgb/metab_corr"
sgb_root   <- dirname(base_dir)
module_dir <- file.path(base_dir, "paper_result/00_RobustDA_MAG_set")
out_tab    <- file.path(module_dir, "tables")
out_data   <- file.path(module_dir, "data")
dir.create(out_tab,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)

da_file <- file.path(sgb_root, "sig_mag/sgb_DA_results/tables/Table_SGB_differential_bins.tsv")
gtdb_file <- file.path(base_dir, "gtdbtk.bac120.summary.tsv")
if (!file.exists(gtdb_file)) {
  gtdb_file <- file.path(sgb_root, "compare_genomics_all/gtdbtk.bac120.summary.tsv")
}

da <- read_tsv(da_file, show_col_types = FALSE) %>%
  rename(MAG = SGB) %>%
  mutate(Robust_DA = as.logical(Robust_DA)) %>%
  filter(Robust_DA, Robust_direction %in% c("Normal", "PoorPlaque"))

parse_gtdb_phylum <- function(classification) {
  if (is.na(classification) || !nzchar(classification)) return(NA_character_)
  m <- regexpr("p__[^;]+", classification, perl = TRUE)
  if (m < 1L) return(NA_character_)
  sub("^p__", "", regmatches(classification, m))
}
parse_gtdb_genus <- function(classification) {
  if (is.na(classification) || !nzchar(classification)) return(NA_character_)
  m <- regexpr("g__[^;]+", classification, perl = TRUE)
  if (m < 1L) return(NA_character_)
  sub("^g__", "", regmatches(classification, m))
}

gtdb <- read_tsv(gtdb_file, show_col_types = FALSE) %>%
  transmute(
    MAG = user_genome,
    classification,
    phylum = vapply(classification, parse_gtdb_phylum, character(1)),
    genus  = vapply(classification, parse_gtdb_genus, character(1))
  ) %>%
  mutate(
    phylum = ifelse(is.na(phylum) | !nzchar(phylum), "Unclassified", phylum),
    genus  = ifelse(is.na(genus)  | !nzchar(genus),  "Unclassified", genus)
  )

screen_df <- da %>%
  left_join(gtdb, by = "MAG") %>%
  mutate(
    enrichment = Robust_direction,
    genome = MAG,
    species = coalesce(species, "Unknown"),
    species_plain = str_replace_all(species, "_", " ")
  ) %>%
  arrange(enrichment, desc(abs(log2FC)))

write_tsv(screen_df, file.path(out_tab, "RobustDA_MAGs.tsv"))
writeLines(screen_df$MAG[screen_df$enrichment == "Normal"],
           file.path(out_data, "Normal_MAGs.txt"))
writeLines(screen_df$MAG[screen_df$enrichment == "PoorPlaque"],
           file.path(out_data, "PoorPlaque_MAGs.txt"))

readme <- c(
  "Robust_DA MAG set for paper_result KO comparative genomics",
  "=========================================================",
  paste0("Source: ", da_file),
  "Criterion: Robust_DA == TRUE (Wilcoxon_FDR < 0.05 ∩ CLRMC_FDR < 0.05)",
  "Direction: Robust_direction (Normal / PoorPlaque)",
  paste0("Total MAGs: ", nrow(screen_df)),
  paste0("  Normal-enriched: ", sum(screen_df$enrichment == "Normal")),
  paste0("  PoorPlaque-enriched: ", sum(screen_df$enrichment == "PoorPlaque")),
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)
writeLines(readme, file.path(module_dir, "README_CRITERIA.txt"))
writeLines(readme, file.path(out_tab, "README_CRITERIA.txt"))

message("[DONE] Robust_DA MAG set: ", nrow(screen_df),
        " (Normal=", sum(screen_df$enrichment == "Normal"),
        ", PoorPlaque=", sum(screen_df$enrichment == "PoorPlaque"), ")")
