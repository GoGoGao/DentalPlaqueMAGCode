## Packaged for revise/; source: project/01this_study_dataset/M2-humann/v2/Figure3/code/Fig3_sankey_top10KO_species.R
## Figures: Fig1G
## Note: top KO-species sankey
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig1G_motility_topKO_species.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript

## ============================================================================
## Fig3 桑基图：Top 10 差异 KO 的主要贡献物种 独立脚本
## ----------------------------------------------------------------------------
## 输入：
## - Figure3/tables/Table_differential_KO_full.tsv
## - Figure3/tables/Table_species_contribution_rank.tsv
## - ../../M1-metaphalan/v2/Figure1/tables/Table_S_differential_species.csv
## 输出：
## - Fig3_sankey_top10_KO_species.pdf
## 运行示例：
##   Rscript Figure3/code/Fig3_sankey_top10KO_species.R \\
##     out_dir=Figure3/code/Fig3_sankey_results \\
##     top_n_ko=10 top_n_species=5
## ============================================================================

suppressPackageStartupMessages({
  library(glue)
})

set.seed(42)

options(
  stringsAsFactors = FALSE,
  scipen           = 999,
  warn             = 1
)

args <- commandArgs(trailingOnly = TRUE)
arg_list <- list()
for (a in args) {
  kv <- strsplit(a, "=", fixed = TRUE)[[1]]
  if (length(kv) == 2L) {
    arg_list[[kv[1]]] <- kv[2]
  }
}

get_arg <- function(name, default = NULL) {
  if (!is.null(arg_list[[name]])) {
    return(arg_list[[name]])
  }
  default
}

tables_dir    <- get_arg("tables_dir",    "Figure3/tables")
out_dir       <- get_arg("out_dir",       "Figure3/code/Fig3_sankey_results")
top_n_ko      <- as.integer(get_arg("top_n_ko",      "10"))
top_n_species <- as.integer(get_arg("top_n_species", "5"))

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

message(glue("[INFO] Fig3_sankey | tables_dir = {tables_dir}"))
message(glue("[INFO] Fig3_sankey | out_dir    = {out_dir}"))

source(file.path("Figure3", "Figure3_Fig3_M2_analysis.R"))

diff_ko_path        <- file.path(tables_dir, "Table_differential_KO_full.tsv")
species_contrib_path <- file.path(tables_dir, "Table_species_contribution_rank.tsv")
diff_species_path    <- file.path("..", "..", "M1-metaphalan", "v2", "Figure1", "tables", "Table_S_differential_species.csv")

res <- plot_sankey_top_ko_species(
  diff_ko_path         = diff_ko_path,
  species_contrib_path = species_contrib_path,
  diff_species_path    = diff_species_path,
  panels_dir           = out_dir,
  top_n_ko             = top_n_ko,
  top_n_species        = top_n_species
)

message("[INFO] Fig3 sankey Top10 KO -> species finished.")

