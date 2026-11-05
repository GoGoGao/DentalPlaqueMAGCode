## Packaged for revise/; source: project/01this_study_dataset/M2-humann/v2/Figure3/code/Fig3_KO_pathway_enrichment.R
## Figures: Fig1F
## Note: KO enrichment stats
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig1F_KO_pathway_enrichment_stats.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript

## ============================================================================
## Fig3 KO Pathway 富集分析 独立脚本
## ----------------------------------------------------------------------------
## 输入：
## - Figure3/tables/Table_differential_KO_full.tsv（必需）
## - Figure3/tables/Table_spearman_marker_KO.tsv（可选，仅当 use_marker_ko=TRUE 时）
## 输出：
## - Fig3_pathway_enrichment_Normal_up.tsv / .pdf（所有差异KO）
## - Fig3_pathway_enrichment_PoorPlaque_up.tsv / .pdf（所有差异KO）
## - Fig3_pathway_enrichment_Normal_marker_up.tsv / .pdf（marker KO，可选）
## - Fig3_pathway_enrichment_PoorPlaque_marker_up.tsv / .pdf（marker KO，可选）
## 运行示例：
##   # 所有差异KO的富集
##   Rscript Figure3/code/Fig3_KO_pathway_enrichment.R \\
##     out_dir=Figure3/code/Fig3_enrich_results
##   # 仅marker KO的富集
##   Rscript Figure3/code/Fig3_KO_pathway_enrichment.R \\
##     out_dir=Figure3/code/Fig3_enrich_results \\
##     use_marker_ko=TRUE marker_fdr_cutoff=0.05 marker_rho_cutoff=0.5
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
    val <- arg_list[[name]]
    # 布尔值转换
    if (val %in% c("TRUE", "True", "true", "T", "1")) {
      return(TRUE)
    }
    if (val %in% c("FALSE", "False", "false", "F", "0")) {
      return(FALSE)
    }
    # 数值转换
    if (grepl("^[0-9.]+$", val)) {
      if (grepl("\\.", val)) {
        return(as.numeric(val))
      } else {
        return(as.integer(val))
      }
    }
    return(val)
  }
  default
}

tables_dir        <- get_arg("tables_dir",        "Figure3/tables")
out_dir           <- get_arg("out_dir",           "Figure3/code/Fig3_enrich_results")
use_marker_ko     <- get_arg("use_marker_ko",     FALSE)
marker_fdr_cutoff <- get_arg("marker_fdr_cutoff", 0.05)
marker_rho_cutoff <- get_arg("marker_rho_cutoff", 0.5)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

message(glue("[INFO] Fig3_enrich | tables_dir        = {tables_dir}"))
message(glue("[INFO] Fig3_enrich | out_dir           = {out_dir}"))
message(glue("[INFO] Fig3_enrich | use_marker_ko     = {use_marker_ko}"))
if (use_marker_ko) {
  message(glue("[INFO] Fig3_enrich | marker_fdr_cutoff = {marker_fdr_cutoff}"))
  message(glue("[INFO] Fig3_enrich | marker_rho_cutoff = {marker_rho_cutoff}"))
}

source(file.path("Figure3", "Figure3_Fig3_M2_analysis.R"))

diff_ko_path   <- file.path(tables_dir, "Table_differential_KO_full.tsv")
spearman_path  <- if (use_marker_ko) {
  file.path(tables_dir, "Table_spearman_marker_KO.tsv")
} else {
  NULL
}

res <- run_ko_pathway_enrichment(
  diff_ko_path      = diff_ko_path,
  spearman_path     = spearman_path,
  panels_dir        = out_dir,
  use_marker_ko     = use_marker_ko,
  marker_fdr_cutoff = marker_fdr_cutoff,
  marker_rho_cutoff = marker_rho_cutoff
)

message("[INFO] Fig3 KO pathway enrichment finished.")

