## Packaged for revise/; source: project/01this_study_dataset/M2-humann/v2/Figure3/code/Fig3_M2_pipeline.R
## Figures: Fig1D;Fig1E;Fig1F;Fig1G
## Note: KO analysis pipeline
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig1_KO_pipeline.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript

## ============================================================================
## Fig3 M2 分析总入口脚本（放在 Figure3/code 目录下）
## ----------------------------------------------------------------------------
## 用途：
## - 统一从命令行调用 Figure3 M2 功能层全部分析：
##   - Fig3A：Top 差异 KO 物种贡献热图
##   - Fig3A-all：全部 Robust_DA KO 物种贡献热图
##   - Fig3B：Marker 物种 × DA KO 相关性热图
##   - KO Pathway 富集分析（两组分别上调 KO）
##   - 桑基图：Top 10 差异 KO 的主要贡献物种
## 运行方式（在项目根目录）：
##   Rscript Figure3/code/Fig3_M2_pipeline.R
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

here_root <- getwd()

message(glue("[INFO] Project root: {here_root}"))

source(file.path("Figure3", "Figure3_Fig3_M2_analysis.R"))

tables_dir <- file.path("Figure3", "tables")
panels_dir <- file.path("Figure3", "code", "Fig3_M2_results")

if (!dir.exists(panels_dir)) {
  dir.create(panels_dir, recursive = TRUE)
}

message(glue("[INFO] Tables dir : {tables_dir}"))
message(glue("[INFO] Panels dir : {panels_dir}"))

res <- run_figure3_analysis(
  tables_dir  = tables_dir,
  panels_dir  = panels_dir,
  top_n_ko    = 60L,
  min_contrib = 0.05,
  fdr_cutoff  = 0.05,
  rho_cutoff  = 0.5
)

message("[INFO] Fig3 M2 analysis completed.")

