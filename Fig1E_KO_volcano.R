## Packaged for revise/; source: project/01this_study_dataset/M2-humann/v2/Figure3/code/Fig3_KO_volcano.R
## Figures: Fig1E
## Note: KO volcano
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig1E_KO_volcano.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript

## ============================================================================
## Fig3 KO volcano plot（差异 KO 火山图）
## ----------------------------------------------------------------------------
## 输入：
## - Figure3/tables/Table_differential_KO_full.tsv
## 输出（到 out_dir）：
## - Fig3_KO_volcano.pdf
## - Fig3_KO_volcano.svg
##
## 运行示例（项目根目录）：
##   Rscript Figure3/code/Fig3_KO_volcano.R \
##     tables_dir=Figure3/tables \
##     out_dir=Figure3/code/Fig3_volcano_results
## ============================================================================

suppressPackageStartupMessages({
  library(glue)
  library(tidyverse)
  library(ggrepel)
  library(Cairo)
})

set.seed(42)

options(
  stringsAsFactors = FALSE,
  scipen           = 999,
  warn             = 1
)

## 解析命令行参数 --------------------------------------------------------------
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

tables_dir <- get_arg("tables_dir", "Figure3/tables")
out_dir    <- get_arg("out_dir",    "Figure3/code/Fig3_volcano_results")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

message(glue("[INFO] Fig3_KO_volcano | tables_dir = {tables_dir}"))
message(glue("[INFO] Fig3_KO_volcano | out_dir    = {out_dir}"))

## 读取差异 KO 结果 ------------------------------------------------------------
diff_ko_path <- file.path(tables_dir, "Table_differential_KO_full.tsv")
diff_ko <- readr::read_tsv(diff_ko_path, show_col_types = FALSE)

required_cols <- c(
  "feature", "ko_id", "function_name",
  "log2FC", "Wilcoxon_FDR", "CLRMC_FDR",
  "Robust_DA", "Robust_direction"
)
missing_cols <- setdiff(required_cols, colnames(diff_ko))
if (length(missing_cols) > 0L) {
  stop(glue("缺少必要列：{paste(missing_cols, collapse = ', ')}"))
}

## 分组配色（与 v2.py / M1 保持一致） -----------------------------------------
COL <- c(
  Normal     = "#2E8B8B",
  PoorPlaque = "#B85450",
  NS         = "grey70"
)

## 构建火山图数据 --------------------------------------------------------------

volcano_df <- diff_ko |>
  dplyr::mutate(
    neg_log10_FDR = -log10(Wilcoxon_FDR),
    # 标记显著性：Robust_DA == TRUE 视为显著差异KO
    sig = dplyr::case_when(
      Robust_DA & Robust_direction == "Normal"     ~ "Normal",
      Robust_DA & Robust_direction == "PoorPlaque" ~ "PoorPlaque",
      TRUE                                         ~ "NS"
    ),
    label = function_name
  )

## 选择需要标注的 KO（可根据需要调节阈值） -----------------------------------
top_label <- volcano_df |>
  dplyr::filter(sig != "NS") |>
  dplyr::arrange(dplyr::desc(neg_log10_FDR * abs(log2FC))) |>
  dplyr::slice_head(n = 20)

## 绘图 ------------------------------------------------------------------------

p <- ggplot(volcano_df, aes(x = log2FC, y = neg_log10_FDR)) +
  geom_point(
    aes(color = sig),
    alpha   = 0.8,
    size    = 1.3,
    stroke  = 0
  ) +
  scale_color_manual(
    values = COL,
    breaks = c("Normal", "PoorPlaque", "NS"),
    labels = c("Normal-enriched KO", "PoorPlaque-enriched KO", "Not significant")
  ) +
  geom_vline(xintercept = 0, color = "#CCCCCC", linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), color = "#CCCCCC", linetype = "dotted", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data        = top_label,
    aes(label   = label, color = sig),
    size        = 2.6,
    max.overlaps = 100,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.size  = 0.2,
    show.legend   = FALSE
  ) +
  labs(
    x     = expression(log[2]("PoorPlaque / Normal")),
    y     = expression(-log[10]("Wilcoxon FDR")),
    color = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "#EEEEEE", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 8),
    legend.position  = "top",
    legend.text      = element_text(size = 8),
    plot.title       = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin      = margin(5.5, 5.5, 5.5, 5.5)
  ) +
  ggtitle("Differential KO volcano plot")

## 导出 PDF + SVG（CNS级别，矢量可编辑） ------------------------------------

pdf_file <- file.path(out_dir, "Fig3_KO_volcano.pdf")
svg_file <- file.path(out_dir, "Fig3_KO_volcano.svg")

Cairo::CairoPDF(
  file   = pdf_file,
  width  = 4.0,
  height = 3.5,
  family = "sans"
)
print(p)
dev.off()

Cairo::CairoSVG(
  file   = svg_file,
  width  = 4.0,
  height = 3.5,
  family = "sans"
)
print(p)
dev.off()

message("[INFO] KO volcano plot saved: ")
message(glue("  - {pdf_file}"))
message(glue("  - {svg_file}"))

