## Packaged for revise/; source: project/01this_study_dataset/M2-humann/v2/Figure3/code/Fig3_enrich_bubble.R
## Figures: Fig1F
## Note: KO KEGG enrichment bubble
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig1F_KO_pathway_enrichment.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript

## ============================================================================
## Fig3 KO Pathway enrichment bubble plot（美化版气泡图）
## ----------------------------------------------------------------------------
## 输入：
## - Figure3/code/Fig3_enrich_results/Fig3_pathway_enrichment_Normal_up.tsv
## - Figure3/code/Fig3_enrich_results/Fig3_pathway_enrichment_PoorPlaque_up.tsv
## 输出（到 out_dir）：
## - Fig3_enrich_bubble.pdf
## - Fig3_enrich_bubble.svg
##
## 运行示例（项目根目录）：
##   Rscript Figure3/code/Fig3_enrich_bubble.R \
##     enrich_dir=Figure3/code/Fig3_enrich_results \
##     out_dir=Figure3/code/Fig3_enrich_results
## ============================================================================

suppressPackageStartupMessages({
  library(glue)
  library(tidyverse)
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

enrich_dir <- get_arg("enrich_dir", "Figure3/code/Fig3_enrich_results")
out_dir    <- get_arg("out_dir",    enrich_dir)
top_n      <- as.integer(get_arg("top_n", "15"))

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

message(glue("[INFO] Fig3_enrich_bubble | enrich_dir = {enrich_dir}"))
message(glue("[INFO] Fig3_enrich_bubble | out_dir    = {out_dir}"))

## 读取富集结果 ---------------------------------------------------------------

norm_path <- file.path(enrich_dir, "Fig3_pathway_enrichment_Normal_up.tsv")
poor_path <- file.path(enrich_dir, "Fig3_pathway_enrichment_PoorPlaque_up.tsv")

normal <- readr::read_tsv(norm_path, show_col_types = FALSE)
poor   <- readr::read_tsv(poor_path, show_col_types = FALSE)

required_cols <- c("category", "ID", "Description", "p.adjust", "qvalue", "GeneRatio", "Count")

# [BUGFIX] 循环元素本身就是 tibble，不能再用 [[1]] 取第一列
for (df in list(normal, poor)) {
  miss <- setdiff(required_cols, names(df))
  if (length(miss) > 0L) {
    stop(glue("缺少必要列：{paste(miss, collapse = ', ')}"))
  }
}

## 统一整理数据 ---------------------------------------------------------------

normal$Group <- "Normal"
poor$Group   <- "PoorPlaque"

all_df <- dplyr::bind_rows(normal, poor) |>
  dplyr::mutate(
    Group         = factor(Group, levels = c("Normal", "PoorPlaque")),
    neg_log10_FDR = -log10(p.adjust),
    # 提取 GeneRatio 的数值部分： \"x/y\" -> x/y
    GeneRatio_num = purrr::map_dbl(GeneRatio, function(x) {
      if (is.na(x)) return(NA_real_)
      sp <- strsplit(as.character(x), "/", fixed = TRUE)[[1]]
      if (length(sp) == 2L) {
        as.numeric(sp[1]) / as.numeric(sp[2])
      } else {
        NA_real_
      }
    })
  )

## 选择展示的 pathway（每组 Top N） -----------------------------------------

top_df <- all_df |>
  dplyr::group_by(Group) |>
  dplyr::arrange(dplyr::desc(neg_log10_FDR)) |>
  dplyr::slice_head(n = top_n) |>
  dplyr::ungroup()

# 为了拼图美观，将Description按neg_log10_FDR排序，统一为因子
top_df <- top_df |>
  dplyr::arrange(Group, neg_log10_FDR) |>
  dplyr::mutate(
    Description = stringr::str_trunc(Description, width = 50),
    Description = factor(Description, levels = unique(Description))
  )

## 分组配色与主题（与 KO 火山图一致的组色） ----------------------------------

COL_GROUP <- c(
  Normal     = "#2E8B8B",
  PoorPlaque = "#B85450"
)

## 绘制气泡图 -----------------------------------------------------------------

p <- ggplot(top_df, aes(x = Group, y = Description)) +
  geom_point(
    aes(size = Count, fill = neg_log10_FDR),
    shape  = 21,
    color  = "black",
    stroke = 0.2,
    alpha  = 0.9
  ) +
  scale_size_continuous(
    range = c(2.0, 7.0),
    name  = "KO count"
  ) +
  scale_fill_gradientn(
    colours = rev(RColorBrewer::brewer.pal(9, "YlOrRd")),
    name    = expression(-log[10]("FDR"))
  ) +
  facet_grid(
    category ~ .,
    scales  = "free_y",
    space   = "free_y"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major   = element_line(color = "#F0F0F0", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "white", color = "grey60", linewidth = 0.4),
    strip.text.y       = element_text(size = 9, face = "bold"),
    axis.text.x        = element_text(size = 9, face = "bold"),
    axis.text.y        = element_text(size = 7),
    legend.position    = "right",
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 8),
    plot.title         = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin        = margin(5.5, 5.5, 5.5, 5.5)
  ) +
  ggtitle("KEGG pathway enrichment for differential KOs")

## 导出 PDF + SVG -------------------------------------------------------------

pdf_file <- file.path(out_dir, "Fig3_enrich_bubble.pdf")
svg_file <- file.path(out_dir, "Fig3_enrich_bubble.svg")

Cairo::CairoPDF(
  file   = pdf_file,
  width  = 5.0,
  height = 5.0,
  family = "sans"
)
print(p)
dev.off()

Cairo::CairoSVG(
  file   = svg_file,
  width  = 5.0,
  height = 5.0,
  family = "sans"
)
print(p)
dev.off()

message("[INFO] Enrichment bubble plot saved: ")
message(glue("  - {pdf_file}"))
message(glue("  - {svg_file}"))

