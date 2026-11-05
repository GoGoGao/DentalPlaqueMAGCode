## Packaged for revise/; source: project/01this_study_dataset/M2-humann/v2/Figure3/code/Fig3_cell_motility_sankey.R
## Figures: Fig1G
## Note: motility KO species contribution
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig1G_motility_species_contribution.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript

## ============================================================================
## Fig3 | Cell motility pathways (ko02040, ko02030) KO→species 桑基图
## ----------------------------------------------------------------------------
## 需求：
## - 基于 PoorPlaque 组富集结果中 Cell motility 两个 pathway：
##   - ko02040 Flagellar assembly
##   - ko02030 Bacterial chemotaxis
## - 找出这些 pathway 中涉及的 KO（geneID 列中的 K 号）
## - 从 submatrix_PoorPlaque.tsv 中计算 KO 在 PoorPlaque 组内的物种贡献度
##   - 每个 KO：按贡献度选 Top 10 物种，其余合并为“Others”
##   - 去除 unclassified 物种
## - 使用桑基图（KO → species）可视化
##   - 差异物种（来自 M1 的 Table_S_differential_species.csv）用红色文字标记
## - 输出：高清 PDF + JPG（Cairo 矢量 / 位图，后续可拼图发表）
##
## 运行示例（项目根目录）：
##   Rscript Figure3/code/Fig3_cell_motility_sankey.R \
##     enrich_path=Figure3/code/Fig3_enrich_results/Fig3_pathway_enrichment_PoorPlaque_up.tsv \
##     submatrix_path=Figure3/tables/submatrix_PoorPlaque.tsv \
##     diff_species_path=../../M1-metaphalan/v2/Figure1/tables/Table_S_differential_species.csv \
##     out_dir=Figure3/code/Fig3_sankey_results
## ============================================================================

suppressPackageStartupMessages({
  library(glue)
  library(tidyverse)
  library(ggalluvial)
  library(Cairo)
})

set.seed(42)

options(
  stringsAsFactors = FALSE,
  scipen           = 999,
  warn             = 1
)

## ---------- 参数解析 ----------------------------------------------------------

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

enrich_path       <- get_arg("enrich_path", "Figure3/code/Fig3_enrich_results/Fig3_pathway_enrichment_PoorPlaque_up.tsv")
submatrix_path    <- get_arg("submatrix_path", "Figure3/tables/submatrix_PoorPlaque.tsv")
diff_species_path <- get_arg("diff_species_path", "../../M1-metaphalan/v2/Figure1/tables/Table_S_differential_species.csv")
out_dir           <- get_arg("out_dir", "Figure3/code/Fig3_sankey_results")
top_n_species     <- as.integer(get_arg("top_n_species", "10"))

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

message(glue("[INFO] Cell motility sankey | enrich_path       = {enrich_path}"))
message(glue("[INFO] Cell motility sankey | submatrix_path    = {submatrix_path}"))
message(glue("[INFO] Cell motility sankey | diff_species_path = {diff_species_path}"))
message(glue("[INFO] Cell motility sankey | out_dir           = {out_dir}"))

## ---------- 工具函数 ----------------------------------------------------------

extract_species_short <- function(species_vec) {
  if (is.null(species_vec)) {
    return(character())
  }
  sub(".*\\.s__", "", species_vec)
}

is_unclassified_species <- function(x) {
  is.na(x) |
    grepl("^unclassified$", x, ignore.case = TRUE) |
    grepl("^unclassified[ _]", x, ignore.case = TRUE)
}

## ---------- 读取富集结果，筛选 Cell motility 两个 pathway --------------------

enrich_df <- readr::read_tsv(enrich_path, show_col_types = FALSE)

required_enrich_cols <- c("category", "subcategory", "ID", "Description", "geneID")
miss_enrich <- setdiff(required_enrich_cols, colnames(enrich_df))
if (length(miss_enrich) > 0L) {
  stop(glue("富集结果缺少必要列：{paste(miss_enrich, collapse = ', ')}"))
}

cell_motility <- enrich_df |>
  dplyr::filter(
    .data$category == "Cellular Processes",
    .data$subcategory == "Cell motility",
    .data$ID %in% c("ko02040", "ko02030")
  )

if (nrow(cell_motility) == 0L) {
  stop("[ERROR] 在富集结果中未找到 Cell motility (ko02040/ko02030) pathway")
}

ko_pathway_df <- cell_motility |>
  dplyr::mutate(
    KO_list = strsplit(.data$geneID, "/", fixed = TRUE)
  ) |>
  tidyr::unnest(cols = "KO_list") |>
  dplyr::transmute(
    ko_id        = KO_list,
    pathway_id   = .data$ID,
    pathway_desc = .data$Description
  )

unique_kos <- unique(ko_pathway_df$ko_id)
message(glue("[INFO] Cell motility 相关 KO 数量：{length(unique_kos)}"))

## ---------- 读取 PoorPlaque submatrix，计算 KO×物种贡献 -----------------------

subm <- readr::read_tsv(submatrix_path, show_col_types = FALSE)
if (!"Features" %in% colnames(subm)) {
  stop("submatrix 文件必须包含 'Features' 列")
}

sample_cols <- setdiff(colnames(subm), "Features")

compute_ko_species_contrib <- function(ko_id) {
  prefix <- paste0(ko_id, ":")
  ko_df <- subm[grepl(paste0("^", prefix), subm$Features), , drop = FALSE]
  if (nrow(ko_df) == 0L) {
    return(NULL)
  }

  # 物种行：含有 '|'
  species_df <- ko_df[grepl("\\|", ko_df$Features), , drop = FALSE]
  if (nrow(species_df) == 0L) {
    return(NULL)
  }

  species_df <- species_df |>
    dplyr::mutate(
      species_full = sub(".*\\|", "", .data$Features),
      species_short = extract_species_short(.data$species_full)
    ) |>
    dplyr::filter(!is_unclassified_species(.data$species_full))

  if (nrow(species_df) == 0L) {
    return(NULL)
  }

  species_sum <- species_df |>
    dplyr::mutate(
      abundance = rowSums(dplyr::across(dplyr::all_of(sample_cols)))
    ) |>
    dplyr::select(.data$species_full, .data$species_short, .data$abundance)

  total_abund <- sum(species_sum$abundance, na.rm = TRUE)
  if (total_abund <= 0) {
    return(NULL)
  }

  species_sum |>
    dplyr::mutate(
      ko_id      = ko_id,
      contrib    = .data$abundance / total_abund
    ) |>
    dplyr::select(.data$ko_id, .data$species_full, .data$species_short, .data$contrib)
}

contrib_list <- lapply(unique_kos, compute_ko_species_contrib)
contrib_list <- contrib_list[!vapply(contrib_list, is.null, logical(1))]

if (length(contrib_list) == 0L) {
  stop("[ERROR] 未能从 submatrix_PoorPlaque 中计算出任何 KO×物种贡献")
}

contrib_df <- dplyr::bind_rows(contrib_list)

## ---------- 结合 pathway 信息，筛选 Top N 物种并合并 Others -------------------

contrib_df <- contrib_df |>
  dplyr::left_join(ko_pathway_df, by = "ko_id")

top_df <- contrib_df |>
  dplyr::group_by(.data$ko_id) |>
  dplyr::arrange(dplyr::desc(.data$contrib), .by_group = TRUE) |>
  dplyr::mutate(rank = dplyr::row_number()) |>
  dplyr::mutate(
    Species_label = dplyr::if_else(.data$rank <= top_n_species, .data$species_short, "Others")
  ) |>
  dplyr::group_by(.data$ko_id, .data$pathway_id, .data$pathway_desc, .data$Species_label) |>
  dplyr::summarise(
    contrib = sum(.data$contrib, na.rm = TRUE),
    .groups = "drop"
  )

## ---------- 标记差异物种（来自 M1） -----------------------------------------

diff_species <- readr::read_csv(diff_species_path, show_col_types = FALSE)
if (!"species" %in% colnames(diff_species)) {
  stop("差异物种表必须包含 'species' 列")
}

marker_species <- diff_species |>
  dplyr::filter(.data$Robust_DA) |>
  dplyr::pull(.data$species) |>
  unique()

top_df <- top_df |>
  dplyr::mutate(
    is_diff_species = .data$Species_label %in% marker_species &
      .data$Species_label != "Others"
  )

## ---------- 构建桑基图数据（wide 接口） -------------------------------------

pathway_short <- function(id, desc) {
  if (id == "ko02040") {
    return("Flagellar assembly")
  } else if (id == "ko02030") {
    return("Bacterial chemotaxis")
  }
  desc
}

sankey_df <- top_df |>
  dplyr::mutate(
    Pathway     = vapply(seq_len(n()), function(i) pathway_short(pathway_id[i], pathway_desc[i]), character(1)),
    KO_label    = glue("{ko_id} ({Pathway})"),
    Species_lab = Species_label
  )

diff_species_short <- marker_species

## ---------- 绘制桑基图（axis1/axis2 wide 模式） -----------------------------

p <- ggplot(
  sankey_df,
  aes(
    axis1 = KO_label,
    axis2 = Species_lab,
    y     = contrib
  )
) +
  ggalluvial::geom_alluvium(
    aes(fill = Pathway),
    alpha   = 0.7,
    colour  = "grey70",
    linewidth = 0.2
  ) +
  ggalluvial::geom_stratum(
    colour   = "grey30",
    linewidth = 0.3
  ) +
  # 物种轴标签：差异物种标红，其余黑色
  geom_text(
    stat = "stratum",
    aes(
      label = after_stat(stratum),
      color = after_stat(stratum %in% diff_species_short & stratum != "Others")
    ),
    size = 3.0
  ) +
  scale_fill_brewer(palette = "Set2", name = "Pathway") +
  scale_color_manual(
    values = c(`TRUE` = "red", `FALSE` = "black"),
    guide  = "none"
  ) +
  labs(
    x = NULL,
    y = "Relative contribution",
    title = "Cell motility pathways (ko02040 & ko02030)\nKO → species contributions in PoorPlaque group"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(size = 9, face = "bold"),
    axis.text.y     = element_text(size = 8),
    plot.title      = element_text(size = 11, face = "bold", hjust = 0),
    legend.position = "right",
    legend.title    = element_text(size = 9),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(5.5, 20, 5.5, 5.5)
  )

## ---------- 导出高清 PDF + JPG -----------------------------------------------

pdf_file <- file.path(out_dir, "Fig3_cell_motility_sankey.pdf")
jpg_file <- file.path(out_dir, "Fig3_cell_motility_sankey.jpg")

Cairo::CairoPDF(
  file   = pdf_file,
  width  = 6.5,
  height = 4.0,
  family = "sans"
)
print(p)
dev.off()

Cairo::CairoJPEG(
  file   = jpg_file,
  width  = 6.5,
  height = 4.0,
  units  = "in",
  dpi    = 600,
  family = "sans"
)
print(p)
dev.off()

message("[INFO] Cell motility KO→species sankey saved:")
message(glue("  - {pdf_file}"))
message(glue("  - {jpg_file}"))

