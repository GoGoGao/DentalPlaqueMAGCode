## Packaged for revise/; source: project/01this_study_dataset/sgb/compare_genomics_all/step4_VF_ARG.R
## Figures: FigS5
## Note: VF/ARG counts
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "FigS5_VF_ARG.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

# Step 4: VFDB / CARD 比较基因组分析（优化版）
#
# 分析要点:
#   - 每个 VF 大类 / 每个 ARG：在 dual-sig MAG 子集中，按 enrichment 组做 MAG 水平 presence（>0）的 Fisher 检验
#   - 堆叠柱：展示各组内「至少携带该特征的 MAG 比例」（prevalence），而非原始 count 加和
#   - 热图 fig4E：仅保留组间有信息量的特征（Fisher padj<0.15 或 组特异检出 或 |Δprevalence|>=0.2）
#   - fig4F/fig4G：按 -log10(padj) 展示最判别性的 VF 大类与 ARG（单侧富集方向）
#
# 输出目录: 05_VF_ARG/

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  library(ComplexHeatmap)
  library(circlize)
  library(Cairo)
  library(patchwork)
  library(scales)
})

options(stringsAsFactors = FALSE)
set.seed(42)

work_dir <- file.path(.pkg_root, "01this_study_dataset/sgb/compare_genomics_all"
source(file.path(work_dir, "export_helper.R"))

out_dir  <- file.path(work_dir, "05_VF_ARG")
fig_dir  <- file.path(work_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

GROUP_COLORS <- c("Normal" = "#00A087", "PoorPlaque" = "#E64B35")

# ---- 1. Load Data ----
screen_df <- read.table(
  file.path(work_dir, "01_MAG_screening/dual_sig_MAGs.tsv"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE
)
normal_mags     <- screen_df$genome[screen_df$enrichment == "Normal"]
poorplaque_mags <- screen_df$genome[screen_df$enrichment == "PoorPlaque"]
mag_group       <- setNames(screen_df$enrichment, screen_df$genome)

read_matrix <- function(path) {
  read.table(path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE,
             stringsAsFactors = FALSE, quote = "", fill = TRUE, comment.char = "")
}

vfdb_mat <- read_matrix(file.path(work_dir, "total.vfdb.xls"))
card_mat <- read_matrix(file.path(work_dir, "total.card.xls"))

all_mags   <- c(normal_mags, poorplaque_mags)
vfdb_sub   <- vfdb_mat[, intersect(all_mags, colnames(vfdb_mat)), drop = FALSE]
card_sub   <- card_mat[, intersect(all_mags, colnames(card_mat)), drop = FALSE]

vf_bin  <- (vfdb_sub > 0) * 1L
arg_bin <- (card_sub > 0) * 1L

n_N <- length(intersect(normal_mags,     colnames(vf_bin)))
n_P <- length(intersect(poorplaque_mags, colnames(vf_bin)))

# ---- 2. Fisher 检验：逐特征（与差异 KO 相同的 MAG×组 2×2 逻辑）----
fisher_binary_rows <- function(bin_mat, mags_N, mags_PP, label) {
  mN <- intersect(mags_N, colnames(bin_mat))
  mP <- intersect(mags_PP, colnames(bin_mat))
  nn <- length(mN)
  np <- length(mP)
  feats <- rownames(bin_mat)
  res <- lapply(seq_along(feats), function(i) {
    r <- feats[i]
    a <- sum(bin_mat[r, mN, drop = TRUE] > 0)
    b <- nn - a
    c <- sum(bin_mat[r, mP, drop = TRUE] > 0)
    d <- np - c
    tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    ft <- suppressWarnings(fisher.test(tab, alternative = "two.sided"))
    prop_N <- a / max(nn, 1L)
    prop_P <- c / max(np, 1L)
    data.frame(
      feature = r,
      n_present_Normal     = a,
      n_present_PoorPlaque = c,
      n_MAG_Normal         = nn,
      n_MAG_PoorPlaque     = np,
      prevalence_Normal     = prop_N,
      prevalence_PoorPlaque = prop_P,
      delta_prev            = prop_P - prop_N,
      odds_ratio            = unname(ft$estimate),
      pval                  = ft$p.value,
      stringsAsFactors      = FALSE
    )
  })
  df <- bind_rows(res) %>%
    mutate(
      padj = p.adjust(pval, method = "BH"),
      enriched_in = case_when(
        padj < 0.05 & prevalence_Normal > prevalence_PoorPlaque ~ "Normal",
        padj < 0.05 & prevalence_PoorPlaque > prevalence_Normal ~ "PoorPlaque",
        TRUE ~ "NS"
      ),
      feature_type = label
    ) %>%
    arrange(pval)
  df
}

vf_fisher  <- fisher_binary_rows(vf_bin,  normal_mags, poorplaque_mags, "VF")
arg_fisher <- fisher_binary_rows(arg_bin, normal_mags, poorplaque_mags, "ARG")

write.table(vf_fisher,  file.path(out_dir, "VF_category_MAG_fisher.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
write.table(arg_fisher, file.path(out_dir, "ARG_gene_MAG_fisher.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

# 高价值结果表：显著或组特异或强流行率差
flag_valuable <- function(df) {
  df %>%
    mutate(
      group_exclusive = (n_present_Normal > 0 & n_present_PoorPlaque == 0) |
        (n_present_PoorPlaque > 0 & n_present_Normal == 0),
      strong_delta    = abs(delta_prev) >= 0.2
    ) %>%
    filter(padj < 0.15 | group_exclusive | strong_delta) %>%
    arrange(padj)
}

vf_valuable  <- flag_valuable(vf_fisher)
arg_valuable <- flag_valuable(arg_fisher)
write.table(vf_valuable,  file.path(out_dir, "VF_high_value_features.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
write.table(arg_valuable, file.path(out_dir, "ARG_high_value_features.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 3. per-MAG 计数（箱线图）----
vf_count_df <- tibble(
  genome   = colnames(vf_bin),
  VF_count = colSums(vf_bin > 0),
  group    = factor(mag_group[genome], levels = c("Normal", "PoorPlaque"))
)

arg_count_df <- tibble(
  genome    = colnames(arg_bin),
  ARG_count = colSums(arg_bin > 0),
  group     = factor(mag_group[genome], levels = c("Normal", "PoorPlaque"))
)

write.table(vf_count_df,  file.path(out_dir, "fig4A_VF_count_per_MAG.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(arg_count_df, file.path(out_dir, "fig4B_ARG_count_per_MAG.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 4A–4B 箱线图 ----
wilcox_safe <- function(data, v, g) {
  x <- data[[v]][data[[g]] == "Normal"]
  y <- data[[v]][data[[g]] == "PoorPlaque"]
  if (length(unique(x)) == 1L && length(unique(y)) == 1L && unique(x) == unique(y)) {
    return(list(p.value = 1))
  }
  suppressWarnings(wilcox.test(x, y))
}

p_vf  <- wilcox_safe(vf_count_df,  "VF_count",  "group")$p.value
p_arg <- wilcox_safe(arg_count_df, "ARG_count", "group")$p.value

p4A <- ggboxplot(vf_count_df, x = "group", y = "VF_count",
                 fill = "group",
                 add = "jitter", add.params = list(size = 2, alpha = 0.7, color = "grey30"),
                 width = 0.5, outlier.shape = NA) +
  scale_fill_manual(values = GROUP_COLORS, drop = FALSE) +
  labs(
    x = NULL, y = "VF categories (presence count)",
    title = "A  VF repertoire size per MAG",
    subtitle = sprintf("Wilcoxon P = %.3g", p_vf)
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x     = element_text(size = 11, face = "bold"),
    legend.position = "none",
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 10, color = "grey35")
  )

p4B <- ggboxplot(arg_count_df, x = "group", y = "ARG_count",
                 fill = "group",
                 add = "jitter", add.params = list(size = 2, alpha = 0.7, color = "grey30"),
                 width = 0.5, outlier.shape = NA) +
  scale_fill_manual(values = GROUP_COLORS, drop = FALSE) +
  labs(
    x = NULL, y = "ARG families (presence count)",
    title = "B  ARG repertoire size per MAG",
    subtitle = sprintf("Wilcoxon P = %.3g", p_arg)
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    axis.text.x     = element_text(size = 11, face = "bold"),
    legend.position = "none",
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 10, color = "grey35")
  )

# ---- 4C–4D：组内 MAG 流行率（prevalence）堆叠 ----
vf_prev_long <- vf_fisher %>%
  select(VF_category = feature, prevalence_Normal, prevalence_PoorPlaque) %>%
  pivot_longer(-VF_category, names_to = "grp", values_to = "prev") %>%
  mutate(
    group = factor(ifelse(grepl("Normal", grp), "Normal", "PoorPlaque"),
                   levels = c("Normal", "PoorPlaque"))
  )

vf_bar_df <- vf_prev_long %>%
  group_by(group) %>%
  mutate(prop = prev / sum(prev) * 100) %>%
  ungroup()

write.table(vf_bar_df, file.path(out_dir, "fig4C_VF_prevalence_stack_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

n_vf_cats <- n_distinct(vf_bar_df$VF_category)
vf_pal <- colorRampPalette(c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
                              "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
                              "#9C755F", "#BAB0AC"))(n_vf_cats)
names(vf_pal) <- sort(unique(vf_bar_df$VF_category))

p4C <- ggplot(vf_bar_df, aes(x = group, y = prop, fill = VF_category)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = vf_pal, name = "VF category") +
  scale_x_discrete(labels = c("Normal" = "Normal-\nenriched", "PoorPlaque" = "PoorPlaque-\nenriched")) +
  labs(x = NULL, y = "Share of total prevalence (%)",
       title = "C  VF: relative contribution of categories",
       subtitle = "Within each group, category weight ~ MAG prevalence") +
  theme_classic(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.subtitle   = element_text(size = 9, color = "grey40"),
    legend.text     = element_text(size = 7),
    legend.key.size = unit(0.35, "cm")
  )

arg_prev_long <- arg_fisher %>%
  filter(n_present_Normal + n_present_PoorPlaque > 0) %>%
  select(ARG = feature, prevalence_Normal, prevalence_PoorPlaque) %>%
  pivot_longer(-ARG, names_to = "grp", values_to = "prev") %>%
  mutate(
    group = factor(ifelse(grepl("Normal", grp), "Normal", "PoorPlaque"),
                   levels = c("Normal", "PoorPlaque"))
  )

arg_bar_df <- arg_prev_long %>%
  group_by(group) %>%
  mutate(prop = prev / sum(prev) * 100) %>%
  ungroup()

write.table(arg_bar_df, file.path(out_dir, "fig4D_ARG_prevalence_stack_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

n_arg <- n_distinct(arg_bar_df$ARG)
arg_pal <- colorRampPalette(c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                               "#F39B7F", "#8491B4", "#91D1C2", "#DC0000"))(max(n_arg, 1))
names(arg_pal) <- sort(unique(arg_bar_df$ARG))

p4D <- ggplot(arg_bar_df, aes(x = group, y = prop, fill = ARG)) +
  geom_col(width = 0.62, color = "white", linewidth = 0.25) +
  scale_fill_manual(values = arg_pal, name = "ARG") +
  scale_x_discrete(labels = c("Normal" = "Normal-\nenriched", "PoorPlaque" = "PoorPlaque-\nenriched")) +
  labs(x = NULL, y = "Share of total prevalence (%)",
       title = "D  ARG: relative contribution",
       subtitle = "Rare ARGs pooled in legend; see ARG_high_value_features.tsv") +
  theme_classic(base_size = 11) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.subtitle   = element_text(size = 9, color = "grey40"),
    legend.text     = element_text(size = 6.5),
    legend.key.size = unit(0.3, "cm")
  )

# ---- 4F–4G：判别性条形图（Fisher padj）----
top_disc <- function(df, k = 12) {
  df %>%
    filter(n_present_Normal + n_present_PoorPlaque > 0) %>%
    mutate(
      neglog = -log10(pmax(padj, 1e-300)),
      dir    = ifelse(delta_prev >= 0, "PoorPlaque", "Normal")
    ) %>%
    arrange(padj) %>%
    slice_head(n = k) %>%
    mutate(feature = reorder(feature, -neglog))
}

plot_disc <- function(df, ttl) {
  udir <- unique(as.character(df$dir))
  cols <- GROUP_COLORS[names(GROUP_COLORS) %in% udir]
  ggplot(df, aes(x = neglog, y = feature, fill = dir)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = cols, drop = FALSE, name = "Higher prevalence") +
    labs(x = expression(-log[10](italic(P)[adj])), y = NULL, title = ttl) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "bottom"
    )
}

vf_top  <- top_disc(vf_fisher,  12)
arg_top <- top_disc(arg_fisher, min(15, nrow(arg_fisher)))

write.table(vf_top,  file.path(out_dir, "fig4F_VF_discriminative_bar_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(arg_top, file.path(out_dir, "fig4G_ARG_discriminative_bar_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p4F <- if (nrow(vf_top) > 0) plot_disc(vf_top, "F  Most discriminative VF categories") else NULL
p4G <- if (nrow(arg_top) > 0) plot_disc(arg_top, "G  Most discriminative ARGs") else NULL

# ---- 4E：精简热图 ----
pick_heatmap_rows <- function(vf_f, arg_f, vf_b, arg_b) {
  vf_keep <- vf_f %>%
    filter(
      padj < 0.15 | abs(delta_prev) >= 0.2 |
        (n_present_Normal > 0 & n_present_PoorPlaque == 0) |
        (n_present_PoorPlaque > 0 & n_present_Normal == 0)
    ) %>%
    pull(feature)
  arg_keep <- arg_f %>%
    filter(
      padj < 0.2 | abs(delta_prev) >= 0.15 |
        (n_present_Normal > 0 & n_present_PoorPlaque == 0) |
        (n_present_PoorPlaque > 0 & n_present_Normal == 0)
    ) %>%
    pull(feature)
  vf_keep  <- intersect(vf_keep,  rownames(vf_b))
  arg_keep <- intersect(arg_keep, rownames(arg_b))
  list(vf = vf_keep, arg = arg_keep)
}

picked <- pick_heatmap_rows(vf_fisher, arg_fisher, vf_bin, arg_bin)
if (length(picked$vf) + length(picked$arg) == 0) {
  picked$vf  <- rownames(vf_bin)
  picked$arg <- rownames(arg_bin)
}

vf_h  <- vf_bin[picked$vf, , drop = FALSE]
arg_h <- arg_bin[picked$arg, , drop = FALSE]

vf_l  <- as.data.frame(vf_h)  %>% rownames_to_column("Feature") %>% mutate(Type = "VF")
arg_l <- as.data.frame(arg_h) %>% rownames_to_column("Feature") %>% mutate(Type = "ARG")
combined_long <- bind_rows(vf_l, arg_l)
combined_mat <- combined_long %>%
  select(-Type) %>%
  column_to_rownames("Feature") %>%
  as.matrix()

col_order <- c(intersect(normal_mags, colnames(combined_mat)),
               intersect(poorplaque_mags, colnames(combined_mat)))
combined_mat <- combined_mat[, col_order, drop = FALSE]

write.table(
  data.frame(Feature = rownames(combined_mat), combined_mat, check.names = FALSE),
  file.path(out_dir, "fig4E_heatmap_data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  tibble(Feature = c(picked$vf, picked$arg), Type = c(rep("VF", length(picked$vf)), rep("ARG", length(picked$arg)))),
  file.path(out_dir, "fig4E_heatmap_row_selection.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

col_anno_df <- tibble(
  Group = c(rep("Normal", length(intersect(normal_mags, col_order))),
            rep("PoorPlaque", length(intersect(poorplaque_mags, col_order))))
)

ha_col <- HeatmapAnnotation(
  Group = col_anno_df$Group,
  col   = list(Group = GROUP_COLORS),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
)

feat_type <- combined_long$Type
names(feat_type) <- combined_long$Feature
feat_type <- feat_type[rownames(combined_mat)]

ha_row <- rowAnnotation(
  Type = feat_type,
  col  = list(Type = c("VF" = "#E64B35", "ARG" = "#4DBBD5")),
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
)

ht <- Heatmap(
  combined_mat,
  name            = "Presence",
  col             = colorRamp2(c(0, 1), c("white", "#2C3E50")),
  cluster_columns = FALSE,
  cluster_rows    = TRUE,
  column_split    = col_anno_df$Group,
  column_gap      = unit(3, "mm"),
  top_annotation  = ha_col,
  left_annotation = ha_row,
  row_names_gp    = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 6.5),
  heatmap_legend_param = list(title = "Presence", title_gp = gpar(fontsize = 9, fontface = "bold"))
)

draw_ht4e <- function() draw(ht, merge_legend = TRUE, heatmap_legend_side = "right")
for (ext in c("pdf", "svg", "png")) {
  f <- file.path(out_dir, paste0("fig4E_VF_ARG_heatmap.", ext))
  tryCatch({
    switch(ext,
      pdf = { Cairo::CairoPDF(f, 14, 9);                    draw_ht4e(); invisible(dev.off()) },
      svg = { Cairo::CairoSVG(f, 14, 9);                    draw_ht4e(); invisible(dev.off()) },
      png = { Cairo::CairoPNG(f, 14, 9, units = "in", dpi = 300); draw_ht4e(); invisible(dev.off()) }
    )
    message("[OK] fig4E ", ext)
  }, error = function(e) message("[ERR] fig4E ", ext, ": ", e$message))
}

# ---- 中文补充表（基础组间差异：负担 Wilcoxon + 特征 Fisher）----
zh_from_vf <- vf_fisher %>%
  transmute(
    `VF功能大类` = feature,
    `Normal组检出MAG数` = n_present_Normal,
    `PoorPlaque组检出MAG数` = n_present_PoorPlaque,
    `Normal组MAG总数` = n_MAG_Normal,
    `PoorPlaque组MAG总数` = n_MAG_PoorPlaque,
    `Normal组流行率` = round(prevalence_Normal, 4),
    `PoorPlaque组流行率` = round(prevalence_PoorPlaque, 4),
    `流行率差_PP减N` = round(delta_prev, 4),
    `OR` = round(odds_ratio, 4),
    `P值` = pval,
    `FDR_BH` = padj,
    `显著富集组_FDR0.05` = enriched_in
  )
write.table(zh_from_vf,
            file.path(out_dir, "补充表S1_VFDB组间Fisher检验_中文.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

zh_from_arg <- arg_fisher %>%
  transmute(
    `CARD条目` = feature,
    `Normal组检出MAG数` = n_present_Normal,
    `PoorPlaque组检出MAG数` = n_present_PoorPlaque,
    `Normal组MAG总数` = n_MAG_Normal,
    `PoorPlaque组MAG总数` = n_MAG_PoorPlaque,
    `Normal组流行率` = round(prevalence_Normal, 4),
    `PoorPlaque组流行率` = round(prevalence_PoorPlaque, 4),
    `流行率差_PP减N` = round(delta_prev, 4),
    `OR` = round(odds_ratio, 4),
    `P值` = pval,
    `FDR_BH` = padj,
    `显著富集组_FDR0.05` = enriched_in
  )
write.table(zh_from_arg,
            file.path(out_dir, "补充表S2_CARD组间Fisher检验_中文.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

burden_tbl <- tibble(
  `指标` = c(
    "VF功能大类数（每MAG检出类别数）",
    "CARD抗性相关条目数（每MAG检出数）"
  ),
  `统计方法` = c("Wilcoxon秩和检验", "Wilcoxon秩和检验"),
  `分组` = c("Normal-enriched vs PoorPlaque-enriched MAG", "同上"),
  `P值` = c(p_vf, p_arg),
  `说明` = c(
    "基于 total.vfdb.xls 列联；仅含双显著差异MAG",
    "基于 total.card.xls；CARD 矩阵未覆盖的 MAG 不计入该侧分母"
  )
)
write.table(burden_tbl,
            file.path(out_dir, "补充表S0_VF与ARG负担组间比较_中文.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")

# ---- Export ggplot panels ----
save_plot(p4A, "fig4A_VF_boxplot",  out_dir, w = 4.8, h = 5.2)
save_plot(p4B, "fig4B_ARG_boxplot", out_dir, w = 4.8, h = 5.2)
save_plot(p4C, "fig4C_VF_bar",      out_dir, w = 6.5, h = 5.8)
save_plot(p4D, "fig4D_ARG_bar",     out_dir, w = 7.5, h = 6.0)
if (!is.null(p4F)) save_plot(p4F, "fig4F_VF_discriminative",  out_dir, w = 7.5, h = 6.5)
if (!is.null(p4G)) save_plot(p4G, "fig4G_ARG_discriminative", out_dir, w = 7.5, h = 6.5)

p4_comb <- (p4A | p4B) / (p4C | p4D) + plot_layout(heights = c(1, 1.05))
save_plot(p4_comb, "Fig4_VFARG_combined", fig_dir, w = 12.5, h = 11)

if (!is.null(p4F) && !is.null(p4G)) {
  p4_fg <- p4F | p4G
  save_plot(p4_fg, "Fig4_FG_discriminative", fig_dir, w = 14, h = 7)
}

cat("[DONE] Step 4 complete. Outputs in:", out_dir, "\n")
