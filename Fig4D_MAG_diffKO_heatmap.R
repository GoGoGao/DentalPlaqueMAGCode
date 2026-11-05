## Packaged for revise/; source: project/01this_study_dataset/sgb/metab_corr/paper_result/03_RobustDA_diffKO_heatmap/scripts/03_diffKO_fisher_heatmap.R
## Figures: Fig4D
## Note: differential KO heatmap
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig4D_MAG_diffKO_heatmap.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# Differential KO (Fisher) + presence/absence heatmap for Robust_DA MAGs
# Adapted from compare_genomics_all/step3_KEGG_enrichment.R (Fisher) + step6_heatmap.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(Cairo)
})

options(stringsAsFactors = FALSE)
set.seed(42)

base_dir   <- file.path(.pkg_root, "01this_study_dataset/sgb/metab_corr"
comp_dir   <- file.path(dirname(base_dir), "compare_genomics_all")
module_dir <- file.path(base_dir, "paper_result/03_RobustDA_diffKO_heatmap")
out_fig    <- file.path(module_dir, "figures")
out_tab    <- file.path(module_dir, "tables")
out_data   <- file.path(module_dir, "data")
dir.create(out_fig,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
source(file.path(comp_dir, "phylum_palette.R"))

GROUP_COLORS <- c("Normal" = "#00A087", "PoorPlaque" = "#E64B35")

# ---- 1. MAG set ----
screen_df <- read.table(
  file.path(base_dir, "paper_result/00_RobustDA_MAG_set/tables/RobustDA_MAGs.tsv"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)
normal_mags     <- screen_df$genome[screen_df$enrichment == "Normal"]
poorplaque_mags <- screen_df$genome[screen_df$enrichment == "PoorPlaque"]
mag_group  <- setNames(screen_df$enrichment, screen_df$genome)
mag_phylum <- setNames(as.character(screen_df$phylum), screen_df$genome)
mag_genus  <- setNames(as.character(screen_df$genus),  screen_df$genome)
all_mags   <- c(normal_mags, poorplaque_mags)
message("[INFO] Robust_DA MAGs: ", length(all_mags),
        " (N=", length(normal_mags), ", PP=", length(poorplaque_mags), ")")

# ---- 2. KO matrix + Fisher ----
ko_mat <- read.table(
  file.path(comp_dir, "total.KO.xls"),
  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE
)
ko_bin <- (ko_mat > 0) * 1L

all_mags_in <- intersect(all_mags, colnames(ko_bin))
n_normal <- length(intersect(normal_mags, colnames(ko_bin)))
n_poor   <- length(intersect(poorplaque_mags, colnames(ko_bin)))
ko_sub <- ko_bin[, all_mags_in, drop = FALSE]
normal_sub <- ko_sub[, intersect(normal_mags, colnames(ko_sub)), drop = FALSE]
poor_sub   <- ko_sub[, intersect(poorplaque_mags, colnames(ko_sub)), drop = FALSE]

n_norm_present <- rowSums(normal_sub > 0)
n_poor_present <- rowSums(poor_sub > 0)
n_norm_absent  <- n_normal - n_norm_present
n_poor_absent  <- n_poor   - n_poor_present

fisher_pvals <- numeric(nrow(ko_sub))
odds_ratio   <- numeric(nrow(ko_sub))
for (i in seq_len(nrow(ko_sub))) {
  mat2x2 <- matrix(c(n_norm_present[i], n_norm_absent[i],
                     n_poor_present[i], n_poor_absent[i]),
                   nrow = 2, byrow = TRUE)
  ft <- fisher.test(mat2x2, alternative = "two.sided")
  fisher_pvals[i] <- ft$p.value
  odds_ratio[i]   <- unname(ft$estimate)
}

fisher_df <- data.frame(
  KO = rownames(ko_sub),
  n_Normal = as.integer(n_norm_present),
  n_PoorPlaque = as.integer(n_poor_present),
  pval = fisher_pvals,
  padj = p.adjust(fisher_pvals, method = "BH"),
  odds_ratio = odds_ratio,
  stringsAsFactors = FALSE
) %>%
  mutate(
    prop_Normal = n_Normal / n_normal,
    prop_PoorPlaque = n_PoorPlaque / n_poor,
    enriched_in = case_when(
      padj < 0.05 & prop_Normal > prop_PoorPlaque ~ "Normal",
      padj < 0.05 & prop_PoorPlaque > prop_Normal ~ "PoorPlaque",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(padj)

write.table(fisher_df, file.path(out_tab, "fisher_KO_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fisher_df %>% filter(enriched_in != "NS"),
            file.path(out_tab, "fisher_KO_significant.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("[INFO] Diff KO padj<0.05: ", sum(fisher_df$padj < 0.05),
        " | Normal=", sum(fisher_df$enriched_in == "Normal"),
        " | PoorPlaque=", sum(fisher_df$enriched_in == "PoorPlaque"))

# ---- 3. Select heatmap KOs: top 25 per side by padj ----
k_per_side <- 25L
ko_N <- fisher_df %>% filter(enriched_in == "Normal") %>%
  arrange(padj, pval) %>% slice_head(n = k_per_side)
ko_P <- fisher_df %>% filter(enriched_in == "PoorPlaque") %>%
  arrange(padj, pval) %>% slice_head(n = k_per_side)
manifest <- bind_rows(
  ko_N %>% mutate(rank_in_side = row_number(), side = "Normal"),
  ko_P %>% mutate(rank_in_side = row_number(), side = "PoorPlaque")
) %>% distinct(KO, .keep_all = TRUE)
write.table(manifest, file.path(out_tab, "fig6_KO_selection_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

top_ko <- unique(c(ko_N$KO, ko_P$KO))
top_ko_avail <- intersect(top_ko, rownames(ko_sub))
message("[INFO] Heatmap KO rows: ", length(top_ko_avail))

# ---- 4. VF/ARG counts for annotation ----
vfdb_mat <- read.table(file.path(comp_dir, "total.vfdb.xls"),
                       sep = "\t", header = TRUE, row.names = 1,
                       check.names = FALSE, stringsAsFactors = FALSE)
card_mat <- read.table(file.path(comp_dir, "total.card.xls"),
                       sep = "\t", header = TRUE, row.names = 1,
                       check.names = FALSE, stringsAsFactors = FALSE)

shorten_ko_label <- function(ko_id) {
  parts <- strsplit(ko_id, ":")[[1]]
  k_num <- parts[1]
  if (length(parts) > 1) {
    gene_name <- strsplit(trimws(parts[2]), ",")[[1]][1]
    paste0(k_num, " [", trimws(gene_name), "]")
  } else k_num
}

ko_heat <- ((ko_sub[top_ko_avail, all_mags[all_mags %in% colnames(ko_sub)], drop = FALSE]) > 0) * 1L
rownames(ko_heat) <- sapply(rownames(ko_heat), shorten_ko_label)

sort_by_genus <- function(mags) {
  g <- mag_genus[mags]
  g[is.na(g)] <- "ZZZ"
  mags[order(g)]
}
col_order <- c(sort_by_genus(intersect(normal_mags, colnames(ko_heat))),
               sort_by_genus(intersect(poorplaque_mags, colnames(ko_heat))))
ko_heat <- ko_heat[, col_order, drop = FALSE]

vf_counts  <- colSums(vfdb_mat[, intersect(col_order, colnames(vfdb_mat)), drop = FALSE] > 0)
arg_counts <- colSums(card_mat[, intersect(col_order, colnames(card_mat)), drop = FALSE] > 0)
vf_vec  <- vf_counts[col_order];  vf_vec[is.na(vf_vec)] <- 0
arg_vec <- arg_counts[col_order]; arg_vec[is.na(arg_vec)] <- 0

ht_data <- data.frame(MAG = col_order, as.data.frame(t(ko_heat)),
                      VF_count = as.numeric(vf_vec),
                      ARG_count = as.numeric(arg_vec),
                      check.names = FALSE)
write.table(ht_data, file.path(out_data, "fig6_heatmap_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 5. ComplexHeatmap ----
group_vec  <- mag_group[col_order]
phylum_vec <- as.character(mag_phylum[col_order])
phylum_vec[is.na(phylum_vec) | !nzchar(phylum_vec)] <- "Unclassified"
phylum_levels <- sort(unique(phylum_vec))
PHYLUM_COLORS <- build_phylum_palette(phylum_levels)
write.table(
  data.frame(phylum_latin = phylum_levels,
             fill_hex = unname(PHYLUM_COLORS[phylum_levels])),
  file.path(out_tab, "fig6_phylum_color_map.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

ko_enrich <- setNames(
  fisher_df$enriched_in[match(top_ko_avail, fisher_df$KO)],
  sapply(top_ko_avail, shorten_ko_label)
)
ko_enrich[is.na(ko_enrich)] <- "NS"

ha_row_anno <- rowAnnotation(
  KO_enrich = ko_enrich[rownames(ko_heat)],
  col = list(KO_enrich = c("Normal" = "#00A087", "PoorPlaque" = "#E64B35", "NS" = "#BBBBBB")),
  show_annotation_name = TRUE,
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  show_legend = TRUE
)

ha_left <- HeatmapAnnotation(
  which = "column",
  Group = group_vec,
  Phylum = factor(phylum_vec, levels = phylum_levels),
  col = list(Group = GROUP_COLORS, Phylum = PHYLUM_COLORS),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 11, fontface = "bold"),
  gap = unit(1.5, "mm"),
  gp = gpar(col = NA)
)

ha_right <- HeatmapAnnotation(
  which = "column",
  VF_count = anno_barplot(as.numeric(vf_vec), bar_width = 0.85,
                          gp = gpar(fill = "#E64B35", col = NA),
                          height = unit(2.0, "cm")),
  ARG_count = anno_barplot(as.numeric(arg_vec), bar_width = 0.85,
                           gp = gpar(fill = "#4DBBD5", col = NA),
                           height = unit(1.5, "cm")),
  show_annotation_name = TRUE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
)

col_split <- factor(group_vec, levels = c("Normal", "PoorPlaque"))

ht <- Heatmap(
  ko_heat,
  name = "KO presence",
  col = colorRamp2(c(0, 1), c("grey96", "#2C3E50")),
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 9),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 60,
  column_split = col_split,
  column_gap = unit(4, "mm"),
  column_title_gp = gpar(fontsize = 13, fontface = "bold",
                         col = c("#00A087", "#E64B35")),
  top_annotation = ha_left,
  bottom_annotation = ha_right,
  left_annotation = ha_row_anno,
  row_title = paste0("Differential KOs (top-by-padj, n=", nrow(ko_heat), ")"),
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "KO presence\n(0=absent, 1=present)",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    legend_height = unit(3, "cm")
  )
)

w_ht <- max(14, 0.28 * ncol(ko_heat) + 4)
h_ht <- max(12, 0.18 * nrow(ko_heat) + 4)

draw_ht <- function() {
  draw(ht, merge_legend = TRUE, heatmap_legend_side = "right",
       annotation_legend_side = "right")
}

for (ext in c("pdf", "svg", "png", "jpg")) {
  f <- file.path(out_fig, paste0("Fig_diffKO_heatmap.", ext))
  tryCatch({
    if (ext == "pdf") {
      Cairo::CairoPDF(f, width = w_ht, height = h_ht); draw_ht(); invisible(dev.off())
    } else if (ext == "svg") {
      Cairo::CairoSVG(f, width = w_ht, height = h_ht); draw_ht(); invisible(dev.off())
    } else if (ext == "png") {
      Cairo::CairoPNG(f, width = w_ht, height = h_ht, units = "in", dpi = 300)
      draw_ht(); invisible(dev.off())
    } else {
      Cairo::CairoJPEG(f, width = w_ht, height = h_ht, units = "in", dpi = 300, quality = 95)
      draw_ht(); invisible(dev.off())
    }
    message("[OK] ", basename(f))
  }, error = function(e) message("[ERR] ", ext, ": ", e$message))
}

writeLines(
  c(
    "差异 KO + 热图（Robust_DA MAG）",
    "MAG 集合: paper_result/00_RobustDA_MAG_set (Robust_DA == TRUE)",
    "矩阵: compare_genomics_all/total.KO.xls 二值化",
    "检验: Fisher 精确检验 + BH padj; padj<0.05 且组间比例更高一侧为 enriched_in",
    "热图行: Normal/PoorPlaque 各 top25 padj KO（合并去重）",
    paste0("n_MAG=", length(all_mags_in), "; n_diffKO=", sum(fisher_df$enriched_in != "NS"),
           "; n_heatmap_rows=", length(top_ko_avail))
  ),
  file.path(module_dir, "README_CRITERIA.txt")
)

message("[DONE] Diff KO heatmap → ", module_dir)
