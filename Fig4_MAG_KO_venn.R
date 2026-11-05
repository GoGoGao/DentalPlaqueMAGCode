## Packaged for revise/; source: project/01this_study_dataset/sgb/metab_corr/paper_result/02_RobustDA_KO_venn/scripts/02_KO_venn.R
## Figures: Fig4
## Note: KO venn
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig4_MAG_KO_venn.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# KO-level Venn: Normal-enriched vs PoorPlaque-enriched Robust_DA MAGs
# Adapted from compare_genomics_all/step2_venn.R (KO panel only)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggVennDiagram)
  library(UpSetR)
  library(Cairo)
})

options(stringsAsFactors = FALSE)
set.seed(42)

base_dir   <- file.path(.pkg_root, "01this_study_dataset/sgb/metab_corr"
comp_dir   <- file.path(dirname(base_dir), "compare_genomics_all")
module_dir <- file.path(base_dir, "paper_result/02_RobustDA_KO_venn")
out_fig    <- file.path(module_dir, "figures")
out_tab    <- file.path(module_dir, "tables")
out_data   <- file.path(module_dir, "data")
dir.create(out_fig,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
source(file.path(base_dir, "paper_result/_shared/export_utils.R"))

# Group palette (aligned with Panel A / RobustDA figures)
GROUP_COLORS <- c("Normal" = "#2E8B8B", "PoorPlaque" = "#B85450")
REGION_FILL  <- c(
  "1"   = "#8FC4C0",  # Normal only
  "2"   = "#E0A8A4",  # PoorPlaque only
  "1/2" = "#C9B6D9"   # Shared
)
EDGE_COLORS  <- c("1" = "#2E8B8B", "2" = "#B85450")

screen_df <- read.table(
  file.path(base_dir, "paper_result/00_RobustDA_MAG_set/tables/RobustDA_MAGs.tsv"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)
normal_mags     <- screen_df$genome[screen_df$enrichment == "Normal"]
poorplaque_mags <- screen_df$genome[screen_df$enrichment == "PoorPlaque"]
message("[INFO] Robust_DA MAGs: Normal=", length(normal_mags),
        " PoorPlaque=", length(poorplaque_mags))

ko_mat <- read.table(
  file.path(comp_dir, "total.KO.xls"),
  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE
)

extract_gene_set <- function(mat, mags) {
  available_mags <- intersect(mags, colnames(mat))
  if (length(available_mags) == 0) return(character(0))
  sub_mat <- mat[, available_mags, drop = FALSE]
  rownames(sub_mat)[rowSums(sub_mat > 0) > 0]
}

ko_Normal     <- extract_gene_set(ko_mat, normal_mags)
ko_PoorPlaque <- extract_gene_set(ko_mat, poorplaque_mags)
message("[INFO] KO union — Normal:", length(ko_Normal),
        " | PoorPlaque:", length(ko_PoorPlaque))

shared  <- intersect(ko_Normal, ko_PoorPlaque)
only_N  <- setdiff(ko_Normal, ko_PoorPlaque)
only_PP <- setdiff(ko_PoorPlaque, ko_Normal)

union_n <- length(only_N) + length(shared) + length(only_PP)
stats_df <- data.frame(
  Category = c("Normal_only", "Shared", "PoorPlaque_only"),
  Count    = c(length(only_N), length(shared), length(only_PP)),
  Percent  = round(100 * c(length(only_N), length(shared), length(only_PP)) / union_n, 2),
  stringsAsFactors = FALSE
)
write.table(stats_df, file.path(out_tab, "fig2_venn_stats_KO.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(stats_df, file.path(out_data, "Fig_KO_venn_region_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

safe_write <- function(vec, cat_name, fname) {
  df <- if (length(vec) == 0)
    data.frame(Feature = character(0), Category = character(0))
  else
    data.frame(Feature = vec, Category = cat_name)
  write.table(df, fname, sep = "\t", quote = FALSE, row.names = FALSE)
}
safe_write(only_N,  "Normal_only",     file.path(out_data, "KO_Normal_only.tsv"))
safe_write(only_PP, "PoorPlaque_only", file.path(out_data, "KO_PoorPlaque_only.tsv"))
safe_write(shared,  "Shared",          file.path(out_data, "KO_shared.tsv"))

write.table(
  data.frame(
    MAG_set = "Robust_DA",
    n_Normal_MAG = length(normal_mags),
    n_PoorPlaque_MAG = length(poorplaque_mags),
    KO_Normal = length(ko_Normal),
    KO_PoorPlaque = length(ko_PoorPlaque),
    Shared = length(shared),
    Normal_only = length(only_N),
    PoorPlaque_only = length(only_PP)
  ),
  file.path(out_tab, "fig2A_upset_KO_data.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Custom 2-set Venn: region fills by group; labels = count + % of union
venn_obj <- Venn(list(Normal = ko_Normal, PoorPlaque = ko_PoorPlaque))
vd <- process_data(venn_obj)
reg_edge  <- as.data.frame(venn_regionedge(vd))
set_edge  <- as.data.frame(venn_setedge(vd))
set_lab   <- as.data.frame(venn_setlabel(vd))
reg_lab   <- as.data.frame(venn_regionlabel(vd))
reg_lab$Percent <- 100 * reg_lab$count / sum(reg_lab$count)
reg_lab$lab <- sprintf("%d\n(%.1f%%)", reg_lab$count, reg_lab$Percent)

venn_ko <- ggplot() +
  geom_polygon(
    data = reg_edge, aes(X, Y, fill = id, group = id),
    color = NA, alpha = 0.92
  ) +
  geom_path(
    data = set_edge, aes(X, Y, color = id, group = id),
    linewidth = 1.35, lineend = "round"
  ) +
  geom_text(
    data = set_lab, aes(X, Y, label = name, color = id),
    size = 5.2, fontface = "bold", hjust = 1, nudge_x = -0.15
  ) +
  geom_text(
    data = reg_lab, aes(X, Y, label = lab),
    size = 4.8, lineheight = 0.95, fontface = "bold", color = "#1F1F1F"
  ) +
  scale_fill_manual(values = REGION_FILL, guide = "none") +
  scale_color_manual(values = EDGE_COLORS, guide = "none") +
  coord_equal(clip = "off") +
  scale_x_continuous(expand = expansion(mult = c(0.22, 0.08))) +
  scale_y_continuous(expand = expansion(mult = 0.08)) +
  labs(
    title = "KO gene sets (Robust_DA MAGs)",
    subtitle = sprintf("Union n=%d; labels = count (%% of union)", union_n)
  ) +
  theme_void(base_size = 14) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "black",
                              margin = margin(b = 4)),
    plot.subtitle = element_text(size = 11, colour = "grey35", hjust = 0.5,
                                 margin = margin(b = 8)),
    plot.margin = margin(10, 14, 10, 18)
  )

save_gg(venn_ko, file.path(out_fig, "Fig_KO_venn"), w = 5.8, h = 5.2)

png(file.path(out_fig, "Fig_KO_upset.png"), width = 2800, height = 2000, res = 300)
upset(fromList(list(Normal = ko_Normal, PoorPlaque = ko_PoorPlaque)),
      sets = c("Normal", "PoorPlaque"), order.by = "freq",
      sets.bar.color = unname(GROUP_COLORS[c("Normal", "PoorPlaque")]),
      main.bar.color = "#555555",
      text.scale = 1.8, keep.order = TRUE,
      mb.ratio = c(0.6, 0.4))
invisible(dev.off())
file.copy(file.path(out_fig, "Fig_KO_upset.png"),
          file.path(out_fig, "Fig_KO_upset.jpg"), overwrite = TRUE)

message("[DONE] KO Venn → ", module_dir)
