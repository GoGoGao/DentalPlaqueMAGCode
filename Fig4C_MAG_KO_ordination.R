## Packaged for revise/; source: project/01this_study_dataset/sgb/metab_corr/paper_result/04_RobustDA_KO_ordination/scripts/04_KO_PCoA.R
## Figures: Fig4C
## Note: MAG KO Jaccard PCoA
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig4C_MAG_KO_ordination.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# KO-level PCoA (Jaccard) + PERMANOVA for Robust_DA MAGs
# Adapted from compare_genomics_all/step5_ordination.R (Fig5A KO only)

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(Cairo)
})

options(stringsAsFactors = FALSE)
set.seed(42)

base_dir   <- file.path(.pkg_root, "01this_study_dataset/sgb/metab_corr"
comp_dir   <- file.path(dirname(base_dir), "compare_genomics_all")
module_dir <- file.path(base_dir, "paper_result/04_RobustDA_KO_ordination")
out_fig    <- file.path(module_dir, "figures")
out_tab    <- file.path(module_dir, "tables")
out_data   <- file.path(module_dir, "data")
dir.create(out_fig,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
source(file.path(base_dir, "paper_result/_shared/export_utils.R"))

GROUP_COLORS <- c("Normal" = "#00A087", "PoorPlaque" = "#E64B35")
GROUP_SHAPES <- c("Normal" = 21, "PoorPlaque" = 24)

screen_df <- read.table(
  file.path(base_dir, "paper_result/00_RobustDA_MAG_set/tables/RobustDA_MAGs.tsv"),
  sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE
)
normal_mags     <- screen_df$genome[screen_df$enrichment == "Normal"]
poorplaque_mags <- screen_df$genome[screen_df$enrichment == "PoorPlaque"]
mag_group <- setNames(screen_df$enrichment, screen_df$genome)
all_mags  <- c(normal_mags, poorplaque_mags)
message("[INFO] Robust_DA MAGs: ", length(all_mags))

ko_mat <- read.table(
  file.path(comp_dir, "total.KO.xls"),
  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE
)

build_pa_matrix <- function(mat, mags) {
  avail <- intersect(mags, colnames(mat))
  sub <- mat[, avail, drop = FALSE]
  sub <- sub[rowSums(sub) > 0, , drop = FALSE]
  t((sub > 0) * 1L)
}

ko_pa <- build_pa_matrix(ko_mat, all_mags)
ko_pa <- ko_pa[intersect(rownames(ko_pa), all_mags), , drop = FALSE]
message("[INFO] KO PA: ", nrow(ko_pa), " MAGs x ", ncol(ko_pa), " KOs")

group_df <- data.frame(
  genome = rownames(ko_pa),
  group  = factor(mag_group[rownames(ko_pa)], levels = c("Normal", "PoorPlaque"))
)
dist_mat <- vegdist(ko_pa, method = "jaccard", binary = TRUE)
pcoa_res <- cmdscale(dist_mat, k = min(nrow(ko_pa) - 1, 10), eig = TRUE)
eig_vals <- pcoa_res$eig
eig_vals[eig_vals < 0] <- 0
var_exp  <- eig_vals / sum(eig_vals) * 100

coord_df <- data.frame(
  genome = rownames(ko_pa),
  PCoA1  = pcoa_res$points[, 1],
  PCoA2  = pcoa_res$points[, 2],
  group  = group_df$group,
  stringsAsFactors = FALSE
)
write.table(coord_df, file.path(out_data, "pcoa_coordinates_KO.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

perm_res <- adonis2(dist_mat ~ group, data = group_df, permutations = 999)
R2   <- round(perm_res$R2[1], 4)
pval <- perm_res$`Pr(>F)`[1]
perm_df <- as.data.frame(perm_res)
perm_df$Term <- rownames(perm_df)
write.table(perm_df, file.path(out_tab, "permanova_KO.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# binary PA matrix for reproducibility (MAG x KO may be large — save summary + full)
write.table(
  data.frame(MAG = rownames(ko_pa), group = group_df$group, n_KO = rowSums(ko_pa)),
  file.path(out_tab, "MAG_KO_richness.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

var1 <- round(var_exp[1], 1)
var2 <- round(var_exp[2], 1)
pval_label <- ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))
stat_label <- sprintf("PERMANOVA: R²=%.3f, P=%s", R2, pval_label)

p5A <- ggplot(coord_df, aes(x = PCoA1, y = PCoA2, fill = group, shape = group)) +
  stat_ellipse(aes(color = group), type = "t", level = 0.95,
               linewidth = 0.9, linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 4.5, color = "grey20", stroke = 0.5, alpha = 0.88) +
  scale_fill_manual(values = GROUP_COLORS,
                    labels = c("Normal" = "Normal-enriched",
                               "PoorPlaque" = "PoorPlaque-enriched"),
                    name = "Group") +
  scale_color_manual(values = GROUP_COLORS) +
  scale_shape_manual(values = GROUP_SHAPES,
                     labels = c("Normal" = "Normal-enriched",
                                "PoorPlaque" = "PoorPlaque-enriched"),
                     name = "Group") +
  annotate("text", x = Inf, y = Inf, label = stat_label,
           hjust = 1.05, vjust = 1.5, size = 4.5, color = "grey20",
           fontface = "italic") +
  labs(
    x = paste0("PCoA1 (", var1, "%)"),
    y = paste0("PCoA2 (", var2, "%)"),
    title = "PCoA: KO functional profile (Jaccard)",
    subtitle = paste0("Robust_DA MAGs (n=", nrow(ko_pa), ")")
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border     = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
    axis.text        = element_text(size = 12, color = "grey10"),
    axis.title       = element_text(size = 13, color = "grey10"),
    axis.line        = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 12, face = "bold"),
    legend.text      = element_text(size = 11),
    plot.title       = element_text(size = 15, face = "bold"),
    plot.subtitle    = element_text(size = 12, colour = "grey30"),
    plot.margin      = margin(10, 14, 8, 8)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 4.5)),
         shape = guide_legend(override.aes = list(size = 4.5)))

save_gg(p5A, file.path(out_fig, "Fig_KO_PCoA"), w = 7.5, h = 6)

message("[DONE] KO PCoA → ", module_dir)
message("[INFO] PERMANOVA R2=", R2, " p=", pval)
