## Packaged for revise/; source: project/01this_study_dataset/sgb/metab_corr/paper_result/14_metab_PCA_ordination/scripts/14_metab_PCA.R
## Figures: Fig2A
## Note: metabolome PCA
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig2A_metab_PCA.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# =============================================================================
# Module 14: Overall metabolome PCA (A/B cohort)
# - Input: select metabolites (metab_abund_AB_select), group labels, Plaque_Index
# - PCA on z-scored metabolite matrix (samples × features)
# - PERMANOVA (Euclidean on scaled data); scree + loadings + score plot
# - Each figure → PDF/SVG/PNG/JPG + TSV; composite + standalone panels
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(vegan)
  library(svglite)
})

options(stringsAsFactors = FALSE)
set.seed(42)

base_dir   <- file.path(.pkg_root, "01this_study_dataset/sgb/metab_corr"
module_dir <- file.path(base_dir, "paper_result/14_metab_PCA_ordination")
out_fig    <- file.path(module_dir, "figures")
out_tab    <- file.path(module_dir, "tables")
out_data   <- file.path(module_dir, "data")
dir.create(out_fig,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
source(file.path(base_dir, "paper_result/_shared/export_utils.R"))

GROUP_COLORS <- c(Normal = "#00A087", PoorPlaque = "#E64B35")
GROUP_SHAPES <- c(Normal = 21, PoorPlaque = 24)
BASE_PT <- 13
N_PC_SCREE <- 10L
N_LOAD_TOP <- 15L  # top |loading| per PC for bar plot

save_pub <- function(plot, stem, width_mm = 160, height_mm = 140, dpi = 600) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  print(plot); invisible(dev.off())
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h)
  print(plot); invisible(dev.off())
  png(paste0(stem, ".png"), width = w, height = h, units = "in", res = dpi, type = "cairo")
  print(plot); invisible(dev.off())
  jpeg(paste0(stem, ".jpg"), width = w, height = h, units = "in", res = dpi, quality = 95, type = "cairo")
  print(plot); invisible(dev.off())
}

mirror_tsv <- function(df, stem_name) {
  write_tsv(df, file.path(out_fig, paste0(stem_name, ".tsv")))
  write_tsv(df, file.path(out_data, paste0(stem_name, ".tsv")))
}

## -------------------- [1] Load data --------------------
message("[1/5] Load groups, Plaque_Index, metabolite abundance ...")
group <- read_tsv(file.path(base_dir, "this.study.group.tsv"),
                  col_names = c("sample", "group"), show_col_types = FALSE) %>%
  mutate(group = factor(group, levels = c("Normal", "PoorPlaque")))
samples <- group$sample

env <- read_tsv(file.path(base_dir, "paper_result/env.tsv"), show_col_types = FALSE) %>%
  rename(variable = 1) %>%
  filter(variable == "Plaque_Index")
pi_vec <- setNames(as.numeric(unlist(env[1, samples])), samples)

metab <- read_tsv(file.path(base_dir, "output/tables/metab_abund_AB_select.tsv"),
                  show_col_types = FALSE)
metab <- as.data.frame(metab)
rownames(metab) <- metab$metab_id
metab$metab_id <- NULL
metab <- metab[, samples, drop = FALSE]
metab[] <- lapply(metab, as.numeric)

# drop zero-variance metabolites before scaling
row_sd <- apply(metab, 1, sd, na.rm = TRUE)
keep <- is.finite(row_sd) & row_sd > 0
metab <- metab[keep, , drop = FALSE]
message("  Samples: ", length(samples),
        " (Normal=", sum(group$group == "Normal"),
        ", PoorPlaque=", sum(group$group == "PoorPlaque"), ")")
message("  Metabolites used (non-zero variance): ", nrow(metab))

sel <- read_tsv(file.path(base_dir, "select.metab.tsv"), show_col_types = FALSE) %>%
  rename(metab_id = `Metab ID`) %>%
  mutate(
    Metabolite = ifelse(is.na(Metabolite) | Metabolite %in% c("-", ""), metab_id, Metabolite)
  ) %>%
  select(metab_id, Metabolite) %>%
  distinct(metab_id, .keep_all = TRUE)

## -------------------- [2] PCA + PERMANOVA --------------------
message("[2/5] Run PCA (scaled) + PERMANOVA (Euclidean) ...")
# samples × metabolites
X <- t(as.matrix(metab))
pca <- prcomp(X, center = TRUE, scale. = TRUE)

var_exp <- (pca$sdev^2) / sum(pca$sdev^2) * 100
n_pc <- min(N_PC_SCREE, length(var_exp))

coord_df <- tibble(
  sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = if (ncol(pca$x) >= 3) pca$x[, 3] else NA_real_,
  group = group$group[match(sample, group$sample)],
  Plaque_Index = pi_vec[sample]
)
write_tsv(coord_df, file.path(out_tab, "pca_coordinates.tsv"))
write_tsv(coord_df, file.path(out_data, "Fig_metab_PCA_score.tsv"))

var_df <- tibble(
  PC = paste0("PC", seq_along(var_exp)),
  PC_index = seq_along(var_exp),
  variance_pct = as.numeric(var_exp),
  cumulative_pct = cumsum(as.numeric(var_exp))
)
write_tsv(var_df, file.path(out_tab, "pca_variance_explained.tsv"))
write_tsv(var_df, file.path(out_data, "Fig_metab_PCA_scree.tsv"))

# Loadings (rotation): metabolites × PCs
load_mat <- pca$rotation
load_long <- as.data.frame(load_mat[, 1:min(5, ncol(load_mat)), drop = FALSE]) %>%
  rownames_to_column("metab_id") %>%
  left_join(sel, by = "metab_id") %>%
  mutate(Metabolite = coalesce(Metabolite, metab_id))
write_tsv(load_long, file.path(out_tab, "pca_loadings_PC1to5.tsv"))
write_tsv(load_long, file.path(out_data, "pca_loadings_PC1to5.tsv"))

# PERMANOVA on Euclidean distance of scaled matrix
X_scaled <- scale(X, center = TRUE, scale = TRUE)
dist_mat <- dist(X_scaled, method = "euclidean")
group_df <- data.frame(
  sample = rownames(X_scaled),
  group = factor(group$group[match(rownames(X_scaled), group$sample)],
                 levels = c("Normal", "PoorPlaque")),
  row.names = rownames(X_scaled)
)
perm_res <- adonis2(dist_mat ~ group, data = group_df, permutations = 999)
R2   <- as.numeric(perm_res$R2[1])
pval <- perm_res$`Pr(>F)`[1]
perm_df <- as.data.frame(perm_res) %>%
  rownames_to_column("Term")
write_tsv(perm_df, file.path(out_tab, "permanova_metab_euclidean.tsv"))
write_tsv(perm_df, file.path(out_data, "permanova_metab_euclidean.tsv"))

# PC1 ~ Plaque_Index Spearman
ct <- cor.test(coord_df$PC1, coord_df$Plaque_Index, method = "spearman", exact = FALSE)
pc1_pi <- tibble(
  metric = c("Spearman_rho", "P_value", "n"),
  value = c(unname(ct$estimate), ct$p.value, nrow(coord_df))
)
write_tsv(pc1_pi, file.path(out_tab, "PC1_vs_PlaqueIndex_spearman.tsv"))
write_tsv(pc1_pi, file.path(out_data, "PC1_vs_PlaqueIndex_spearman.tsv"))

pval_label <- ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))
stat_label <- sprintf("PERMANOVA: R\u00b2=%.3f, P=%s", R2, pval_label)
message("  ", stat_label)
message(sprintf("  PC1 vs Plaque_Index: rho=%.3f, P=%.4g",
                unname(ct$estimate), ct$p.value))

## -------------------- [3] Score plot --------------------
message("[3/5] Draw score / scree / loadings / PC1–PI panels ...")
var1 <- round(var_exp[1], 1)
var2 <- round(var_exp[2], 1)

p_score <- ggplot(coord_df, aes(x = PC1, y = PC2, fill = group, shape = group)) +
  stat_ellipse(aes(color = group), type = "t", level = 0.95,
               linewidth = 0.9, linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 4.2, color = "grey20", stroke = 0.5, alpha = 0.9) +
  scale_fill_manual(values = GROUP_COLORS, name = "Group") +
  scale_color_manual(values = GROUP_COLORS) +
  scale_shape_manual(values = GROUP_SHAPES, name = "Group") +
  annotate("text", x = Inf, y = Inf, label = stat_label,
           hjust = 1.05, vjust = 1.5, size = 3.8, color = "grey20",
           fontface = "italic") +
  labs(
    x = paste0("PC1 (", var1, "%)"),
    y = paste0("PC2 (", var2, "%)"),
    title = "PCA: overall metabolome profile",
    subtitle = paste0("Select metabolites (n=", nrow(metab),
                      "), samples A/B (n=", nrow(coord_df), "); z-scored")
  ) +
  theme_classic(base_size = BASE_PT) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border     = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
    axis.line        = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", size = BASE_PT + 1),
    plot.subtitle    = element_text(colour = "grey30", size = BASE_PT - 1),
    plot.margin      = margin(8, 12, 6, 6)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 4, shape = c(21, 24))),
         shape = "none")

stem_score <- file.path(out_fig, "Fig_metab_PCA_score")
save_pub(p_score, stem_score, width_mm = 170, height_mm = 145)
mirror_tsv(coord_df, "Fig_metab_PCA_score")

## -------------------- [4] Scree plot --------------------
scree_df <- var_df %>%
  filter(PC_index <= n_pc) %>%
  mutate(PC = factor(PC, levels = PC))

p_scree <- ggplot(scree_df, aes(x = PC, y = variance_pct)) +
  geom_col(fill = "#3C5488", width = 0.72, alpha = 0.9) +
  geom_line(aes(x = as.numeric(PC), y = cumulative_pct),
            color = "#E64B35", linewidth = 0.9) +
  geom_point(aes(x = as.numeric(PC), y = cumulative_pct),
             color = "#E64B35", size = 2.4) +
  scale_y_continuous(
    name = "Variance explained (%)",
    sec.axis = sec_axis(~ ., name = "Cumulative (%)")
  ) +
  labs(
    x = NULL,
    title = "PCA scree plot",
    subtitle = paste0("Top ", n_pc, " principal components")
  ) +
  theme_classic(base_size = BASE_PT) +
  theme(
    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
    axis.line = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey30", size = BASE_PT - 1),
    axis.title.y.right = element_text(color = "#E64B35"),
    axis.text.y.right = element_text(color = "#E64B35")
  )

stem_scree <- file.path(out_fig, "Fig_metab_PCA_scree")
save_pub(p_scree, stem_scree, width_mm = 150, height_mm = 120)
mirror_tsv(scree_df, "Fig_metab_PCA_scree")

## -------------------- [5] Top loadings (PC1 & PC2) --------------------
pick_top_load <- function(pc_col, n = N_LOAD_TOP) {
  load_long %>%
    transmute(
      metab_id, Metabolite,
      loading = .data[[pc_col]],
      abs_loading = abs(.data[[pc_col]]),
      PC = pc_col
    ) %>%
    arrange(desc(abs_loading)) %>%
    slice_head(n = n) %>%
    mutate(
      label = str_trunc(Metabolite, 42),
      label = factor(label, levels = rev(unique(label))),
      direction = ifelse(loading >= 0, "Positive", "Negative")
    )
}

load_plot_df <- bind_rows(pick_top_load("PC1"), pick_top_load("PC2"))
write_tsv(load_plot_df, file.path(out_tab, "Fig_metab_PCA_loadings.tsv"))
mirror_tsv(load_plot_df, "Fig_metab_PCA_loadings")

p_load <- ggplot(load_plot_df, aes(x = loading, y = label, fill = direction)) +
  geom_col(width = 0.72, alpha = 0.9) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
  facet_wrap(~ PC, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(Positive = "#E64B35", Negative = "#4DBBD5"),
                    name = "Loading") +
  labs(
    x = "PCA loading",
    y = NULL,
    title = "Top metabolite loadings",
    subtitle = paste0("Top ", N_LOAD_TOP, " |loading| metabolites per PC")
  ) +
  theme_classic(base_size = BASE_PT - 1) +
  theme(
    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.55),
    axis.line = element_blank(),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey30"),
    legend.position = "bottom"
  )

stem_load <- file.path(out_fig, "Fig_metab_PCA_loadings")
save_pub(p_load, stem_load, width_mm = 220, height_mm = 160)

## -------------------- PC1 vs Plaque Index --------------------
pi_lab <- sprintf("Spearman \u03c1=%.3f, P=%s",
                  unname(ct$estimate),
                  ifelse(ct$p.value < 0.001, "<0.001", sprintf("%.3f", ct$p.value)))

p_pi <- ggplot(coord_df, aes(x = Plaque_Index, y = PC1, fill = group, shape = group)) +
  geom_smooth(method = "lm", formula = y ~ x, color = "grey35",
              fill = "grey80", linewidth = 0.8, alpha = 0.35,
              inherit.aes = FALSE, aes(x = Plaque_Index, y = PC1)) +
  geom_point(size = 3.8, color = "grey20", stroke = 0.45, alpha = 0.9) +
  scale_fill_manual(values = GROUP_COLORS, name = "Group") +
  scale_shape_manual(values = GROUP_SHAPES, name = "Group") +
  annotate("text", x = Inf, y = Inf, label = pi_lab,
           hjust = 1.05, vjust = 1.5, size = 3.6, fontface = "italic") +
  labs(
    x = "Plaque Index",
    y = paste0("PC1 (", var1, "%)"),
    title = "PC1 versus Plaque Index",
    subtitle = "Overall metabolome principal component vs clinical plaque score"
  ) +
  theme_classic(base_size = BASE_PT) +
  theme(
    panel.border = element_rect(color = "grey30", fill = NA, linewidth = 0.6),
    axis.line = element_blank(),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey30", size = BASE_PT - 1),
    legend.position = "right"
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3.8, shape = c(21, 24))),
         shape = "none")

stem_pi <- file.path(out_fig, "Fig_metab_PCA_PC1_PI")
save_pub(p_pi, stem_pi, width_mm = 160, height_mm = 135)
mirror_tsv(coord_df %>% select(sample, group, Plaque_Index, PC1, PC2),
           "Fig_metab_PCA_PC1_PI")

## -------------------- Composite --------------------
message("[4/5] Assemble composite figure ...")
p_comp <- (p_score | p_scree) / (p_load | p_pi) +
  plot_annotation(
    tag_levels = "A",
    title = "Overall metabolome PCA",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.tag = element_text(face = "bold", size = 14)
    )
  ) +
  plot_layout(heights = c(1, 1.15))

stem_comp <- file.path(out_fig, "Fig_metab_PCA_composite")
save_pub(p_comp, stem_comp, width_mm = 280, height_mm = 240)
comp_tsv <- tibble(
  panel = c("A", "B", "C", "D"),
  figure = c("Fig_metab_PCA_score", "Fig_metab_PCA_scree",
             "Fig_metab_PCA_loadings", "Fig_metab_PCA_PC1_PI"),
  description = c(
    "PCA score plot (PC1–PC2) with 95% ellipse + PERMANOVA",
    "Scree + cumulative variance of top PCs",
    "Top |loading| metabolites on PC1 and PC2",
    "PC1 vs Plaque Index (Spearman)"
  ),
  n_metabolites = nrow(metab),
  n_samples = nrow(coord_df),
  permanova_R2 = R2,
  permanova_P = pval,
  PC1_var_pct = var_exp[1],
  PC2_var_pct = var_exp[2]
)
write_tsv(comp_tsv, paste0(stem_comp, ".tsv"))
file.copy(paste0(stem_comp, ".tsv"),
          file.path(out_data, "Fig_metab_PCA_composite.tsv"),
          overwrite = TRUE)

## -------------------- Summary --------------------
message("[5/5] Write analysis summary ...")
summary_df <- tibble(
  item = c(
    "n_samples", "n_Normal", "n_PoorPlaque", "n_metabolites",
    "PC1_variance_pct", "PC2_variance_pct", "PC1_PC2_cumulative_pct",
    "permanova_R2", "permanova_P",
    "PC1_PlaqueIndex_Spearman_rho", "PC1_PlaqueIndex_P"
  ),
  value = c(
    nrow(coord_df),
    sum(coord_df$group == "Normal"),
    sum(coord_df$group == "PoorPlaque"),
    nrow(metab),
    round(var_exp[1], 4),
    round(var_exp[2], 4),
    round(sum(var_exp[1:2]), 4),
    round(R2, 4),
    pval,
    round(unname(ct$estimate), 4),
    ct$p.value
  )
)
write_tsv(summary_df, file.path(out_tab, "pca_analysis_summary.tsv"))
write_tsv(summary_df, file.path(out_data, "pca_analysis_summary.tsv"))

message("[DONE] Module 14 metabolome PCA → ", module_dir)
message("[INFO] ", stat_label)
