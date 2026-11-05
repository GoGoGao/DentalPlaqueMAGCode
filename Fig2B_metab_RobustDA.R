## Packaged for revise/; source: project/01this_study_dataset/sgb/metab_corr/paper_result/05_RobustDA_metab_DA/scripts/05_RobustDA_metab_bar_volcano.R
## Figures: Fig2B
## Note: Robust_DA metabolite volcano/bar
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig2B_metab_RobustDA.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# Panel B SCI-ready: Robust_DA metabolite bar+scatter + volcano
# - drop metab_xxx from axis labels
# - color metabolite labels by select.metab.tsv col23 = HMDB Class
# - fonts 2× vs original panel_B (base ~7 → 14)
# - volcano of differential metabolites
# Output: paper_result/05_RobustDA_metab_DA/ + sync output/figures/panel_B/

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggtext)
  library(ggrepel)
  library(scales)
  library(svglite)
})

base_dir   <- file.path(.pkg_root, "01this_study_dataset/sgb/metab_corr"
module_dir <- file.path(base_dir, "paper_result/05_RobustDA_metab_DA")
out_fig    <- file.path(module_dir, "figures")
out_tab    <- file.path(module_dir, "tables")
out_data   <- file.path(module_dir, "data")
sync_fig   <- file.path(base_dir, "output/figures/panel_B")
dir.create(out_fig,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
dir.create(sync_fig, recursive = TRUE, showWarnings = FALSE)
source(file.path(base_dir, "paper_result/_shared/export_utils.R"))

COL    <- c(Normal = "#2E8B8B", PoorPlaque = "#B85450", NS = "#B0B0B0")
COL_PT <- c(Normal = "#A8D0D0", PoorPlaque = "#E0B4B2")
STRIP_BG <- c(Normal = "#D6EEEE", PoorPlaque = "#F6D9D7")
BASE_PT <- 14  # 2× original ~7
N_TOP_SIDE <- 20L
has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)

save_pub <- function(plot, stem, width_mm = 180, height_mm = 220, dpi = 600) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h); print(plot); invisible(dev.off())
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h); print(plot); invisible(dev.off())
  png(paste0(stem, ".png"), width = w, height = h, units = "in", res = dpi, type = "cairo")
  print(plot); invisible(dev.off())
  jpeg(paste0(stem, ".jpg"), width = w, height = h, units = "in", res = dpi, quality = 95, type = "cairo")
  print(plot); invisible(dev.off())
}

# Distinct palette for HMDB Class
build_class_palette <- function(classes) {
  classes <- unique(as.character(classes))
  classes <- classes[!is.na(classes) & nzchar(classes)]
  n <- length(classes)
  # high-contrast multi-hue palette
  base <- c(
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
    "#E18727", "#20854E", "#7876B1", "#6F99AD", "#FFDC91",
    "#EE4C97", "#0072B5", "#BC3C29", "#E18727", "#20854E",
    "#7876B1", "#6F99AD", "#FFDC91", "#EE4C97", "#631879",
    "#008280", "#BB0021", "#5F559B", "#A20056", "#808180"
  )
  if (n <= length(base)) {
    pal <- base[seq_len(n)]
  } else {
    pal <- scales::hue_pal(h = c(0, 360) + 15, c = 90, l = 50)(n)
  }
  names(pal) <- classes
  pal
}

normalize_class <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x) | x %in% c("-", "NA", "null")] <- "Unclassified"
  x
}

message("[1/4] Load groups, abundance, RobustDA metabolite table ...")
group <- read_tsv(file.path(base_dir, "this.study.group.tsv"),
                  col_names = c("sample", "group"), show_col_types = FALSE) %>%
  mutate(group = factor(group, levels = c("Normal", "PoorPlaque")))
samples <- group$sample

metab <- read_tsv(file.path(base_dir, "output/tables/metab_abund_AB_select.tsv"),
                  show_col_types = FALSE)
metab <- as.data.frame(metab)
rownames(metab) <- metab$metab_id
metab$metab_id <- NULL
metab <- metab[, samples, drop = FALSE]
metab[] <- lapply(metab, as.numeric)

# select.metab.tsv col23 = HMDB Class
sel_raw <- read_tsv(file.path(base_dir, "select.metab.tsv"), show_col_types = FALSE)
hmdb_col <- names(sel_raw)[23]
message("  Classification column (col23): ", hmdb_col)

sel <- sel_raw %>%
  rename(metab_id = `Metab ID`) %>%
  mutate(
    Metabolite = ifelse(is.na(Metabolite) | Metabolite %in% c("-", ""), metab_id, Metabolite),
    HMDB_Class = normalize_class(.data[[hmdb_col]])
  )

da_file <- file.path(base_dir, "output/tables/metabolite_RobustDA.tsv")
if (!file.exists(da_file)) stop("Missing ", da_file, " — run scripts/04_RobustDA_network.R first")
metab_da <- read_tsv(da_file, show_col_types = FALSE) %>%
  left_join(sel %>% select(metab_id, HMDB_Class, Metabolite_sel = Metabolite), by = "metab_id") %>%
  mutate(
    Metabolite = coalesce(Metabolite, Metabolite_sel, metab_id),
    HMDB_Class = normalize_class(HMDB_Class),
    # effect size for volcano: signed log2 fold-change of means
    log2FC = log2((mean_PoorPlaque + 1e-6) / (mean_Normal + 1e-6)),
    neglog10_FDR = -log10(pmax(Wilcoxon_FDR, 1e-300))
  )

metab_rob <- metab_da %>% filter(Robust_DA)
message("  Robust_DA metabolites: ", nrow(metab_rob),
        " (Normal=", sum(metab_rob$Robust_direction == "Normal"),
        ", PoorPlaque=", sum(metab_rob$Robust_direction == "PoorPlaque"), ")")

write_tsv(
  metab_da %>%
    select(metab_id, Metabolite, HMDB_Class, mean_Normal, mean_PoorPlaque, delta, log2FC,
           Wilcoxon_p, Wilcoxon_FDR, CLRMC_p, CLRMC_FDR, Robust_DA, Robust_direction),
  file.path(out_tab, "metabolite_RobustDA_with_HMDB_Class.tsv")
)
write_tsv(
  metab_rob %>%
    select(metab_id, Metabolite, HMDB_Class, mean_Normal, mean_PoorPlaque, delta, log2FC,
           Wilcoxon_FDR, CLRMC_FDR, Robust_direction) %>%
    arrange(Robust_direction, desc(abs(delta))),
  file.path(out_tab, "metabolite_RobustDA_only_with_HMDB_Class.tsv")
)

## -------------------- Panel B bar + scatter --------------------
message("[2/4] Bar + scatter (top ", N_TOP_SIDE, " / side by |delta|) ...")
top_metab <- metab_rob %>%
  arrange(Robust_direction, desc(abs(delta))) %>%
  group_by(Robust_direction) %>%
  slice_head(n = N_TOP_SIDE) %>%
  ungroup() %>%
  mutate(
    # axis label: metabolite name only (no metab_xxx)
    lab_plain = make.unique(str_trunc(Metabolite, 40, side = "right"), sep = " ·")
  )

class_levels <- sort(unique(top_metab$HMDB_Class))
CLASS_PAL <- build_class_palette(class_levels)
# markdown y labels colored by HMDB Class
top_metab <- top_metab %>%
  mutate(
    lab_md = sprintf(
      "<span style='color:%s'>%s</span>",
      CLASS_PAL[HMDB_Class],
      str_replace_all(lab_plain, c("&" = "&amp;", "<" = "&lt;", ">" = "&gt;"))
    )
  )

metab_long <- metab[top_metab$metab_id, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("metab_id") %>%
  pivot_longer(-metab_id, names_to = "sample", values_to = "intensity") %>%
  left_join(group, by = "sample") %>%
  left_join(
    top_metab %>% select(metab_id, Metabolite, lab_plain, lab_md, Robust_direction,
                         delta, HMDB_Class),
    by = "metab_id"
  )

metab_sum <- metab_long %>%
  group_by(metab_id, Metabolite, lab_plain, lab_md, Robust_direction, group,
           delta, HMDB_Class) %>%
  summarise(mean = mean(intensity), se = sd(intensity) / sqrt(n()), .groups = "drop")

ord_labs <- top_metab %>%
  arrange(Robust_direction, desc(abs(delta))) %>%
  pull(lab_md)
metab_long$lab_md <- factor(metab_long$lab_md, levels = rev(ord_labs))
metab_sum$lab_md  <- factor(metab_sum$lab_md,  levels = rev(ord_labs))
metab_sum$Robust_direction <- factor(metab_sum$Robust_direction, levels = c("Normal", "PoorPlaque"))
metab_long$Robust_direction <- factor(metab_long$Robust_direction, levels = c("Normal", "PoorPlaque"))
metab_sum$HMDB_Class <- factor(metab_sum$HMDB_Class, levels = class_levels)

write_tsv(metab_long, file.path(out_data, "Fig_panelB_metab_points.tsv"))
write_tsv(metab_sum,  file.path(out_data, "Fig_panelB_metab_summary.tsv"))
write_tsv(
  top_metab %>%
    select(metab_id, Metabolite, lab_plain, HMDB_Class, Robust_direction, delta, log2FC,
           Wilcoxon_FDR, CLRMC_FDR),
  file.path(out_tab, "Fig_panelB_metab_list.tsv")
)
write_tsv(
  tibble(HMDB_Class = names(CLASS_PAL), fill_hex = unname(CLASS_PAL)),
  file.path(out_tab, "Fig_panelB_HMDB_Class_palette.tsv")
)

# left-side class ticks (one per metabolite) for HMDB Class legend
class_tick <- metab_sum %>%
  distinct(lab_md, HMDB_Class) %>%
  mutate(HMDB_Class = factor(HMDB_Class, levels = class_levels))

use_ggnewscale <- requireNamespace("ggnewscale", quietly = TRUE)

pB <- ggplot() +
  geom_col(
    data = metab_sum, aes(x = lab_md, y = mean, fill = group),
    position = position_dodge(width = 0.75), width = 0.7, alpha = 0.85
  ) +
  geom_errorbar(
    data = metab_sum,
    aes(x = lab_md, ymin = mean - se, ymax = mean + se, group = group),
    position = position_dodge(width = 0.75), width = 0.25,
    linewidth = 0.45, color = "grey30"
  ) +
  geom_point(
    data = metab_long, aes(x = lab_md, y = intensity, color = group),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
    size = 1.4, alpha = 0.75, stroke = 0
  ) +
  scale_fill_manual(values = COL[c("Normal", "PoorPlaque")], name = "Group") +
  scale_color_manual(values = COL_PT, guide = "none")

if (use_ggnewscale) {
  pB <- pB +
    ggnewscale::new_scale_fill() +
    geom_point(
      data = class_tick, aes(x = lab_md, y = 0, fill = HMDB_Class),
      shape = 22, size = 2.8, color = "grey20", stroke = 0.2
    ) +
    scale_fill_manual(values = CLASS_PAL, name = "HMDB Class")
} else {
  # fallback: class color only via y-axis markdown (no separate class legend layer)
  pB <- pB
}

pB <- pB +
  coord_flip(clip = "off") +
  labs(
    title = "Robust_DA metabolites (Wilcoxon ∩ CLR-MC, FDR<0.05)",
    subtitle = paste0("Top ", nrow(top_metab), " / ", nrow(metab_rob),
                      " by |Δ|; metabolite names colored by HMDB Class (select.metab col23)"),
    x = NULL, y = "Intensity"
  )

if (has_ggh4x) {
  pB <- pB + ggh4x::facet_grid2(
    Robust_direction ~ ., scales = "free_y", space = "free_y",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(fill = unname(STRIP_BG[c("Normal", "PoorPlaque")])),
      text_y = ggh4x::elem_list_text(
        colour = unname(COL[c("Normal", "PoorPlaque")]),
        face = "bold", size = BASE_PT - 1
      )
    )
  )
} else {
  pB <- pB +
    facet_grid(Robust_direction ~ ., scales = "free_y", space = "free_y") +
    theme(strip.text.y = element_text(size = BASE_PT - 1, face = "bold"))
}

pB <- pB +
  theme_classic(base_size = BASE_PT) +
  theme(
    axis.text.y = element_markdown(size = BASE_PT - 4, lineheight = 1.05),
    axis.text.x = element_text(size = BASE_PT - 2),
    axis.title  = element_text(size = BASE_PT),
    legend.position = "bottom",
    legend.title = element_text(size = BASE_PT - 2, face = "bold"),
    legend.text  = element_text(size = BASE_PT - 4),
    legend.key.size = unit(0.35, "cm"),
    plot.title = element_text(size = BASE_PT + 2, face = "bold"),
    plot.subtitle = element_text(size = BASE_PT - 2, colour = "grey30"),
    strip.background = element_rect(fill = "grey95", colour = NA),
    plot.margin = margin(8, 12, 8, 8)
  )

stem_b <- file.path(out_fig, "Fig_panelB_RobustDA_metab_bar_scatter")
message("[3/4] Export bar+scatter → ", stem_b)
save_pub(pB, stem_b, width_mm = 200, height_mm = 280)

## -------------------- Volcano --------------------
message("[4/4] Volcano plot ...")
vol_df <- metab_da %>%
  mutate(
    sig_group = factor(
      case_when(
        Robust_DA & Robust_direction == "Normal" ~ "Normal",
        Robust_DA & Robust_direction == "PoorPlaque" ~ "PoorPlaque",
        TRUE ~ "NS"
      ),
      levels = c("PoorPlaque", "Normal", "NS")
    )
  )

# Balanced top labels per side (score = |log2FC| × −log10(FDR)); keep short
N_LAB_PER_SIDE <- 5L
lab_tbl <- metab_rob %>%
  mutate(
    score = abs(log2FC) * neglog10_FDR,
    label_show = str_trunc(Metabolite, 20, side = "right")
  ) %>%
  group_by(Robust_direction) %>%
  arrange(desc(score), Wilcoxon_FDR) %>%
  slice_head(n = N_LAB_PER_SIDE) %>%
  ungroup() %>%
  select(metab_id, label_show, Robust_direction, log2FC, neglog10_FDR, score)

vol_df <- vol_df %>%
  left_join(lab_tbl %>% select(metab_id, label_show), by = "metab_id")

vol_lab <- vol_df %>%
  filter(!is.na(label_show)) %>%
  mutate(
    # push Normal (down) labels left, PoorPlaque (up) labels right
    nudge_x = ifelse(log2FC < 0, -0.55, 0.55)
  )

write_tsv(
  vol_df %>%
    select(metab_id, Metabolite, HMDB_Class, log2FC, delta, Wilcoxon_FDR, CLRMC_FDR,
           neglog10_FDR, Robust_DA, Robust_direction, sig_group, label_show),
  file.path(out_data, "Fig_panelB_volcano_data.tsv")
)
write_tsv(lab_tbl, file.path(out_tab, "Fig_panelB_volcano_labeled_metabolites.tsv"))

x_pad <- max(0.6, 0.12 * diff(range(vol_df$log2FC, na.rm = TRUE)))
y_max <- max(vol_df$neglog10_FDR, na.rm = TRUE) * 1.08

# Volcano colored by HMDB Class (Robust only)
pV <- ggplot(vol_df, aes(x = log2FC, y = neglog10_FDR)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_point(data = filter(vol_df, sig_group == "NS"),
             color = COL[["NS"]], size = 1.6, alpha = 0.45) +
  geom_point(data = filter(vol_df, sig_group != "NS"),
             aes(color = HMDB_Class), size = 2.4, alpha = 0.9) +
  scale_color_manual(values = build_class_palette(unique(metab_rob$HMDB_Class)),
                     name = "HMDB Class (Robust_DA)") +
  ggrepel::geom_label_repel(
    data = vol_lab,
    aes(x = log2FC, y = neglog10_FDR, label = label_show),
    inherit.aes = FALSE,
    color = "grey15",
    fill = alpha("white", 0.92),
    label.size = 0.15,
    label.padding = unit(0.12, "lines"),
    size = 3.0,
    max.overlaps = Inf,
    force = 2.0, force_pull = 0.35,
    box.padding = 0.5, point.padding = 0.4,
    min.segment.length = 0.02,
    segment.size = 0.35, segment.color = "grey45", segment.alpha = 0.85,
    nudge_x = vol_lab$nudge_x,
    xlim = c(min(vol_df$log2FC) - x_pad, max(vol_df$log2FC) + x_pad),
    ylim = c(0, y_max),
    seed = 42, max.time = 2.5, max.iter = 30000,
    show.legend = FALSE
  ) +
  coord_cartesian(xlim = range(vol_df$log2FC) + c(-x_pad * 0.35, x_pad * 0.35),
                  ylim = c(0, y_max), clip = "off") +
  labs(
    title = "Volcano: metabolite differential abundance",
    subtitle = "x = log2(mean_PoorPlaque / mean_Normal); y = −log10(Wilcoxon FDR); colored = Robust_DA (Wilcoxon∩CLR-MC)",
    x = expression(log[2]~FC~(PoorPlaque / Normal)),
    y = expression(-log[10](Wilcoxon~FDR))
  ) +
  theme_classic(base_size = BASE_PT) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = BASE_PT - 3, face = "bold"),
    legend.text  = element_text(size = BASE_PT - 5),
    legend.key.size = unit(0.32, "cm"),
    plot.title = element_text(size = BASE_PT + 2, face = "bold"),
    plot.subtitle = element_text(size = BASE_PT - 3, colour = "grey30"),
    axis.title = element_text(size = BASE_PT),
    axis.text  = element_text(size = BASE_PT - 2),
    plot.margin = margin(8, 14, 8, 14)
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

stem_v <- file.path(out_fig, "Fig_panelB_RobustDA_metab_volcano")
save_pub(pV, stem_v, width_mm = 240, height_mm = 170)

# Direction counts for on-plot annotation (color = enrichment group)
n_Normal_vol <- sum(metab_rob$Robust_direction == "Normal")
n_PP_vol     <- sum(metab_rob$Robust_direction == "PoorPlaque")
x_lim_v2 <- range(vol_df$log2FC, na.rm = TRUE) + c(-x_pad * 0.45, x_pad * 0.45)
count_ann <- tibble(
  side = c("Normal", "PoorPlaque"),
  n = c(n_Normal_vol, n_PP_vol),
  label = c(
    sprintf("Normal n = %d", n_Normal_vol),
    sprintf("PoorPlaque n = %d", n_PP_vol)
  ),
  # left = Normal (down in PP / up in Normal); right = PoorPlaque
  x = c(x_lim_v2[1] + 0.04 * diff(x_lim_v2),
        x_lim_v2[2] - 0.04 * diff(x_lim_v2)),
  y = rep(y_max * 0.98, 2),
  hjust = c(0, 1),
  color = unname(COL[c("Normal", "PoorPlaque")])
)
write_tsv(count_ann, file.path(out_data, "Fig_panelB_volcano_direction_counts.tsv"))
write_tsv(count_ann, file.path(out_tab, "Fig_panelB_volcano_direction_counts.tsv"))

# Direction-colored volcano: white label boxes, outward nudge, more canvas
pV2 <- ggplot(vol_df, aes(x = log2FC, y = neglog10_FDR)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_point(aes(color = sig_group, size = sig_group, alpha = sig_group)) +
  scale_color_manual(values = COL, name = "Robust_DA") +
  scale_size_manual(values = c(Normal = 2.6, PoorPlaque = 2.6, NS = 1.4), guide = "none") +
  scale_alpha_manual(values = c(Normal = 0.9, PoorPlaque = 0.9, NS = 0.35), guide = "none") +
  ggrepel::geom_label_repel(
    data = vol_lab,
    aes(x = log2FC, y = neglog10_FDR, label = label_show, color = sig_group),
    inherit.aes = FALSE,
    fill = alpha("white", 0.93),
    label.size = 0.18,
    label.padding = unit(0.15, "lines"),
    size = 3.1,
    max.overlaps = Inf,
    force = 2.4, force_pull = 0.3,
    box.padding = 0.6, point.padding = 0.45,
    min.segment.length = 0.01,
    segment.size = 0.4, segment.color = "grey40", segment.alpha = 0.9,
    nudge_x = vol_lab$nudge_x,
    direction = "both",
    xlim = c(min(vol_df$log2FC) - x_pad, max(vol_df$log2FC) + x_pad),
    ylim = c(0, y_max),
    seed = 42, max.time = 3, max.iter = 40000,
    show.legend = FALSE
  ) +
  # up/down counts: color matches enrichment group
  geom_text(
    data = count_ann,
    aes(x = x, y = y, label = label, hjust = hjust, colour = side),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 4.2,
    show.legend = FALSE
  ) +
  coord_cartesian(xlim = x_lim_v2, ylim = c(0, y_max), clip = "off") +
  labs(
    title = "Volcano: Robust_DA metabolites",
    subtitle = sprintf(
      "Robust_DA n=%d; labels = top %d / side by |log2FC|×−log10(FDR)",
      nrow(metab_rob),
      N_LAB_PER_SIDE
    ),
    x = expression(log[2]~FC~(PoorPlaque / Normal)),
    y = expression(-log[10](Wilcoxon~FDR))
  ) +
  theme_classic(base_size = BASE_PT) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = BASE_PT - 1, face = "bold"),
    legend.text  = element_text(size = BASE_PT - 2),
    plot.title = element_text(size = BASE_PT + 2, face = "bold"),
    plot.subtitle = element_text(size = BASE_PT - 2, colour = "grey30"),
    axis.title = element_text(size = BASE_PT),
    axis.text  = element_text(size = BASE_PT - 2),
    plot.margin = margin(12, 22, 12, 22)
  )

stem_v2 <- file.path(out_fig, "Fig_panelB_RobustDA_metab_volcano_direction")
save_pub(pV2, stem_v2, width_mm = 205, height_mm = 170)
# also sync count TSV next to figure outputs
file.copy(file.path(out_data, "Fig_panelB_volcano_direction_counts.tsv"),
          file.path(out_fig, "Fig_panelB_volcano_direction_counts.tsv"), overwrite = TRUE)

# sync bar figure + TSVs to legacy panel_B
for (ext in c("pdf", "png", "jpg", "svg")) {
  file.copy(paste0(stem_b, ".", ext),
            file.path(sync_fig, paste0("Fig_panelB_RobustDA_metab_bar_scatter.", ext)),
            overwrite = TRUE)
}
file.copy(file.path(out_data, "Fig_panelB_metab_points.tsv"),
          file.path(sync_fig, "Fig_panelB_metab_points.tsv"), overwrite = TRUE)
file.copy(file.path(out_data, "Fig_panelB_metab_summary.tsv"),
          file.path(sync_fig, "Fig_panelB_metab_summary.tsv"), overwrite = TRUE)

message("[DONE] Panel B SCI figures → ", module_dir)
message("  Bar+scatter: ", stem_b)
message("  Volcano (HMDB Class): ", stem_v)
message("  Volcano (direction):  ", stem_v2)
