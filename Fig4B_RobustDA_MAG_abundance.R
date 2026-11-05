## Packaged for revise/; source: project/01this_study_dataset/sgb/metab_corr/paper_result/01_RobustDA_MAG_abundance/scripts/Fig_panelA_RobustDA_bar_scatter.R
## Figures: Fig4B
## Note: Robust_DA MAG abundance
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig4B_RobustDA_MAG_abundance.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# Fig panel_A: Robust_DA MAG abundance (bar = mean±SE + sample jitter)
# SCI-ready: 2× fonts; species Latin italic; MAG ID upright (ggtext markdown)
# Output: paper_result/01_RobustDA_MAG_abundance/{figures,tables,data}
#         + sync to metab_corr/output/figures/panel_A/

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggtext)
  library(scales)
  library(svglite)
})

base_dir   <- file.path(.pkg_root, "01this_study_dataset/sgb/metab_corr"
study_dir  <- dirname(dirname(base_dir))  # 01this_study_dataset
sig_dir    <- file.path(dirname(base_dir), "sig_mag/sgb_DA_results/tables")
module_dir <- file.path(base_dir, "paper_result/01_RobustDA_MAG_abundance")
out_fig    <- file.path(module_dir, "figures")
out_tab    <- file.path(module_dir, "tables")
out_data   <- file.path(module_dir, "data")
sync_fig   <- file.path(base_dir, "output/figures/panel_A")
# MetaPhlAn Robust_DA species (relative: ../../../M1-metaphalan/... from paper_result)
mp_csv <- file.path(
  study_dir, "M1-metaphalan/v2/Figure1/tables/Table_S_differential_species.csv"
)
dir.create(out_fig,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_tab,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
dir.create(sync_fig, recursive = TRUE, showWarnings = FALSE)

COL    <- c(Normal = "#2E8B8B", PoorPlaque = "#B85450")
COL_PT <- c(Normal = "#A8D0D0", PoorPlaque = "#E0B4B2")
STRIP_BG <- c(Normal = "#D6EEEE", PoorPlaque = "#F6D9D7")
# Highlight label color follows enrichment group color
HL_COLOR <- COL  # Normal / PoorPlaque

# GTDB/MetaPhlAn name harmonization → genus + species epithet
binomial_core <- function(x) {
  s <- str_trim(str_replace_all(as.character(x), "_", " "))
  toks <- str_split(s, "\\s+", simplify = TRUE)
  toks <- toks[nzchar(toks)]
  if (length(toks) >= 2 && grepl("^[A-Z]$", toks[2])) {
    toks <- c(toks[1], toks[-(1:2)])
  }
  if (length(toks) >= 3 && grepl("^[A-Z]$", toks[length(toks)])) {
    toks <- toks[-length(toks)]
  }
  genus <- tolower(toks[1])
  epithet <- if (length(toks) >= 2) tolower(toks[2]) else ""
  paste(genus, epithet)
}

# ---- fonts: original ~7 pt base → 2× ----
BASE_PT <- 14
has_ggh4x <- requireNamespace("ggh4x", quietly = TRUE)

theme_set(
  theme_classic(base_size = BASE_PT, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.5, colour = "black"),
      axis.ticks = element_line(linewidth = 0.5, colour = "black"),
      axis.title = element_text(size = BASE_PT),
      axis.text  = element_text(size = BASE_PT - 2, colour = "black"),
      legend.title = element_text(size = BASE_PT - 1),
      legend.text  = element_text(size = BASE_PT - 2),
      strip.text   = element_text(size = BASE_PT - 1, face = "bold"),
      plot.title   = element_text(size = BASE_PT + 2, face = "bold"),
      plot.subtitle = element_text(size = BASE_PT - 1, colour = "grey30"),
      panel.grid = element_blank()
    )
)

save_pub <- function(plot, stem, width_mm = 180, height_mm = 280, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  print(plot); invisible(dev.off())
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h)
  print(plot); invisible(dev.off())
  png(paste0(stem, ".png"), width = w, height = h, units = "in", res = dpi, type = "cairo")
  print(plot); invisible(dev.off())
  jpeg(paste0(stem, ".jpg"), width = w, height = h, units = "in", res = dpi, quality = 95)
  print(plot); invisible(dev.off())
}

# Species Latin italic via markdown; MAG ID upright; Unknown not italic
# highlight=TRUE → color by enrichment group (Normal teal / PoorPlaque red-brown)
make_y_lab <- function(species, mag, highlight = FALSE, direction = NA_character_) {
  sp <- str_replace_all(species %||% "Unknown", "_", " ")
  sp <- str_trim(sp)
  if (is.na(sp) || !nzchar(sp) || tolower(sp) == "unknown") {
    lab <- paste0("Unknown (", mag, ")")
  } else {
    sp_esc <- str_replace_all(sp, "([*_])", "\\\\\\1")
    lab <- paste0("*", sp_esc, "* (", mag, ")")
  }
  if (isTRUE(highlight)) {
    col <- HL_COLOR[[as.character(direction)]]
    if (is.null(col) || is.na(col)) col <- COL[["PoorPlaque"]]
    paste0("<span style='color:", col, "'>", lab, "</span>")
  } else {
    lab
  }
}
`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a[1])) a else b

message("[1/3] Load inputs ...")
group <- read_tsv(file.path(base_dir, "this.study.group.tsv"),
                  col_names = c("sample", "group"), show_col_types = FALSE) %>%
  mutate(group = factor(group, levels = c("Normal", "PoorPlaque")))
samples <- group$sample

mag_ab <- read_tsv(file.path(base_dir, "bin_abundance_table.tsv"), show_col_types = FALSE) %>%
  rename(MAG = `Genomic bins`)
mag_ab <- as.data.frame(mag_ab)
rownames(mag_ab) <- mag_ab$MAG
mag_ab$MAG <- NULL
mag_ab <- mag_ab[, samples, drop = FALSE]
mag_ab[] <- lapply(mag_ab, as.numeric)

da <- read_tsv(file.path(sig_dir, "Table_SGB_differential_bins.tsv"), show_col_types = FALSE) %>%
  rename(MAG = SGB) %>%
  mutate(Robust_DA = as.logical(Robust_DA))

anno_path <- file.path(base_dir, "output/tables/RobustDA_MAG_annotation.tsv")
if (file.exists(anno_path)) {
  mag_anno <- read_tsv(anno_path, show_col_types = FALSE)
} else {
  gtdb <- read_tsv(file.path(base_dir, "gtdbtk.bac120.summary.tsv"), show_col_types = FALSE) %>%
    transmute(
      MAG = user_genome,
      species_gtdb = str_extract(classification, "s__[^;]*") %>% str_remove("^s__")
    )
  mag_anno <- da %>%
    filter(Robust_DA) %>%
    left_join(gtdb, by = "MAG") %>%
    mutate(
      species = coalesce(species, species_gtdb, "Unknown"),
      focus_direction = Robust_direction,
      short_label = paste0(str_replace_all(species, "_", " "), " (", MAG, ")")
    )
}

robust_ids <- da$MAG[da$Robust_DA %in% TRUE]

# MetaPhlAn Robust_DA species (direction = Normal / PoorPlaque)
stopifnot(file.exists(mp_csv))
mp_robust <- read_csv(mp_csv, show_col_types = FALSE) %>%
  mutate(
    Robust_DA = as.logical(Robust_DA),
    species_plain_mp = str_replace_all(species, "_", " "),
    core = vapply(species, binomial_core, character(1))
  ) %>%
  filter(Robust_DA %in% TRUE, Robust_direction %in% c("Normal", "PoorPlaque"))

message("  MetaPhlAn Robust_DA species: ", nrow(mp_robust),
        " (Normal=", sum(mp_robust$Robust_direction == "Normal"),
        ", PoorPlaque=", sum(mp_robust$Robust_direction == "PoorPlaque"), ")")

mag_anno <- mag_anno %>%
  filter(MAG %in% rownames(mag_ab), MAG %in% robust_ids) %>%
  mutate(
    species_plain = str_replace_all(coalesce(species, "Unknown"), "_", " "),
    focus_direction = coalesce(focus_direction, Robust_direction),
    core = vapply(species_plain, binomial_core, character(1))
  )

# Match MAG ↔ MetaPhlAn by binomial core + same Robust_direction
mp_by_dir <- split(mp_robust, mp_robust$Robust_direction)
match_one <- function(core, direction, mag_plain) {
  pool <- mp_by_dir[[as.character(direction)]]
  empty <- list(hit = FALSE, mp_species = NA_character_, mp_log2FC = NA_real_,
                match_how = NA_character_)
  if (is.null(pool) || is.na(core) || !str_detect(core, "\\S+\\s+\\S+")) return(empty)
  idx <- which(pool$core == core)
  if (!length(idx)) return(empty)
  i <- idx[[1]]
  mp_plain <- str_replace_all(pool$species[[i]], "_", " ")
  how <- if (identical(tolower(mag_plain), tolower(mp_plain))) "exact" else "binomial_core"
  list(
    hit = TRUE,
    mp_species = pool$species[[i]],
    mp_log2FC = as.numeric(pool$log2FC[[i]]),
    match_how = how
  )
}

hits <- Map(match_one, mag_anno$core, mag_anno$focus_direction, mag_anno$species_plain)
mag_anno <- mag_anno %>%
  mutate(
    is_highlight = vapply(hits, `[[`, logical(1), "hit"),
    mp_species = vapply(hits, `[[`, character(1), "mp_species"),
    mp_log2FC = vapply(hits, `[[`, numeric(1), "mp_log2FC"),
    match_how = vapply(hits, `[[`, character(1), "match_how"),
    note_zh = ifelse(
      is_highlight,
      paste0("MetaPhlAn Robust_DA same direction: ",
             str_replace_all(mp_species, "_", " ")),
      NA_character_
    ),
    y_lab = mapply(make_y_lab, species_plain, MAG, is_highlight,
                   focus_direction, USE.NAMES = FALSE),
    short_label = paste0(species_plain, " (", MAG, ")")
  )

message("  Highlight (MetaPhlAn∩same direction): ",
        sum(mag_anno$is_highlight), " MAGs → ",
        paste(unique(mag_anno$species_plain[mag_anno$is_highlight]), collapse = "; "))

# overlap lookup table
write_tsv(
  mag_anno %>%
    filter(is_highlight) %>%
    select(MAG, species_plain, focus_direction, log2FC,
           mp_species, mp_log2FC, match_how, note_zh, y_lab) %>%
    arrange(focus_direction, desc(abs(log2FC))),
  file.path(out_tab, "Fig_panelA_MetaPhlAn_overlap.tsv")
)
# also keep unmatched MetaPhlAn Robust_DA for transparency
matched_mp <- unique(na.omit(mag_anno$mp_species))
write_tsv(
  mp_robust %>%
    mutate(matched_in_PanelA = species %in% matched_mp) %>%
    select(species, species_plain_mp, log2FC, Robust_direction,
           Wilcoxon_FDR, CLRMC_FDR, matched_in_PanelA, core) %>%
    arrange(Robust_direction, desc(abs(log2FC))),
  file.path(out_tab, "Fig_panelA_MetaPhlAn_RobustDA_species.tsv")
)

focus_mags <- mag_anno$MAG
message("  Robust_DA MAGs: ", length(focus_mags),
        " (Normal=", sum(mag_anno$focus_direction == "Normal"),
        ", PoorPlaque=", sum(mag_anno$focus_direction == "PoorPlaque"), ")")

message("[2/3] Summarise & plot ...")
mag_long <- mag_ab[focus_mags, , drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("MAG") %>%
  pivot_longer(-MAG, names_to = "sample", values_to = "abund") %>%
  left_join(group, by = "sample") %>%
  left_join(
    mag_anno %>% select(MAG, short_label, y_lab, focus_direction, log2FC, species_plain),
    by = "MAG"
  )

mag_sum <- mag_long %>%
  group_by(MAG, short_label, y_lab, focus_direction, group, log2FC, species_plain) %>%
  summarise(
    mean = mean(abund, na.rm = TRUE),
    se = sd(abund, na.rm = TRUE) / sqrt(n()),
    sd = sd(abund, na.rm = TRUE),
    .groups = "drop"
  )

ord_labs <- mag_anno %>%
  arrange(focus_direction, desc(abs(log2FC))) %>%
  pull(y_lab)
mag_long$y_lab <- factor(mag_long$y_lab, levels = rev(ord_labs))
mag_sum$y_lab  <- factor(mag_sum$y_lab,  levels = rev(ord_labs))

# plotting TSVs (plain short_label for readability + markdown y_lab)
write_tsv(mag_long, file.path(out_data, "Fig_panelA_data_points.tsv"))
write_tsv(mag_sum,  file.path(out_data, "Fig_panelA_data_summary.tsv"))
keep_cols <- c("MAG", "species_plain", "short_label", "y_lab", "focus_direction",
               "is_highlight", "note_zh", "mp_species", "mp_log2FC", "match_how",
               "log2FC", "Wilcoxon_FDR", "CLRMC_FDR", "Robust_DA", "Robust_direction")
write_tsv(
  mag_anno %>%
    select(any_of(keep_cols)) %>%
    arrange(focus_direction, desc(abs(log2FC))),
  file.path(out_tab, "Fig_panelA_MAG_list.tsv")
)
write_tsv(
  mag_anno %>%
    filter(is_highlight) %>%
    select(MAG, species_plain, note_zh, focus_direction, log2FC,
           mp_species, mp_log2FC, match_how, y_lab) %>%
    arrange(desc(abs(log2FC))),
  file.path(out_tab, "Fig_panelA_highlight_taxa.tsv")
)

# Facet labels colored by group (markdown fallback; ggh4x strip colors if available)
strip_labs <- c(
  Normal = paste0("<span style='color:", COL[["Normal"]], "'><b>Normal</b></span>"),
  PoorPlaque = paste0("<span style='color:", COL[["PoorPlaque"]], "'><b>PoorPlaque</b></span>")
)
mag_sum$focus_direction <- factor(mag_sum$focus_direction, levels = c("Normal", "PoorPlaque"))
mag_long$focus_direction <- factor(mag_long$focus_direction, levels = c("Normal", "PoorPlaque"))

pA <- ggplot() +
  geom_col(
    data = mag_sum, aes(x = y_lab, y = mean, fill = group),
    position = position_dodge(width = 0.75), width = 0.7, alpha = 0.85
  ) +
  geom_errorbar(
    data = mag_sum,
    aes(x = y_lab, ymin = pmax(mean - se, 0), ymax = mean + se, group = group),
    position = position_dodge(width = 0.75), width = 0.25,
    linewidth = 0.45, color = "grey30"
  ) +
  geom_point(
    data = mag_long, aes(x = y_lab, y = abund, color = group),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
    size = 1.4, alpha = 0.75, stroke = 0
  ) +
  scale_fill_manual(values = COL, name = "Group") +
  scale_color_manual(values = COL_PT, name = "Group") +
  coord_flip() +
  labs(
    title = "Robust_DA MAG abundance (Wilcoxon ∩ CLR-MC)",
    subtitle = "Bars = mean ± SE; points = samples; colored labels = MetaPhlAn Robust_DA (same direction)",
    x = NULL, y = "Relative abundance"
  )

if (has_ggh4x) {
  pA <- pA + ggh4x::facet_grid2(
    focus_direction ~ .,
    scales = "free_y", space = "free_y",
    strip = ggh4x::strip_themed(
      background_y = ggh4x::elem_list_rect(
        fill = unname(STRIP_BG[c("Normal", "PoorPlaque")])
      ),
      text_y = ggh4x::elem_list_text(
        colour = unname(COL[c("Normal", "PoorPlaque")]),
        face = "bold", size = BASE_PT - 1
      )
    )
  )
} else {
  pA <- pA +
    facet_grid(
      focus_direction ~ .,
      scales = "free_y", space = "free_y",
      labeller = labeller(focus_direction = strip_labs)
    ) +
    theme(strip.text.y = element_markdown(size = BASE_PT - 1))
}

pA <- pA +
  theme(
    axis.text.y = element_markdown(size = BASE_PT - 4, colour = "black", lineheight = 1.05),
    axis.text.x = element_text(size = BASE_PT - 2),
    legend.position = "bottom",
    legend.text = element_text(size = BASE_PT - 2),
    legend.title = element_text(size = BASE_PT - 1),
    strip.background = element_rect(fill = "grey95", colour = NA),
    plot.margin = margin(8, 10, 8, 8)
  )

stem <- file.path(out_fig, "Fig_panelA_RobustDA_bar_scatter")
message("[3/3] Export figures → ", stem)
save_pub(pA, stem, width_mm = 180, height_mm = 300)

# sync to legacy output path + TSVs
for (ext in c("pdf", "png", "jpg", "svg")) {
  file.copy(paste0(stem, ".", ext),
            file.path(sync_fig, paste0("Fig_panelA_RobustDA_bar_scatter.", ext)),
            overwrite = TRUE)
}
file.copy(file.path(out_data, "Fig_panelA_data_points.tsv"),
          file.path(sync_fig, "Fig_panelA_data_points.tsv"), overwrite = TRUE)
file.copy(file.path(out_data, "Fig_panelA_data_summary.tsv"),
          file.path(sync_fig, "Fig_panelA_data_summary.tsv"), overwrite = TRUE)

message("[DONE] Panel A optimized figure written.")
message("  figures: ", out_fig)
message("  synced : ", sync_fig)
