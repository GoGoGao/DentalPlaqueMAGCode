## Packaged for revise/; source: project/writing/16_map00650_hubMAG_KO_presence/scripts/16_map00650_hubMAG_KO_presence.R
## Figures: Fig6B;Fig6C
## Note: Δprop lollipop + pathview
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig6BC_map00650_MAG_KO_presence.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# 16: Hub MAG genomic capacity on map00650 (Butanoate metabolism)
# Claim: MAG3 (Normal hub) vs MAG95 (PoorPlaque hub) encode complementary
#        map00650 modules that echo Fig.3 community KO directions;
#        MAG83/MAG36 supply PoorPlaque-side community KO capacity.
# Panels:
#   A — presence heatmap (key KOs × MAG3/MAG95/MAG83/MAG36)
#   B — pathview presence overlay on the same Fig.3 KO node set (4 MAGs)
#   C — Panel-C Spearman rho: hub MAGs × map00650 metabolites
#   Composite A | (B / C)
# Font: Times New Roman, plain (not bold)
# Exports: PDF / SVG / PNG / JPG + TSV per panel and each subfigure

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(grid)
  library(png)
  library(pathview)
})

set.seed(42)
options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

## ---------- paths ------------------------------------------------------------
base_dir   <- file.path(.pkg_root, "01this_study_dataset/metab_corr"
mod06      <- file.path(base_dir, "paper_result/06_RobustDA_cooccurrence_network")
mod15      <- file.path(.pkg_root, "writing/15_map00650_KO_metab_pathway"
module_dir <- file.path(.pkg_root, "writing/16_map00650_hubMAG_KO_presence"
out_fig    <- file.path(module_dir, "figures")
out_tab    <- file.path(module_dir, "tables")
out_data   <- file.path(module_dir, "data")
cache_dir  <- file.path(module_dir, "cache")
for (d in c(out_fig, out_tab, out_data, cache_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
source(file.path(base_dir, "paper_result/_shared/export_utils.R"))

## ---------- constants --------------------------------------------------------
PATHWAY_ID   <- "map00650"
PATHWAY_NUM  <- "00650"
PATHWAY_NAME <- "Butanoate metabolism"
FONT         <- "Times New Roman"
BASE_SIZE    <- 9
TXT_COL      <- "black"
COL_PP       <- "#E64B35"
COL_NM       <- "#00A087"
COL_ABS      <- "#F0F0F0"
COL_PRES     <- "#3C5488"
FIG_MM       <- 180
FIG_IN       <- FIG_MM / 25.4

FOCUS_MAGS <- c("MAG3", "MAG95", "MAG83", "MAG36")
MAG_LABELS <- c(
  MAG3  = "Rothia dentocariosa (MAG3)",
  MAG95 = "Prevotella intermedia (MAG95)",
  MAG83 = "Unknown (MAG83)",
  MAG36 = "F0058 sp905372605 (MAG36)"
)
MAG_SHORT <- c(
  MAG3  = "MAG3\nRothia\n(Normal)",
  MAG95 = "MAG95\nPrevotella\n(PoorPlaque)",
  MAG83 = "MAG83\nUnknown\n(PoorPlaque)",
  MAG36 = "MAG36\nF0058\n(PoorPlaque)"
)
MAG_DIR <- c(
  MAG3 = "Normal", MAG95 = "PoorPlaque",
  MAG83 = "PoorPlaque", MAG36 = "PoorPlaque"
)

# Fig.3 Robust_DA KO nodes (same set as pathview gene layer)
FIG3_KOS <- c(
  "K00074", "K03821", "K01715", "K01641", "K18122", "K01035", "K01034",
  "K00823", "K00244", "K07250", "K00247", "K00245", "K00246"
)

# Hub-discriminative map00650 KOs (genomic contrast MAG3 vs MAG95 + carriers)
HUB_KOS <- c(
  "K07250", "K00135", "K00004",                 # MAG3 GABA / BDH module
  "K00174", "K00175", "K23351",                 # PP MAG Fisher markers
  "K00929", "K03737", "K00656",                 # MAG95 butyrate / pyruvate
  "K01580", "K00074",                           # MAG83 gad + hbd
  "K01034", "K01035", "K18122"                  # MAG36 ato / 4HB-CoA-T
)

KO_SHORT <- c(
  K00074 = "3-HB-CoA DH (hbd)",
  K03821 = "PHA polymerase",
  K01715 = "enoyl-CoA hydratase",
  K01641 = "HMG-CoA synthase",
  K18122 = "4-HB CoA-T",
  K01035 = "atoA (beta)",
  K01034 = "atoD (alpha)",
  K00823 = "GABA AT (puuE)",
  K00244 = "SDH flavoprotein",
  K07250 = "GABA AT (gabT)",
  K00247 = "SDH subunit D",
  K00245 = "SDH Fe-S",
  K00246 = "SDH subunit C",
  K00135 = "SSA DH (gabD)",
  K00004 = "BDH / butB",
  K00174 = "korA / oorA",
  K00175 = "korB / oorB",
  K23351 = "gcdC / mmdC",
  K00929 = "butyrate kinase (buk)",
  K03737 = "pyruvate:ferredoxin OR (por)",
  K00656 = "pyruvate formate-lyase",
  K01580 = "Glu decarboxylase (gad)"
)

KO_MODULE <- c(
  K07250 = "GABA shunt (Normal)",
  K00135 = "GABA shunt (Normal)",
  K00004 = "GABA / butanol (Normal)",
  K00823 = "GABA shunt (PP community)",
  K00244 = "SDH / FRD (Normal)",
  K00245 = "SDH / FRD (Normal)",
  K00246 = "SDH / FRD (Normal)",
  K00247 = "SDH / FRD (Normal)",
  K00174 = "2OG anaerobic (PP)",
  K00175 = "2OG anaerobic (PP)",
  K23351 = "2OG anaerobic (PP)",
  K00929 = "Butyrate (PP)",
  K03737 = "Pyruvate anaerobic (PP)",
  K00656 = "Pyruvate anaerobic (PP)",
  K01580 = "GABA production (PP)",
  K00074 = "Butanoate CoA (PP community)",
  K03821 = "Butanoate CoA (PP community)",
  K01715 = "Butanoate CoA (PP community)",
  K01641 = "Butanoate CoA (PP community)",
  K18122 = "Butanoate CoA (PP community)",
  K01035 = "Butanoate CoA (PP community)",
  K01034 = "Butanoate CoA (PP community)"
)

theme_set(
  theme_classic(base_size = BASE_SIZE, base_family = FONT) +
    theme(
      text = element_text(family = FONT, colour = TXT_COL, face = "plain", size = BASE_SIZE),
      plot.title = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL),
      plot.subtitle = element_text(family = FONT, face = "plain", size = BASE_SIZE - 1, colour = TXT_COL),
      axis.title = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL),
      axis.text = element_text(family = FONT, face = "plain", size = BASE_SIZE - 0.5, colour = TXT_COL),
      legend.title = element_text(family = FONT, face = "plain", size = BASE_SIZE - 0.5, colour = TXT_COL),
      legend.text = element_text(family = FONT, face = "plain", size = BASE_SIZE - 0.5, colour = TXT_COL),
      plot.tag = element_text(family = FONT, face = "plain", size = BASE_SIZE + 1, colour = TXT_COL),
      strip.text = element_text(family = FONT, face = "plain", size = BASE_SIZE - 0.5, colour = TXT_COL),
      strip.background = element_blank(),
      panel.grid = element_blank()
    )
)

## ---------- helpers ----------------------------------------------------------
parse_ko_id <- function(x) {
  x <- as.character(x)
  ifelse(str_detect(x, "^K\\d+"), str_extract(x, "^K\\d+"), x)
}

save_raster_multi <- function(png_in, stem, dpi = 300) {
  img <- png::readPNG(png_in)
  h <- dim(img)[1]; w <- dim(img)[2]
  asp <- h / w
  width_in <- FIG_IN * 0.55
  height_in <- width_in * asp
  for (ext in c("png", "jpg", "pdf", "svg")) {
    f <- paste0(stem, ".", ext)
    tryCatch({
      if (ext == "png") {
        png::writePNG(img, f)
      } else if (ext == "jpg") {
        grDevices::jpeg(f, width = w, height = h, units = "px", res = dpi, quality = 95)
        grid::grid.raster(img)
        grDevices::dev.off()
      } else if (ext == "pdf") {
        grDevices::cairo_pdf(f, width = width_in, height = height_in)
        grid::grid.raster(img)
        grDevices::dev.off()
      } else if (ext == "svg") {
        svglite::svglite(f, width = width_in, height = height_in)
        grid::grid.raster(img)
        grDevices::dev.off()
      }
      message("[OK] ", basename(f))
    }, error = function(e) message("[ERR] ", basename(f), ": ", e$message))
  }
  list(w = width_in, h = height_in, px_w = w, px_h = h)
}

## ---------- [1] Load KO matrix + Fig.3 node annotation -----------------------
message("[1] Build hub-MAG × map00650 key-KO presence matrix ...")

ko_nodes <- read_tsv(file.path(mod15, "tables/map00650_KO_nodes.tsv"), show_col_types = FALSE) %>%
  mutate(KO = parse_ko_id(KO))

ko_raw <- read_tsv(file.path(base_dir, "total.KO.xls"), show_col_types = FALSE)
ko_ids <- parse_ko_id(ko_raw$KO)
ko_full_names <- setNames(as.character(ko_raw$KO), ko_ids)
stopifnot(all(FOCUS_MAGS %in% colnames(ko_raw)))

# Union KO set: Fig.3 nodes + hub markers (unique, ordered)
key_kos <- unique(c(FIG3_KOS, HUB_KOS))
missing_in_mat <- setdiff(key_kos, ko_ids)
if (length(missing_in_mat)) {
  message("  KOs absent from total.KO.xls (treated as 0): ",
          paste(missing_in_mat, collapse = ", "))
}

pres_mat <- sapply(FOCUS_MAGS, function(m) {
  v <- integer(length(key_kos))
  names(v) <- key_kos
  for (k in key_kos) {
    if (k %in% ko_ids) {
      v[k] <- as.integer(ko_raw[[m]][match(k, ko_ids)] > 0)
    } else {
      v[k] <- 0L
    }
  }
  v
})
# ensure matrix rownames
if (is.null(dim(pres_mat))) {
  pres_mat <- matrix(pres_mat, ncol = length(FOCUS_MAGS),
                     dimnames = list(key_kos, FOCUS_MAGS))
}

fig3_dir <- setNames(ko_nodes$Robust_direction, ko_nodes$KO)
fig3_lfc <- setNames(ko_nodes$log2FC, ko_nodes$KO)

count_mat <- sapply(FOCUS_MAGS, function(m) {
  v <- integer(length(key_kos))
  names(v) <- key_kos
  for (k in key_kos) {
    if (k %in% ko_ids) {
      v[k] <- as.integer(ko_raw[[m]][match(k, ko_ids)])
    } else {
      v[k] <- 0L
    }
  }
  v
})

pres_long <- as.data.frame(pres_mat, stringsAsFactors = FALSE) %>%
  rownames_to_column("KO") %>%
  pivot_longer(-KO, names_to = "MAG", values_to = "present") %>%
  mutate(present = as.integer(present)) %>%
  left_join(
    as.data.frame(count_mat, stringsAsFactors = FALSE) %>%
      rownames_to_column("KO") %>%
      pivot_longer(-KO, names_to = "MAG", values_to = "gene_count"),
    by = c("KO", "MAG")
  ) %>%
  mutate(
    MAG = factor(MAG, levels = FOCUS_MAGS),
    MAG_label = unname(MAG_LABELS[as.character(MAG)]),
    MAG_short = unname(MAG_SHORT[as.character(MAG)]),
    MAG_direction = unname(MAG_DIR[as.character(MAG)]),
    KO_short = ifelse(KO %in% names(KO_SHORT), unname(KO_SHORT[KO]), KO),
    module = ifelse(KO %in% names(KO_MODULE), unname(KO_MODULE[KO]), "Other"),
    in_Fig3_pathview = KO %in% FIG3_KOS,
    community_direction = ifelse(KO %in% names(fig3_dir), unname(fig3_dir[KO]), NA_character_),
    community_log2FC = ifelse(KO %in% names(fig3_lfc), unname(fig3_lfc[KO]), NA_real_),
    presence_label = ifelse(present == 1L, "Present", "Absent")
  )

# row order for heatmap: Fig.3 by log2FC then hub-only markers
row_order <- {
  fig3_ord <- ko_nodes %>% arrange(desc(log2FC)) %>% pull(KO)
  hub_only <- setdiff(HUB_KOS, FIG3_KOS)
  unique(c(fig3_ord, hub_only))
}
row_order <- row_order[row_order %in% key_kos]

pres_long <- pres_long %>%
  mutate(
    KO_axis = paste0(KO, "  ", KO_short),
    KO_axis = factor(KO_axis, levels = rev(
      paste0(row_order, "  ", ifelse(row_order %in% names(KO_SHORT),
                                     unname(KO_SHORT[row_order]), row_order))
    )),
    module = factor(
      module,
      levels = c(
        "Butanoate CoA (PP community)",
        "GABA shunt (PP community)",
        "SDH / FRD (Normal)",
        "GABA shunt (Normal)",
        "GABA / butanol (Normal)",
        "2OG anaerobic (PP)",
        "Butyrate (PP)",
        "Pyruvate anaerobic (PP)",
        "GABA production (PP)",
        "Other"
      )
    )
  )

write_tsv(pres_long, file.path(out_tab, "hubMAG_map00650_KO_presence_long.tsv"))
write_tsv(pres_long, file.path(out_data, "hubMAG_map00650_KO_presence_long.tsv"))

pres_wide <- as.data.frame(pres_mat) %>%
  rownames_to_column("KO") %>%
  mutate(
    KO_short = ifelse(KO %in% names(KO_SHORT), unname(KO_SHORT[KO]), KO),
    module = ifelse(KO %in% names(KO_MODULE), unname(KO_MODULE[KO]), "Other"),
    in_Fig3_pathview = KO %in% FIG3_KOS,
    community_direction = ifelse(KO %in% names(fig3_dir), unname(fig3_dir[KO]), NA_character_),
    community_log2FC = ifelse(KO %in% names(fig3_lfc), unname(fig3_lfc[KO]), NA_real_)
  ) %>%
  relocate(KO, KO_short, module, in_Fig3_pathview, community_direction, community_log2FC)
write_tsv(pres_wide, file.path(out_tab, "hubMAG_map00650_KO_presence_wide.tsv"))
write_tsv(pres_wide, file.path(out_data, "hubMAG_map00650_KO_presence_wide.tsv"))

# summary: n present per MAG among Fig.3 / hub sets
sum_tbl <- pres_long %>%
  group_by(MAG, MAG_label, MAG_direction) %>%
  summarise(
    n_key_KO = n_distinct(KO),
    n_present = sum(present == 1L),
    n_Fig3_nodes = sum(in_Fig3_pathview),
    n_Fig3_present = sum(in_Fig3_pathview & present == 1L),
    n_PP_community_present = sum(community_direction == "PoorPlaque" & present == 1L, na.rm = TRUE),
    n_N_community_present = sum(community_direction == "Normal" & present == 1L, na.rm = TRUE),
    .groups = "drop"
  )
write_tsv(sum_tbl, file.path(out_tab, "hubMAG_map00650_presence_summary.tsv"))
write_tsv(sum_tbl, file.path(out_data, "hubMAG_map00650_presence_summary.tsv"))
message("  Presence summary:")
print(sum_tbl)

## ---------- [2] Panel A: presence heatmap ------------------------------------
message("[2] Plot presence heatmap ...")

# left annotation strip: community direction for Fig.3 KOs
ann_df <- pres_long %>%
  distinct(KO, KO_axis, community_direction, in_Fig3_pathview, module) %>%
  mutate(
    ann_col = case_when(
      community_direction == "PoorPlaque" ~ COL_PP,
      community_direction == "Normal" ~ COL_NM,
      TRUE ~ "grey70"
    )
  )

pA <- ggplot(pres_long, aes(x = MAG, y = KO_axis, fill = factor(present))) +
  geom_tile(colour = "white", linewidth = 0.35, width = 0.92, height = 0.92) +
  scale_fill_manual(
    values = c("0" = COL_ABS, "1" = COL_PRES),
    labels = c("0" = "Absent", "1" = "Present"),
    name = NULL
  ) +
  scale_x_discrete(labels = MAG_SHORT, expand = c(0, 0.05)) +
  scale_y_discrete(expand = c(0, 0.05)) +
  labs(
    x = NULL, y = NULL,
    title = "map00650 key KO presence in hub / carrier MAGs",
    subtitle = "Rows: Fig.3 Robust_DA nodes (top) + hub-discriminative KOs; teal/red marks = community DA direction"
  ) +
  theme(
    axis.text.x = element_text(size = BASE_SIZE - 1, lineheight = 0.9, colour = TXT_COL),
    axis.text.y = element_text(size = BASE_SIZE - 1.2, colour = TXT_COL),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "top",
    legend.margin = margin(0, 0, -4, 0),
    legend.key.size = unit(0.28, "cm"),
    plot.margin = margin(2, 4, 2, 2),
    plot.title = element_text(size = BASE_SIZE, colour = TXT_COL),
    plot.subtitle = element_text(size = BASE_SIZE - 1.5, colour = "grey30")
  )

# community-direction side tiles (as separate narrow plot)
ann_plot_df <- ann_df %>%
  mutate(x = "DA")
pA_ann <- ggplot(ann_plot_df, aes(x = x, y = KO_axis, fill = community_direction)) +
  geom_tile(colour = "white", linewidth = 0.2, width = 0.9, height = 0.92) +
  scale_fill_manual(
    values = c("PoorPlaque" = COL_PP, "Normal" = COL_NM),
    na.value = "grey80",
    name = "Community\nRobust_DA",
    breaks = c("PoorPlaque", "Normal")
  ) +
  scale_y_discrete(expand = c(0, 0.05)) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.28, "cm"),
    plot.margin = margin(28, 0, 2, 0)
  )

pA_comb <- pA_ann + pA + plot_layout(widths = c(0.08, 1))

stem_A <- file.path(out_fig, "Fig_hubMAG_map00650_presence_heatmap")
save_gg(pA_comb, stem_A, w = FIG_IN * 0.72, h = FIG_IN * 0.95, dpi = 300)
write_tsv(
  pres_long %>%
    select(KO, KO_short, module, in_Fig3_pathview, community_direction, community_log2FC,
           MAG, MAG_label, MAG_direction, present, gene_count),
  paste0(stem_A, ".tsv")
)
write_tsv(
  pres_long %>%
    select(KO, KO_short, module, in_Fig3_pathview, community_direction, community_log2FC,
           MAG, MAG_label, MAG_direction, present, gene_count),
  file.path(out_data, "Fig_hubMAG_map00650_presence_heatmap.tsv")
)

## ---------- [3] Panel B: pathview presence on Fig.3 KO nodes -----------------
message("[3] pathview multi-MAG presence on Fig.3 KO node set ...")

# gene.data matrix: rows = Fig.3 KOs, cols = MAGs; values 0/1
# Use only FIG3_KOS so node set matches Fig.3 pathview gene layer.
gene_mat <- sapply(FOCUS_MAGS, function(m) {
  v <- as.numeric(pres_mat[FIG3_KOS, m])
  names(v) <- FIG3_KOS
  v
})
# pathview prefers non-zero contrast; map Absent=0 Present=1, limit 0-1
# Replace 0 with a tiny NA-like display via na.col for zeros:
# keep 0/1; mid grey / high blue-ish via custom low/high
gene_mat_pv <- gene_mat
# pathview drops all-zero genes from colouring; keep them as 0

old_wd <- getwd()
setwd(cache_dir)
# copy XML cache from module 15 if available to avoid re-download
xml15 <- file.path(mod15, "cache/ko00650.xml")
if (file.exists(xml15) && !file.exists(file.path(cache_dir, "ko00650.xml"))) {
  file.copy(xml15, file.path(cache_dir, "ko00650.xml"), overwrite = FALSE)
}

pv_list <- list()
for (m in FOCUS_MAGS) {
  gdat <- setNames(as.numeric(gene_mat_pv[, m]), rownames(gene_mat_pv))
  # present=1, absent kept as NA so white (not mid-grey for zero)
  gdat[gdat == 0] <- NA_real_
  suffix <- paste0("pres_", m)
  pv_list[[m]] <- tryCatch(
    pathview(
      gene.data   = gdat,
      pathway.id  = PATHWAY_NUM,
      species     = "ko",
      gene.idtype = "KEGG",
      keys.align  = "y",
      limit       = list(gene = c(0, 1)),
      bins        = list(gene = 2),
      low         = list(gene = "grey90"),
      mid         = list(gene = COL_PRES),
      high        = list(gene = COL_PRES),
      na.col      = "white",
      same.layer  = FALSE,
      kegg.native = TRUE,
      out.suffix  = suffix
    ),
    error = function(e) {
      message("[ERR] pathview ", m, ": ", e$message)
      NULL
    }
  )
}
setwd(old_wd)

# Collect per-MAG pathview PNGs and export
pv_pngs <- list()
for (m in FOCUS_MAGS) {
  cand <- list.files(cache_dir, pattern = paste0("00650.*pres_", m, ".*\\.png$"),
                     full.names = TRUE)
  if (!length(cand)) {
    cand <- list.files(cache_dir, pattern = paste0("pres_", m, ".*\\.png$"),
                       full.names = TRUE)
  }
  if (length(cand)) {
    pv_pngs[[m]] <- cand[[1]]
    stem_m <- file.path(out_fig, paste0("Fig_hubMAG_map00650_pathview_", m))
    save_raster_multi(cand[[1]], stem_m, dpi = 300)
    write_tsv(
      tibble(
        MAG = m,
        MAG_label = MAG_LABELS[[m]],
        MAG_direction = MAG_DIR[[m]],
        KO = FIG3_KOS,
        KO_short = unname(KO_SHORT[FIG3_KOS]),
        community_direction = unname(fig3_dir[FIG3_KOS]),
        community_log2FC = unname(fig3_lfc[FIG3_KOS]),
        present = as.integer(gene_mat[FIG3_KOS, m])
      ),
      paste0(stem_m, ".tsv")
    )
  } else {
    message("[WARN] no pathview PNG for ", m)
  }
}

# Combined 2×2 pathview raster panel
message("  Assemble 2x2 pathview cartoon panel ...")
if (length(pv_pngs) == 4L) {
  imgs <- lapply(FOCUS_MAGS, function(m) png::readPNG(pv_pngs[[m]]))
  names(imgs) <- FOCUS_MAGS
  # draw with grid
  stem_B <- file.path(out_fig, "Fig_hubMAG_map00650_pathview_grid")
  # approximate aspect from first image
  asp <- dim(imgs[[1]])[1] / dim(imgs[[1]])[2]
  w_in <- FIG_IN * 0.95
  h_in <- w_in * asp * 0.55  # 2 rows, slightly compressed
  for (ext in c("pdf", "svg", "png", "jpg")) {
    f <- paste0(stem_B, ".", ext)
    tryCatch({
      if (ext == "pdf") {
        grDevices::cairo_pdf(f, width = w_in, height = h_in, family = FONT)
      } else if (ext == "svg") {
        svglite::svglite(f, width = w_in, height = h_in)
      } else if (ext == "png") {
        grDevices::png(f, width = w_in, height = h_in, units = "in", res = 300)
      } else {
        grDevices::jpeg(f, width = w_in, height = h_in, units = "in", res = 300,
                        quality = 95)
      }
      grid::grid.newpage()
      # title
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(
        nrow = 3, ncol = 2,
        heights = grid::unit(c(0.07, 0.465, 0.465), "npc"),
        widths = grid::unit(c(0.5, 0.5), "npc")
      )))
      grid::grid.text(
        "Fig.3 map00650 KO nodes: genomic presence in hub / carrier MAGs",
        gp = grid::gpar(fontfamily = FONT, fontsize = BASE_SIZE + 1, col = TXT_COL),
        vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
      )
      positions <- list(
        MAG3 = c(2, 1), MAG95 = c(2, 2),
        MAG83 = c(3, 1), MAG36 = c(3, 2)
      )
      for (m in FOCUS_MAGS) {
        r <- positions[[m]][1]; c <- positions[[m]][2]
        grid::pushViewport(grid::viewport(layout.pos.row = r, layout.pos.col = c))
        grid::grid.raster(imgs[[m]], width = 0.96, height = 0.88, interpolate = TRUE)
        grid::grid.text(
          MAG_LABELS[[m]],
          x = 0.02, y = 0.98, just = c("left", "top"),
          gp = grid::gpar(fontfamily = FONT, fontsize = BASE_SIZE - 1, col = TXT_COL)
        )
        grid::popViewport()
      }
      grDevices::dev.off()
      message("[OK] ", basename(f))
    }, error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
      message("[ERR] ", basename(f), ": ", e$message)
    })
  }
  write_tsv(
    tibble(
      MAG = rep(FOCUS_MAGS, each = length(FIG3_KOS)),
      MAG_label = unname(MAG_LABELS[rep(FOCUS_MAGS, each = length(FIG3_KOS))]),
      MAG_direction = unname(MAG_DIR[rep(FOCUS_MAGS, each = length(FIG3_KOS))]),
      KO = rep(FIG3_KOS, times = length(FOCUS_MAGS)),
      KO_short = unname(KO_SHORT[rep(FIG3_KOS, times = length(FOCUS_MAGS))]),
      community_direction = unname(fig3_dir[rep(FIG3_KOS, times = length(FOCUS_MAGS))]),
      community_log2FC = unname(fig3_lfc[rep(FIG3_KOS, times = length(FOCUS_MAGS))]),
      present = as.integer(as.vector(gene_mat[FIG3_KOS, FOCUS_MAGS]))
    ),
    paste0(stem_B, ".tsv")
  )
  write_tsv(
    tibble(
      MAG = rep(FOCUS_MAGS, each = length(FIG3_KOS)),
      KO = rep(FIG3_KOS, times = length(FOCUS_MAGS)),
      present = as.integer(as.vector(gene_mat[FIG3_KOS, FOCUS_MAGS]))
    ),
    file.path(out_data, "Fig_hubMAG_map00650_pathview_grid.tsv")
  )
}

# Also export pathview gene input matrix
write_tsv(
  as.data.frame(gene_mat) %>% rownames_to_column("KO"),
  file.path(out_tab, "hubMAG_Fig3KO_presence_matrix.tsv")
)
write_tsv(
  as.data.frame(gene_mat) %>% rownames_to_column("KO"),
  file.path(out_data, "hubMAG_Fig3KO_presence_matrix.tsv")
)

## ---------- [3b] Extended pathview: Fig.3 nodes + hub markers on map --------
message("[3b] Extended pathview (Fig.3 + hub-marker KOs on map00650) ...")

ext_kos <- unique(c(FIG3_KOS, HUB_KOS))
gene_mat_ext <- sapply(FOCUS_MAGS, function(m) {
  v <- as.numeric(pres_mat[ext_kos, m])
  names(v) <- ext_kos
  v
})
write_tsv(
  as.data.frame(gene_mat_ext) %>% rownames_to_column("KO"),
  file.path(out_tab, "hubMAG_extendedKO_presence_matrix.tsv")
)
write_tsv(
  as.data.frame(gene_mat_ext) %>% rownames_to_column("KO"),
  file.path(out_data, "hubMAG_extendedKO_presence_matrix.tsv")
)

old_wd <- getwd()
setwd(cache_dir)
pv_ext_pngs <- list()
for (m in FOCUS_MAGS) {
  gdat <- setNames(as.numeric(gene_mat_ext[, m]), rownames(gene_mat_ext))
  gdat[gdat == 0] <- NA_real_
  suffix <- paste0("presExt_", m)
  tryCatch(
    pathview(
      gene.data   = gdat,
      pathway.id  = PATHWAY_NUM,
      species     = "ko",
      gene.idtype = "KEGG",
      keys.align  = "y",
      limit       = list(gene = c(0, 1)),
      bins        = list(gene = 2),
      low         = list(gene = "grey90"),
      mid         = list(gene = COL_PRES),
      high        = list(gene = COL_PRES),
      na.col      = "white",
      same.layer  = FALSE,
      kegg.native = TRUE,
      out.suffix  = suffix
    ),
    error = function(e) message("[ERR] pathview ext ", m, ": ", e$message)
  )
  cand <- list.files(cache_dir, pattern = paste0("00650.*presExt_", m, ".*\\.png$"),
                     full.names = TRUE)
  if (length(cand)) {
    pv_ext_pngs[[m]] <- cand[[1]]
    stem_m <- file.path(out_fig, paste0("Fig_hubMAG_map00650_pathviewExt_", m))
    save_raster_multi(cand[[1]], stem_m, dpi = 300)
    write_tsv(
      tibble(
        MAG = m,
        MAG_label = MAG_LABELS[[m]],
        MAG_direction = MAG_DIR[[m]],
        KO = ext_kos,
        KO_short = ifelse(ext_kos %in% names(KO_SHORT), unname(KO_SHORT[ext_kos]), ext_kos),
        in_Fig3_pathview = ext_kos %in% FIG3_KOS,
        present = as.integer(gene_mat_ext[ext_kos, m])
      ),
      paste0(stem_m, ".tsv")
    )
  }
}
setwd(old_wd)

if (length(pv_ext_pngs) == 4L) {
  imgs <- lapply(FOCUS_MAGS, function(m) png::readPNG(pv_ext_pngs[[m]]))
  names(imgs) <- FOCUS_MAGS
  stem_B2 <- file.path(out_fig, "Fig_hubMAG_map00650_pathviewExt_grid")
  asp <- dim(imgs[[1]])[1] / dim(imgs[[1]])[2]
  w_in <- FIG_IN * 0.95
  h_in <- w_in * asp * 0.55
  for (ext in c("pdf", "svg", "png", "jpg")) {
    f <- paste0(stem_B2, ".", ext)
    tryCatch({
      if (ext == "pdf") {
        grDevices::cairo_pdf(f, width = w_in, height = h_in, family = FONT)
      } else if (ext == "svg") {
        svglite::svglite(f, width = w_in, height = h_in)
      } else if (ext == "png") {
        grDevices::png(f, width = w_in, height = h_in, units = "in", res = 300)
      } else {
        grDevices::jpeg(f, width = w_in, height = h_in, units = "in", res = 300,
                        quality = 95)
      }
      grid::grid.newpage()
      grid::pushViewport(grid::viewport(layout = grid::grid.layout(
        nrow = 3, ncol = 2,
        heights = grid::unit(c(0.07, 0.465, 0.465), "npc"),
        widths = grid::unit(c(0.5, 0.5), "npc")
      )))
      grid::grid.text(
        "map00650 genomic presence: Fig.3 Robust_DA nodes + hub-marker KOs",
        gp = grid::gpar(fontfamily = FONT, fontsize = BASE_SIZE + 1, col = TXT_COL),
        vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:2)
      )
      positions <- list(
        MAG3 = c(2, 1), MAG95 = c(2, 2),
        MAG83 = c(3, 1), MAG36 = c(3, 2)
      )
      for (m in FOCUS_MAGS) {
        r <- positions[[m]][1]; c <- positions[[m]][2]
        grid::pushViewport(grid::viewport(layout.pos.row = r, layout.pos.col = c))
        grid::grid.raster(imgs[[m]], width = 0.96, height = 0.88, interpolate = TRUE)
        grid::grid.text(
          MAG_LABELS[[m]],
          x = 0.02, y = 0.98, just = c("left", "top"),
          gp = grid::gpar(fontfamily = FONT, fontsize = BASE_SIZE - 1, col = TXT_COL)
        )
        grid::popViewport()
      }
      grDevices::dev.off()
      message("[OK] ", basename(f))
    }, error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
      message("[ERR] ", basename(f), ": ", e$message)
    })
  }
  write_tsv(
    tibble(
      MAG = rep(FOCUS_MAGS, each = length(ext_kos)),
      MAG_label = unname(MAG_LABELS[rep(FOCUS_MAGS, each = length(ext_kos))]),
      KO = rep(ext_kos, times = length(FOCUS_MAGS)),
      in_Fig3_pathview = rep(ext_kos %in% FIG3_KOS, times = length(FOCUS_MAGS)),
      present = as.integer(as.vector(gene_mat_ext[ext_kos, FOCUS_MAGS]))
    ),
    paste0(stem_B2, ".tsv")
  )
  write_tsv(
    tibble(
      MAG = rep(FOCUS_MAGS, each = length(ext_kos)),
      KO = rep(ext_kos, times = length(FOCUS_MAGS)),
      present = as.integer(as.vector(gene_mat_ext[ext_kos, FOCUS_MAGS]))
    ),
    file.path(out_data, "Fig_hubMAG_map00650_pathviewExt_grid.tsv")
  )
}

## ---------- [4] Panel C: MAG × map00650 metabolite network rho ---------------
message("[4] Plot hub-MAG × map00650 metabolite association heatmap ...")

m650_metab <- c(
  metab_17037 = "2-Oxoglutarate",
  metab_16879 = "Fumarate",
  metab_679   = "4-Aminobutanoate",
  metab_28050 = "Succinate",
  metab_27177 = "L-Glutamate",
  metab_542   = "L-Glutamate",
  metab_26768 = "Pyruvate",
  metab_13490 = "Butanal"
)

edges <- read_tsv(file.path(mod06, "data/Fig_panelC_edges.tsv"), show_col_types = FALSE) %>%
  filter(MAG %in% FOCUS_MAGS, metab_id %in% names(m650_metab)) %>%
  mutate(
    metab_name = unname(m650_metab[metab_id]),
    # collapse duplicate glutamate features for display: keep both IDs in TSV,
    # but for plot use metab_name; if both edges exist for same MAG-name, keep max |rho|
  )

edge_plot <- edges %>%
  group_by(MAG, metab_name) %>%
  slice_max(order_by = abs(rho), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    MAG = factor(MAG, levels = FOCUS_MAGS),
    metab_name = factor(
      metab_name,
      levels = c("2-Oxoglutarate", "4-Aminobutanoate", "L-Glutamate",
                 "Pyruvate", "Succinate", "Fumarate", "Butanal")
    )
  )

# complete grid for absent edges (show as NA / light)
grid_df <- expand_grid(
  MAG = factor(FOCUS_MAGS, levels = FOCUS_MAGS),
  metab_name = factor(
    c("2-Oxoglutarate", "4-Aminobutanoate", "L-Glutamate",
      "Pyruvate", "Succinate", "Fumarate", "Butanal"),
    levels = c("2-Oxoglutarate", "4-Aminobutanoate", "L-Glutamate",
               "Pyruvate", "Succinate", "Fumarate", "Butanal")
  )
) %>%
  left_join(edge_plot %>% select(MAG, metab_name, rho, FDR, edge_sign, metab_id),
            by = c("MAG", "metab_name"))

write_tsv(
  edges %>%
    transmute(MAG, species_plain, focus_direction, metab_id, Metabolite,
              map00650_name = unname(m650_metab[metab_id]),
              rho, p, FDR, edge_sign, same_direction),
  file.path(out_tab, "hubMAG_map00650_metab_edges.tsv")
)
write_tsv(
  grid_df,
  file.path(out_data, "hubMAG_map00650_metab_rho_grid.tsv")
)

pC <- ggplot(grid_df, aes(x = MAG, y = metab_name, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.4, width = 0.92, height = 0.92) +
  geom_text(
    aes(label = ifelse(is.na(rho), "", sprintf("%+.2f", rho))),
    family = FONT, size = 2.6, colour = TXT_COL
  ) +
  scale_fill_gradient2(
    low = COL_NM, mid = "white", high = COL_PP,
    midpoint = 0, na.value = "grey95",
    limits = c(-1, 1),
    name = expression(Spearman ~ rho)
  ) +
  scale_x_discrete(labels = MAG_SHORT, expand = c(0, 0.05)) +
  scale_y_discrete(expand = c(0, 0.05)) +
  labs(
    x = NULL, y = NULL,
    title = "Panel-C significant edges: hub MAGs × map00650 metabolites",
    subtitle = "|ρ| ≥ 0.80, FDR ≤ 0.05; blank = no significant edge"
  ) +
  theme(
    axis.text.x = element_text(size = BASE_SIZE - 1, lineheight = 0.9),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    legend.key.height = unit(0.35, "cm"),
    legend.key.width = unit(0.25, "cm"),
    plot.margin = margin(2, 4, 2, 2),
    plot.subtitle = element_text(size = BASE_SIZE - 1.5, colour = "grey30")
  )

stem_C <- file.path(out_fig, "Fig_hubMAG_map00650_metab_rho")
save_gg(pC, stem_C, w = FIG_IN * 0.72, h = FIG_IN * 0.42, dpi = 300)
write_tsv(grid_df, paste0(stem_C, ".tsv"))
write_tsv(grid_df, file.path(out_data, "Fig_hubMAG_map00650_metab_rho.tsv"))

## ---------- [5] Module cartoon (compact genomic logic) -----------------------
message("[5] Plot module cartoon ...")

cartoon <- tribble(
  ~MAG, ~module_lab, ~kos, ~role,
  "MAG3",  "GABA shunt",           "gabT (K07250)\ngabD (K00135)\nBDH (K00004)", "Normal hub",
  "MAG95", "2OG / butyrate",       "korA/B (K00174/75)\nmmdC (K23351)\nbuk (K00929)\npor (K03737)", "PoorPlaque hub",
  "MAG83", "GABA production +\nhbd", "gad (K01580)\nhbd (K00074)", "PP carrier",
  "MAG36", "ato / 4HB-CoA-T",      "atoD/A (K01034/35)\nabfT (K18122)", "PP carrier"
) %>%
  mutate(
    MAG = factor(MAG, levels = FOCUS_MAGS),
    x = as.numeric(MAG),
    fill = ifelse(MAG_DIR[as.character(MAG)] == "Normal", COL_NM, COL_PP)
  )

pCart <- ggplot(cartoon, aes(x = x, y = 1)) +
  geom_tile(aes(fill = fill), width = 0.9, height = 0.85, colour = "white") +
  geom_text(aes(label = paste0(MAG, "\n", role)), y = 1.28,
            family = FONT, size = 2.7, colour = "white", lineheight = 0.95) +
  geom_text(aes(label = module_lab), y = 1.05,
            family = FONT, size = 2.6, colour = "white", fontface = "plain") +
  geom_text(aes(label = kos), y = 0.78,
            family = FONT, size = 2.2, colour = "white", lineheight = 0.92) +
  scale_fill_identity() +
  scale_x_continuous(breaks = 1:4, labels = FOCUS_MAGS, expand = c(0.02, 0.02)) +
  coord_cartesian(ylim = c(0.45, 1.45)) +
  labs(
    x = NULL, y = NULL,
    title = "Complementary map00650 modules in hub / carrier MAGs"
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    plot.margin = margin(4, 4, 2, 4)
  )

stem_cart <- file.path(out_fig, "Fig_hubMAG_map00650_module_cartoon")
save_gg(pCart, stem_cart, w = FIG_IN * 0.95, h = FIG_IN * 0.28, dpi = 300)
write_tsv(cartoon %>% select(MAG, module_lab, kos, role), paste0(stem_cart, ".tsv"))
write_tsv(cartoon %>% select(MAG, module_lab, kos, role),
          file.path(out_data, "Fig_hubMAG_map00650_module_cartoon.tsv"))

## ---------- [6] Composite ----------------------------------------------------
message("[6] Assemble composite ...")

# Use heatmap + rho + cartoon (pathview grid kept as standalone B for clarity)
p_comp <- (pA_comb / pC / pCart) +
  plot_layout(heights = c(1.55, 0.72, 0.38)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(family = FONT, face = "plain", size = BASE_SIZE + 1, colour = TXT_COL)
  )

stem_comp <- file.path(out_fig, "Fig_hubMAG_map00650_composite")
save_gg(p_comp, stem_comp, w = FIG_IN * 0.78, h = FIG_IN * 1.35, dpi = 300)

# composite source TSV = presence long + edge grid keys
write_tsv(
  bind_rows(
    pres_long %>%
      transmute(panel = "A_presence", feature = KO, feature_name = KO_short,
                MAG, value = as.numeric(present), detail = module),
    grid_df %>%
      transmute(panel = "C_rho", feature = as.character(metab_name),
                feature_name = as.character(metab_name),
                MAG = as.character(MAG), value = rho,
                detail = ifelse(is.na(edge_sign), "no_sig_edge", edge_sign)),
    cartoon %>%
      transmute(panel = "D_cartoon", feature = MAG, feature_name = module_lab,
                MAG, value = NA_real_, detail = kos)
  ),
  paste0(stem_comp, ".tsv")
)
write_tsv(
  bind_rows(
    pres_long %>%
      transmute(panel = "A_presence", feature = KO, MAG, value = as.numeric(present)),
    grid_df %>%
      transmute(panel = "C_rho", feature = as.character(metab_name),
                MAG = as.character(MAG), value = rho)
  ),
  file.path(out_data, "Fig_hubMAG_map00650_composite.tsv")
)

# Also a second composite that includes pathview grid if available
if (length(pv_pngs) == 4L && file.exists(paste0(file.path(out_fig, "Fig_hubMAG_map00650_pathview_grid"), ".png"))) {
  # export a "full" composite note file pointing to panels
  write_tsv(
    tibble(
      panel = c("A", "B", "C", "D"),
      file = c(
        "Fig_hubMAG_map00650_presence_heatmap",
        "Fig_hubMAG_map00650_pathview_grid",
        "Fig_hubMAG_map00650_metab_rho",
        "Fig_hubMAG_map00650_module_cartoon"
      ),
      content = c(
        "Key KO presence heatmap (Fig.3 nodes + hub markers)",
        "pathview presence on same Fig.3 KO node set (2x2 MAG grid)",
        "Significant Panel-C rho for map00650 metabolites",
        "Complementary genomic module cartoon"
      )
    ),
    file.path(out_tab, "composite_panel_index.tsv")
  )
}

## ---------- [7] Contrast stats (MAG3 vs MAG95 on key set) --------------------
message("[7] MAG3 vs MAG95 KO contingency on key set ...")

m3 <- pres_mat[key_kos, "MAG3"]
m95 <- pres_mat[key_kos, "MAG95"]
tab <- table(
  MAG3 = factor(m3, levels = c(0, 1), labels = c("absent", "present")),
  MAG95 = factor(m95, levels = c(0, 1), labels = c("absent", "present"))
)
ft <- fisher.test(tab)
contrast <- tibble(
  n_KO = length(key_kos),
  n_MAG3_only = sum(m3 == 1 & m95 == 0),
  n_MAG95_only = sum(m3 == 0 & m95 == 1),
  n_shared = sum(m3 == 1 & m95 == 1),
  n_neither = sum(m3 == 0 & m95 == 0),
  Fisher_OR = unname(ft$estimate),
  Fisher_p = ft$p.value,
  MAG3_only_KOs = paste(names(m3)[m3 == 1 & m95 == 0], collapse = ";"),
  MAG95_only_KOs = paste(names(m95)[m3 == 0 & m95 == 1], collapse = ";"),
  shared_KOs = paste(names(m3)[m3 == 1 & m95 == 1], collapse = ";")
)
write_tsv(contrast, file.path(out_tab, "MAG3_vs_MAG95_keyKO_contrast.tsv"))
write_tsv(contrast, file.path(out_data, "MAG3_vs_MAG95_keyKO_contrast.tsv"))
print(contrast)

message("Done. Figures in: ", out_fig)
message("Tables in: ", out_tab)
