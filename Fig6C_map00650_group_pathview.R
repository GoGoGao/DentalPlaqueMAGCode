## Packaged for revise/; source: project/writing/16_map00650_hubMAG_KO_presence/scripts/16b_map00650_group_pathview.R
## Figures: Fig6C
## Note: group pathview
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig6C_map00650_group_pathview.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# 16b: Group-level pathview of map00650 KO presence
#   Normal- vs PoorPlaque-enriched Robust_DA MAGs
# Panels:
#   A — presence proportion in Normal MAGs
#   B — presence proportion in PoorPlaque MAGs
#   C — Δprop (PP − Normal); teal=Normal-biased, red=PP-biased
#   D — lollipop of top |Δprop| / FDR-significant KOs
#   Composite: A|B|C (top) + D (bottom)
# Font: Times New Roman, plain
# KO profile: sgb/compare_genomics/total.KO.xls

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(png)
  library(pathview)
})

set.seed(42)
options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

## ---------- paths ------------------------------------------------------------
ko_path    <- file.path(.pkg_root, "01this_study_dataset/sgb/compare_genomics/total.KO.xls"
mag_path   <- file.path(.pkg_root, "01this_study_dataset/metab_corr/paper_result/00_RobustDA_MAG_set/tables/RobustDA_MAGs.tsv"
ko2pw_path <- file.path(.pkg_root, "01this_study_dataset/sgb/compare_genomics_all/00_data/kegg_ko_pathway.tab"
fig3_path  <- file.path(.pkg_root, "writing/15_map00650_KO_metab_pathway/tables/map00650_KO_nodes.tsv"
module_dir <- file.path(.pkg_root, "writing/16_map00650_hubMAG_KO_presence"
out_fig    <- file.path(module_dir, "figures")
out_tab    <- file.path(module_dir, "tables")
out_data   <- file.path(module_dir, "data")
cache_dir  <- file.path(module_dir, "cache")
for (d in c(out_fig, out_tab, out_data, cache_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
source(file.path(.pkg_root, "01this_study_dataset/metab_corr/paper_result/_shared/export_utils.R")

## ---------- constants --------------------------------------------------------
PATHWAY_NUM <- "00650"
PATHWAY_ID  <- "map00650"
FONT        <- "Times New Roman"
BASE_SIZE   <- 9
TXT_COL     <- "black"
COL_PP      <- "#E64B35"
COL_NM      <- "#00A087"
FIG_IN      <- 180 / 25.4

theme_set(
  theme_classic(base_size = BASE_SIZE, base_family = FONT) +
    theme(
      text = element_text(family = FONT, colour = TXT_COL, face = "plain", size = BASE_SIZE),
      plot.title = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL),
      plot.subtitle = element_text(family = FONT, face = "plain", size = BASE_SIZE - 1, colour = "grey30"),
      axis.title = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL),
      axis.text = element_text(family = FONT, face = "plain", size = BASE_SIZE - 0.5, colour = TXT_COL),
      legend.title = element_text(family = FONT, face = "plain", size = BASE_SIZE - 0.5, colour = TXT_COL),
      legend.text = element_text(family = FONT, face = "plain", size = BASE_SIZE - 0.5, colour = TXT_COL),
      plot.tag = element_text(family = FONT, face = "plain", size = BASE_SIZE + 1, colour = TXT_COL),
      panel.grid = element_blank()
    )
)

parse_ko <- function(x) str_extract(as.character(x), "K\\d+")

save_raster_multi <- function(png_in, stem, width_in = FIG_IN * 0.55, dpi = 300) {
  img <- png::readPNG(png_in)
  h <- dim(img)[1]; w <- dim(img)[2]
  height_in <- width_in * (h / w)
  for (ext in c("png", "jpg", "pdf", "svg")) {
    f <- paste0(stem, ".", ext)
    tryCatch({
      if (ext == "png") {
        file.copy(png_in, f, overwrite = TRUE)
      } else if (ext == "jpg") {
        grDevices::jpeg(f, width = w, height = h, units = "px", quality = 95)
        grid::grid.raster(img); grDevices::dev.off()
      } else if (ext == "pdf") {
        grDevices::cairo_pdf(f, width = width_in, height = height_in)
        grid::grid.raster(img); grDevices::dev.off()
      } else {
        svglite::svglite(f, width = width_in, height = height_in)
        grid::grid.raster(img); grDevices::dev.off()
      }
      message("[OK] ", basename(f))
    }, error = function(e) message("[ERR] ", basename(f), ": ", e$message))
  }
  list(w = width_in, h = height_in)
}

bh_fdr <- function(p) {
  p <- as.numeric(p)
  n <- length(p)
  o <- order(p)
  r <- rep(NA_real_, n)
  # BH: q_(i) = min_{j>=i} (n/j * p_(j)) on ascending ranks
  r[o] <- pmin(1, rev(cummin(rev(p[o] * n / seq_len(n)))))
  r
}

## ---------- [1] Group presence on map00650 -----------------------------------
message("[1] Compute group-level map00650 KO presence ...")

mags <- read_tsv(mag_path, show_col_types = FALSE)
dir_col <- if ("Robust_direction" %in% names(mags)) "Robust_direction" else "enrichment"
nm_mags <- mags %>% filter(.data[[dir_col]] == "Normal") %>% pull(MAG)
pp_mags <- mags %>% filter(.data[[dir_col]] == "PoorPlaque") %>% pull(MAG)

ko_raw <- read_tsv(ko_path, show_col_types = FALSE)
ko_ids <- parse_ko(ko_raw$KO)
ko_names <- setNames(sub("^[^:]+:", "", as.character(ko_raw$KO)), ko_ids)

nm_mags <- intersect(nm_mags, colnames(ko_raw))
pp_mags <- intersect(pp_mags, colnames(ko_raw))
stopifnot(length(nm_mags) > 0, length(pp_mags) > 0)

ko2pw <- read_tsv(ko2pw_path, col_names = c("path", "ko"), show_col_types = FALSE) %>%
  mutate(
    pathway_id = str_replace(str_replace(path, "^path:", ""), "^ko", "map"),
    KO = str_replace(ko, "^ko:", "")
  ) %>%
  filter(pathway_id == PATHWAY_ID) %>%
  distinct(KO)
m650 <- intersect(ko2pw$KO, ko_ids)
message(sprintf("  Robust_DA MAGs: Normal=%d, PoorPlaque=%d; map00650 KOs=%d",
                length(nm_mags), length(pp_mags), length(m650)))

fisher_one <- function(n_pp, N_pp, n_nm, N_nm) {
  tab <- matrix(c(n_pp, N_pp - n_pp, n_nm, N_nm - n_nm), nrow = 2, byrow = TRUE)
  ft <- tryCatch(fisher.test(tab), error = function(e) list(estimate = NA_real_, p.value = 1))
  list(OR = unname(ft$estimate), p = ft$p.value)
}

rows <- lapply(m650, function(k) {
  idx <- which(ko_ids == k)[1]
  n_nm <- sum(ko_raw[idx, nm_mags] > 0, na.rm = TRUE)
  n_pp <- sum(ko_raw[idx, pp_mags] > 0, na.rm = TRUE)
  ft <- fisher_one(n_pp, length(pp_mags), n_nm, length(nm_mags))
  tibble(
    KO = k,
    KO_name = unname(ko_names[k]),
    n_Normal = n_nm, N_Normal = length(nm_mags),
    prop_Normal = n_nm / length(nm_mags),
    n_PoorPlaque = n_pp, N_PoorPlaque = length(pp_mags),
    prop_PoorPlaque = n_pp / length(pp_mags),
    delta_prop = (n_pp / length(pp_mags)) - (n_nm / length(nm_mags)),
    OR_PP_vs_N = ft$OR,
    Fisher_p = ft$p
  )
})
pres <- bind_rows(rows) %>%
  mutate(Fisher_FDR = bh_fdr(Fisher_p))

if (file.exists(fig3_path)) {
  fig3 <- read_tsv(fig3_path, show_col_types = FALSE) %>%
    mutate(KO = parse_ko(KO)) %>%
    distinct(KO, .keep_all = TRUE) %>%
    select(KO, community_DA = Robust_direction, community_log2FC = log2FC)
  pres <- left_join(pres, fig3, by = "KO")
} else {
  pres <- pres %>% mutate(community_DA = NA_character_, community_log2FC = NA_real_)
}

pres <- pres %>% arrange(desc(delta_prop))
write_tsv(pres, file.path(out_tab, "map00650_KO_presence_by_MAG_group.tsv"))
write_tsv(pres, file.path(out_data, "map00650_KO_presence_by_MAG_group.tsv"))

sig <- pres %>% filter(Fisher_FDR < 0.05)
message(sprintf("  FDR<0.05 differential KOs: %d", nrow(sig)))

## ---------- [2] pathview overlays -------------------------------------------
message("[2] Render group pathview panels ...")

# copy XML cache if available
xml_src <- file.path(cache_dir, "ko00650.xml")
xml15 <- file.path(.pkg_root, "writing/15_map00650_KO_metab_pathway/cache/ko00650.xml"
if (!file.exists(xml_src) && file.exists(xml15)) file.copy(xml15, xml_src)

run_pv <- function(gene_vec, suffix, low, mid, high, lim) {
  # drop all-NA
  g <- gene_vec[!is.na(gene_vec)]
  if (!length(g)) return(NULL)
  old <- getwd(); setwd(cache_dir)
  on.exit(setwd(old), add = TRUE)
  tryCatch(
    pathview(
      gene.data = g,
      pathway.id = PATHWAY_NUM,
      species = "ko",
      gene.idtype = "KEGG",
      keys.align = "y",
      limit = list(gene = lim),
      bins = list(gene = 10),
      low = list(gene = low),
      mid = list(gene = mid),
      high = list(gene = high),
      na.col = "white",
      same.layer = FALSE,
      kegg.native = TRUE,
      out.suffix = suffix
    ),
    error = function(e) {
      message("[ERR] pathview ", suffix, ": ", e$message)
      NULL
    }
  )
}

find_pv_png <- function(suffix) {
  cand <- list.files(cache_dir, pattern = paste0("00650.*", suffix, ".*\\.png$"),
                     full.names = TRUE)
  if (!length(cand)) cand <- list.files(cache_dir, pattern = paste0(suffix, ".*\\.png$"),
                                        full.names = TRUE)
  if (length(cand)) cand[[1]] else NA_character_
}

# A: Normal presence prop (0-1), white→teal
g_nm <- setNames(pres$prop_Normal, pres$KO)
g_nm[g_nm == 0] <- NA_real_
pv_nm <- run_pv(g_nm, "grpNormal", low = "grey95", mid = "#7FCDC0", high = COL_NM, lim = c(0, 1))

# B: PP presence prop (0-1), white→red
g_pp <- setNames(pres$prop_PoorPlaque, pres$KO)
g_pp[g_pp == 0] <- NA_real_
pv_pp <- run_pv(g_pp, "grpPP", low = "grey95", mid = "#F4A582", high = COL_PP, lim = c(0, 1))

# C: delta_prop, symmetric around 0
g_d <- setNames(pres$delta_prop, pres$KO)
# keep zeros as mid (shared / no bias); only truly unused pathway KOs absent from both → still 0
# For readability: set KOs with both props==0 to NA (not on any Robust_DA MAG)
both0 <- pres$prop_Normal == 0 & pres$prop_PoorPlaque == 0
g_d[both0] <- NA_real_
dlim <- max(0.55, ceiling(max(abs(g_d), na.rm = TRUE) * 20) / 20)
pv_d <- run_pv(g_d, "grpDelta", low = COL_NM, mid = "grey90", high = COL_PP, lim = c(-dlim, dlim))

export_pv <- function(suffix, stem_name, tsv_df) {
  png <- find_pv_png(suffix)
  if (is.na(png) || !file.exists(png)) {
    message("[WARN] missing pathview PNG for ", suffix)
    return(invisible(NULL))
  }
  stem <- file.path(out_fig, stem_name)
  save_raster_multi(png, stem, width_in = FIG_IN * 0.52)
  write_tsv(tsv_df, paste0(stem, ".tsv"))
  write_tsv(tsv_df, file.path(out_data, paste0(stem_name, ".tsv")))
  png
}

png_A <- export_pv(
  "grpNormal", "Fig_map00650_group_pathview_Normal",
  pres %>% transmute(KO, KO_name, color_value = prop_Normal, n = n_Normal, N = N_Normal,
                     panel = "Normal_presence_prop")
)
png_B <- export_pv(
  "grpPP", "Fig_map00650_group_pathview_PoorPlaque",
  pres %>% transmute(KO, KO_name, color_value = prop_PoorPlaque, n = n_PoorPlaque, N = N_PoorPlaque,
                     panel = "PoorPlaque_presence_prop")
)
png_C <- export_pv(
  "grpDelta", "Fig_map00650_group_pathview_delta",
  pres %>% transmute(KO, KO_name, color_value = delta_prop,
                     prop_Normal, prop_PoorPlaque, Fisher_p, Fisher_FDR, community_DA,
                     panel = "delta_prop_PP_minus_Normal")
)

## ---------- [3] 1×3 pathview grid -------------------------------------------
message("[3] Assemble A|B|C pathview grid ...")

pngs <- c(Normal = png_A, PoorPlaque = png_B, Delta = png_C)
pngs <- pngs[!is.na(pngs) & file.exists(pngs)]
stopifnot(length(pngs) == 3L)

stem_grid <- file.path(out_fig, "Fig_map00650_group_pathview_ABC")
imgs <- lapply(pngs, png::readPNG)
asp <- dim(imgs[[1]])[1] / dim(imgs[[1]])[2]
w_in <- FIG_IN * 1.15
h_in <- (w_in / 3) * asp * 1.12
titles <- c(
  Normal = sprintf("Normal MAGs (n=%d)\npresence proportion", length(nm_mags)),
  PoorPlaque = sprintf("PoorPlaque MAGs (n=%d)\npresence proportion", length(pp_mags)),
  Delta = sprintf("Δprop (PP − Normal)\nlimit ±%.2f", dlim)
)

for (ext in c("pdf", "svg", "png", "jpg")) {
  f <- paste0(stem_grid, ".", ext)
  tryCatch({
    if (ext == "pdf") grDevices::cairo_pdf(f, width = w_in, height = h_in, family = FONT)
    else if (ext == "svg") svglite::svglite(f, width = w_in, height = h_in)
    else if (ext == "png") grDevices::png(f, width = w_in, height = h_in, units = "in", res = 300)
    else grDevices::jpeg(f, width = w_in, height = h_in, units = "in", res = 300, quality = 95)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
      nrow = 2, ncol = 3,
      heights = grid::unit(c(0.10, 0.90), "npc"),
      widths = grid::unit(rep(1, 3), "null")
    )))
    grid::grid.text(
      "map00650 KO presence in Robust_DA MAGs by enrichment group",
      gp = grid::gpar(fontfamily = FONT, fontsize = BASE_SIZE + 1, col = TXT_COL),
      vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:3)
    )
    for (i in seq_along(imgs)) {
      grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = i))
      grid::grid.raster(imgs[[i]], width = 0.98, height = 0.88, interpolate = TRUE)
      grid::grid.text(
        paste0(LETTERS[i], "  ", titles[[names(imgs)[i]]]),
        x = 0.02, y = 0.99, just = c("left", "top"),
        gp = grid::gpar(fontfamily = FONT, fontsize = BASE_SIZE - 0.5, col = TXT_COL, lineheight = 0.95)
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
  bind_rows(
    pres %>% transmute(panel = "A_Normal", KO, value = prop_Normal, Fisher_FDR, community_DA),
    pres %>% transmute(panel = "B_PoorPlaque", KO, value = prop_PoorPlaque, Fisher_FDR, community_DA),
    pres %>% transmute(panel = "C_delta", KO, value = delta_prop, Fisher_FDR, community_DA)
  ),
  paste0(stem_grid, ".tsv")
)
write_tsv(
  bind_rows(
    pres %>% transmute(panel = "A_Normal", KO, value = prop_Normal),
    pres %>% transmute(panel = "B_PoorPlaque", KO, value = prop_PoorPlaque),
    pres %>% transmute(panel = "C_delta", KO, value = delta_prop)
  ),
  file.path(out_data, "Fig_map00650_group_pathview_ABC.tsv")
)

## ---------- [4] Companion lollipop ------------------------------------------
message("[4] Plot differential presence lollipop ...")

# Short names aligned with Figure 3B (map00650 KO lollipop); extras in same style
KO_NAME_SHORT <- c(
  # Fig.3B identical
  K00074 = "3-HB-CoA DH",
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
  # group-differential extras (Fig.3B style)
  K23351 = "methylmalonyl-CoA decarboxylase",
  K00174 = "2OG:ferredoxin OR alpha",
  K00175 = "2OG:ferredoxin OR beta",
  K03737 = "pyruvate:ferredoxin OR",
  K00929 = "butyrate kinase",
  K00135 = "SSA dehydrogenase",
  K00004 = "butanediol DH",
  K01652 = "acetolactate synthase L",
  K01653 = "acetolactate synthase S",
  K11258 = "acetolactate synthase ilvM",
  K00239 = "SDH/FRD flavoprotein",
  K00240 = "SDH/FRD Fe-S",
  K00242 = "SDH/FRD subunit D",
  K00248 = "butyryl-CoA DH",
  K01575 = "acetolactate decarboxylase",
  K01782 = "3-hydroxyacyl-CoA DH",
  K18369 = "alcohol DH (adh2)",
  K23259 = "alcohol DH",
  K01580 = "Glu decarboxylase"
)

# Full KEGG-style names for Fig.3B overlap (axis keeps short; TSV stores both)
KEGG_KO_NAME <- c(
  K00074 = "3-hydroxybutyryl-CoA dehydrogenase [EC:1.1.1.157]",
  K03821 = "poly[(R)-3-hydroxyalkanoate] polymerase subunit PhaC [EC:2.3.1.304]",
  K01715 = "short-chain-enoyl-CoA hydratase [EC:4.2.1.150]",
  K01641 = "hydroxymethylglutaryl-CoA synthase [EC:2.3.3.10]",
  K18122 = "4-hydroxybutyrate CoA-transferase [EC:2.8.3.-]",
  K01035 = "acetate CoA/acetoacetate CoA-transferase beta subunit [EC:2.8.3.8 2.8.3.9]",
  K01034 = "acetate CoA/acetoacetate CoA-transferase alpha subunit [EC:2.8.3.8 2.8.3.9]",
  K00823 = "4-aminobutyrate aminotransferase [EC:2.6.1.19]",
  K00244 = "succinate dehydrogenase flavoprotein subunit [EC:1.3.5.1]",
  K07250 = "4-aminobutyrate aminotransferase / (S)-3-amino-2-methylpropionate transaminase / 5-aminovalerate transaminase [EC:2.6.1.19 2.6.1.22 2.6.1.48]",
  K00247 = "succinate dehydrogenase subunit D",
  K00245 = "succinate dehydrogenase iron-sulfur subunit [EC:1.3.5.1]",
  K00246 = "succinate dehydrogenase subunit C"
)

lol <- pres %>%
  filter(Fisher_p < 0.05 | abs(delta_prop) >= 0.20) %>%
  mutate(
    KO_name = ifelse(KO %in% names(KEGG_KO_NAME), unname(KEGG_KO_NAME[KO]), KO_name),
    name_short = ifelse(
      KO %in% names(KO_NAME_SHORT),
      unname(KO_NAME_SHORT[KO]),
      str_trunc(
        str_trim(str_replace(KO_name, "\\s*\\[EC:[^\\]]+\\]\\s*$", "")),
        28L, side = "right", ellipsis = "..."
      )
    ),
    # Fig.3B axis format: "KO  short-name"
    label = paste0(KO, "  ", name_short),
    sig = case_when(
      Fisher_FDR < 0.05 ~ "FDR < 0.05",
      Fisher_p < 0.05 ~ "p < 0.05",
      TRUE ~ "|Δprop| ≥ 0.20"
    ),
    bias = ifelse(delta_prop >= 0, "PoorPlaque-biased", "Normal-biased")
  ) %>%
  arrange(delta_prop) %>%
  mutate(label = factor(label, levels = label))

pD <- ggplot(lol, aes(x = delta_prop, y = label, colour = bias)) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = delta_prop, yend = label), linewidth = 0.7) +
  geom_point(aes(shape = sig), size = 2.3) +
  scale_colour_manual(values = c("Normal-biased" = COL_NM, "PoorPlaque-biased" = COL_PP), name = NULL) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "p < 0.05" = 17, "|Δprop| ≥ 0.20" = 15), name = NULL) +
  scale_x_continuous(
    limits = c(-0.70, 0.70),
    breaks = c(-0.6, 0, 0.6),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  labs(
    x = "Δ presence proportion\n(PoorPlaque − Normal)",
    y = NULL,
    title = "Group-differential map00650 KO presence in Robust_DA MAGs",
    subtitle = NULL
  ) +
  theme(
    text = element_text(family = FONT, colour = TXT_COL, face = "plain"),
    plot.title = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL,
                              margin = margin(b = 2)),
    axis.title.x = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL,
                                lineheight = 0.95, margin = margin(t = 4)),
    axis.text = element_text(family = FONT, face = "plain", size = BASE_SIZE - 0.5, colour = TXT_COL),
    axis.text.x = element_text(family = FONT, face = "plain", size = BASE_SIZE - 0.5, colour = TXT_COL),
    legend.position = "top",
    legend.box = "vertical",
    legend.margin = margin(0, 0, -2, 0),
    legend.key.size = unit(0.28, "cm"),
    legend.text = element_text(family = FONT, face = "plain", size = BASE_SIZE - 0.5, colour = TXT_COL),
    plot.margin = margin(6, 10, 6, 4),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3)
  )

stem_D <- file.path(out_fig, "Fig_map00650_group_presence_lollipop")
# moderately tall aspect; x breaks thinned so labels are not cramped
save_gg(pD, stem_D, w = FIG_IN * 0.60, h = FIG_IN * 0.78, dpi = 300)
write_tsv(
  lol %>% select(KO, KO_name, name_short, prop_Normal, prop_PoorPlaque, delta_prop,
                 Fisher_p, Fisher_FDR, community_DA, bias, sig),
  paste0(stem_D, ".tsv")
)
write_tsv(
  lol %>% select(KO, KO_name, name_short, prop_Normal, prop_PoorPlaque, delta_prop,
                 Fisher_p, Fisher_FDR, community_DA, bias, sig),
  file.path(out_data, "Fig_map00650_group_presence_lollipop.tsv")
)

## ---------- [5] Full composite (ABC raster + D) ------------------------------
message("[5] Composite ABC + D ...")

# Render ABC grid png already; combine with ggplot D via patchwork + inset raster
# Simpler: write a tall composite with grid
stem_comp <- file.path(out_fig, "Fig_map00650_group_pathview_composite")
grid_png <- paste0(stem_grid, ".png")
# save D as temp raster via ggsave png then stack
tmp_d <- file.path(cache_dir, "tmp_lollipop.png")
ggsave(tmp_d, pD, width = FIG_IN * 0.60, height = FIG_IN * 0.78, dpi = 300, bg = "white")

img_abc <- png::readPNG(grid_png)
img_d <- png::readPNG(tmp_d)
w_in2 <- FIG_IN * 1.15
h_in2 <- w_in2 * 0.95

for (ext in c("pdf", "svg", "png", "jpg")) {
  f <- paste0(stem_comp, ".", ext)
  tryCatch({
    if (ext == "pdf") grDevices::cairo_pdf(f, width = w_in2, height = h_in2, family = FONT)
    else if (ext == "svg") svglite::svglite(f, width = w_in2, height = h_in2)
    else if (ext == "png") grDevices::png(f, width = w_in2, height = h_in2, units = "in", res = 300)
    else grDevices::jpeg(f, width = w_in2, height = h_in2, units = "in", res = 300, quality = 95)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
      2, 1, heights = grid::unit(c(0.58, 0.42), "npc")
    )))
    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::grid.raster(img_abc, width = 1, height = 1)
    grid::popViewport()
    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    grid::grid.raster(img_d, width = 0.92, height = 0.96)
    grid::grid.text("D", x = 0.02, y = 0.98, just = c("left", "top"),
                    gp = grid::gpar(fontfamily = FONT, fontsize = BASE_SIZE + 1, col = TXT_COL))
    grid::popViewport()
    grDevices::dev.off()
    message("[OK] ", basename(f))
  }, error = function(e) {
    try(grDevices::dev.off(), silent = TRUE)
    message("[ERR] ", basename(f), ": ", e$message)
  })
}

write_tsv(
  bind_rows(
    pres %>% transmute(panel = "A", KO, value = prop_Normal, Fisher_FDR),
    pres %>% transmute(panel = "B", KO, value = prop_PoorPlaque, Fisher_FDR),
    pres %>% transmute(panel = "C", KO, value = delta_prop, Fisher_FDR),
    lol %>% transmute(panel = "D", KO, value = delta_prop, Fisher_FDR)
  ),
  paste0(stem_comp, ".tsv")
)
write_tsv(
  bind_rows(
    pres %>% transmute(panel = "A", KO, value = prop_Normal),
    pres %>% transmute(panel = "B", KO, value = prop_PoorPlaque),
    pres %>% transmute(panel = "C", KO, value = delta_prop),
    lol %>% transmute(panel = "D", KO, value = delta_prop)
  ),
  file.path(out_data, "Fig_map00650_group_pathview_composite.tsv")
)

# panel index
write_tsv(
  tibble(
    panel = c("A", "B", "C", "D", "ABC", "Composite"),
    file = c(
      "Fig_map00650_group_pathview_Normal",
      "Fig_map00650_group_pathview_PoorPlaque",
      "Fig_map00650_group_pathview_delta",
      "Fig_map00650_group_presence_lollipop",
      "Fig_map00650_group_pathview_ABC",
      "Fig_map00650_group_pathview_composite"
    ),
    content = c(
      "Presence proportion in Normal Robust_DA MAGs",
      "Presence proportion in PoorPlaque Robust_DA MAGs",
      "Δprop (PoorPlaque − Normal)",
      "Lollipop of group-differential KOs",
      "A|B|C pathview strip",
      "ABC + D composite"
    )
  ),
  file.path(out_tab, "group_pathview_panel_index.tsv")
)

message("Done.")
message("  Figures: ", out_fig)
message("  KO source: ", ko_path)
