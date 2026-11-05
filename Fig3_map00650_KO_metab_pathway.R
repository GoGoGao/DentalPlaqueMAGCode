## Packaged for revise/; source: project/01this_study_dataset/sgb/metab_corr/paper_result/15_map00650_KO_metab_pathway/scripts/15_map00650_KO_metab_pathway.R
## Figures: Fig3A;Fig3B;Fig3C
## Note: map00650 pathview + lollipops
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig3_map00650_KO_metab_pathway.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# 15: Cross-omics KEGG map00650 (Butanoate metabolism) — KO + metabolite overlay
# - pathview native KEGG map colored by KO direction (+1/−1) and metabolite log2FC
# - companion node panels (KO presence proportions; metabolite log2FC)
# - composite figure for manuscript
# Font: Times New Roman, plain (not bold)
# Exports: PDF / SVG / PNG / JPG + TSV per panel

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(purrr)
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
base_dir   <- file.path(.pkg_root, "01this_study_dataset/sgb/metab_corr"
comp_dir   <- file.path(dirname(base_dir), "compare_genomics_all")
mod03      <- file.path(base_dir, "paper_result/03_RobustDA_diffKO_heatmap")
mod07      <- file.path(base_dir, "paper_result/07_RobustDA_metab_KEGG_enrich")
mod12      <- file.path(base_dir, "paper_result/12_RobustDA_diffKO_KEGG_enrich")
humann_dir <- file.path(.pkg_root, "01this_study_dataset/M2-humann"
module_dir <- file.path(base_dir, "paper_result/15_map00650_KO_metab_pathway")
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
# Unified black text size for the whole figure (axes / titles / tags / legend)
BASE_SIZE    <- 10
TXT_COL      <- "black"
COL_PP       <- "#E64B35"
COL_NM       <- "#00A087"
COL_NS       <- "#D9D9D9"
GROUP_COLORS <- c("Normal" = COL_NM, "PoorPlaque" = COL_PP)
FN_TRUNC     <- 34L
# Nature-style square canvas: 170 mm × 170 mm
FIG_MM       <- 170
FIG_IN       <- FIG_MM / 25.4
# A (top) : B|C (bottom)
COMP_HEIGHTS <- c(1.65, 1)

# KEGG REST NAME (queried 2026-08-08) — override HUMAnN/metabolome synonyms
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
KEGG_CPD_NAME <- c(
  C00026 = "2-Oxoglutarate",
  C00122 = "Fumarate",
  C00334 = "4-Aminobutanoate",
  C00042 = "Succinate",
  C00025 = "L-Glutamate",
  C00022 = "Pyruvate",
  C01412 = "Butanal"
)

theme_set(
  theme_classic(base_size = BASE_SIZE, base_family = FONT) +
    theme(
      text = element_text(family = FONT, colour = TXT_COL, face = "plain", size = BASE_SIZE),
      plot.title = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL),
      plot.subtitle = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL),
      axis.title = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL),
      axis.text = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL),
      legend.title = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL),
      legend.text = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL),
      plot.tag = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL),
      strip.text = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL),
      strip.background = element_blank()
    )
)

parse_ko_id <- function(x) sub(":.*", "", trimws(as.character(x)))

encode_b64 <- function(path) {
  if (requireNamespace("base64enc", quietly = TRUE)) {
    return(base64enc::base64encode(path))
  }
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    raw <- readBin(path, what = "raw", n = file.info(path)$size)
    return(jsonlite::base64_enc(raw))
  }
  out <- system2("base64", args = c("-w", "0", path), stdout = TRUE)
  paste(out, collapse = "")
}

save_raster_multi <- function(png_path, stem, w_in = NULL, h_in = NULL, dpi = 300) {
  stopifnot(file.exists(png_path))
  img <- png::readPNG(png_path)
  nr <- dim(img)[1]
  nc <- dim(img)[2]
  if (is.null(w_in) || is.null(h_in)) {
    w_in <- max(nc / dpi, 4)
    h_in <- max(nr / dpi, 3)
  }
  grob <- grid::rasterGrob(img, interpolate = TRUE)

  file.copy(png_path, paste0(stem, ".png"), overwrite = TRUE)
  message("[OK] ", basename(paste0(stem, ".png")))

  grDevices::jpeg(paste0(stem, ".jpg"), width = w_in, height = h_in, units = "in",
                  res = dpi, quality = 95, type = "cairo")
  grid::grid.newpage(); grid::grid.draw(grob)
  invisible(grDevices::dev.off())
  message("[OK] ", basename(paste0(stem, ".jpg")))

  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w_in, height = h_in, family = FONT)
  grid::grid.newpage(); grid::grid.draw(grob)
  invisible(grDevices::dev.off())
  message("[OK] ", basename(paste0(stem, ".pdf")))

  b64 <- encode_b64(png_path)
  svg <- sprintf(
    paste0(
      '<?xml version="1.0" encoding="UTF-8"?>\n',
      '<svg xmlns="http://www.w3.org/2000/svg" ',
      'xmlns:xlink="http://www.w3.org/1999/xlink" ',
      'width="%.2fin" height="%.2fin" viewBox="0 0 %d %d">\n',
      '  <image width="%d" height="%d" xlink:href="data:image/png;base64,%s"/>\n',
      '</svg>\n'
    ),
    w_in, h_in, nc, nr, nc, nr, b64
  )
  writeLines(svg, paste0(stem, ".svg"), useBytes = TRUE)
  message("[OK] ", basename(paste0(stem, ".svg")))
  invisible(list(w = w_in, h = h_in, ncol = nc, nrow = nr))
}

## ---------- [1] HUMAnN Robust_DA KO nodes on map00650 -----------------------
message("[1] Build map00650 HUMAnN Robust_DA KO node table ...")

ko2pw <- read_tsv(
  file.path(comp_dir, "00_data/kegg_ko_pathway.tab"),
  col_names = c("path", "ko"),
  show_col_types = FALSE
) %>%
  transmute(
    pathway_id = str_replace(path, "^path:", ""),
    KO = str_replace(ko, "^ko:", "")
  ) %>%
  mutate(pathway_id = str_replace(pathway_id, "^ko", "map")) %>%
  filter(pathway_id == PATHWAY_ID) %>%
  distinct()

humann_ko <- read_tsv(
  file.path(humann_dir, "v2/Figure3/tables/Table_differential_KO_full.tsv"),
  show_col_types = FALSE
) %>%
  mutate(
    Robust_DA = as.logical(Robust_DA),
    KO = as.character(ko_id),
    KO_name = as.character(function_name)
  ) %>%
  filter(
    Robust_DA,
    KO %in% ko2pw$KO,
    !str_detect(KO_name, regex("unclassified|NO_NAME", ignore_case = TRUE))
  ) %>%
  mutate(
    pathway_id = PATHWAY_ID,
    pathway_name = PATHWAY_NAME,
    omics = "KO",
    source = "HUMAnN",
    color_value = as.numeric(log2FC),
    padj = pmin(Wilcoxon_FDR, CLRMC_FDR, na.rm = TRUE),
    # prefer current KEGG NAME; fall back to HUMAnN function_name
    KO_name = ifelse(KO %in% names(KEGG_KO_NAME), unname(KEGG_KO_NAME[KO]), KO_name)
  ) %>%
  arrange(desc(log2FC))

ko_nodes <- humann_ko
write_tsv(ko_nodes, file.path(out_tab, "map00650_KO_nodes.tsv"))
write_tsv(ko_nodes, file.path(out_data, "map00650_KO_nodes.tsv"))
write_tsv(ko_nodes, file.path(out_fig, "Fig_map00650_KO_lollipop.tsv"))
message(sprintf("  HUMAnN Robust_DA KOs on %s: %d", PATHWAY_ID, nrow(ko_nodes)))

# archive MAG Fisher presence hits for reference (not used in panel B)
fisher_arch <- read_tsv(
  file.path(mod03, "tables/fisher_KO_significant.tsv"),
  show_col_types = FALSE
) %>%
  mutate(
    KO_full = KO,
    KO = parse_ko_id(KO),
    KO_name = ifelse(str_detect(KO_full, ":"), sub("^[^:]+:", "", KO_full), KO)
  ) %>%
  filter(KO %in% ko2pw$KO) %>%
  mutate(source = "MAG_Fisher_presence")
write_tsv(fisher_arch, file.path(out_tab, "map00650_KO_MAG_Fisher_archive.tsv"))

## ---------- [2] Metabolite nodes on map00650 --------------------------------
message("[2] Build map00650 metabolite node table ...")

cpd_map <- read_tsv(
  file.path(mod07, "data/kegg_cpd_pathway_map_filtered.tsv"),
  show_col_types = FALSE
) %>%
  filter(path == PATHWAY_ID) %>%
  distinct(cpd)

metab <- read_tsv(
  file.path(mod07, "tables/metab_RobustDA_with_KEGG_CompoundID.tsv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(KEGG_Compound_ID), nzchar(KEGG_Compound_ID)) %>%
  filter(KEGG_Compound_ID %in% cpd_map$cpd)

# aggregate duplicate Compound IDs (e.g. C00025)
metab_agg <- metab %>%
  group_by(KEGG_Compound_ID) %>%
  summarise(
    metab_ids = paste(unique(metab_id), collapse = ";"),
    Metabolite = paste(unique(Metabolite), collapse = " / "),
    Robust_direction = {
      tab <- table(Robust_direction)
      names(tab)[which.max(tab)]
    },
    log2FC = mean(log2FC, na.rm = TRUE),
    Wilcoxon_FDR = min(Wilcoxon_FDR, na.rm = TRUE),
    CLRMC_FDR = min(CLRMC_FDR, na.rm = TRUE),
    n_features = n(),
    .groups = "drop"
  ) %>%
  mutate(
    pathway_id = PATHWAY_ID,
    pathway_name = PATHWAY_NAME,
    omics = "Metabolite",
    PP_direction = ifelse(
      Robust_direction == "PoorPlaque", "Up in PoorPlaque", "Down in PoorPlaque"
    ),
    color_value = log2FC,
    Metabolite = ifelse(
      KEGG_Compound_ID %in% names(KEGG_CPD_NAME),
      unname(KEGG_CPD_NAME[KEGG_Compound_ID]),
      Metabolite
    )
  ) %>%
  arrange(desc(log2FC))

write_tsv(metab, file.path(out_tab, "map00650_metab_features_raw.tsv"))
write_tsv(metab_agg, file.path(out_tab, "map00650_metab_nodes.tsv"))
write_tsv(metab_agg, file.path(out_data, "map00650_metab_nodes.tsv"))
write_tsv(metab_agg, file.path(out_fig, "Fig_map00650_metab_lollipop.tsv"))
message(sprintf("  Robust_DA compounds on %s: %d (from %d features)",
                PATHWAY_ID, nrow(metab_agg), nrow(metab)))

## ---------- [3] Enrichment summary + unified node table ---------------------
message("[3] Assemble enrichment summary & unified node table ...")

# HUMAnN pathway enrichment (ko00650)
humann_enrich_files <- c(
  PoorPlaque = file.path(
    humann_dir, "v2/Figure3/code/Fig3_enrich_results/Fig3_pathway_enrichment_PoorPlaque_up.tsv"
  ),
  Normal = file.path(
    humann_dir, "v2/Figure3/code/Fig3_enrich_results/Fig3_pathway_enrichment_Normal_up.tsv"
  )
)
humann_enrich <- map_dfr(names(humann_enrich_files), function(grp) {
  f <- humann_enrich_files[[grp]]
  if (!file.exists(f)) return(tibble())
  read_tsv(f, show_col_types = FALSE) %>%
    filter(ID == "ko00650" | ID == PATHWAY_ID) %>%
    mutate(
      omics = "KO_HUMAnN",
      enrich_group = grp,
      ID = PATHWAY_ID,
      sig_p05 = pvalue < 0.05,
      sig_padj05 = p.adjust < 0.05
    )
})

ko_enrich_mag <- read_tsv(
  file.path(mod12, "tables/KEGG_enrichment_all_sets.tsv"),
  show_col_types = FALSE
) %>%
  filter(ID == PATHWAY_ID) %>%
  mutate(omics = "KO_MAG_Fisher")

metab_enrich <- read_tsv(
  file.path(mod07, "tables/Fig_metab_KEGG_enrichment_all_sets.tsv"),
  show_col_types = FALSE
) %>%
  filter(ID == PATHWAY_ID) %>%
  mutate(
    omics = "Metabolite",
    enrich_group = .data[["set"]],
    sig_p05 = pvalue < 0.05,
    sig_padj05 = if ("sig_padj05" %in% names(.)) sig_padj05 else p.adjust < 0.05
  )

pick_enrich_cols <- function(df) {
  if (nrow(df) == 0L) return(tibble())
  df %>%
    transmute(
      omics, enrich_group,
      ID, Description,
      GeneRatio, BgRatio, Count,
      pvalue, p.adjust,
      qvalue = if ("qvalue" %in% names(df)) qvalue else NA_real_,
      geneID, sig_p05, sig_padj05
    )
}

enrich_sum <- bind_rows(
  pick_enrich_cols(humann_enrich),
  pick_enrich_cols(ko_enrich_mag),
  pick_enrich_cols(metab_enrich)
)

write_tsv(enrich_sum, file.path(out_tab, "map00650_enrichment_summary.tsv"))
write_tsv(enrich_sum, file.path(out_data, "map00650_enrichment_summary.tsv"))

node_summary <- bind_rows(
  ko_nodes %>%
    transmute(
      pathway_id, pathway_name, omics = "KO_HUMAnN",
      node_id = KO, node_name = KO_name,
      direction = Robust_direction,
      score = color_value,
      padj,
      detail = sprintf("log2FC=%.3f; Wilcoxon_FDR=%.3g; CLRMC_FDR=%.3g",
                       log2FC, Wilcoxon_FDR, CLRMC_FDR)
    ),
  metab_agg %>%
    transmute(
      pathway_id, pathway_name, omics,
      node_id = KEGG_Compound_ID, node_name = Metabolite,
      direction = Robust_direction,
      score = color_value,
      padj = pmin(Wilcoxon_FDR, CLRMC_FDR, na.rm = TRUE),
      detail = sprintf("log2FC=%.3f; n_features=%d", log2FC, n_features)
    )
)
write_tsv(node_summary, file.path(out_tab, "map00650_node_summary.tsv"))
write_tsv(node_summary, file.path(out_data, "map00650_node_summary.tsv"))
write_tsv(node_summary, file.path(out_fig, "Fig_map00650_composite.tsv"))

## ---------- [4] pathview overlay --------------------------------------------
message("[4] Render pathview KO+metabolite overlay ...")

gene.data <- setNames(ko_nodes$color_value, ko_nodes$KO)
cpd.data  <- setNames(metab_agg$color_value, metab_agg$KEGG_Compound_ID)
gene_lim  <- max(1, ceiling(max(abs(gene.data)) * 2) / 2)
cpd_lim   <- max(0.3, ceiling(max(abs(cpd.data)) * 10) / 10)

# pathview writes into getwd(); use cache then copy
old_wd <- getwd()
setwd(cache_dir)
pv <- tryCatch(
  pathview(
    gene.data   = gene.data,
    cpd.data    = cpd.data,
    pathway.id  = PATHWAY_NUM,
    species     = "ko",
    gene.idtype = "KEGG",
    cpd.idtype  = "kegg",
    keys.align  = "y",
    limit       = list(gene = gene_lim, cpd = cpd_lim),
    bins        = list(gene = 10, cpd = 10),
    low         = list(gene = COL_NM, cpd = COL_NM),
    mid         = list(gene = "grey90", cpd = "grey90"),
    high        = list(gene = COL_PP, cpd = COL_PP),
    na.col      = "white",
    same.layer  = FALSE,
    kegg.native = TRUE,
    out.suffix  = "KO_metab"
  ),
  error = function(e) {
    message("[ERR] pathview failed: ", e$message)
    NULL
  }
)
setwd(old_wd)

pv_png <- file.path(cache_dir, "ko00650.KO_metab.png")
if (!file.exists(pv_png)) {
  # sometimes naming differs
  cand <- list.files(cache_dir, pattern = "00650.*KO_metab.*\\.png$", full.names = TRUE)
  if (length(cand)) pv_png <- cand[[1]]
}
stopifnot(file.exists(pv_png))

# export pathview plot data if present
if (!is.null(pv) && is.list(pv)) {
  plot_df <- tryCatch({
    if (!is.null(pv$plot.data.gene) || !is.null(pv$plot.data.cpd)) {
      bind_rows(
        if (!is.null(pv$plot.data.gene)) {
          as_tibble(pv$plot.data.gene, rownames = "row_id") %>%
            mutate(layer = "gene")
        } else NULL,
        if (!is.null(pv$plot.data.cpd)) {
          as_tibble(pv$plot.data.cpd, rownames = "row_id") %>%
            mutate(layer = "compound")
        } else NULL
      )
    } else NULL
  }, error = function(e) NULL)
  if (!is.null(plot_df) && nrow(plot_df) > 0) {
    write_tsv(plot_df, file.path(out_data, "map00650_pathview_plot_data.tsv"))
    write_tsv(plot_df, file.path(out_tab, "map00650_pathview_plot_data.tsv"))
  }
}

stem_pv <- file.path(out_fig, "Fig_map00650_pathview")
pv_dim <- save_raster_multi(pv_png, stem_pv, dpi = 300)
# plotting TSV for pathview panel = node color inputs
write_tsv(
  bind_rows(
    tibble(node_id = names(gene.data), omics = "KO", color_value = as.numeric(gene.data)),
    tibble(node_id = names(cpd.data), omics = "Metabolite", color_value = as.numeric(cpd.data))
  ),
  paste0(stem_pv, ".tsv")
)
file.copy(file.path(cache_dir, "ko00650.xml"),
          file.path(out_data, "ko00650.xml"), overwrite = TRUE)

## ---------- [5] Companion: HUMAnN KO log2FC lollipop (panel B) ---------------
message("[5] Plot HUMAnN KO log2FC lollipop ...")

# Single-line axis: "KO  short-name" (full KEGG NAME kept in TSV)
KO_NAME_SHORT <- c(
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
  K00246 = "SDH subunit C"
)

ko_plot <- ko_nodes %>%
  mutate(
    name_short = unname(KO_NAME_SHORT[KO]),
    name_short = ifelse(
      is.na(name_short),
      str_trunc(
        str_trim(str_replace(KO_name, "\\s*\\[EC:[^\\]]+\\]\\s*$", "")),
        28L, side = "right", ellipsis = "..."
      ),
      name_short
    ),
    label = paste0(KO, "  ", name_short),
    label = factor(label, levels = label[order(log2FC)])
  )

p_ko <- ggplot(ko_plot, aes(x = log2FC, y = label, colour = Robust_direction)) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_segment(
    aes(x = 0, xend = log2FC, yend = label),
    linewidth = 0.7
  ) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = GROUP_COLORS, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.10))) +
  labs(
    x = expression(log[2] * "FC (PoorPlaque / Normal)"),
    y = NULL,
    title = "HUMAnN Robust_DA KOs on map00650"
  ) +
  theme(
    legend.position = "top",
    legend.margin = margin(0, 0, -2, 0),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = BASE_SIZE, colour = TXT_COL),
    plot.margin = margin(2, 4, 2, 2),
    plot.title = element_text(size = BASE_SIZE, colour = TXT_COL),
    axis.text = element_text(size = BASE_SIZE, colour = TXT_COL),
    axis.title = element_text(size = BASE_SIZE, colour = TXT_COL),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3)
  )

stem_ko <- file.path(out_fig, "Fig_map00650_KO_lollipop")
save_gg(p_ko, stem_ko, w = FIG_IN * 0.52, h = FIG_IN * 0.48, dpi = 300)
write_tsv(
  ko_plot %>%
    select(KO, KO_name, name_short, Robust_direction, log2FC, Wilcoxon_FDR, CLRMC_FDR,
           mean_Normal, mean_PoorPlaque, padj),
  paste0(stem_ko, ".tsv")
)
write_tsv(
  ko_plot %>%
    select(KO, KO_name, name_short, Robust_direction, log2FC, Wilcoxon_FDR, CLRMC_FDR,
           mean_Normal, mean_PoorPlaque, padj),
  file.path(out_data, "Fig_map00650_KO_lollipop.tsv")
)

## ---------- [6] Companion: metabolite log2FC lollipop ------------------------
message("[6] Plot metabolite log2FC panel ...")

metab_plot <- metab_agg %>%
  mutate(
    label = paste0(KEGG_Compound_ID, "  ", Metabolite),
    label = factor(label, levels = label[order(log2FC)])
  )

p_metab <- ggplot(metab_plot, aes(x = log2FC, y = label, colour = Robust_direction)) +
  geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.4) +
  geom_segment(
    aes(x = 0, xend = log2FC, yend = label),
    linewidth = 0.7
  ) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = GROUP_COLORS, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.10))) +
  labs(
    x = expression(log[2] * "FC (PoorPlaque / Normal)"),
    y = NULL,
    title = "Robust_DA metabolites on map00650"
  ) +
  theme(
    legend.position = "top",
    legend.margin = margin(0, 0, -2, 0),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = BASE_SIZE, colour = TXT_COL),
    plot.margin = margin(2, 4, 2, 2),
    plot.title = element_text(size = BASE_SIZE, colour = TXT_COL),
    axis.text = element_text(size = BASE_SIZE, colour = TXT_COL),
    axis.title = element_text(size = BASE_SIZE, colour = TXT_COL),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3)
  )

stem_metab <- file.path(out_fig, "Fig_map00650_metab_lollipop")
save_gg(p_metab, stem_metab, w = FIG_IN * 0.5, h = FIG_IN * 0.42, dpi = 300)
# TSV already written above; refresh from plot order
write_tsv(
  metab_plot %>%
    select(KEGG_Compound_ID, Metabolite, Robust_direction, log2FC,
           Wilcoxon_FDR, CLRMC_FDR, n_features, metab_ids),
  paste0(stem_metab, ".tsv")
)

## ---------- [7] Composite ----------------------------------------------------
message("[7] Build composite figure ...")

img <- png::readPNG(pv_png)

# B/C for bottom row (shared legend on B); unified black text size
p_ko_b <- p_ko +
  labs(title = NULL) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(0, "pt"),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = BASE_SIZE, colour = TXT_COL),
    plot.margin = margin(0, 2, 0, 0),
    axis.text = element_text(size = BASE_SIZE, colour = TXT_COL),
    axis.title = element_text(size = BASE_SIZE, colour = TXT_COL, margin = margin(2, 0, 0, 0))
  )
p_metab_c <- p_metab +
  labs(title = NULL) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 2),
    axis.text = element_text(size = BASE_SIZE, colour = TXT_COL),
    axis.title = element_text(size = BASE_SIZE, colour = TXT_COL, margin = margin(2, 0, 0, 0))
  )

# Grid assembly: A full-bleed top; B|C bottom — reliable vs patchwork letterboxing
draw_composite <- function() {
  grid::grid.newpage()
  txt_gp <- grid::gpar(fontsize = BASE_SIZE, fontfamily = FONT, fontface = "plain", col = TXT_COL)
  # compact title strip (same size, all black)
  grid::pushViewport(grid::viewport(x = 0.5, y = 0.978, width = 0.98, height = 0.03))
  grid::grid.text("Butanoate metabolism (map00650 / ko00650)", gp = txt_gp)
  grid::popViewport()
  grid::pushViewport(grid::viewport(x = 0.5, y = 0.952, width = 0.98, height = 0.022))
  grid::grid.text(
    "A: pathview (HUMAnN KO + metabolite log2FC); B/C: node log2FC (KEGG NAME)",
    gp = txt_gp
  )
  grid::popViewport()

  # content area under title (tight; no inter-panel gap)
  content_bottom <- 0.012
  content_top <- 0.938
  content_h <- content_top - content_bottom
  a_h <- COMP_HEIGHTS[1] / sum(COMP_HEIGHTS) * content_h
  bc_h <- COMP_HEIGHTS[2] / sum(COMP_HEIGHTS) * content_h
  a_y <- content_bottom + bc_h + a_h / 2
  bc_y <- content_bottom + bc_h / 2

  grid::pushViewport(grid::viewport(x = 0.5, y = a_y, width = 0.99, height = a_h))
  grid::grid.raster(img, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
  grid::grid.text(
    "A", x = unit(0.008, "npc"), y = unit(0.99, "npc"),
    just = c("left", "top"), gp = txt_gp
  )
  grid::popViewport()

  grid::pushViewport(grid::viewport(x = 0.5, y = bc_y, width = 0.99, height = bc_h))
  grid::grid.draw(gridExtra::arrangeGrob(
    gridExtra::arrangeGrob(
      ggplot2::ggplotGrob(p_ko_b),
      top = grid::textGrob("B", x = 0.01, hjust = 0, gp = txt_gp)
    ),
    gridExtra::arrangeGrob(
      ggplot2::ggplotGrob(p_metab_c),
      top = grid::textGrob("C", x = 0.01, hjust = 0, gp = txt_gp)
    ),
    ncol = 2, widths = c(1.3, 1)
  ))
  grid::popViewport()
}

save_composite_multi <- function(stem, w_in, h_in, dpi = 300) {
  # PNG
  grDevices::png(paste0(stem, ".png"), width = w_in, height = h_in, units = "in",
                 res = dpi, type = "cairo", bg = "white")
  draw_composite()
  invisible(grDevices::dev.off())
  message("[OK] ", basename(paste0(stem, ".png")))
  # JPG
  grDevices::jpeg(paste0(stem, ".jpg"), width = w_in, height = h_in, units = "in",
                  res = dpi, quality = 95, type = "cairo", bg = "white")
  draw_composite()
  invisible(grDevices::dev.off())
  message("[OK] ", basename(paste0(stem, ".jpg")))
  # PDF
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w_in, height = h_in, family = FONT, bg = "white")
  draw_composite()
  invisible(grDevices::dev.off())
  message("[OK] ", basename(paste0(stem, ".pdf")))
  # SVG via Cairo
  if (isTRUE(capabilities("cairo"))) {
    grDevices::svg(paste0(stem, ".svg"), width = w_in, height = h_in, family = FONT, bg = "white")
    draw_composite()
    invisible(grDevices::dev.off())
    message("[OK] ", basename(paste0(stem, ".svg")))
  }
}

if (!requireNamespace("gridExtra", quietly = TRUE)) {
  stop("Package gridExtra is required for composite assembly")
}

stem_comp <- file.path(out_fig, "Fig_map00650_composite")
save_composite_multi(stem_comp, w_in = FIG_IN, h_in = FIG_IN, dpi = 300)
write_tsv(node_summary, paste0(stem_comp, ".tsv"))
write_tsv(node_summary, file.path(out_data, "Fig_map00650_composite.tsv"))

## also export B|C as its own subfigure
p_side_only <- p_ko_b + p_metab_c +
  plot_layout(nrow = 1, widths = c(1.25, 1), guides = "collect") +
  plot_annotation(
    title = "map00650 node-level log2FC (KEGG NAME)",
    tag_levels = list(c("B", "C")),
    theme = theme(
      plot.title = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL),
      plot.tag = element_text(family = FONT, face = "plain", size = BASE_SIZE, colour = TXT_COL),
      plot.margin = margin(2, 3, 2, 3)
    )
  ) &
  theme(
    legend.position = "bottom",
    text = element_text(size = BASE_SIZE, colour = TXT_COL)
  )
stem_side <- file.path(out_fig, "Fig_map00650_node_panels")
save_gg(p_side_only, stem_side, w = FIG_IN, h = FIG_IN * 0.48, dpi = 300)
write_tsv(node_summary, paste0(stem_side, ".tsv"))

## ---------- [8] README -------------------------------------------------------
readme <- c(
  "# 15_map00650_KO_metab_pathway",
  "",
  "Cross-omics visualization of **Butanoate metabolism** (`map00650` / `ko00650`).",
  "",
  "## Claim",
  "HUMAnN community KO abundance and metabolome converge on butanoate metabolism:",
  "PoorPlaque-up KOs (e.g. 3-hydroxybutyryl-CoA dehydrogenase, short-chain-enoyl-CoA hydratase)",
  "and PoorPlaque-up metabolites (2-Oxoglutarate, Fumarate, Succinate, 4-Aminobutanoate, L-Glutamate, Pyruvate),",
  "with Normal-enriched succinate dehydrogenase (frd) subunits / gabT (K07250) and Butanal;",
  "note puuE (K00823) is PoorPlaque-up (opposite to gabT).",
  "",
  "## Panels",
  "- `Fig_map00650_pathview` — native KEGG map; KO boxes = HUMAnN log2FC; circles = metabolite log2FC",
  "- `Fig_map00650_KO_lollipop` — HUMAnN Robust_DA KO log2FC on the pathway (panel B)",
  "- `Fig_map00650_metab_lollipop` — Robust_DA metabolite log2FC on the pathway (panel C)",
  "- `Fig_map00650_composite` — A (top, larger) pathview + B|C (bottom row) node panels; 170×170 mm",
  "- `Fig_map00650_node_panels` — B|C only",
  "",
  "## Inputs",
  "- HUMAnN: `M2-humann/v2/Figure3/tables/Table_differential_KO_full.tsv`",
  "- HUMAnN enrichment: `Fig3_pathway_enrichment_*_up.tsv`",
  "- Module 07: metabolite Robust_DA + Compound IDs; cpd→pathway map",
  "- Module 12: MAG Fisher KO enrichment (archived in enrichment summary)",
  "- `compare_genomics_all/00_data/kegg_ko_pathway.tab`",
  "",
  "## Notes",
  "- Panel B uses HUMAnN community KO abundance log2FC (PoorPlaque/Normal), Robust_DA only.",
  "- Node labels use current KEGG NAME (REST, 2026-08-08); axis text omits trailing [EC:...].",
  "- Pathview KO colouring matches panel B (continuous log2FC).",
  "- Duplicate Compound IDs (e.g. C00025) are mean-aggregated for pathview.",
  "- Colours: Normal `#00A087`; PoorPlaque `#E64B35`.",
  "- Font: Times New Roman, plain; composite canvas 170 mm × 170 mm.",
  "",
  sprintf("- pathview KO colour limit: ±%.2f; compound limit: ±%.2f", gene_lim, cpd_lim),
  sprintf("- n HUMAnN KO nodes: %d; n metabolite nodes: %d", nrow(ko_nodes), nrow(metab_agg))
)
writeLines(readme, file.path(module_dir, "README.md"))

message("[DONE] Module 15 map00650 KO+metabolite pathway figure")
message("  figures: ", out_fig)
