## Packaged for revise/; source: project/review/R1.4/scripts/17_FigS_bio8_KOmetab_lm_scatter.R
## Figures: FigS3
## Note: curated high-|rho| LM scatters
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "FigS3_KOmetab_lm_scatter.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# Fig.S: 8 biologically curated KO–metabolite LM scatters (2×4).
# X-axis: log10(KO+ε) scale but label without "log10" wording; Y: full metabolite name.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(svglite)
})

ROOT <- file.path(.pkg_root, ""
BASE <- file.path(ROOT, "01this_study_dataset/sgb/metab_corr")
HUMANN <- file.path(ROOT, "01this_study_dataset/M2-humann")
OUT <- file.path(ROOT, "review/R1.4/Fig2D_KOmetab_noPI")
TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")
FIG_SUB <- file.path(FIG, "FigS_top10_KOmetab_lm_scatter_panels")
dir.create(FIG_SUB, recursive = TRUE, showWarnings = FALSE)

COL <- c(Normal = "#1F6F78", PoorPlaque = "#A8483E")
CI_FILL <- "#9AA4AD"
LINE_COL <- "#1A1A1A"
EPS <- 0.1
FONT <- if (any(grepl("^Liberation Serif$", system("fc-list : family", intern = TRUE)))) {
  "Liberation Serif"
} else if (any(grepl("Times New Roman", system("fc-list : family", intern = TRUE), fixed = TRUE))) {
  "Times New Roman"
} else {
  "Nimbus Roman"
}
BASE_PT <- 16  # was 10; +60%

dual_write <- function(df, stem, dir = TAB) {
  write_tsv(df, file.path(dir, paste0(stem, ".tsv")), na = "")
  write_csv(df, file.path(dir, paste0(stem, ".csv")), na = "")
}

save_pub <- function(plot, stem, dir, width_in, height_in, dpi = 300) {
  sysf <- list(sans = FONT, serif = FONT, mono = FONT)
  svglite::svglite(file.path(dir, paste0(stem, ".svg")), width = width_in, height = height_in, system_fonts = sysf)
  print(plot); invisible(dev.off())
  grDevices::cairo_pdf(file.path(dir, paste0(stem, ".pdf")), width = width_in, height = height_in, family = FONT)
  print(plot); invisible(dev.off())
  png(file.path(dir, paste0(stem, ".png")), width = width_in, height = height_in, units = "in",
      res = dpi, type = "cairo", family = FONT)
  print(plot); invisible(dev.off())
  jpeg(file.path(dir, paste0(stem, ".jpg")), width = width_in, height = height_in, units = "in",
       res = dpi, quality = 95, type = "cairo", family = FONT)
  print(plot); invisible(dev.off())
}

theme_pub <- function(base_size = BASE_PT) {
  # use + (not %+replace%) so title/margin inheritance is not wiped
  theme_classic(base_size = base_size, base_family = FONT) +
    theme(
      text = element_text(family = FONT, colour = "black", face = "plain"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      axis.title = element_text(family = FONT, colour = "black", size = base_size - 1),
      axis.title.x = element_text(margin = margin(t = 6), size = base_size - 2, lineheight = 1.05),
      axis.title.y = element_text(margin = margin(r = 6), angle = 90, vjust = 0.5,
                                 size = base_size - 2, lineheight = 1.05),
      axis.text = element_text(family = FONT, colour = "black", size = base_size - 4),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.3),
      legend.title = element_blank(),
      legend.text = element_text(family = FONT, size = base_size - 2),
      legend.position = "bottom",
      legend.margin = margin(2, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.key = element_blank(),
      legend.key.size = unit(0.4, "cm"),
      plot.margin = margin(12, 14, 8, 14),
      plot.tag = element_text(family = FONT, colour = "black", size = base_size + 1),
      plot.tag.position = c(0.01, 0.99)
    )
}

short_gene <- function(x, n = 18) {
  x <- gsub("\\s+", " ", as.character(x))
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "…"), x)
}

wrap_ylab <- function(x, width = 16) {
  stringr::str_wrap(as.character(x), width = width)
}

# Curated 8 pairs (order = panel A–H)
pairs_sel <- tribble(
  ~KO, ~metabolite,
  "K15371", "metab_28000",
  "K15371", "metab_17260",
  "K10006", "metab_542",
  "K05989", "metab_542",
  "K15923", "metab_19574",
  "K15923", "metab_28000",
  "K02257", "metab_16879",
  "K02257", "metab_16880"
) %>%
  mutate(panel = LETTERS[row_number()], rank = row_number())

message("font: ", FONT)
message("[1/4] Load edges + abundance ...")
group <- read_tsv(file.path(BASE, "this.study.group.tsv"),
                  col_names = c("sample", "group"), show_col_types = FALSE) %>%
  mutate(group = factor(group, levels = c("Normal", "PoorPlaque")))
samples <- group$sample

edges_all <- read_tsv(file.path(TAB, "Table_Fig2D_KOmetab_edges_noPI.tsv"), show_col_types = FALSE) %>%
  mutate(
    abs_rho = abs(rho),
    KO = ifelse(str_starts(from, "K"), from, to),
    metabolite = ifelse(str_starts(from, "K"), to, from),
    metab_lab = ifelse(str_starts(from, "K"), as.character(label_to), as.character(label_from))
  )

edges <- pairs_sel %>%
  left_join(
    edges_all %>% select(KO, metabolite, rho, p, FDR, edge_sign, enrich_from, enrich_to, metab_lab, abs_rho),
    by = c("KO", "metabolite")
  )
stopifnot(!any(is.na(edges$rho)))

diff_ko_path <- file.path(HUMANN, "v2/Figure3/tables/Table_differential_KO_full.tsv")
ko_ann <- read_tsv(diff_ko_path, show_col_types = FALSE) %>%
  transmute(KO = ko_id, function_name) %>% distinct(KO, .keep_all = TRUE)
edges <- edges %>%
  left_join(ko_ann, by = "KO") %>%
  mutate(
    ko_func = ifelse(is.na(function_name), KO, as.character(function_name)),
    ko_lab = KO,
    metab_short = metab_lab,  # keep full name for tables; axis uses full name
    prevalence_filtered = TRUE
  )

dual_write(edges %>% select(
  rank, panel, KO, metabolite, ko_lab, metab_lab, metab_short, ko_func,
  rho, abs_rho, p, FDR, edge_sign, enrich_from, enrich_to, prevalence_filtered
), "Table_FigS_top10_KOmetab_pairs")

metab_path <- file.path(BASE, "output/tables/metab_abund_AB_select.tsv")
metab <- read_tsv(metab_path, show_col_types = FALSE) %>% as.data.frame()
rownames(metab) <- metab$metab_id
metab$metab_id <- NULL
metab <- metab[, samples, drop = FALSE]

ko_file <- file.path(HUMANN, "comb_genefamilies-cpm.uniref90_ko.normed.named.tsv")
ko_raw <- read_tsv(ko_file, show_col_types = FALSE, comment = "")
feat_col <- names(ko_raw)[1]
col_map <- tibble(col = setdiff(names(ko_raw), feat_col)) %>%
  mutate(sample = str_extract(col, "^[AB]\\d+")) %>%
  filter(sample %in% samples)

needed_ko <- unique(edges$KO)
ko_ab <- ko_raw %>%
  rename(feature = !!feat_col) %>%
  filter(!str_detect(feature, "\\|")) %>%
  mutate(ko_id = str_extract(feature, "K\\d{5}")) %>%
  filter(ko_id %in% needed_ko) %>%
  distinct(ko_id, .keep_all = TRUE) %>%
  select(ko_id, all_of(col_map$col)) %>%
  as.data.frame()
rownames(ko_ab) <- ko_ab$ko_id
ko_ab$ko_id <- NULL
colnames(ko_ab) <- col_map$sample[match(colnames(ko_ab), col_map$col)]
ko_ab <- ko_ab[, samples, drop = FALSE]
stopifnot(all(needed_ko %in% rownames(ko_ab)))
stopifnot(all(edges$metabolite %in% rownames(metab)))

message("[2/4] Pairwise LM ...")
pair_list <- list()
stat_rows <- list()
for (i in seq_len(nrow(edges))) {
  r <- edges[i, ]
  df <- tibble(
    sample = samples,
    group = group$group[match(samples, group$sample)],
    KO_ab = as.numeric(ko_ab[r$KO, samples]),
    Metab_ab = as.numeric(metab[r$metabolite, samples]),
    panel = r$panel,
    rank = r$rank,
    KO = r$KO,
    metabolite = r$metabolite,
    metab_lab = r$metab_lab,
    metab_short = r$metab_lab,
    spearman_rho = r$rho,
    spearman_FDR = r$FDR,
    prevalence_filtered = TRUE
  ) %>%
    mutate(KO_ab_log10 = log10(KO_ab + EPS))

  fit <- lm(Metab_ab ~ KO_ab_log10, data = df)
  s <- summary(fit)
  ct <- cor.test(df$KO_ab, df$Metab_ab, method = "spearman", exact = FALSE)
  pear <- cor.test(df$KO_ab, df$Metab_ab, method = "pearson")
  fit0 <- lm(Metab_ab ~ KO_ab, data = df)
  s0 <- summary(fit0)
  gene_short <- short_gene(r$ko_func)

  stat_rows[[i]] <- tibble(
    rank = r$rank, panel = r$panel, KO = r$KO, metabolite = r$metabolite,
    metab_lab = r$metab_lab, metab_short = r$metab_lab,
    ko_func = r$ko_func, gene_short = gene_short, n = nrow(df),
    network_rho = r$rho, network_FDR = r$FDR,
    spearman_rho = unname(ct$estimate), spearman_p = ct$p.value,
    pearson_r = unname(pear$estimate), pearson_p = pear$p.value,
    lm_raw_r2 = s0$r.squared, lm_raw_slope_p = s0$coefficients[2, 4],
    lm_log_r2 = s$r.squared, lm_log_slope = unname(coef(fit)[2]),
    lm_log_intercept = unname(coef(fit)[1]), lm_log_slope_p = s$coefficients[2, 4],
    prevalence_filtered = TRUE, log_eps = EPS,
    # no "log10" wording on axis; tick labels already show raw KO abundance
    # KO + gene on x; stats go to caption (outside data area) to avoid occlusion
    xlab = sprintf("%s\n(%s)", r$KO, gene_short),
    ylab = wrap_ylab(r$metab_lab, width = 14),
    stats_lab = sprintf("ρ = %+.3f · R² = %.2f", r$rho, s$r.squared)
  )
  pair_list[[i]] <- df
  dual_write(df, sprintf("FigS_top10_KOmetab_lm_scatter_%s_plotdata", r$panel))
}

plot_long <- bind_rows(pair_list)
stats_df <- bind_rows(stat_rows)
dual_write(plot_long, "FigS_top10_KOmetab_lm_scatter_panel_plotdata")
dual_write(stats_df, "Table_FigS_top10_KOmetab_lm_stats")

log_breaks_for <- function(x) {
  rng <- range(x, na.rm = TRUE)
  raw_cand <- c(0, 0.1, 0.3, 1, 3, 10, 30, 100)
  br <- log10(raw_cand + EPS)
  br <- br[br >= rng[1] - 1e-6 & br <= rng[2] + 1e-6]
  if (length(br) < 3) pretty(rng, n = 4) else unique(br)
}
log_lab <- function(x) {
  raw <- pmax(10^x - EPS, 0)
  ifelse(abs(raw) < 1e-8, "0",
         ifelse(raw < 1, sprintf("%.1f", raw), sprintf("%.0f", raw)))
}

make_one <- function(df, st) {
  br <- log_breaks_for(df$KO_ab_log10)
  if (length(br) > 5) {
    br <- br[unique(round(seq(1, length(br), length.out = 5)))]
  }
  ggplot(df, aes(x = KO_ab_log10, y = Metab_ab)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, level = 0.95,
                colour = LINE_COL, fill = CI_FILL, alpha = 0.30, linewidth = 0.6) +
    geom_point(aes(fill = group), shape = 21, size = 2.1, stroke = 0.25,
               colour = "white", alpha = 0.95) +
    scale_fill_manual(values = COL, name = NULL) +
    scale_x_continuous(breaks = br, labels = log_lab, expand = expansion(mult = c(0.05, 0.06))) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.08))) +
    labs(
      title = NULL, subtitle = NULL,
      caption = st$stats_lab[1],
      x = st$xlab[1], y = st$ylab[1]
    ) +
    theme_pub() +
    theme(
      legend.position = "none",
      plot.caption = element_text(
        family = FONT, colour = "#222222", size = BASE_PT - 3,
        hjust = 0.5, margin = margin(t = 4, b = 0)
      )
    )
}

message("[3/4] Draw singles (stats as caption; no in-panel annotate) ...")
panels <- edges$panel
plots <- map(seq_along(panels), function(i) {
  pn <- panels[i]
  p <- make_one(pair_list[[i]], stats_df %>% filter(panel == pn)) + labs(tag = pn)
  save_pub(p, sprintf("FigS_top10_KOmetab_lm_scatter_%s", pn), FIG_SUB, 3.8, 3.7)
  p
})

message("[4/4] Draw 2×4 panel ...")
panel <- wrap_plots(plots, ncol = 4, nrow = 2) +
  plot_layout(guides = "collect", widths = rep(1, 4), heights = c(1, 1)) &
  theme(
    legend.position = "bottom",
    legend.text = element_text(family = FONT, size = 12),
    plot.margin = margin(10, 12, 6, 12)
  )
panel <- panel +
  plot_annotation(
    title = "Curated KO–metabolite associations: LM fit with 95% CI",
    subtitle = "Eight pairs with |ρ|≥0.9 and KO nonzero ≥80%; grey band = LM 95% CI",
    theme = theme(
      plot.title = element_text(family = FONT, size = 16, hjust = 0.5, colour = "black",
                                margin = margin(b = 2)),
      plot.subtitle = element_text(family = FONT, size = 11, hjust = 0.5, colour = "#555555",
                                   margin = margin(b = 8))
    )
  )
# taller/wider canvas: room for wrapped y, 2-line x, caption stats
save_pub(panel, "FigS_top10_KOmetab_lm_scatter_panel", FIG, 16.0, 9.2)

writeLines(c(
  "FigS bio8 curated KO–metabolite LM scatters",
  "Pairs: K15371–NAc-Tyr/NAc-Glu; K10006–Glu; K05989–Glu; K15923–NAc-Phe/NAc-Tyr; K02257–Fumarate/Malate",
  "X uses log10(KO+0.1) for fitting/display ticks show raw KO; axis title has no 'log10' wording",
  "Y-axis: full metabolite library names"
), file.path(FIG_SUB, "README_axis_note.txt"))

message("DONE")
print(stats_df %>% select(panel, KO, metabolite, metab_lab, network_rho, lm_log_r2, xlab))
message("  ", file.path(FIG, "FigS_top10_KOmetab_lm_scatter_panel.jpg"))
