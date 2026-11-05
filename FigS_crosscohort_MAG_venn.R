## Packaged for revise/; source: project/02dataset/input/dataset_venn_analysis.R
## Figures: FigS8B
## Note: MAG Venn across cohorts
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "FigS_crosscohort_MAG_venn.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

# ==============================================================================
#  Dataset-level dRep MAG species composition overlap analysis
#  -- Four CNS-level publication figures, unified design
#
#  Fig 1 : Venn diagram  (three-set SGB overlap)
#  Fig 2 : Alluvial plot (Venn region -> Phylum -> Genus)
#  Fig 3 : KEGG KO richness boxplot  (medium-/high-quality MAGs only)
#  Fig 4 : Region-enriched KO heatmap  (Fisher exact + FDR)
#
#  Usage
#  -----
#    Rscript dataset_venn_analysis.R \
#      --cdb       Cdb.csv \
#      --wdb       Wdb.csv \
#      --group     group.tsv \
#      --gtdb      gtdbtk_bac120_summary.tsv \
#      --ko        total_KO.xls \
#      --outdir    results_dataset_venn \
#      --datasets  "ThisStudy,PRJNA766357,PRJNA597465" \
#      --min_score 50
# ==============================================================================

suppressPackageStartupMessages(library(methods))

# ==============================================================================
#  0.  PARSE CLI ARGUMENTS
# ==============================================================================

parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  d <- list(
    cdb       = "Cdb.csv",
    wdb       = "Wdb.csv",
    group     = "group.tsv",
    gtdb      = "gtdbtk_bac120_summary.tsv",
    ko        = "total_KO.xls",
    outdir    = "results_dataset_venn",
    datasets  = "ThisStudy,PRJNA766357,PRJNA597465",
    min_score = "50"
  )
  i <- 1
  while (i <= length(a)) {
    k <- sub("^--", "", a[i])
    if (k %in% names(d) && i + 1 <= length(a)) { d[[k]] <- a[i + 1]; i <- i + 2 }
    else i <- i + 1
  }
  d$datasets  <- trimws(unlist(strsplit(d$datasets, ",")))
  d$min_score <- as.numeric(d$min_score)
  d
}

ARGS      <- parse_args()
OUTDIR    <- ARGS$outdir
DATASETS  <- ARGS$datasets
MIN_SCORE <- ARGS$min_score
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

msg <- function(...)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ==============================================================================
#  1.  PACKAGES
# ==============================================================================

need <- c("ggplot2", "dplyr", "tidyr", "patchwork", "ggalluvial",
          "scales", "openxlsx", "grid", "showtext")
for (p in need)
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org", quiet = TRUE)

suppressPackageStartupMessages({
  library(ggplot2);    library(dplyr);     library(tidyr)
  library(patchwork);  library(ggalluvial)
  library(scales);     library(openxlsx);  library(grid)
})

# ==============================================================================
#  2.  FONT SETUP  (Times New Roman -> fallback serif)
# ==============================================================================

FONT <- "serif"
tryCatch({
  if (requireNamespace("showtext", quietly = TRUE)) {
    showtext::showtext_auto()
    tnr_paths <- Sys.glob(c(
      "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf",
      "/usr/share/fonts/truetype/times-new-roman/Times_New_Roman.ttf",
      "/usr/share/fonts/TTF/times.ttf",
      "/usr/local/share/fonts/Times_New_Roman.ttf",
      "C:/Windows/Fonts/times.ttf"
    ))
    if (length(tnr_paths) > 0) {
      showtext::font_add("TNR", regular = tnr_paths[1])
      FONT <- "TNR"
      msg("Font: Times New Roman loaded from ", tnr_paths[1])
    } else {
      msg("Times New Roman not found; using default serif.")
    }
  }
}, error = function(e) NULL)

# ==============================================================================
#  3.  UNIFIED DESIGN TOKENS
# ==============================================================================

# -- Unicode symbols as variables (avoid \uXXXX inside backticks) --
SYM_CAP  <- "\u2229"   # intersection symbol
SYM_ARR  <- "\u2192"   # right arrow
SYM_DASH <- "\u2013"   # en dash
SYM_MUL  <- "\u00d7"   # multiplication sign
SYM_GEQ  <- "\u2265"   # greater or equal
SYM_APO  <- "\u2019"   # right single quote

# -- Canvas (identical for all four figures) --
FIG_W_IN  <- 7.2          # ~183 mm (Nature double column)
FIG_H_IN  <- 6.0          # ~152 mm
DPI       <- 600

# -- Font sizes (pt) -- all figures share these --
FS_TITLE   <- 11
FS_AXIS    <- 9
FS_TICK    <- 7.5
FS_LEGEND  <- 7
FS_ANNOT   <- 6.5

# -- Seven-region names --
RN <- list(
  r1 = "only_ThisStudy",
  r2 = "only_PRJNA766357",
  r3 = "only_PRJNA597465",
  r4 = paste0("ThisStudy ", SYM_CAP, " PRJNA766357"),
  r5 = paste0("ThisStudy ", SYM_CAP, " PRJNA597465"),
  r6 = paste0("PRJNA766357 ", SYM_CAP, " PRJNA597465"),
  r7 = "All_shared"
)

REGION_ORDER <- unlist(RN, use.names = FALSE)

# -- Colour palette (JAMA-inspired, colour-blind friendly) --
PAL <- setNames(
  c("#3B4992", "#EE0000", "#008B45", "#631879", "#E07233", "#008280", "#BB0021"),
  REGION_ORDER
)

# -- Region display labels --
REGION_LABELS <- setNames(
  c("ThisStudy\nonly",
    "PRJNA766357\nonly",
    "PRJNA597465\nonly",
    paste0("ThisStudy ", SYM_CAP, "\nPRJNA766357"),
    paste0("ThisStudy ", SYM_CAP, "\nPRJNA597465"),
    paste0("PRJNA766357 ", SYM_CAP, "\nPRJNA597465"),
    "All\nshared"),
  REGION_ORDER
)

# -- Taxonomy palette (12 phyla + Other) --
TAX_PAL <- c(
  "#3B4992", "#EE0000", "#008B45", "#631879", "#E07233",
  "#008280", "#BB0021", "#5F559B", "#A20056", "#D8A903",
  "#1B9E77", "#808180", "#CCCCCC"
)

# -- Master ggplot2 theme --
theme_pub <- function(base_size = FS_AXIS) {
  theme_classic(base_size = base_size, base_family = FONT) %+replace%
    theme(
      line                = element_line(linewidth = 0.35, colour = "black"),
      axis.line           = element_line(linewidth = 0.35),
      axis.ticks          = element_line(linewidth = 0.25),
      axis.ticks.length   = unit(1.2, "mm"),
      axis.text           = element_text(size = FS_TICK,  colour = "black"),
      axis.title          = element_text(size = FS_AXIS,  face = "bold"),
      legend.background   = element_blank(),
      legend.key          = element_blank(),
      legend.key.size     = unit(3, "mm"),
      legend.text         = element_text(size = FS_LEGEND),
      legend.title        = element_text(size = FS_LEGEND, face = "bold"),
      panel.background    = element_blank(),
      panel.grid          = element_blank(),
      plot.background     = element_blank(),
      plot.title          = element_text(size = FS_TITLE, face = "bold",
                                         hjust = 0.5, margin = margin(b = 3)),
      plot.subtitle       = element_text(size = FS_ANNOT, hjust = 0.5,
                                         colour = "grey35", margin = margin(b = 5)),
      plot.margin         = margin(8, 8, 6, 6),
      strip.background    = element_rect(fill = "grey94", colour = NA),
      strip.text          = element_text(size = FS_TICK, face = "bold")
    )
}

# ==============================================================================
#  4.  HELPER FUNCTIONS
# ==============================================================================

strip_ext <- function(x) {
  x <- basename(as.character(x))
  x <- sub("\\.gz$",  "", x, ignore.case = TRUE)
  sub("\\.(fa|fasta|fna|fas)$", "", x, ignore.case = TRUE)
}

extract_sample <- function(g) {
  s <- strip_ext(g)
  s <- sub("\\.bin\\..*$", "", s, ignore.case = TRUE)
  s <- sub("_bin_.*$",    "", s, ignore.case = TRUE)
  sub("\\..*$", "", s)
}

read_auto <- function(path) {
  ext <- tolower(tools::file_ext(path))
  sep <- if (ext %in% c("tsv", "txt", "xls")) "\t" else ","
  tryCatch(
    read.table(path, sep = sep, header = TRUE, quote = "\"",
               comment.char = "", check.names = FALSE, fill = TRUE,
               stringsAsFactors = FALSE, encoding = "UTF-8"),
    error = function(e)
      read.table(path, sep = "\t", header = TRUE, quote = "\"",
                 comment.char = "", check.names = FALSE, fill = TRUE,
                 stringsAsFactors = FALSE, encoding = "UTF-8")
  )
}

save_fig <- function(obj, prefix, w = FIG_W_IN, h = FIG_H_IN) {
  pp <- file.path(OUTDIR, paste0(prefix, c(".pdf", ".svg", ".jpg")))
  if (is.function(obj)) {
    cairo_pdf(pp[1], width = w, height = h); obj(); dev.off()
    svg(pp[2],       width = w, height = h); obj(); dev.off()
    jpeg(pp[3], width = w, height = h, units = "in", res = DPI, quality = 100)
    obj(); dev.off()
  } else {
    ggsave(pp[1], obj, width = w, height = h, device = cairo_pdf, dpi = DPI)
    ggsave(pp[2], obj, width = w, height = h, device = svg)
    ggsave(pp[3], obj, width = w, height = h, device = "jpeg",
           dpi = DPI, quality = 100)
  }
  msg("  -> ", prefix, "  (.pdf .svg .jpg)")
}

write_tsv <- function(df, name) {
  p <- file.path(OUTDIR, name)
  write.table(df, p, sep = "\t", quote = FALSE,
              row.names = FALSE, fileEncoding = "UTF-8")
  msg("  -> ", name)
}

parse_taxonomy <- function(cls) {
  ranks <- c("domain","phylum","class","order","family","genus","species")
  parts <- strsplit(as.character(cls), ";")
  mat <- do.call(rbind, lapply(parts, function(p) {
    p <- trimws(sub("^[a-z]__", "", trimws(p)))
    p[!nzchar(p)] <- NA_character_
    c(p, rep(NA_character_, max(0, 7 - length(p))))[1:7]
  }))
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  names(df) <- ranks
  df
}

# ==============================================================================
#  5.  READ DATA
# ==============================================================================

msg("Reading input files ...")
cdb  <- read_auto(ARGS$cdb)
wdb  <- read_auto(ARGS$wdb)
grp  <- read_auto(ARGS$group)
gtdb <- read_auto(ARGS$gtdb)
ko   <- read_auto(ARGS$ko)

# ==============================================================================
#  6.  MAG -> DATASET MAPPING
# ==============================================================================

cdb2 <- data.frame(
  genome  = as.character(cdb$genome),
  cluster = as.character(cdb$secondary_cluster),
  stringsAsFactors = FALSE
) %>%
  mutate(genome_id = strip_ext(genome),
         Sample    = extract_sample(genome))

grp2 <- data.frame(
  Sample  = trimws(as.character(grp$Sample)),
  Dataset = trimws(as.character(grp$Dataset)),
  stringsAsFactors = FALSE
) %>% distinct(Sample, .keep_all = TRUE)

merged <- left_join(cdb2, grp2, by = "Sample") %>%
  filter(!is.na(Dataset), nzchar(Dataset))

# ==============================================================================
#  7.  VENN REGIONS
# ==============================================================================

ds_sets <- lapply(setNames(DATASETS, DATASETS),
                  function(d) unique(merged$cluster[merged$Dataset == d]))
A <- ds_sets[[1]]; B <- ds_sets[[2]]; C <- ds_sets[[3]]

regions <- list()
regions[[RN$r1]] <- setdiff(A, union(B, C))
regions[[RN$r2]] <- setdiff(B, union(A, C))
regions[[RN$r3]] <- setdiff(C, union(A, B))
regions[[RN$r4]] <- setdiff(intersect(A, B), C)
regions[[RN$r5]] <- setdiff(intersect(A, C), B)
regions[[RN$r6]] <- setdiff(intersect(B, C), A)
regions[[RN$r7]] <- Reduce(intersect, list(A, B, C))

summary_df <- data.frame(
  Venn_region = names(regions),
  SGB_count   = vapply(regions, length, integer(1)),
  stringsAsFactors = FALSE
)

# ==============================================================================
#  8.  REPRESENTATIVE MAG  (Wdb highest score)
# ==============================================================================

cluster_first <- tapply(cdb2$genome, cdb2$cluster, `[`, 1)

wdb2 <- wdb %>%
  mutate(genome_id = strip_ext(genome),
         w_cluster = as.character(cluster),
         score     = as.numeric(score))

reps <- data.frame(cluster = unique(cdb2$cluster), stringsAsFactors = FALSE) %>%
  mutate(rep_mag    = unname(cluster_first[cluster]),
         rep_source = "fallback_first",
         rep_score  = NA_real_)

for (i in seq_len(nrow(reps))) {
  sub_w <- wdb2 %>% filter(w_cluster == reps$cluster[i])
  if (nrow(sub_w) > 0) {
    best <- sub_w %>% slice_max(score, n = 1, with_ties = FALSE)
    reps$rep_mag[i]    <- best$genome[1]
    reps$rep_source[i] <- "Wdb_best_score"
    reps$rep_score[i]  <- best$score[1]
  }
}

# ==============================================================================
#  9.  TAXONOMY
# ==============================================================================

gtdb2 <- gtdb %>%
  mutate(genome_id = strip_ext(user_genome)) %>%
  select(genome_id, classification,
         closest_genome_ani, closest_genome_af, msa_percent)

# ==============================================================================
#  10. BUILD LONG TABLE
# ==============================================================================

long_rows <- list()
for (rname in names(regions)) {
  for (cl in regions[[rname]]) {
    rr  <- reps %>% filter(cluster == cl)
    rid <- strip_ext(rr$rep_mag[1])
    tr  <- gtdb2 %>% filter(genome_id == rid)
    if (nrow(tr) == 0) {
      all_ids <- cdb2$genome_id[cdb2$cluster == cl]
      tr <- gtdb2 %>% filter(genome_id %in% all_ids) %>% slice(1)
    }
    long_rows[[length(long_rows) + 1]] <- data.frame(
      Venn_region        = rname,
      secondary_cluster  = cl,
      Representative_MAG = rr$rep_mag[1],
      Rep_source         = rr$rep_source[1],
      dRep_score         = rr$rep_score[1],
      Dataset_presence   = paste(
        DATASETS[vapply(DATASETS, function(d) cl %in% ds_sets[[d]], logical(1))],
        collapse = ";"
      ),
      classification     = if (nrow(tr) > 0) tr$classification[1] else "",
      closest_genome_ani = if (nrow(tr) > 0) as.character(tr$closest_genome_ani[1]) else NA_character_,
      closest_genome_af  = if (nrow(tr) > 0) as.character(tr$closest_genome_af[1])  else NA_character_,
      msa_percent        = if (nrow(tr) > 0) as.character(tr$msa_percent[1])        else NA_character_,
      stringsAsFactors   = FALSE
    )
  }
}

long_table <- bind_rows(long_rows)
tax_parsed <- parse_taxonomy(long_table$classification)
long_table <- bind_cols(long_table, tax_parsed)
long_table$Venn_region <- factor(long_table$Venn_region, levels = REGION_ORDER)

msg("Long table: ", nrow(long_table), " SGBs")

# -- Supplementary Tables 1 & 2 --
write_tsv(summary_df, "STable1_venn_region_summary.tsv")
write_tsv(long_table, "STable2_SGB_long_table.tsv")

# #############################################################################
#  FIGURE 1  --  VENN DIAGRAM  (grid graphics)
# #############################################################################

msg("Fig 1: Venn diagram ...")

draw_venn <- function() {
  counts <- vapply(regions, length, integer(1))
  total  <- sum(counts)
  pct    <- counts / max(total, 1) * 100
  lab    <- function(k) sprintf("%d\n(%.1f%%)", counts[[k]], pct[[k]])

  grid.newpage()
  pushViewport(viewport(width = 0.94, height = 0.94))
  grid.rect(gp = gpar(fill = "white", col = NA))

  # -- Title --
  grid.text("Species-level genome bin overlap across datasets",
            x = 0.50, y = 0.96,
            gp = gpar(fontsize = FS_TITLE, fontface = "bold", fontfamily = FONT))
  grid.text(sprintf("Total non-redundant SGBs = %d", total),
            x = 0.50, y = 0.915,
            gp = gpar(fontsize = FS_ANNOT, fontfamily = FONT, col = "grey35"))

  # -- Three circles --
  fill_c   <- c("#3B4992", "#EE0000", "#008B45")
  border_c <- c("#2A3570", "#AA0000", "#006030")
  cx <- c(0.385, 0.615, 0.500)
  cy <- c(0.560, 0.560, 0.385)
  r  <- 0.225
  for (i in 1:3)
    grid.circle(cx[i], cy[i], r,
                gp = gpar(fill = adjustcolor(fill_c[i], 0.28),
                          col = border_c[i], lwd = 1.6))

  # -- Dataset labels --
  lx <- c(0.21, 0.79, 0.50)
  ly <- c(0.81, 0.81, 0.115)
  for (i in 1:3)
    grid.text(DATASETS[i], lx[i], ly[i],
              gp = gpar(fontsize = FS_AXIS, fontface = "bold",
                        fontfamily = FONT, col = border_c[i]))

  # -- Count labels --
  pos <- list()
  pos[[RN$r1]] <- c(0.28, 0.59)
  pos[[RN$r2]] <- c(0.72, 0.59)
  pos[[RN$r3]] <- c(0.50, 0.23)
  pos[[RN$r4]] <- c(0.50, 0.635)
  pos[[RN$r5]] <- c(0.40, 0.42)
  pos[[RN$r6]] <- c(0.60, 0.42)
  pos[[RN$r7]] <- c(0.50, 0.505)

  for (k in names(pos)) {
    fg <- if (k == RN$r7) "#BB0021" else "black"
    fs <- if (k == RN$r7) FS_TICK + 1.5 else FS_TICK
    grid.text(lab(k), pos[[k]][1], pos[[k]][2],
              gp = gpar(fontsize = fs, fontface = "bold",
                        fontfamily = FONT, col = fg, lineheight = 0.82))
  }

  # -- Footnote --
  grid.text("Percentages relative to all non-redundant SGBs.",
            x = 0.50, y = 0.035,
            gp = gpar(fontsize = FS_ANNOT - 0.5, fontfamily = FONT, col = "grey50"))
  popViewport()
}

save_fig(draw_venn, "Fig1_Venn_SGB_overlap")

# #############################################################################
#  FIGURE 2  --  ALLUVIAL:  Venn region -> Phylum -> Genus
# #############################################################################

msg("Fig 2: Alluvial plot ...")

top_n_label <- function(x, n) {
  tab <- sort(table(x[!is.na(x) & nzchar(x)]), decreasing = TRUE)
  top <- names(tab)[seq_len(min(n, length(tab)))]
  ifelse(!is.na(x) & x %in% top, x, "Other")
}

allu_df <- long_table %>%
  filter(!is.na(phylum), nzchar(phylum)) %>%
  mutate(Phylum = top_n_label(phylum, 8),
         Genus  = top_n_label(genus,  10)) %>%
  count(Venn_region, Phylum, Genus, name = "Freq")

phy_lvl <- c(setdiff(sort(unique(allu_df$Phylum)), "Other"), "Other")
allu_df$Phylum <- factor(allu_df$Phylum, levels = phy_lvl)
phy_pal <- setNames(rep_len(TAX_PAL, length(phy_lvl)), phy_lvl)
phy_pal["Other"] <- "#CCCCCC"

sub2 <- paste0("Alluvial:  Venn region  ", SYM_ARR, "  Phylum  ", SYM_ARR, "  Genus")

fig2 <- ggplot(
  allu_df,
  aes(axis1 = Venn_region, axis2 = Phylum, axis3 = Genus, y = Freq)
) +
  geom_alluvium(aes(fill = Phylum), width = 1/5, alpha = 0.65,
                curve_type = "sigmoid", decreasing = FALSE) +
  geom_stratum(width = 1/5, fill = "grey97", colour = "grey50",
               linewidth = 0.25) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = FS_ANNOT / ggplot2::.pt, family = FONT, lineheight = 0.82) +
  scale_x_discrete(limits = c("Venn region", "Phylum", "Genus"),
                   expand = c(0.15, 0.05)) +
  scale_fill_manual(values = phy_pal) +
  labs(title    = "Taxonomic flow of SGBs across Venn regions",
       subtitle = sub2,
       y = "Number of SGBs") +
  theme_pub() +
  theme(legend.position  = "right",
        legend.key.size  = unit(3, "mm"),
        axis.title.x     = element_blank(),
        axis.text.x      = element_text(size = FS_AXIS, face = "bold"))

save_fig(fig2, "Fig2_Alluvial_Region_Phylum_Genus")

# -- S-Table 3 --
tax_counts <- long_table %>%
  group_by(Venn_region, phylum, class, order, family, genus) %>%
  summarise(n_SGB = n(), .groups = "drop") %>%
  arrange(Venn_region, desc(n_SGB))
write_tsv(tax_counts, "STable3_taxonomy_composition.tsv")

# #############################################################################
#  FIGURE 3  --  KEGG KO RICHNESS  (medium-/high-quality MAGs)
# #############################################################################

msg("Fig 3: KEGG KO richness (dRep score >= ", MIN_SCORE, ") ...")

ko_col   <- names(ko)[1]
ko_ids   <- as.character(ko[[ko_col]])
g_cols   <- setdiff(names(ko), ko_col)
ko_mat   <- as.matrix(
  data.frame(lapply(ko[g_cols], function(x) suppressWarnings(as.numeric(x))),
             check.names = FALSE)
)
rownames(ko_mat) <- ko_ids
colnames(ko_mat) <- g_cols
ko_pa <- ko_mat > 0
ko_pa[is.na(ko_pa)] <- FALSE

lt_ko <- long_table %>%
  mutate(rep_id = strip_ext(Representative_MAG)) %>%
  filter(!is.na(dRep_score), dRep_score >= MIN_SCORE,
         rep_id %in% colnames(ko_pa))

msg("  Passed filter: ", nrow(lt_ko), " / ", nrow(long_table), " SGBs")

richness          <- colSums(ko_pa[, lt_ko$rep_id, drop = FALSE])
lt_ko$KO_richness <- as.numeric(richness)

kegg_sum <- lt_ko %>%
  group_by(Venn_region) %>%
  summarise(n = n(), mean = mean(KO_richness), sd = sd(KO_richness),
            median = median(KO_richness), .groups = "drop") %>%
  mutate(sd = ifelse(is.na(sd), 0, sd))

kw <- kruskal.test(KO_richness ~ Venn_region, data = lt_ko)
kw_lab <- paste0("Kruskal", SYM_DASH, "Wallis  H = ",
                 sprintf("%.1f", kw$statistic), ",  P ",
                 ifelse(kw$p.value < 0.001, "< 0.001",
                        sprintf("= %.3f", kw$p.value)))

n_lab <- lt_ko %>% count(Venn_region, name = "nn")

qual_note <- paste0("Medium-/high-quality MAGs  (dRep score ", SYM_GEQ, " ",
                    MIN_SCORE, ")")
sub3 <- paste0("Diamond = mean;  box = IQR;  whiskers = 1.5", SYM_MUL, "IQR")

fig3 <- ggplot(lt_ko, aes(x = Venn_region, y = KO_richness)) +
  geom_boxplot(aes(fill = Venn_region), width = 0.55,
               outlier.shape = NA, linewidth = 0.3,
               show.legend = FALSE, alpha = 0.82) +
  geom_jitter(aes(colour = Venn_region), width = 0.16, size = 0.6,
              alpha = 0.50, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 2.2,
               colour = "black", show.legend = FALSE) +
  geom_text(data = n_lab, aes(x = Venn_region, y = -Inf,
                               label = paste0("n=", nn)),
            vjust = 1.6, size = FS_ANNOT / ggplot2::.pt,
            family = FONT, inherit.aes = FALSE) +
  scale_fill_manual(values = PAL) +
  scale_colour_manual(values = PAL) +
  scale_x_discrete(labels = REGION_LABELS) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  coord_cartesian(clip = "off") +
  annotate("text", x = Inf, y = Inf, label = kw_lab,
           hjust = 1.05, vjust = 1.6,
           size = FS_ANNOT / ggplot2::.pt, family = FONT, colour = "grey30") +
  annotate("text", x = Inf, y = Inf, label = qual_note,
           hjust = 1.05, vjust = 3.3,
           size = FS_ANNOT / ggplot2::.pt, family = FONT,
           fontface = "italic", colour = "grey45") +
  labs(x        = NULL,
       y        = "KO richness per representative MAG",
       title    = "KEGG functional repertoire across Venn regions",
       subtitle = sub3) +
  theme_pub() +
  theme(axis.text.x = element_text(size = FS_TICK, lineheight = 0.85))

save_fig(fig3, "Fig3_KEGG_KO_richness_boxplot")

# -- S-Tables 4 & 5 --
write_tsv(kegg_sum, "STable4_KEGG_richness_summary.tsv")
ko_detail <- lt_ko %>%
  select(Venn_region, secondary_cluster, Representative_MAG,
         dRep_score, classification, KO_richness)
write_tsv(ko_detail, "STable5_KEGG_richness_per_SGB.tsv")

# #############################################################################
#  FIGURE 4  --  ENRICHED KO HEATMAP  (Fisher + FDR)
# #############################################################################

msg("Fig 4: Enriched KO heatmap ...")

pa_rep <- ko_pa[, lt_ko$rep_id, drop = FALSE]
colnames(pa_rep) <- lt_ko$secondary_cluster

enrich_list <- list()
for (r in levels(droplevels(lt_ko$Venn_region))) {
  in_idx  <- which(lt_ko$Venn_region == r)
  out_idx <- which(lt_ko$Venn_region != r)
  n_in <- length(in_idx); n_out <- length(out_idx)
  if (n_in < 1 || n_out < 1) next

  in_prev  <- rowMeans(pa_rep[, in_idx,  drop = FALSE])
  out_prev <- rowMeans(pa_rep[, out_idx, drop = FALSE])
  in_pos   <- rowSums( pa_rep[, in_idx,  drop = FALSE])
  out_pos  <- rowSums( pa_rep[, out_idx, drop = FALSE])

  pval <- vapply(seq_along(ko_ids), function(j) {
    m <- matrix(c(in_pos[j], n_in  - in_pos[j],
                  out_pos[j], n_out - out_pos[j]), nrow = 2, byrow = TRUE)
    fisher.test(m, alternative = "greater")$p.value
  }, numeric(1))

  enrich_list[[r]] <- data.frame(
    Venn_region = r,
    KO          = ko_ids,
    in_prev     = in_prev,
    out_prev    = out_prev,
    prev_diff   = in_prev - out_prev,
    log2_OR     = log2(((in_pos + 0.5) / (n_in  - in_pos + 0.5)) /
                       ((out_pos + 0.5) / (n_out - out_pos + 0.5))),
    p_value     = pval,
    stringsAsFactors = FALSE
  )
}

enrich <- bind_rows(enrich_list)
enrich$FDR <- ave(enrich$p_value, enrich$Venn_region,
                  FUN = function(x) p.adjust(x, "BH"))

enrich_sorted <- enrich %>% arrange(Venn_region, FDR, desc(prev_diff))
write_tsv(enrich_sorted, "STable6_KEGG_enrichment_full.tsv")

# -- Select top 5 per region --
heat_sel <- enrich %>%
  filter(prev_diff > 0) %>%
  group_by(Venn_region) %>%
  arrange(FDR, desc(prev_diff)) %>%
  slice_head(n = 5) %>%
  ungroup()

if (nrow(heat_sel) > 40) heat_sel <- heat_sel[1:40, ]

heat_sel <- heat_sel %>%
  mutate(
    KO_id    = sub(":.*", "", KO),
    KO_desc  = trimws(sub(",.*", "", sub("^[^:]+:", "", KO))),
    KO_label = paste0(KO_id, " (", KO_desc, ")"),
    sig      = case_when(FDR < 0.001 ~ "***",
                         FDR < 0.01  ~ "**",
                         FDR < 0.05  ~ "*",
                         TRUE        ~ ""),
    Venn_region = factor(Venn_region, levels = REGION_ORDER)
  )

heat_sel$KO_label <- factor(heat_sel$KO_label,
                            levels = rev(unique(heat_sel$KO_label)))

heat_h <- max(FIG_H_IN, 0.18 * nrow(heat_sel) + 2.8)

sub4 <- paste0("Fisher", SYM_APO, "s exact test (one-sided);   ",
               "* FDR < 0.05    ** < 0.01    *** < 0.001")

fig4 <- ggplot(heat_sel, aes(x = Venn_region, y = KO_label, fill = prev_diff)) +
  geom_tile(colour = "white", linewidth = 0.45) +
  geom_text(aes(label = sig), size = FS_ANNOT / ggplot2::.pt,
            vjust = 0.72, fontface = "bold") +
  scale_fill_gradient2(low      = "#2166AC",
                       mid      = "white",
                       high     = "#B2182B",
                       midpoint = 0,
                       name     = "Prevalence\ndifference",
                       limits   = c(-1, 1),
                       breaks   = c(-1, -0.5, 0, 0.5, 1)) +
  scale_x_discrete(labels = REGION_LABELS, position = "top") +
  labs(title    = "Region-enriched KEGG orthologs",
       subtitle = sub4,
       x = NULL, y = NULL) +
  theme_pub() +
  theme(
    axis.text.x       = element_text(size = FS_TICK, lineheight = 0.85,
                                     angle = 45, hjust = 0),
    axis.text.y       = element_text(size = FS_ANNOT),
    legend.key.height = unit(10, "mm"),
    legend.key.width  = unit(3, "mm")
  )

save_fig(fig4, "Fig4_KEGG_enriched_KO_heatmap", w = FIG_W_IN, h = heat_h)

# -- S-Table 7 --
heat_export <- heat_sel %>%
  select(Venn_region, KO_id, KO_desc, in_prev, out_prev,
         prev_diff, log2_OR, p_value, FDR)
write_tsv(heat_export, "STable7_KEGG_enriched_KO_detail.tsv")

# #############################################################################
#  CONSOLIDATED EXCEL WORKBOOK
# #############################################################################

msg("Writing Excel workbook ...")

wb <- createWorkbook()
hs <- createStyle(fontName = "Arial", fontSize = 9, textDecoration = "bold",
                  fgFill = "#D9E2F3", border = "Bottom", borderColour = "#4472C4")

add_ws <- function(nm, df) {
  sn <- substr(nm, 1, 31)
  addWorksheet(wb, sn)
  writeData(wb, sn, df)
  addStyle(wb, sn, hs, rows = 1, cols = seq_len(ncol(df)), gridExpand = TRUE)
  setColWidths(wb, sn, cols = seq_len(ncol(df)), widths = "auto")
}

add_ws("S1_Region_Summary",        summary_df)
add_ws("S2_SGB_Long_Table",        long_table)
add_ws("S3_Taxonomy_Composition",  tax_counts)
add_ws("S4_KEGG_Richness_Summary", kegg_sum)
add_ws("S5_KEGG_Richness_Per_SGB", ko_detail)
add_ws("S6_KEGG_Enrichment_Full",  enrich_sorted)
add_ws("S7_Enriched_KO_Detail",    heat_export)

saveWorkbook(wb, file.path(OUTDIR, "Supplementary_Tables.xlsx"), overwrite = TRUE)
msg("  -> Supplementary_Tables.xlsx")

# ==============================================================================
#  METHODS / RESULTS  MARKDOWN
# ==============================================================================

md <- c(
  "# Dataset-level dRep MAG Species Composition Overlap",
  "",
  "## Methods",
  "",
  sprintf(paste0(
    "Species-level genome bins (SGBs) were defined by dRep secondary_cluster ",
    "(ANI >= 95%%). Representative MAGs were selected from Wdb.csv by highest ",
    "dRep score. KEGG functional analyses included only medium-/high-quality ",
    "MAGs (dRep score >= %g). Taxonomy was assigned via GTDB-Tk bac120. ",
    "Region-enriched KOs were identified by one-sided Fisher's exact test ",
    "with Benjamini-Hochberg correction."
  ), MIN_SCORE),
  "",
  "## Results",
  "",
  paste(sprintf("- %s: %d SGBs", summary_df$Venn_region, summary_df$SGB_count),
        collapse = "\n"),
  sprintf("\nKEGG analysis: %d SGBs passed quality filter.", nrow(lt_ko)),
  sprintf("Kruskal-Wallis: H = %.2f, P = %.4g", kw$statistic, kw$p.value),
  "",
  "## Figures",
  "",
  "| Figure | Description |",
  "|--------|-------------|",
  "| Fig 1 | Venn diagram: SGB overlap across three datasets |",
  "| Fig 2 | Alluvial plot: Venn region -> Phylum -> Genus |",
  sprintf("| Fig 3 | KO richness boxplot (dRep score >= %g) |", MIN_SCORE),
  "| Fig 4 | Enriched KO heatmap (Fisher + FDR) |",
  "",
  "## Supplementary Tables",
  "",
  "| Table | Description |",
  "|-------|-------------|",
  "| S1 | Venn region summary counts |",
  "| S2 | SGB-level long table with taxonomy |",
  "| S3 | Taxonomy composition per region |",
  "| S4 | KEGG KO richness summary per region |",
  "| S5 | Per-SGB KO richness detail |",
  "| S6 | Full KEGG enrichment (Fisher + FDR) |",
  "| S7 | Heatmap selected KO detail |"
)
writeLines(md, file.path(OUTDIR, "Methods_Results.md"))

msg("=== All done. ===")
