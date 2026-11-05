## Packaged for revise/; source: project/01this_study_dataset/sgb/metab_corr/paper_result/07_RobustDA_metab_KEGG_enrich/scripts/07_metab_KEGG_enrich_bubble.R
## Figures: Fig2C
## Note: metabolite KEGG enrichment
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig2C_metab_KEGG_enrich.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# 07 Robust_DA differential metabolites → KEGG pathway enrichment (by Compound ID)
# - KEGG Compound ID from select.metab.tsv column 16
# - One comparison: PoorPlaque-up vs PoorPlaque-down metabolites
#   (PoorPlaque-down ≡ Normal-up; PoorPlaque-up ≡ Normal-down)
# - Bubble style from Fig3_enrich_bubble.claude.R
# Output: paper_result/07_RobustDA_metab_KEGG_enrich/

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(glue)
  library(clusterProfiler)
  library(KEGGREST)
  library(ggnewscale)
})

set.seed(42)
options(stringsAsFactors = FALSE, scipen = 999, warn = 1)

base_dir   <- file.path(.pkg_root, "01this_study_dataset/sgb/metab_corr"
module_dir <- file.path(base_dir, "paper_result/07_RobustDA_metab_KEGG_enrich")
out_fig    <- file.path(module_dir, "figures")
out_tab    <- file.path(module_dir, "tables")
out_data   <- file.path(module_dir, "data")
cache_dir  <- file.path(module_dir, "cache")
for (d in c(out_fig, out_tab, out_data, cache_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
source(file.path(base_dir, "paper_result/_shared/export_utils.R"))

# X-axis focused on PoorPlaque direction
X_UP   <- "Up in PoorPlaque"
X_DOWN <- "Down in PoorPlaque"
COL_GROUP <- setNames(c("#B85450", "#2E8B8B"), c(X_UP, X_DOWN))
COL_SIDE  <- c(COL_GROUP, Shared = "#8A8A8A")
GRADIENT_COLS <- c("#FDEBD0", "#F0A868", "#D4572A", "#9C2010", "#5C0A04")
BASE_SIZE <- 22                      # doubled from 11 for publication readability
FONT_FAMILY <- "Times New Roman"
PATH_WRAP_CHARS <- 38L               # wrap pathway names longer than this to 2 lines
TOP_N <- 10L                 # max pathways kept per direction after relaxed cutoff
P_CUT_DISPLAY <- 0.05        # relaxed annotation threshold (nominal p)
MIN_GS <- 2L
# mega overview maps dilute metabolite enrichment signal
EXCLUDE_MAPS <- c("map01100", "map01110", "map01120", "map01130")

## ---------- helpers -----------------------------------------------------------
save_pub <- function(plot, stem, width_in = 8, height_in = 6, dpi = 300) {
  save_gg(plot, stem, w = width_in, h = height_in, dpi = dpi)
}

# Split long pathway names into exactly two lines at a whitespace near midpoint
wrap_two_lines <- function(x, width = PATH_WRAP_CHARS) {
  vapply(as.character(x), function(s) {
    if (is.na(s) || !nzchar(s)) return(s)
    if (nchar(s) <= width) return(s)
    mid <- ceiling(nchar(s) / 2)
    m <- gregexpr("\\s+", s, perl = TRUE)[[1]]
    if (m[1] < 0) {
      return(paste0(substr(s, 1, mid), "\n", substr(s, mid + 1L, nchar(s))))
    }
    idx <- which.min(abs(m - mid))
    br <- m[idx]
    len <- attr(m, "match.length")[idx]
    paste0(
      str_trim(substr(s, 1, br - 1L)),
      "\n",
      str_trim(substr(s, br + len, nchar(s)))
    )
  }, character(1), USE.NAMES = FALSE)
}

normalize_cpd <- function(x) {
  x <- as.character(x)
  x <- str_extract(x, "C\\d+")
  x
}

parse_generatio <- function(x) {
  sp <- strsplit(as.character(x), "/", fixed = TRUE)[[1]]
  if (length(sp) == 2L) as.numeric(sp[1]) / as.numeric(sp[2]) else NA_real_
}

cache_rds <- function(path, builder) {
  if (file.exists(path)) {
    message("[cache] load ", basename(path))
    return(readRDS(path))
  }
  obj <- builder()
  saveRDS(obj, path)
  message("[cache] save ", basename(path))
  obj
}

parse_br08901 <- function() {
  raw <- keggGet("br:br08901")[[1]]
  lines <- unlist(strsplit(paste(raw, collapse = "\n"), "\n", fixed = TRUE))
  cat_cur <- NA_character_
  sub_cur <- NA_character_
  out <- list()
  for (ln in lines) {
    if (grepl("^A", ln)) {
      cat_cur <- sub("^A", "", ln)
      sub_cur <- NA_character_
    } else if (grepl("^B", ln)) {
      sub_cur <- str_trim(sub("^B", "", ln))
    } else if (grepl("^C", ln)) {
      m <- str_match(ln, "^C\\s+(\\d{5})\\s+(.*)$")
      if (!is.na(m[1, 1])) {
        map_id <- paste0("map", m[1, 2])
        out[[length(out) + 1L]] <- tibble(
          ID = map_id,
          Description = str_trim(m[1, 3]),
          category = ifelse(is.na(cat_cur) || !nzchar(cat_cur), "Unclassified", cat_cur),
          subcategory = ifelse(is.na(sub_cur) || !nzchar(sub_cur), "Unclassified", sub_cur)
        )
      }
    }
  }
  bind_rows(out) %>% distinct(ID, .keep_all = TRUE)
}

get_cpd_pathway_map <- function() {
  link <- keggLink("pathway", "cpd")
  tibble(
    cpd = sub("^cpd:", "", names(link)),
    path = sub("^path:", "", unname(link))
  ) %>%
    filter(str_detect(path, "^map\\d+"), !path %in% EXCLUDE_MAPS)
}

annotate_enrich_df <- function(ek_df, path_anno) {
  if (is.null(ek_df) || nrow(ek_df) == 0L) return(ek_df)
  ek_df %>%
    as_tibble() %>%
    mutate(ID = as.character(ID)) %>%
    left_join(path_anno %>% select(ID, category, subcategory, Desc_anno = Description),
              by = "ID") %>%
    mutate(
      Description = ifelse(!is.na(Desc_anno) & nzchar(Desc_anno), Desc_anno, Description),
      category = ifelse(is.na(category) | !nzchar(category), "Unclassified", category),
      subcategory = ifelse(is.na(subcategory) | !nzchar(subcategory), "Unclassified", subcategory)
    ) %>%
    select(category, subcategory, ID, Description, GeneRatio, BgRatio,
           pvalue, p.adjust, qvalue, geneID, Count, everything(), -Desc_anno) %>%
    arrange(p.adjust, pvalue)
}

run_enrich <- function(gene_ids, universe, term2gene, term2name, path_anno,
                       p_cut = 0.05, q_cut = 0.2) {
  gene_ids <- unique(gene_ids[!is.na(gene_ids)])
  gene_ids <- intersect(gene_ids, unique(term2gene$gene))
  uni <- unique(universe[!is.na(universe)])
  uni <- intersect(uni, unique(term2gene$gene))
  if (length(gene_ids) < MIN_GS) {
    return(tibble())
  }
  ek <- tryCatch(
    enricher(
      gene          = gene_ids,
      universe      = uni,
      TERM2GENE     = term2gene,
      TERM2NAME     = term2name,
      pvalueCutoff  = 1,
      qvalueCutoff  = 1,
      minGSSize     = MIN_GS,
      maxGSSize     = 500,
      pAdjustMethod = "BH"
    ),
    error = function(e) {
      message("[WARN] enricher failed: ", e$message)
      NULL
    }
  )
  if (is.null(ek) || nrow(as.data.frame(ek)) == 0L) return(tibble())
  annotate_enrich_df(as.data.frame(ek), path_anno) %>%
    mutate(
      gene_set_size = length(gene_ids),
      universe_size = length(uni),
      sig_padj05 = p.adjust < p_cut,
      sig_q02 = qvalue < q_cut
    )
}

## ---------- select pathways: relaxed p + keep group-specific hits ------------
# Prefer pathways passing nominal p < P_CUT_DISPLAY in each group.
# Keep all Normal hits (usually few) and top_n PoorPlaque / other-group hits,
# prioritizing group-specific pathways so both sides appear on the plot.
select_group_enrich_paths <- function(normal_df, poor_df,
                                      p_cut = P_CUT_DISPLAY,
                                      top_n = TOP_N,
                                      x_levels = c("Normal", "PoorPlaque")) {
  req <- c("category", "subcategory", "ID", "Description",
           "pvalue", "p.adjust", "qvalue", "GeneRatio", "Count")
  for (df in list(normal_df, poor_df)) {
    miss <- setdiff(req, names(df))
    if (length(miss) > 0L) stop("缺少必要列：", paste(miss, collapse = ", "))
  }

  n_pass <- normal_df %>% filter(pvalue < p_cut) %>% arrange(pvalue, p.adjust)
  p_pass <- poor_df   %>% filter(pvalue < p_cut) %>% arrange(pvalue, p.adjust)

  # fallback if one side has zero hits under p_cut: take top_n by pvalue
  if (nrow(n_pass) == 0L && nrow(normal_df) > 0L) {
    n_pass <- normal_df %>% arrange(pvalue, p.adjust) %>% slice_head(n = min(3L, top_n))
    message(glue("[WARN] {x_levels[1]}: no pathway with p<{p_cut}; fallback top {nrow(n_pass)} by p"))
  }
  if (nrow(p_pass) == 0L && nrow(poor_df) > 0L) {
    p_pass <- poor_df %>% arrange(pvalue, p.adjust) %>% slice_head(n = top_n)
    message(glue("[WARN] {x_levels[2]}: no pathway with p<{p_cut}; fallback top {nrow(p_pass)} by p"))
  }

  n_ids <- n_pass$ID
  p_ids <- p_pass$ID
  shared_ids <- intersect(n_ids, p_ids)
  n_only <- setdiff(n_ids, p_ids)
  p_only <- setdiff(p_ids, n_ids)

  # Cap each side at top_n, but always keep ≥1 pathway per side when available.
  # Priority within each side: group-specific first, then shared.
  pick_side <- function(only_ids, shared, all_ordered, n_keep) {
    prefer <- unique(c(only_ids, shared))
    # preserve p-value order from all_ordered
    prefer <- all_ordered[all_ordered %in% prefer]
    if (length(prefer) == 0L) prefer <- all_ordered
    n_keep <- max(1L, as.integer(n_keep))
    prefer[seq_len(min(length(prefer), n_keep))]
  }
  keep_n <- pick_side(n_only, shared_ids, n_ids, top_n)
  keep_p <- pick_side(p_only, shared_ids, p_ids, top_n)

  side_of <- function(id) {
    in_n <- id %in% n_ids
    in_p <- id %in% p_ids
    if (in_n && in_p) "Shared" else if (in_n) x_levels[1] else x_levels[2]
  }

  # Bubble only where that group itself passes the display cutoff
  top_df <- bind_rows(
    n_pass %>% filter(ID %in% keep_n) %>% mutate(Group = x_levels[1]),
    p_pass %>% filter(ID %in% keep_p) %>% mutate(Group = x_levels[2])
  ) %>%
    mutate(
      Group = factor(Group, levels = x_levels),
      enrich_side = vapply(ID, side_of, character(1)),
      enrich_side = factor(enrich_side, levels = c(x_levels[1], x_levels[2], "Shared")),
      neg_log10_p = -log10(pmax(pvalue, 1e-300)),
      neg_log10_FDR = -log10(pmax(p.adjust, 1e-300)),
      GeneRatio_num = map_dbl(GeneRatio, parse_generatio),
      pass_p = pvalue < p_cut,
      pass_padj05 = p.adjust < 0.05
    )

  message(glue(
    "[select] p<{p_cut}: {x_levels[1]}={length(n_ids)} ({length(n_only)} specific), ",
    "{x_levels[2]}={length(p_ids)} ({length(p_only)} specific), shared={length(shared_ids)}; ",
    "display {x_levels[1]}={length(keep_n)}, {x_levels[2]}={length(keep_p)}, ",
    "pathways={length(unique(c(keep_n, keep_p)))} (bubbles={nrow(top_df)})"
  ))
  top_df
}

## ---------- bubble plot (redeveloped from Fig3_enrich_bubble.claude.R) --------
plot_enrich_bubble <- function(normal_df, poor_df,
                               p_cut = P_CUT_DISPLAY,
                               top_n = TOP_N,
                               count_legend = "Compound count",
                               x_levels = c("Normal", "PoorPlaque"),
                               title = NULL) {
  if ((is.null(normal_df) || nrow(normal_df) == 0L) &&
      (is.null(poor_df) || nrow(poor_df) == 0L)) {
    return(list(plot = NULL, data = tibble(), height = 5.5, width = 8.5))
  }

  top_df <- select_group_enrich_paths(
    normal_df, poor_df, p_cut = p_cut, top_n = top_n, x_levels = x_levels
  )
  if (nrow(top_df) == 0L) {
    return(list(plot = NULL, data = tibble(), height = 5.5, width = 8.5))
  }

  GAP <- 1.6
  ROW_STEP <- 1.65  # extra vertical room for 2-line Times labels at 2× font
  # order: category by mean significance; within category put
  # Normal-specific → Shared → PoorPlaque-specific, then by max -log10(p)
  desc_meta <- top_df %>%
    group_by(category, subcategory, Description, ID) %>%
    summarise(
      max_nlp = max(neg_log10_p),
      enrich_side = {
        sides <- unique(as.character(enrich_side))
        if (length(sides) > 1L) "Shared" else sides[1]
      },
      .groups = "drop"
    ) %>%
    mutate(
      side_rank = case_when(
        enrich_side == x_levels[1] ~ 1L,
        enrich_side == "Shared" ~ 2L,
        TRUE ~ 3L
      )
    )

  cat_order <- desc_meta %>%
    group_by(category) %>%
    summarise(mean_nlp = mean(max_nlp), .groups = "drop") %>%
    arrange(desc(mean_nlp)) %>%
    pull(category)

  # Keep subcategory blocks contiguous within each category (for strip + legend).
  desc_meta <- desc_meta %>%
    mutate(category = factor(category, levels = cat_order)) %>%
    arrange(category, subcategory, side_rank, max_nlp)

  y_pos <- numeric(nrow(desc_meta))
  cursor <- 0
  for (cat in levels(desc_meta$category)) {
    rows <- which(desc_meta$category == cat)
    y_pos[rows] <- cursor + seq_along(rows) * ROW_STEP
    cursor <- max(y_pos[rows]) + GAP
  }
  desc_meta$y_pos <- y_pos
  desc_meta$label_wrap <- wrap_two_lines(desc_meta$Description)
  desc_meta$enrich_side <- factor(
    desc_meta$enrich_side, levels = c(x_levels[1], x_levels[2], "Shared")
  )

  top_df <- top_df %>%
    select(-any_of(c("enrich_side", "y_pos"))) %>%
    left_join(desc_meta %>% select(Description, y_pos, enrich_side),
              by = "Description")

  n_cats <- length(cat_order)
  CAT_PALETTE <- setNames(
    c("#4E9E9A", "#D47B6A", "#7B9EC9", "#C4A35A",
      "#8BBF7A", "#B07BAF", "#E08C5A", "#6BAED6")[seq_len(n_cats)],
    cat_order
  )
  # Legend order follows top→bottom appearance on the plot
  subcat_levels <- unique(as.character(desc_meta$subcategory))
  SUBCAT_PALETTE <- setNames(
    colorRampPalette(c("#A8D5D1", "#F2C4BE", "#B8C9E0",
                       "#E8D5A0", "#C5E0B4", "#D9C0D8",
                       "#F5D0A9", "#C9B8E0"))(max(length(subcat_levels), 1L)),
    subcat_levels
  )

  shade_df <- desc_meta %>%
    group_by(category) %>%
    mutate(row_idx = row_number()) %>%
    ungroup() %>%
    filter(row_idx %% 2 == 0) %>%
    mutate(ymin = y_pos - ROW_STEP / 2, ymax = y_pos + ROW_STEP / 2)

  cat_bands <- desc_meta %>%
    group_by(category) %>%
    summarise(ymin = min(y_pos) - ROW_STEP / 2, ymax = max(y_pos) + ROW_STEP / 2,
              ymid = mean(y_pos), .groups = "drop") %>%
    mutate(
      cat_fill = CAT_PALETTE[as.character(category)],
      category = factor(category, levels = cat_order),
      cat_lab = wrap_two_lines(as.character(category), width = 18L)
    )

  # Contiguous runs of the same subcategory (one strip band per run)
  subcat_bands <- desc_meta %>%
    mutate(
      sub_run = cumsum(
        c(TRUE, as.character(subcategory)[-1] != as.character(subcategory)[-n()])
      )
    ) %>%
    group_by(sub_run, subcategory) %>%
    summarise(ymin = min(y_pos) - ROW_STEP / 2, ymax = max(y_pos) + ROW_STEP / 2,
              ymid = mean(y_pos), .groups = "drop") %>%
    mutate(
      subcategory = factor(as.character(subcategory), levels = subcat_levels),
      sub_fill = SUBCAT_PALETTE[as.character(subcategory)]
    )

  # enrich-side tick on far right of panel (group ownership)
  side_bands <- desc_meta %>%
    mutate(ymin = y_pos - 0.45, ymax = y_pos + 0.45,
           side_fill = COL_SIDE[as.character(enrich_side)])

  boundary_df <- cat_bands %>%
    arrange(desc(ymax)) %>%
    slice(-1) %>%
    mutate(boundary_y = ymax + GAP / 2)

  dumbbell_df <- top_df %>%
    group_by(Description) %>%
    filter(n() == 2) %>%
    ungroup()

  y_breaks <- desc_meta$y_pos
  y_labels <- desc_meta$label_wrap

  X_SUB_R <- 0.48
  X_SUB_L <- 0.38
  X_CAT_R <- 0.36
  X_CAT_L <- 0.14
  X_SIDE_L <- 2.48
  X_SIDE_R <- 2.58

  count_range <- range(top_df$Count, na.rm = TRUE)
  pretty_breaks_count <- scales::pretty_breaks(n = 4)(seq(count_range[1], count_range[2]))
  group_col <- COL_GROUP[x_levels]
  names(group_col) <- x_levels
  fill_lim <- range(top_df$neg_log10_p, na.rm = TRUE)
  # keep pale-but-visible fills for weaker (e.g. Normal) hits
  fill_lim[1] <- max(0.8, fill_lim[1] - 0.15)

  p <- ggplot(top_df, aes(x = Group, y = y_pos)) +
    ## ① alternating row shade
    geom_rect(
      data = shade_df,
      aes(ymin = ymin, ymax = ymax, xmin = 0.5, xmax = 2.5),
      inherit.aes = FALSE, fill = "#00000010", color = NA
    )

  if (nrow(boundary_df) > 0L) {
    p <- p + geom_hline(
      data = boundary_df, aes(yintercept = boundary_y),
      inherit.aes = FALSE, color = "#AAAAAA", linewidth = 0.5
    )
  }
  if (nrow(dumbbell_df) > 0L) {
    p <- p + geom_line(
      data = dumbbell_df, aes(group = Description),
      color = "#CCCCCC", linewidth = 0.45, alpha = 0.85
    )
  }

  p <- p +
    ## ④ halo
    geom_point(
      aes(size = Count * 1.8, color = Group),
      shape = 16, alpha = 0.12, show.legend = FALSE
    ) +
    ## ⑤ bubbles (stroke = group colour; fill = -log10 p)
    geom_point(
      data = filter(top_df, Group == x_levels[1]),
      aes(size = Count, fill = neg_log10_p),
      shape = 21, color = group_col[1], stroke = 1.15, alpha = 0.95
    ) +
    geom_point(
      data = filter(top_df, Group == x_levels[2]),
      aes(size = Count, fill = neg_log10_p),
      shape = 21, color = group_col[2], stroke = 1.15, alpha = 0.95
    ) +
    scale_color_manual(values = group_col, guide = "none") +
    scale_size_continuous(
      range = c(3.5, 13.5), name = count_legend, breaks = pretty_breaks_count
    ) +
    scale_fill_gradientn(
      colours = GRADIENT_COLS,
      limits = fill_lim,
      oob = scales::squish,
      name = expression(-log[10](italic(P))),
      guide = guide_colorbar(
        barwidth = unit(0.55, "cm"), barheight = unit(3.4, "cm"),
        ticks = TRUE, ticks.colour = "white",
        frame.colour = "#999999", title.position = "top",
        order = 1
      )
    ) +
    ## ⑥ subcategory strip + legend (new fill scale; no on-plot text)
    new_scale_fill() +
    geom_rect(
      data = subcat_bands,
      aes(ymin = ymin, ymax = ymax, xmin = X_SUB_L, xmax = X_SUB_R,
          fill = subcategory),
      inherit.aes = FALSE,
      color = "white", linewidth = 0.3
    ) +
    scale_fill_manual(
      name = "KEGG subcategory",
      values = SUBCAT_PALETTE,
      breaks = subcat_levels,
      drop = FALSE,
      guide = guide_legend(
        title.position = "top",
        order = 3,
        ncol = 1,
        override.aes = list(color = "white", linewidth = 0.3),
        keywidth = unit(0.55, "cm"),
        keyheight = unit(0.55, "cm")
      )
    ) +
    ## ⑦ category strip (primary class; text stays on strip)
    geom_rect(
      data = cat_bands,
      aes(ymin = ymin, ymax = ymax, xmin = X_CAT_L, xmax = X_CAT_R),
      fill = cat_bands$cat_fill, inherit.aes = FALSE,
      color = "white", linewidth = 0.3
    ) +
    geom_text(
      data = cat_bands,
      aes(x = (X_CAT_L + X_CAT_R) / 2, y = ymid, label = cat_lab),
      inherit.aes = FALSE, angle = 90, hjust = 0.5, vjust = 0.5,
      family = FONT_FAMILY, lineheight = 0.92,
      size = BASE_SIZE * 0.22, fontface = "bold", color = "white"
    ) +
    ## ⑧ enrich-side ownership ticks (who is enriched)
    geom_rect(
      data = side_bands,
      aes(ymin = ymin, ymax = ymax, xmin = X_SIDE_L, xmax = X_SIDE_R),
      fill = side_bands$side_fill, inherit.aes = FALSE, color = NA
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.42, 0.10))) +
    scale_y_continuous(
      breaks = y_breaks, labels = y_labels,
      expand = expansion(add = c(0.55, 0.55))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = NULL, y = NULL,
      title = title,
      subtitle = glue(
        "Display cutoff: nominal P < {p_cut}; left strips = category / subcategory ",
        "(legend); right ticks = enrichment ownership"
      )
    ) +
    theme_classic(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
    theme(
      text = element_text(family = FONT_FAMILY),
      panel.border = element_rect(color = "#888888", fill = NA, linewidth = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks.x = element_line(color = "#AAAAAA", linewidth = 0.35),
      axis.ticks.y = element_line(color = "#AAAAAA", linewidth = 0.35),
      axis.ticks.length = unit(3, "pt"),
      axis.text.x = element_text(
        family = FONT_FAMILY, size = BASE_SIZE, face = "bold", color = "black",
        angle = 30, hjust = 1, vjust = 1, margin = margin(t = 4)
      ),
      axis.text.y = element_text(
        family = FONT_FAMILY, size = BASE_SIZE * 0.80, color = "black",
        hjust = 1, lineheight = 0.95, margin = margin(r = 4)
      ),
      plot.title = element_text(
        family = FONT_FAMILY, face = "bold", size = BASE_SIZE * 1.05, color = "black"
      ),
      plot.subtitle = element_text(
        family = FONT_FAMILY, size = BASE_SIZE * 0.62, color = "#666666"
      ),
      legend.position = "right",
      legend.title = element_text(
        family = FONT_FAMILY, face = "bold", color = "black", size = BASE_SIZE * 0.78
      ),
      legend.text = element_text(
        family = FONT_FAMILY, color = "black", size = BASE_SIZE * 0.66
      ),
      legend.key = element_blank(),
      legend.spacing.y = unit(8, "pt"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.box = "vertical",
      legend.box.just = "left",
      plot.margin = margin(16, 20, 16, 28, unit = "pt"),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(
      size = guide_legend(
        title.position = "top",
        override.aes = list(fill = "#F0A868", color = "#555555", stroke = 0.8),
        order = 2
      )
    )

  total_y <- max(desc_meta$y_pos) + ROW_STEP
  fig_h <- max(10.5, total_y * 0.42 + 2.8)
  fig_w <- 14.5

  list(plot = p, data = top_df, meta = desc_meta, height = fig_h, width = fig_w)
}

## ---------- single-set bubble (one column; relaxed p cutoff) ------------------
plot_single_enrich_bubble <- function(df, group_label,
                                      p_cut = P_CUT_DISPLAY,
                                      top_n = TOP_N,
                                      count_legend = "Compound count") {
  if (is.null(df) || nrow(df) == 0L) {
    return(list(plot = NULL, data = tibble(), height = 5, width = 7))
  }
  plot_df <- df %>%
    filter(pvalue < p_cut) %>%
    arrange(pvalue, p.adjust)
  if (nrow(plot_df) == 0L) {
    plot_df <- df %>% arrange(pvalue, p.adjust) %>% slice_head(n = top_n)
  } else {
    plot_df <- plot_df %>% slice_head(n = top_n)
  }
  plot_df <- plot_df %>%
    mutate(
      Group = group_label,
      neg_log10_p = -log10(pmax(pvalue, 1e-300)),
      neg_log10_FDR = -log10(pmax(p.adjust, 1e-300)),
      GeneRatio_num = map_dbl(GeneRatio, parse_generatio),
      Description = factor(Description, levels = rev(unique(Description)))
    )

  col <- COL_GROUP[[group_label]]
  if (is.null(col) || is.na(col)) col <- "#555555"

  p <- ggplot(plot_df, aes(x = GeneRatio_num, y = Description)) +
    geom_point(
      aes(size = Count, fill = neg_log10_p),
      shape = 21, color = col, stroke = 1.05, alpha = 0.95
    ) +
    scale_size_continuous(range = c(2.5, 10), name = count_legend) +
    scale_fill_gradientn(
      colours = GRADIENT_COLS,
      name = expression(-log[10](italic(P)))
    ) +
    scale_x_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      x = "Compound ratio", y = NULL,
      title = glue("KEGG enrichment: {group_label}"),
      subtitle = glue("Pathways with nominal P < {p_cut} (top {top_n})")
    ) +
    theme_classic(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
    theme(
      text = element_text(family = FONT_FAMILY),
      panel.border = element_rect(color = "#888888", fill = NA, linewidth = 0.5),
      axis.line = element_blank(),
      plot.title = element_text(family = FONT_FAMILY, face = "bold", size = BASE_SIZE * 1.05),
      plot.subtitle = element_text(family = FONT_FAMILY, size = BASE_SIZE * 0.75, color = "#666666"),
      axis.text.y = element_text(family = FONT_FAMILY, size = BASE_SIZE * 0.85, color = "black",
                                lineheight = 0.95),
      legend.title = element_text(family = FONT_FAMILY, face = "bold", size = BASE_SIZE * 0.88)
    )

  fig_h <- max(4.5, nrow(plot_df) * 0.34 + 1.4)
  list(plot = p, data = plot_df, height = fig_h, width = 7.2)
}

## ============================================================================
## [1] Load metabolites + KEGG Compound IDs
## ============================================================================
message("[1] Load Robust_DA metabolites & select.metab.tsv col16 (KEGG Compound ID) ...")

sel_raw <- read_tsv(file.path(base_dir, "select.metab.tsv"), show_col_types = FALSE)
cpd_col <- names(sel_raw)[16]
message("  Column 16: ", cpd_col)

sel <- sel_raw %>%
  transmute(
    metab_id = `Metab ID`,
    Metabolite_anno = Metabolite,
    KEGG_Compound_ID = normalize_cpd(.data[[cpd_col]]),
    KEGG_Pathway_ID = `KEGG Pathway ID`,
    KEGG_Pathway_First = `KEGG Pathway First Category`,
    KEGG_Pathway_Second = `KEGG Pathway Second Category`
  )

da_file <- file.path(base_dir, "paper_result/05_RobustDA_metab_DA/tables/metabolite_RobustDA_only_with_HMDB_Class.tsv")
if (!file.exists(da_file)) {
  da_file <- file.path(base_dir, "output/tables/metabolite_RobustDA.tsv")
}
metab_da <- read_tsv(da_file, show_col_types = FALSE) %>%
  filter(if ("Robust_DA" %in% names(.)) Robust_DA else TRUE) %>%
  left_join(sel, by = "metab_id") %>%
  mutate(
    Metabolite = coalesce(Metabolite, Metabolite_anno, metab_id),
    # log2FC = log2(PoorPlaque / Normal); Normal-high => log2FC < 0
    log2FC = if ("log2FC" %in% names(.)) log2FC else
      log2((mean_PoorPlaque + 1e-6) / (mean_Normal + 1e-6)),
    Robust_direction = as.character(Robust_direction)
  )

# Two non-redundant sets (two-group DA):
#   PoorPlaque_up   = higher in PoorPlaque (= Normal down)
#   PoorPlaque_down = lower in PoorPlaque  (= Normal up)
sets <- list(
  PoorPlaque_up = metab_da %>%
    filter(Robust_direction == "PoorPlaque", !is.na(KEGG_Compound_ID)) %>%
    pull(KEGG_Compound_ID) %>% unique(),
  PoorPlaque_down = metab_da %>%
    filter(Robust_direction == "Normal", !is.na(KEGG_Compound_ID)) %>%
    pull(KEGG_Compound_ID) %>% unique()
)

input_map <- metab_da %>%
  filter(!is.na(KEGG_Compound_ID)) %>%
  mutate(
    PP_direction = ifelse(Robust_direction == "PoorPlaque", "Up in PoorPlaque", "Down in PoorPlaque")
  ) %>%
  select(metab_id, Metabolite, Robust_direction, PP_direction, log2FC,
         KEGG_Compound_ID, Wilcoxon_FDR, CLRMC_FDR) %>%
  distinct()
write_tsv(input_map, file.path(out_tab, "metab_RobustDA_with_KEGG_CompoundID.tsv"))

set_summary <- tibble(
  set = names(sets),
  n_compound = map_int(sets, length),
  definition = c(
    "Robust_direction=PoorPlaque: up in PoorPlaque (= down in Normal)",
    "Robust_direction=Normal: down in PoorPlaque (= up in Normal)"
  ),
  x_label = c(X_UP, X_DOWN)
)
write_tsv(set_summary, file.path(out_tab, "enrichment_gene_set_summary.tsv"))
message("  Gene-set sizes:\n", paste(sprintf("    %s: %d", set_summary$set, set_summary$n_compound), collapse = "\n"))

## ============================================================================
## [2] KEGG compound→pathway map + hierarchy
## ============================================================================
message("[2] Build KEGG compound–pathway mapping & pathway hierarchy ...")

map_df <- cache_rds(file.path(cache_dir, "kegg_cpd_pathway_map.rds"), get_cpd_pathway_map)
path_anno <- cache_rds(file.path(cache_dir, "kegg_br08901_pathway_anno.rds"), parse_br08901)

# pathway names from list as fallback
path_names <- cache_rds(file.path(cache_dir, "kegg_pathway_names.rds"), function() {
  kn <- keggList("pathway")
  ids <- sub("^path:", "", names(kn))
  ids <- ifelse(str_detect(ids, "^map"), ids, paste0("map", ids))
  tibble(ID = ids, Description = as.character(unname(kn)))
})

path_anno <- path_anno %>%
  left_join(path_names %>% select(ID, Desc_list = Description), by = "ID") %>%
  mutate(Description = ifelse(!is.na(Description) & nzchar(Description),
                              Description, coalesce(Desc_list, ID))) %>%
  select(-Desc_list)

term2gene <- map_df %>% transmute(term = path, gene = cpd)
term2name <- path_anno %>% transmute(term = ID, name = Description)
# add any map missing from br08901
miss_paths <- setdiff(unique(map_df$path), term2name$term)
if (length(miss_paths) > 0L) {
  extra <- path_names %>% filter(ID %in% miss_paths) %>% transmute(term = ID, name = Description)
  term2name <- bind_rows(term2name, extra) %>% distinct(term, .keep_all = TRUE)
  path_anno <- bind_rows(
    path_anno,
    tibble(ID = miss_paths, Description = term2name$name[match(miss_paths, term2name$term)],
           category = "Unclassified", subcategory = "Unclassified")
  ) %>% distinct(ID, .keep_all = TRUE)
}

universe <- unique(sel$KEGG_Compound_ID)
universe <- universe[!is.na(universe)]
universe <- intersect(universe, unique(term2gene$gene))
message("  Mapped pathways: ", n_distinct(map_df$path),
        "; universe compounds (detected ∩ mapped): ", length(universe))

write_tsv(map_df, file.path(out_data, "kegg_cpd_pathway_map_filtered.tsv"))
write_tsv(path_anno, file.path(out_data, "kegg_pathway_category_anno.tsv"))

## ============================================================================
## [3] Enrichment: PoorPlaque up vs down only
## ============================================================================
message("[3] Run KEGG enrichment for PoorPlaque-up / PoorPlaque-down ...")

enrich_list <- list()
for (nm in names(sets)) {
  message("  → ", nm, " (n=", length(sets[[nm]]), ")")
  ek_df <- run_enrich(sets[[nm]], universe, term2gene, term2name, path_anno)
  enrich_list[[nm]] <- ek_df
  out_tsv <- file.path(out_tab, glue("Fig_metab_KEGG_enrichment_{nm}.tsv"))
  if (nrow(ek_df) == 0L) {
    empty <- tibble(
      category = character(), subcategory = character(), ID = character(),
      Description = character(), GeneRatio = character(), BgRatio = character(),
      pvalue = double(), p.adjust = double(), qvalue = double(),
      geneID = character(), Count = integer()
    )
    write_tsv(empty, out_tsv)
    message("    no enrichment hits")
  } else {
    write_tsv(ek_df, out_tsv)
    sig_padj <- ek_df %>% filter(p.adjust < 0.05)
    sig_p <- ek_df %>% filter(pvalue < P_CUT_DISPLAY)
    write_tsv(sig_padj, file.path(out_tab, glue("Fig_metab_KEGG_enrichment_{nm}_padj05.tsv")))
    write_tsv(sig_p, file.path(out_tab, glue("Fig_metab_KEGG_enrichment_{nm}_p05.tsv")))
    message(glue(
      "    pathways: {nrow(ek_df)} (p<{P_CUT_DISPLAY}: {nrow(sig_p)}; padj<0.05: {nrow(sig_padj)})"
    ))
  }
}

enrich_all <- bind_rows(lapply(names(enrich_list), function(nm) {
  if (nrow(enrich_list[[nm]]) == 0L) return(tibble())
  mutate(enrich_list[[nm]], set = nm, .before = 1)
}))
write_tsv(enrich_all, file.path(out_tab, "Fig_metab_KEGG_enrichment_all_sets.tsv"))

## ============================================================================
## [4] ONE bubble: Up in PoorPlaque | Down in PoorPlaque
## ============================================================================
message("[4] Draw single comparison bubble (PoorPlaque up vs down) ...")

# clear previous multi-panel exports in figures/
old_figs <- list.files(out_fig, full.names = TRUE)
if (length(old_figs) > 0L) {
  unlink(old_figs)
  message("  cleaned ", length(old_figs), " old files in figures/")
}
old_data_figs <- list.files(out_data, pattern = "^Fig_metab_KEGG_enrich_bubble_", full.names = TRUE)
if (length(old_data_figs) > 0L) unlink(old_data_figs)

fig_res <- plot_enrich_bubble(
  enrich_list$PoorPlaque_up,
  enrich_list$PoorPlaque_down,
  p_cut = P_CUT_DISPLAY,
  top_n = TOP_N,
  count_legend = "Compound count",
  x_levels = c(X_UP, X_DOWN),
  title = "KEGG enrichment of PoorPlaque differential metabolites"
)
if (is.null(fig_res$plot)) stop("No pathways to plot after selection")

stem <- file.path(out_fig, "Fig_metab_KEGG_enrich_bubble")
save_pub(fig_res$plot, stem, width_in = fig_res$width, height_in = fig_res$height)
write_tsv(fig_res$data, paste0(stem, ".tsv"))
write_tsv(fig_res$data, file.path(out_data, "Fig_metab_KEGG_enrich_bubble.tsv"))
if (!is.null(fig_res$meta)) {
  meta_out <- fig_res$meta
  if ("label_wrap" %in% names(meta_out)) {
    meta_out$label_wrap <- gsub("\n", " | ", meta_out$label_wrap, fixed = TRUE)
  }
  write_tsv(meta_out, file.path(out_data, "Fig_metab_KEGG_enrich_bubble_pathway_meta.tsv"))
}

## ============================================================================
## Done
## ============================================================================
message("[DONE] One figure: ", stem, ".{pdf,svg,png,jpg,tsv}")
message("  tables/: PoorPlaque_up / PoorPlaque_down enrichment TSVs")
message("  data/: plot TSV + KEGG mapping")
