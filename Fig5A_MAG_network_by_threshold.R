## Packaged for revise/; source: project/review/R1.4/Fig5_MAG_network_rho_thresholds/scripts/02_plot_Fig5_by_absR_prev50.R
## Figures: Fig5A
## Note: networks by |r| threshold
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig5A_MAG_network_by_threshold.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# Figure 5 MAG co-occurrence networks after within-group prevalence filter (≥50%),
# then |r| threshold (0.8 / 0.9). Keeps ALL qualifying edges and incident nodes.
#
# Usage:
#   Rscript 02_plot_Fig5_by_absR_prev50.R 0.8
#   Rscript 02_plot_Fig5_by_absR_prev50.R 0.9

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(scales)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
RHO_MIN <- if (length(args) >= 1) as.numeric(args[[1]]) else 0.8
PREV_MIN <- if (length(args) >= 2) as.numeric(args[[2]]) else 0.5
stopifnot(is.finite(RHO_MIN), RHO_MIN > 0, RHO_MIN < 1)
stopifnot(is.finite(PREV_MIN), PREV_MIN > 0, PREV_MIN <= 1)

ROOT <- file.path(.pkg_root, ""
tag  <- sprintf("rho%03d_prev%02d", as.integer(round(RHO_MIN * 100)), as.integer(round(PREV_MIN * 100)))
OUT  <- file.path(ROOT, "review/R1.4", paste0("Fig5_MAG_network_", tag))
TAB  <- file.path(OUT, "tables")
FIG  <- file.path(OUT, "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

message("=== Fig.5 |r|>=", RHO_MIN, " + within-group prev>=", PREV_MIN, " -> ", OUT)

dual_write <- function(df, stem) {
  write_tsv(df, file.path(TAB, paste0(stem, ".tsv")))
  write_csv(df, file.path(TAB, paste0(stem, ".csv")))
}

save_plot <- function(p, stem, width, height) {
  for (ext in c("pdf", "svg", "png", "jpg")) {
    ggsave(file.path(FIG, paste0(stem, ".", ext)), p,
           width = width, height = height, dpi = 300, bg = "white")
  }
  message("  fig ", stem)
}

# ---- inputs ----
mag <- read_tsv(file.path(ROOT, "01this_study_dataset/input/bin_abundance_table.tsv"),
                show_col_types = FALSE) %>%
  column_to_rownames(var = names(.)[1]) %>%
  as.matrix()

grp <- read_tsv(
  file.path(ROOT, "01this_study_dataset/sgb/species_network/this.study.group.tsv"),
  col_names = c("Sample", "Group"), show_col_types = FALSE
)
stopifnot(all(grp$Sample %in% colnames(mag)))
mag <- mag[, grp$Sample, drop = FALSE]

gtdb <- read_tsv(
  file.path(ROOT, "01this_study_dataset/input/gtdbtk.bac120.summary.tsv"),
  show_col_types = FALSE
) %>%
  transmute(
    MAG = user_genome,
    classification,
    phylum = str_replace(str_extract(classification, "p__[^;]+"), "^p__", ""),
    genus  = str_replace(str_extract(classification, "g__[^;]+"), "^g__", ""),
    species = str_replace(str_extract(classification, "s__.*$"), "^s__", "")
  ) %>%
  mutate(
    genus = if_else(is.na(genus) | genus == "", "Unassigned", genus),
    phylum = if_else(is.na(phylum) | phylum == "", "Unassigned", phylum),
    species_short = if_else(is.na(species) | species == "", MAG, str_replace_all(species, "_", " "))
  )

da <- read_tsv(
  file.path(ROOT, "01this_study_dataset/sgb/metab_corr/paper_result/00_RobustDA_MAG_set/tables/RobustDA_MAGs.tsv"),
  show_col_types = FALSE
) %>%
  select(MAG, enrichment = Robust_direction, species_plain) %>%
  distinct(MAG, .keep_all = TRUE)

# ---- prevalence audit ----
prev_rows <- list()
for (gname in c("Normal", "PoorPlaque")) {
  smp <- grp$Sample[grp$Group == gname]
  X <- mag[, smp, drop = FALSE]
  n_s <- ncol(X)
  n_nz <- rowSums(X > 0)
  prev_rows[[gname]] <- tibble(
    MAG = rownames(X),
    group = gname,
    n_samples = n_s,
    n_nonzero = as.integer(n_nz),
    prevalence = n_nz / n_s,
    pass_prev50 = n_nz >= ceiling(PREV_MIN * n_s),
    n_zero = as.integer(n_s - n_nz),
    median_abundance = apply(X, 1, median),
    mean_abundance = rowMeans(X)
  )
}
prev_by_group <- bind_rows(prev_rows)
dual_write(prev_by_group, paste0("Table_Fig5_", tag, "_MAG_prevalence_by_group"))

prev_overall <- tibble(
  MAG = rownames(mag),
  n_samples = ncol(mag),
  n_nonzero = as.integer(rowSums(mag > 0)),
  prevalence = rowSums(mag > 0) / ncol(mag),
  pass_prev50_overall = rowSums(mag > 0) >= ceiling(PREV_MIN * ncol(mag))
) %>%
  left_join(gtdb %>% select(MAG, genus, species_short), by = "MAG")
dual_write(prev_overall, paste0("Table_Fig5_", tag, "_MAG_prevalence_overall"))

message("Prevalence filter (within-group >=", PREV_MIN, "):")
for (gname in c("Normal", "PoorPlaque")) {
  d <- prev_by_group %>% filter(group == gname)
  message(sprintf("  %s: pass %d / %d; fail %d; prev==0: %d",
                  gname, sum(d$pass_prev50), nrow(d), sum(!d$pass_prev50),
                  sum(d$n_nonzero == 0)))
}

# ---- rebuild Spearman edges on prevalence-filtered MAG set ----
spearman_edges <- function(X, names, r_min = 0.5, p_max = 0.05) {
  # X: features x samples
  n_f <- nrow(X); n_s <- ncol(X)
  if (n_s < 4 || n_f < 2) return(tibble())
  ranks <- t(apply(X, 1, rank))
  sd <- apply(ranks, 1, sd)
  valid <- which(sd > 0)
  if (length(valid) < 2) return(tibble())
  R <- cor(t(ranks[valid, , drop = FALSE]), method = "pearson")
  dof <- n_s - 2
  tstat <- R * sqrt(dof / pmax(1 - R^2, 1e-12))
  P <- 2 * pt(-abs(tstat), df = dof)
  idx <- valid
  out <- list()
  k <- 0L
  for (i in seq_along(idx)) {
    for (j in seq_along(idx)) {
      if (j <= i) next
      r <- R[i, j]; p <- P[i, j]
      if (is.finite(r) && is.finite(p) && abs(r) >= r_min && p < p_max) {
        k <- k + 1L
        out[[k]] <- tibble(
          source = names[idx[i]],
          target = names[idx[j]],
          correlation = r,
          p = p,
          abs_r = abs(r)
        )
      }
    }
  }
  if (k == 0) tibble() else bind_rows(out)
}

edge_list <- list()
kept_mags <- list()
for (gname in c("Normal", "PoorPlaque")) {
  smp <- grp$Sample[grp$Group == gname]
  keep <- prev_by_group %>% filter(group == gname, pass_prev50) %>% pull(MAG)
  kept_mags[[gname]] <- keep
  X <- mag[keep, smp, drop = FALSE]
  ed <- spearman_edges(X, keep, r_min = 0.5, p_max = 0.05) %>%
    filter(abs_r >= RHO_MIN) %>%
    mutate(
      group = gname,
      sign = if_else(correlation >= 0, "positive", "negative"),
      prevalence_filter = PREV_MIN,
      abs_r_threshold = RHO_MIN
    )
  edge_list[[gname]] <- ed
  message(sprintf("  %s retained MAGs=%d; |r|>=%s edges=%d",
                  gname, length(keep), RHO_MIN, nrow(ed)))
}

edges_N <- edge_list[["Normal"]]
edges_P <- edge_list[["PoorPlaque"]]
stopifnot(nrow(edges_N) > 0, nrow(edges_P) > 0)
dual_write(edges_N, paste0("Table_Fig5_", tag, "_edges_Normal"))
dual_write(edges_P, paste0("Table_Fig5_", tag, "_edges_PoorPlaque"))
dual_write(bind_rows(edges_N, edges_P), paste0("Table_Fig5_", tag, "_edges_all"))

# compare vs unfiltered threshold edges if available
for (gname in c("Normal", "PoorPlaque")) {
  unf_path <- file.path(
    ROOT, "review/R1.4",
    sprintf("Fig5_MAG_network_rho%03d/tables/Table_Fig5_rho%03d_edges_%s.tsv",
            as.integer(round(RHO_MIN * 100)), as.integer(round(RHO_MIN * 100)), gname)
  )
  if (file.exists(unf_path)) {
    unf <- read_tsv(unf_path, show_col_types = FALSE)
    message(sprintf("  vs unfiltered %s: edges %d -> %d (delta %d)",
                    gname, nrow(unf), nrow(edge_list[[gname]]),
                    nrow(edge_list[[gname]]) - nrow(unf)))
  }
}

build_igraph <- function(edge_df) {
  verts <- tibble(name = sort(unique(c(edge_df$source, edge_df$target)))) %>%
    left_join(gtdb, by = c("name" = "MAG")) %>%
    left_join(da, by = c("name" = "MAG")) %>%
    mutate(
      enrichment = replace_na(enrichment, "ns"),
      label_hub = coalesce(species_plain, species_short, name)
    )
  g <- graph_from_data_frame(
    d = edge_df %>% select(from = source, to = target, correlation, sign, abs_r),
    directed = FALSE, vertices = verts
  )
  stopifnot(ecount(g) == nrow(edge_df), vcount(g) == nrow(verts))
  set.seed(42)
  if (ecount(g) >= 1 && vcount(g) >= 3) {
    comm <- tryCatch(
      cluster_louvain(g, weights = abs(E(g)$correlation)),
      error = function(e) components(g)
    )
    if (inherits(comm, "communities")) {
      V(g)$module <- as.character(membership(comm)[V(g)$name])
    } else {
      V(g)$module <- as.character(comm$membership[V(g)$name])
    }
  } else {
    V(g)$module <- "1"
  }
  V(g)$degree_v <- degree(g)
  V(g)$closeness_v <- if (vcount(g) > 1) closeness(g, normalized = TRUE) else 0
  V(g)$betweenness_v <- betweenness(g, normalized = FALSE)
  g
}

net_stats <- function(g, group) {
  e_cor <- E(g)$correlation
  tibble(
    group = group,
    abs_r_threshold = RHO_MIN,
    prevalence_min = PREV_MIN,
    n_MAGs_pass_prev = length(kept_mags[[group]]),
    n_nodes = vcount(g),
    n_edges = ecount(g),
    n_positive = sum(e_cor >= 0),
    n_negative = sum(e_cor < 0),
    connectance = ecount(g) / (vcount(g) * (vcount(g) - 1) / 2),
    mean_degree = mean(degree(g)),
    median_degree = median(degree(g)),
    mean_closeness = mean(V(g)$closeness_v),
    mean_betweenness = mean(V(g)$betweenness_v),
    n_modules = length(unique(V(g)$module)),
    median_abs_r = median(abs(e_cor))
  )
}

plot_one_network <- function(g, group_label, size_limits, seed = 123) {
  set.seed(seed)
  lay <- create_layout(g, layout = "fr")
  hubs <- lay %>% arrange(desc(degree_v)) %>% slice_head(n = min(10L, nrow(lay))) %>% pull(name)
  lay <- lay %>% mutate(
    hub_label = if_else(
      name %in% hubs,
      paste0(name, "\n", str_trunc(coalesce(label_hub, name), 28)),
      NA_character_
    )
  )
  mod_levels <- sort(unique(as.character(lay$module)))
  mod_cols <- hue_pal()(max(length(mod_levels), 1))
  names(mod_cols) <- mod_levels

  p <- ggraph(lay) +
    geom_edge_link(aes(color = sign, edge_alpha = abs_r, edge_width = abs_r),
                   lineend = "round", show.legend = TRUE) +
    scale_edge_color_manual(values = c(positive = "#377eb8", negative = "#e41a1c"),
                            name = "Correlation") +
    scale_edge_alpha(range = c(0.15, 0.85), guide = "none") +
    scale_edge_width(range = c(0.2, 1.4), guide = "none") +
    geom_node_point(aes(size = degree_v, fill = module),
                    color = "grey20", shape = 21, stroke = 0.25, alpha = 0.92) +
    scale_size_continuous(range = c(1.8, 9), limits = size_limits, name = "Degree") +
    scale_fill_manual(values = mod_cols, name = "Module") +
    geom_text_repel(
      data = lay %>% filter(!is.na(hub_label)),
      aes(x = x, y = y, label = hub_label), inherit.aes = FALSE,
      size = 2.4, max.overlaps = 50, segment.alpha = 0.35, segment.size = 0.2,
      box.padding = 0.35, point.padding = 0.15, lineheight = 0.9,
      color = "black", family = "serif"
    ) +
    theme_void(base_family = "serif") +
    theme(
      legend.position = "right", legend.box = "vertical",
      legend.text = element_text(size = 8), legend.title = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey35"),
      plot.margin = margin(8, 8, 8, 8)
    ) +
    guides(edge_color = guide_legend(order = 1), size = guide_legend(order = 2),
           fill = guide_legend(order = 3, override.aes = list(size = 3.5))) +
    labs(
      title = paste0("MAG co-occurrence network: ", group_label),
      subtitle = paste0(
        "within-group prev>=", PREV_MIN, "; |r|>=", RHO_MIN,
        " (all qualifying edges) | ",
        vcount(g), " nodes | ", ecount(g), " edges | ",
        length(mod_levels), " modules"
      )
    )
  list(plot = p, layout = lay)
}

export_node_edge <- function(g, lay, group) {
  nodes <- as_tibble(lay) %>%
    transmute(
      group = !!group, abs_r_threshold = RHO_MIN, prevalence_min = PREV_MIN,
      MAG = name, x, y, degree = degree_v, closeness = closeness_v,
      betweenness = betweenness_v, module, enrichment, phylum, genus,
      species = coalesce(species_plain, species_short),
      is_top10_hub = name %in% (arrange(lay, desc(degree_v)) %>% slice_head(n = 10) %>% pull(name))
    )
  edges <- as_data_frame(g, what = "edges") %>%
    transmute(group = !!group, abs_r_threshold = RHO_MIN, prevalence_min = PREV_MIN,
              source = from, target = to, correlation, abs_r, sign)
  stopifnot(nrow(edges) == ecount(g))
  dual_write(nodes, paste0("Table_Fig5_", tag, "_nodes_", group))
  dual_write(edges, paste0("Fig_Fig5_", tag, "_network_", group, "_edges_plotdata"))
  dual_write(nodes, paste0("Fig_Fig5_", tag, "_network_", group, "_nodes_plotdata"))
  nodes
}

gN <- build_igraph(edges_N)
gP <- build_igraph(edges_P)
size_limits <- c(1, max(c(degree(gN), degree(gP), 1)))

summ <- bind_rows(net_stats(gN, "Normal"), net_stats(gP, "PoorPlaque"))
dual_write(summ, paste0("Table_Fig5_", tag, "_network_summary"))
print(summ)

hub_tbl <- bind_rows(
  tibble(group = "Normal", MAG = names(sort(degree(gN), decreasing = TRUE))) %>%
    mutate(degree_rank = row_number(), degree = as.integer(degree(gN)[MAG])) %>%
    filter(degree_rank <= 20) %>% left_join(gtdb, by = "MAG") %>% left_join(da, by = "MAG"),
  tibble(group = "PoorPlaque", MAG = names(sort(degree(gP), decreasing = TRUE))) %>%
    mutate(degree_rank = row_number(), degree = as.integer(degree(gP)[MAG])) %>%
    filter(degree_rank <= 20) %>% left_join(gtdb, by = "MAG") %>% left_join(da, by = "MAG")
) %>% mutate(abs_r_threshold = RHO_MIN, prevalence_min = PREV_MIN)
dual_write(hub_tbl, paste0("Table_Fig5_", tag, "_hubs_top20"))

resN <- plot_one_network(gN, "Normal", size_limits, 123)
resP <- plot_one_network(gP, "PoorPlaque", size_limits, 456)
nodesN <- export_node_edge(gN, resN$layout, "Normal")
nodesP <- export_node_edge(gP, resP$layout, "PoorPlaque")

stem <- paste0("Fig_Fig5_", tag, "_network")
save_plot(resN$plot, paste0(stem, "_Normal"), 7.2, 6.2)
save_plot(resP$plot, paste0(stem, "_PoorPlaque"), 7.2, 6.2)

if (requireNamespace("patchwork", quietly = TRUE)) {
  p_comb <- (resN$plot + theme(legend.position = "none")) +
    (resP$plot + theme(legend.position = "right")) +
    patchwork::plot_annotation(
      title = paste0("Figure 5. MAG co-occurrence (|r| >= ", RHO_MIN,
                     "; within-group prev >= ", PREV_MIN, ")"),
      subtitle = "All qualifying edges and incident nodes retained",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 13, family = "serif"),
        plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey30", family = "serif")
      )
    )
  save_plot(p_comb, paste0(stem, "_combined"), 13.5, 6.4)
}

node_long <- bind_rows(nodesN, nodesP) %>%
  select(group, MAG, degree, closeness, betweenness) %>%
  pivot_longer(c(degree, closeness, betweenness), names_to = "metric", values_to = "value")
dual_write(node_long, paste0("Fig_Fig5_", tag, "_node_metrics_plotdata"))
p_box <- ggplot(node_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.55, outlier.size = 0.8, alpha = 0.85, colour = "grey25") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(Normal = "#4C72B0", PoorPlaque = "#C44E52")) +
  theme_bw(base_family = "serif", base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        plot.title = element_text(hjust = 0.5, size = 11)) +
  labs(title = paste0("Node topology (|r|>=", RHO_MIN, "; prev>=", PREV_MIN, ")"),
       subtitle = "Descriptive only", x = NULL, y = NULL)
save_plot(p_box, paste0("Fig_Fig5_", tag, "_node_metrics_box"), 8.2, 3.4)

# failed-prev MAG list for Normal (main effect)
failN <- prev_by_group %>% filter(group == "Normal", !pass_prev50) %>%
  arrange(prevalence, MAG) %>%
  left_join(gtdb %>% select(MAG, genus, species_short), by = "MAG")
dual_write(failN, paste0("Table_Fig5_", tag, "_MAGs_fail_prev_Normal"))

cap <- paste0(
  "Figure 5. MAG co-occurrence networks after within-group prevalence filtering ",
  "(nonzero in >=", PREV_MIN * 100, "% of group samples) and |r| >= ", RHO_MIN, ".\n\n",
  "MAGs were retained within each group if detected (abundance > 0) in at least ",
  "ceil(", PREV_MIN, "*n) samples (n = 10/group). Spearman networks were rebuilt on the filtered ",
  "feature set (|r| >= 0.5, P < 0.05), then edges with |r| >= ", RHO_MIN, " were kept in full ",
  "(no display pruning). ",
  "Normal: ", summ$n_MAGs_pass_prev[1], " MAGs pass prev; network ",
  summ$n_nodes[1], " nodes / ", summ$n_edges[1], " edges (connectance ",
  sprintf("%.3f", summ$connectance[1]), "). ",
  "PoorPlaque: ", summ$n_MAGs_pass_prev[2], " MAGs pass prev; network ",
  summ$n_nodes[2], " nodes / ", summ$n_edges[2], " edges (connectance ",
  sprintf("%.3f", summ$connectance[2]), "). ",
  "Exploratory; co-occurrence ≠ interaction.\n"
)
writeLines(cap, file.path(FIG, paste0(stem, "_CAPTION.txt")))
writeLines(cap, file.path(TAB, paste0(stem, "_CAPTION.txt")))

md <- paste0(
  "# Fig.5 |r|>=", RHO_MIN, " + within-group prev>=", PREV_MIN, "\n\n",
  "## Prevalence\n",
  "- Normal pass: ", summ$n_MAGs_pass_prev[1], "/125; fail listed in ",
  "`Table_Fig5_", tag, "_MAGs_fail_prev_Normal`\n",
  "- PoorPlaque pass: ", summ$n_MAGs_pass_prev[2], "/125\n\n",
  "## Network\n",
  "| group | pass prev | nodes | edges | connectance | mean deg |\n",
  "|---|---:|---:|---:|---:|---:|\n",
  "| Normal | ", summ$n_MAGs_pass_prev[1], " | ", summ$n_nodes[1], " | ",
  summ$n_edges[1], " | ", sprintf("%.4f", summ$connectance[1]), " | ",
  sprintf("%.2f", summ$mean_degree[1]), " |\n",
  "| PoorPlaque | ", summ$n_MAGs_pass_prev[2], " | ", summ$n_nodes[2], " | ",
  summ$n_edges[2], " | ", sprintf("%.4f", summ$connectance[2]), " | ",
  sprintf("%.2f", summ$mean_degree[2]), " |\n\n",
  "## Output\n`", OUT, "`\n"
)
writeLines(md, file.path(OUT, paste0("Fig5_", tag, "_结论.md")))
message("Done -> ", OUT)
