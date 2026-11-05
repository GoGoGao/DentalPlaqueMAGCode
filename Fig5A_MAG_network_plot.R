## Packaged for revise/; source: project/review/R1.4/Fig5_MAG_network_robust/scripts/01_plot_Fig5_robust_networks.R
## Figures: Fig5A
## Note: robust MAG networks
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig5A_MAG_network_plot.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# Rebuild Figure 5 MAG co-occurrence networks using resampling-robust edges only.
# Robust = bootstrap support ≥ 0.8 AND LOO support ≥ 0.8 (from Fig5_MAG_network_stability).

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(scales)
  library(ggrepel)
})

ROOT <- file.path(.pkg_root, ""
STAB <- file.path(ROOT, "review/R1.4/Fig5_MAG_network_stability/tables")
OUT  <- file.path(ROOT, "review/R1.4/Fig5_MAG_network_robust")
TAB  <- file.path(OUT, "tables")
FIG  <- file.path(OUT, "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

dual_write <- function(df, stem) {
  write_tsv(df, file.path(TAB, paste0(stem, ".tsv")))
  write_csv(df, file.path(TAB, paste0(stem, ".csv")))
}

save_plot <- function(p, stem, width, height, plot_df = NULL) {
  if (!is.null(plot_df)) dual_write(plot_df, paste0(stem, "_plotdata"))
  for (ext in c("pdf", "svg", "png", "jpg")) {
    ggsave(
      filename = file.path(FIG, paste0(stem, ".", ext)),
      plot = p, width = width, height = height, dpi = 300, bg = "white"
    )
  }
  message("  fig ", stem)
}

# ---- taxonomy / DA ----
gtdb <- read_tsv(
  file.path(ROOT, "01this_study_dataset/input/gtdbtk.bac120.summary.tsv"),
  show_col_types = FALSE
) %>%
  transmute(
    MAG = user_genome,
    classification = classification,
    phylum = str_replace(str_extract(classification, "p__[^;]+"), "^p__", ""),
    genus  = str_replace(str_extract(classification, "g__[^;]+"), "^g__", ""),
    species = str_replace(str_extract(classification, "s__.*$"), "^s__", "")
  ) %>%
  mutate(
    genus = if_else(is.na(genus) | genus == "", "Unassigned", genus),
    phylum = if_else(is.na(phylum) | phylum == "", "Unassigned", phylum),
    species_short = if_else(
      is.na(species) | species == "",
      MAG,
      str_replace_all(species, "_", " ")
    )
  )

da <- read_tsv(
  file.path(
    ROOT,
    "01this_study_dataset/sgb/metab_corr/paper_result/00_RobustDA_MAG_set/tables/RobustDA_MAGs.tsv"
  ),
  show_col_types = FALSE
) %>%
  select(MAG, enrichment = Robust_direction, species_plain) %>%
  distinct(MAG, .keep_all = TRUE)

# ---- robust edges ----
load_robust <- function(group) {
  read_tsv(file.path(STAB, paste0("Table_Fig5_edge_support_", group, ".tsv")),
           show_col_types = FALSE) %>%
    filter(in_observed_network == 1, robust_edge == 1) %>%
    transmute(
      group = group,
      source, target,
      correlation = r_obs,
      abs_r = abs_r_obs,
      p = p_obs,
      bootstrap_support, loo_support,
      sign = if_else(r_obs >= 0, "positive", "negative"),
      direct = "undirected",
      cor = if_else(r_obs >= 0, "+", "-")
    )
}

edges_N <- load_robust("Normal")
edges_P <- load_robust("PoorPlaque")
dual_write(edges_N, "Table_Fig5_robust_edges_Normal")
dual_write(edges_P, "Table_Fig5_robust_edges_PoorPlaque")
dual_write(bind_rows(edges_N, edges_P), "Table_Fig5_robust_edges_all")

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
    directed = FALSE,
    vertices = verts
  )
  set.seed(42)
  comm <- cluster_louvain(g, weights = abs(E(g)$correlation))
  V(g)$module <- as.character(membership(comm)[V(g)$name])
  V(g)$degree_v <- degree(g)
  V(g)$closeness_v <- closeness(g, normalized = TRUE)
  V(g)$betweenness_v <- betweenness(g, normalized = FALSE)
  g
}

net_stats <- function(g, group) {
  e_cor <- E(g)$correlation
  tibble(
    group = group,
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
    median_abs_r = median(abs(e_cor)),
    median_bootstrap_support = NA_real_
  )
}

plot_one_network <- function(g, group_label, seed = 123) {
  set.seed(seed)
  lay <- create_layout(g, layout = "fr")
  top_n <- min(10L, nrow(lay))
  hubs <- lay %>% arrange(desc(degree_v)) %>% slice_head(n = top_n) %>% pull(name)
  lay <- lay %>%
    mutate(
      hub_label = if_else(
        name %in% hubs,
        paste0(name, "\n", str_trunc(coalesce(label_hub, name), 28)),
        NA_character_
      )
    )

  mod_levels <- sort(unique(lay$module))
  mod_cols <- hue_pal()(length(mod_levels))
  names(mod_cols) <- mod_levels

  p <- ggraph(lay) +
    geom_edge_link(
      aes(color = sign, edge_alpha = abs_r, edge_width = abs_r),
      lineend = "round", show.legend = TRUE
    ) +
    scale_edge_color_manual(
      values = c(positive = "#377eb8", negative = "#e41a1c"),
      name = "Correlation"
    ) +
    scale_edge_alpha(range = c(0.15, 0.85), guide = "none") +
    scale_edge_width(range = c(0.2, 1.4), guide = "none") +
    geom_node_point(
      aes(size = degree_v, fill = module),
      color = "grey20", shape = 21, stroke = 0.25, alpha = 0.92
    ) +
    scale_size_continuous(range = c(1.8, 9), name = "Degree") +
    scale_fill_manual(values = mod_cols, name = "Module") +
    geom_text_repel(
      data = lay %>% filter(!is.na(hub_label)),
      aes(x = x, y = y, label = hub_label),
      inherit.aes = FALSE,
      size = 2.4, max.overlaps = 40,
      segment.alpha = 0.35, segment.size = 0.2,
      box.padding = 0.35, point.padding = 0.15,
      lineheight = 0.9, color = "black",
      family = "serif"
    ) +
    theme_void(base_family = "serif") +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey35"),
      plot.margin = margin(8, 8, 8, 8)
    ) +
    guides(
      edge_color = guide_legend(order = 1),
      size = guide_legend(order = 2),
      fill = guide_legend(order = 3, override.aes = list(size = 3.5))
    ) +
    labs(
      title = paste0("Robust MAG co-occurrence network: ", group_label),
      subtitle = paste0(
        vcount(g), " nodes | ", ecount(g), " robust edges | ",
        length(mod_levels), " modules (boot≥0.8 & LOO≥0.8)"
      )
    )

  list(plot = p, layout = lay, hubs = hubs)
}

export_node_edge_tables <- function(g, lay, group) {
  nodes <- as_tibble(lay) %>%
    transmute(
      group = !!group,
      MAG = name,
      x, y,
      degree = degree_v,
      closeness = closeness_v,
      betweenness = betweenness_v,
      module,
      enrichment,
      phylum, genus,
      species = coalesce(species_plain, species_short),
      is_top10_hub = name %in% (
        arrange(lay, desc(degree_v)) %>% slice_head(n = 10) %>% pull(name)
      )
    )
  edges <- as_data_frame(g, what = "edges") %>%
    transmute(
      group = !!group,
      source = from, target = to,
      correlation, abs_r, sign
    )
  dual_write(nodes, paste0("Table_Fig5_robust_nodes_", group))
  dual_write(edges, paste0("Fig_Fig5_robust_network_", group, "_edges_plotdata"))
  dual_write(nodes, paste0("Fig_Fig5_robust_network_", group, "_nodes_plotdata"))
  nodes
}

message("=== Build robust networks ===")
gN <- build_igraph(edges_N)
gP <- build_igraph(edges_P)

# attach median bootstrap from edge tables for summary
summ <- bind_rows(net_stats(gN, "Normal"), net_stats(gP, "PoorPlaque")) %>%
  mutate(
    median_bootstrap_support = c(
      median(edges_N$bootstrap_support),
      median(edges_P$bootstrap_support)
    )
  )
dual_write(summ, "Table_Fig5_robust_network_summary")
print(summ)

# hubs
hub_tbl <- bind_rows(
  tibble(group = "Normal", MAG = names(sort(degree(gN), decreasing = TRUE))) %>%
    mutate(degree_rank = row_number(), degree = as.integer(degree(gN)[MAG])) %>%
    filter(degree_rank <= 20) %>%
    left_join(gtdb, by = "MAG") %>%
    left_join(da, by = "MAG"),
  tibble(group = "PoorPlaque", MAG = names(sort(degree(gP), decreasing = TRUE))) %>%
    mutate(degree_rank = row_number(), degree = as.integer(degree(gP)[MAG])) %>%
    filter(degree_rank <= 20) %>%
    left_join(gtdb, by = "MAG") %>%
    left_join(da, by = "MAG")
)
dual_write(hub_tbl, "Table_Fig5_robust_hubs_top20")

shared <- intersect(
  hub_tbl %>% filter(group == "Normal", degree_rank <= 20) %>% pull(MAG),
  hub_tbl %>% filter(group == "PoorPlaque", degree_rank <= 20) %>% pull(MAG)
)
message("Shared top20 hubs: ", length(shared), " -> ", paste(shared, collapse = ", "))

message("=== Plot networks ===")
resN <- plot_one_network(gN, "Normal", seed = 123)
resP <- plot_one_network(gP, "PoorPlaque", seed = 456)
nodesN <- export_node_edge_tables(gN, resN$layout, "Normal")
nodesP <- export_node_edge_tables(gP, resP$layout, "PoorPlaque")

save_plot(resN$plot, "Fig_Fig5_robust_network_Normal", width = 7.2, height = 6.2,
          plot_df = NULL)
save_plot(resP$plot, "Fig_Fig5_robust_network_PoorPlaque", width = 7.2, height = 6.2,
          plot_df = NULL)

# side-by-side combined (drop legends on left to save space via patchwork if available)
if (requireNamespace("patchwork", quietly = TRUE)) {
  p_comb <- (resN$plot + theme(legend.position = "none")) +
    (resP$plot + theme(legend.position = "right")) +
    patchwork::plot_annotation(
      title = "Figure 5. Resampling-robust MAG co-occurrence networks",
      subtitle = "Edges with bootstrap support ≥ 0.8 and leave-one-out support ≥ 0.8",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 13, family = "serif"),
        plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey30", family = "serif")
      )
    )
  save_plot(p_comb, "Fig_Fig5_robust_network_combined", width = 13.5, height = 6.4)
} else {
  # fallback: save individually only
  message("patchwork not installed; skip combined panel")
}

# Fig 5B-style descriptive node metrics (no inferential Wilcoxon)
node_long <- bind_rows(nodesN, nodesP) %>%
  select(group, MAG, degree, closeness, betweenness) %>%
  pivot_longer(c(degree, closeness, betweenness), names_to = "metric", values_to = "value")
dual_write(node_long, "Fig_Fig5_robust_node_metrics_plotdata")
dual_write(node_long, "Table_Fig5_robust_node_metrics_long")

p_box <- ggplot(node_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(width = 0.55, outlier.size = 0.8, alpha = 0.85, colour = "grey25") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(Normal = "#4C72B0", PoorPlaque = "#C44E52")) +
  theme_bw(base_family = "serif", base_size = 10) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    plot.title = element_text(hjust = 0.5, size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 8, colour = "grey35")
  ) +
  labs(
    title = "Node topology in resampling-robust MAG networks",
    subtitle = "Descriptive distributions only (nodes are not independent replicates)",
    x = NULL, y = NULL
  )
save_plot(p_box, "Fig_Fig5_robust_node_metrics_box", width = 8.2, height = 3.4,
          plot_df = node_long)

# caption + conclusion
cap <- paste0(
  "Figure 5 (robust-edge revision). MAG co-occurrence networks rebuilt from resampling-stable edges only.\n\n",
  "Group-specific Spearman MAG–MAG edges were retained if bootstrap support ≥ 0.8 and leave-one-out support ≥ 0.8 ",
  "(participant-level resampling; B = 500; n = 10 per group). ",
  "Normal: ", summ$n_nodes[1], " nodes, ", summ$n_edges[1], " robust edges (",
  summ$n_positive[1], " positive / ", summ$n_negative[1], " negative; connectance = ",
  sprintf("%.3f", summ$connectance[1]), "; mean degree = ", sprintf("%.2f", summ$mean_degree[1]), "). ",
  "PoorPlaque: ", summ$n_nodes[2], " nodes, ", summ$n_edges[2], " robust edges (",
  summ$n_positive[2], " positive / ", summ$n_negative[2], " negative; connectance = ",
  sprintf("%.3f", summ$connectance[2]), "; mean degree = ", sprintf("%.2f", summ$mean_degree[2]), "). ",
  "Node size ∝ degree; colour = Louvain module; edge colour = positive (blue) / negative (red). ",
  "Labels mark top-10 degree hubs. Co-occurrence denotes statistical association, not direct microbial interaction; ",
  "analyses remain exploratory given small n.\n"
)
writeLines(cap, file.path(FIG, "Fig_Fig5_robust_network_CAPTION.txt"))
writeLines(cap, file.path(TAB, "Fig_Fig5_robust_network_CAPTION.txt"))

md <- paste0(
  "# Fig.5 robust-edge MAG co-occurrence networks\n\n",
  "## 边定义\n",
  "- 来自 `Fig5_MAG_network_stability`：观察边中 **bootstrap≥0.8 且 LOO≥0.8**\n",
  "- Normal **", summ$n_edges[1], "** / PoorPlaque **", summ$n_edges[2], "** 稳健边\n\n",
  "## 拓扑摘要\n",
  "| group | nodes | edges | +/− | connectance | mean degree |\n",
  "|---|---:|---:|---|---:|---:|\n",
  "| Normal | ", summ$n_nodes[1], " | ", summ$n_edges[1], " | ",
  summ$n_positive[1], "/", summ$n_negative[1], " | ",
  sprintf("%.4f", summ$connectance[1]), " | ", sprintf("%.2f", summ$mean_degree[1]), " |\n",
  "| PoorPlaque | ", summ$n_nodes[2], " | ", summ$n_edges[2], " | ",
  summ$n_positive[2], "/", summ$n_negative[2], " | ",
  sprintf("%.4f", summ$connectance[2]), " | ", sprintf("%.2f", summ$mean_degree[2]), " |\n\n",
  "## Hub\n",
  "- Shared top-20 hubs: ", length(shared),
  if (length(shared)) paste0(" (", paste(shared, collapse = ", "), ")") else "", "\n",
  "- 详见 `Table_Fig5_robust_hubs_top20`\n\n",
  "## 相对原稿\n",
  "- 原稿全网 771 / 1426 边；稳健网边数大幅收缩（约 28–32%）。\n",
  "- PoorPlaque 稳健边仍更多，但密度差远小于全网；勿再写 markedly denser / 节点 Wilcoxon 推断。\n",
  "- 共现 ≠ 互作；exploratory。\n\n",
  "## 输出\n`", OUT, "`\n"
)
writeLines(md, file.path(OUT, "Fig5_robust_network_结论.md"))
message("Done -> ", OUT)
