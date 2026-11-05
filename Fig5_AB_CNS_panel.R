## Packaged for revise/; source: project/review/R1.4/Fig5_MAG_network_rho080_prev50/scripts/04_plot_Fig5_AB_CNS_panel.R
## Figures: Fig5A;Fig5B
## Note: Fig5 AB collage
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig5_AB_CNS_panel.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# CNS-style combined Figure 5 (A networks + B topology) for rho080_prev50.
# Liberation Serif (= Times metrics), larger labels, legends over in-plot text.

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(ggrepel)
  library(ggpubr)
  library(patchwork)
  library(cowplot)
  library(showtext)
  library(sysfonts)
})

OUT <- file.path(.pkg_root, "review/R1.4/Fig5_MAG_network_rho080_prev50"
TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")
STAB <- file.path(.pkg_root, "review/R1.4/Fig5_MAG_network_stability/tables"

# Times-compatible serif
font_add("TimesNewRoman",
         regular = "/usr/share/fonts/liberation/LiberationSerif-Regular.ttf",
         bold    = "/usr/share/fonts/liberation/LiberationSerif-Bold.ttf",
         italic  = "/usr/share/fonts/liberation/LiberationSerif-Italic.ttf",
         bolditalic = "/usr/share/fonts/liberation/LiberationSerif-BoldItalic.ttf")
showtext_auto()
showtext_opts(dpi = 320)
FONT <- "TimesNewRoman"

# Coordinated CNS palette
COL_N <- "#1F6F6B"
COL_P <- "#B83A3A"
COL_POS <- "#2F5D8C"
COL_NEG <- "#C46B4A"
COL_UNST <- "#B0B0B0"  # unstable: bootstrap/LOO support < 0.8
COL_NODE_N <- "#5BA8A2"
COL_NODE_P <- "#D48A8A"
COL_INK <- "#222222"

dual_write <- function(df, stem) {
  write_tsv(df, file.path(TAB, paste0(stem, ".tsv")))
  write_csv(df, file.path(TAB, paste0(stem, ".csv")))
}

# Save vector/raster without jpeg device (avoids colour corruption); JPG via png later
save_png_pdf_svg <- function(p, stem, w, h) {
  ggsave(file.path(FIG, paste0(stem, ".pdf")), p, width = w, height = h,
         device = cairo_pdf, bg = "white")
  ggsave(file.path(FIG, paste0(stem, ".svg")), p, width = w, height = h, bg = "white")
  ggsave(file.path(FIG, paste0(stem, ".png")), p, width = w, height = h,
         dpi = 320, bg = "white", limitsize = FALSE)
  message("saved ", stem, " (pdf/svg/png)")
}

edge_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "|")

load_group <- function(group) {
  edges <- read_csv(file.path(TAB, paste0("Table_Fig5_rho080_prev50_edges_", group, ".csv")),
                    show_col_types = FALSE)
  nodes0 <- read_csv(file.path(TAB, paste0("Table_Fig5_rho080_prev50_nodes_", group, ".csv")),
                     show_col_types = FALSE)
  support <- read_tsv(
    file.path(STAB, paste0("Table_Fig5_edge_support_prev50_", group, ".tsv")),
    show_col_types = FALSE
  ) %>%
    mutate(ekey = edge_key(source, target)) %>%
    select(ekey, bootstrap_support, loo_support, robust_edge)

  edges <- edges %>%
    mutate(ekey = edge_key(source, target)) %>%
    left_join(support, by = "ekey")
  stopifnot(!any(is.na(edges$robust_edge)))

  # Draw unstable first so Positive/Negative sit on top
  edges <- edges %>%
    mutate(
      edge_type = case_when(
        robust_edge == 0 ~ "Unstable",
        correlation >= 0 ~ "Positive",
        TRUE ~ "Negative"
      ),
      edge_type = factor(edge_type, levels = c("Unstable", "Positive", "Negative"))
    ) %>%
    arrange(edge_type)

  dual_write(
    edges %>% transmute(
      group, source, target, correlation, abs_r, p,
      bootstrap_support, loo_support, robust_edge, edge_type = as.character(edge_type)
    ),
    paste0("Fig_Fig5_rho080_prev50_CNS_network_", group, "_edges_plotdata")
  )

  g <- graph_from_data_frame(
    d = edges %>% transmute(
      from = source, to = target, correlation, abs_r,
      bootstrap_support, loo_support, robust_edge,
      edge_type = as.character(edge_type)
    ),
    directed = FALSE,
    vertices = nodes0 %>% transmute(name = MAG, species = coalesce(species, MAG))
  )
  V(g)$degree_v <- degree(g)
  V(g)$closeness_v <- closeness(g, normalized = TRUE)
  V(g)$betweenness_v <- betweenness(g, normalized = FALSE)
  V(g)$cen_degree <- degree(g) / pmax(vcount(g) - 1, 1)
  g
}

gN <- load_group("Normal")
gP <- load_group("PoorPlaque")
size_lim <- c(1, max(c(degree(gN), degree(gP))))

short_sp <- function(x) {
  x <- str_replace_all(x, "_", " ")
  # Drop trailing MAG IDs if species field fell back to MAG
  x <- if_else(str_detect(x, "^MAG\\d+$"), x, x)
  parts <- str_split_fixed(x, "\\s+", 3)
  # Prefer Genus + species/epithet; keep Genus-only if nothing else
  out <- ifelse(parts[, 2] == "", parts[, 1], paste(parts[, 1], parts[, 2]))
  # Soft-wrap long binomials for hub labels (Genus\nepithet)
  ifelse(nchar(out) > 18 & parts[, 2] != "",
         paste(parts[, 1], parts[, 2], sep = "\n"),
         out)
}

plot_net <- function(g, group, node_fill, seed, show_guides = TRUE) {
  set.seed(seed)
  lay <- create_layout(g, layout = "fr") %>%
    mutate(
      species_short = short_sp(species),
      hub_rank = rank(-degree_v, ties.method = "first"),
      # Core hubs labelled by species name (not MAG ID)
      hub_lab = if_else(hub_rank <= 5, species_short, NA_character_)
    )

  dual_write(
    lay %>% transmute(
      group = !!group, MAG = name, x, y, degree = degree_v,
      closeness = closeness_v, betweenness = betweenness_v, cen_degree,
      species = species_short, hub_rank, hub_lab
    ),
    paste0("Fig_Fig5_rho080_prev50_CNS_network_", group, "_layout_plotdata")
  )

  p <- ggraph(lay) +
    # Robust Positive / Negative first
    geom_edge_link(
      aes(filter = edge_type != "Unstable", colour = edge_type,
          edge_width = abs_r, edge_alpha = abs_r),
      lineend = "round", show.legend = show_guides
    ) +
    # Unstable on top (grey; fixed alpha so grey stays visible among hubs)
    geom_edge_link(
      aes(filter = edge_type == "Unstable", colour = edge_type,
          edge_width = abs_r),
      edge_alpha = 0.70, lineend = "round", show.legend = show_guides
    ) +
    scale_edge_colour_manual(
      values = c(Positive = COL_POS, Negative = COL_NEG, Unstable = COL_UNST),
      breaks = c("Positive", "Negative", "Unstable"),
      name = "Edge"
    ) +
    scale_edge_width(range = c(0.3, 1.25), guide = "none") +
    scale_edge_alpha(range = c(0.22, 0.65), guide = "none") +
    geom_node_point(
      aes(size = degree_v),
      fill = node_fill, colour = "#333333", shape = 21,
      stroke = 0.4, alpha = 0.95
    ) +
    scale_size_continuous(
      name = "Degree", range = c(2.4, 10.5), limits = size_lim,
      breaks = c(5, 10, 15)
    ) +
    geom_text_repel(
      data = filter(lay, !is.na(hub_lab)),
      aes(x = x, y = y, label = hub_lab),
      inherit.aes = FALSE,
      family = FONT,
      size = 4.2,
      colour = COL_INK,
      fontface = "italic",
      lineheight = 0.92,
      max.overlaps = 40,
      box.padding = 0.55,
      point.padding = 0.35,
      segment.color = "#777777",
      segment.size = 0.35,
      segment.alpha = 0.75,
      min.segment.length = 0,
      force = 1.6,
      seed = 1
    ) +
    labs(title = group) +
    theme_void(base_family = FONT) +
    theme(
      text = element_text(family = FONT, colour = COL_INK),
      plot.title = element_text(family = FONT, size = 18, hjust = 0.5,
                                margin = margin(b = 8), colour = COL_INK),
      legend.title = element_text(family = FONT, size = 13),
      legend.text = element_text(family = FONT, size = 12),
      legend.position = if (show_guides) "right" else "none",
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(8, 6, 8, 6)
    ) +
    guides(
      edge_colour = guide_legend(order = 1, override.aes = list(edge_width = 1.4, edge_alpha = 1)),
      size = guide_legend(order = 2)
    )
  list(plot = p, layout = lay)
}

resN <- plot_net(gN, "Normal", COL_NODE_N, 123, show_guides = FALSE)
resP <- plot_net(gP, "PoorPlaque", COL_NODE_P, 456, show_guides = TRUE)

hub_key <- bind_rows(
  resN$layout %>% filter(!is.na(hub_lab)) %>%
    transmute(panel = "Normal", MAG = name, degree = degree_v,
              species = str_replace_all(hub_lab, "\n", " ")),
  resP$layout %>% filter(!is.na(hub_lab)) %>%
    transmute(panel = "PoorPlaque", MAG = name, degree = degree_v,
              species = str_replace_all(hub_lab, "\n", " "))
) %>% arrange(panel, desc(degree))
dual_write(hub_key, "Table_Fig5_rho080_prev50_CNS_hub_label_key")

# Shared legend extracted from PoorPlaque plot
leg <- cowplot::get_legend(resP$plot + theme(legend.position = "right"))
p_nets <- cowplot::plot_grid(
  resN$plot + theme(legend.position = "none"),
  resP$plot + theme(legend.position = "none"),
  ncol = 2, align = "hv", rel_widths = c(1, 1)
)
p_a <- cowplot::plot_grid(p_nets, leg, ncol = 2, rel_widths = c(1, 0.16))

# Panel B
node_long <- bind_rows(
  as_tibble(resN$layout) %>% transmute(
    MAG = name, group = "Normal",
    betweenness = betweenness_v, cen_degree, closeness = closeness_v, degree = degree_v
  ),
  as_tibble(resP$layout) %>% transmute(
    MAG = name, group = "PoorPlaque",
    betweenness = betweenness_v, cen_degree, closeness = closeness_v, degree = degree_v
  )
) %>%
  pivot_longer(c(betweenness, cen_degree, closeness, degree),
               names_to = "property", values_to = "value") %>%
  mutate(
    group = factor(group, levels = c("Normal", "PoorPlaque")),
    property = factor(
      property,
      levels = c("betweenness", "cen_degree", "closeness", "degree"),
      labels = c("Betweenness", "Degree centrality", "Closeness", "Degree")
    )
  )
dual_write(node_long, "Fig_Fig5_rho080_prev50_CNS_panelB_plotdata")

sig_lab <- function(p) dplyr::case_when(
  p < 1e-3 ~ "***", p < 1e-2 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns"
)
wilcox_df <- node_long %>%
  group_by(property) %>%
  summarise(p = wilcox.test(value ~ group)$p.value,
            y = max(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(label = sig_lab(p), group1 = "Normal", group2 = "PoorPlaque",
         y.position = y * 1.10)
dual_write(wilcox_df %>% transmute(property, wilcox_p = p, label),
           "Table_Fig5_rho080_prev50_CNS_panelB_wilcoxon")

p_b <- ggplot(node_long, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.40, colour = NA, width = 0.9) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.95,
               colour = COL_INK, linewidth = 0.35) +
  geom_jitter(aes(colour = group), width = 0.09, size = 1.0,
              alpha = 0.30, show.legend = FALSE) +
  facet_wrap(~ property, scales = "free_y", nrow = 1) +
  stat_pvalue_manual(
    wilcox_df, label = "label", tip.length = 0.012,
    size = 5.2, family = FONT, bracket.size = 0.4
  ) +
  scale_fill_manual(values = c(Normal = COL_N, PoorPlaque = COL_P), name = NULL) +
  scale_colour_manual(values = c(Normal = COL_N, PoorPlaque = COL_P), guide = "none") +
  labs(y = "Node property value") +
  coord_cartesian(clip = "off") +
  theme_bw(base_family = FONT, base_size = 14) +
  theme(
    text = element_text(family = FONT, colour = COL_INK),
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = FONT, size = 14),
    axis.text = element_text(family = FONT, size = 12, colour = COL_INK),
    strip.background = element_blank(),
    strip.text = element_text(family = FONT, size = 14, colour = COL_INK),
    panel.grid.major = element_line(colour = "#EDEDED", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "#C8C8C8", linewidth = 0.45),
    legend.position = "right",
    legend.text = element_text(family = FONT, size = 12),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(6, 10, 6, 10)
  )

fig <- cowplot::plot_grid(
  p_a, p_b,
  ncol = 1,
  rel_heights = c(1.35, 1.0),
  labels = c("A", "B"),
  label_fontfamily = FONT,
  label_size = 20,
  label_fontface = "plain",
  label_colour = COL_INK,
  label_x = 0.005, label_y = 0.995,
  hjust = 0, vjust = 1
) + theme(plot.background = element_rect(fill = "white", colour = NA))

stem <- "Fig_Fig5_rho080_prev50_AB_CNS"
save_png_pdf_svg(fig, stem, w = 13.0, h = 10.4)
save_png_pdf_svg(p_a, "Fig_Fig5_rho080_prev50_network_combined_CNS", w = 12.4, h = 5.6)
save_png_pdf_svg(p_b, "Fig_Fig5_rho080_prev50_panelB_CNS", w = 11.8, h = 4.1)

# JPG from PNG (avoids R jpeg colour issues)
for (s in c(stem, "Fig_Fig5_rho080_prev50_network_combined_CNS", "Fig_Fig5_rho080_prev50_panelB_CNS")) {
  png_path <- file.path(FIG, paste0(s, ".png"))
  jpg_path <- file.path(FIG, paste0(s, ".jpg"))
  cmd <- sprintf(
    "python3 -c \"from PIL import Image; Image.open('%s').convert('RGB').save('%s', quality=95, subsampling=0)\"",
    png_path, jpg_path
  )
  system(cmd)
}
message("jpg exported from png")

cap <- paste0(
  "Fig. 5 | Prevalence-filtered MAG co-occurrence networks and node topology.\n",
  "a, Group-specific MAG–MAG Spearman co-occurrence networks for Normal and PoorPlaque ",
  "after within-group prevalence filtering (≥50% of participants; n = 10 per group), ",
  "retaining edges with |r| ≥ 0.8 and P < 0.05. Nodes are MAGs (teal, Normal; rose, PoorPlaque); ",
  "node size scales with degree. Edges are coloured by sign among resampling-robust associations ",
  "(slate blue, positive; terracotta, negative) or shown in grey if unstable ",
  "(bootstrap support < 0.8 or leave-one-out support < 0.8; B = 500). ",
  "Labels mark the five highest-degree hubs per network (species names, italic). ",
  "Normal: 85 nodes, 199 edges (189 robust, 10 unstable); ",
  "PoorPlaque: 123 nodes, 428 edges (421 robust, 7 unstable). ",
  "b, Node-level topology from the networks in a: betweenness, degree centrality, closeness and degree. ",
  "Violin plots with inset boxplots and MAG points; fill denotes group. ",
  "Asterisks indicate two-sided Wilcoxon rank-sum tests between groups ",
  "(***P < 0.001; ns, P ≥ 0.05); tests are descriptive because MAG nodes within a network ",
  "are not independent replicates. Co-occurrence denotes statistical association, ",
  "not direct microbial interaction. Source data are provided as a Source Data file.\n"
)
writeLines(cap, file.path(FIG, paste0(stem, "_CAPTION.txt")))
writeLines(cap, file.path(TAB, paste0(stem, "_CAPTION.txt")))
print(hub_key)
message("Done")
