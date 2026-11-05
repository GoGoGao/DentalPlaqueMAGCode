## Packaged for revise/; source: project/review/R1.4/Fig5_MAG_network_rho080_prev50/scripts/03_plot_Fig5B_node_topology.R
## Figures: Fig5B
## Note: node topology violins
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig5B_node_topology.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# Figure 5B-style node topology comparison for rho080_prev50 MAG networks.
# Violin + box + jitter; Wilcoxon P (descriptive; nodes are not independent replicates).

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
})

OUT <- file.path(.pkg_root, "review/R1.4/Fig5_MAG_network_rho080_prev50"
TAB <- file.path(OUT, "tables")
FIG <- file.path(OUT, "figures")
NET_FIG <- file.path(FIG, "Fig_Fig5_rho080_prev50_network_combined.pdf")

group_colors <- c(Normal = "#00A087", PoorPlaque = "#E64B35")

dual_write <- function(df, stem) {
  write_tsv(df, file.path(TAB, paste0(stem, ".tsv")))
  write_csv(df, file.path(TAB, paste0(stem, ".csv")))
}

nodes <- bind_rows(
  read_csv(file.path(TAB, "Table_Fig5_rho080_prev50_nodes_Normal.csv"), show_col_types = FALSE),
  read_csv(file.path(TAB, "Table_Fig5_rho080_prev50_nodes_PoorPlaque.csv"), show_col_types = FALSE)
)

# Degree centrality = degree / (n-1) within each network (standard definition).
# Keep raw degree as a separate panel to match Fig. 5 legend (betweenness, degree centrality, closeness, degree).
node_long <- nodes %>%
  group_by(group) %>%
  mutate(cen_degree = degree / pmax(n() - 1, 1)) %>%
  ungroup() %>%
  select(MAG, group, betweenness, cen_degree, closeness, degree) %>%
  pivot_longer(
    cols = c(betweenness, cen_degree, closeness, degree),
    names_to = "property",
    values_to = "value"
  ) %>%
  mutate(
    group = factor(group, levels = c("Normal", "PoorPlaque")),
    property = factor(
      property,
      levels = c("betweenness", "cen_degree", "closeness", "degree"),
      labels = c("Betweenness", "Degree centrality", "Closeness", "Degree")
    )
  )

dual_write(node_long, "Fig_Fig5_rho080_prev50_panelB_node_topology_plotdata")
dual_write(node_long, "Table_Fig5_rho080_prev50_panelB_node_topology_long")

# Summary + Wilcoxon (reported as in original Fig.5B; interpretive caution in caption)
wilcox_rows <- node_long %>%
  group_by(property) %>%
  summarise(
    n_Normal = sum(group == "Normal"),
    n_PoorPlaque = sum(group == "PoorPlaque"),
    median_Normal = median(value[group == "Normal"]),
    median_PoorPlaque = median(value[group == "PoorPlaque"]),
    mean_Normal = mean(value[group == "Normal"]),
    mean_PoorPlaque = mean(value[group == "PoorPlaque"]),
    wilcox_p = wilcox.test(value ~ group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(wilcox_p_label = sprintf("P = %.2g", wilcox_p))
dual_write(wilcox_rows, "Table_Fig5_rho080_prev50_panelB_wilcoxon")
print(wilcox_rows)

nN <- n_distinct(nodes$MAG[nodes$group == "Normal"])
nP <- n_distinct(nodes$MAG[nodes$group == "PoorPlaque"])

p_b <- ggplot(node_long, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.55, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.9, color = "black") +
  geom_jitter(
    aes(color = group),
    width = 0.12, alpha = 0.25, size = 0.55, show.legend = FALSE
  ) +
  facet_wrap(~ property, scales = "free_y", nrow = 1) +
  stat_compare_means(
    comparisons = list(c("Normal", "PoorPlaque")),
    method = "wilcox.test",
    label = "p.format",
    hide.ns = FALSE,
    vjust = -0.2,
    tip.length = 0.015,
    size = 3.1
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  coord_cartesian(clip = "off") +
  theme_bw(base_family = "serif", base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 8.5, colour = "grey35"),
    plot.margin = margin(10, 10, 8, 10)
  ) +
  labs(
    y = "Node property value",
    fill = "Group",
    title = "Node-level topology of |r| >= 0.8 MAG networks (prev >= 50%)",
    subtitle = sprintf(
      "Violin + boxplots; points = individual MAGs. Teal, Normal (n = %d); red, PoorPlaque (n = %d). Wilcoxon rank-sum P shown.",
      nN, nP
    )
  )

stem <- "Fig_Fig5_rho080_prev50_panelB_node_topology"
for (ext in c("pdf", "svg", "png", "jpg")) {
  ggsave(file.path(FIG, paste0(stem, ".", ext)), p_b,
         width = 10.2, height = 4.6, dpi = 300, bg = "white")
}
message("saved ", stem)

# Optional A+B composite if network PDF exists
if (file.exists(NET_FIG) && requireNamespace("patchwork", quietly = TRUE)) {
  # Rebuild lightweight network thumbnails from existing node/edge plotdata is heavy;
  # instead stack panelB under a note that (A) is the combined network figure.
  # Create labeled panelB as standalone Fig5B; write AB caption pointing to both files.
  message("panel B standalone; combine with (A) via caption / downstream layout")
}

cap <- paste0(
  "Figure 5B. Node-level topology of prevalence-filtered MAG co-occurrence networks (|r| ≥ 0.8).\n\n",
  "Networks correspond to Fig. 5A rebuilt after within-group prevalence ≥ 50% and |r| ≥ 0.8 ",
  "(all qualifying edges retained). Panels show betweenness, degree centrality (degree/(n−1)), ",
  "closeness and degree. Violin + boxplots; points, individual MAGs. Teal, Normal (n = ", nN,
  "); red, PoorPlaque (n = ", nP, "). Wilcoxon rank-sum P values are shown for visual comparison ",
  "with the original Figure 5 style; because MAG nodes within a network are not independent ",
  "biological replicates, these P values are descriptive and are not used for participant-level inference.\n"
)
writeLines(cap, file.path(FIG, paste0(stem, "_CAPTION.txt")))
writeLines(cap, file.path(TAB, paste0(stem, "_CAPTION.txt")))

# Also try patchwork AB: read network combined png as grob is awkward in R without magick.
# Produce a simple two-row PDF using cowplot if network jpg available.
net_jpg <- file.path(FIG, "Fig_Fig5_rho080_prev50_network_combined.jpg")
if (file.exists(net_jpg) && requireNamespace("cowplot", quietly = TRUE)) {
  p_a <- cowplot::ggdraw() + cowplot::draw_image(net_jpg)
  p_ab <- cowplot::plot_grid(
    p_a, p_b,
    ncol = 1, rel_heights = c(1.15, 0.85),
    labels = c("A", "B"), label_fontfamily = "serif", label_size = 14
  )
  ab_stem <- "Fig_Fig5_rho080_prev50_AB"
  for (ext in c("pdf", "png", "jpg")) {
    cowplot::save_plot(
      file.path(FIG, paste0(ab_stem, ".", ext)),
      p_ab, base_width = 12.5, base_height = 10.5, dpi = 300
    )
  }
  # svg via ggsave if possible
  ggsave(file.path(FIG, paste0(ab_stem, ".svg")), p_ab, width = 12.5, height = 10.5, bg = "white")
  writeLines(
    paste0(
      "Figure 5. Prevalence-filtered MAG co-occurrence networks (|r| ≥ 0.8).\n\n",
      "(A) MAG co-occurrence networks built separately in Normal (left) and PoorPlaque (right) ",
      "from MAG relative abundances after within-group prevalence ≥ 50%, retaining edges with ",
      "Spearman |r| ≥ 0.8 (P < 0.05 at the primary |r| ≥ 0.5 screen). Nodes, MAGs (size = degree); ",
      "edges, positive (blue) or negative (red) correlations. Labels, top-degree hubs.\n",
      "(B) Node-level topology of the networks in (A) (betweenness, degree centrality, closeness and degree). ",
      "Violin + boxplots; points, individual MAGs. Teal, Normal (n = ", nN, "); red, PoorPlaque (n = ", nP,
      "). Wilcoxon rank-sum P values shown (descriptive only; nodes are not independent replicates).\n"
    ),
    file.path(FIG, paste0(ab_stem, "_CAPTION.txt"))
  )
  message("saved ", ab_stem)
}

message("Done")
