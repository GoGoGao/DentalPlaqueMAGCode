## Packaged for revise/; source: project/review/R1.4/scripts/08_Fig2D_network_R_layouts.R
## Figures: Fig2E
## Note: network layouts
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "Fig2E_KOmetab_network_layouts.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

#!/usr/bin/env Rscript
# Fig.2D no-PI KO–metabolite network: multiple R layouts for selection (no node labels).
# Fonts: Times New Roman via cairo; SVG/PDF editable text (no showtext).

suppressPackageStartupMessages({
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
  library(graphlayouts)
})

ROOT <- file.path(.pkg_root, ""
OUT  <- file.path(ROOT, "review/R1.4/Fig2D_KOmetab_noPI")
TAB  <- file.path(OUT, "tables")
FIG  <- file.path(OUT, "figures")
CAND <- file.path(FIG, "layout_candidates_R")
dir.create(CAND, showWarnings = FALSE, recursive = TRUE)

RHO_THR <- 0.85
# Times New Roman often absent on Linux; Liberation Serif is metric-compatible
FONT <- if (any(grepl("Times New Roman", system("fc-list : family", intern = TRUE), fixed = TRUE))) {
  "Times New Roman"
} else if (any(grepl("^Liberation Serif$", system("fc-list : family", intern = TRUE)))) {
  "Liberation Serif"
} else {
  "Nimbus Roman"
}
message("font family: ", FONT)

COL <- list(
  KO = "#2C4A6E",
  Metabolite = "#D97845",
  pos = "#B83B3B",
  neg = "#3A6EA5",
  unstable = "#A0A0A0"
)

dual_write <- function(df, stem) {
  write_tsv(df, file.path(TAB, paste0(stem, ".tsv")), na = "")
  write_csv(df, file.path(TAB, paste0(stem, ".csv")), na = "")
}

save_fig <- function(p, stem, w = 7.2, h = 6.8) {
  # PDF/SVG: cairo, editable text; also PNG/JPG
  ggsave(file.path(CAND, paste0(stem, ".pdf")), p, width = w, height = h,
         device = cairo_pdf, bg = "white")
  ggsave(file.path(CAND, paste0(stem, ".svg")), p, width = w, height = h,
         device = svglite::svglite, bg = "white")
  ggsave(file.path(CAND, paste0(stem, ".png")), p, width = w, height = h,
         dpi = 300, bg = "white")
  # JPG via magick if available, else png copy route
  png_path <- file.path(CAND, paste0(stem, ".png"))
  jpg_path <- file.path(CAND, paste0(stem, ".jpg"))
  if (requireNamespace("magick", quietly = TRUE)) {
    magick::image_write(magick::image_read(png_path), path = jpg_path, format = "jpeg", quality = 95)
  } else {
    # fallback: raster via jpeg device
    jpeg(jpg_path, width = w * 300, height = h * 300, res = 300, quality = 95)
    grid::grid.raster(png::readPNG(png_path))
    dev.off()
  }
  message("saved ", stem)
}

# ---------- data ----------
edges_all <- read_tsv(file.path(TAB, "Table_Fig2D_KOmetab_edges_noPI.tsv"), show_col_types = FALSE)
sup <- read_tsv(file.path(TAB, "Table_Fig2D_noPI_edge_support_module13edges.tsv"), show_col_types = FALSE)

# normalize key like Python
norm_key_ko <- function(a, b) {
  ifelse(startsWith(as.character(a), "K"), paste(a, b, sep = "||"), paste(b, a, sep = "||"))
}
edges_all <- edges_all %>%
  mutate(key = norm_key_ko(from, to))
sup <- sup %>%
  mutate(key = paste(KO, metabolite, sep = "||"))

display <- edges_all %>%
  filter(abs(rho) >= RHO_THR) %>%
  left_join(sup %>% select(key, bootstrap_support, loo_support, robust_edge), by = "key") %>%
  mutate(robust_edge = ifelse(is.na(robust_edge), 0L, as.integer(robust_edge)))

n_rob <- sum(display$robust_edge == 1)
n_unst <- sum(display$robust_edge == 0)

g0 <- graph_from_data_frame(
  display %>% select(from, to, rho, robust_edge),
  directed = FALSE
)
# node type
V(g0)$type <- ifelse(startsWith(V(g0)$name, "K"), "KO", "Metabolite")
V(g0)$deg <- degree(g0)
E(g0)$sign <- ifelse(E(g0)$rho > 0, "pos", "neg")
E(g0)$stable <- ifelse(E(g0)$robust_edge == 1, "stable", "unstable")

subtitle <- sprintf(
  paste0("|ρ|≥%.2f: %d edges · %d nodes (stable %d; unstable %d) · ",
         "full |ρ|≥0.80: %d edges"),
  RHO_THR, ecount(g0), vcount(g0), n_rob, n_unst, nrow(edges_all)
)

# ---------- layout helpers ----------
layout_fr_tight <- function(g) {
  layout_with_fr(g, niter = 3000, grid = "nogrid")
}
layout_kk <- function(g) layout_with_kk(g)
layout_drl <- function(g) layout_with_drl(g)
layout_graphopt <- function(g) layout_with_graphopt(g, niter = 500, charge = 0.05, mass = 30, spring.length = 0.5)
layout_lgl <- function(g) layout_with_lgl(g)
layout_mds <- function(g) layout_with_mds(g)
layout_dh <- function(g) layout_with_dh(g)
layout_stress <- function(g) layout_with_stress(g, weights = NA, iter = 500)
layout_backbone <- function(g) {
  # backbone layout emphasizes communities; fallback to stress
  tryCatch(layout_as_backbone(g, keep = 0.4)$xy,
           error = function(e) layout_with_stress(g, weights = NA))
}

# Pack components: stress on each, then circle-pack centers (igraph + manual)
layout_component_pack <- function(g) {
  comps <- components(g)
  n_comp <- comps$no
  xy <- matrix(0, nrow = vcount(g), ncol = 2)
  radii <- numeric(n_comp)
  centers <- matrix(0, nrow = n_comp, ncol = 2)
  local <- vector("list", n_comp)
  sizes <- as.integer(table(comps$membership))
  ord <- order(sizes, decreasing = TRUE)
  for (ii in seq_along(ord)) {
    cid <- ord[ii]
    idx <- which(comps$membership == cid)
    sub <- induced_subgraph(g, idx)
    if (vcount(sub) == 1) {
      loc <- matrix(c(0, 0), nrow = 1)
    } else if (vcount(sub) <= 4) {
      loc <- layout_in_circle(sub)
    } else {
      loc <- tryCatch(layout_with_stress(sub, weights = NA),
                      error = function(e) layout_with_fr(sub, niter = 800))
    }
    loc <- scale(loc, center = TRUE, scale = FALSE)
    # target radius ~ sqrt(n)
    ext <- max(apply(loc, 2, function(z) diff(range(z))) / 2, 0.15)
    target_r <- 0.25 * sqrt(vcount(sub)) + if (ii == 1) 0.8 else 0.15
    if (ii == 1) target_r <- target_r * 1.6
    loc <- loc * (target_r / ext)
    local[[cid]] <- list(idx = idx, loc = loc, r = target_r)
    radii[cid] <- target_r
  }
  # greedy pack largest first
  placed <- list()
  for (cid in ord) {
    r <- radii[cid] + 0.05
    if (length(placed) == 0) {
      centers[cid, ] <- c(0, 0)
      placed[[length(placed) + 1]] <- c(0, 0, r)
      next
    }
    best <- NULL
    best_s <- Inf
    for (pl in placed) {
      n_ang <- max(12L, as.integer(16 * (pl[3] + r)))
      for (k in seq_len(n_ang)) {
        ang <- 2 * pi * (k - 1) / n_ang
        x <- pl[1] + (pl[3] + r) * cos(ang)
        y <- pl[2] + (pl[3] + r) * sin(ang)
        ok <- TRUE
        for (q in placed) {
          if ((x - q[1])^2 + (y - q[2])^2 < (q[3] + r - 1e-6)^2) {
            ok <- FALSE
            break
          }
        }
        if (!ok) next
        sc <- x * x + y * y
        if (sc < best_s) {
          best_s <- sc
          best <- c(x, y)
        }
      }
    }
    if (is.null(best)) {
      xmax <- max(vapply(placed, function(q) q[1] + q[3], 0))
      best <- c(xmax + r, 0)
    }
    centers[cid, ] <- best
    placed[[length(placed) + 1]] <- c(best[1], best[2], r)
  }
  for (cid in seq_len(n_comp)) {
    idx <- local[[cid]]$idx
    xy[idx, ] <- sweep(local[[cid]]$loc, 2, centers[cid, ], "+")
  }
  # square max-fill
  xy <- scale(xy, center = TRUE, scale = FALSE)
  m <- max(abs(xy))
  if (m > 0) xy <- xy / m
  xy
}

layouts <- list(
  stress        = layout_stress,
  fr_tight      = layout_fr_tight,
  kk            = layout_kk,
  graphopt      = layout_graphopt,
  drl           = layout_drl,
  mds           = layout_mds,
  backbone      = layout_backbone,
  component_pack = layout_component_pack
)

# ---------- plot one layout ----------
plot_one <- function(g, xy, layout_name) {
  # attach coords
  V(g)$x <- xy[, 1]
  V(g)$y <- xy[, 2]
  tg <- as_tbl_graph(g)

  # edge aesthetic data
  edf <- igraph::as_data_frame(g, what = "edges") %>%
    mutate(
      edge_col = case_when(
        stable == "unstable" ~ COL$unstable,
        sign == "pos" ~ COL$pos,
        TRUE ~ COL$neg
      ),
      edge_lty = ifelse(stable == "unstable", "22", "solid"),
      edge_alpha = ifelse(stable == "unstable", 0.85, 0.32),
      edge_width = ifelse(
        stable == "unstable", 0.55,
        0.25 + 1.2 * pmin(1, (abs(rho) - RHO_THR) / 0.15)
      )
    )

  ndf <- igraph::as_data_frame(g, what = "vertices") %>%
    mutate(
      node_col = ifelse(type == "KO", COL$KO, COL$Metabolite),
      node_size = 1.6 + 0.22 * deg,
      node_size = ifelse(deg >= 8, pmax(node_size, 5.5), node_size)
    )

  # rebuild tbl with layout for ggraph
  tg <- tg %>%
    activate(nodes) %>%
    mutate(x = xy[, 1], y = xy[, 2],
           node_col = ndf$node_col, node_size = ndf$node_size)

  # Manual ggplot for full control (avoid ggraph edge aes quirks with mixed lty)
  # Use ggraph create_layout then ggplot layers
  lay <- create_layout(tg, layout = "manual", x = xy[, 1], y = xy[, 2])

  p <- ggraph(lay) +
    geom_edge_link(
      aes(edge_colour = I(edf$edge_col),
          edge_width = I(edf$edge_width),
          edge_alpha = I(edf$edge_alpha),
          edge_linetype = I(edf$edge_lty)),
      lineend = "round"
    ) +
    geom_node_point(
      aes(fill = type, size = node_size),
      shape = 21, colour = "white", stroke = 0.15
    ) +
    scale_fill_manual(values = c(KO = COL$KO, Metabolite = COL$Metabolite), name = NULL) +
    scale_size_identity() +
    coord_fixed() +
    labs(
      title = "Robust_DA KO–metabolite associations (no Plaque Index)",
      subtitle = paste0(subtitle, "\nLayout: ", layout_name),
      caption = NULL
    ) +
    theme_void(base_family = FONT, base_size = 11) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, family = FONT),
      plot.subtitle = element_text(hjust = 0.5, size = 7.5, colour = "#444444",
                                   family = FONT, lineheight = 1.15),
      legend.position = "bottom",
      legend.text = element_text(size = 8, family = FONT),
      plot.margin = margin(6, 6, 6, 6),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.caption = element_text(hjust = 0.5, size = 7.5, family = FONT, colour = "#333333")
    ) +
    guides(fill = guide_legend(override.aes = list(size = 3, stroke = 0.2))) +
    labs(caption = "Edges: solid red = stable +ρ; solid blue = stable −ρ; dashed grey = unstable")
  p
}

# ---------- run all layouts ----------
message("=== R multi-layout Fig.2D |rho|>=", RHO_THR, " ===")
set.seed(42)
all_layouts_df <- list()

for (nm in names(layouts)) {
  message("layout: ", nm)
  xy <- tryCatch(layouts[[nm]](g0), error = function(e) {
    message("  FAIL ", nm, ": ", conditionMessage(e))
    NULL
  })
  if (is.null(xy)) next
  if (nrow(xy) != vcount(g0)) {
    message("  bad dims for ", nm)
    next
  }
  # center
  xy <- scale(xy, center = TRUE, scale = FALSE)
  colnames(xy) <- c("x", "y")

  p <- plot_one(g0, xy, nm)
  stem <- paste0("Fig2D_KOmetab_network_noPI__", nm)
  save_fig(p, stem, w = 7.0, h = 6.6)

  df <- data.frame(
    id = V(g0)$name,
    x = xy[, 1],
    y = xy[, 2],
    type = V(g0)$type,
    degree = V(g0)$deg,
    layout = nm,
    stringsAsFactors = FALSE
  )
  dual_write(df, paste0("Fig2D_KOmetab_network_layout_R__", nm))
  all_layouts_df[[nm]] <- df
}

# Combined overview: 2x4 facet of point+edge miniatures is heavy; instead save index table
idx <- data.frame(
  layout = names(all_layouts_df),
  n_nodes = vcount(g0),
  n_edges = ecount(g0),
  file_stem = paste0("Fig2D_KOmetab_network_noPI__", names(all_layouts_df)),
  stringsAsFactors = FALSE
)
dual_write(idx, "Fig2D_KOmetab_network_R_layout_index")

# Default pick: stress (usually clean for association nets) → overwrite main figure
pick <- if ("stress" %in% names(all_layouts_df)) "stress" else names(all_layouts_df)[[1]]
message("default pick → ", pick)
src <- file.path(CAND, paste0("Fig2D_KOmetab_network_noPI__", pick))
for (ext in c("pdf", "svg", "png", "jpg")) {
  f <- paste0(src, ".", ext)
  if (file.exists(f)) {
    file.copy(f, file.path(FIG, paste0("Fig2D_KOmetab_network_noPI.", ext)), overwrite = TRUE)
  }
}
# also copy layout coords of pick to main layout name
dual_write(all_layouts_df[[pick]] %>% select(-layout), "Fig2D_KOmetab_network_layout")
dual_write(display, "Fig2D_KOmetab_network_display_edges")

# Write a short CHOICE note
writeLines(c(
  "Fig.2D R layout candidates (no node labels)",
  paste0("Directory: ", CAND),
  paste0("Default copied to main figure: ", pick),
  "Layouts:",
  paste0("  - ", idx$layout, " → ", idx$file_stem, ".jpg"),
  "",
  "Reply with preferred layout name to set as main Fig2D_KOmetab_network_noPI.*"
), file.path(CAND, "README_choose_layout.txt"))

message("DONE; candidates in ", CAND)
message("main figure updated from layout: ", pick)
