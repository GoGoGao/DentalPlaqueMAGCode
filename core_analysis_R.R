#!/usr/bin/env Rscript
# ══════════════════════════════════════════════════════════════════════════════
#  Core Analysis Functions for Oral Microbiome Metagenomic Study
# ══════════════════════════════════════════════════════════════════════════════
#
#  Reusable statistical functions for:
#    - MAG-level gene set comparison (Venn / Fisher exact test)
#    - KO differential prevalence (Fisher exact test, BH correction)
#    - VF / ARG burden comparison (Wilcoxon rank-sum)
#    - Functional profile ordination (PCoA + PERMANOVA)
#    - KO–clinical indicator correlation (Spearman / point-biserial)
#    - KEGG pathway enrichment (clusterProfiler)
#
#  Dependencies:
#    tidyverse, vegan, ComplexHeatmap, circlize, clusterProfiler, Cairo
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(clusterProfiler)
})


# ══════════════════════════════════════════════════════════════════════════════
#  1. Gene Set Extraction & Venn Comparison
# ══════════════════════════════════════════════════════════════════════════════

#' Extract union of detected gene features from a binary/count matrix
#'
#' @param mat Matrix with rows = features, cols = genomes/MAGs.
#' @param genome_ids Character vector of genome IDs to include.
#' @return Character vector of feature names detected in >= 1 genome.
extract_gene_set <- function(mat, genome_ids) {
  avail <- intersect(genome_ids, colnames(mat))
  if (length(avail) == 0L) return(character(0))
  sub <- mat[, avail, drop = FALSE]
  rownames(sub)[rowSums(sub > 0) > 0]
}


#' Compute Venn set statistics between two groups
#'
#' @param set_a Character vector of features in group A.
#' @param set_b Character vector of features in group B.
#' @return List with shared, only_A, only_B, and count summary.
venn_set_stats <- function(set_a, set_b) {
  shared <- intersect(set_a, set_b)
  only_a <- setdiff(set_a, set_b)
  only_b <- setdiff(set_b, set_a)

  list(
    shared = shared,
    only_A = only_a,
    only_B = only_b,
    counts = data.frame(
      Category = c("A_only", "Shared", "B_only"),
      Count    = c(length(only_a), length(shared), length(only_b))
    )
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  2. Fisher Exact Test for Differential Feature Prevalence
# ══════════════════════════════════════════════════════════════════════════════

#' Per-feature Fisher exact test comparing binary prevalence between groups
#'
#' Constructs a 2x2 contingency table (present/absent × group) for each row
#' of a binary matrix, applies Fisher's exact test, and corrects with BH.
#'
#' @param bin_mat Binary matrix (features × genomes), 1 = detected.
#' @param ids_a Character vector of genome IDs in group A.
#' @param ids_b Character vector of genome IDs in group B.
#' @param group_a_name Label for group A (default "GroupA").
#' @param group_b_name Label for group B (default "GroupB").
#' @return data.frame with per-feature prevalence, odds ratio, p, padj, enrichment.
fisher_binary_test <- function(bin_mat, ids_a, ids_b,
                               group_a_name = "GroupA",
                               group_b_name = "GroupB") {
  mA <- intersect(ids_a, colnames(bin_mat))
  mB <- intersect(ids_b, colnames(bin_mat))
  nA <- length(mA)
  nB <- length(mB)

  feats <- rownames(bin_mat)
  res <- lapply(feats, function(f) {
    a <- sum(bin_mat[f, mA, drop = TRUE] > 0)
    b <- nA - a
    cc <- sum(bin_mat[f, mB, drop = TRUE] > 0)
    d <- nB - cc
    tab <- matrix(c(a, b, cc, d), nrow = 2, byrow = TRUE)
    ft <- suppressWarnings(fisher.test(tab, alternative = "two.sided"))
    data.frame(
      feature            = f,
      n_present_A        = a,
      n_present_B        = cc,
      prevalence_A       = a / max(nA, 1L),
      prevalence_B       = cc / max(nB, 1L),
      delta_prevalence   = cc / max(nB, 1L) - a / max(nA, 1L),
      odds_ratio         = unname(ft$estimate),
      p_value            = ft$p.value,
      stringsAsFactors   = FALSE
    )
  })

  df <- bind_rows(res) %>%
    mutate(
      padj = p.adjust(p_value, method = "BH"),
      enriched_in = case_when(
        padj < 0.05 & prevalence_A > prevalence_B ~ group_a_name,
        padj < 0.05 & prevalence_B > prevalence_A ~ group_b_name,
        TRUE ~ "NS"
      )
    ) %>%
    arrange(p_value)

  colnames(df) <- gsub("_A$", paste0("_", group_a_name), colnames(df))
  colnames(df) <- gsub("_B$", paste0("_", group_b_name), colnames(df))

  df
}


# ══════════════════════════════════════════════════════════════════════════════
#  3. Per-Genome Feature Burden Comparison (Wilcoxon)
# ══════════════════════════════════════════════════════════════════════════════

#' Compare per-genome feature counts between two groups (Wilcoxon rank-sum)
#'
#' @param bin_mat Binary feature matrix (features × genomes).
#' @param group_map Named vector: genome_id -> group_label.
#' @param group_a Label for group A.
#' @param group_b Label for group B.
#' @return List with test_stat, p_value, and per-genome count table.
compare_feature_burden <- function(bin_mat, group_map,
                                   group_a = "GroupA",
                                   group_b = "GroupB") {
  counts <- colSums(bin_mat > 0)
  count_df <- tibble(
    genome = names(counts),
    count  = unname(counts),
    group  = group_map[names(counts)]
  ) %>%
    filter(group %in% c(group_a, group_b))

  x <- count_df$count[count_df$group == group_a]
  y <- count_df$count[count_df$group == group_b]
  wt <- suppressWarnings(wilcox.test(x, y))

  list(
    p_value  = wt$p.value,
    W_stat   = unname(wt$statistic),
    count_df = count_df
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  4. PCoA + PERMANOVA on Binary Functional Profiles
# ══════════════════════════════════════════════════════════════════════════════

#' PCoA ordination with PERMANOVA on presence/absence functional profiles
#'
#' @param pa_mat Presence/absence matrix (genomes × features).
#' @param group_map Named vector: genome_id -> group_label.
#' @param dist_method Distance method for vegan::vegdist (default "jaccard").
#' @param n_perm Number of permutations for PERMANOVA.
#' @return List with coordinates, variance explained, PERMANOVA results.
run_pcoa_permanova <- function(pa_mat, group_map,
                               dist_method = "jaccard",
                               n_perm = 999) {
  genomes <- intersect(rownames(pa_mat), names(group_map))
  pa_sub <- pa_mat[genomes, , drop = FALSE]
  group_df <- data.frame(
    genome = genomes,
    group  = factor(group_map[genomes])
  )

  dist_mat <- vegdist(pa_sub, method = dist_method, binary = TRUE)
  pcoa_res <- cmdscale(dist_mat, k = min(nrow(pa_sub) - 1, 10), eig = TRUE)

  eig_vals <- pcoa_res$eig
  eig_vals[eig_vals < 0] <- 0
  var_exp <- eig_vals / sum(eig_vals) * 100

  perm_res <- adonis2(dist_mat ~ group, data = group_df, permutations = n_perm)

  list(
    coords       = data.frame(
      genome = genomes,
      PCoA1  = pcoa_res$points[, 1],
      PCoA2  = pcoa_res$points[, 2],
      group  = group_df$group
    ),
    var_explained = c(PC1 = var_exp[1], PC2 = var_exp[2]),
    permanova     = list(
      R2      = perm_res$R2[1],
      F_stat  = perm_res$F[1],
      p_value = perm_res$`Pr(>F)`[1]
    )
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  5. KO–Clinical Indicator Correlation
# ══════════════════════════════════════════════════════════════════════════════

#' Correlation analysis between functional features and clinical indicators
#'
#' Uses Spearman rank correlation for continuous clinical variables and
#' point-biserial (Pearson on 0/1 coded) for binary categorical variables.
#' P-values are corrected by Benjamini-Hochberg.
#'
#' @param feature_mat Samples × features numeric matrix.
#' @param clinical_df Samples × clinical variables data.frame.
#' @param continuous_vars Character vector of continuous variable names.
#' @param categorical_vars Character vector of binary variable names (coded 0/1).
#' @return data.frame with feature, clinical_variable, method, correlation,
#'         CI_lower, CI_upper, p_value, FDR.
correlate_features_clinical <- function(feature_mat, clinical_df,
                                        continuous_vars = NULL,
                                        categorical_vars = NULL) {
  all_vars <- c(continuous_vars, categorical_vars)
  feature_ids <- colnames(feature_mat)

  results <- list()
  idx <- 1L
  for (feat in feature_ids) {
    y <- as.numeric(feature_mat[, feat])
    for (clin in all_vars) {
      x <- as.numeric(clinical_df[[clin]])
      ok <- complete.cases(x, y)
      if (sum(ok) < 5L) next
      xx <- x[ok]; yy <- y[ok]

      if (clin %in% categorical_vars) {
        ct <- cor.test(xx, yy, method = "pearson")
        method <- "point-biserial"
      } else {
        ct <- cor.test(xx, yy, method = "spearman", exact = FALSE)
        method <- "Spearman"
      }

      # Fisher-z 95% CI
      n <- sum(ok)
      r_val <- ct$estimate
      z <- atanh(r_val)
      z_se <- 1 / sqrt(n - 3)
      ci_lo <- tanh(z - 1.96 * z_se)
      ci_hi <- tanh(z + 1.96 * z_se)

      results[[idx]] <- data.frame(
        feature           = feat,
        clinical_variable = clin,
        method            = method,
        correlation       = as.numeric(r_val),
        CI_lower          = ci_lo,
        CI_upper          = ci_hi,
        p_value           = ct$p.value,
        n                 = n,
        stringsAsFactors  = FALSE
      )
      idx <- idx + 1L
    }
  }

  df <- bind_rows(results)
  if (nrow(df) > 0) {
    df$FDR <- p.adjust(df$p_value, method = "BH")
  }

  df
}


# ══════════════════════════════════════════════════════════════════════════════
#  6. KEGG Pathway Enrichment
# ══════════════════════════════════════════════════════════════════════════════

#' KEGG pathway enrichment for a list of KO identifiers
#'
#' @param ko_ids Character vector of KO identifiers (e.g., "K00001").
#' @param padj_cutoff Adjusted p-value threshold.
#' @return data.frame of enrichment results, or NULL if none significant.
run_kegg_enrichment <- function(ko_ids, padj_cutoff = 0.05) {
  if (length(ko_ids) == 0L) {
    message("[WARN] Empty KO list, skipping enrichment.")
    return(NULL)
  }

  res <- tryCatch(
    enrichKEGG(
      gene         = ko_ids,
      organism     = "ko",
      keyType      = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    ),
    error = function(e) {
      message("[WARN] enrichKEGG failed: ", e$message)
      NULL
    }
  )

  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)

  df <- as.data.frame(res) %>%
    mutate(
      GeneRatio_num = map_dbl(GeneRatio, function(x) {
        sp <- strsplit(as.character(x), "/")[[1]]
        if (length(sp) == 2) as.numeric(sp[1]) / as.numeric(sp[2]) else NA
      }),
      neg_log10_padj = -log10(p.adjust)
    )

  df
}


#' Parse KO IDs from annotation strings (e.g., "K00001:gene_name" -> "K00001")
#'
#' @param x Character vector of annotated KO strings.
#' @return Character vector of clean KO identifiers.
parse_ko_id <- function(x) {
  sub(":.*", "", trimws(x))
}


# ══════════════════════════════════════════════════════════════════════════════
#  7. Utility Functions
# ══════════════════════════════════════════════════════════════════════════════

#' Build a presence/absence matrix from a count matrix (transposed)
#'
#' @param mat Count matrix (features × genomes).
#' @param genome_ids Character vector of genome IDs to include.
#' @return Binary matrix (genomes × features) suitable for PCoA.
build_pa_matrix <- function(mat, genome_ids) {
  avail <- intersect(genome_ids, colnames(mat))
  sub <- mat[, avail, drop = FALSE]
  sub <- sub[rowSums(sub) > 0, , drop = FALSE]
  t((sub > 0) * 1L)
}


#' Save a ggplot object in PDF, SVG, and PNG formats
#'
#' @param plot_obj ggplot object.
#' @param prefix Output file path prefix (without extension).
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @param dpi Resolution for PNG output.
save_plot_multi <- function(plot_obj, prefix, width, height, dpi = 300) {
  for (ext in c("pdf", "svg", "png")) {
    fpath <- paste0(prefix, ".", ext)
    if (ext == "png") {
      ggsave(fpath, plot_obj, width = width, height = height, dpi = dpi)
    } else {
      ggsave(fpath, plot_obj, width = width, height = height, device = ext)
    }
  }
}
