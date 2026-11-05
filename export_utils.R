## Packaged for revise/; source: project/01this_study_dataset/sgb/metab_corr/paper_result/_shared/export_utils.R
## Figures: shared
## Note: figure export helpers
args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else NA_character_
.pkg_root <- Sys.getenv("PLAQUE_PKG_ROOT", unset = "")
if (!nzchar(.pkg_root)) {
  if (!is.na(this_file)) .pkg_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = FALSE)
  else .pkg_root <- normalizePath("..", mustWork = FALSE)
}
SHARED <- file.path(.pkg_root, "shared_input")
OUTDIR <- file.path(.pkg_root, "output", tools::file_path_sans_ext(basename(if (!is.na(this_file)) this_file else "export_utils.R")))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
## --- original script below (host paths may remain; override with SHARED/OUTDIR as needed) ---

# Shared export helpers for paper_result SCI figures
# PDF / SVG / PNG / JPG + optional TSV handled by callers

save_gg <- function(p, stem, w = 8, h = 6, dpi = 300) {
  dir.create(dirname(stem), showWarnings = FALSE, recursive = TRUE)
  for (ext in c("pdf", "svg", "png", "jpg")) {
    f <- paste0(stem, ".", ext)
    tryCatch({
      if (ext == "pdf") {
        suppressMessages(ggplot2::ggsave(f, plot = p, width = w, height = h, dpi = dpi,
                         device = grDevices::cairo_pdf, bg = "white"))
      } else if (ext == "jpg") {
        grDevices::jpeg(f, width = w, height = h, units = "in", res = dpi,
                        quality = 95, bg = "white", type = "cairo")
        print(p)
        invisible(grDevices::dev.off())
      } else {
        suppressMessages(ggplot2::ggsave(f, plot = p, width = w, height = h, dpi = dpi,
                         device = ext, bg = "white"))
      }
      message("[OK] ", basename(f))
    }, error = function(e) message("[ERR] ", basename(f), ": ", e$message))
  }
}
