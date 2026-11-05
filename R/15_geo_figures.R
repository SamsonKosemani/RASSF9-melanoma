# ==============================================================================
# 15_geo_figures.R
#
# Rebuild every GEO validation figure and table, ONE GENE PER FILE.
#
# RASSF9 and RASSF10 never share a figure or a table, so either can be promoted
# to a main figure or demoted to supplementary without any rework. Same rule as
# the TCGA outputs.
#
# This script RECOMPUTES NOTHING. It reads the result CSVs written by scripts 13
# and 14 and redraws from them, so the figures cannot drift from the numbers
# that were tested. Run 13 and 14 first.
#
# DESIGN NOTES
#   Colour. The blue/vermillion pair (#0072B2 / #D55E00) was checked with a
#   six-point palette validator: lightness band, chroma floor, colour-vision
#   separation, normal-vision floor, and contrast against the page. It passes
#   all five, with protanope dE 21.9. The teal/coral pair used previously
#   failed the chroma floor and scored dE 7.1 under protanopia -- near-identical
#   greys for a red-green colourblind reader.
#
#   Form follows the data's job. Two unordered groups (nevus vs melanoma) get a
#   categorical pair. Progression stages are ORDERED, so they get a sequential
#   single-hue ramp, which survives greyscale printing and encodes the ordering
#   in lightness rather than hue.
#
#   Identity is never colour alone: every group is labelled on the axis, and the
#   effect-size panel labels each cohort directly. A reader who sees no colour
#   at all loses nothing.
#
# Output: figures/geo/geo_rassf9_*.{pdf,png}
#         figures/geo/geo_rassf10_*.{pdf,png}
#         results/geo/table_geo_rassf9_*.csv, table_geo_rassf10_*.csv
#         results/logs/15_geo_figures.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(patchwork)
  library(scales)
})

LOG <- file.path(PATHS$logs, "15_geo_figures.log")
cat("GEO figures rebuild —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

RES_GEO <- file.path(PATHS$results, "geo")
FIG_GEO <- file.path(PATHS$figures, "geo")
dir.create(RES_GEO, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_GEO, recursive = TRUE, showWarnings = FALSE)

GENES <- c("RASSF9", "RASSF10")

# ==============================================================================
# 0. Clear the old mixed-gene outputs
# ==============================================================================
hdr("0. Removing superseded figures")

old <- list.files(FIG_GEO, full.names = TRUE)
if (length(old)) {
  file.remove(old)
  log_line("  removed ", length(old), " old figure file(s):", logfile = LOG)
  for (o in basename(old)) log_line("    ", o, logfile = LOG)
} else {
  log_line("  none present", logfile = LOG)
}

# The analysis CSVs from 13 and 14 are the inputs here, so they stay. Only the
# per-gene table copies this script writes are refreshed.
old_tab <- list.files(RES_GEO, pattern = "^table_geo_", full.names = TRUE)
if (length(old_tab)) {
  file.remove(old_tab)
  log_line("  removed ", length(old_tab), " old per-gene table(s)", logfile = LOG)
}

# ==============================================================================
# 1. Load results
# ==============================================================================
hdr("1. Inputs")

need <- function(p) {
  if (!file.exists(p)) stop("Missing input: ", p,
                            "\nRun scripts 13 and 14 first.")
  p
}

vals  <- readr::read_csv(need(file.path(RES_GEO, "nevus_vs_melanoma_values.csv")))
tests <- readr::read_csv(need(file.path(RES_GEO, "nevus_vs_melanoma_tests.csv")))
prog  <- readr::read_csv(need(file.path(RES_GEO, "gse4587_values.csv")))

log_line("  nevus/melanoma values: ", nrow(vals), " samples", logfile = LOG)
log_line("  nevus/melanoma tests:  ", nrow(tests), " rows", logfile = LOG)
log_line("  GSE4587 values:        ", nrow(prog), " rows", logfile = LOG)

vals$group <- factor(vals$group, levels = c("Nevus", "Melanoma"))

# Effect-size CI for the forest panel, recomputed from the stored d and n only
# (no re-testing -- the same arithmetic script 14 already logged).
tests <- tests |>
  dplyr::mutate(label = paste0(dataset, "\n(n = ", n_nevus, " / ", n_melanoma, ")"))

# ==============================================================================
# 2. Panel builders
# ==============================================================================

# --- lesion comparison, one cohort -------------------------------------------
panel_lesion <- function(gene, ds) {
  if (!gene %in% names(vals)) return(NULL)
  pd <- vals |> dplyr::filter(dataset == ds, !is.na(.data[[gene]]))
  if (!nrow(pd)) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, size = 2.8, colour = PAL$grey_pt,
                 lineheight = 1.2,
                 label = paste0(gene, " is not present\nin the ", ds,
                                " matrices")) +
        labs(subtitle = ds) +
        theme_void() +
        theme(plot.subtitle = element_text(size = FIG$base_size,
                                           colour = PAL$grey_pt, hjust = 0.5))
    )
  }

  r <- tests |> dplyr::filter(dataset == ds, gene == !!gene)
  rng <- range(pd[[gene]], na.rm = TRUE)
  pad <- diff(rng) * 0.14
  if (!is.finite(pad) || pad == 0) pad <- 1

  n_lab <- pd |> dplyr::count(group) |>
    dplyr::mutate(lab = paste0(group, "\nn = ", n))
  labs_v <- setNames(n_lab$lab, n_lab$group)

  p <- ggplot(pd, aes(x = group, y = .data[[gene]], fill = group)) +
    # Distribution shape, then the raw points, then the summary on top: the
    # reader sees the actual sample, not only a box.
    geom_violin(trim = TRUE, linewidth = 0.3, colour = NA, alpha = 0.22,
                scale = "width") +
    geom_jitter(width = 0.11, height = 0, size = 0.7, alpha = 0.55,
                colour = PAL$grey_pt, stroke = 0) +
    geom_boxplot(width = 0.19, outlier.shape = NA, linewidth = 0.45,
                 colour = PAL$ink, fill = "white", alpha = 0.92) +
    stat_summary(fun = median, geom = "point", shape = 21, size = 1.9,
                 stroke = 0.45, colour = PAL$ink,
                 aes(fill = group), show.legend = FALSE) +
    scale_fill_manual(values = PAL_LESION, guide = "none") +
    scale_x_discrete(labels = labs_v) +
    labs(x = NULL, y = bquote(bold(italic(.(gene))~"expression")),
         subtitle = ds) +
    theme_pub()

  if (nrow(r)) {
    ybar <- rng[2] + pad * 0.45
    p <- p +
      annotate("segment", x = 1, xend = 2, y = ybar, yend = ybar,
               linewidth = 0.45, colour = PAL$ink) +
      annotate("segment", x = 1, xend = 1, y = ybar, yend = ybar - pad * 0.12,
               linewidth = 0.45, colour = PAL$ink) +
      annotate("segment", x = 2, xend = 2, y = ybar, yend = ybar - pad * 0.12,
               linewidth = 0.45, colour = PAL$ink) +
      annotate("text", x = 1.5, y = ybar + pad * 0.1,
               label = fmt_p_plot(r$p_value), parse = TRUE, size = 2.9,
               vjust = 0, colour = PAL$ink) +
      coord_cartesian(ylim = c(rng[1] - pad * 0.25, rng[2] + pad * 1.15))
  }
  p
}

# --- effect sizes across cohorts ---------------------------------------------
panel_effect <- function(gene) {
  r <- tests |> dplyr::filter(gene == !!gene)
  if (!nrow(r)) return(NULL)
  r$label <- factor(r$label, levels = rev(r$label))

  ggplot(r, aes(x = cohens_d, y = label)) +
    geom_vline(xintercept = 0, linetype = "22", colour = PAL$grey_pt,
               linewidth = 0.4) +
    geom_errorbar(aes(xmin = d_lower, xmax = d_upper), orientation = "y",
                  width = 0.12, linewidth = 0.5, colour = PAL$ink) +
    geom_point(size = 2.4, shape = 21, fill = PAL$melanoma, colour = PAL$ink,
               stroke = 0.45) +
    # Direct value labels: identity and magnitude are legible with no colour.
    geom_text(aes(label = sprintf("%.2f", cohens_d)), vjust = -1.05,
              size = 2.6, colour = PAL$ink) +
    scale_x_continuous(n.breaks = 5) +
    labs(x = "Cohen's d (melanoma − nevus)", y = NULL,
         subtitle = "Effect size, 95% CI") +
    theme_pub(grid = TRUE)
}

# --- GSE4587 progression, one probe ------------------------------------------
STAGE_ORDER <- c("melanocytes", "benign nevi", "atypical nevi",
                 "melanoma in situ", "VGP melanoma", "MGP melanoma",
                 "lymph node metastasis", "metastatic melanoma culture")

panel_progression <- function(probe) {
  pd <- prog |>
    dplyr::filter(probe == !!probe, group %in% STAGE_ORDER) |>
    dplyr::mutate(group = factor(group, levels = STAGE_ORDER))
  if (!nrow(pd)) return(NULL)

  cols <- setNames(pal_sequential(nlevels(pd$group)), levels(pd$group))

  ggplot(pd, aes(x = group, y = value)) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.52, linewidth = 0.3,
                 colour = PAL$ink, fatten = 0) +
    geom_point(aes(fill = group), shape = 21, size = 2.1, stroke = 0.4,
               colour = PAL$ink,
               position = position_jitter(width = 0.07, height = 0)) +
    scale_fill_manual(values = cols, guide = "none") +
    labs(x = NULL, y = expression(bold(log[2]~"intensity")),
         subtitle = probe) +
    theme_pub(grid = TRUE) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

# ==============================================================================
# 3. Save helper
# ==============================================================================
save_geo <- function(plot, name, width, height) {
  for (fmt in c("pdf", "png")) {
    f <- file.path(FIG_GEO, paste0(name, ".", fmt))
    if (file.exists(f)) file.remove(f)
    args <- list(filename = f, plot = plot, width = width, height = height,
                 units = "mm")
    if (fmt == "pdf") args$device <- grDevices::cairo_pdf
    else { args$dpi <- 600; args$bg <- "white" }
    do.call(ggplot2::ggsave, args)
  }
  log_line("  wrote: ", name, ".{pdf,png}  (", width, " x ", height, " mm)",
           logfile = LOG)
}

# ==============================================================================
# 4. Nevus vs melanoma, one figure per gene
# ==============================================================================
hdr("4. Nevus vs melanoma")

for (g in GENES) {
  pA <- panel_lesion(g, "GSE112509")
  pB <- panel_lesion(g, "GSE98394")
  pC <- panel_effect(g)

  if (is.null(pA) && is.null(pB)) {
    log_line("  ", g, ": no data, skipped", logfile = LOG)
    next
  }

  parts <- Filter(Negate(is.null), list(pA, pB, pC))
  fig <- wrap_plots(parts, nrow = 1, widths = c(1, 1, 1.25)[seq_along(parts)]) +
    plot_annotation(
      title = bquote(bold(italic(.(g))~"in melanocytic nevi versus primary melanoma")),
      caption = "Independent patient cohorts. Wilcoxon rank-sum, two-sided. Points are individual samples.",
      tag_levels = "A") &
    theme(plot.tag = element_text(size = FIG$base_size + 3, face = "bold"),
          plot.title = element_text(size = FIG$base_size + 1, face = "bold"),
          plot.caption = element_text(size = FIG$base_size - 1,
                                      colour = PAL$grey_pt, hjust = 0))

  save_geo(fig, paste0("geo_", tolower(g), "_nevus_vs_melanoma"),
           width = FIG$width_double, height = 78)

  readr::write_csv(tests |> dplyr::filter(gene == g),
                   file.path(RES_GEO,
                             paste0("table_geo_", tolower(g),
                                    "_nevus_vs_melanoma.csv")))
}

# ==============================================================================
# 5. GSE4587 progression, one figure per gene
# ==============================================================================
hdr("5. GSE4587 progression series")

probe_map <- list(
  RASSF9  = c("RASSF9_210335_at", "RASSF9_241470_x_at"),
  RASSF10 = c("RASSF10_238755_at")
)

for (g in GENES) {
  ps <- Filter(Negate(is.null), lapply(probe_map[[g]], panel_progression))
  if (!length(ps)) {
    log_line("  ", g, ": no probes found in gse4587_values.csv", logfile = LOG)
    next
  }

  cap <- if (g == "RASSF9") {
    paste("1-3 samples per group: descriptive only, no test performed.",
          "Probes shown separately; 241470_x_at is a cross-hybridising",
          "probe and less reliable. Bars are group means.")
  } else {
    paste("1-3 samples per group: descriptive only, no test performed.",
          "Bars are group means.")
  }

  fig <- wrap_plots(ps, ncol = 1) +
    plot_annotation(
      title = bquote(bold(italic(.(g))~"across the melanocytic progression series (GSE4587)")),
      caption = cap,
      tag_levels = if (length(ps) > 1) "A" else NULL) &
    theme(plot.tag = element_text(size = FIG$base_size + 3, face = "bold"),
          plot.title = element_text(size = FIG$base_size + 1, face = "bold"),
          plot.caption = element_text(size = FIG$base_size - 1,
                                      colour = PAL$grey_pt, hjust = 0))

  save_geo(fig, paste0("geo_", tolower(g), "_progression_gse4587"),
           width = FIG$width_1_5,
           height = if (length(ps) > 1) 150 else 85)

  readr::write_csv(prog |> dplyr::filter(probe %in% probe_map[[g]]),
                   file.path(RES_GEO,
                             paste0("table_geo_", tolower(g),
                                    "_progression_gse4587.csv")))
}

# ==============================================================================
# 6. Inventory
# ==============================================================================
hdr("6. What is now in figures/geo")

for (f in sort(list.files(FIG_GEO))) log_line("  ", f, logfile = LOG)
log_line("\n  per-gene tables in results/geo:", logfile = LOG)
for (f in sort(list.files(RES_GEO, pattern = "^table_geo_"))) {
  log_line("  ", f, logfile = LOG)
}

log_line("\n  Every figure covers ONE gene, so either can become a main or a",
         logfile = LOG)
log_line("  supplementary figure without rework.", logfile = LOG)
log_line("\n  NOTE: the TCGA figures (scripts 04-10) still carry the previous",
         logfile = LOG)
log_line("  palette. Re-run those scripts to pick up the validated colours and",
         logfile = LOG)
log_line("  keep one house style across the paper.", logfile = LOG)

report_warnings(LOG)
write_session_info("15_geo_figures")
log_line("\nFull log: ", LOG, logfile = LOG)
