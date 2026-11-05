# ==============================================================================
# 18_depthmatched_gse230574.R
#
# Script 17 left one comparison unsafe. RASSF9 is sporadically detected in the
# primary nevus cells (2-3% of cells in one donor) and in NO cultured cell at
# all -- but the two runs differ in sequencing depth and in capture quality
# (ACTB 82-85% in the nevus run versus 100% in the culture run), so the
# contrast could be technical rather than biological.
#
# This script removes that confound three ways, in increasing strength:
#
#   (1) DEPTH MATCHING by binomial thinning. Every cell is computationally
#       downsampled to a common library size: a cell with k molecules of the
#       gene and L total is re-drawn as Binomial(k, T/L). This is the standard
#       way to ask "what would this cell look like if it had been sequenced as
#       shallowly as the worst group?" -- and it removes depth as an
#       explanation by construction, not by argument.
#
#   (2) UNIT OF ANALYSIS at the donor/line level, not the cell. Seven of the
#       twelve culture tags are 501Mel derivatives -- the same parental line
#       with different labels or treatments. Counting them as seven melanoma
#       replicates would be fake replication. Collapsed properly there are
#       4 nevus donors, 2 melanocyte cultures and 4 independent melanoma lines.
#
#   (3) AN EXPLICIT UPPER BOUND for groups with zero counts. "Not detected" is
#       not a number. The exact binomial upper bound on abundance, given the
#       total molecules sequenced, is -- and it is the honest way to state a
#       negative result.
#
# No p value is reported for cell-level comparisons. At the line level the
# groups are 4 versus 4 and a Wilcoxon is computed, but it is labelled
# underpowered and is NOT the basis of any claim.
#
# Input:  results/geo/table_geo_rassf9_singlecell_gse230574.csv  (and rassf10)
# Output: results/geo/table_geo_<gene>_depthmatched_gse230574.csv
#         results/geo/gse230574_upper_bounds.csv
#         figures/geo/geo_<gene>_depthmatched_gse230574.{pdf,png}
#         results/logs/18_depthmatched.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(patchwork))

LOG <- file.path(PATHS$logs, "18_depthmatched.log")
cat("GSE230574 depth-matched —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

RES_GEO <- file.path(PATHS$results, "geo")
FIG_GEO <- file.path(PATHS$figures, "geo")
GENES   <- c("RASSF9", "RASSF10")

set.seed(123)        # the thinning is stochastic; the seed makes it reproducible
N_THIN    <- 20      # thinning repeats, averaged -- one draw is noisy
MIN_FRAC  <- 0.01

# ==============================================================================
# 1. Load the per-cell tables written by script 17
# ==============================================================================
hdr("1. Input")

load_gene <- function(g) {
  f <- file.path(RES_GEO, paste0("table_geo_", tolower(g),
                                 "_singlecell_gse230574.csv"))
  if (!file.exists(f)) {
    stop("Missing ", basename(f), "\n  Run R/17_singlecell_gse230574.R first.")
  }
  tb <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE)
  names(tb)[names(tb) == g] <- "count"
  tb |> dplyr::select(run, tag, cell, lib_size, count) |> dplyr::mutate(gene = g)
}

dat <- purrr::map_dfr(GENES, load_gene)
log_line("  cells loaded: ", nrow(dat) / length(GENES), " per gene",
         logfile = LOG)

# ---- Biological grouping -----------------------------------------------------
# The 501Mel derivatives are collapsed to one biological line. mcherry_high /
# low / off are sorted fractions of the same culture and vemurafenib is a
# treatment arm; CDKN2A is an engineered variant of it. Treating them as
# independent melanoma samples would inflate n from 4 to 10 with no new biology.
classify <- function(tag) {
  dplyr::case_when(
    grepl("^SampleTag", tag)        ~ "Nevus (primary)",
    grepl("^melanocyte", tag)       ~ "Melanocyte (cultured)",
    TRUE                            ~ "Melanoma line")
}
line_of <- function(tag) {
  dplyr::case_when(
    grepl("^501Mel", tag)     ~ "501Mel",
    grepl("^melanocytes$", tag) ~ "NHM",
    grepl("^melanocyte_CDKN2A$", tag) ~ "NHM_CDKN2A",
    TRUE                      ~ tag)
}

dat <- dat |>
  dplyr::mutate(class = classify(tag), line = line_of(tag))

log_line("", logfile = LOG)
log_line("  biological grouping (cells):", logfile = LOG)
grp <- dat |> dplyr::filter(gene == GENES[1]) |>
  dplyr::count(class, line, name = "n_cells") |> dplyr::arrange(class, line)
for (i in seq_len(nrow(grp))) {
  log_line(sprintf("    %-24s %-14s %6d", grp$class[i], grp$line[i],
                   grp$n_cells[i]), logfile = LOG)
}
log_line("", logfile = LOG)
log_line("  independent units: ",
         dplyr::n_distinct(grp$line[grp$class == "Nevus (primary)"]),
         " nevus donors, ",
         dplyr::n_distinct(grp$line[grp$class == "Melanocyte (cultured)"]),
         " melanocyte cultures, ",
         dplyr::n_distinct(grp$line[grp$class == "Melanoma line"]),
         " melanoma lines", logfile = LOG)

# ==============================================================================
# 2. Choose the thinning target
# ==============================================================================
# The target is the depth the SHALLOWEST group being compared can actually
# support, so that no group is given an advantage. Cells below the target are
# dropped (they cannot be upsampled), and how many are lost is reported per
# group, because heavy loss in one group is itself a caveat.
hdr("2. Depth-matching target")

one <- dat |> dplyr::filter(gene == GENES[1])
med_by_line <- one |>
  dplyr::group_by(class, line) |>
  dplyr::summarise(median_lib = median(lib_size), n = dplyr::n(),
                   .groups = "drop")

TARGET <- floor(min(med_by_line$median_lib))
log_line("  median library size per unit:", logfile = LOG)
for (i in seq_len(nrow(med_by_line))) {
  log_line(sprintf("    %-24s %-14s %8d", med_by_line$class[i],
                   med_by_line$line[i], round(med_by_line$median_lib[i])),
           logfile = LOG)
}
log_line("", logfile = LOG)
log_line("  thinning target = ", TARGET,
         " molecules per cell (the shallowest unit's median)", logfile = LOG)

keep <- dat |> dplyr::filter(lib_size >= TARGET)
loss <- dat |>
  dplyr::filter(gene == GENES[1]) |>
  dplyr::group_by(class, line) |>
  dplyr::summarise(n = dplyr::n(), n_keep = sum(lib_size >= TARGET),
                   .groups = "drop") |>
  dplyr::mutate(pct_keep = 100 * n_keep / n)
log_line("", logfile = LOG)
log_line("  cells retained at the target depth:", logfile = LOG)
for (i in seq_len(nrow(loss))) {
  log_line(sprintf("    %-24s %-14s %6d of %6d (%5.1f%%)", loss$class[i],
                   loss$line[i], loss$n_keep[i], loss$n[i], loss$pct_keep[i]),
           logfile = LOG)
}

# ==============================================================================
# 3. Binomial thinning
# ==============================================================================
hdr("3. Binomial thinning")

log_line("  each retained cell re-drawn as Binomial(count, ", TARGET,
         "/lib_size),", logfile = LOG)
log_line("  averaged over ", N_THIN, " independent draws", logfile = LOG)

thin_once <- function(count, lib) {
  stats::rbinom(length(count), size = count, prob = pmin(TARGET / lib, 1))
}

thinned <- keep |>
  dplyr::mutate(
    det_rate = {
      hits <- replicate(N_THIN, thin_once(count, lib_size) > 0)
      if (is.matrix(hits)) rowMeans(hits) else as.numeric(hits)
    })

# ==============================================================================
# 4. Per-unit results, depth matched
# ==============================================================================
hdr("4. Depth-matched detection per unit")

res <- thinned |>
  dplyr::group_by(gene, class, line) |>
  dplyr::summarise(
    n_cells        = dplyr::n(),
    total_mol      = sum(lib_size),
    raw_pos        = mean(count > 0),
    thin_pos       = mean(det_rate),
    gene_molecules = sum(count),
    .groups = "drop") |>
  dplyr::arrange(gene, class, dplyr::desc(thin_pos))

for (g in GENES) {
  r <- res |> dplyr::filter(gene == g)
  log_line("", logfile = LOG)
  log_line("  ", g, logfile = LOG)
  log_line(sprintf("    %-24s %-14s %7s %9s %10s %10s",
                   "class", "unit", "cells", "%pos_raw", "%pos_thin",
                   "molecules"), logfile = LOG)
  for (i in seq_len(nrow(r))) {
    log_line(sprintf("    %-24s %-14s %7d %8.2f%% %9.2f%% %10d",
                     r$class[i], r$line[i], r$n_cells[i],
                     100 * r$raw_pos[i], 100 * r$thin_pos[i],
                     r$gene_molecules[i]), logfile = LOG)
  }

  by_class <- r |>
    dplyr::group_by(class) |>
    dplyr::summarise(units = dplyr::n(),
                     mean_thin = mean(thin_pos),
                     molecules = sum(gene_molecules),
                     total_mol = sum(total_mol), .groups = "drop")
  log_line("", logfile = LOG)
  log_line("    collapsed by class, AFTER depth matching:", logfile = LOG)
  for (i in seq_len(nrow(by_class))) {
    log_line(sprintf("      %-24s %d unit(s)  %6.2f%% of cells positive  (%d molecules in %d)",
                     by_class$class[i], by_class$units[i],
                     100 * by_class$mean_thin[i], by_class$molecules[i],
                     by_class$total_mol[i]), logfile = LOG)
  }
}

readr::write_csv(res, file.path(RES_GEO, "gse230574_depthmatched.csv"))
log_line("", logfile = LOG)
log_line("  wrote: results/geo/gse230574_depthmatched.csv", logfile = LOG)

# ==============================================================================
# 5. Upper bounds where nothing was detected
# ==============================================================================
# "Not detected" is a statement about the experiment, not the gene. The exact
# binomial upper bound turns it into a quantity: given this many molecules
# sequenced and this many hits, abundance is at most X. That is what belongs in
# a manuscript sentence about a negative result.
hdr("5. Upper bounds on abundance")

ub <- res |>
  dplyr::group_by(gene, class) |>
  dplyr::summarise(hits = sum(gene_molecules), total = sum(total_mol),
                   .groups = "drop") |>
  dplyr::rowwise() |>
  dplyr::mutate(
    frac_hat  = hits / total,
    upper_95  = stats::binom.test(hits, total)$conf.int[2],
    tpm_hat   = 1e6 * frac_hat,
    tpm_upper = 1e6 * upper_95) |>
  dplyr::ungroup()

log_line(sprintf("  %-10s %-24s %8s %12s %10s %12s",
                 "gene", "class", "hits", "molecules", "est_TPM", "upper_TPM"),
         logfile = LOG)
for (i in seq_len(nrow(ub))) {
  log_line(sprintf("  %-10s %-24s %8d %12.0f %10.3f %12.3f",
                   ub$gene[i], ub$class[i], ub$hits[i], ub$total[i],
                   ub$tpm_hat[i], ub$tpm_upper[i]), logfile = LOG)
}
readr::write_csv(ub, file.path(RES_GEO, "gse230574_upper_bounds.csv"))
log_line("", logfile = LOG)
log_line("  wrote: results/geo/gse230574_upper_bounds.csv", logfile = LOG)
log_line("  Read these as: 'in <class>, abundance is at most <upper_TPM> TPM",
         logfile = LOG)
log_line("  (95% binomial upper bound on <total> molecules sequenced)'.",
         logfile = LOG)

# ==============================================================================
# 6. Line-level test, reported as underpowered
# ==============================================================================
hdr("6. Line-level comparison")

for (g in GENES) {
  r <- res |> dplyr::filter(gene == g)
  a <- r$thin_pos[r$class == "Nevus (primary)"]
  b <- r$thin_pos[r$class %in% c("Melanocyte (cultured)", "Melanoma line")]
  log_line("", logfile = LOG)
  log_line("  ", g, ": nevus donors n = ", length(a),
           " vs cultured units n = ", length(b), logfile = LOG)
  if (length(a) >= 2 && length(b) >= 2) {
    wt <- suppressWarnings(stats::wilcox.test(a, b))
    log_line(sprintf("    nevus mean %.3f%% | cultured mean %.3f%% | W = %.1f, p = %s",
                     100 * mean(a), 100 * mean(b), wt$statistic,
                     fmt_p(wt$p.value)), logfile = LOG)
    log_line("    UNDERPOWERED. With ", length(a), " vs ", length(b),
             " units this p value is a description,", logfile = LOG)
    log_line("    not evidence. The upper bounds in section 5 are the result.",
             logfile = LOG)
  } else {
    log_line("    too few units for any test", logfile = LOG)
  }
}

# ==============================================================================
# 7. Figures
# ==============================================================================
hdr("7. Figures")

# PAL_CLASS comes from 00_setup.R and was selected by colour-vision simulation.
# Both panels run horizontally: unit names read at full size with no rotation,
# which is what rotated labels were costing in the earlier version.
make_fig <- function(g) {
  r <- res |> dplyr::filter(gene == g)
  if (!nrow(r)) return(NULL)
  r <- r |> dplyr::mutate(class = factor(class, levels = names(PAL_CLASS))) |>
    dplyr::arrange(class, thin_pos)
  r$line <- factor(r$line, levels = r$line)

  pa <- ggplot(r, aes(y = line, x = 100 * thin_pos, fill = class)) +
    geom_col(width = 0.68, colour = PAL$ink, linewidth = 0.25) +
    geom_vline(xintercept = 100 * MIN_FRAC, linetype = "22",
               colour = PAL$grey_pt, linewidth = 0.4) +
    facet_grid(class ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = PAL_CLASS, guide = "none", drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.20)), n.breaks = 5) +
    coord_cartesian(clip = "off") +
    labs(y = NULL, x = "Cells positive (%)") +
    theme_journal(grid = FALSE) +
    theme(panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey92"),
          panel.spacing.y = unit(5, "pt"),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 0,
                                           size = FIG$base_size - 1),
          axis.text.y = element_text(size = FIG$base_size - 1, hjust = 1),
          plot.margin = margin(4, 10, 4, 4, unit = "pt"))

  # The upper bounds span two to three orders of magnitude, so a linear axis
  # would flatten the cultured classes to invisible slivers and hide the very
  # comparison the panel exists to make. Log axis, with the value printed on
  # each bar so nothing depends on reading a log scale correctly. The label is
  # placed INSIDE the bar when the bar is long enough to hold it, which removes
  # any chance of text running off the panel.
  u <- ub |> dplyr::filter(gene == g) |>
    dplyr::mutate(class = factor(class, levels = rev(names(PAL_CLASS))),
                  lbl = sprintf("%.3f", tpm_upper),
                  inside = tpm_upper >= 0.5 * max(tpm_upper))
  pb <- ggplot(u, aes(y = class, x = tpm_upper, fill = class)) +
    geom_col(width = 0.6, colour = PAL$ink, linewidth = 0.25) +
    geom_text(data = dplyr::filter(u, inside),
              aes(label = lbl), hjust = 1.15, colour = "white",
              size = (FIG$base_size - 1.5) / .pt, fontface = "bold") +
    geom_text(data = dplyr::filter(u, !inside),
              aes(label = lbl), hjust = -0.18, colour = PAL$ink,
              size = (FIG$base_size - 1.5) / .pt, fontface = "bold") +
    scale_fill_manual(values = PAL_CLASS, guide = "none", drop = FALSE) +
    scale_x_log10(expand = expansion(mult = c(0.02, 0.32)),
                  labels = scales::label_number(accuracy = 0.001)) +
    coord_cartesian(clip = "off") +
    labs(y = NULL, x = "Upper bound on abundance (TPM, 95%)") +
    theme_journal(grid = FALSE) +
    theme(panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey92"),
          axis.text.y = element_text(size = FIG$base_size - 0.5, hjust = 1),
          axis.text.x = element_text(size = FIG$base_size - 1.5),
          plot.margin = margin(4, 10, 4, 4, unit = "pt"))

  (pa / pb) + patchwork::plot_layout(heights = c(2.2, 1))
}

for (g in GENES) {
  p <- make_fig(g)
  if (is.null(p)) { log_line("  ", g, ": skipped", logfile = LOG); next }

  fig <- p + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = FIG$base_size + 3, face = "bold",
                                  colour = PAL$ink))

  stem <- paste0("geo_", tolower(g), "_depthmatched_gse230574")
  save_fig_to(fig, stem, FIG_GEO, width = FIG$width_double, height = 145)

  uu <- ub |> dplyr::filter(gene == g)
  add_legend(
    stem,
    paste0("Figure. ", g, " in GSE230574 after depth matching."),
    paste0(
      "(A) Percentage of cells in which ", g, " was detected, after every ",
      "cell was computationally downsampled to ", TARGET, " molecules by ",
      "binomial thinning and averaged over ", N_THIN, " independent draws, so ",
      "that sequencing depth is matched across groups by construction. The ",
      "501Mel parental line, its CDKN2A and vemurafenib arms and its three ",
      "mCherry sort fractions are collapsed into a single unit, since they ",
      "are one biological line. Dashed line, the ", 100 * MIN_FRAC,
      "% detection floor. (B) Exact binomial 95% upper bound on ", g,
      " abundance for each class, given the total molecules sequenced (",
      paste(sprintf("%s, %.3f TPM", uu$class, uu$tpm_upper), collapse = "; "),
      "); the axis is logarithmic and each value is printed on its bar. ",
      "Units of analysis: 4 nevus donors, 2 melanocyte cultures and 4 ",
      "independent melanoma lines. Cells within a culture are not independent ",
      "biological replicates, so the bounds rather than a group test are ",
      "reported."))

  log_line("  wrote: ", stem, ".{pdf,png,svg}", logfile = LOG)

  readr::write_csv(res |> dplyr::filter(gene == g),
                   file.path(RES_GEO, paste0("table_geo_", tolower(g),
                                             "_depthmatched_gse230574.csv")))
}

report_warnings(logfile = LOG)
write_session_info("18_depthmatched")
log_line("\nFull log: ", LOG, logfile = NULL)
