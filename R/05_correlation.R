# ==============================================================================
# 05_correlation.R
#
# Question 3: does expression scale quantitatively with UV mutational burden?
#
# Pre-specified (docs/analysis-plan.md section 5):
#   - Spearman is PRIMARY (no distributional assumption, robust to the heavy
#     right skew of SBS7 and the floor effect in expression)
#   - Pearson on log10(SBS7+1) is SECONDARY and reported alongside
#   - Pearson on RAW SBS7 is also computed, because the earlier analysis found
#     a significant Spearman alongside a null Pearson. That discrepancy is a
#     scale artefact, not a contradiction, and this script demonstrates it
#     rather than leaving it unexplained.
#   - Secondary proxies: SBS7 fraction, total mutational burden
#
# OUTPUT: RASSF9  -> figures/figure2_rassf9_correlation.*   + table_main_*
#         RASSF10 -> figures/figureS2_rassf10_correlation.* + table_supp_*
#
# Input:  data/processed/analysis_cohort.csv
# Output: results/table_main_05_rassf9_correlation.csv
#         results/table_supp_05_rassf10_correlation.csv
#         results/logs/05_correlation.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(patchwork)
  library(scales)
})

LOG <- file.path(PATHS$logs, "05_correlation.log")
cat("Correlation analysis —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

MAIN_GENE <- "RASSF9"
SUPP_GENE <- "RASSF10"

# ==============================================================================
# 1. Load and re-derive groups (for point colouring only)
# ==============================================================================
hdr("1. Cohort")

coh <- readr::read_csv(file.path(PATHS$processed, "analysis_cohort.csv"))
tert <- quantile(coh$SBS7_total, c(1/3, 2/3), na.rm = TRUE)

coh <- coh |>
  dplyr::mutate(
    uv_group = dplyr::case_when(
      SBS7_total <= tert[1] ~ "UV-low",
      SBS7_total >= tert[2] ~ "UV-high",
      TRUE                  ~ "Middle"
    ),
    uv_group = factor(uv_group, levels = c("UV-low", "Middle", "UV-high"))
  )

log_line("Cohort: n = ", nrow(coh), logfile = LOG)

# ==============================================================================
# 2. Correlations
# ==============================================================================
hdr("2. Correlation tests")

cor_row <- function(d, xvar, gene, xlabel, method) {
  ycol <- paste0(gene, "_log2tpm")
  ok <- stats::complete.cases(d[[xvar]], d[[ycol]])
  x <- d[[xvar]][ok]; y <- d[[ycol]][ok]
  ct <- suppressWarnings(cor.test(x, y, method = method))

  tibble::tibble(
    gene      = gene,
    predictor = xlabel,
    variable  = xvar,
    method    = ifelse(method == "spearman", "Spearman", "Pearson"),
    n         = length(x),
    estimate  = unname(ct$estimate),
    conf_low  = if (!is.null(ct$conf.int)) ct$conf.int[1] else NA_real_,
    conf_high = if (!is.null(ct$conf.int)) ct$conf.int[2] else NA_real_,
    p_value   = ct$p.value
  )
}

# Spearman is invariant to monotone transformation, so SBS7_total and
# log10_SBS7 give identical rho. Pearson is not, which is the whole point.
spec <- list(
  list("SBS7_total",    "UV burden (SBS7)",        "spearman"),
  list("SBS7_total",    "UV burden (SBS7), raw",   "pearson"),
  list("log10_SBS7",    "UV burden (SBS7), log10", "pearson"),
  list("SBS7_fraction", "SBS7 fraction",           "spearman"),
  list("SBS7_fraction", "SBS7 fraction",           "pearson"),
  list("TMB_sig",       "Total mutational burden", "spearman"),
  list("log10_TMB",     "Total burden, log10",     "pearson")
)

rows <- list()
for (gene in c(MAIN_GENE, SUPP_GENE)) {
  log_line("\n", gene, ":", logfile = LOG)
  for (s in spec) {
    r <- cor_row(coh, s[[1]], gene, s[[2]], s[[3]])
    rows[[length(rows) + 1]] <- r
    log_line(sprintf("  %-26s %-9s n=%d  estimate=%+.3f  p=%s",
                     r$predictor, r$method, r$n, r$estimate, fmt_p(r$p_value)),
             logfile = LOG)
  }
}
cor_tbl <- dplyr::bind_rows(rows)

# --- Explain the Spearman/Pearson divergence explicitly -----------------------
hdr("3. Why Spearman and Pearson differ")

for (gene in c(MAIN_GENE, SUPP_GENE)) {
  sp  <- cor_tbl |> dplyr::filter(gene == !!gene, variable == "SBS7_total",
                                  method == "Spearman")
  pr  <- cor_tbl |> dplyr::filter(gene == !!gene, variable == "SBS7_total",
                                  method == "Pearson")
  plg <- cor_tbl |> dplyr::filter(gene == !!gene, variable == "log10_SBS7",
                                  method == "Pearson")
  log_line(sprintf(
    "  %-8s Spearman rho = %+.3f (p=%s) | Pearson on raw = %+.3f (p=%s) | Pearson on log10 = %+.3f (p=%s)",
    gene, sp$estimate, fmt_p(sp$p_value), pr$estimate, fmt_p(pr$p_value),
    plg$estimate, fmt_p(plg$p_value)), logfile = LOG)
}
log_line("\n  SBS7 spans three orders of magnitude (5 to 23,393). Pearson on the",
         logfile = LOG)
log_line("  raw scale is dominated by a handful of hypermutated tumours and is",
         logfile = LOG)
log_line("  therefore uninformative here; on the log scale it should align with",
         logfile = LOG)
log_line("  Spearman. Report Spearman as primary and the log-scale Pearson as",
         logfile = LOG)
log_line("  support; do not present the raw-scale Pearson as a contradiction.",
         logfile = LOG)

# ==============================================================================
# 4. Figures
# ==============================================================================
hdr("4. Figures")

purge_fig_prefix("figure2")
purge_fig_prefix("figureS2")

scatter_panel <- function(gene, xvar, xlab_expr, log_x = FALSE, subtitle) {
  ycol <- paste0(gene, "_log2tpm")
  r <- cor_tbl |> dplyr::filter(gene == !!gene, variable == !!xvar,
                                method == "Spearman")
  d <- coh[stats::complete.cases(coh[[xvar]], coh[[ycol]]), ]

  lab <- sprintf("rho == %+.3f", r$estimate)
  plab <- if (r$p_value < 0.0001) "italic(P) < 0.0001" else
    sprintf("italic(P) == %s", formatC(r$p_value, format = "f",
                                       digits = ifelse(r$p_value < 0.01, 4, 3)))

  p <- ggplot(d, aes(x = .data[[xvar]], y = .data[[ycol]])) +
    geom_point(aes(colour = uv_group), size = 0.9, alpha = 0.65, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = PAL$ink, fill = "grey70", linewidth = 0.6,
                alpha = 0.25) +
    scale_colour_manual(values = PAL_UV, name = NULL) +
    labs(x = xlab_expr,
         y = bquote(bold(italic(.(gene))~log[2]*"(TPM+1)")),
         subtitle = subtitle) +
    theme_pub() +
    theme(legend.position = "none")

  if (log_x) {
    p <- p + scale_x_log10(labels = label_number(big.mark = ","))
  }

  # Annotate in the top-right corner, in data coordinates.
  xr <- range(d[[xvar]], na.rm = TRUE)
  yr <- range(d[[ycol]], na.rm = TRUE)
  xpos <- if (log_x) 10^(log10(xr[1]) + 0.97 * diff(log10(xr))) else
    xr[1] + 0.97 * diff(xr)

  p +
    annotate("text", x = xpos, y = yr[2] * 1.02, label = lab, parse = TRUE,
             hjust = 1, vjust = 1, size = 2.9, colour = PAL$ink) +
    annotate("text", x = xpos, y = yr[2] * 0.90, label = plab, parse = TRUE,
             hjust = 1, vjust = 1, size = 2.9, colour = PAL$ink)
}

build_fig <- function(gene) {
  pA <- scatter_panel(gene, "SBS7_total",
                      expression(bold("SBS7 mutations")), log_x = TRUE,
                      subtitle = "UV burden (primary proxy)")
  pB <- scatter_panel(gene, "SBS7_fraction",
                      expression(bold("SBS7 fraction of total burden")),
                      log_x = FALSE,
                      subtitle = "UV burden, relative")
  pC <- scatter_panel(gene, "TMB_sig",
                      expression(bold("Total mutations")), log_x = TRUE,
                      subtitle = "Total mutational burden")

  (pA | pB | pC) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = FIG$base_size + 3, face = "bold"))
}

save_fig(build_fig(MAIN_GENE), "figure2_rassf9_correlation",
         width = FIG$width_double, height = 68)

save_fig(build_fig(SUPP_GENE), "figureS2_rassf10_correlation",
         width = FIG$width_double, height = 68)

# ==============================================================================
# 5. Results tables
# ==============================================================================
hdr("5. Results tables")

main_tbl <- cor_tbl |> dplyr::filter(gene == MAIN_GENE)
supp_tbl <- cor_tbl |> dplyr::filter(gene == SUPP_GENE)

readr::write_csv(main_tbl,
                 file.path(PATHS$results, "table_main_05_rassf9_correlation.csv"))
readr::write_csv(supp_tbl,
                 file.path(PATHS$results, "table_supp_05_rassf10_correlation.csv"))

log_line("\nMAIN (", MAIN_GENE, "):", logfile = LOG)
cat(capture.output(print(as.data.frame(main_tbl), row.names = FALSE, digits = 3)),
    sep = "\n", file = LOG, append = TRUE)
print(as.data.frame(main_tbl), digits = 3)

log_line("\nSUPPLEMENTARY (", SUPP_GENE, "):", logfile = LOG)
cat(capture.output(print(as.data.frame(supp_tbl), row.names = FALSE, digits = 3)),
    sep = "\n", file = LOG, append = TRUE)
print(as.data.frame(supp_tbl), digits = 3)

hdr("6. Warnings")
report_warnings(LOG)

write_session_info("05_correlation")

hdr("DONE")
log_line("Figure 2 (main) : figures/figure2_rassf9_correlation.*", logfile = LOG)
log_line("Figure S2 (supp): figures/figureS2_rassf10_correlation.*", logfile = LOG)
log_line("\nNote: a correlation of this magnitude explains a few percent of ",
         "variance at most. The adjusted model in script 06 is what decides ",
         "how this is reported.", logfile = LOG)
log_line("Full log: ", LOG, logfile = LOG)
