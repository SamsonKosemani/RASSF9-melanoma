# ==============================================================================
# 10_mediation_exploratory.R
#
#                        *** EXPLORATORY ANALYSIS ***
#
# This analysis was NOT pre-specified. It was prompted by a pattern noticed in
# the script 08 results: RASSF9 promoter methylation correlates with UV burden
# (rho = +0.173, p = 0.002) and strongly with expression (rho = -0.351,
# p = 1.1e-10), while the direct UV-expression association is weak. That is the
# shape an indirect effect would take.
#
# Everything here must be labelled exploratory in the manuscript. It generates
# a hypothesis; it does not test one.
#
# WHAT THIS CANNOT ESTABLISH
#   1. Direction. These data are cross-sectional and measured at one moment in
#      one tumour. "UV damage -> methylation -> silencing" fits the data, but so
#      does "methylation state -> altered repair/expression", and so does a
#      common cause driving both. Statistical mediation does not order events
#      in time.
#   2. No unmeasured confounding. Mediation estimates are unbiased only if
#      nothing unmeasured affects both mediator and outcome. In bulk tumour
#      data that assumption is heroic: cell-type composition, differentiation
#      state and global methylation drift all plausibly affect both.
#   3. Causality of any kind. Nothing here is an intervention.
#
# The honest claim this supports: "the association between UV burden and RASSF9
# expression is statistically consistent with an indirect path through promoter
# methylation" -- and nothing stronger.
#
# RASSF10 is included as a specificity check, not as a parallel finding: script
# 08 showed its promoter methylation does not track UV burden at all, so the
# indirect effect should be null there. If it is, that supports the RASSF9
# result being specific rather than an artefact of the method.
#
# Model (all paths adjusted for CPE purity and sample type):
#   a  : UV        -> methylation
#   b  : methylation -> expression, holding UV
#   c  : UV        -> expression            (total)
#   c' : UV        -> expression, holding methylation   (direct)
#   ab : indirect effect, with percentile bootstrap CI
#
# Output: figures/figureE1_rassf9_mediation.*  (E = exploratory)
#         results/table_exploratory_10_mediation.csv
#         results/logs/10_mediation.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(patchwork)
  library(scales)
})

LOG <- file.path(PATHS$logs, "10_mediation.log")
cat("Exploratory mediation —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

GENES   <- c("RASSF9", "RASSF10")
PURITY  <- "CPE"
N_BOOT  <- 5000

log_line("*** EXPLORATORY -- not pre-specified. See header for limits. ***",
         logfile = LOG)

# ==============================================================================
# 1. Data
# ==============================================================================
hdr("1. Data")

coh    <- readr::read_csv(file.path(PATHS$processed, "analysis_cohort.csv"))
betas  <- readr::read_csv(file.path(PATHS$processed, "methylation_betas.csv"))
probes <- readr::read_csv(file.path(PATHS$processed, "methylation_probes.csv"))

d <- dplyr::inner_join(coh, betas, by = "patient") |>
  dplyr::mutate(
    sample_type = dplyr::case_when(
      grepl("^Primary", sample_type, ignore.case = TRUE)  ~ "Primary",
      grepl("Metastatic", sample_type, ignore.case = TRUE) ~ "Metastatic",
      TRUE ~ NA_character_),
    sample_type = factor(sample_type, levels = c("Primary", "Metastatic"))
  )

for (g in GENES) {
  pr <- probes |> dplyr::filter(gene == g, promoter, on_array)
  pr <- pr[pr$probe %in% names(d), ]
  d[[paste0(g, "_prom_beta")]] <- rowMeans(
    as.matrix(d[, pr$probe, drop = FALSE]), na.rm = TRUE)
}

log_line("Cohort: n = ", nrow(d), logfile = LOG)

# ==============================================================================
# 2. Path estimation
# ==============================================================================
hdr("2. Path models")

# Standardising X, M and Y makes a, b and the indirect effect comparable
# across genes and keeps the bootstrap on a sensible scale.
z <- function(x) as.numeric(scale(x))

paths <- function(dat) {
  fit_m  <- lm(M ~ X + P + S, data = dat)                 # a path
  fit_y  <- lm(Y ~ X + M + P + S, data = dat)             # b and c' paths
  fit_t  <- lm(Y ~ X + P + S, data = dat)                 # c path (total)
  c(a  = unname(coef(fit_m)["X"]),
    b  = unname(coef(fit_y)["M"]),
    cp = unname(coef(fit_y)["X"]),
    c  = unname(coef(fit_t)["X"]),
    ab = unname(coef(fit_m)["X"] * coef(fit_y)["M"]))
}

rows <- list()
boot_store <- list()

for (g in GENES) {
  dd <- d[, c(paste0(g, "_log2tpm"), "log10_SBS7", paste0(g, "_prom_beta"),
              PURITY, "sample_type")]
  names(dd) <- c("Y", "X", "M", "P", "S")
  dd <- dd[stats::complete.cases(dd), ]
  dd$Y <- z(dd$Y); dd$X <- z(dd$X); dd$M <- z(dd$M); dd$P <- z(dd$P)

  est <- paths(dd)

  # --- p values for the individual paths ---
  fm <- summary(lm(M ~ X + P + S, data = dd))$coefficients
  fy <- summary(lm(Y ~ X + M + P + S, data = dd))$coefficients
  ft <- summary(lm(Y ~ X + P + S, data = dd))$coefficients

  # --- percentile bootstrap for the indirect effect ---
  bs <- matrix(NA_real_, nrow = N_BOOT, ncol = 5,
               dimnames = list(NULL, c("a", "b", "cp", "c", "ab")))
  for (i in seq_len(N_BOOT)) {
    idx <- sample.int(nrow(dd), replace = TRUE)
    bs[i, ] <- tryCatch(paths(dd[idx, ]), error = function(e) rep(NA_real_, 5))
  }
  ab_ci <- quantile(bs[, "ab"], c(0.025, 0.975), na.rm = TRUE)

  # Proportion mediated is only interpretable when total and indirect effects
  # share a sign; otherwise it is a meaningless ratio.
  prop_med <- if (sign(est["ab"]) == sign(est["c"]) && abs(est["c"]) > 1e-8) {
    est["ab"] / est["c"]
  } else NA_real_

  boot_store[[g]] <- bs

  rows[[g]] <- tibble::tibble(
    gene = g, n = nrow(dd),
    a = est["a"],  a_p  = fm["X", "Pr(>|t|)"],
    b = est["b"],  b_p  = fy["M", "Pr(>|t|)"],
    c_total = est["c"], c_p = ft["X", "Pr(>|t|)"],
    c_direct = est["cp"], cp_p = fy["X", "Pr(>|t|)"],
    ab = est["ab"], ab_ci_low = ab_ci[1], ab_ci_high = ab_ci[2],
    ab_significant = !(ab_ci[1] <= 0 && ab_ci[2] >= 0),
    prop_mediated = prop_med
  )

  r <- rows[[g]]
  log_line("\n  ", g, " (n = ", r$n, ", standardised units):", logfile = LOG)
  log_line(sprintf("    a  UV -> methylation        %+.4f  (p = %s)",
                   r$a, fmt_p(r$a_p)), logfile = LOG)
  log_line(sprintf("    b  methylation -> expression %+.4f  (p = %s)",
                   r$b, fmt_p(r$b_p)), logfile = LOG)
  log_line(sprintf("    c  UV -> expression, total   %+.4f  (p = %s)",
                   r$c_total, fmt_p(r$c_p)), logfile = LOG)
  log_line(sprintf("    c' UV -> expression, direct  %+.4f  (p = %s)",
                   r$c_direct, fmt_p(r$cp_p)), logfile = LOG)
  log_line(sprintf("    ab indirect effect           %+.4f  [%+.4f, %+.4f]  %s",
                   r$ab, r$ab_ci_low, r$ab_ci_high,
                   ifelse(r$ab_significant, "CI excludes zero",
                          "CI includes zero")), logfile = LOG)
  if (!is.na(r$prop_mediated)) {
    log_line(sprintf("    proportion of total effect via methylation: %.1f%%",
                     100 * r$prop_mediated), logfile = LOG)
  } else {
    log_line("    proportion mediated: not interpretable (indirect and total ",
             "effects differ in sign, or the total effect is ~0)", logfile = LOG)
  }
}
med_tbl <- dplyr::bind_rows(rows)

# ==============================================================================
# 3. Reading the result
# ==============================================================================
hdr("3. How to read this")

r9  <- med_tbl |> dplyr::filter(gene == "RASSF9")
r10 <- med_tbl |> dplyr::filter(gene == "RASSF10")

log_line("  Specificity check: the indirect path should be present for RASSF9",
         logfile = LOG)
log_line("  and absent for RASSF10, whose promoter methylation does not track",
         logfile = LOG)
log_line("  UV burden (script 08: rho = -0.001).", logfile = LOG)
log_line(sprintf("    RASSF9  indirect effect CI excludes zero: %s",
                 r9$ab_significant), logfile = LOG)
log_line(sprintf("    RASSF10 indirect effect CI excludes zero: %s",
                 r10$ab_significant), logfile = LOG)

if (isTRUE(r9$ab_significant) && !isTRUE(r10$ab_significant)) {
  log_line("\n  Pattern as predicted. This STRENGTHENS the plausibility of the",
           logfile = LOG)
  log_line("  indirect path for RASSF9 but still does not establish direction",
           logfile = LOG)
  log_line("  or causality -- see the header. Report as exploratory.",
           logfile = LOG)
} else if (isTRUE(r9$ab_significant) && isTRUE(r10$ab_significant)) {
  log_line("\n  Both genes show an indirect effect. The specificity check did",
           logfile = LOG)
  log_line("  NOT hold, which weakens the RASSF9 interpretation: a method that",
           logfile = LOG)
  log_line("  finds mediation everywhere is finding structure, not biology.",
           logfile = LOG)
} else {
  log_line("\n  No indirect effect for RASSF9. The pattern noticed in script 08",
           logfile = LOG)
  log_line("  does not survive formal estimation. Report nothing, or report it",
           logfile = LOG)
  log_line("  as a negative exploratory result.", logfile = LOG)
}

log_line("\n  Whatever the numbers, the wording in the manuscript should be:",
         logfile = LOG)
log_line("  'statistically consistent with an indirect path' -- never",
         logfile = LOG)
log_line("  'mediated by', 'drives', 'leads to' or 'causes'.", logfile = LOG)

# ==============================================================================
# 4. Figure
# ==============================================================================
hdr("4. Figure")

purge_fig_prefix("figureE1")

panel_paths <- function(g) {
  r <- med_tbl |> dplyr::filter(gene == g)
  bs <- boot_store[[g]]
  pd <- tibble::tibble(
    path = factor(c("a: UV->meth", "b: meth->expr", "c: total",
                    "c': direct", "ab: indirect"),
                  levels = rev(c("a: UV->meth", "b: meth->expr", "c: total",
                                 "c': direct", "ab: indirect"))),
    est = c(r$a, r$b, r$c_total, r$c_direct, r$ab),
    low = c(quantile(bs[, "a"],  0.025, na.rm = TRUE),
            quantile(bs[, "b"],  0.025, na.rm = TRUE),
            quantile(bs[, "c"],  0.025, na.rm = TRUE),
            quantile(bs[, "cp"], 0.025, na.rm = TRUE), r$ab_ci_low),
    high = c(quantile(bs[, "a"],  0.975, na.rm = TRUE),
             quantile(bs[, "b"],  0.975, na.rm = TRUE),
             quantile(bs[, "c"],  0.975, na.rm = TRUE),
             quantile(bs[, "cp"], 0.975, na.rm = TRUE), r$ab_ci_high)
  ) |>
    dplyr::mutate(excl = low > 0 | high < 0)

  ggplot(pd, aes(x = est, y = path)) +
    geom_vline(xintercept = 0, linetype = "22", colour = PAL$grey_pt,
               linewidth = 0.4) +
    geom_errorbar(aes(xmin = low, xmax = high), orientation = "y",
                  width = 0.14, linewidth = 0.5, colour = PAL$ink) +
    geom_point(aes(fill = excl), size = 2.3, shape = 21, colour = PAL$ink,
               stroke = 0.4) +
    scale_fill_manual(values = c("TRUE" = PAL$uv_high, "FALSE" = "white"),
                      guide = "none") +
    labs(x = "Standardised effect", y = NULL,
         subtitle = "Path estimates, bootstrap CI") +
    theme_pub(grid = TRUE)
}

panel_boot <- function(g) {
  r <- med_tbl |> dplyr::filter(gene == g)
  bs <- as.data.frame(boot_store[[g]])
  ggplot(bs, aes(x = ab)) +
    geom_histogram(bins = 50, fill = PAL$hist, colour = "white",
                   linewidth = 0.1) +
    geom_vline(xintercept = 0, linetype = "22", colour = PAL$ink,
               linewidth = 0.45) +
    geom_vline(xintercept = c(r$ab_ci_low, r$ab_ci_high), linetype = "solid",
               colour = PAL$uv_high, linewidth = 0.45) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = "Indirect effect (ab)", y = "Bootstrap samples",
         subtitle = paste0(format(N_BOOT, big.mark = ","), " resamples")) +
    theme_pub()
}

panel_apath <- function(g) {
  bcol <- paste0(g, "_prom_beta")
  pd <- d[stats::complete.cases(d[[bcol]], d$log10_SBS7), ]
  ct <- suppressWarnings(cor.test(pd[[bcol]], pd$log10_SBS7,
                                  method = "spearman"))
  ggplot(pd, aes(x = log10_SBS7, y = .data[[bcol]])) +
    geom_point(size = 0.9, alpha = 0.6, colour = PAL$accent, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = PAL$ink,
                fill = "grey70", linewidth = 0.6, alpha = 0.25) +
    annotate("text", x = max(pd$log10_SBS7) * 0.99,
             y = max(pd[[bcol]]) * 0.99,
             label = sprintf("rho == %+.3f", ct$estimate), parse = TRUE,
             hjust = 1, vjust = 1, size = 2.9, colour = PAL$ink) +
    labs(x = expression(bold(log[10]*"(SBS7 + 1)")),
         y = "Promoter methylation (beta)",
         subtitle = "a path: UV vs methylation") +
    theme_pub()
}

for (g in GENES) {
  tag <- if (g == "RASSF9") "figureE1_rassf9_mediation" else
    "figureE2_rassf10_mediation"
  purge_fig_prefix(substr(tag, 1, 9))
  fig <- (panel_paths(g) | panel_boot(g) | panel_apath(g)) +
    plot_annotation(
      tag_levels = "A",
      caption = "Exploratory analysis; not pre-specified. See Methods.") &
    theme(plot.tag = element_text(size = FIG$base_size + 3, face = "bold"),
          plot.caption = element_text(size = FIG$base_size - 1,
                                      colour = PAL$grey_pt, hjust = 0))
  save_fig(fig, tag, width = FIG$width_double, height = 75)
}

# ==============================================================================
# 5. Table
# ==============================================================================
hdr("5. Results table")

readr::write_csv(med_tbl,
                 file.path(PATHS$results, "table_exploratory_10_mediation.csv"))
cat(capture.output(print(as.data.frame(med_tbl), row.names = FALSE, digits = 3)),
    sep = "\n", file = LOG, append = TRUE)
print(as.data.frame(med_tbl), digits = 3)

hdr("6. Warnings")
report_warnings(LOG)

write_session_info("10_mediation")

hdr("DONE")
log_line("Files are named 'exploratory' and 'figureE' deliberately, so that ",
         "this analysis cannot be mistaken for a pre-specified result when ",
         "the manuscript is assembled.", logfile = LOG)
log_line("Full log: ", LOG, logfile = LOG)
