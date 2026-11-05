# ==============================================================================
# 06_adjusted_model.R
#
# Question 4 (DECISIVE): does the association between UV mutational burden and
# expression survive adjustment for tumour purity and sample type?
#
# Pre-specified model (docs/analysis-plan.md section 5):
#     log2(expr TPM+1) ~ log10(SBS7+1) + purity + sample_type
#
# Purity: CPE (consensus estimate, complete coverage) is PRIMARY per the
# amended D10; ESTIMATE, LUMP and IHC are pre-specified sensitivity analyses.
# ABSOLUTE is unavailable for SKCM.
#
# Additional sensitivity, declared here and justified by script 05: RASSF9's
# univariate association is rank-based (Spearman significant, Pearson on the
# log scale null). A linear model is therefore not the natural test for it, so
# the same model is refit on rank-inverse-normal transformed expression. This
# is reported as a sensitivity analysis, NOT substituted for the primary model.
#
# Per the analysis plan: if the association does not survive adjustment, that
# is reported as a negative result in the main text.
#
# Input:  data/processed/analysis_cohort.csv
# Output: figures/figure3_rassf9_adjusted.*, figureS3_rassf10_adjusted.*
#         results/table_main_06_rassf9_adjusted.csv
#         results/table_supp_06_rassf10_adjusted.csv
#         results/logs/06_adjusted.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(patchwork)
  library(scales)
})

LOG <- file.path(PATHS$logs, "06_adjusted.log")
cat("Adjusted models —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

MAIN_GENE <- "RASSF9"
SUPP_GENE <- "RASSF10"
PURITY_PRIMARY <- "CPE"
PURITY_SENS    <- c("ESTIMATE", "LUMP", "IHC")

# ==============================================================================
# 1. Prepare data
# ==============================================================================
hdr("1. Model data")

coh <- readr::read_csv(file.path(PATHS$processed, "analysis_cohort.csv"))

coh <- coh |>
  dplyr::mutate(
    sample_type = dplyr::case_when(
      grepl("^Primary", sample_type, ignore.case = TRUE) ~ "Primary",
      grepl("Metastatic", sample_type, ignore.case = TRUE) ~ "Metastatic",
      TRUE ~ NA_character_
    ),
    sample_type = factor(sample_type, levels = c("Primary", "Metastatic"))
  )

log_line("Cohort: n = ", nrow(coh), logfile = LOG)
log_line("Sample type (reference = Primary):", logfile = LOG)
st <- as.data.frame(table(coh$sample_type, useNA = "ifany"))
names(st) <- c("sample_type", "n")
cat(capture.output(print(st, row.names = FALSE)), sep = "\n",
    file = LOG, append = TRUE)
print(st)

for (p in c(PURITY_PRIMARY, PURITY_SENS)) {
  if (p %in% names(coh)) {
    log_line(sprintf("  %-9s available %d/%d, median %.3f", p,
                     sum(!is.na(coh[[p]])), nrow(coh),
                     median(coh[[p]], na.rm = TRUE)), logfile = LOG)
  } else {
    log_line("  [FLAG] purity column '", p, "' absent from cohort file",
             logfile = LOG)
  }
}

# ==============================================================================
# 2. Is purity itself related to UV burden?  (the reason for adjusting)
# ==============================================================================
hdr("2. Confounding structure")

for (p in c(PURITY_PRIMARY, PURITY_SENS)) {
  if (!p %in% names(coh)) next
  ct <- suppressWarnings(cor.test(coh$log10_SBS7, coh[[p]], method = "spearman"))
  log_line(sprintf("  log10(SBS7) vs %-9s : rho = %+.3f, p = %s", p,
                   ct$estimate, fmt_p(ct$p.value)), logfile = LOG)
}

pt <- suppressWarnings(wilcox.test(log10_SBS7 ~ sample_type, data = coh,
                                   exact = FALSE))
log_line(sprintf("  log10(SBS7) by sample type: p = %s", fmt_p(pt$p.value)),
         logfile = LOG)

for (g in c(MAIN_GENE, SUPP_GENE)) {
  col <- paste0(g, "_log2tpm")
  ct <- suppressWarnings(cor.test(coh[[col]], coh[[PURITY_PRIMARY]],
                                  method = "spearman"))
  log_line(sprintf("  %-8s vs %s purity: rho = %+.3f, p = %s", g,
                   PURITY_PRIMARY, ct$estimate, fmt_p(ct$p.value)),
           logfile = LOG)
}

log_line("\n  If purity correlates with BOTH the exposure and the outcome, it",
         logfile = LOG)
log_line("  is a confounder and the adjusted estimate is the one to believe.",
         logfile = LOG)

# ==============================================================================
# 3. Fit models
# ==============================================================================
hdr("3. Model fitting")

rank_inv_norm <- function(x) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  qnorm((r - 0.5) / sum(!is.na(x)))
}

# Variance inflation without extra dependencies.
vif_manual <- function(fit) {
  X <- model.matrix(fit)[, -1, drop = FALSE]
  if (ncol(X) < 2) return(setNames(rep(NA_real_, ncol(X)), colnames(X)))
  sapply(seq_len(ncol(X)), function(j) {
    r2 <- summary(lm(X[, j] ~ X[, -j, drop = FALSE]))$r.squared
    1 / (1 - r2)
  }) |> setNames(colnames(X))
}

fit_model <- function(d, gene, purity, label, transform = FALSE) {
  ycol <- paste0(gene, "_log2tpm")
  dd <- d[, c(ycol, "log10_SBS7", purity, "sample_type")]
  names(dd)[1:3] <- c("y", "uv", "purity")
  dd <- dd[stats::complete.cases(dd), ]
  if (nrow(dd) < 30) return(NULL)
  if (transform) dd$y <- rank_inv_norm(dd$y)

  fit <- lm(y ~ uv + purity + sample_type, data = dd)
  s <- summary(fit)
  co <- s$coefficients

  # Standardised slope for the UV term, for cross-model comparability.
  beta_std <- co["uv", "Estimate"] * sd(dd$uv) / sd(dd$y)
  ci <- confint(fit)["uv", ]

  list(
    fit = fit,
    row = tibble::tibble(
      gene        = gene,
      model       = label,
      purity_used = purity,
      n           = nrow(dd),
      uv_beta     = co["uv", "Estimate"],
      uv_se       = co["uv", "Std. Error"],
      uv_ci_low   = ci[1],
      uv_ci_high  = ci[2],
      uv_beta_std = beta_std,
      uv_p        = co["uv", "Pr(>|t|)"],
      purity_beta = co["purity", "Estimate"],
      purity_p    = co["purity", "Pr(>|t|)"],
      type_beta   = co[grep("^sample_type", rownames(co))[1], "Estimate"],
      type_p      = co[grep("^sample_type", rownames(co))[1], "Pr(>|t|)"],
      adj_r2      = s$adj.r.squared,
      max_vif     = max(vif_manual(fit), na.rm = TRUE)
    )
  )
}

all_rows <- list()
fits <- list()

for (gene in c(MAIN_GENE, SUPP_GENE)) {
  log_line("\n", gene, ":", logfile = LOG)

  # --- unadjusted, for direct comparison ---
  ycol <- paste0(gene, "_log2tpm")
  du <- coh[stats::complete.cases(coh[[ycol]], coh$log10_SBS7), ]
  fu <- lm(as.formula(paste0(ycol, " ~ log10_SBS7")), data = du)
  su <- summary(fu)$coefficients
  ciu <- confint(fu)["log10_SBS7", ]
  all_rows[[length(all_rows) + 1]] <- tibble::tibble(
    gene = gene, model = "Unadjusted", purity_used = NA_character_,
    n = nrow(du),
    uv_beta = su["log10_SBS7", "Estimate"], uv_se = su["log10_SBS7", "Std. Error"],
    uv_ci_low = ciu[1], uv_ci_high = ciu[2],
    uv_beta_std = su["log10_SBS7", "Estimate"] * sd(du$log10_SBS7) / sd(du[[ycol]]),
    uv_p = su["log10_SBS7", "Pr(>|t|)"],
    purity_beta = NA_real_, purity_p = NA_real_,
    type_beta = NA_real_, type_p = NA_real_,
    adj_r2 = summary(fu)$adj.r.squared, max_vif = NA_real_
  )
  log_line(sprintf("  %-28s n=%d  beta=%+.4f [%+.4f, %+.4f]  p=%s  adjR2=%.4f",
                   "Unadjusted", nrow(du), su["log10_SBS7", "Estimate"],
                   ciu[1], ciu[2], fmt_p(su["log10_SBS7", "Pr(>|t|)"]),
                   summary(fu)$adj.r.squared), logfile = LOG)

  # --- primary adjusted model ---
  m <- fit_model(coh, gene, PURITY_PRIMARY,
                 paste0("Adjusted (", PURITY_PRIMARY, ")"))
  if (!is.null(m)) {
    all_rows[[length(all_rows) + 1]] <- m$row
    fits[[paste(gene, "primary")]] <- m$fit
    log_line(sprintf("  %-28s n=%d  beta=%+.4f [%+.4f, %+.4f]  p=%s  adjR2=%.4f  maxVIF=%.2f",
                     m$row$model, m$row$n, m$row$uv_beta, m$row$uv_ci_low,
                     m$row$uv_ci_high, fmt_p(m$row$uv_p), m$row$adj_r2,
                     m$row$max_vif), logfile = LOG)
  }

  # --- purity sensitivity ---
  for (p in PURITY_SENS) {
    if (!p %in% names(coh)) next
    m2 <- fit_model(coh, gene, p, paste0("Adjusted (", p, ")"))
    if (!is.null(m2)) {
      all_rows[[length(all_rows) + 1]] <- m2$row
      log_line(sprintf("  %-28s n=%d  beta=%+.4f [%+.4f, %+.4f]  p=%s",
                       m2$row$model, m2$row$n, m2$row$uv_beta, m2$row$uv_ci_low,
                       m2$row$uv_ci_high, fmt_p(m2$row$uv_p)), logfile = LOG)
    }
  }

  # --- rank-transformed sensitivity ---
  m3 <- fit_model(coh, gene, PURITY_PRIMARY, "Adjusted, rank-normal outcome",
                  transform = TRUE)
  if (!is.null(m3)) {
    all_rows[[length(all_rows) + 1]] <- m3$row
    log_line(sprintf("  %-28s n=%d  beta=%+.4f [%+.4f, %+.4f]  p=%s",
                     m3$row$model, m3$row$n, m3$row$uv_beta, m3$row$uv_ci_low,
                     m3$row$uv_ci_high, fmt_p(m3$row$uv_p)), logfile = LOG)
  }
}

mod_tbl <- dplyr::bind_rows(all_rows)

# ==============================================================================
# 4. Verdict
# ==============================================================================
hdr("4. Verdict per the analysis plan")

for (gene in c(MAIN_GENE, SUPP_GENE)) {
  r <- mod_tbl |> dplyr::filter(gene == !!gene,
                                model == paste0("Adjusted (", PURITY_PRIMARY, ")"))
  if (nrow(r) == 0) next
  verdict <- if (r$uv_p < 0.05) {
    "SURVIVES adjustment -- report as a finding"
  } else {
    "does NOT survive adjustment -- report as a negative result in the main text"
  }
  log_line(sprintf("  %-8s adjusted beta = %+.4f, p = %s  ->  %s",
                   gene, r$uv_beta, fmt_p(r$uv_p), verdict), logfile = LOG)
}

# ==============================================================================
# 5. Figures
# ==============================================================================
hdr("5. Figures")

purge_fig_prefix("figure3")
purge_fig_prefix("figureS3")

# --- Panel A: UV coefficient across model specifications ----------------------
panel_coefs <- function(gene) {
  # Short labels: the long model names were wide enough to collide with the
  # neighbouring panel.
  short <- c(
    "Unadjusted"                    = "Unadjusted",
    "Adjusted (CPE)"                = "+ CPE",
    "Adjusted (ESTIMATE)"           = "+ ESTIMATE",
    "Adjusted (LUMP)"               = "+ LUMP",
    "Adjusted (IHC)"                = "+ IHC",
    "Adjusted, rank-normal outcome" = "+ CPE, rank"
  )

  r <- mod_tbl |>
    dplyr::filter(gene == !!gene) |>
    dplyr::mutate(
      lab = dplyr::coalesce(unname(short[model]), model),
      lab = factor(lab, levels = rev(unname(short[unique(model)]))),
      sig = uv_p < 0.05
    )

  ggplot(r, aes(x = uv_beta, y = lab)) +
    geom_vline(xintercept = 0, linetype = "22", colour = PAL$grey_pt,
               linewidth = 0.4) +
    geom_errorbar(aes(xmin = uv_ci_low, xmax = uv_ci_high), orientation = "y",
                  width = 0.14, linewidth = 0.5, colour = PAL$ink) +
    geom_point(aes(fill = sig), size = 2.3, shape = 21, colour = PAL$ink,
               stroke = 0.4) +
    scale_fill_manual(values = c("TRUE" = PAL$uv_high, "FALSE" = "white"),
                      guide = "none") +
    scale_x_continuous(n.breaks = 4) +
    labs(x = expression(bold("UV coefficient"~beta)),
         y = NULL, subtitle = "Model specifications") +
    theme_pub(grid = TRUE)
}

# --- Panel B: partial residual plot -------------------------------------------
panel_partial <- function(gene) {
  ycol <- paste0(gene, "_log2tpm")
  dd <- coh[, c(ycol, "log10_SBS7", PURITY_PRIMARY, "sample_type")]
  names(dd) <- c("y", "uv", "purity", "sample_type")
  dd <- dd[stats::complete.cases(dd), ]

  # Residualise outcome AND exposure on the covariates, then plot one against
  # the other: the slope of this line is exactly the adjusted coefficient.
  ry <- resid(lm(y  ~ purity + sample_type, data = dd))
  rx <- resid(lm(uv ~ purity + sample_type, data = dd))
  pd <- tibble::tibble(rx = rx, ry = ry)

  r <- mod_tbl |> dplyr::filter(gene == !!gene,
                                model == paste0("Adjusted (", PURITY_PRIMARY, ")"))
  lab <- sprintf("beta == %+.3f", r$uv_beta)
  plab <- if (r$uv_p < 0.0001) "italic(P) < 0.0001" else
    sprintf("italic(P) == %s", formatC(r$uv_p, format = "f",
                                       digits = ifelse(r$uv_p < 0.01, 4, 3)))

  ggplot(pd, aes(x = rx, y = ry)) +
    geom_point(size = 0.9, alpha = 0.55, colour = PAL$hist, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = PAL$ink,
                fill = "grey70", linewidth = 0.6, alpha = 0.25) +
    annotate("text", x = max(pd$rx) * 0.98, y = max(pd$ry) * 0.98, label = lab,
             parse = TRUE, hjust = 1, vjust = 1, size = 2.9, colour = PAL$ink) +
    annotate("text", x = max(pd$rx) * 0.98, y = max(pd$ry) * 0.84, label = plab,
             parse = TRUE, hjust = 1, vjust = 1, size = 2.9, colour = PAL$ink) +
    labs(x = expression(bold("UV burden, adjusted")),
         y = bquote(bold(italic(.(gene))*", adjusted")),
         subtitle = "Partial association") +
    theme_pub()
}

# --- Panel C: purity against UV burden ----------------------------------------
panel_confound <- function() {
  ggplot(coh, aes(x = log10_SBS7, y = .data[[PURITY_PRIMARY]])) +
    geom_point(size = 0.9, alpha = 0.55, colour = PAL$accent, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = PAL$ink,
                fill = "grey70", linewidth = 0.6, alpha = 0.25) +
    labs(x = expression(bold(log[10]*"(SBS7 + 1)")),
         y = paste0("Purity (", PURITY_PRIMARY, ")"),
         subtitle = "Confounder check") +
    theme_pub()
}

build_fig <- function(gene) {
  # The coefficient panel carries the row labels, so it gets extra width.
  (panel_coefs(gene) | panel_partial(gene) | panel_confound()) +
    plot_layout(widths = c(1.25, 1, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = FIG$base_size + 3, face = "bold"))
}

save_fig(build_fig(MAIN_GENE), "figure3_rassf9_adjusted",
         width = FIG$width_double, height = 85)
save_fig(build_fig(SUPP_GENE), "figureS3_rassf10_adjusted",
         width = FIG$width_double, height = 85)

# ==============================================================================
# 6. Tables
# ==============================================================================
hdr("6. Results tables")

main_tbl <- mod_tbl |> dplyr::filter(gene == MAIN_GENE)
supp_tbl <- mod_tbl |> dplyr::filter(gene == SUPP_GENE)

readr::write_csv(main_tbl,
                 file.path(PATHS$results, "table_main_06_rassf9_adjusted.csv"))
readr::write_csv(supp_tbl,
                 file.path(PATHS$results, "table_supp_06_rassf10_adjusted.csv"))

log_line("\nMAIN (", MAIN_GENE, "):", logfile = LOG)
cat(capture.output(print(as.data.frame(main_tbl), row.names = FALSE, digits = 3)),
    sep = "\n", file = LOG, append = TRUE)
print(as.data.frame(main_tbl), digits = 3)

log_line("\nSUPPLEMENTARY (", SUPP_GENE, "):", logfile = LOG)
cat(capture.output(print(as.data.frame(supp_tbl), row.names = FALSE, digits = 3)),
    sep = "\n", file = LOG, append = TRUE)
print(as.data.frame(supp_tbl), digits = 3)

hdr("7. Warnings")
report_warnings(LOG)

write_session_info("06_adjusted")

hdr("DONE")
log_line("Figure 3 (main) : figures/figure3_rassf9_adjusted.*", logfile = LOG)
log_line("Figure S3 (supp): figures/figureS3_rassf10_adjusted.*", logfile = LOG)
log_line("Full log: ", LOG, logfile = LOG)
