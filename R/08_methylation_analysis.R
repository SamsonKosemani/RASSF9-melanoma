# ==============================================================================
# 08_methylation_analysis.R
#
# Question 5: does promoter methylation inversely correlate with expression?
#
# This tests the mechanism behind what scripts 04-06 found. It also tests
# Gustavo's stated premise directly: if RASSF10 is epigenetically silenced in
# metastatic melanoma, its promoter should be MORE methylated in metastases --
# and script 06 already found RASSF10 expression significantly LOWER in
# metastases (beta = -0.347, p < 0.0001), so there is a concrete prediction.
#
# PRE-SPECIFIED (D12-D14, fixed before any association was examined):
#   D12 Promoter summary = MEAN beta across promoter-window probes per sample.
#       Per-probe results are reported alongside, so the mean cannot hide a
#       probe pointing the other way.
#   D13 Primary test = Spearman, promoter-mean beta vs log2(TPM+1). Per-probe
#       Spearman with Benjamini-Hochberg correction WITHIN each gene.
#   D14 Adjusted model = log2(TPM+1) ~ promoter beta + purity (CPE) +
#       sample type. Methylation is measured in the same tissue as expression,
#       so purity confounding applies here as it did in script 06.
#
# Secondary, also pre-specified: promoter methylation by sample type
# (primary vs metastatic), and promoter methylation vs UV mutational burden.
#
# Input:  data/processed/methylation_betas.csv, methylation_probes.csv,
#         analysis_cohort.csv
# Output: figures/figure4_rassf9_methylation.*, figureS4_rassf10_methylation.*
#         results/table_main_08_rassf9_methylation.csv
#         results/table_supp_08_rassf10_methylation.csv
#         results/logs/08_methylation.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(patchwork)
  library(scales)
})

LOG <- file.path(PATHS$logs, "08_methylation.log")
cat("Methylation analysis —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

MAIN_GENE <- "RASSF9"
SUPP_GENE <- "RASSF10"
PURITY    <- "CPE"

# ==============================================================================
# 1. Assemble
# ==============================================================================
hdr("1. Data")

probes <- readr::read_csv(file.path(PATHS$processed, "methylation_probes.csv"))
betas  <- readr::read_csv(file.path(PATHS$processed, "methylation_betas.csv"))
coh    <- readr::read_csv(file.path(PATHS$processed, "analysis_cohort.csv"))

coh <- coh |>
  dplyr::mutate(
    sample_type = dplyr::case_when(
      grepl("^Primary", sample_type, ignore.case = TRUE)  ~ "Primary",
      grepl("Metastatic", sample_type, ignore.case = TRUE) ~ "Metastatic",
      TRUE ~ NA_character_),
    sample_type = factor(sample_type, levels = c("Primary", "Metastatic"))
  )

d <- dplyr::inner_join(coh, betas, by = "patient")
log_line("Patients with expression + methylation: ", nrow(d), logfile = LOG)

prom_list <- list()
for (g in c(MAIN_GENE, SUPP_GENE)) {
  pr <- probes |>
    dplyr::filter(gene == g, promoter, on_array) |>
    dplyr::arrange(dist_tss)
  pr <- pr[pr$probe %in% names(d), ]
  prom_list[[g]] <- pr
  log_line("  ", g, ": ", nrow(pr), " promoter probes available", logfile = LOG)

  # D12: promoter mean, computed per sample across available probes.
  mat <- as.matrix(d[, pr$probe, drop = FALSE])
  d[[paste0(g, "_prom_beta")]] <- rowMeans(mat, na.rm = TRUE)

  log_line(sprintf("    promoter mean beta: median [IQR] = %s, range [%.3f, %.3f]",
                   fmt_median_iqr(d[[paste0(g, "_prom_beta")]], 3),
                   min(d[[paste0(g, "_prom_beta")]], na.rm = TRUE),
                   max(d[[paste0(g, "_prom_beta")]], na.rm = TRUE)),
           logfile = LOG)
}

# ==============================================================================
# 2. Per-probe correlations  (D13)
# ==============================================================================
hdr("2. Per-probe methylation vs expression")

probe_rows <- list()
for (g in c(MAIN_GENE, SUPP_GENE)) {
  ycol <- paste0(g, "_log2tpm")
  pr <- prom_list[[g]]
  out <- list()
  for (k in seq_len(nrow(pr))) {
    p <- pr$probe[k]
    ok <- stats::complete.cases(d[[p]], d[[ycol]])
    ct <- suppressWarnings(cor.test(d[[p]][ok], d[[ycol]][ok],
                                    method = "spearman"))
    out[[k]] <- tibble::tibble(
      gene = g, probe = p, dist_tss = pr$dist_tss[k], n = sum(ok),
      median_beta = median(d[[p]], na.rm = TRUE),
      rho = unname(ct$estimate), p_value = ct$p.value
    )
  }
  o <- dplyr::bind_rows(out)
  o$p_adj <- p.adjust(o$p_value, method = "BH")   # BH within gene
  probe_rows[[g]] <- o

  log_line("\n  ", g, " (BH-corrected within gene):", logfile = LOG)
  for (k in seq_len(nrow(o))) {
    log_line(sprintf("    %-12s %+6d bp  median beta=%.3f  rho=%+.3f  p=%s  p.adj=%s",
                     o$probe[k], o$dist_tss[k], o$median_beta[k], o$rho[k],
                     fmt_p(o$p_value[k]), fmt_p(o$p_adj[k])), logfile = LOG)
  }
  log_line(sprintf("    %d of %d probes significant after BH correction",
                   sum(o$p_adj < 0.05), nrow(o)), logfile = LOG)
}
probe_tbl <- dplyr::bind_rows(probe_rows)

# ==============================================================================
# 3. Promoter-mean association  (D13, D14)
# ==============================================================================
hdr("3. Promoter mean: unadjusted and adjusted")

summary_rows <- list()
for (g in c(MAIN_GENE, SUPP_GENE)) {
  ycol <- paste0(g, "_log2tpm")
  bcol <- paste0(g, "_prom_beta")

  ok <- stats::complete.cases(d[[bcol]], d[[ycol]])
  sp <- suppressWarnings(cor.test(d[[bcol]][ok], d[[ycol]][ok],
                                  method = "spearman"))

  dd <- d[, c(ycol, bcol, PURITY, "sample_type", "log10_SBS7")]
  names(dd) <- c("y", "beta", "purity", "sample_type", "uv")
  dd <- dd[stats::complete.cases(dd), ]

  fit_u <- lm(y ~ beta, data = dd)
  fit_a <- lm(y ~ beta + purity + sample_type, data = dd)
  cu <- summary(fit_u)$coefficients
  ca <- summary(fit_a)$coefficients
  ciu <- confint(fit_u)["beta", ]
  cia <- confint(fit_a)["beta", ]

  # Secondary: is the promoter more methylated in metastases?
  wt <- suppressWarnings(wilcox.test(beta ~ sample_type, data = dd,
                                     exact = FALSE))
  med_p <- median(dd$beta[dd$sample_type == "Primary"], na.rm = TRUE)
  med_m <- median(dd$beta[dd$sample_type == "Metastatic"], na.rm = TRUE)

  # Secondary: does methylation track UV burden?
  uv_ct <- suppressWarnings(cor.test(dd$beta, dd$uv, method = "spearman"))

  summary_rows[[g]] <- tibble::tibble(
    gene = g, n = nrow(dd),
    spearman_rho = unname(sp$estimate), spearman_p = sp$p.value,
    lm_beta_unadj = cu["beta", "Estimate"],
    lm_unadj_low = ciu[1], lm_unadj_high = ciu[2],
    lm_unadj_p = cu["beta", "Pr(>|t|)"],
    lm_beta_adj = ca["beta", "Estimate"],
    lm_adj_low = cia[1], lm_adj_high = cia[2],
    lm_adj_p = ca["beta", "Pr(>|t|)"],
    adj_r2 = summary(fit_a)$adj.r.squared,
    median_beta_primary = med_p, median_beta_metastatic = med_m,
    beta_by_type_p = wt$p.value,
    beta_vs_uv_rho = unname(uv_ct$estimate), beta_vs_uv_p = uv_ct$p.value
  )

  s <- summary_rows[[g]]
  log_line("\n  ", g, ":", logfile = LOG)
  log_line(sprintf("    Spearman (promoter mean vs expression): rho = %+.3f, p = %s",
                   s$spearman_rho, fmt_p(s$spearman_p)), logfile = LOG)
  log_line(sprintf("    Linear, unadjusted: beta = %+.3f [%+.3f, %+.3f], p = %s",
                   s$lm_beta_unadj, s$lm_unadj_low, s$lm_unadj_high,
                   fmt_p(s$lm_unadj_p)), logfile = LOG)
  log_line(sprintf("    Linear, adjusted  : beta = %+.3f [%+.3f, %+.3f], p = %s, adjR2 = %.3f",
                   s$lm_beta_adj, s$lm_adj_low, s$lm_adj_high,
                   fmt_p(s$lm_adj_p), s$adj_r2), logfile = LOG)
  log_line(sprintf("    Promoter beta by sample type: primary %.3f vs metastatic %.3f, p = %s",
                   s$median_beta_primary, s$median_beta_metastatic,
                   fmt_p(s$beta_by_type_p)), logfile = LOG)
  log_line(sprintf("    Promoter beta vs UV burden: rho = %+.3f, p = %s",
                   s$beta_vs_uv_rho, fmt_p(s$beta_vs_uv_p)), logfile = LOG)
}
sum_tbl <- dplyr::bind_rows(summary_rows)

hdr("4. Interpretation guide")
log_line("  A NEGATIVE methylation-expression correlation supports silencing.",
         logfile = LOG)
log_line("  HIGHER promoter beta in metastases would support Gustavo's premise",
         logfile = LOG)
log_line("  that RASSF10 is epigenetically silenced in metastatic melanoma.",
         logfile = LOG)
log_line("  Note median beta values: a promoter sitting near beta = 0 is",
         logfile = LOG)
log_line("  unmethylated in most tumours, which limits how much expression",
         logfile = LOG)
log_line("  variation methylation could possibly explain.", logfile = LOG)

# ==============================================================================
# 5. Figures
# ==============================================================================
hdr("5. Figures")

purge_fig_prefix("figure4")
purge_fig_prefix("figureS4")

panel_probewise <- function(g) {
  o <- probe_rows[[g]] |> dplyr::mutate(sig = p_adj < 0.05)
  ggplot(o, aes(x = dist_tss, y = rho)) +
    geom_hline(yintercept = 0, linetype = "22", colour = PAL$grey_pt,
               linewidth = 0.4) +
    geom_segment(aes(xend = dist_tss, yend = 0), linewidth = 0.4,
                 colour = PAL$grey_pt) +
    geom_point(aes(fill = sig), size = 2.2, shape = 21, colour = PAL$ink,
               stroke = 0.4) +
    scale_fill_manual(values = c("TRUE" = PAL$uv_high, "FALSE" = "white"),
                      guide = "none") +
    labs(x = "Distance from TSS (bp)", y = expression(bold("Spearman"~rho)),
         subtitle = "Per-probe association") +
    theme_pub(grid = TRUE)
}

panel_scatter <- function(g) {
  ycol <- paste0(g, "_log2tpm"); bcol <- paste0(g, "_prom_beta")
  s <- sum_tbl |> dplyr::filter(gene == g)
  lab  <- sprintf("rho == %+.3f", s$spearman_rho)
  plab <- if (s$spearman_p < 0.0001) "italic(P) < 0.0001" else
    sprintf("italic(P) == %s", formatC(s$spearman_p, format = "f",
                                       digits = ifelse(s$spearman_p < 0.01, 4, 3)))
  pd <- d[stats::complete.cases(d[[bcol]], d[[ycol]]), ]

  ggplot(pd, aes(x = .data[[bcol]], y = .data[[ycol]])) +
    geom_point(size = 0.9, alpha = 0.6, colour = PAL$hist, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = PAL$ink,
                fill = "grey70", linewidth = 0.6, alpha = 0.25) +
    annotate("text", x = max(pd[[bcol]]) * 0.98, y = max(pd[[ycol]]) * 0.98,
             label = lab, parse = TRUE, hjust = 1, vjust = 1, size = 2.9,
             colour = PAL$ink) +
    annotate("text", x = max(pd[[bcol]]) * 0.98, y = max(pd[[ycol]]) * 0.86,
             label = plab, parse = TRUE, hjust = 1, vjust = 1, size = 2.9,
             colour = PAL$ink) +
    labs(x = "Promoter methylation (mean beta)",
         y = bquote(bold(italic(.(g))~log[2]*"(TPM+1)")),
         subtitle = "Methylation vs expression") +
    theme_pub()
}

panel_bytype <- function(g) {
  bcol <- paste0(g, "_prom_beta")
  s <- sum_tbl |> dplyr::filter(gene == g)
  pd <- d[!is.na(d[[bcol]]) & !is.na(d$sample_type), ]
  ymax <- max(pd[[bcol]], na.rm = TRUE)

  ggplot(pd, aes(x = sample_type, y = .data[[bcol]], fill = sample_type)) +
    geom_violin(trim = TRUE, linewidth = 0.3, colour = NA, alpha = 0.30,
                scale = "width") +
    geom_jitter(width = 0.11, size = 0.5, alpha = 0.40, colour = PAL$grey_pt,
                stroke = 0) +
    geom_boxplot(width = 0.20, outlier.shape = NA, linewidth = 0.45,
                 colour = PAL$ink, fill = "white", alpha = 0.95) +
    annotate("segment", x = 1, xend = 2, y = ymax * 1.05, yend = ymax * 1.05,
             linewidth = 0.45, colour = PAL$ink) +
    annotate("text", x = 1.5, y = ymax * 1.08,
             label = fmt_p_plot(s$beta_by_type_p), parse = TRUE, size = 2.9,
             colour = PAL$ink, vjust = 0) +
    scale_fill_manual(values = c("Primary" = PAL$uv_low,
                                 "Metastatic" = PAL$uv_high), guide = "none") +
    coord_cartesian(ylim = c(0, ymax * 1.20)) +
    labs(x = NULL, y = "Promoter methylation (beta)",
         subtitle = "By tumour site") +
    theme_pub()
}

build_fig <- function(g) {
  (panel_probewise(g) | panel_scatter(g) | panel_bytype(g)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = FIG$base_size + 3, face = "bold"))
}

save_fig(build_fig(MAIN_GENE), "figure4_rassf9_methylation",
         width = FIG$width_double, height = 80)
save_fig(build_fig(SUPP_GENE), "figureS4_rassf10_methylation",
         width = FIG$width_double, height = 80)

# ==============================================================================
# 6. Tables
# ==============================================================================
hdr("6. Results tables")

readr::write_csv(
  dplyr::bind_rows(
    sum_tbl |> dplyr::filter(gene == MAIN_GENE) |> dplyr::mutate(level = "promoter mean"),
    probe_tbl |> dplyr::filter(gene == MAIN_GENE) |> dplyr::mutate(level = "per probe")
  ),
  file.path(PATHS$results, "table_main_08_rassf9_methylation.csv"))

readr::write_csv(
  dplyr::bind_rows(
    sum_tbl |> dplyr::filter(gene == SUPP_GENE) |> dplyr::mutate(level = "promoter mean"),
    probe_tbl |> dplyr::filter(gene == SUPP_GENE) |> dplyr::mutate(level = "per probe")
  ),
  file.path(PATHS$results, "table_supp_08_rassf10_methylation.csv"))

log_line("Wrote: results/table_main_08_rassf9_methylation.csv", logfile = LOG)
log_line("Wrote: results/table_supp_08_rassf10_methylation.csv", logfile = LOG)

hdr("7. Warnings")
report_warnings(LOG)

write_session_info("08_methylation")

hdr("DONE")
log_line("Figure 4 (main) : figures/figure4_rassf9_methylation.*", logfile = LOG)
log_line("Figure S4 (supp): figures/figureS4_rassf10_methylation.*", logfile = LOG)
log_line("Full log: ", LOG, logfile = LOG)
