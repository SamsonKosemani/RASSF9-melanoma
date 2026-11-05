# ==============================================================================
# 04_group_comparison.R
#
# Question 2: does RASSF9 expression differ between melanomas with high versus
# low UV mutational burden? RASSF10 is analysed in parallel but kept entirely
# separate in the outputs.
#
# OUTPUT STRUCTURE (fixed here so the manuscript never needs re-splitting):
#   RASSF9  -> figures/figure1_rassf9_uv_groups.*   + results/table_main_*.csv
#   RASSF10 -> figures/figureS1_rassf10_uv_groups.* + results/table_supp_*.csv
#
# Pre-specified (docs/analysis-plan.md section 5):
#   - groups = outer tertiles of SBS7_total; middle tertile not analysed
#   - test   = Wilcoxon rank-sum, two-sided, alpha = 0.05
#   - effect = Cohen's d
#   - sensitivity: full cohort incl. SBS7 = 0; signature accuracy >= 0.80
#
# Input:  data/processed/analysis_cohort.csv, full_cohort.csv
# Output: see above, plus results/logs/04_groups.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(patchwork)
  library(scales)
})

LOG <- file.path(PATHS$logs, "04_groups.log")
cat("Group comparison —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

MAIN_GENE <- "RASSF9"    # the spine of the manuscript
SUPP_GENE <- "RASSF10"   # parallel observation, supplementary

# ==============================================================================
# 1. Load
# ==============================================================================
hdr("1. Cohort")

coh  <- readr::read_csv(file.path(PATHS$processed, "analysis_cohort.csv"))
full <- readr::read_csv(file.path(PATHS$processed, "full_cohort.csv"))

log_line("Primary cohort: n = ", nrow(coh), logfile = LOG)
log_line("Full cohort:    n = ", nrow(full), logfile = LOG)

# ==============================================================================
# 2. UV stratification
# ==============================================================================
hdr("2. UV stratification")

assign_uv_group <- function(d) {
  tert <- quantile(d$SBS7_total, c(1/3, 2/3), na.rm = TRUE)
  out <- d |>
    dplyr::mutate(
      uv_group = dplyr::case_when(
        SBS7_total <= tert[1] ~ "UV-low",
        SBS7_total >= tert[2] ~ "UV-high",
        TRUE                  ~ "Middle"
      ),
      uv_group = factor(uv_group, levels = c("UV-low", "Middle", "UV-high"))
    )
  attr(out, "cutoffs") <- tert
  out
}

coh  <- assign_uv_group(coh)
full <- assign_uv_group(full)
cut_primary <- attr(coh, "cutoffs")

log_line("Tertile cutoffs (primary cohort): lower = ", round(cut_primary[1]),
         ", upper = ", round(cut_primary[2]), logfile = LOG)

grp_n <- as.data.frame(table(coh$uv_group))
names(grp_n) <- c("group", "n")
cat(capture.output(print(grp_n, row.names = FALSE)),
    sep = "\n", file = LOG, append = TRUE)
print(grp_n)

# ==============================================================================
# 3. Tests
# ==============================================================================
hdr("3. Wilcoxon rank-sum tests")

cohens_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  s  <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  (mean(x) - mean(y)) / s
}

# Hedges & Olkin large-sample standard error for d.
cohens_d_ci <- function(d, n1, n2, level = 0.95) {
  se <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))
  z  <- qnorm(1 - (1 - level) / 2)
  c(lower = d - z * se, upper = d + z * se, se = se)
}

run_comparison <- function(d, gene, label) {
  col <- paste0(gene, "_log2tpm")
  if (!col %in% names(d)) return(NULL)

  lo <- d[[col]][d$uv_group == "UV-low"];  lo <- lo[!is.na(lo)]
  hi <- d[[col]][d$uv_group == "UV-high"]; hi <- hi[!is.na(hi)]
  if (length(lo) < 3 || length(hi) < 3) return(NULL)

  w  <- suppressWarnings(wilcox.test(hi, lo, exact = FALSE))
  dd <- cohens_d(hi, lo)
  ci <- cohens_d_ci(dd, length(hi), length(lo))

  res <- tibble::tibble(
    analysis    = label,
    gene        = gene,
    n_low       = length(lo),
    n_high      = length(hi),
    median_low  = median(lo),
    median_high = median(hi),
    W           = unname(w$statistic),
    p_value     = w$p.value,
    cohens_d    = dd,
    d_lower     = unname(ci["lower"]),
    d_upper     = unname(ci["upper"]),
    direction   = ifelse(median(hi) < median(lo), "lower in UV-high",
                         "higher in UV-high")
  )

  log_line(sprintf(
    "  %-22s %-8s n=%d/%d  median %.3f vs %.3f  W=%.0f  p=%s  d=%.3f [%.2f, %.2f]  (%s)",
    label, gene, length(lo), length(hi), median(lo), median(hi),
    w$statistic, fmt_p(w$p.value), dd, ci["lower"], ci["upper"],
    res$direction), logfile = LOG)

  res
}

coh_acc <- coh |> dplyr::filter(accuracy >= 0.80) |> assign_uv_group()

cohort_sets <- list(
  "Primary"        = coh,
  "S1 full cohort" = full,
  "S2 accuracy"    = coh_acc
)

results <- list()
for (gene in c(MAIN_GENE, SUPP_GENE)) {
  log_line("\n", gene, ":", logfile = LOG)
  for (lab in names(cohort_sets)) {
    results[[length(results) + 1]] <- run_comparison(cohort_sets[[lab]], gene, lab)
  }
}
res_tbl <- dplyr::bind_rows(results)

# ==============================================================================
# 4. Shared plotting pieces
# ==============================================================================
hdr("4. Figures")

purge_fig_prefix("figure1")
purge_fig_prefix("figureS1")

pri <- coh |> dplyr::filter(uv_group != "Middle") |> droplevels()

grp_labels <- pri |>
  dplyr::count(uv_group) |>
  dplyr::mutate(lab = paste0(uv_group, "\n(n = ", n, ")"))
lab_vec <- setNames(grp_labels$lab, grp_labels$uv_group)

# --- UV burden distribution (context panel, shared by both figures) -----------
panel_distribution <- function() {
  ggplot(coh, aes(x = log10_SBS7)) +
    geom_histogram(bins = 32, fill = PAL$hist, colour = "white",
                   linewidth = 0.15) +
    geom_vline(xintercept = log10(cut_primary + 1), linetype = "22",
               colour = PAL$ink, linewidth = 0.5) +
    annotate("text", x = log10(cut_primary[1] + 1), y = Inf,
             label = round(cut_primary[1]), vjust = 1.4, hjust = 1.15,
             size = 2.7, colour = PAL$ink, fontface = "bold") +
    annotate("text", x = log10(cut_primary[2] + 1), y = Inf,
             label = round(cut_primary[2]), vjust = 1.4, hjust = -0.15,
             size = 2.7, colour = PAL$ink, fontface = "bold") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.04))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = expression(bold(log[10]*"(SBS7 mutations + 1)")), y = "Tumours",
         subtitle = paste0("UV mutational burden (n = ", nrow(coh), ")")) +
    theme_pub()
}

panel_strata <- function() {
  ggplot(coh, aes(x = uv_group, y = SBS7_total, fill = uv_group)) +
    geom_boxplot(outlier.size = 0.4, outlier.colour = PAL$grey_pt,
                 linewidth = 0.4, width = 0.62, colour = PAL$ink) +
    scale_y_log10(labels = label_number(big.mark = ","),
                  expand = expansion(mult = c(0.05, 0.08))) +
    scale_fill_manual(values = PAL_UV, guide = "none") +
    labs(x = NULL, y = "SBS7 mutations", subtitle = "Tertile strata") +
    theme_pub(grid = TRUE)
}

# --- Expression by UV group ---------------------------------------------------
panel_expression <- function(gene) {
  col <- paste0(gene, "_log2tpm")
  r <- res_tbl |> dplyr::filter(analysis == "Primary", gene == !!gene)
  ymax <- max(pri[[col]], na.rm = TRUE)
  ybar <- ymax * 1.05

  ggplot(pri, aes(x = uv_group, y = .data[[col]], fill = uv_group)) +
    geom_violin(trim = TRUE, linewidth = 0.3, colour = NA, alpha = 0.30,
                scale = "width") +
    geom_jitter(width = 0.11, size = 0.5, alpha = 0.40, colour = PAL$grey_pt,
                stroke = 0) +
    geom_boxplot(width = 0.20, outlier.shape = NA, linewidth = 0.45,
                 colour = PAL$ink, fill = "white", alpha = 0.95) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 1.5,
                 fill = "white", colour = PAL$ink, stroke = 0.4) +
    annotate("segment", x = 1, xend = 2, y = ybar, yend = ybar,
             linewidth = 0.45, colour = PAL$ink) +
    annotate("segment", x = 1, xend = 1, y = ybar, yend = ybar * 0.975,
             linewidth = 0.45, colour = PAL$ink) +
    annotate("segment", x = 2, xend = 2, y = ybar, yend = ybar * 0.975,
             linewidth = 0.45, colour = PAL$ink) +
    annotate("text", x = 1.5, y = ybar * 1.02, label = fmt_p_plot(r$p_value),
             parse = TRUE, size = 2.9, colour = PAL$ink, vjust = 0) +
    scale_fill_manual(values = PAL_UV, guide = "none") +
    scale_x_discrete(labels = lab_vec) +
    coord_cartesian(ylim = c(0, ymax * 1.20)) +
    labs(x = NULL, y = bquote(bold(italic(.(gene))~log[2]*"(TPM+1)")),
         subtitle = bquote(italic(.(gene))~"expression")) +
    theme_pub()
}

# --- Effect size across cohort definitions ------------------------------------
panel_effects <- function(gene) {
  r <- res_tbl |>
    dplyr::filter(gene == !!gene) |>
    dplyr::mutate(analysis = factor(analysis,
                                    levels = rev(names(cohort_sets))))
  ggplot(r, aes(x = cohens_d, y = analysis)) +
    geom_vline(xintercept = 0, linetype = "22", colour = PAL$grey_pt,
               linewidth = 0.4) +
    geom_errorbar(aes(xmin = d_lower, xmax = d_upper), orientation = "y",
                  width = 0.14, linewidth = 0.5, colour = PAL$ink) +
    geom_point(size = 2.2, shape = 21, fill = PAL$uv_high, colour = PAL$ink,
               stroke = 0.4) +
    labs(x = "Cohen's d (UV-high \u2212 UV-low)", y = NULL,
         subtitle = "Effect size, 95% CI") +
    theme_pub(grid = TRUE)
}

# ==============================================================================
# 5. Figure 1 — RASSF9 (main)
# ==============================================================================
fig1 <- (panel_distribution() | panel_strata()) /
        (panel_expression(MAIN_GENE) | panel_effects(MAIN_GENE)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = FIG$base_size + 3, face = "bold"))

save_fig(fig1, "figure1_rassf9_uv_groups",
         width = FIG$width_double, height = 135)

# ==============================================================================
# 6. Supplementary Figure S1 — RASSF10
# ==============================================================================
figS1 <- (panel_expression(SUPP_GENE) | panel_effects(SUPP_GENE)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = FIG$base_size + 3, face = "bold"))

save_fig(figS1, "figureS1_rassf10_uv_groups",
         width = FIG$width_1_5, height = 72)

# ==============================================================================
# 7. Results tables, split by gene
# ==============================================================================
hdr("5. Results tables")

main_tbl <- res_tbl |> dplyr::filter(gene == MAIN_GENE)
supp_tbl <- res_tbl |> dplyr::filter(gene == SUPP_GENE)

readr::write_csv(main_tbl,
                 file.path(PATHS$results, "table_main_04_rassf9_groups.csv"))
readr::write_csv(supp_tbl,
                 file.path(PATHS$results, "table_supp_04_rassf10_groups.csv"))

log_line("\nMAIN (", MAIN_GENE, "):", logfile = LOG)
cat(capture.output(print(as.data.frame(main_tbl), row.names = FALSE)),
    sep = "\n", file = LOG, append = TRUE)
print(as.data.frame(main_tbl))

log_line("\nSUPPLEMENTARY (", SUPP_GENE, "):", logfile = LOG)
cat(capture.output(print(as.data.frame(supp_tbl), row.names = FALSE)),
    sep = "\n", file = LOG, append = TRUE)
print(as.data.frame(supp_tbl))

log_line("\nWrote: results/table_main_04_rassf9_groups.csv", logfile = LOG)
log_line("Wrote: results/table_supp_04_rassf10_groups.csv", logfile = LOG)

hdr("6. Warnings")
report_warnings(LOG)

write_session_info("04_groups")

hdr("DONE")
log_line("Figure 1 (main)      : figures/figure1_rassf9_uv_groups.*", logfile = LOG)
log_line("Figure S1 (supp)     : figures/figureS1_rassf10_uv_groups.*", logfile = LOG)
log_line("Full log: ", LOG, logfile = LOG)
