# ==============================================================================
# 09_aggressiveness.R
#
# Question 6: is lower expression associated with more aggressive melanoma?
#
# IMPORTANT CONTEXT. Gustavo recalls an earlier analysis showing that lower
# RASSF9 was associated with more aggressive disease. No such analysis exists
# in the project record -- what existed was the UV-burden association and a
# melanoma-vs-normal comparison, neither of which tested aggressiveness. This
# script runs the analysis properly for the first time. It is included because
# it was asked for, NOT because a result is anticipated, and a null result here
# is a perfectly good answer that should be reported as such.
#
# PRE-SPECIFIED (D15, fixed before any of these associations were examined):
#   Four aggressiveness measures, each tested once per gene:
#     A1 tumour site        primary vs metastatic      Wilcoxon
#     A2 AJCC stage         ordered groups             Kruskal-Wallis + trend
#     A3 Breslow thickness  continuous                 Spearman
#     A4 overall survival   continuous expression      Cox, unadjusted + adjusted
#   Benjamini-Hochberg correction across the four tests WITHIN each gene,
#   because these are four looks at one question, not four questions.
#   Adjusted survival model includes CPE purity and sample type, as elsewhere.
#
# Input:  data/processed/analysis_cohort.csv, data/raw/skcm_star_counts.rds
# Output: figures/figure5_rassf9_aggressiveness.*, figureS5_rassf10_*
#         results/table_main_09_rassf9_aggressiveness.csv
#         results/table_supp_09_rassf10_aggressiveness.csv
#         results/logs/09_aggressiveness.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(survival)
  library(patchwork)
  library(scales)
})
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))

LOG <- file.path(PATHS$logs, "09_aggressiveness.log")
cat("Aggressiveness analysis —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

MAIN_GENE <- "RASSF9"
SUPP_GENE <- "RASSF10"
PURITY    <- "CPE"

# ==============================================================================
# 1. Clinical variables
# ==============================================================================
hdr("1. Clinical data")

coh <- readr::read_csv(file.path(PATHS$processed, "analysis_cohort.csv"))
se  <- readRDS(file.path(PATHS$raw, "skcm_star_counts.rds"))
cd  <- as.data.frame(colData(se))

safe_chr <- function(x) {
  if (is.null(x)) return(NA_character_)
  if (is.list(x)) x <- vapply(x, function(e)
    if (length(e)) as.character(e)[1] else NA_character_, character(1))
  as.character(x)
}
safe_num <- function(x) suppressWarnings(as.numeric(safe_chr(x)))

pick <- function(patterns, numeric = FALSE) {
  for (p in patterns) {
    hit <- grep(p, names(cd), ignore.case = TRUE, value = TRUE)
    if (length(hit)) {
      log_line("    using column: ", hit[1], logfile = LOG)
      return(if (numeric) safe_num(cd[[hit[1]]]) else safe_chr(cd[[hit[1]]]))
    }
  }
  log_line("    [FLAG] none of: ", paste(patterns, collapse = ", "),
           logfile = LOG)
  rep(if (numeric) NA_real_ else NA_character_, nrow(cd))
}

log_line("  stage:", logfile = LOG)
stage_raw <- pick(c("^ajcc_pathologic_stage$", "ajcc.*stage", "tumor_stage"))
log_line("  Breslow:", logfile = LOG)
breslow   <- pick(c("breslow"), numeric = TRUE)
log_line("  vital status:", logfile = LOG)
vital     <- pick(c("^vital_status$", "vital"))
log_line("  days to death:", logfile = LOG)
d_death   <- pick(c("days_to_death"), numeric = TRUE)
log_line("  days to last follow-up:", logfile = LOG)
d_follow  <- pick(c("days_to_last_follow", "days_to_last_known"), numeric = TRUE)

clin <- tibble::tibble(
  patient   = tcga_patient(safe_chr(if ("barcode" %in% names(cd)) cd$barcode
                                    else colnames(se))),
  stage_raw = stage_raw,
  breslow   = breslow,
  vital     = vital,
  d_death   = d_death,
  d_follow  = d_follow
) |>
  dplyr::group_by(patient) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    # Collapse sub-stages: I/IA/IB -> I, and so on.
    stage = dplyr::case_when(
      grepl("^Stage\\s*IV", stage_raw)        ~ "IV",
      grepl("^Stage\\s*III", stage_raw)       ~ "III",
      grepl("^Stage\\s*II", stage_raw)        ~ "II",
      grepl("^Stage\\s*I([^IV]|$)", stage_raw) ~ "I",
      grepl("^Stage\\s*0", stage_raw)         ~ "0",
      TRUE                                    ~ NA_character_
    ),
    stage = factor(stage, levels = c("0", "I", "II", "III", "IV")),
    os_event = dplyr::case_when(
      grepl("Dead", vital, ignore.case = TRUE)  ~ 1L,
      grepl("Alive", vital, ignore.case = TRUE) ~ 0L,
      TRUE ~ NA_integer_),
    os_time = dplyr::coalesce(d_death, d_follow)
  )

d <- dplyr::inner_join(coh, clin, by = "patient") |>
  dplyr::mutate(
    sample_type = dplyr::case_when(
      grepl("^Primary", sample_type, ignore.case = TRUE)  ~ "Primary",
      grepl("Metastatic", sample_type, ignore.case = TRUE) ~ "Metastatic",
      TRUE ~ NA_character_),
    sample_type = factor(sample_type, levels = c("Primary", "Metastatic"))
  )

log_line("\n  Cohort with clinical data: ", nrow(d), logfile = LOG)
log_line("  Stage available:   ", sum(!is.na(d$stage)), logfile = LOG)
log_line("  Breslow available: ", sum(!is.na(d$breslow)), logfile = LOG)
log_line("  Survival usable:   ",
         sum(!is.na(d$os_time) & !is.na(d$os_event) & d$os_time > 0),
         " (", sum(d$os_event == 1, na.rm = TRUE), " deaths)", logfile = LOG)

if (sum(!is.na(d$stage)) > 0) {
  stg <- as.data.frame(table(d$stage, useNA = "ifany"))
  names(stg) <- c("stage", "n")
  cat(capture.output(print(stg, row.names = FALSE)), sep = "\n",
      file = LOG, append = TRUE)
  print(stg)
}

# ==============================================================================
# 2. Tests
# ==============================================================================
hdr("2. Association tests (D15)")

results <- list()

for (g in c(MAIN_GENE, SUPP_GENE)) {
  ycol <- paste0(g, "_log2tpm")
  log_line("\n  ", g, ":", logfile = LOG)
  rows <- list()

  # --- A1 tumour site ---
  dd <- d[!is.na(d$sample_type) & !is.na(d[[ycol]]), ]
  w <- suppressWarnings(wilcox.test(dd[[ycol]] ~ dd$sample_type, exact = FALSE))
  rows[[1]] <- tibble::tibble(
    gene = g, test = "A1 tumour site (primary vs metastatic)",
    n = nrow(dd), statistic = unname(w$statistic),
    estimate = median(dd[[ycol]][dd$sample_type == "Metastatic"], na.rm = TRUE) -
               median(dd[[ycol]][dd$sample_type == "Primary"], na.rm = TRUE),
    p_value = w$p.value)
  log_line(sprintf("    A1 site     n=%d  median diff (met - primary) = %+.3f  p = %s",
                   nrow(dd), rows[[1]]$estimate, fmt_p(w$p.value)), logfile = LOG)

  # --- A2 stage ---
  ds <- d[!is.na(d$stage) & !is.na(d[[ycol]]), ]
  if (nlevels(droplevels(ds$stage)) > 1) {
    kw <- kruskal.test(ds[[ycol]] ~ droplevels(ds$stage))
    trend <- suppressWarnings(cor.test(as.integer(droplevels(ds$stage)),
                                       ds[[ycol]], method = "spearman"))
    rows[[2]] <- tibble::tibble(
      gene = g, test = "A2 AJCC stage (Kruskal-Wallis)", n = nrow(ds),
      statistic = unname(kw$statistic), estimate = unname(trend$estimate),
      p_value = kw$p.value)
    log_line(sprintf("    A2 stage    n=%d  KW p = %s   | monotone trend rho = %+.3f, p = %s",
                     nrow(ds), fmt_p(kw$p.value), trend$estimate,
                     fmt_p(trend$p.value)), logfile = LOG)
  }

  # --- A3 Breslow ---
  db <- d[!is.na(d$breslow) & !is.na(d[[ycol]]), ]
  if (nrow(db) > 20) {
    ct <- suppressWarnings(cor.test(db$breslow, db[[ycol]], method = "spearman"))
    rows[[3]] <- tibble::tibble(
      gene = g, test = "A3 Breslow thickness (Spearman)", n = nrow(db),
      statistic = NA_real_, estimate = unname(ct$estimate), p_value = ct$p.value)
    log_line(sprintf("    A3 Breslow  n=%d  rho = %+.3f  p = %s",
                     nrow(db), ct$estimate, fmt_p(ct$p.value)), logfile = LOG)
  } else {
    log_line("    A3 Breslow  insufficient data (n = ", nrow(db), ")",
             logfile = LOG)
  }

  # --- A4 survival ---
  dv <- d[!is.na(d$os_time) & d$os_time > 0 & !is.na(d$os_event) &
            !is.na(d[[ycol]]), ]
  if (nrow(dv) > 30 && sum(dv$os_event) > 10) {
    dv$expr <- dv[[ycol]]
    cx <- coxph(Surv(os_time, os_event) ~ expr, data = dv)
    sc <- summary(cx)
    cxa <- coxph(Surv(os_time, os_event) ~ expr + get(PURITY) + sample_type,
                 data = dv)
    sca <- summary(cxa)

    rows[[4]] <- tibble::tibble(
      gene = g, test = "A4 overall survival (Cox, per log2 unit)", n = nrow(dv),
      statistic = unname(sc$coefficients[1, "z"]),
      estimate = unname(sc$coefficients[1, "exp(coef)"]),
      p_value = unname(sc$coefficients[1, "Pr(>|z|)"]))
    log_line(sprintf("    A4 survival n=%d (%d deaths)  HR = %.3f [%.3f, %.3f]  p = %s",
                     nrow(dv), sum(dv$os_event),
                     sc$conf.int[1, "exp(coef)"], sc$conf.int[1, "lower .95"],
                     sc$conf.int[1, "upper .95"],
                     fmt_p(sc$coefficients[1, "Pr(>|z|)"])), logfile = LOG)
    log_line(sprintf("       adjusted (purity + site): HR = %.3f [%.3f, %.3f]  p = %s",
                     sca$conf.int[1, "exp(coef)"], sca$conf.int[1, "lower .95"],
                     sca$conf.int[1, "upper .95"],
                     fmt_p(sca$coefficients[1, "Pr(>|z|)"])), logfile = LOG)
  } else {
    log_line("    A4 survival insufficient data", logfile = LOG)
  }

  r <- dplyr::bind_rows(rows)
  r$p_adj <- p.adjust(r$p_value, method = "BH")   # BH within gene, across A1-A4
  results[[g]] <- r

  log_line(sprintf("    -> %d of %d tests significant after BH correction",
                   sum(r$p_adj < 0.05, na.rm = TRUE), nrow(r)), logfile = LOG)
}
res_tbl <- dplyr::bind_rows(results)

# ==============================================================================
# 3. Verdict
# ==============================================================================
hdr("3. Verdict")

for (g in c(MAIN_GENE, SUPP_GENE)) {
  r <- results[[g]]
  sig <- r[which(r$p_adj < 0.05), ]
  if (nrow(sig) == 0) {
    log_line("  ", g, ": NO aggressiveness measure associated with expression ",
             "after correction. Report as a null result.", logfile = LOG)
  } else {
    log_line("  ", g, ": associated with ", nrow(sig), " measure(s):",
             logfile = LOG)
    for (k in seq_len(nrow(sig))) {
      log_line("    - ", sig$test[k], " (p.adj = ", fmt_p(sig$p_adj[k]), ")",
               logfile = LOG)
    }
  }
}

log_line("\n  Reminder: tumour site is NOT a clean aggressiveness measure -- a ",
         logfile = LOG)
log_line("  metastatic sample indicates the disease metastasised, but the ",
         logfile = LOG)
log_line("  expression difference may reflect the tissue the lesion sits in ",
         logfile = LOG)
log_line("  rather than tumour behaviour. Stage, Breslow and survival are the ",
         logfile = LOG)
log_line("  measures that bear on aggressiveness directly.", logfile = LOG)

# ==============================================================================
# 4. Figures
# ==============================================================================
hdr("4. Figures")

purge_fig_prefix("figure5")
purge_fig_prefix("figureS5")

panel_site <- function(g) {
  ycol <- paste0(g, "_log2tpm")
  r <- results[[g]] |> dplyr::filter(grepl("^A1", test))
  pd <- d[!is.na(d$sample_type) & !is.na(d[[ycol]]), ]
  ymax <- max(pd[[ycol]], na.rm = TRUE)
  ggplot(pd, aes(x = sample_type, y = .data[[ycol]], fill = sample_type)) +
    geom_violin(trim = TRUE, linewidth = 0.3, colour = NA, alpha = 0.3,
                scale = "width") +
    geom_jitter(width = 0.11, size = 0.45, alpha = 0.35, colour = PAL$grey_pt,
                stroke = 0) +
    geom_boxplot(width = 0.2, outlier.shape = NA, linewidth = 0.45,
                 colour = PAL$ink, fill = "white") +
    annotate("segment", x = 1, xend = 2, y = ymax * 1.05, yend = ymax * 1.05,
             linewidth = 0.45, colour = PAL$ink) +
    annotate("text", x = 1.5, y = ymax * 1.08, label = fmt_p_plot(r$p_value),
             parse = TRUE, size = 2.9, vjust = 0, colour = PAL$ink) +
    scale_fill_manual(values = c("Primary" = PAL$uv_low,
                                 "Metastatic" = PAL$uv_high), guide = "none") +
    coord_cartesian(ylim = c(0, ymax * 1.2)) +
    labs(x = NULL, y = bquote(bold(italic(.(g))~log[2]*"(TPM+1)")),
         subtitle = "Tumour site") +
    theme_pub()
}

panel_stage <- function(g) {
  ycol <- paste0(g, "_log2tpm")
  pd <- d[!is.na(d$stage) & !is.na(d[[ycol]]), ]
  pd$stage <- droplevels(pd$stage)
  r <- results[[g]] |> dplyr::filter(grepl("^A2", test))
  ggplot(pd, aes(x = stage, y = .data[[ycol]], fill = stage)) +
    geom_boxplot(outlier.size = 0.4, outlier.colour = PAL$grey_pt,
                 linewidth = 0.4, width = 0.62, colour = PAL$ink) +
    scale_fill_brewer(palette = "YlOrRd", guide = "none") +
    labs(x = "AJCC stage", y = bquote(bold(italic(.(g))~log[2]*"(TPM+1)")),
         subtitle = if (nrow(r)) paste0("Stage (KW p = ",
                                        formatC(r$p_value, format = "f",
                                                digits = 3), ")") else "Stage") +
    theme_pub(grid = TRUE)
}

panel_breslow <- function(g) {
  ycol <- paste0(g, "_log2tpm")
  pd <- d[!is.na(d$breslow) & !is.na(d[[ycol]]), ]
  r <- results[[g]] |> dplyr::filter(grepl("^A3", test))
  p <- ggplot(pd, aes(x = breslow, y = .data[[ycol]])) +
    geom_point(size = 0.9, alpha = 0.6, colour = PAL$hist, stroke = 0) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = PAL$ink,
                fill = "grey70", linewidth = 0.6, alpha = 0.25) +
    scale_x_log10() +
    labs(x = "Breslow thickness (mm)",
         y = bquote(bold(italic(.(g))~log[2]*"(TPM+1)")),
         subtitle = "Tumour thickness") +
    theme_pub()
  if (nrow(r)) {
    p <- p + annotate("text", x = max(pd$breslow, na.rm = TRUE),
                      y = max(pd[[ycol]], na.rm = TRUE),
                      label = sprintf("rho == %+.3f", r$estimate), parse = TRUE,
                      hjust = 1, vjust = 1, size = 2.9, colour = PAL$ink)
  }
  p
}

panel_km <- function(g) {
  ycol <- paste0(g, "_log2tpm")
  dv <- d[!is.na(d$os_time) & d$os_time > 0 & !is.na(d$os_event) &
            !is.na(d[[ycol]]), ]
  dv$grp <- factor(ifelse(dv[[ycol]] > median(dv[[ycol]], na.rm = TRUE),
                          "High", "Low"), levels = c("Low", "High"))
  sf <- survfit(Surv(os_time / 365.25, os_event) ~ grp, data = dv)
  sd <- survdiff(Surv(os_time, os_event) ~ grp, data = dv)
  pval <- 1 - pchisq(sd$chisq, length(sd$n) - 1)

  km <- tibble::tibble(
    time = sf$time, surv = sf$surv,
    grp = rep(gsub("grp=", "", names(sf$strata)), sf$strata)
  )
  km$grp <- factor(km$grp, levels = c("Low", "High"))

  ggplot(km, aes(x = time, y = surv, colour = grp)) +
    geom_step(linewidth = 0.6) +
    scale_colour_manual(values = c("Low" = PAL$uv_high, "High" = PAL$uv_low),
                        name = NULL) +
    annotate("text", x = max(km$time) * 0.95, y = 0.95,
             label = fmt_p_plot(pval), parse = TRUE, hjust = 1, size = 2.9,
             colour = PAL$ink) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Years", y = "Overall survival",
         subtitle = bquote("Median split,"~italic(.(g)))) +
    theme_pub(grid = TRUE) +
    theme(legend.position = c(0.12, 0.2))
}

build_fig <- function(g) {
  (panel_site(g) | panel_stage(g)) / (panel_breslow(g) | panel_km(g)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = FIG$base_size + 3, face = "bold"))
}

save_fig(build_fig(MAIN_GENE), "figure5_rassf9_aggressiveness",
         width = FIG$width_double, height = 140)
save_fig(build_fig(SUPP_GENE), "figureS5_rassf10_aggressiveness",
         width = FIG$width_double, height = 140)

# ==============================================================================
# 5. Tables
# ==============================================================================
hdr("5. Results tables")

readr::write_csv(res_tbl |> dplyr::filter(gene == MAIN_GENE),
                 file.path(PATHS$results, "table_main_09_rassf9_aggressiveness.csv"))
readr::write_csv(res_tbl |> dplyr::filter(gene == SUPP_GENE),
                 file.path(PATHS$results, "table_supp_09_rassf10_aggressiveness.csv"))

cat(capture.output(print(as.data.frame(res_tbl), row.names = FALSE, digits = 3)),
    sep = "\n", file = LOG, append = TRUE)
print(as.data.frame(res_tbl), digits = 3)

hdr("6. Warnings")
report_warnings(LOG)

write_session_info("09_aggressiveness")

hdr("DONE")
log_line("Figure 5 (main) : figures/figure5_rassf9_aggressiveness.*", logfile = LOG)
log_line("Figure S5 (supp): figures/figureS5_rassf10_aggressiveness.*", logfile = LOG)
log_line("\nWhatever this shows, it is now an analysis that exists and can be ",
         "described accurately to Gustavo.", logfile = LOG)
log_line("Full log: ", LOG, logfile = LOG)
