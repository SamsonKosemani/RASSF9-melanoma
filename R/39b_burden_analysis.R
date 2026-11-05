# ==============================================================================
# 39b_burden_analysis.R
#
# Mutational burden against RASSF9 and RASSF10 expression.
#
# THE REQUEST, and how it is handled. Gustavo has asked whether tumours with
# higher mutational burden show lower RASSF9, which would support the proposal
# that burden drives progression by repressing RASSF9. He has also asked for all
# negative results. A request that names a desired direction is exactly the
# situation in which an analysis drifts -- try total burden, then UV burden,
# then the fraction, then a tertile split, and one of them clears p < 0.05 by
# chance. So the primary test and the whole family are fixed below, before any
# number is computed, the family is corrected for multiplicity, and every
# variant is reported whatever it shows.
#
# THE CONFOUND THAT MATTERS MOST, and it is one this project created today.
# Three platforms now agree that RASSF9 in bulk tumour tissue is largely NOT
# melanocyte-derived, and in this same cohort purity depresses RASSF9
# expression (rho = -0.204, p = 3.05e-05). Mutational burden estimates are
# themselves purity-dependent, because variant allele frequencies fall as
# normal tissue dilutes the sample. A negative burden-expression correlation is
# therefore precisely what composition alone would produce, with no repression
# anywhere in the mechanism. Section 4 tests the premises of that explanation
# individually, as was done for the methylation result, before any association
# is interpreted.
#
# WHAT CANNOT BE CLAIMED. These are cross-sectional bulk measurements of
# resected tumours. They can show whether an association exists and whether it
# survives adjustment. They cannot establish that burden causes repression,
# that repression causes progression, or the direction of any arrow. The script
# reports association only, and says so in its own output.
#
# PRE-SPECIFIED, before the numbers:
#   PRIMARY    Spearman correlation of RASSF9 expression with TOTAL burden
#              (TMB_sig) across all 411 samples.
#   FAMILY     4 burden definitions (total, UV/SBS7, non-UV = total - SBS7,
#              UV fraction) x 2 genes = 8 tests, Benjamini-Hochberg corrected
#              across the family.
#   CRITERION  ASSOCIATION PRESENT  = BH p < 0.05 for the primary test AND the
#                                     adjusted model keeps sign at p < 0.05
#              PURITY-EXPLAINED     = unadjusted significant, adjusted effect
#                                     loses more than 60% or significance, AND
#                                     purity relates to both burden and
#                                     expression in the directions required
#              ABSENT               = no significant association
#   ROBUSTNESS Spearman beside Pearson on log10 burden, plus recomputation with
#              the top-decile hypermutators removed. The original 317-sample
#              analysis gave rho = -0.130 (p = 0.021) against r = -0.033
#              (p = 0.564); that gap usually means a few extreme tumours carry
#              the rank statistic, and it should be diagnosed, not quoted.
#
# Output: results/tcga/burden_family.csv
#         results/tcga/burden_adjusted.csv
#         results/tcga/burden_strata.csv
#         figures/tcga/tcga_burden_rassf9.{pdf,png,svg}
#         results/logs/39b_burden_analysis.log
# ==============================================================================

source(here::here("R", "00_setup.R"))

LOG <- file.path(PATHS$logs, "39b_burden_analysis.log")
cat("Mutational burden vs RASSF9/RASSF10 —", format(Sys.time()), "\n\n", file = LOG)
hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

PROC <- PATHS$processed
RES  <- file.path(PATHS$results, "tcga")
FIG_TCGA <- file.path(PATHS$figures, "tcga")
invisible(lapply(c(RES, FIG_TCGA), dir.create, recursive = TRUE, showWarnings = FALSE))

GENES    <- c("RASSF9", "RASSF10")
RETAIN   <- 0.60
TOP_FRAC <- 0.10

pick <- function(slot, fb) {
  if (!exists("PAL", inherits = TRUE)) return(fb)
  x <- get("PAL", inherits = TRUE)[[slot]]
  if (!is.null(x) && is.character(x) && nzchar(x[1])) x[1] else fb
}
COL_A <- pick("nevus", "#0095D9"); COL_B <- pick("melanoma", "#F2680C")
COL_GREY <- pick("grey_pt", "#8E8E8E"); COL_INK <- pick("ink", "#1A1A1A")

# ==============================================================================
# 1. Data
# ==============================================================================
hdr("1. Data")

f <- file.path(PROC, "full_cohort.csv")
if (!file.exists(f)) stop("Missing ", f)
D <- readr::read_csv(f, show_col_types = FALSE, progress = FALSE)
log_line("  full_cohort.csv: ", nrow(D), " samples", logfile = LOG)

need <- c("TMB_sig", "SBS7_total", "SBS7_fraction", "CPE", "sample_type",
          "RASSF9_log2tpm", "RASSF10_log2tpm")
miss <- setdiff(need, names(D))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

# Non-UV burden separates a UV-specific effect from a general mutation load.
D$nonUV <- D$TMB_sig - D$SBS7_total
neg <- sum(D$nonUV < 0, na.rm = TRUE)
if (neg) {
  log_line("  [FLAG] ", neg, " samples have SBS7_total > TMB_sig, which should",
           logfile = LOG)
  log_line("  not happen; those are set to NA rather than clamped to zero.",
           logfile = LOG)
  D$nonUV[D$nonUV < 0] <- NA_real_
}

BURDENS <- list(
  "Total burden"   = "TMB_sig",
  "UV (SBS7)"      = "SBS7_total",
  "Non-UV burden"  = "nonUV",
  "UV fraction"    = "SBS7_fraction")

for (nm in names(BURDENS)) {
  v <- D[[BURDENS[[nm]]]]
  log_line(sprintf("    %-15s median %10.1f | IQR %9.1f-%-9.1f | max %9.1f | NA %d",
                   nm, stats::median(v, na.rm = TRUE),
                   stats::quantile(v, 0.25, na.rm = TRUE),
                   stats::quantile(v, 0.75, na.rm = TRUE),
                   max(v, na.rm = TRUE), sum(is.na(v))), logfile = LOG)
}
log_line("  purity (CPE): median ", sprintf("%.3f", stats::median(D$CPE, na.rm = TRUE)),
         ", missing ", sum(is.na(D$CPE)), logfile = LOG)
log_line("  sample types: ",
         paste(names(table(D$sample_type)), as.integer(table(D$sample_type)),
               sep = "=", collapse = ", "), logfile = LOG)
log_line("  NOTE: 411 is the ceiling. Burden exists only for samples present in",
         logfile = LOG)
log_line("  the SigProfiler table, so the 469-sample expression/methylation",
         logfile = LOG)
log_line("  rebuild does not extend this analysis.", logfile = LOG)

# ==============================================================================
# 2. The pre-specified family
# ==============================================================================
hdr("2. The family of tests, all eight reported")

lg <- function(x) log10(x + 1)

rows <- list()
for (g in GENES) {
  y <- D[[paste0(g, "_log2tpm")]]
  for (nm in names(BURDENS)) {
    x <- D[[BURDENS[[nm]]]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 30) next
    sp <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman"))
    xl <- if (nm == "UV fraction") x[ok] else lg(x[ok])
    pe <- suppressWarnings(stats::cor.test(xl, y[ok], method = "pearson"))
    rows[[length(rows) + 1L]] <- tibble::tibble(
      gene = g, burden = nm, n = sum(ok),
      rho = unname(sp$estimate), rho_p = sp$p.value,
      r = unname(pe$estimate), r_p = pe$p.value)
  }
}
FAM <- dplyr::bind_rows(rows)
FAM$rho_bh <- stats::p.adjust(FAM$rho_p, method = "BH")
FAM$r_bh   <- stats::p.adjust(FAM$r_p,   method = "BH")

log_line(sprintf("  %-8s %-15s %5s %8s %11s %11s %9s %11s %11s", "gene",
                 "burden", "n", "rho", "p", "BH p", "r", "p", "BH p"),
         logfile = LOG)
for (i in seq_len(nrow(FAM)))
  log_line(sprintf("  %-8s %-15s %5d %+8.3f %11s %11s %+9.3f %11s %11s",
                   FAM$gene[i], FAM$burden[i], FAM$n[i], FAM$rho[i],
                   fmt_p(FAM$rho_p[i]), fmt_p(FAM$rho_bh[i]), FAM$r[i],
                   fmt_p(FAM$r_p[i]), fmt_p(FAM$r_bh[i])), logfile = LOG)
readr::write_csv(FAM, file.path(RES, "burden_family.csv"))

log_line("", logfile = LOG)
log_line("  Spearman uses raw burden; Pearson uses log10(burden + 1), except",
         logfile = LOG)
log_line("  the UV fraction which is already bounded. BH correction is across",
         logfile = LOG)
log_line("  all eight tests, as pre-specified.", logfile = LOG)

# ==============================================================================
# 3. Robustness: are hypermutators carrying the rank statistic?
# ==============================================================================
hdr("3. Robustness to extreme burden")

log_line("  Burden spans 8 to 26,500 mutations, so a handful of hypermutated",
         logfile = LOG)
log_line("  tumours can dominate a rank correlation while contributing little",
         logfile = LOG)
log_line("  to a linear one. Each test is recomputed with the top ",
         round(100 * TOP_FRAC), "% of", logfile = LOG)
log_line("  burden removed.", logfile = LOG)

rb <- list()
for (g in GENES) {
  y <- D[[paste0(g, "_log2tpm")]]
  for (nm in names(BURDENS)) {
    x <- D[[BURDENS[[nm]]]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 40) next
    cut <- stats::quantile(x[ok], 1 - TOP_FRAC, na.rm = TRUE)
    keep <- ok & x <= cut
    s0 <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman"))
    s1 <- suppressWarnings(stats::cor.test(x[keep], y[keep], method = "spearman"))
    log_line(sprintf("  %-8s %-15s full n = %3d rho = %+.3f (%s) | trimmed n = %3d rho = %+.3f (%s)",
                     g, nm, sum(ok), s0$estimate, fmt_p(s0$p.value),
                     sum(keep), s1$estimate, fmt_p(s1$p.value)), logfile = LOG)
    rb[[length(rb) + 1L]] <- tibble::tibble(
      gene = g, burden = nm, rho_full = unname(s0$estimate), p_full = s0$p.value,
      rho_trim = unname(s1$estimate), p_trim = s1$p.value, n_trim = sum(keep))
  }
}
RB <- dplyr::bind_rows(rb)

# ==============================================================================
# 4. The purity confound, premises first
# ==============================================================================
hdr("4. The purity-composition confound")

log_line("  For purity to manufacture a negative burden-expression correlation,",
         logfile = LOG)
log_line("  all three must hold:", logfile = LOG)
log_line("    (a) purity NEGATIVELY correlated with expression", logfile = LOG)
log_line("    (b) purity POSITIVELY correlated with burden", logfile = LOG)
log_line("    (c) and the burden-expression association must then weaken when",
         logfile = LOG)
log_line("        purity is held constant", logfile = LOG)
log_line("", logfile = LOG)

prem <- list()
sc <- function(a, b, lbl) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 30) return(NULL)
  ct <- suppressWarnings(stats::cor.test(a[ok], b[ok], method = "spearman"))
  log_line(sprintf("    %-48s n = %3d  rho = %+.3f  p = %s", lbl, sum(ok),
                   ct$estimate, fmt_p(ct$p.value)), logfile = LOG)
  tibble::tibble(comparison = lbl, n = sum(ok), rho = unname(ct$estimate),
                 p = ct$p.value)
}
for (g in GENES)
  prem[[length(prem) + 1L]] <- sc(D$CPE, D[[paste0(g, "_log2tpm")]],
    paste0("(a) purity vs ", g, " expression"))
for (nm in names(BURDENS))
  prem[[length(prem) + 1L]] <- sc(D$CPE, D[[BURDENS[[nm]]]],
    paste0("(b) purity vs ", nm))
PREM <- dplyr::bind_rows(prem)

# ==============================================================================
# 5. Adjustment
# ==============================================================================
hdr("5. Adjusted models")

fit1 <- function(g, bn, extra, lbl) {
  x <- D[[BURDENS[[bn]]]]
  xb <- if (bn == "UV fraction") x else lg(x)
  dd <- data.frame(y = D[[paste0(g, "_log2tpm")]], b = xb)
  for (nm in names(extra)) dd[[nm]] <- extra[[nm]]
  dd <- dd[stats::complete.cases(dd), , drop = FALSE]
  if (nrow(dd) < 40) return(NULL)
  keep <- vapply(setdiff(names(dd), c("y", "b")),
                 function(cn) dplyr::n_distinct(dd[[cn]]) > 1, logical(1))
  use <- c("b", names(keep)[keep])
  fit <- try(stats::lm(stats::as.formula(paste("y ~", paste(use, collapse = " + "))),
                       data = dd), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  cf <- summary(fit)$coefficients
  if (!"b" %in% rownames(cf)) return(NULL)
  ci <- try(stats::confint(fit, "b"), silent = TRUE)
  tibble::tibble(gene = g, burden = bn, model = lbl, n = stats::nobs(fit),
                 beta = cf["b", 1],
                 lo = if (inherits(ci, "try-error")) NA_real_ else ci[1],
                 hi = if (inherits(ci, "try-error")) NA_real_ else ci[2],
                 p = cf["b", 4])
}

ad <- list()
for (g in GENES) for (bn in names(BURDENS)) {
  mods <- list(
    list("M0 unadjusted", list()),
    list("M1 + sample type", list(sample_type = D$sample_type)),
    list("M2 + purity", list(CPE = D$CPE)),
    list("M3 + purity + sample type",
         list(CPE = D$CPE, sample_type = D$sample_type)))
  for (m in mods) {
    r <- fit1(g, bn, m[[2]], m[[1]])
    if (!is.null(r)) ad[[length(ad) + 1L]] <- r
  }
}
AD <- dplyr::bind_rows(ad)
readr::write_csv(AD, file.path(RES, "burden_adjusted.csv"))

for (g in GENES) for (bn in names(BURDENS)) {
  s <- AD |> dplyr::filter(gene == g, burden == bn)
  if (!nrow(s)) next
  b0 <- s$beta[s$model == "M0 unadjusted"]
  log_line("\n  ", g, " / ", bn, ":", logfile = LOG)
  for (i in seq_len(nrow(s))) {
    ret <- if (length(b0) && is.finite(b0) && abs(b0) >= 0.02)
      sprintf("%4.0f%% retained", 100 * s$beta[i] / b0) else "  baseline ~ 0"
    log_line(sprintf("    %-28s n = %3d  beta = %+.4f [%+.4f, %+.4f]  p = %-10s %s",
                     s$model[i], s$n[i], s$beta[i], s$lo[i], s$hi[i],
                     fmt_p(s$p[i]), ret), logfile = LOG)
  }
}

# ==============================================================================
# 6. Strata
# ==============================================================================
hdr("6. Strata")

st <- list()
strat_run <- function(sub, lbl) {
  for (g in GENES) {
    y <- sub[[paste0(g, "_log2tpm")]]
    for (nm in c("Total burden", "UV (SBS7)")) {
      x <- sub[[BURDENS[[nm]]]]
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 30) next
      ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman"))
      log_line(sprintf("    %-8s %-15s %-26s n = %3d  rho = %+.3f  p = %s",
                       g, nm, lbl, sum(ok), ct$estimate, fmt_p(ct$p.value)),
               logfile = LOG)
      st[[length(st) + 1L]] <<- tibble::tibble(
        gene = g, burden = nm, stratum = lbl, n = sum(ok),
        rho = unname(ct$estimate), p = ct$p.value)
    }
  }
}
for (k in sort(unique(D$sample_type))) {
  s <- D[D$sample_type == k, ]
  if (nrow(s) >= 30) strat_run(s, paste0("sample type: ", k))
}
if ("has_uv_signature" %in% names(D)) {
  strat_run(D[which(D$has_uv_signature), ],  "UV signature present")
  strat_run(D[which(!D$has_uv_signature), ], "no UV signature")
}
ST <- dplyr::bind_rows(st)
if (nrow(ST)) readr::write_csv(ST, file.path(RES, "burden_strata.csv"))

# ==============================================================================
# 7. Verdict
# ==============================================================================
hdr("7. Verdict against the pre-set criterion")

prim <- FAM |> dplyr::filter(gene == "RASSF9", burden == "Total burden")
for (g in GENES) {
  p0 <- FAM |> dplyr::filter(gene == g, burden == "Total burden")
  a3 <- AD |> dplyr::filter(gene == g, burden == "Total burden",
                            model == "M3 + purity + sample type")
  a0 <- AD |> dplyr::filter(gene == g, burden == "Total burden",
                            model == "M0 unadjusted")
  if (!nrow(p0) || !nrow(a3) || !nrow(a0)) next
  ret <- if (abs(a0$beta[1]) >= 0.02) a3$beta[1] / a0$beta[1] else NA_real_
  sig_un <- p0$rho_bh[1] < 0.05
  sig_adj <- a3$p[1] < 0.05 && sign(a3$beta[1]) == sign(a0$beta[1])

  verdict <- if (sig_un && sig_adj) "ASSOCIATION PRESENT"
  else if (sig_un && !sig_adj) "PURITY-EXPLAINED or lost on adjustment"
  else "ABSENT"
  log_line("  ", g, " vs total burden: ", verdict, logfile = LOG)
  log_line(sprintf("    rho = %+.3f, BH p = %s | adjusted beta %+.4f (p = %s)%s",
                   p0$rho[1], fmt_p(p0$rho_bh[1]), a3$beta[1], fmt_p(a3$p[1]),
                   if (is.finite(ret)) sprintf(", %.0f%% retained", 100 * ret) else ""),
           logfile = LOG)
}
log_line("", logfile = LOG)
log_line("  WHAT THIS CANNOT SAY, and it should be stated in any write-up:",
         logfile = LOG)
log_line("  these are cross-sectional bulk measurements of resected tumours.",
         logfile = LOG)
log_line("  An association, if present, does not establish that burden causes",
         logfile = LOG)
log_line("  repression of RASSF9, that repression drives progression, or the",
         logfile = LOG)
log_line("  direction of any arrow. A tumour that is more mutated and has less",
         logfile = LOG)
log_line("  RASSF9 is consistent with repression, with selection, with", logfile = LOG)
log_line("  differing cell composition, and with reverse causation equally.",
         logfile = LOG)
log_line("", logfile = LOG)
log_line("  FURTHER: RASSF9 in bulk tumour is largely non-melanocytic across",
         logfile = LOG)
log_line("  three platforms, so any burden association measured in bulk may",
         logfile = LOG)
log_line("  reflect stromal and vascular content rather than melanocyte",
         logfile = LOG)
log_line("  biology, whatever the purity adjustment shows.", logfile = LOG)

# ==============================================================================
# 8. Figure
# ==============================================================================
hdr("8. Figure")

reset_legends("# TCGA mutational burden figure legends\n")

pa <- ggplot(D, aes(x = lg(TMB_sig), y = RASSF9_log2tpm)) +
  geom_point(size = 0.8, alpha = 0.5, colour = COL_GREY) +
  geom_smooth(method = "loess", formula = y ~ x, se = TRUE,
              colour = COL_B, fill = COL_B, alpha = 0.15, linewidth = 0.6) +
  labs(x = expression(bold(log[10]*"(total burden + 1)")),
       y = expression(bold("RASSF9  "*log[2]*"(TPM+1)"))) +
  theme_journal() +
  theme(plot.margin = margin(4, 12, 4, 4, unit = "pt"))

pb <- ggplot(FAM, aes(x = rho, y = interaction(burden, gene, sep = " — "),
                      colour = gene)) +
  geom_vline(xintercept = 0, linetype = "22", colour = COL_GREY, linewidth = 0.4) +
  geom_point(aes(shape = rho_bh < 0.05), size = 2) +
  scale_colour_manual(values = c(RASSF9 = COL_B, RASSF10 = COL_A), name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "BH p < 0.05", `FALSE` = "ns"),
                     name = NULL) +
  labs(x = "Spearman rho", y = NULL) +
  theme_journal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = FIG$base_size - 3),
        legend.key.size = unit(8, "pt"),
        plot.margin = margin(4, 12, 4, 4, unit = "pt"))

fig <- patchwork::wrap_plots(list(pa, pb), nrow = 1) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = FIG$base_size + 2, face = "bold",
                                colour = COL_INK))
save_fig_to(fig, "tcga_burden_rassf9", FIG_TCGA,
            width = FIG$width_double, height = 88)
log_line("  wrote: tcga_burden_rassf9", logfile = LOG)

add_legend("tcga_burden_rassf9",
  "Figure. Mutational burden and RASSF9 expression in TCGA-SKCM.",
  paste0("(A) RASSF9 expression against total mutational burden across ",
         nrow(D), " cutaneous melanomas, with a loess fit and 95% confidence ",
         "band; burden is plotted on a log scale because it spans three orders ",
         "of magnitude. (B) Spearman correlations for all eight pre-specified ",
         "tests — four burden definitions (total, UV-associated SBS7, non-UV ",
         "remainder, and UV fraction) for each gene — with filled points ",
         "significant after Benjamini-Hochberg correction across the family. ",
         "The non-UV term separates a UV-specific effect from a general ",
         "mutation-load effect. Associations are reported without causal ",
         "interpretation; these are cross-sectional bulk measurements."))

report_warnings(logfile = LOG)
write_session_info("39b_burden_analysis")
log_line("\nFull log: ", LOG, logfile = NULL)
