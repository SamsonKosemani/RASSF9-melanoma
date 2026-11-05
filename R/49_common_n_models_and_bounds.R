## =============================================================================
## 49_common_n_models_and_bounds.R
##
## Addresses the three points raised by Prof. Amarante-Mendes (26 Oct 2026).
##
## POINT 1 -- the p-value floor belongs to the normal approximation, not the
##            exact test
##   Script 48 stated that no Mann-Whitney p below about 1e-36 is attainable
##   for n1 = 56 against n2 = 940. That bound is a property of the NORMAL
##   APPROXIMATION: |z| cannot exceed (n1*n2/2 - 0.5)/sigma_U = 12.59, and the
##   normal tail at that z is ~1e-36.
##
##   The EXACT Mann-Whitney has a far smaller floor. Under complete separation
##   the exact two-sided p is 2 / choose(N, n1) = 2 / choose(996, 56), which is
##   of order 1e-92. Both are computed below rather than asserted.
##
##   This changes no reported value. The two significant comparisons give
##   p = 8.7e-17 and 1.1e-16 under the normal approximation, far above either
##   floor. It changes only what the Methods may claim.
##
## POINT 2a -- compare nested models on the same lines
##   Script 46 fitted M0 and M1 on 736 lines and M2 and M3 on 611, because copy
##   number is missing for 125. Percentages retained across models fitted on
##   different samples confound adjustment with sample change. All four models
##   are refitted here on the common complete-case set.
##
## POINT 2b -- subtypes among the 39 lines with paired data
##   Reported below, with expression and methylation summarised per subtype.
##
## POINT 3 -- wording
##   No analysis change. Handled in the note: methylation-associated repression
##   rather than demonstrated silencing, and the open question of whether UV
##   inducibility is affected is stated explicitly.
##
## INPUTS   results/ccle/46_line_level_data.csv
## OUTPUTS  results/ccle/49_models_common_n_rassf9.csv  /  ..._rassf10.csv
##          results/ccle/49_subtypes_paired_lines.csv
##          results/ccle/49_pvalue_bounds.csv
##          results/logs/49_common_n_models_and_bounds.log
## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------
if (file.exists("R/00_setup.R")) source("R/00_setup.R")

suppressPackageStartupMessages({
  library(data.table)
})

ROOT    <- if (exists("PROJECT_ROOT")) PROJECT_ROOT else getwd()
RES_DIR <- file.path(ROOT, "results", "ccle")
LOG_DIR <- file.path(ROOT, "results", "logs")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_PATH <- file.path(LOG_DIR, "49_common_n_models_and_bounds.log")
.logcon  <- file(LOG_PATH, open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = .logcon) }
WARN <- character(0)
warn <- function(...) { m <- paste0(...); WARN <<- c(WARN, m); say("  [WARN] ", m) }
on.exit({
  if (length(WARN)) { say(""); say("WARNINGS RAISED (", length(WARN), "):"); for (w in WARN) say("  - ", w) }
  else say("  no warnings raised")
  say(""); say("Full log: ", LOG_PATH); close(.logcon)
}, add = TRUE)

set.seed(1234)
GENES <- c("RASSF9", "RASSF10")

f_ll <- file.path(RES_DIR, "46_line_level_data.csv")
if (!file.exists(f_ll)) stop("missing ", basename(f_ll), " -- run 46 first", call. = FALSE)
D <- fread(f_ll, showProgress = FALSE)
say("")
say("== 0. Input ====================================================")
say("  ", basename(f_ll), ": ", nrow(D), " lines")

## ---- 1. The two floors ------------------------------------------------------
say("")
say("== 1. Mann-Whitney p-value floors ==============================")
say("  Both are computed, not asserted. n1 = 56 melanoma, n2 = 940 others.")
say("")

n1 <- 56; n2 <- 940; N <- n1 + n2
sigma <- sqrt(n1 * n2 * (N + 1) / 12)
z_max <- (n1 * n2 / 2 - 0.5) / sigma
log10_norm <- (log(2) + stats::pnorm(-z_max, log.p = TRUE)) / log(10)
log10_exact <- (log(2) - lchoose(N, n1)) / log(10)

say(sprintf("  NORMAL APPROXIMATION"))
say(sprintf("    sigma_U = %.1f | max |z| = %.3f", sigma, z_max))
say(sprintf("    floor on the two-sided p  ~ 1e%.1f", log10_norm))
say(sprintf("    This is the bound that applies to the tests actually reported,"))
say(sprintf("    because those use the normal approximation with tie correction."))
say("")
say(sprintf("  EXACT TEST"))
say(sprintf("    number of distinct group assignments: choose(%d, %d) = 1e%.1f",
            N, n1, lchoose(N, n1) / log(10)))
say(sprintf("    floor on the two-sided p  = 2 / choose(N, n1) ~ 1e%.1f", log10_exact))
say(sprintf("    The exact floor is about 1e%.0f times smaller than the normal one.",
            abs(log10_exact) - abs(log10_norm)))
say("")
say("  WHAT THIS CORRECTS: script 48 described ~1e-36 as the smallest")
say("  attainable p for this comparison. That holds for the normal")
say("  approximation only. The exact test can in principle go far lower. No")
say("  reported value changes: the two significant comparisons are 8.7e-17 and")
say("  1.1e-16, far above both floors, and the exact test is not used here")
say("  because with n = 996 it is computationally impractical and the normal")
say("  approximation is the standard choice at this sample size.")

BND <- data.table(
  method = c("Normal approximation (used)", "Exact Mann-Whitney (not used)"),
  n1 = n1, n2 = n2,
  statistic = c(sprintf("max |z| = %.3f", z_max),
                sprintf("choose(%d,%d) = 1e%.1f", N, n1, lchoose(N, n1) / log(10))),
  log10_floor = c(log10_norm, log10_exact))
fwrite(BND, file.path(RES_DIR, "49_pvalue_bounds.csv"))

## ---- 2. Subtypes among the paired melanoma lines ----------------------------
say("")
say("== 2. Melanoma subtypes among the lines with paired data =======")

D[, is_mel := !is.na(disease) & grepl("melanoma", disease, ignore.case = TRUE)]
D[, is_uveal := !is.na(subtype) &
     grepl("uveal|ocular", paste(subtype, disease), ignore.case = TRUE)]

MEL <- D[is_mel == TRUE & is_uveal == FALSE]
PAIR <- MEL[is.finite(RASSF9_core) & is.finite(RASSF9)]

say("  melanoma lines (non-uveal) in the matched file: ", nrow(MEL))
say("  with both methylation and expression          : ", nrow(PAIR))

ST <- PAIR[, .(n = .N,
               pct_det_R9  = 100 * mean(RASSF9  >= log2(2)),
               median_R9   = median(RASSF9),
               max_R9      = max(RASSF9),
               median_core = median(RASSF9_core),
               min_core    = min(RASSF9_core),
               max_core    = max(RASSF9_core)),
           by = subtype][order(-n)]

say("")
say(sprintf("  %-38s %4s %8s %9s %8s %11s %13s",
            "subtype", "n", "%TPM>1", "medianExp", "maxExp", "medianCore", "core range"))
for (i in seq_len(nrow(ST))) {
  s <- ST[i]
  say(sprintf("  %-38s %4d %7.1f%% %9.3f %8.3f %11.3f %6.3f-%.3f",
              s$subtype, s$n, s$pct_det_R9, s$median_R9, s$max_R9,
              s$median_core, s$min_core, s$max_core))
}
say("")
say("  The subtype labels are DepMap Oncotree values. Lines labelled simply")
say("  'Melanoma' are of unstated site and are NOT assumed to be cutaneous.")
if (!any(grepl("acral", ST$subtype, ignore.case = TRUE)))
  say("  No acral line appears here: none is present in the methylation panel.")

fwrite(ST, file.path(RES_DIR, "49_subtypes_paired_lines.csv"))

## ---- 3. Nested models on a common set ---------------------------------------
say("")
say("== 3. Nested models refitted on identical lines =================")
say("  Script 46 fitted M0 and M1 on 736 lines and M2 and M3 on 611, because")
say("  copy number is missing for some. A percentage retained computed across")
say("  models fitted on different samples confounds adjustment with sample")
say("  change. All four models are refitted here on the complete-case set.")

for (g in GENES) {
  cn   <- paste0(g, "_core")
  cncn <- paste0(g, "_cn")
  if (!all(c(cn, cncn, g) %in% names(D))) { warn("columns missing for ", g); next }

  d <- D[is.finite(get(cn)) & is.finite(get(g)) & !is.na(lineage)]
  d[, lin := factor(lineage)]
  keep <- d[, .N, by = lin][N >= 3, lin]
  d <- d[lin %in% keep]
  d[, lin := droplevels(lin)]

  n_all <- nrow(d)
  dc <- d[is.finite(get(cncn)) & get(cncn) > 0]
  dc[, lin := droplevels(lin)]
  ## a factor level with a single observation cannot be fitted
  keep2 <- dc[, .N, by = lin][N >= 3, lin]
  dc <- dc[lin %in% keep2]
  dc[, lin := droplevels(lin)]

  say("")
  say("  ", g, ": lines with methylation + expression + lineage = ", n_all)
  say("         of those with usable copy number (complete case) = ", nrow(dc))
  if (nrow(dc) < 100) { warn(g, ": too few complete cases"); next }

  mods <- list(
    "M0 unadjusted"                  = as.formula(paste(g, "~", cn)),
    "M1 + lineage"                   = as.formula(paste(g, "~", cn, "+ lin")),
    "M2 + log2 copy number"          = as.formula(paste(g, "~", cn, "+ log2(", cncn, ")")),
    "M3 + lineage + log2 copy number"= as.formula(paste(g, "~", cn, "+ lin + log2(", cncn, ")"))
  )

  say("")
  say(sprintf("    %-32s %5s %10s %22s %11s %9s",
              "model", "n", "beta", "95% CI", "p", "retained"))
  rows <- list(); b0 <- NA_real_
  for (nm in names(mods)) {
    fit <- try(stats::lm(mods[[nm]], data = dc), silent = TRUE)
    if (inherits(fit, "try-error")) { say("    ", nm, ": failed"); next }
    co <- summary(fit)$coefficients
    if (!(cn %in% rownames(co))) { say("    ", nm, ": core term dropped"); next }
    est <- co[cn, 1]; se <- co[cn, 2]; pv <- co[cn, 4]
    ci  <- est + c(-1.96, 1.96) * se
    if (is.na(b0)) b0 <- est
    say(sprintf("    %-32s %5d %+10.3f  [%+.3f, %+.3f] %11.2e %8.0f%%",
                nm, stats::nobs(fit), est, ci[1], ci[2], pv, 100 * est / b0))
    rows[[length(rows) + 1]] <- data.table(gene = g, model = nm, n = stats::nobs(fit),
                                           beta = est, lo = ci[1], hi = ci[2], p = pv,
                                           pct_retained = 100 * est / b0)
  }
  if (length(rows)) {
    RR <- rbindlist(rows)
    fwrite(RR, file.path(RES_DIR, sprintf("49_models_common_n_%s.csv", tolower(g))))
    if (length(unique(RR$n)) > 1)
      warn(g, ": models still differ in n (", paste(unique(RR$n), collapse = ", "),
           ") -- a lineage level was dropped between models")
    else
      say("    all four models fitted on the same ", RR$n[1], " lines")
  }

  ## the comparison script 46 could not make: does restricting to complete
  ## cases change the unadjusted estimate?
  f0 <- stats::lm(as.formula(paste(g, "~", cn)), data = d)
  c0 <- summary(f0)$coefficients[cn, ]
  say(sprintf("    for reference, M0 on all %d lines: beta = %+.3f, p = %.2e",
              n_all, c0[1], c0[4]))
  say("    (a difference between this and M0 above reflects the restriction to")
  say("     lines with copy number, not any adjustment)")
}

## ---- 4. What changes in the note --------------------------------------------
say("")
say("== 4. Changes required in the analysis note =====================")
say("  1. Methods: attribute the ~1e-36 floor to the NORMAL APPROXIMATION and")
say("     state the exact-test floor separately, with the reason the normal")
say("     approximation is used at this sample size.")
say("  2. Table 8: replace with the four models above, all on the same lines,")
say("     and report the all-lines M0 separately for reference.")
say("  3. Table 2 or Results 1: add the subtype breakdown of the 39 paired")
say("     lines from section 2.")
say("  4. Throughout: methylation-ASSOCIATED repression, not demonstrated")
say("     silencing; and state that whether promoter methylation constrains UV")
say("     inducibility is untested here.")
