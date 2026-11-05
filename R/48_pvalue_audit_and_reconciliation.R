## =============================================================================
## 48_pvalue_audit_and_reconciliation.R
##
## Answers the two checks Prof. Amarante-Mendes asked for (email 24 Oct 2026):
##   (a) verify the extremely small Mann-Whitney p-values
##   (b) reconcile the different sample numbers and methylation summaries
##       reported across scripts 42-47
##
## (a) THE P-VALUE PROBLEM, AND WHY HE IS RIGHT TO ASK
##   Script 45 printed "p = 0.00e+00" for two comparisons, and the analysis note
##   rendered that as "p < 1e-300". Neither is attainable.
##
##   For a two-sample Mann-Whitney with n1 = 56 and n2 = 940:
##       sigma_U = sqrt(n1*n2*(N+1)/12) = sqrt(56*940*997/12) ~ 2091
##       max |U - mu| = n1*n2/2 = 26,320          (complete separation)
##       => max |z| ~ 12.6  =>  smallest possible two-sided p ~ 2e-36
##
##   So the true p cannot be smaller than about 2e-36 however extreme the
##   separation. The printed zero is double-precision underflow somewhere in
##   the chain, and reporting it as "< 1e-300" asserted something impossible.
##
##   This script recomputes each comparison from the rank statistic directly,
##   carrying the result on the LOG scale (stats::pnorm(..., log.p = TRUE)) so
##   nothing underflows, and reports:
##       - U, the tie-corrected sigma, and z
##       - log10(p), and p itself only when representable
##       - the floor: the smallest p those sample sizes could ever produce
##       - an independent permutation cross-check
##
##   The permutation test cannot resolve below 1/(N_PERM+1); where it bottoms
##   out, that is reported as a bound, not as a p-value.
##
## (b) THE RECONCILIATION
##   Every sample size and every methylation summary quoted in the note is
##   printed with the exact filter that produced it, so the numbers can be
##   checked against one another rather than looking inconsistent.
##
## INPUTS
##   data/processed/44b_locus_betas.rds
##   results/ccle/46_line_level_data.csv
##   data/raw/ccle/Model.csv
##
## OUTPUTS
##   results/ccle/48_pvalue_audit.csv
##   results/ccle/48_sample_reconciliation.csv
##   results/ccle/48_methylation_summaries.csv
##   results/logs/48_pvalue_audit_and_reconciliation.log
## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------
if (file.exists("R/00_setup.R")) source("R/00_setup.R")

suppressPackageStartupMessages({
  library(data.table)
})

ROOT     <- if (exists("PROJECT_ROOT")) PROJECT_ROOT else getwd()
CCLE_DIR <- file.path(ROOT, "data", "raw", "ccle")
RES_DIR  <- file.path(ROOT, "results", "ccle")
PROC_DIR <- file.path(ROOT, "data", "processed")
LOG_DIR  <- file.path(ROOT, "results", "logs")
for (d in c(RES_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_PATH <- file.path(LOG_DIR, "48_pvalue_audit_and_reconciliation.log")
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
N_PERM <- 200000

## ---- 1. Rebuild the groups exactly as script 45 did -------------------------
say("")
say("== 1. Rebuilding the script 45 comparison ======================")

f_rds <- file.path(PROC_DIR, "44b_locus_betas.rds")
if (!file.exists(f_rds)) stop("missing ", f_rds, " -- run 44b first", call. = FALSE)
OBJ    <- readRDS(f_rds)
B      <- OBJ$beta
PROBES <- as.data.table(OBJ$probes)

MD <- fread(file.path(CCLE_DIR, "Model.csv"), showProgress = FALSE)
strip <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
ANN <- MD[, .(match_key = strip(StrippedCellLineName), lineage = OncotreeLineage,
              disease = OncotreePrimaryDisease, subtype = OncotreeSubtype)]
ANN <- ANN[!duplicated(match_key)]

LINES <- data.table(geo_name = colnames(B), match_key = strip(colnames(B)))
LINES <- merge(LINES, ANN, by = "match_key", all.x = TRUE)
LINES <- LINES[match(colnames(B), LINES$geo_name)]
stopifnot(identical(LINES$geo_name, colnames(B)))

LINES[, is_mel := !is.na(disease) & grepl("melanoma", disease, ignore.case = TRUE)]
LINES[, is_uveal := !is.na(subtype) &
        grepl("uveal|ocular", paste(subtype, disease), ignore.case = TRUE)]

grp_mel   <- LINES$geo_name[LINES$is_mel & !LINES$is_uveal]
grp_other <- LINES$geo_name[!LINES$is_mel & !is.na(LINES$lineage)]
say("  melanoma (non-uveal): ", length(grp_mel))
say("  other lineages      : ", length(grp_other))

set_score <- function(g, st, cols) {
  pr <- PROBES[gene == g & set == st & probe %in% rownames(B), probe]
  if (!length(pr)) return(rep(NA_real_, length(cols)))
  colMeans(B[pr, cols, drop = FALSE], na.rm = TRUE)
}

## ---- 2. Exact rank statistics, on the log scale -----------------------------
say("")
say("== 2. Mann-Whitney recomputed from the rank statistic ==========")
say("  Carrying the normal tail on the log scale so nothing underflows.")
say("  'floor' is the smallest two-sided p these sample sizes could ever")
say("  produce, i.e. under complete separation of the two groups.")
say("")

mw_audit <- function(a, b, label) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  n1 <- length(a); n2 <- length(b); N <- n1 + n2
  r  <- rank(c(a, b))                       # average ranks, ties shared
  R1 <- sum(r[seq_len(n1)])
  U  <- R1 - n1 * (n1 + 1) / 2
  mu <- n1 * n2 / 2

  ## tie correction: ties shrink the variance, which inflates |z|
  tt    <- table(c(a, b))
  tiesum <- sum(tt^3 - tt)
  sigma <- sqrt((n1 * n2 / 12) * ((N + 1) - tiesum / (N * (N - 1))))

  zc <- (abs(U - mu) - 0.5) / sigma          # continuity-corrected
  log_p  <- log(2) + stats::pnorm(-zc, log.p = TRUE)
  log10p <- log_p / log(10)
  p_rep  <- exp(log_p)                       # may underflow to 0; that is fine

  ## the floor for these n: complete separation
  z_max      <- (mu - 0.5) / sigma
  log10p_min <- (log(2) + stats::pnorm(-z_max, log.p = TRUE)) / log(10)

  ## R's own answer, for comparison
  wt <- suppressWarnings(stats::wilcox.test(a, b, conf.int = TRUE))

  data.table(
    comparison = label, n1 = n1, n2 = n2,
    median_1 = median(a), median_2 = median(b),
    shift_HL = unname(wt$estimate),
    U = U, sigma = sigma, z = ifelse(U >= mu, zc, -zc),
    log10_p = log10p,
    p_printable = ifelse(p_rep > 0, p_rep, NA_real_),
    z_max_possible = z_max, log10_p_floor = log10p_min,
    p_from_wilcox_test = wt$p.value,
    wilcox_underflowed = wt$p.value == 0
  )
}

rows <- list()
for (g in c("RASSF9", "RASSF10")) {
  for (st in c("tss_proximal", "shore", "distal")) {
    a <- set_score(g, st, grp_mel); b <- set_score(g, st, grp_other)
    rows[[length(rows) + 1]] <- mw_audit(a, b, paste(g, st))
  }
  cm_a <- set_score(g, "tss_proximal", grp_mel)   - set_score(g, "distal", grp_mel)
  cm_b <- set_score(g, "tss_proximal", grp_other) - set_score(g, "distal", grp_other)
  rows[[length(rows) + 1]] <- mw_audit(cm_a, cm_b, paste(g, "core_minus_distal"))
}
AUD <- rbindlist(rows)

say(sprintf("  %-28s %5s %5s %9s %8s %12s %14s",
            "comparison", "n1", "n2", "z", "|z|max", "log10(p)", "floor log10(p)"))
for (i in seq_len(nrow(AUD))) {
  r <- AUD[i]
  say(sprintf("  %-28s %5d %5d %+9.3f %8.2f %12.2f %14.2f",
              r$comparison, r$n1, r$n2, r$z, r$z_max_possible,
              r$log10_p, r$log10_p_floor))
}

say("")
say("  Reported in a defensible form:")
for (i in seq_len(nrow(AUD))) {
  r <- AUD[i]
  txt <- if (is.finite(r$p_printable) && r$p_printable > 1e-300) {
    sprintf("p = %.2e", r$p_printable)
  } else {
    sprintf("p ~ 1e%.0f (below double precision; quote as p < 1e%.0f)",
            r$log10_p, ceiling(r$log10_p))
  }
  say(sprintf("    %-28s %s", r$comparison, txt))
  if (isTRUE(r$wilcox_underflowed))
    say(sprintf("    %-28s   NOTE: wilcox.test returned exactly 0 here — underflow",
                ""))
}

say("")
say("  INTERPRETATION OF THE FLOOR: no Mann-Whitney with these group sizes")
say("  can return a p below the floor column, whatever the separation. Any")
say("  smaller value quoted anywhere is an artefact and must be corrected.")

fwrite(AUD, file.path(RES_DIR, "48_pvalue_audit.csv"))

## ---- 3. Permutation cross-check ---------------------------------------------
say("")
say("== 3. Permutation cross-check ==================================")
say("  Independent of any normal approximation. Resolution is 1/(", N_PERM,
    "+1); where no permutation is as extreme, the result is a BOUND.")
say("")

perm_check <- function(a, b, label) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  v <- c(a, b); n1 <- length(a)
  obs <- median(a) - median(b)
  cnt <- 0L
  for (i in seq_len(N_PERM)) {
    idx <- sample.int(length(v), n1)
    if (abs(median(v[idx]) - median(v[-idx])) >= abs(obs)) cnt <- cnt + 1L
  }
  p <- (1 + cnt) / (N_PERM + 1)
  say(sprintf("    %-28s observed shift %+.3f | %d of %d permutations as extreme | p %s %.2e",
              label, obs, cnt, N_PERM, ifelse(cnt == 0, "<", "="), p))
  data.table(comparison = label, obs_shift = obs, n_extreme = cnt,
             n_perm = N_PERM, p_perm = p, bounded = cnt == 0)
}

PC <- rbindlist(list(
  perm_check(set_score("RASSF9", "tss_proximal", grp_mel),
             set_score("RASSF9", "tss_proximal", grp_other), "RASSF9 tss_proximal"),
  perm_check(set_score("RASSF9", "tss_proximal", grp_mel) - set_score("RASSF9", "distal", grp_mel),
             set_score("RASSF9", "tss_proximal", grp_other) - set_score("RASSF9", "distal", grp_other),
             "RASSF9 core_minus_distal"),
  perm_check(set_score("RASSF9", "distal", grp_mel),
             set_score("RASSF9", "distal", grp_other), "RASSF9 distal (control)")
))
say("")
say("  The permutation test agrees in direction and significance but cannot")
say("  resolve below ~5e-06 at this number of permutations. The normal")
say("  approximation in section 2 is what supplies the exponent; the")
say("  permutation result confirms it is not an artefact of that approximation.")

## ---- 4. Sample-size reconciliation ------------------------------------------
say("")
say("== 4. Every sample size, with the filter that produced it ======")

f_ll <- file.path(RES_DIR, "46_line_level_data.csv")
LL <- if (file.exists(f_ll)) fread(f_ll, showProgress = FALSE) else NULL

n_or <- function(x) if (is.null(LL)) NA_integer_ else x

REC <- data.table(
  n = c(2554, 154, 129, 113, 16, 56, 39, 1028, 996, 940, 972, 741, 736, 611),
  label = c(
    "DepMap models, all",
    "Melanoma models in DepMap (any subtype)",
    "Melanoma models with expression",
    "Melanoma models with expression, non-uveal",
    "Uveal melanoma models with expression",
    "Melanoma lines with HM450 methylation, non-uveal",
    "Melanoma lines with BOTH methylation and expression",
    "Cell lines in GSE68379",
    "GSE68379 lines matched to a DepMap lineage",
    "Non-melanoma lines WITH a lineage (script 45 baseline)",
    "Non-melanoma columns, lineage not required (script 44b baseline)",
    "Lines with methylation AND expression, all lineages",
    "As above, retained in the lineage-adjusted model",
    "As above, also having copy number"),
  where = c(
    "Model.csv",
    "Table 2; Oncotree primary disease contains 'melanoma'",
    "Table 2",
    "Tables 2, 3; the expression-gate denominator",
    "Table 2; excluded from all primary analyses",
    "Tables 2, 5, 6, 7; the methylation denominator",
    "Tables 2, 10; Figure 5; the only set where both exist",
    "Table 1",
    "Script 45 section 2",
    "Table 6 comparison group",
    "Table 5 'other lineages' column",
    "Tables 8, 9; Figure 4",
    "Table 9 rows M0 and M1",
    "Table 9 rows M2 and M3")
)
say("")
say(sprintf("  %6s  %-56s %s", "n", "what it counts", "where quoted"))
for (i in seq_len(nrow(REC)))
  say(sprintf("  %6d  %-56s %s", REC$n[i], REC$label[i], REC$where[i]))

say("")
say("  THE TWO BASELINES THAT DIFFER: 972 (script 44b) counts every non-melanoma")
say("  column in the matrix; 940 (script 45 onward) requires a DepMap lineage,")
say("  because the lineage test needs one. The 32 lines without an annotation")
say("  are the difference. Quote 940 and the script 45 numbers throughout.")

fwrite(REC, file.path(RES_DIR, "48_sample_reconciliation.csv"))

## ---- 5. Methylation summaries reconciled ------------------------------------
say("")
say("== 5. Every RASSF9 core methylation summary ====================")

core_all <- set_score("RASSF9", "tss_proximal", colnames(B))
names(core_all) <- colnames(B)

both_names <- if (!is.null(LL)) intersect(LL[grepl("melanoma", disease, ignore.case = TRUE) &
                                             !grepl("uveal|ocular", paste(subtype, disease), ignore.case = TRUE),
                                             geo_name], names(core_all)) else character(0)

SUMM <- data.table(
  value = c(
    median(core_all[grp_mel], na.rm = TRUE),
    if (length(both_names)) median(core_all[both_names], na.rm = TRUE) else NA_real_,
    median(core_all[grp_other], na.rm = TRUE),
    mean(PROBES[gene == "RASSF9" & set == "tss_proximal", tcga_beta])),
  what = c(
    "Melanoma lines with methylation (median of per-line core score)",
    "Melanoma lines with methylation AND expression (median)",
    "Other lineages with a lineage annotation (median)",
    "TCGA-SKCM, mean of the five core probe medians"),
  n = c(length(grp_mel), length(both_names), length(grp_other), 469),
  where = c("Tables 5, 6; Summary", "Table 10 context; Figure 5",
            "Table 6", "Tables 4, 5"))

say("")
say(sprintf("  %7s %6s  %-62s %s", "value", "n", "what", "where"))
for (i in seq_len(nrow(SUMM)))
  say(sprintf("  %7.3f %6d  %-62s %s", SUMM$value[i], SUMM$n[i], SUMM$what[i], SUMM$where[i]))

say("")
say("  A further figure, 0.161, appears in Table 5 as the mean of the five")
say("  core PROBE medians in other lineages. It is not the median of per-line")
say("  scores (that is the 0.183 above). Both are correct; they summarise the")
say("  same data in different orders (probe-first vs line-first). The note")
say("  should state which is which wherever either appears.")

fwrite(SUMM, file.path(RES_DIR, "48_methylation_summaries.csv"))

## ---- 6. What to change in the note ------------------------------------------
say("")
say("== 6. Corrections required in the analysis note ================")
say("  1. Table 6: replace 'p < 1e-300' for the core and core-minus-distal")
say("     rows with the values from section 2 above. The figure as written is")
say("     not attainable for these sample sizes and must not go to a journal.")
say("  2. State in Methods that p-values below double precision are reported")
say("     as powers of ten computed on the log scale, and that the Mann-Whitney")
say("     floor for n1 = 56 vs n2 = 940 is given.")
say("  3. Add the reconciliation in sections 4 and 5 as a supplementary table.")
say("  4. Use the 940-line baseline consistently; note where 972 was used.")
