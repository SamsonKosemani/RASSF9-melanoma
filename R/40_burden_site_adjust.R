# ==============================================================================
# 40_burden_site_adjust.R
#
# Is the RASSF9-UV burden association a site effect in disguise?
#
# THE PROBLEM. Script 39b found RASSF9 falls specifically with UV-associated
# mutations (rho = -0.174, BH p = 0.001) and not with the non-UV remainder
# (+0.093), surviving purity adjustment at 100% of effect. But acral and
# mucosal melanomas are both UV-low AND biologically distinct, so a UV-burden
# gradient could partly be an anatomical-site gradient with nothing to do with
# UV. full_cohort.csv carries only sample_type, not site, so the question has
# been open.
#
# A SUBTLETY THAT WOULD OTHERWISE PRODUCE A CONFIDENT WRONG ANSWER. TCGA-SKCM
# is 311 metastatic and 100 primary samples. For a metastatic sample,
# site_of_resection_or_biopsy records where the METASTASIS was taken -- often a
# lymph node -- not where the melanoma arose. Only tissue_or_organ_of_origin or
# primary_site speaks to acral versus cutaneous. This script prefers those and
# reports which column it used, with its full value table, so the choice can be
# checked rather than trusted.
#
# NO DOWNLOAD IS ATTEMPTED FIRST. The expression object already carries clinical
# colData, so this looks there before reaching for the GDC. If no usable site
# column exists, the script stops and says what to download rather than
# proceeding on a worse variable.
#
# PRE-SPECIFIED, before the numbers:
#   SITE-INDEPENDENT = the UV effect retains more than 60% of its magnitude
#                      with site in the model at p < 0.05, AND the association
#                      holds within the cutaneous-only subset
#   SITE-CONFOUNDED  = the effect loses more than 60%, or significance, once
#                      site is included
#   UNRESOLVED       = site data too coarse to decide -- which is a real
#                      possibility, since TCGA records many cases as "Skin,
#                      NOS", and would itself be the honest finding
#
# Output: results/tcga/burden_site.csv
#         results/tcga/burden_site_models.csv
#         results/logs/40_burden_site_adjust.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(SummarizedExperiment))

LOG <- file.path(PATHS$logs, "40_burden_site_adjust.log")
cat("Burden vs site confound —", format(Sys.time()), "\n\n", file = LOG)
hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

PROC <- PATHS$processed
RES  <- file.path(PATHS$results, "tcga")
dir.create(RES, recursive = TRUE, showWarnings = FALSE)
RETAIN <- 0.60

# ==============================================================================
# 1. Cohort
# ==============================================================================
hdr("1. Cohort")

D <- readr::read_csv(file.path(PROC, "full_cohort.csv"), show_col_types = FALSE,
                     progress = FALSE)
D$nonUV <- D$TMB_sig - D$SBS7_total
D$nonUV[D$nonUV < 0] <- NA_real_
log_line("  ", nrow(D), " samples with burden, purity and expression",
         logfile = LOG)

# ==============================================================================
# 2. Find a site column in the expression object's clinical data
# ==============================================================================
hdr("2. Looking for anatomical site")

RNA_F <- file.path(PATHS$raw, "skcm_star_counts.rds")
if (!file.exists(RNA_F)) stop("Missing ", RNA_F)
log_line("  loading colData from ", basename(RNA_F), " -- about a minute",
         logfile = LOG)
SE <- readRDS(RNA_F)
CD <- as.data.frame(colData(SE))
rm(SE); gc(verbose = FALSE)
log_line("  clinical columns available: ", ncol(CD), logfile = LOG)

site_like <- grep("site|organ|origin|anatom|location|primary",
                  names(CD), value = TRUE, ignore.case = TRUE)
log_line("  site-like columns: ",
         if (length(site_like)) paste(site_like, collapse = ", ") else "NONE",
         logfile = LOG)

for (cn in site_like) {
  v <- as.character(CD[[cn]])
  u <- sort(table(v[!is.na(v)]), decreasing = TRUE)
  log_line("\n    ", cn, " — ", length(u), " distinct, ", sum(is.na(v)),
           " missing", logfile = LOG)
  for (i in seq_len(min(length(u), 14)))
    log_line(sprintf("      %-54s %4d", names(u)[i], as.integer(u[i])),
             logfile = LOG)
  if (length(u) > 14) log_line("      ... and ", length(u) - 14, " more",
                               logfile = LOG)
}

# Priority: origin beats resection site, because for a metastatic sample the
# resection site is where the metastasis was taken, not where the tumour arose.
PREF <- c("tissue_or_organ_of_origin", "primary_site", "site_of_resection_or_biopsy")
chosen <- NULL
for (cn in PREF) if (cn %in% names(CD)) {
  v <- as.character(CD[[cn]])
  if (dplyr::n_distinct(v[!is.na(v)]) > 1) { chosen <- cn; break }
}
if (is.null(chosen) && length(site_like)) {
  cand <- site_like[vapply(site_like, function(cn)
    dplyr::n_distinct(CD[[cn]][!is.na(CD[[cn]])]) > 1, logical(1))]
  if (length(cand)) chosen <- cand[1]
}
if (is.null(chosen)) {
  log_line("\n  [STOP] no usable site column in the expression object.",
           logfile = LOG)
  log_line("  Fetch it from the GDC instead:", logfile = LOG)
  log_line('    cl <- TCGAbiolinks::GDCquery_clinic("TCGA-SKCM", "clinical")',
           logfile = LOG)
  log_line('    readr::write_csv(cl, file.path(PATHS$raw, "skcm_clinical.csv"))',
           logfile = LOG)
  log_line("  then rerun this script.", logfile = LOG)
  quit(save = "no", status = 0)
}
log_line("\n  USING: ", chosen, logfile = LOG)
if (chosen == "site_of_resection_or_biopsy")
  log_line("  [FLAG] this is the RESECTION site. For metastatic samples it is",
           logfile = LOG)
if (chosen == "site_of_resection_or_biopsy")
  log_line("  where the metastasis was taken, not the primary site. Treat any",
           logfile = LOG)
if (chosen == "site_of_resection_or_biopsy")
  log_line("  conclusion from it as weak.", logfile = LOG)

SITE <- tibble::tibble(patient = substr(rownames(CD), 1, 12),
                       site_raw = as.character(CD[[chosen]])) |>
  dplyr::distinct(patient, .keep_all = TRUE)

# ==============================================================================
# 3. Classify site
# ==============================================================================
hdr("3. Classifying site")

s <- tolower(SITE$site_raw)
SITE$site <- dplyr::case_when(
  is.na(s)                                                    ~ "Unknown",
  grepl("nail|finger|toe|hand|foot|feet|palm|sole|heel|acral", s) ~ "Acral",
  grepl("mucos|vulva|vagin|anal|rect|nasal|sinus|oral|mouth|conjunctiv|penis|urethra", s) ~ "Mucosal",
  grepl("choroid|uvea|ciliary|iris", s)                       ~ "Uveal",
  grepl("lymph", s)                                           ~ "Lymph node",
  # "Skin, NOS" is 361 of 473 samples and says NOTHING about acral versus
  # non-acral, which is the whole question. An earlier version matched "skin"
  # and filed all of them as Cutaneous, producing a 329-sample "cutaneous
  # only" subset that may well contain acral cases and a SITE-INDEPENDENT
  # verdict that overstated what had been tested. It is uninformative here.
  grepl("^skin, ?nos$|^skin$", s)                             ~ "Skin NOS",
  grepl("skin|trunk|back|shoulder|arm|leg|thigh|scalp|neck|face|ear|chest|abdomen|thorax|limb|extremit", s) ~ "Cutaneous",
  grepl("nos|not reported|unknown|other", s)                  ~ "Unspecified",
  TRUE                                                        ~ "Other")

D <- dplyr::left_join(D, SITE, by = "patient")
tb <- sort(table(D$site, useNA = "ifany"), decreasing = TRUE)
for (i in seq_along(tb))
  log_line(sprintf("    %-14s %4d", names(tb)[i], as.integer(tb[i])),
           logfile = LOG)

# Only a site that distinguishes acral/mucosal from specified cutaneous counts
# as informative. "Skin NOS" does not.
inf <- sum(D$site %in% c("Acral", "Mucosal", "Cutaneous", "Uveal"), na.rm = TRUE)
log_line("  acral cases identified: ", sum(D$site == "Acral", na.rm = TRUE),
         " | mucosal: ", sum(D$site == "Mucosal", na.rm = TRUE),
         " | Skin NOS (uninformative): ", sum(D$site == "Skin NOS", na.rm = TRUE),
         logfile = LOG)
if (sum(D$site %in% c("Acral", "Mucosal"), na.rm = TRUE) == 0) {
  log_line("  [FLAG] NO acral or mucosal case is identifiable from this column.",
           logfile = LOG)
  log_line("  The site adjustment below therefore CANNOT address the acral",
           logfile = LOG)
  log_line("  confound; it only separates skin from lymph node and other. Any",
           logfile = LOG)
  log_line("  verdict must say so. The informative test is the one restricted",
           logfile = LOG)
  log_line("  to UV-POSITIVE tumours in 39b (rho -0.120, p = 0.033, n = 320),",
           logfile = LOG)
  log_line("  which is a gradient by degree of UV damage among UV-exposed",
           logfile = LOG)
  log_line("  melanomas rather than a contrast between exposed and shielded",
           logfile = LOG)
  log_line("  sites.", logfile = LOG)
}
log_line("  informative for the acral/mucosal question: ", inf, " of ", nrow(D),
         sprintf(" (%.0f%%)", 100 * inf / nrow(D)), logfile = LOG)
if (inf < 0.5 * nrow(D)) {
  log_line("  [FLAG] fewer than half the samples carry an informative site.",
           logfile = LOG)
  log_line("  TCGA records many melanomas as 'Skin, NOS', so the adjustment",
           logfile = LOG)
  log_line("  below is limited and the verdict may be UNRESOLVED. That is a",
           logfile = LOG)
  log_line("  real answer about the metadata, not a failure of the analysis.",
           logfile = LOG)
}

log_line("\n  site against UV signature -- the premise of the worry:",
         logfile = LOG)
if ("has_uv_signature" %in% names(D)) {
  ct <- table(D$site, D$has_uv_signature)
  log_line(sprintf("    %-14s %10s %10s %10s", "site", "no UV", "UV", "% UV"),
           logfile = LOG)
  for (r in rownames(ct)) {
    a <- ct[r, "FALSE"]; b <- ct[r, "TRUE"]
    log_line(sprintf("    %-14s %10d %10d %9.1f%%", r, a, b,
                     100 * b / max(a + b, 1)), logfile = LOG)
  }
  log_line("  If acral and mucosal are overwhelmingly UV-negative while",
           logfile = LOG)
  log_line("  cutaneous are UV-positive, site and UV burden are entangled and",
           logfile = LOG)
  log_line("  the adjustment below is the whole question.", logfile = LOG)
}

log_line("\n  RASSF9 and RASSF10 expression by site:", logfile = LOG)
BY <- D |> dplyr::group_by(site) |>
  dplyr::summarise(n = dplyr::n(),
                   RASSF9 = stats::median(RASSF9_log2tpm, na.rm = TRUE),
                   RASSF10 = stats::median(RASSF10_log2tpm, na.rm = TRUE),
                   SBS7 = stats::median(SBS7_total, na.rm = TRUE),
                   .groups = "drop") |> dplyr::arrange(dplyr::desc(n))
for (i in seq_len(nrow(BY)))
  log_line(sprintf("    %-14s n = %4d | RASSF9 %6.3f | RASSF10 %6.3f | median SBS7 %8.1f",
                   BY$site[i], BY$n[i], BY$RASSF9[i], BY$RASSF10[i], BY$SBS7[i]),
           logfile = LOG)
readr::write_csv(BY, file.path(RES, "burden_site.csv"))

for (g in c("RASSF9", "RASSF10")) {
  sub <- D[D$site %in% names(which(table(D$site) >= 10)), ]
  if (dplyr::n_distinct(sub$site) < 2) next
  kw <- suppressWarnings(stats::kruskal.test(sub[[paste0(g, "_log2tpm")]],
                                             factor(sub$site)))
  log_line(sprintf("    %s differs by site: Kruskal-Wallis p = %s", g,
                   fmt_p(kw$p.value)), logfile = LOG)
}

# ==============================================================================
# 4. Does the UV association survive site adjustment?
# ==============================================================================
hdr("4. Adjustment for site")

lg <- function(x) log10(x + 1)
BURD <- list("UV (SBS7)" = "SBS7_total", "UV fraction" = "SBS7_fraction",
             "Total burden" = "TMB_sig", "Non-UV burden" = "nonUV")

fit1 <- function(g, bn, dat, extra, lbl) {
  x <- dat[[BURD[[bn]]]]
  dd <- data.frame(y = dat[[paste0(g, "_log2tpm")]],
                   b = if (bn == "UV fraction") x else lg(x))
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

# Collapse rare categories so the model is not fitting one-sample levels.
D$site_f <- as.character(D$site)
small <- names(which(table(D$site_f) < 10))
D$site_f[D$site_f %in% small] <- "Other/small"
D$site_f <- factor(D$site_f)
log_line("  site factor levels used: ", paste(levels(D$site_f), collapse = ", "),
         logfile = LOG)

rows <- list()
for (g in c("RASSF9", "RASSF10")) for (bn in names(BURD)) {
  mods <- list(
    list("M0 unadjusted", list()),
    list("M1 + purity + type", list(CPE = D$CPE, sample_type = D$sample_type)),
    list("M2 + site", list(site = D$site_f)),
    list("M3 + site + purity + type",
         list(site = D$site_f, CPE = D$CPE, sample_type = D$sample_type)))
  for (m in mods) {
    r <- fit1(g, bn, D, m[[2]], m[[1]])
    if (!is.null(r)) rows[[length(rows) + 1L]] <- r
  }
}
AD <- dplyr::bind_rows(rows)

for (g in c("RASSF9", "RASSF10")) for (bn in names(BURD)) {
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
# 5. Cutaneous only -- the cleanest version of the test
# ==============================================================================
hdr("5. Within cutaneous melanomas only")

log_line("  Restricting to cutaneous cases removes acral and mucosal entirely,",
         logfile = LOG)
log_line("  so any surviving gradient cannot be a site effect between those",
         logfile = LOG)
log_line("  groups. It costs sample size.", logfile = LOG)

# Skin NOS is included here only because excluding it would leave almost
# nothing; the subset is therefore "skin of unspecified type", NOT verified
# non-acral skin, and the log says so.
CUT <- D[D$site %in% c("Cutaneous", "Skin NOS"), ]
log_line("  NOTE: this subset is skin of unspecified type. It is not verified",
         logfile = LOG)
log_line("  to exclude acral melanoma, because the metadata cannot identify",
         logfile = LOG)
log_line("  acral cases at all.", logfile = LOG)
log_line("  cutaneous samples: ", nrow(CUT), logfile = LOG)
cut_rows <- list()
if (nrow(CUT) >= 40) {
  for (g in c("RASSF9", "RASSF10")) for (bn in names(BURD)) {
    x <- CUT[[BURD[[bn]]]]; y <- CUT[[paste0(g, "_log2tpm")]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 30) next
    ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman"))
    r <- fit1(g, bn, CUT, list(CPE = CUT$CPE, sample_type = CUT$sample_type),
              "cutaneous, adjusted")
    log_line(sprintf("    %-8s %-14s n = %3d  rho = %+.3f  p = %-10s | adj beta %s",
                     g, bn, sum(ok), ct$estimate, fmt_p(ct$p.value),
                     if (is.null(r)) "n/a" else
                       sprintf("%+.4f (p = %s)", r$beta, fmt_p(r$p))),
             logfile = LOG)
    cut_rows[[length(cut_rows) + 1L]] <- tibble::tibble(
      gene = g, burden = bn, subset = "Cutaneous only", n = sum(ok),
      rho = unname(ct$estimate), rho_p = ct$p.value,
      adj_beta = if (is.null(r)) NA_real_ else r$beta,
      adj_p = if (is.null(r)) NA_real_ else r$p)
  }
} else {
  log_line("  too few cutaneous samples for a separate analysis", logfile = LOG)
}
CUTR <- dplyr::bind_rows(cut_rows)
readr::write_csv(dplyr::bind_rows(AD |> dplyr::mutate(subset = "All"),
                                  CUTR |> dplyr::mutate(model = "cutaneous")),
                 file.path(RES, "burden_site_models.csv"))

# ==============================================================================
# 6. Verdict
# ==============================================================================
hdr("6. Verdict against the pre-set criterion")

log_line("  SITE-INDEPENDENT = UV effect keeps > ", round(100 * RETAIN),
         "% with site in the model at", logfile = LOG)
log_line("                     p < 0.05, AND holds within cutaneous only",
         logfile = LOG)
log_line("  SITE-CONFOUNDED  = loses more than that, or significance",
         logfile = LOG)
log_line("  UNRESOLVED       = site data too coarse to decide", logfile = LOG)
log_line("", logfile = LOG)

for (g in c("RASSF9", "RASSF10")) {
  s <- AD |> dplyr::filter(gene == g, burden == "UV (SBS7)")
  b0 <- s$beta[s$model == "M0 unadjusted"]
  b3 <- s |> dplyr::filter(model == "M3 + site + purity + type")
  if (!length(b0) || !nrow(b3) || !is.finite(b0) || abs(b0) < 0.02) next
  ret <- b3$beta[1] / b0
  cu <- CUTR |> dplyr::filter(gene == g, burden == "UV (SBS7)")
  cut_ok <- nrow(cu) > 0 && is.finite(cu$rho_p[1]) && cu$rho_p[1] < 0.05 &&
            sign(cu$rho[1]) == sign(b0)

  no_acral <- sum(D$site %in% c("Acral", "Mucosal"), na.rm = TRUE) == 0
  verdict <- if (no_acral) {
    "UNRESOLVED for acral/mucosal -- no such case is identifiable in the metadata; see the UV-positive-only test in 39b"
  } else if (inf < 0.5 * nrow(D) && !isTRUE(cut_ok)) {
    "UNRESOLVED -- site metadata too coarse"
  } else if (ret >= RETAIN && b3$p[1] < 0.05 && isTRUE(cut_ok)) {
    "SITE-INDEPENDENT"
  } else if (ret < RETAIN || b3$p[1] >= 0.05) {
    "SITE-CONFOUNDED -- the UV gradient is partly anatomical"
  } else "PARTIAL -- unresolved"

  log_line("  ", g, " vs UV burden: ", verdict, logfile = LOG)
  log_line(sprintf("    unadjusted %+.4f -> with site %+.4f (%.0f%% retained, p = %s)",
                   b0, b3$beta[1], 100 * ret, fmt_p(b3$p[1])), logfile = LOG)
  if (nrow(cu))
    log_line(sprintf("    cutaneous only: n = %d, rho = %+.3f, p = %s",
                     cu$n[1], cu$rho[1], fmt_p(cu$rho_p[1])), logfile = LOG)
  log_line("", logfile = LOG)
}

log_line("  A caution on over-adjustment: if UV burden is largely a CONSEQUENCE",
         logfile = LOG)
log_line("  of where a melanoma arose, then site sits on the causal path rather",
         logfile = LOG)
log_line("  than beside it, and adjusting for it removes real signal along with",
         logfile = LOG)
log_line("  the confound. Surviving adjustment is therefore strong evidence;",
         logfile = LOG)
log_line("  failing it is not proof that the association is spurious. The",
         logfile = LOG)
log_line("  cutaneous-only analysis is the cleaner test for that reason.",
         logfile = LOG)

report_warnings(logfile = LOG)
write_session_info("40_burden_site_adjust")
log_line("\nFull log: ", LOG, logfile = NULL)
