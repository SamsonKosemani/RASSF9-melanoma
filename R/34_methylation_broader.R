# ==============================================================================
# 34_methylation_broader.R
#
# Does the RASSF9 promoter methylation - expression association hold in the
# broader TCGA-SKCM cohort?
#
# WHY. The primary analysis was restricted to 320 cutaneous melanomas with a
# detectable UV signature, because the question then was about UV burden. The
# methylation result (Spearman rho = -0.351, p = 1.1e-10, all 7 promoter probes
# significant) was the strongest human finding in that analysis, and a finding
# that strong should not depend on a filter chosen for a different purpose.
# Prof. Amarante-Mendes asked for it to be validated more broadly.
#
# WHAT "BROADER" MEANS HERE. Methylation is on disk for 470 samples, well
# beyond the 320, so no new download is needed. Three nested tiers are run:
#
#   T1  the original 320 (UV-signature-positive cutaneous) -- the comparator
#   T2  all cutaneous melanomas with both expression and methylation, including
#       the 91 with no detectable UV signature
#   T3  every sample with both measurements, whatever the cohort filters
#
# Uveal melanoma, which would be the ideal non-UV negative control, is NOT
# included: it is a separate GDC project (TCGA-UVM) and is not in the cohort
# files, so it would need a fresh download. Worth doing separately; not
# silently approximated here.
#
# PRE-SPECIFIED CRITERION, fixed before the numbers:
#   VALIDATED  = in T2, the promoter-mean association is in the same direction
#                and significant at alpha = 0.05, and the adjusted coefficient
#                keeps its sign.
#   WEAKENED   = same direction but no longer significant.
#   NOT SUPPORTED = direction reverses, or the effect vanishes.
# Writing this down first matters because a larger n almost always lowers a p
# value, so "still significant" is a weak test on its own. The informative
# comparison is whether the EFFECT SIZE holds, not whether p stays below 0.05.
#
# THREE CHECKS THE ORIGINAL ANALYSIS DID NOT RUN, added because a validation
# that only repeats the original test cannot find anything new:
#
#   1. NON-PROMOTER PROBES AS A SPECIFICITY CONTROL. The saved probe table
#      carries 23 probes, of which 18 are promoter. If promoter methylation
#      genuinely silences the gene, promoter probes should correlate negatively
#      with expression while probes elsewhere should not behave the same way.
#      If every probe on the locus correlates equally, the signal is about the
#      locus or the sample, not about the promoter.
#
#   2. WITHIN-STRATUM CONSISTENCY. RASSF9 promoter methylation differs between
#      primary and metastatic samples (0.522 vs 0.466, p = 0.013), and so does
#      expression. An association driven entirely by that split would be a
#      between-group artefact. The test is repeated within each stratum.
#
#   3. THE SAMPLES THE FILTERS EXCLUDED, analysed on their own. If the
#      association holds in the excluded samples too, the filter was irrelevant
#      to it, which is a stronger statement than "it holds in the larger set".
#
# Output: results/tcga/methylation_broader_tiers.csv
#         results/tcga/methylation_broader_probes.csv
#         results/tcga/methylation_broader_strata.csv
#         figures/tcga/tcga_<gene>_methylation_broader.{pdf,png,svg}
#         results/logs/34_methylation_broader.log
# ==============================================================================

source(here::here("R", "00_setup.R"))

LOG <- file.path(PATHS$logs, "34_methylation_broader.log")
cat("TCGA-SKCM methylation, broader cohort —", format(Sys.time()), "\n\n",
    file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

PROC <- PATHS$processed
RES <- file.path(PATHS$results, "tcga")
# NOT named FIG: 00_setup.R defines FIG as the figure-settings list, and a
# character vector of the same name shadows it, so FIG$base_size then fails
# inside every theme call.
FIG_TCGA <- file.path(PATHS$figures, "tcga")
invisible(lapply(c(RES, FIG_TCGA), dir.create, recursive = TRUE,
                 showWarnings = FALSE))

GENES <- c("RASSF9", "RASSF10")

pick <- function(slot, fallback) {
  if (!exists("PAL", inherits = TRUE)) return(fallback)
  x <- get("PAL", inherits = TRUE)[[slot]]
  if (!is.null(x) && is.character(x) && nzchar(x[1])) x[1] else fallback
}
COL_A    <- pick("nevus",    "#0095D9")
COL_B    <- pick("melanoma", "#F2680C")
COL_GREY <- pick("grey_pt",  "#8E8E8E")
COL_INK  <- pick("ink",      "#1A1A1A")

# ==============================================================================
# 1. Inputs
# ==============================================================================
hdr("1. Inputs")

f_probes <- file.path(PROC, "methylation_probes.csv")
f_betas  <- file.path(PROC, "methylation_betas.csv")
f_ana    <- file.path(PROC, "analysis_cohort.csv")
f_full   <- file.path(PROC, "full_cohort.csv")
for (f in c(f_probes, f_betas, f_ana, f_full)) {
  if (!file.exists(f)) stop("Missing ", f)
}

PR <- readr::read_csv(f_probes, show_col_types = FALSE, progress = FALSE)
BE <- readr::read_csv(f_betas,  show_col_types = FALSE, progress = FALSE)
AN <- readr::read_csv(f_ana,    show_col_types = FALSE, progress = FALSE)
FU <- readr::read_csv(f_full,   show_col_types = FALSE, progress = FALSE)

log_line("  probes:          ", nrow(PR), logfile = LOG)
log_line("  methylation rows: ", nrow(BE), logfile = LOG)
log_line("  analysis cohort:  ", nrow(AN), logfile = LOG)
log_line("  full cohort:      ", nrow(FU), logfile = LOG)

for (g in GENES) {
  sub <- PR[PR$gene == g, ]
  log_line(sprintf("  %-8s %d probes, %d promoter (TSS %+d to %+d bp)", g,
                   nrow(sub), sum(sub$promoter),
                   min(sub$dist_tss[sub$promoter]),
                   max(sub$dist_tss[sub$promoter])), logfile = LOG)
}
np <- PR$probe[!PR$promoter]
log_line("  non-promoter probes available as a specificity control: ",
         length(np), logfile = LOG)
if (!length(np)) {
  log_line("  [NOTE] none present; check 1 will be skipped", logfile = LOG)
}

PROBE_COLS <- intersect(PR$probe, names(BE))
log_line("  probes present in the beta table: ", length(PROBE_COLS), " of ",
         nrow(PR), logfile = LOG)

# ==============================================================================
# 2. Tiers
# ==============================================================================
hdr("2. Cohort tiers")

key_be <- if ("patient" %in% names(BE)) "patient" else names(BE)[1]
log_line("  joining methylation on: ", key_be, logfile = LOG)

mk <- function(df, label) {
  j <- dplyr::inner_join(df, BE, by = "patient", suffix = c("", ".meth"))
  j$tier <- label
  j
}

T1 <- mk(AN, "T1 primary (UV-positive cutaneous)")
T2 <- mk(FU, "T2 all cutaneous")

# T3: everything with methylation AND an expression value from either cohort.
EXPR <- dplyr::bind_rows(
  AN |> dplyr::select(dplyr::any_of(c("patient", "sample_type", "CPE",
                                      paste0(GENES, "_log2tpm")))),
  FU |> dplyr::select(dplyr::any_of(c("patient", "sample_type", "CPE",
                                      paste0(GENES, "_log2tpm"))))) |>
  dplyr::distinct(patient, .keep_all = TRUE)
T3 <- mk(EXPR, "T3 all with both measurements")

log_line("  T1: ", nrow(T1), " samples", logfile = LOG)
log_line("  T2: ", nrow(T2), " samples", logfile = LOG)
log_line("  T3: ", nrow(T3), " samples", logfile = LOG)
log_line("  methylation samples with no expression match: ",
         sum(!BE$patient %in% EXPR$patient), logfile = LOG)

# The samples T1 excluded, on their own. If the association holds here too, the
# filter was irrelevant to it.
EXCL <- T2 |> dplyr::filter(!patient %in% T1$patient)
log_line("  samples in T2 but not T1 (no UV signature): ", nrow(EXCL),
         logfile = LOG)

# ==============================================================================
# 3. The association, per tier
# ==============================================================================
hdr("3. Promoter methylation and expression, by tier")

promoter_mean <- function(d, g) {
  p <- intersect(PR$probe[PR$gene == g & PR$promoter], names(d))
  if (!length(p)) return(rep(NA_real_, nrow(d)))
  rowMeans(as.matrix(d[, p, drop = FALSE]), na.rm = TRUE)
}

run_tier <- function(d, g, label) {
  ycol <- paste0(g, "_log2tpm")
  if (!ycol %in% names(d)) return(NULL)
  d$.beta <- promoter_mean(d, g)
  d$.y <- suppressWarnings(as.numeric(d[[ycol]]))
  s <- d[is.finite(d$.beta) & is.finite(d$.y), ]
  if (nrow(s) < 20) {
    log_line("  ", g, " / ", label, ": only ", nrow(s), " samples; skipped",
             logfile = LOG)
    return(NULL)
  }
  ct <- suppressWarnings(stats::cor.test(s$.beta, s$.y, method = "spearman"))

  # Adjusted model, matching the original specification where the covariates
  # are available.
  # A covariate with one level inside a stratum is aliased and makes the whole
  # model unfittable, which is why the stratified adjusted betas came back NA.
  cov <- intersect(c("CPE", "sample_type"), names(s))
  cov <- cov[vapply(cov, function(cn)
    dplyr::n_distinct(s[[cn]][!is.na(s[[cn]])]) > 1, logical(1))]
  fml <- stats::as.formula(paste(".y ~ .beta",
                                 if (length(cov)) paste("+", paste(cov, collapse = " + ")) else ""))
  fit <- try(stats::lm(fml, data = s), silent = TRUE)
  b <- ci_lo <- ci_hi <- p_adj <- NA_real_
  if (!inherits(fit, "try-error")) {
    cf <- summary(fit)$coefficients
    if (".beta" %in% rownames(cf)) {
      b <- cf[".beta", 1]; p_adj <- cf[".beta", 4]
      ci <- try(stats::confint(fit, ".beta"), silent = TRUE)
      if (!inherits(ci, "try-error")) { ci_lo <- ci[1]; ci_hi <- ci[2] }
    }
  }

  log_line(sprintf("  %-8s %-36s n = %3d | rho = %+.3f, p = %-10s | adj beta = %+.3f [%.3f, %.3f], p = %s",
                   g, label, nrow(s), ct$estimate, fmt_p(ct$p.value),
                   b, ci_lo, ci_hi, fmt_p(p_adj)), logfile = LOG)

  tibble::tibble(gene = g, tier = label, n = nrow(s),
                 rho = unname(ct$estimate), rho_p = ct$p.value,
                 adj_beta = b, adj_lo = ci_lo, adj_hi = ci_hi, adj_p = p_adj,
                 median_beta = stats::median(s$.beta, na.rm = TRUE),
                 covariates = paste(cov, collapse = " + "))
}

tiers <- list()
for (g in GENES) {
  for (tt in list(list(T1, "T1 UV-positive cutaneous (original)"),
                  list(T2, "T2 all cutaneous"),
                  list(T3, "T3 all with both"),
                  list(EXCL, "Excluded from T1 only"))) {
    r <- run_tier(tt[[1]], g, tt[[2]])
    if (!is.null(r)) tiers[[length(tiers) + 1L]] <- r
  }
}
TI <- dplyr::bind_rows(tiers)
if (nrow(TI)) readr::write_csv(TI, file.path(RES, "methylation_broader_tiers.csv"))

# ==============================================================================
# 4. Verdict against the pre-specified criterion
# ==============================================================================
hdr("4. Verdict")

log_line("  CRITERION (fixed before the numbers):", logfile = LOG)
log_line("    VALIDATED     = T2 same direction, p < 0.05, adjusted sign kept",
         logfile = LOG)
log_line("    WEAKENED      = same direction, no longer significant",
         logfile = LOG)
log_line("    NOT SUPPORTED = direction reverses or effect vanishes",
         logfile = LOG)
log_line("  Note that a larger n lowers p almost automatically, so the",
         logfile = LOG)
log_line("  informative comparison is the EFFECT SIZE across tiers, not",
         logfile = LOG)
log_line("  whether p stayed below 0.05.", logfile = LOG)

for (g in GENES) {
  a <- TI |> dplyr::filter(gene == g, grepl("^T1", tier))
  b <- TI |> dplyr::filter(gene == g, grepl("^T2", tier))
  if (!nrow(a) || !nrow(b)) next
  same_dir <- sign(b$rho[1]) == sign(a$rho[1])
  sig <- b$rho_p[1] < 0.05
  sign_kept <- is.na(b$adj_beta[1]) || is.na(a$adj_beta[1]) ||
    sign(b$adj_beta[1]) == sign(a$adj_beta[1])
  verdict <- if (same_dir && sig && sign_kept) {
    "VALIDATED"
  } else if (same_dir) {
    "WEAKENED"
  } else {
    "NOT SUPPORTED"
  }
  log_line("", logfile = LOG)
  log_line(sprintf("  %-8s T1 rho %+.3f (n = %d)  ->  T2 rho %+.3f (n = %d)   %s",
                   g, a$rho[1], a$n[1], b$rho[1], b$n[1], verdict),
           logfile = LOG)
  log_line(sprintf("           effect size change: %+.3f (%.0f%% of the original)",
                   b$rho[1] - a$rho[1], 100 * b$rho[1] / a$rho[1]),
           logfile = LOG)
}

# ==============================================================================
# 5. Check 1 -- are promoter probes special?
# ==============================================================================
hdr("5. Specificity: promoter against non-promoter probes")

probe_rows <- list()
for (g in GENES) {
  ycol <- paste0(g, "_log2tpm")
  if (!ycol %in% names(T2)) next
  ps <- PR[PR$gene == g & PR$probe %in% names(T2), ]
  if (!nrow(ps)) next
  y <- suppressWarnings(as.numeric(T2[[ycol]]))
  for (i in seq_len(nrow(ps))) {
    v <- suppressWarnings(as.numeric(T2[[ps$probe[i]]]))
    ok <- is.finite(v) & is.finite(y)
    if (sum(ok) < 20) next
    ct <- suppressWarnings(stats::cor.test(v[ok], y[ok], method = "spearman"))
    probe_rows[[length(probe_rows) + 1L]] <- tibble::tibble(
      gene = g, probe = ps$probe[i], dist_tss = ps$dist_tss[i],
      region = if ("region" %in% names(ps)) ps$region[i] else NA_character_,
      promoter = ps$promoter[i], n = sum(ok),
      rho = unname(ct$estimate), p = ct$p.value,
      median_beta = stats::median(v[ok]))
  }
}
PB <- dplyr::bind_rows(probe_rows)
if (nrow(PB)) {
  PB <- PB |>
    dplyr::group_by(gene) |>
    dplyr::mutate(p_bh = stats::p.adjust(p, method = "BH")) |>
    dplyr::ungroup()
  readr::write_csv(PB, file.path(RES, "methylation_broader_probes.csv"))

  for (g in GENES) {
    sub <- PB |> dplyr::filter(gene == g) |> dplyr::arrange(dist_tss)
    if (!nrow(sub)) next
    log_line("", logfile = LOG)
    log_line("  ", g, " per probe (T2):", logfile = LOG)
    log_line(sprintf("    %-12s %9s %-12s %8s %8s %10s", "probe", "TSS dist",
                     "region", "med beta", "rho", "BH p"), logfile = LOG)
    for (i in seq_len(nrow(sub))) {
      log_line(sprintf("    %-12s %+9d %-12s %8.3f %+8.3f %10s",
                       sub$probe[i], sub$dist_tss[i],
                       ifelse(is.na(sub$region[i]), "-", sub$region[i]),
                       sub$median_beta[i], sub$rho[i], fmt_p(sub$p_bh[i])),
               logfile = LOG)
    }
    pr <- sub |> dplyr::filter(promoter)
    nq <- sub |> dplyr::filter(!promoter)
    if (!nrow(nq)) {
      log_line("    no non-promoter probes for this gene, so there is no",
               logfile = LOG)
      log_line("    external specificity control. The internal one is the",
               logfile = LOG)
      log_line("    island-shore pattern: probes near the TSS should correlate",
               logfile = LOG)
      log_line("    negatively while distal shore probes need not.",
               logfile = LOG)
    }
    if (nrow(pr) && nrow(nq)) {
      log_line(sprintf("    promoter probes: mean rho %+.3f (%d of %d significant)",
                       mean(pr$rho), sum(pr$p_bh < 0.05), nrow(pr)),
               logfile = LOG)
      log_line(sprintf("    non-promoter:    mean rho %+.3f (%d of %d significant)",
                       mean(nq$rho), sum(nq$p_bh < 0.05), nrow(nq)),
               logfile = LOG)
      # Judge the comparison rather than leaving the reader to. Opposite signs
      # are evidence OF specificity, not against it: gene-body methylation
      # commonly tracks expression positively while promoter methylation
      # tracks it negatively. Only a similar negative correlation everywhere
      # would undercut the silencing interpretation.
      if (sign(mean(pr$rho)) != sign(mean(nq$rho))) {
        log_line("    -> OPPOSITE SIGNS. Promoter and non-promoter methylation",
                 logfile = LOG)
        log_line("       relate to expression in opposite directions, which is",
                 logfile = LOG)
        log_line("       the textbook promoter-versus-gene-body pattern and is",
                 logfile = LOG)
        log_line("       evidence that the signal is promoter-specific.",
                 logfile = LOG)
      } else if (abs(mean(nq$rho)) > 0.7 * abs(mean(pr$rho))) {
        log_line("    -> SIMILAR. Methylation anywhere on the locus tracks",
                 logfile = LOG)
        log_line("       expression about as well, so the signal is not",
                 logfile = LOG)
        log_line("       specific to the promoter and the silencing",
                 logfile = LOG)
        log_line("       interpretation weakens considerably.", logfile = LOG)
      } else {
        log_line("    -> Promoter probes correlate substantially more strongly",
                 logfile = LOG)
        log_line("       than the rest of the locus, consistent with promoter",
                 logfile = LOG)
        log_line("       silencing.", logfile = LOG)
      }
    }
  }
}

# ==============================================================================
# 6. Check 2 -- does it hold within primary and within metastatic?
# ==============================================================================
hdr("6. Within-stratum consistency")

log_line("  RASSF9 promoter methylation and expression both differ by tumour",
         logfile = LOG)
log_line("  site, so an association driven entirely by that split would be a",
         logfile = LOG)
log_line("  between-group artefact rather than a within-tumour relationship.",
         logfile = LOG)

strata <- list()
if ("sample_type" %in% names(T2)) {
  for (g in GENES) {
    for (st in sort(unique(stats::na.omit(T2$sample_type)))) {
      r <- run_tier(T2[T2$sample_type == st, ], g, paste0("T2 / ", st))
      if (!is.null(r)) strata[[length(strata) + 1L]] <- r
    }
  }
  ST <- dplyr::bind_rows(strata)
  if (nrow(ST)) {
    readr::write_csv(ST, file.path(RES, "methylation_broader_strata.csv"))
    log_line("", logfile = LOG)
    log_line("  -> the association must hold INSIDE each stratum, not only",
             logfile = LOG)
    log_line("     across them, for the silencing interpretation to stand",
             logfile = LOG)
  }
} else {
  log_line("  sample_type not available; stratified check skipped",
           logfile = LOG)
}

# ==============================================================================
# 7. Figure
# ==============================================================================
hdr("7. Figure")

reset_legends("# TCGA methylation (broader cohort) figure legends\n")

for (g in GENES) {
  sub <- TI |> dplyr::filter(gene == g, !grepl("^Excluded", tier))
  pb  <- PB |> dplyr::filter(gene == g)
  if (!nrow(sub)) next

  pa <- ggplot(sub, aes(x = rho, y = factor(tier, levels = rev(tier)))) +
    geom_vline(xintercept = 0, linetype = "22", colour = COL_GREY,
               linewidth = 0.4) +
    geom_point(aes(size = n), colour = COL_B) +
    scale_size_continuous(range = c(1.8, 3.6), guide = "none") +
    labs(x = "Spearman rho", y = NULL) +
    theme_journal(grid = FALSE) +
    theme(panel.grid.major.x = element_line(linewidth = 0.25,
                                            colour = "grey92"),
          plot.margin = margin(4, 10, 4, 4, unit = "pt"))

  p_list <- list(pa)
  if (nrow(pb)) {
    pbp <- ggplot(pb, aes(x = dist_tss, y = rho)) +
      geom_hline(yintercept = 0, linetype = "22", colour = COL_GREY,
                 linewidth = 0.4) +
      geom_point(aes(colour = promoter, shape = p_bh < 0.05), size = 1.8) +
      scale_colour_manual(values = c(`TRUE` = COL_B, `FALSE` = COL_A),
                          labels = c(`TRUE` = "Promoter",
                                     `FALSE` = "Elsewhere"), name = NULL) +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                         labels = c(`TRUE` = "BH p < 0.05",
                                    `FALSE` = "not significant"), name = NULL) +
      labs(x = "Distance from TSS (bp)", y = "Spearman rho") +
      theme_journal(grid = TRUE) +
      theme(legend.position = "bottom",
            legend.text = element_text(size = FIG$base_size - 3),
            legend.key.size = unit(8, "pt"),
            plot.margin = margin(4, 10, 4, 4, unit = "pt"))
    p_list <- list(pa, pbp)
  }

  fig <- patchwork::wrap_plots(p_list, nrow = 1) +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = FIG$base_size + 2, face = "bold",
                                  colour = COL_INK))
  stem <- paste0("tcga_", tolower(g), "_methylation_broader")
  save_fig_to(fig, stem, FIG_TCGA, width = FIG$width_double, height = 88)
  log_line("  wrote: ", stem, logfile = LOG)

  t2 <- sub |> dplyr::filter(grepl("^T2", tier))
  add_legend(stem,
    paste0("Figure. ", g,
           " promoter methylation and expression across cohort definitions."),
    paste0("(A) Spearman correlation between mean promoter methylation and ",
           g, " expression in TCGA-SKCM, across nested cohort definitions; ",
           "point size is sample number. (B) Per-probe correlation against ",
           "distance from the transcription start site, with promoter and ",
           "non-promoter probes distinguished; filled points are significant ",
           "after Benjamini-Hochberg correction. ",
           if (nrow(t2)) paste0("In the full cutaneous cohort (n = ",
                                t2$n[1], "), rho = ",
                                sprintf("%.3f", t2$rho[1]), ", p = ",
                                fmt_p(t2$rho_p[1]), ". ") else "",
           "Non-promoter probes serve as a specificity control: if methylation ",
           "elsewhere on the locus tracked expression as strongly, the signal ",
           "would not be attributable to the promoter."))
}

# ==============================================================================
# 8. Reading this
# ==============================================================================
hdr("8. Reading this")

log_line("  WHAT WOULD MAKE THIS A STRONG VALIDATION:", logfile = LOG)
log_line("    - the effect size holds across tiers, not just the p value",
         logfile = LOG)
log_line("    - promoter probes behave differently from non-promoter probes",
         logfile = LOG)
log_line("    - the association survives within primary and within metastatic",
         logfile = LOG)
log_line("    - it is present in the samples the original filter excluded",
         logfile = LOG)
log_line("", logfile = LOG)
log_line("  WHAT WOULD UNDERCUT IT:", logfile = LOG)
log_line("    - every probe on the locus correlates equally, which would point",
         logfile = LOG)
log_line("      at a sample-level property (purity, cell composition) rather",
         logfile = LOG)
log_line("      than promoter silencing", logfile = LOG)
log_line("    - the association disappears within strata", logfile = LOG)
log_line("", logfile = LOG)
log_line("  NOT TESTED HERE: uveal melanoma, which would be the ideal non-UV",
         logfile = LOG)
log_line("  control, is a separate GDC project (TCGA-UVM) and is not in these",
         logfile = LOG)
log_line("  cohort files. It needs its own download and script.", logfile = LOG)

report_warnings(logfile = LOG)
write_session_info("34_methylation_broader")
log_line("\nFull log: ", LOG, logfile = NULL)
