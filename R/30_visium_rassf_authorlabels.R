# ==============================================================================
# 30_visium_rassf_authorlabels.R
#
# THE READOUT, on the authors' own spot labels.
#
# Script 29 joined Kreuger et al.'s published per-spot assignments to the GEO
# count matrices: 39,264 of their 39,264 spots matched, exactly, on
# section + 16-mer barcode. So every spot's identity here is THEIRS, taken from
# the fit reported in their paper, not our reproduction of it. That retires the
# largest limitation the Visium arm had.
#
# THE QUESTION. RASSF9 falls sharply between nevus and melanoma in bulk
# cohorts. Two explanations are compatible with that: melanocytes lose RASSF9
# as they transform, or melanoma lesions simply contain less of the stroma that
# carries RASSF9. In bulk tissue the two are inseparable. Here they are not:
# the nevus and the melanoma sit in the same section, in the same stroma, in
# the same patient, so a paired within-lesion comparison answers it.
#
# WHAT RODERICK SLIEKER OBSERVED, unprompted, on 14 Oct 2026: "I do not see
# this gene to be present in melanoma, but rather in the immune cells?" That is
# an independent reading of his own data, and it agrees with what we found in
# GSE230574, where RASSF9 molecules in nevi were ~44% fibroblast and ~18%
# endothelial against ~10% melanocytic. Section 4 answers his question
# quantitatively across all fifteen of their cell types, because an informal
# impression and a measured distribution are not the same thing and the second
# is what belongs in a paper.
#
# DESIGN, fixed before any RASSF9 value was seen:
#
#   GATE (section 3). PRAME and CDK2 are the two genes the authors themselves
#   report as separating their tumour from their nevus. If the join were wrong
#   -- a mismatched barcode, a scrambled section -- those genes would NOT
#   separate. So the gate is run first, on marker genes, and the RASSF9
#   comparison is only run if it passes. The bar is set here, before the
#   result: PRAME higher in tumour in at least 75% of patients.
#
#   PRIMARY (section 5). Pseudobulk, paired within patient. Counts summed
#   across a patient's nevus spots and separately across their melanoma spots,
#   expressed per million. ONE paired observation per patient; n = patients,
#   not spots. Wilcoxon signed-rank, with a paired t on log2 CPM and a sign
#   test alongside. 39,000 spots from 16 patients are not 39,000 independent
#   observations, and a spot-level p value would be pseudo-replication.
#
#   Pairing also disposes of age, sex, Breslow thickness and site: they are
#   constant within a patient and cancel in a within-patient difference.
#
#   SENSITIVITY (section 6). Repeated at purity >= 0.75 using THEIR theta, and
#   with Tumor 1 and Tumor 2 kept separate, since pooling them is our choice
#   and ought to be shown not to matter.
#
#   SECONDARY (section 7). Spot-level zero-inflated negative binomial with
#   patient as a random intercept -- what they used. Direction and coefficient
#   only; the p value is not the result.
#
#   BOTH GENES, reported whatever comes out, on the same footing.
#
# A NULL RESULT IS A RESULT HERE. If RASSF9 does not differ between the two
# melanocytic compartments of the same lesion, that is evidence the bulk fall
# is compositional -- which is what we already demonstrated for RASSF10 -- and
# it is reported as such, with the confidence interval, not buried.
#
# Output: results/geo/gse298774_authorlabels_perpatient.csv
#         results/geo/gse298774_authorlabels_tests.csv
#         results/geo/gse298774_authorlabels_lineage.csv
#         figures/geo/geo_<gene>_visium_paired.{pdf,png,svg}
#         figures/geo/geo_<gene>_visium_lineage.{pdf,png,svg}
#         results/logs/30_visium_rassf_authorlabels.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(Matrix)
  library(patchwork)
})

LOG <- file.path(PATHS$logs, "30_visium_rassf_authorlabels.log")
cat("GSE298774 readout on author labels —", format(Sys.time()), "\n\n",
    file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

EXT_DIR <- file.path(PATHS$raw, "geo", "GSE298774", "extracted")
PROC    <- PATHS$processed
RES_GEO <- file.path(PATHS$results, "geo")
FIG_GEO <- file.path(PATHS$figures, "geo")
invisible(lapply(c(RES_GEO, FIG_GEO), dir.create, recursive = TRUE,
                 showWarnings = FALSE))

GENES     <- c("RASSF9", "RASSF10")
CONTROLS  <- c("PRAME", "CDK2", "MLANA", "PMEL", "PTPRC", "COL1A1")
MIN_SPOTS <- 10
PURITY_HI <- 0.75
GATE_BAR  <- 0.75          # PRAME must be up in tumour in >= 75% of patients

pick <- function(slot, fallback) {
  if (!exists("PAL", inherits = TRUE)) return(fallback)
  x <- get("PAL", inherits = TRUE)[[slot]]
  if (!is.null(x) && is.character(x) && nzchar(x[1])) x[1] else fallback
}
COL_NEV  <- pick("nevus",    "#0095D9")
COL_TUM  <- pick("melanoma", "#F2680C")
COL_GREY <- pick("grey_pt",  "#8E8E8E")
COL_INK  <- pick("ink",      "#1A1A1A")

# ==============================================================================
# 1. Labels
# ==============================================================================
hdr("1. Author labels")

LF <- file.path(PROC, "gse298774_kreuger_labels.rds")
if (!file.exists(LF)) stop("Run R/29_kreuger_metadata.R first (missing ", LF, ")")
KK <- readRDS(LF)
L  <- KK$labels

log_line("  spots with author labels: ", nrow(L), logfile = LOG)
log_line("  sections: ", dplyr::n_distinct(L$section), " | patients: ",
         dplyr::n_distinct(L$patient), logfile = LOG)
log_line("  source: ", KK$source$received, logfile = LOG)

NEVUS_LAB  <- "Nevus"
TUMOUR_LAB <- grep("^Tumor", unique(L$label), value = TRUE)
log_line("  nevus label: ", NEVUS_LAB, " | tumour labels: ",
         paste(TUMOUR_LAB, collapse = ", "), logfile = LOG)

L <- L |>
  dplyr::mutate(
    grp = dplyr::case_when(label == NEVUS_LAB  ~ "nevus",
                           label %in% TUMOUR_LAB ~ "tumour",
                           TRUE ~ NA_character_))

# ==============================================================================
# 2. Expression, from the GEO matrices
# ==============================================================================
hdr("2. Expression")

files <- list.files(EXT_DIR, full.names = TRUE)
sec_of <- function(b) sub("^GSM\\d+_(.+?)_(aligned_fiducials|barcodes|features|matrix|scalefactors_json|tissue_hires_image|tissue_lowres_image|tissue_positions).*$", "\\1", b)
inv <- tibble::tibble(path = files, base = basename(files)) |>
  dplyr::mutate(section = sec_of(base),
                role = dplyr::case_when(
                  grepl("_matrix\\.mtx", base)   ~ "matrix",
                  grepl("_barcodes\\.tsv", base) ~ "barcodes",
                  grepl("_features\\.tsv", base) ~ "features",
                  TRUE ~ "other")) |>
  dplyr::filter(role != "other")
g2 <- function(sec, role) {
  p <- inv$path[inv$section == sec & inv$role == role]
  if (!length(p)) NA_character_ else p[1]
}
core16 <- function(x) sub("-\\d+$", "", x)

WANT <- unique(c(GENES, CONTROLS))
vals <- list()
for (sec in sort(unique(L$section))) {
  mp <- g2(sec, "matrix"); bp <- g2(sec, "barcodes"); fp <- g2(sec, "features")
  if (any(is.na(c(mp, bp, fp)))) {
    log_line("  [FLAG] ", sec, ": files missing; skipped", logfile = LOG); next
  }
  m  <- Matrix::readMM(mp)
  bc <- readr::read_tsv(bp, col_names = FALSE, show_col_types = FALSE,
                        progress = FALSE)[[1]]
  ft <- readr::read_tsv(fp, col_names = FALSE, show_col_types = FALSE,
                        progress = FALSE)
  if (ncol(m) == nrow(ft) && nrow(m) == length(bc)) m <- Matrix::t(m)
  if (nrow(m) != nrow(ft) || ncol(m) != length(bc)) {
    log_line("  [FLAG] ", sec, ": dimension mismatch; skipped", logfile = LOG)
    rm(m); next
  }
  sy  <- toupper(as.character(ft[[min(2, ncol(ft))]]))
  out <- tibble::tibble(spot = paste0(sec, "|", core16(bc)),
                        lib_size = as.numeric(Matrix::colSums(m)))
  for (g in WANT) {
    i <- which(sy == g)
    out[[g]] <- if (!length(i)) {
      NA_real_
    } else if (length(i) > 1) {
      as.numeric(Matrix::colSums(m[i, , drop = FALSE]))
    } else {
      as.numeric(m[i, ])
    }
  }
  vals[[sec]] <- out
  rm(m); gc(verbose = FALSE)
}
E <- dplyr::bind_rows(vals)
d <- L |> dplyr::inner_join(E, by = "spot")
log_line("  spots with labels and counts: ", nrow(d), " of ", nrow(L),
         sprintf(" (%.1f%%)", 100 * nrow(d) / nrow(L)), logfile = LOG)
if (nrow(d) < 0.95 * nrow(L)) {
  log_line("  [FLAG] more than 5% of labelled spots lost in the expression",
           logfile = LOG)
  log_line("  join. Script 29 matched 100%, so something differs here.",
           logfile = LOG)
}

for (g in WANT) {
  if (all(is.na(d[[g]]))) {
    log_line(sprintf("  %-8s NOT on the annotation", g), logfile = LOG)
  } else {
    log_line(sprintf("  %-8s detected in %5.2f%% of spots, %8.0f molecules total",
                     g, 100 * mean(d[[g]] > 0, na.rm = TRUE),
                     sum(d[[g]], na.rm = TRUE)), logfile = LOG)
  }
}

# Helper: pseudobulk CPM by any grouping.
pb_cpm <- function(dat, keys, genes) {
  dat |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(n_spots = dplyr::n(), total_umi = sum(lib_size),
                     dplyr::across(dplyr::all_of(genes),
                                   ~ sum(.x, na.rm = TRUE), .names = "{.col}_n"),
                     .groups = "drop") |>
    dplyr::mutate(dplyr::across(dplyr::ends_with("_n"),
                                ~ 1e6 * .x / pmax(total_umi, 1),
                                .names = "{.col}_cpm"))
}

# ==============================================================================
# 3. THE GATE: do the authors' labels behave as the authors say?
# ==============================================================================
# Run on marker genes, before any RASSF9 number is produced. If the join were
# wrong, PRAME would not separate tumour from nevus.
hdr("3. Gate: marker genes on the author labels")

paired_frame <- function(dat, gene, extra_filter = TRUE) {
  s <- dat |> dplyr::filter(!is.na(grp), extra_filter)
  p <- pb_cpm(s, c("patient", "grp"), gene)
  nm <- paste0(gene, "_n_cpm")
  p |>
    dplyr::select(patient, grp, n_spots, cpm = dplyr::all_of(nm)) |>
    tidyr::pivot_wider(names_from = grp, values_from = c(cpm, n_spots)) |>
    dplyr::filter(!is.na(cpm_nevus), !is.na(cpm_tumour),
                  n_spots_nevus >= MIN_SPOTS, n_spots_tumour >= MIN_SPOTS) |>
    dplyr::mutate(diff = cpm_tumour - cpm_nevus)
}

gate_ok <- TRUE
for (g in CONTROLS) {
  if (all(is.na(d[[g]]))) {
    log_line(sprintf("  %-8s not on the annotation; skipped", g), logfile = LOG)
    next
  }
  w <- paired_frame(d, g)
  up <- sum(w$diff > 0)
  wt <- suppressWarnings(stats::wilcox.test(w$cpm_tumour, w$cpm_nevus,
                                            paired = TRUE))
  log_line(sprintf("  %-8s nevus %8.2f CPM, tumour %8.2f CPM | higher in tumour %2d/%2d | p = %s",
                   g, stats::median(w$cpm_nevus), stats::median(w$cpm_tumour),
                   up, nrow(w), fmt_p(wt$p.value)), logfile = LOG)
  if (g == "PRAME") {
    gate_ok <- (up / nrow(w)) >= GATE_BAR
    log_line("", logfile = LOG)
    log_line(sprintf("  GATE on PRAME: higher in tumour in %d of %d patients (%.0f%%), bar is %.0f%%",
                     up, nrow(w), 100 * up / nrow(w), 100 * GATE_BAR),
             logfile = LOG)
  }
}

if (!gate_ok) {
  log_line("", logfile = LOG)
  log_line("  [STOP] PRAME does not separate the authors' tumour spots from",
           logfile = LOG)
  log_line("  their nevus spots. Their paper reports that it does, so the",
           logfile = LOG)
  log_line("  failure is ours: most likely the barcode join or the section",
           logfile = LOG)
  log_line("  mapping. No RASSF9 comparison is run on a join that cannot",
           logfile = LOG)
  log_line("  reproduce the authors' own result.", logfile = LOG)
  report_warnings(logfile = LOG)
  stop("Marker gate failed; see ", LOG)
}
log_line("  gate passed -- the join reproduces the authors' own marker result",
         logfile = LOG)

# ==============================================================================
# 4. Where the genes actually sit, across all fifteen cell types
# ==============================================================================
# This answers Roderick Slieker's question directly and quantitatively.
hdr("4. Expression by cell type")

lin <- pb_cpm(d, "label", WANT)
lin_long <- lin |>
  dplyr::select(label, n_spots, total_umi, dplyr::ends_with("_n_cpm")) |>
  tidyr::pivot_longer(dplyr::ends_with("_n_cpm"), names_to = "gene",
                      values_to = "cpm") |>
  dplyr::mutate(gene = sub("_n_cpm$", "", gene))
readr::write_csv(lin_long, file.path(RES_GEO,
                                     "gse298774_authorlabels_lineage.csv"))

for (g in GENES) {
  sub <- lin_long |> dplyr::filter(gene == g) |> dplyr::arrange(dplyr::desc(cpm))
  if (all(is.na(sub$cpm))) next
  tot <- d |>
    dplyr::group_by(label) |>
    dplyr::summarise(mol = sum(.data[[g]], na.rm = TRUE), .groups = "drop")
  tot$share <- 100 * tot$mol / sum(tot$mol)
  sub <- dplyr::left_join(sub, tot, by = "label")
  log_line("", logfile = LOG)
  log_line("  ", g, " by author cell type:", logfile = LOG)
  log_line(sprintf("    %-22s %10s %10s %8s %8s", "cell type", "CPM",
                   "molecules", "% of all", "spots"), logfile = LOG)
  for (i in seq_len(nrow(sub))) {
    log_line(sprintf("    %-22s %10.2f %10.0f %7.1f%% %8d", sub$label[i],
                     sub$cpm[i], sub$mol[i], sub$share[i], sub$n_spots[i]),
             logfile = LOG)
  }
  mel <- sub |> dplyr::filter(label %in% c(NEVUS_LAB, TUMOUR_LAB))
  log_line(sprintf("    -> melanocytic compartments (%s) hold %.1f%% of all %s molecules",
                   paste(c(NEVUS_LAB, TUMOUR_LAB), collapse = " + "),
                   sum(mel$share, na.rm = TRUE), g), logfile = LOG)
}

# ==============================================================================
# 5. PRIMARY: pseudobulk, paired within patient
# ==============================================================================
hdr("5. Primary: paired within patient")

run_paired <- function(w, g, label) {
  if (nrow(w) < 4) {
    log_line("  ", g, " / ", label, ": only ", nrow(w),
             " usable patients -- not tested", logfile = LOG)
    return(NULL)
  }
  wt <- suppressWarnings(stats::wilcox.test(w$cpm_tumour, w$cpm_nevus,
                                            paired = TRUE))
  tt <- try(stats::t.test(log2(w$cpm_tumour + 0.1), log2(w$cpm_nevus + 0.1),
                          paired = TRUE), silent = TRUE)
  st <- stats::binom.test(sum(w$diff < 0), nrow(w), 0.5)

  log_line("", logfile = LOG)
  log_line("  ", g, " / ", label, " -- ", nrow(w), " paired patients",
           logfile = LOG)
  log_line(sprintf("    %-10s %10s %10s %11s %7s %7s", "patient", "nevus CPM",
                   "tum CPM", "difference", "n nev", "n tum"), logfile = LOG)
  for (i in order(w$patient)) {
    log_line(sprintf("    %-10s %10.2f %10.2f %+11.2f %7d %7d", w$patient[i],
                     w$cpm_nevus[i], w$cpm_tumour[i], w$diff[i],
                     w$n_spots_nevus[i], w$n_spots_tumour[i]), logfile = LOG)
  }
  log_line("", logfile = LOG)
  log_line(sprintf("    median nevus %.2f CPM, median tumour %.2f CPM",
                   stats::median(w$cpm_nevus), stats::median(w$cpm_tumour)),
           logfile = LOG)
  log_line(sprintf("    lower in tumour in %d of %d patients",
                   sum(w$diff < 0), nrow(w)), logfile = LOG)
  log_line(sprintf("    Wilcoxon signed-rank: V = %.1f, p = %s",
                   wt$statistic, fmt_p(wt$p.value)), logfile = LOG)
  if (!inherits(tt, "try-error")) {
    log_line(sprintf("    paired t on log2 CPM: %+.3f [%.3f, %.3f], p = %s",
                     tt$estimate, tt$conf.int[1], tt$conf.int[2],
                     fmt_p(tt$p.value)), logfile = LOG)
    log_line(sprintf("    -> fold change %.2fx [%.2f, %.2f]",
                     2^tt$estimate, 2^tt$conf.int[1], 2^tt$conf.int[2]),
             logfile = LOG)
  }
  log_line(sprintf("    sign test: %d/%d, p = %s", sum(w$diff < 0), nrow(w),
                   fmt_p(st$p.value)), logfile = LOG)

  list(stats = tibble::tibble(
         gene = g, set = label, n_patients = nrow(w),
         median_nevus_cpm = stats::median(w$cpm_nevus),
         median_tumour_cpm = stats::median(w$cpm_tumour),
         n_lower_in_tumour = sum(w$diff < 0),
         wilcoxon_V = unname(wt$statistic), wilcoxon_p = wt$p.value,
         sign_p = st$p.value,
         log2fc = if (inherits(tt, "try-error")) NA_real_ else unname(tt$estimate),
         ci_lo = if (inherits(tt, "try-error")) NA_real_ else tt$conf.int[1],
         ci_hi = if (inherits(tt, "try-error")) NA_real_ else tt$conf.int[2],
         t_p = if (inherits(tt, "try-error")) NA_real_ else tt$p.value),
       wide = w)
}

tests <- list(); wides <- list(); pbs <- list()
for (g in GENES) {
  if (all(is.na(d[[g]]))) {
    log_line("  ", g, ": not on the annotation; skipped", logfile = LOG); next
  }
  w <- paired_frame(d, g)
  r <- run_paired(w, g, "all labelled spots")
  if (!is.null(r)) {
    tests[[length(tests) + 1L]] <- r$stats
    wides[[g]] <- r$wide
  }
}
pb_all <- pb_cpm(dplyr::filter(d, !is.na(grp)), c("patient", "grp"), WANT)
readr::write_csv(pb_all,
                 file.path(RES_GEO, "gse298774_authorlabels_perpatient.csv"))

# ==============================================================================
# 6. Sensitivity
# ==============================================================================
hdr("6. Sensitivity")

log_line("  (a) high-purity spots only, purity >= ", PURITY_HI,
         " on the authors' theta", logfile = LOG)
for (g in GENES) {
  if (all(is.na(d[[g]]))) next
  w <- paired_frame(d, g, extra_filter = d$purity >= PURITY_HI)
  r <- run_paired(w, g, paste0("purity >= ", PURITY_HI))
  if (!is.null(r)) tests[[length(tests) + 1L]] <- r$stats
}

log_line("", logfile = LOG)
log_line("  (b) Tumor 1 and Tumor 2 separately, since pooling them was our",
         logfile = LOG)
log_line("      choice and should be shown not to drive the answer",
         logfile = LOG)
for (tl in TUMOUR_LAB) {
  dd <- d |>
    dplyr::mutate(grp = dplyr::case_when(label == NEVUS_LAB ~ "nevus",
                                         label == tl ~ "tumour",
                                         TRUE ~ NA_character_))
  for (g in GENES) {
    if (all(is.na(dd[[g]]))) next
    w <- paired_frame(dd, g)
    r <- run_paired(w, g, paste0("Nevus vs ", tl))
    if (!is.null(r)) tests[[length(tests) + 1L]] <- r$stats
  }
}

TT <- dplyr::bind_rows(tests)
if (nrow(TT)) {
  readr::write_csv(TT, file.path(RES_GEO, "gse298774_authorlabels_tests.csv"))
  log_line("", logfile = LOG)
  log_line("  wrote: results/geo/gse298774_authorlabels_tests.csv", logfile = LOG)
}

# ==============================================================================
# 7. Secondary: spot-level mixed model
# ==============================================================================
hdr("7. Secondary: spot-level ZINB with patient random effect")

if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  log_line("  glmmTMB not installed; attempting install", logfile = LOG)
  try(utils::install.packages("glmmTMB", repos = "https://cloud.r-project.org"),
      silent = TRUE)
}
if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  log_line("  glmmTMB unavailable -- skipped. The paired test in section 5",
           logfile = LOG)
  log_line("  stands on its own; this was supporting evidence only.",
           logfile = LOG)
} else {
  for (g in GENES) {
    if (all(is.na(d[[g]]))) next
    dd <- d |>
      dplyr::filter(!is.na(grp), !is.na(.data[[g]])) |>
      dplyr::mutate(y = .data[[g]],
                    grp = factor(grp, levels = c("nevus", "tumour")),
                    off = log(pmax(lib_size, 1)), pid = factor(patient))
    if (dplyr::n_distinct(dd$pid) < 4) next
    fit <- try(glmmTMB::glmmTMB(y ~ grp + offset(off) + (1 | pid),
                                ziformula = ~ 1, family = glmmTMB::nbinom2,
                                data = dd), silent = TRUE)
    if (inherits(fit, "try-error")) {
      log_line("  ", g, ": model did not converge", logfile = LOG); next
    }
    cf  <- summary(fit)$coefficients$cond
    row <- grep("^grptumour$", rownames(cf))
    if (!length(row)) next
    log_line(sprintf("  %-8s beta %+.3f (SE %.3f), fold %.2fx, p = %s   [%d spots, %d patients]",
                     g, cf[row, 1], cf[row, 2], exp(cf[row, 1]),
                     fmt_p(cf[row, 4]), nrow(dd), dplyr::n_distinct(dd$pid)),
             logfile = LOG)
  }
  log_line("", logfile = LOG)
  log_line("  Read the direction and the coefficient here, not the p value:",
           logfile = LOG)
  log_line("  spots within one lesion are not independent. The paired",
           logfile = LOG)
  log_line("  patient-level test is the result.", logfile = LOG)
}

# ==============================================================================
# 8. Figures
# ==============================================================================
hdr("8. Figures")

reset_legends("# GSE298774 figure legends\n")

for (g in GENES) {
  w <- wides[[g]]
  if (is.null(w) || !nrow(w)) next

  long <- dplyr::bind_rows(
    tibble::tibble(patient = w$patient, grp = "Nevus", cpm = w$cpm_nevus),
    tibble::tibble(patient = w$patient, grp = "Melanoma", cpm = w$cpm_tumour)) |>
    dplyr::mutate(grp = factor(grp, levels = c("Nevus", "Melanoma")))

  pa <- ggplot(long, aes(x = grp, y = cpm, group = patient)) +
    geom_line(colour = COL_GREY, linewidth = 0.35, alpha = 0.8) +
    geom_point(aes(colour = grp), size = 1.4) +
    scale_colour_manual(values = c(Nevus = COL_NEV, Melanoma = COL_TUM),
                        guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = bquote(bold(.(g)~"(counts per million)"))) +
    theme_journal(grid = TRUE) +
    theme(plot.margin = margin(4, 10, 4, 4, unit = "pt"))

  pbf <- ggplot(w, aes(x = reorder(patient, diff), y = diff)) +
    geom_hline(yintercept = 0, linetype = "22", colour = COL_GREY,
               linewidth = 0.4) +
    geom_col(aes(fill = diff < 0), width = 0.66, colour = COL_INK,
             linewidth = 0.25) +
    scale_fill_manual(values = c(`TRUE` = COL_NEV, `FALSE` = COL_TUM),
                      guide = "none") +
    coord_flip(clip = "off") +
    labs(x = NULL, y = "Melanoma − nevus (CPM)") +
    theme_journal(grid = FALSE) +
    theme(panel.grid.major.x = element_line(linewidth = 0.25,
                                            colour = "grey92"),
          axis.text.y = element_text(size = FIG$base_size - 1),
          plot.margin = margin(4, 10, 4, 4, unit = "pt"))

  fig <- (pa | pbf) + patchwork::plot_layout(widths = c(1, 1.25)) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = FIG$base_size + 3, face = "bold",
                                  colour = COL_INK))
  stem <- paste0("geo_", tolower(g), "_visium_paired")
  save_fig_to(fig, stem, FIG_GEO, width = FIG$width_double, height = 100)

  st <- TT |> dplyr::filter(gene == g, set == "all labelled spots")
  add_legend(stem,
    paste0("Figure. ", g, " in the nevus and melanoma compartments of the same lesion."),
    paste0("Spatial transcriptomics of nevus-associated melanoma (GSE298774, ",
           "Kreuger et al.). Spot identities are the authors' own published ",
           "assignments from their K = 15 latent Dirichlet allocation, ",
           "obtained from the corresponding author and joined to the GEO ",
           "count matrices on section and barcode (39,264 of 39,264 spots ",
           "matched). Counts were summed across all nevus spots and all ",
           "melanoma spots within each patient and expressed per million, ",
           "giving one paired observation per patient (n = ",
           st$n_patients[1], "). (A) Each line is one patient. (B) The ",
           "within-patient difference. Median ",
           sprintf("%.2f", st$median_nevus_cpm[1]), " CPM in nevus against ",
           sprintf("%.2f", st$median_tumour_cpm[1]), " CPM in melanoma; ",
           "Wilcoxon signed-rank p = ", fmt_p(st$wilcoxon_p[1]), ". Because ",
           "both compartments lie in the same section they share the same ",
           "stroma and the same patient, so this comparison is not subject to ",
           "the tissue-composition confounding that affects bulk ",
           "nevus-versus-melanoma cohorts."))
  log_line("  wrote: ", stem, logfile = LOG)

  # ---- Where the gene sits across cell types --------------------------------
  sub <- lin_long |>
    dplyr::filter(gene == g, !is.na(cpm)) |>
    dplyr::mutate(melano = label %in% c(NEVUS_LAB, TUMOUR_LAB))
  if (nrow(sub)) {
    pl <- ggplot(sub, aes(x = reorder(label, cpm), y = cpm, fill = melano)) +
      geom_col(width = 0.7, colour = COL_INK, linewidth = 0.25) +
      scale_fill_manual(values = c(`TRUE` = COL_TUM, `FALSE` = COL_GREY),
                        guide = "none") +
      coord_flip(clip = "off") +
      labs(x = NULL, y = bquote(bold(.(g)~"(counts per million)"))) +
      theme_journal(grid = FALSE) +
      theme(panel.grid.major.x = element_line(linewidth = 0.25,
                                              colour = "grey92"),
            plot.margin = margin(4, 12, 4, 4, unit = "pt"))
    stem2 <- paste0("geo_", tolower(g), "_visium_lineage")
    save_fig_to(pl, stem2, FIG_GEO, width = FIG$width_single, height = 95)
    add_legend(stem2,
      paste0("Figure. ", g, " across the cell types of nevus-associated melanoma."),
      paste0("Pseudobulk expression per million within each of the fifteen ",
             "spot identities assigned by Kreuger et al. (GSE298774). The two ",
             "melanocytic compartments are highlighted. This addresses where ",
             "the gene is expressed in intact skin, as distinct from whether ",
             "it differs between the melanocytic compartments."))
    log_line("  wrote: ", stem2, logfile = LOG)
  }
}

# ==============================================================================
# 9. Interpretation
# ==============================================================================
hdr("9. What this does and does not establish")

log_line("  WHAT THIS DESIGN REMOVES, relative to bulk cohorts:", logfile = LOG)
log_line("    - tissue composition: both compartments lie in one section and",
         logfile = LOG)
log_line("      share the same stroma", logfile = LOG)
log_line("    - between-patient variation: the comparison is paired", logfile = LOG)
log_line("    - patient covariates: age, sex, Breslow and site are constant",
         logfile = LOG)
log_line("      within a patient and cancel", logfile = LOG)
log_line("    - our own judgement about which spots are which: the labels are",
         logfile = LOG)
log_line("      the authors' published assignments", logfile = LOG)
log_line("    - culture artefact: this is tissue, not a cell line", logfile = LOG)
log_line("", logfile = LOG)
log_line("  WHAT REMAINS A LIMITATION:", logfile = LOG)
log_line("    - a 55 um spot holds several cells, so a 'nevus' spot is",
         logfile = LOG)
log_line("      nevus-DOMINATED, not pure, even above the purity threshold",
         logfile = LOG)
log_line("    - their object retains an uneven fraction of spots per section",
         logfile = LOG)
log_line("      (14.9% to 99.5%), consistent with depth filtering; this costs",
         logfile = LOG)
log_line("      precision but is symmetric within each patient's pair",
         logfile = LOG)
log_line("    - all lesions are WHO class I low-CSD melanoma, median Breslow",
         logfile = LOG)
log_line("      1.05 mm, so this is nevus versus EARLY melanoma and not",
         logfile = LOG)
log_line("      melanoma generally", logfile = LOG)
log_line("    - section 4 describes where the gene sits in intact skin; it is",
         logfile = LOG)
log_line("      descriptive and not a test", logfile = LOG)

report_warnings(logfile = LOG)
write_session_info("30_visium_rassf_authorlabels")
log_line("\nFull log: ", LOG, logfile = NULL)
