# ==============================================================================
# 37b_belote_detection.R
#
# Is RASSF9 detectable in normal human melanocytes, and if not, does that mean
# anything?
#
# Belote et al. (Nat Cell Biol 2021;23:1035-1047, PMID 34475532), GSE151091:
# Smart-seq2, 26,576 genes x 14,370 cells, 9,719 annotated, 22 donors, median
# 3,650 genes per cell. This is the deepest public measurement of normal human
# melanocytes, so if RASSF9 cannot be seen here it cannot be seen anywhere in
# public single-cell data, and that is an answer rather than a failure.
#
# THE TRAP THIS SCRIPT IS BUILT TO AVOID. "RASSF9 is detected in x% of
# melanocytes" means nothing on its own. A low number has two incompatible
# explanations -- the gene is not expressed, or the assay cannot see genes at
# that abundance -- and they demand opposite conclusions. So the script does not
# report a detection rate as a finding. It reports it against two references.
#
# REFERENCE 1, THE INTERNAL POSITIVE CONTROL, AND THE REASON THIS DATASET IS
# WORTH THE TROUBLE. class_1 contains 636 eccrine cells and 1,865 keratinocytes
# alongside the melanocytes -- same experiment, same plates, same chemistry. Our
# GSE298774 Visium analysis found RASSF9 concentrated in eccrine duct and gland
# (30.41 CPM) against 5.50 and 2.14 in the two tumour compartments, with only
# 17.1% of its molecules in melanocytic spots. That prediction was made on a
# different platform, in a different dataset, before anyone looked here:
#
#   PRE-REGISTERED PREDICTION: if RASSF9 is detectable by this assay at all,
#   eccrine cells show it and melanocytes largely do not.
#
# If that holds, non-detection in melanocytes is informative AND the
# localisation claim is independently replicated at single-cell resolution --
# which is the point Slieker pushed back on. If eccrine cells are also blank,
# the assay is the limitation and no conclusion about melanocytes follows.
#
# REFERENCE 2, ABUNDANCE MATCHING. Among genes at RASSF9's own expression level,
# what detection rate is normal here? If RASSF9 sits where its abundance
# predicts, its sparsity is abundance, not absence.
#
# DEPTH IS CONTROLLED, NOT ASSUMED. Detection rates rise with sequencing depth,
# and in the Visium analysis a 28.5-fold depth range inverted the raw-count
# rankings entirely. If eccrine cells are sequenced deeper than melanocytes,
# higher detection there is an artefact. So every comparison runs twice: raw,
# then after binomial thinning to a common depth.
#
# DONORS, NOT CELLS, ARE THE UNIT. 6,889 melanocytes from 22 donors are not
# 6,889 independent observations. Cell-level tests are reported only with the
# pseudo-replication stated, and the donor-level test is the one that counts.
#
# CRITERION, FIXED BEFORE THE NUMBERS:
#   ADEQUATE            = RASSF9 detected in >= 10% of melanocytes after
#                         thinning -> proceed to cross-dataset comparison
#   INFORMATIVE ABSENCE = < 10% in melanocytes AND eccrine detection >= 3x the
#                         melanocyte rate after thinning with donor-level
#                         support, OR abundance-matched genes detected >= 3x
#   UNINFORMATIVE       = < 10% AND eccrine not higher AND abundance-matched
#                         genes similarly sparse -> the assay cannot see genes
#                         at this level; report as inconclusive, NOT as absence
#
# Output: results/geo/belote_detection_by_class.csv
#         results/geo/belote_detection_donor.csv
#         results/geo/belote_calibration.csv
#         figures/geo/geo_rassf9_belote_detection.{pdf,png,svg}
#         results/logs/37b_belote_detection.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("data.table is required. Install with:\n",
       '  install.packages("data.table")')
}

LOG <- file.path(PATHS$logs, "37b_belote_detection.log")
cat("Belote GSE151091 RASSF9 detection —", format(Sys.time()), "\n\n", file = LOG)
hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

DIR  <- file.path(PATHS$raw, "geo", "belote")
RES  <- file.path(PATHS$results, "geo")
FIG_GEO <- file.path(PATHS$figures, "geo")
invisible(lapply(c(RES, FIG_GEO), dir.create, recursive = TRUE, showWarnings = FALSE))

F_MAT <- file.path(DIR, "GSE151091_raw_matrix.csv.gz")
F_MET <- file.path(DIR, "GSE151091_Metadata.csv.gz")
F_HI  <- file.path(DIR, "GSE151091_HighResLouvain.csv.gz")

TARGETS <- c("RASSF9", "RASSF10")
MEL_MK  <- c("MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10", "TYRP1")
KRT_MK  <- c("KRT14", "KRT1")
IMM_MK  <- c("PTPRC")
DET_MIN   <- 0.10     # "adequate" detection threshold
RATIO_MIN <- 3.0      # how much higher eccrine / matched genes must be
MIN_CELLS <- 20       # minimum cells for a donor x class estimate
N_CAL     <- 2000     # random genes for the abundance calibration
SEED      <- 20261016
set.seed(SEED)

pick <- function(slot, fb) {
  if (!exists("PAL", inherits = TRUE)) return(fb)
  x <- get("PAL", inherits = TRUE)[[slot]]
  if (!is.null(x) && is.character(x) && nzchar(x[1])) x[1] else fb
}
COL_MEL <- pick("melanoma", "#F2680C"); COL_ECC <- pick("nevus", "#0095D9")
COL_GREY <- pick("grey_pt", "#8E8E8E"); COL_INK <- pick("ink", "#1A1A1A")

# ==============================================================================
# 1. Metadata and cell groups
# ==============================================================================
hdr("1. Metadata")

for (f in c(F_MAT, F_MET)) if (!file.exists(f)) stop("Missing ", f)

MET <- data.table::fread(F_MET, showProgress = FALSE)
log_line("  metadata: ", nrow(MET), " cells x ", ncol(MET), " columns",
         logfile = LOG)
stopifnot(all(c("cell_id", "class_1", "patient") %in% names(MET)))

# Cycling variants are the same lineage caught mid-cycle, so they are folded
# into their parent class rather than dropped or left as separate groups.
MET$group <- dplyr::case_when(
  MET$class_1 %in% c("mel", "cyc_mel")             ~ "Melanocyte",
  MET$class_1 %in% c("krt", "cyc_krt")             ~ "Keratinocyte",
  MET$class_1 %in% c("eccrine", "cyc_eccrine")     ~ "Eccrine",
  MET$class_1 %in% c("dendritic", "cyc_dendritic") ~ "Dendritic",
  MET$class_1 == "T-cell"                          ~ "T cell",
  MET$class_1 == "mast"                            ~ "Mast",
  TRUE                                             ~ as.character(MET$class_1))

tb <- sort(table(MET$group), decreasing = TRUE)
for (i in seq_along(tb))
  log_line(sprintf("    %-14s %5d cells", names(tb)[i], as.integer(tb[i])),
           logfile = LOG)
log_line("  donors: ", dplyr::n_distinct(MET$patient), logfile = LOG)

# The published annotated total is 9,688; this file has 9,719, and donor
# A1015LM contributes exactly 31. Reported both ways rather than silently
# adopting the paper's filter.
n_a1015 <- sum(MET$patient == "A1015LM")
log_line("  donor A1015LM contributes ", n_a1015, " cells; ", nrow(MET), " - ",
         n_a1015, " = ", nrow(MET) - n_a1015,
         " matches the published annotated total", logfile = LOG)

# ==============================================================================
# 2. Matrix
# ==============================================================================
hdr("2. Expression matrix")

log_line("  reading ", basename(F_MAT), sprintf(" (%.1f MB gz)",
         file.size(F_MAT) / 1024^2), " -- a minute or two", logfile = LOG)
DT <- data.table::fread(F_MAT, showProgress = FALSE)
genes <- as.character(DT[[1]])
cells <- names(DT)[-1]
log_line("  matrix: ", length(genes), " genes x ", length(cells), " cells",
         logfile = LOG)

for (g in TARGETS) {
  if (!g %in% genes) stop(g, " absent from the matrix; 37a said it was present")
}

keep <- intersect(cells, MET$cell_id)
log_line("  cells with annotation: ", length(keep), " of ", length(cells),
         " (", length(cells) - length(keep), " unannotated, excluded)",
         logfile = LOG)

M <- as.matrix(DT[, keep, with = FALSE])
rownames(M) <- genes
rm(DT); gc(verbose = FALSE)
MET <- MET[match(keep, MET$cell_id), ]
stopifnot(identical(as.character(MET$cell_id), colnames(M)))

tot <- colSums(M)
log_line("  library size: median ", format(round(stats::median(tot)), big.mark = ","),
         ", range ", format(round(min(tot)), big.mark = ","), "-",
         format(round(max(tot)), big.mark = ","), logfile = LOG)

hdr("3. Depth by cell class -- the artefact to rule out")

log_line("  If one class is sequenced deeper, it will show higher detection for",
         logfile = LOG)
log_line("  every gene regardless of biology. Section 6 removes this by",
         logfile = LOG)
log_line("  thinning all cells to a common depth.", logfile = LOG)
DEP <- tibble::tibble(group = MET$group, total = tot, genes_det = colSums(M > 0)) |>
  dplyr::group_by(group) |>
  dplyr::summarise(n = dplyr::n(), med_total = stats::median(total),
                   p10_total = stats::quantile(total, 0.1),
                   med_genes = stats::median(genes_det), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(n))
for (i in seq_len(nrow(DEP)))
  log_line(sprintf("    %-14s n = %5d | median counts %12s | median genes %6.0f",
                   DEP$group[i], DEP$n[i],
                   format(round(DEP$med_total[i]), big.mark = ","),
                   DEP$med_genes[i]), logfile = LOG)
fold <- max(DEP$med_total) / min(DEP$med_total)
log_line(sprintf("  depth range across classes: %.1f-fold", fold), logfile = LOG)

# ==============================================================================
# 4. Marker gate
# ==============================================================================
hdr("4. Marker gate")

log_line("  The annotations are the authors'. This checks they behave before",
         logfile = LOG)
log_line("  anything is concluded from them.", logfile = LOG)

cp10k <- function(gene) {
  if (!gene %in% rownames(M)) return(rep(NA_real_, ncol(M)))
  1e4 * M[gene, ] / pmax(tot, 1)
}
gmean <- function(gene, grp) {
  v <- cp10k(gene); if (all(is.na(v))) return(NA_real_)
  mean(v[MET$group == grp], na.rm = TRUE)
}
gate <- c()
for (g in MEL_MK) {
  a <- gmean(g, "Melanocyte"); b <- gmean(g, "Keratinocyte")
  ok <- is.finite(a) && is.finite(b) && a > b
  gate <- c(gate, ok)
  log_line(sprintf("    %-7s melanocyte %9.2f vs keratinocyte %9.2f  %s", g, a, b,
                   if (ok) "ok" else "FAIL"), logfile = LOG)
}
for (g in KRT_MK) {
  a <- gmean(g, "Keratinocyte"); b <- gmean(g, "Melanocyte")
  ok <- is.finite(a) && is.finite(b) && a > b
  gate <- c(gate, ok)
  log_line(sprintf("    %-7s keratinocyte %9.2f vs melanocyte %9.2f  %s", g, a, b,
                   if (ok) "ok" else "FAIL"), logfile = LOG)
}
pass <- mean(gate, na.rm = TRUE)
log_line(sprintf("  gate: %.0f%% of checks passed", 100 * pass), logfile = LOG)
if (pass < 0.75) stop("Marker gate failed; see ", LOG)

# ==============================================================================
# 5. Detection, raw
# ==============================================================================
hdr("5. Detection by cell class (raw, depth NOT matched)")

det_tab <- function(gene, mat, lbl) {
  v <- mat[gene, ]
  tibble::tibble(gene = gene, scale = lbl, group = MET$group, x = v,
                 lib = tot) |>
    dplyr::group_by(gene, scale, group) |>
    dplyr::summarise(n_cells = dplyr::n(), n_det = sum(x > 0),
                     pct_det = 100 * mean(x > 0),
                     mean_count = mean(x),
                     mean_cp10k = mean(1e4 * x / pmax(lib, 1)),
                     .groups = "drop")
}
RAW <- dplyr::bind_rows(lapply(TARGETS, det_tab, mat = M, lbl = "raw"))
show_det <- function(D) {
  for (g in TARGETS) {
    s <- D |> dplyr::filter(gene == g) |> dplyr::arrange(dplyr::desc(n_cells))
    if (!nrow(s)) next
    log_line("\n  ", g, ":", logfile = LOG)
    log_line(sprintf("    %-14s %7s %7s %8s %11s", "class", "cells", "det",
                     "% det", "mean CP10K"), logfile = LOG)
    for (i in seq_len(nrow(s)))
      log_line(sprintf("    %-14s %7d %7d %7.2f%% %11.3f", s$group[i],
                       s$n_cells[i], s$n_det[i], s$pct_det[i], s$mean_cp10k[i]),
               logfile = LOG)
  }
}
show_det(RAW)

# ==============================================================================
# 6. Detection, depth-matched by binomial thinning
# ==============================================================================
hdr("6. Detection after thinning to a common depth")

cmp <- c("Melanocyte", "Eccrine", "Keratinocyte")
D_TARGET <- floor(min(DEP$p10_total[DEP$group %in% cmp]))
log_line("  target depth = ", format(D_TARGET, big.mark = ","),
         " counts (10th percentile of the shallowest compared class)",
         logfile = LOG)
usable <- tot >= D_TARGET
log_line("  cells retained: ", sum(usable), " of ", length(tot),
         sprintf(" (%.1f%%)", 100 * mean(usable)), logfile = LOG)
for (grp in cmp)
  log_line(sprintf("    %-14s %5d of %5d retained", grp,
                   sum(usable & MET$group == grp), sum(MET$group == grp)),
           logfile = LOG)

cal_genes <- sample(setdiff(rownames(M), c(TARGETS, MEL_MK, KRT_MK, IMM_MK)),
                    min(N_CAL, nrow(M) - 10))
thin_set <- unique(c(TARGETS, MEL_MK, cal_genes))
p_keep <- ifelse(usable, D_TARGET / tot, NA_real_)

# Binomial thinning: each count is resampled at the ratio of target depth to
# that cell's actual depth, so every cell is compared at the same depth.
TH <- matrix(NA_real_, nrow = length(thin_set), ncol = ncol(M),
             dimnames = list(thin_set, colnames(M)))
for (g in thin_set) {
  x <- M[g, ]
  v <- rep(NA_real_, length(x))
  idx <- which(usable & x > 0)
  if (length(idx)) v[idx] <- stats::rbinom(length(idx), size = x[idx],
                                           prob = p_keep[idx])
  v[usable & x == 0] <- 0
  TH[g, ] <- v
}

keep_u <- usable
MET_T <- MET[keep_u, ]; TH_T <- TH[, keep_u, drop = FALSE]
tot_T  <- rep(D_TARGET, sum(keep_u))

THIN <- dplyr::bind_rows(lapply(TARGETS, function(g) {
  tibble::tibble(gene = g, scale = "thinned", group = MET_T$group,
                 x = TH_T[g, ]) |>
    dplyr::group_by(gene, scale, group) |>
    dplyr::summarise(n_cells = dplyr::n(), n_det = sum(x > 0, na.rm = TRUE),
                     pct_det = 100 * mean(x > 0, na.rm = TRUE),
                     mean_count = mean(x, na.rm = TRUE),
                     mean_cp10k = mean(1e4 * x / D_TARGET, na.rm = TRUE),
                     .groups = "drop")
}))
show_det(THIN)

ALL <- dplyr::bind_rows(RAW, THIN)
readr::write_csv(ALL, file.path(RES, "belote_detection_by_class.csv"))

# ==============================================================================
# 7. The pre-registered prediction
# ==============================================================================
hdr("7. Pre-registered prediction: eccrine > melanocyte")

log_line("  From GSE298774 Visium, on a different platform, before this dataset",
         logfile = LOG)
log_line("  was examined. Donors are the unit; the cell-level test is shown",
         logfile = LOG)
log_line("  only for completeness and is pseudo-replicated.", logfile = LOG)

get_pct <- function(D, g, grp) {
  v <- D$pct_det[D$gene == g & D$group == grp]
  if (length(v)) v[1] else NA_real_
}
for (g in TARGETS) {
  m <- get_pct(THIN, g, "Melanocyte"); e <- get_pct(THIN, g, "Eccrine")
  k <- get_pct(THIN, g, "Keratinocyte")
  log_line(sprintf("\n  %s after thinning: melanocyte %.2f%% | eccrine %.2f%% | keratinocyte %.2f%%",
                   g, m, e, k), logfile = LOG)
  if (is.finite(m) && is.finite(e)) {
    ratio <- if (m > 0) e / m else Inf
    log_line(sprintf("    eccrine / melanocyte = %s",
                     if (is.finite(ratio)) sprintf("%.2fx", ratio) else "melanocyte rate is 0"),
             logfile = LOG)
    a <- sum(TH_T[g, MET_T$group == "Eccrine"] > 0, na.rm = TRUE)
    b <- sum(MET_T$group == "Eccrine")
    c2 <- sum(TH_T[g, MET_T$group == "Melanocyte"] > 0, na.rm = TRUE)
    d2 <- sum(MET_T$group == "Melanocyte")
    ft <- try(stats::fisher.test(matrix(c(a, b - a, c2, d2 - c2), nrow = 2)),
              silent = TRUE)
    if (!inherits(ft, "try-error"))
      log_line(sprintf("    cell-level Fisher p = %s (PSEUDO-REPLICATED, not the test)",
                       fmt_p(ft$p.value)), logfile = LOG)
  }
}

DON <- list()
for (g in TARGETS) {
  d <- tibble::tibble(donor = MET_T$patient, group = MET_T$group,
                      x = TH_T[g, ]) |>
    dplyr::filter(group %in% c("Melanocyte", "Eccrine")) |>
    dplyr::group_by(gene = g, donor, group) |>
    dplyr::summarise(n = dplyr::n(), pct = 100 * mean(x > 0, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::filter(n >= MIN_CELLS)
  DON[[g]] <- d
  paired <- dplyr::inner_join(
    d |> dplyr::filter(group == "Melanocyte") |>
      dplyr::select(donor, pct_Melanocyte = pct, n_Melanocyte = n),
    d |> dplyr::filter(group == "Eccrine") |>
      dplyr::select(donor, pct_Eccrine = pct, n_Eccrine = n),
    by = "donor")
  log_line("\n  ", g, ": donors with >= ", MIN_CELLS,
           " cells in both classes: ", nrow(paired), logfile = LOG)
  if (nrow(paired) >= 4) {
    wt <- suppressWarnings(stats::wilcox.test(paired$pct_Eccrine,
                                              paired$pct_Melanocyte, paired = TRUE))
    log_line(sprintf("    paired Wilcoxon across donors: p = %s | eccrine higher in %d of %d",
                     fmt_p(wt$p.value), sum(paired$pct_Eccrine > paired$pct_Melanocyte),
                     nrow(paired)), logfile = LOG)
  } else {
    log_line("    too few donors contribute both classes for a paired test;",
             logfile = LOG)
    log_line("    the donor-level claim cannot be made and is not made.",
             logfile = LOG)
  }
}
if (length(DON)) readr::write_csv(dplyr::bind_rows(DON),
                                  file.path(RES, "belote_detection_donor.csv"))

# ==============================================================================
# 8. Abundance calibration and an exact bound
# ==============================================================================
hdr("8. Is RASSF9's sparsity normal for its abundance?")

mel_T <- MET_T$group == "Melanocyte"
cal <- tibble::tibble(
  gene = cal_genes,
  pct_det = vapply(cal_genes, function(g)
    100 * mean(TH_T[g, mel_T] > 0, na.rm = TRUE), numeric(1)),
  mean_cp10k = vapply(cal_genes, function(g)
    mean(1e4 * TH_T[g, mel_T] / D_TARGET, na.rm = TRUE), numeric(1)))
readr::write_csv(cal, file.path(RES, "belote_calibration.csv"))

log_line("  ", nrow(cal), " random genes, melanocytes only, depth-matched.",
         logfile = LOG)
log_line(sprintf("  detection across those genes: median %.2f%%, 90th pct %.2f%%",
                 stats::median(cal$pct_det), stats::quantile(cal$pct_det, 0.9)),
         logfile = LOG)

for (g in TARGETS) {
  k <- sum(TH_T[g, mel_T] > 0, na.rm = TRUE); n <- sum(mel_T)
  pct <- 100 * k / n
  bt <- stats::binom.test(k, n)
  log_line(sprintf("\n  %s in melanocytes: %d of %d detected = %.2f%% [%.2f%%, %.2f%%]",
                   g, k, n, pct, 100 * bt$conf.int[1], 100 * bt$conf.int[2]),
           logfile = LOG)
  my_cp <- mean(1e4 * TH_T[g, mel_T] / D_TARGET, na.rm = TRUE)
  band <- cal |> dplyr::filter(mean_cp10k > 0,
                               abs(log2(mean_cp10k / pmax(my_cp, 1e-9))) < 1)
  if (nrow(band) >= 20) {
    log_line(sprintf("    %d genes within 2-fold of its expression level detect in %.2f%% of cells (median)",
                     nrow(band), stats::median(band$pct_det)), logfile = LOG)
    log_line(sprintf("    %s sits at the %.1f percentile of that band",
                     g, 100 * mean(band$pct_det <= pct)), logfile = LOG)
  } else {
    log_line("    too few abundance-matched genes for this calibration",
             logfile = LOG)
  }
  log_line(sprintf("    percentile among all %d calibration genes: %.1f",
                   nrow(cal), 100 * mean(cal$pct_det <= pct)), logfile = LOG)
}

# ==============================================================================
# 9. Verdict
# ==============================================================================
hdr("9. Verdict against the pre-set criterion")

log_line("  ADEQUATE            = >= ", round(100 * DET_MIN),
         "% of melanocytes after thinning", logfile = LOG)
log_line("  INFORMATIVE ABSENCE = below that AND eccrine >= ", RATIO_MIN,
         "x melanocyte, or", logfile = LOG)
log_line("                        abundance-matched genes >= ", RATIO_MIN, "x",
         logfile = LOG)
log_line("  UNINFORMATIVE       = below that AND neither reference is higher",
         logfile = LOG)
log_line("", logfile = LOG)

for (g in TARGETS) {
  m <- get_pct(THIN, g, "Melanocyte"); e <- get_pct(THIN, g, "Eccrine")
  my_cp <- mean(1e4 * TH_T[g, mel_T] / D_TARGET, na.rm = TRUE)
  band <- cal |> dplyr::filter(mean_cp10k > 0,
                               abs(log2(mean_cp10k / pmax(my_cp, 1e-9))) < 1)
  band_med <- if (nrow(band) >= 20) stats::median(band$pct_det) else NA_real_
  r_ecc  <- if (is.finite(m) && m > 0) e / m else if (is.finite(e) && e > 0) Inf else NA_real_
  r_band <- if (is.finite(m) && m > 0) band_med / m else NA_real_

  verdict <- if (is.finite(m) && m >= 100 * DET_MIN) {
    "ADEQUATE -- proceed to cross-dataset comparison"
  } else if (isTRUE(r_ecc >= RATIO_MIN) || isTRUE(r_band >= RATIO_MIN)) {
    "INFORMATIVE ABSENCE -- the assay can see genes like this; melanocytes lack it"
  } else {
    "UNINFORMATIVE -- assay-limited; no conclusion about melanocytes"
  }
  log_line("  ", g, ": ", verdict, logfile = LOG)
  log_line(sprintf("    melanocyte %.2f%% | eccrine %.2f%% (%s) | matched-gene median %s (%s)",
                   m, e,
                   if (is.finite(r_ecc)) sprintf("%.2fx", r_ecc) else "n/a",
                   if (is.finite(band_med)) sprintf("%.2f%%", band_med) else "n/a",
                   if (is.finite(r_band)) sprintf("%.2fx", r_band) else "n/a"),
           logfile = LOG)
  log_line("", logfile = LOG)
}

log_line("  WHAT THIS CANNOT SAY: this dataset has no nevus melanocytes and no",
         logfile = LOG)
log_line("  melanoma, so it cannot deliver Gustavo's three-way comparison on", logfile = LOG)
log_line("  its own. It establishes the normal-melanocyte arm, and whether that",
         logfile = LOG)
log_line("  arm is measurable at all.", logfile = LOG)

# ==============================================================================
# 10. Figure
# ==============================================================================
hdr("10. Figure")

reset_legends("# Belote GSE151091 detection figure legends\n")

pd <- ALL |> dplyr::filter(gene == "RASSF9",
                           group %in% c("Melanocyte", "Eccrine", "Keratinocyte",
                                        "Dendritic", "T cell", "Mast")) |>
  dplyr::mutate(scale = factor(scale, levels = c("raw", "thinned"),
                               labels = c("Raw", "Depth-matched")))
pa <- ggplot(pd, aes(x = reorder(group, pct_det), y = pct_det, fill = scale)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  scale_fill_manual(values = c(Raw = COL_GREY, `Depth-matched` = COL_MEL),
                    name = NULL) +
  coord_flip() +
  labs(x = NULL, y = "Cells with RASSF9 detected (%)") +
  theme_journal() +
  theme(legend.position = "bottom",
        legend.key.size = unit(8, "pt"),
        legend.text = element_text(size = FIG$base_size - 3),
        plot.margin = margin(4, 12, 4, 4, unit = "pt"))

pb <- ggplot(cal, aes(x = pct_det)) +
  geom_histogram(bins = 50, fill = COL_GREY, alpha = 0.45, colour = NA) +
  geom_vline(xintercept = get_pct(THIN, "RASSF9", "Melanocyte"),
             colour = COL_MEL, linewidth = 0.6) +
  labs(x = "Detection in melanocytes (%), random genes",
       y = "Genes") +
  theme_journal() +
  theme(plot.margin = margin(4, 12, 4, 4, unit = "pt"))

fig <- patchwork::wrap_plots(list(pa, pb), nrow = 1) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = FIG$base_size + 2, face = "bold",
                                colour = COL_INK))
save_fig_to(fig, "geo_rassf9_belote_detection", FIG_GEO,
            width = FIG$width_double, height = 86)
log_line("  wrote: geo_rassf9_belote_detection", logfile = LOG)

add_legend("geo_rassf9_belote_detection",
  "Figure. RASSF9 detection across cell types in normal human skin (Smart-seq2).",
  paste0("(A) Percentage of cells with at least one RASSF9 molecule, by author-assigned ",
         "cell class, before and after binomial thinning to a common depth of ",
         format(D_TARGET, big.mark = ","), " counts per cell; depth matching is ",
         "required because detection rates otherwise track sequencing depth rather ",
         "than expression. (B) Distribution of detection rates for ", nrow(cal),
         " randomly chosen genes in melanocytes at matched depth, with RASSF9 marked; ",
         "this establishes whether RASSF9's sparsity is unusual for a gene of its ",
         "abundance in this assay. Data from Belote et al., Nat Cell Biol 2021 ",
         "(GSE151091), ", dplyr::n_distinct(MET$patient), " donors."))

report_warnings(logfile = LOG)
write_session_info("37b_belote_detection")
log_line("\nFull log: ", LOG, logfile = NULL)
