# ==============================================================================
# 38b_tirosh_detection.R
#
# RASSF9 in malignant melanoma cells, tested inside Tirosh before it is tested
# across datasets.
#
# WHY THE INTERNAL TEST COMES FIRST. The plan was a cross-dataset contrast --
# Belote normal melanocytes against Tirosh malignant cells, justified by shared
# Smart-seq2 chemistry. That comparison is still worth making, but it is
# confounded by study, lab and processing, so it can only ever be suggestive.
# 38a revealed something better: Tirosh contains 65 ENDOTHELIAL cells and 61
# CAFs, and those are exactly the two non-melanocytic types the GSE298774
# Visium deconvolution flagged as RASSF9-high (endothelial 35.98 CPM,
# fibroblasts 23.81, against 2.14 and 5.50 in the two tumour compartments).
#
# So Tirosh carries its own internal positive control, on the same plates and
# the same tumours as the malignant cells, and the prediction again comes from
# a different platform before this dataset was examined:
#
#   PRE-REGISTERED PREDICTION, from GSE298774:
#     RASSF9 is detected more often in endothelial cells and CAFs than in
#     malignant melanoma cells from the same tumours.
#   AND, from the Belote result (immune cells all zero):
#     T cells, B cells, macrophages and NK cells should be LOW -- a negative
#     control with far better donor coverage than the positive one.
#
# WHAT 38a ESTABLISHED, and what it costs:
#   - values are log2(TPM/10 + 1), NOT counts, so binomial thinning is
#     IMPOSSIBLE. Depth is instead controlled by restricting every comparison
#     to cells within a common detected-gene range. This is a weaker control
#     than the thinning used on Belote and is labelled as such throughout.
#   - the first THREE data rows are annotation (tumor; malignant flag;
#     non-malignant cell type), not two. 38a's detector missed row 1 because
#     "tumor" passes a gene-symbol test, and falsely flagged TRAF3IP2-AS1.
#     The row count is therefore hard-coded and then VERIFIED against the
#     labels rather than trusted.
#   - encodings, read off the file: malignant 1=no, 2=yes, 0=unresolved;
#     cell type 1=T, 2=B, 3=Macro, 4=Endo, 5=CAF, 6=NK, 0=unassigned
#   - two gene labels are duplicated and must be collapsed on the LINEAR
#     scale, since summing log values is meaningless
#
# DONOR COVERAGE IS CHECKED BEFORE A REFERENCE IS CHOSEN. On the Belote data I
# built a headline on eccrine cells and only then discovered they came from
# three usable donors. 65 endothelial cells across 19 tumours may be just as
# concentrated. Section 3 settles that before any test is run.
#
# CRITERION, FIXED BEFORE THE NUMBERS:
#   SUPPORTED     = within the common depth band, endothelial+CAF detection is
#                   >= 3x the malignant rate, AND the donor-level paired test
#                   agrees where it is available
#   SUGGESTIVE    = the ratio holds but donor coverage cannot support a test
#   NOT SUPPORTED = the ratio does not hold
#   Separately, immune classes are expected LOW; if they are not, the
#   positive-control logic is unreliable and nothing here should be trusted.
#
# Output: results/geo/tirosh_detection_by_class.csv
#         results/geo/tirosh_donor_level.csv
#         results/geo/tirosh_vs_belote.csv
#         figures/geo/geo_rassf9_tirosh.{pdf,png,svg}
#         results/logs/38b_tirosh_detection.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required")

LOG <- file.path(PATHS$logs, "38b_tirosh_detection.log")
cat("Tirosh GSE72056 RASSF9 detection —", format(Sys.time()), "\n\n", file = LOG)
hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

T_DIR <- file.path(PATHS$raw, "geo", "tirosh")
B_DIR <- file.path(PATHS$raw, "geo", "belote")
RES   <- file.path(PATHS$results, "geo")
FIG_GEO <- file.path(PATHS$figures, "geo")
invisible(lapply(c(RES, FIG_GEO), dir.create, recursive = TRUE, showWarnings = FALSE))

F_MAT <- file.path(T_DIR, "GSE72056_melanoma_single_cell_revised_v2.txt.gz")
TARGETS   <- c("RASSF9", "RASSF10")
MEL_MK    <- c("MLANA", "PMEL", "TYR", "DCT", "TYRP1")
N_ANN     <- 3          # verified below, not assumed
MIN_CELLS <- 15
RATIO_MIN <- 3.0
SEED      <- 20261016
set.seed(SEED)

pick <- function(slot, fb) {
  if (!exists("PAL", inherits = TRUE)) return(fb)
  x <- get("PAL", inherits = TRUE)[[slot]]
  if (!is.null(x) && is.character(x) && nzchar(x[1])) x[1] else fb
}
COL_MAL <- pick("melanoma", "#F2680C"); COL_REF <- pick("nevus", "#0095D9")
COL_GREY <- pick("grey_pt", "#8E8E8E"); COL_INK <- pick("ink", "#1A1A1A")

# ==============================================================================
# 1. Load, strip annotation rows, collapse duplicate genes
# ==============================================================================
hdr("1. Load")

if (!file.exists(F_MAT)) stop("Missing ", F_MAT)
log_line("  reading ", basename(F_MAT), sprintf(" (%.1f MB gz)",
         file.size(F_MAT) / 1024^2), logfile = LOG)
DT <- data.table::fread(F_MAT, sep = "\t", showProgress = FALSE)
lab <- as.character(DT[[1]])
cells <- names(DT)[-1]
log_line("  ", length(lab), " data rows x ", length(cells), " cells",
         logfile = LOG)

# Verify the annotation rows rather than trusting the count.
want <- c("tumor", "malignant", "non-malignant")
got  <- tolower(lab[seq_len(N_ANN)])
ok <- all(vapply(seq_len(N_ANN), function(i) grepl(want[i], got[i], fixed = TRUE),
                 logical(1)))
log_line("  annotation rows: ", paste(lab[seq_len(N_ANN)], collapse = " | "),
         logfile = LOG)
if (!ok) stop("The first ", N_ANN, " rows are not the expected annotation rows.")
log_line("  verified: the first ", N_ANN, " rows are annotation, not genes",
         logfile = LOG)

tumour <- as.integer(unlist(DT[1, -1], use.names = FALSE))
malig  <- as.integer(unlist(DT[2, -1], use.names = FALSE))
ctype  <- as.integer(unlist(DT[3, -1], use.names = FALSE))

glab <- lab[-seq_len(N_ANN)]
G <- as.matrix(DT[-seq_len(N_ANN), -1, with = FALSE])
rownames(G) <- glab
rm(DT); gc(verbose = FALSE)

dups <- unique(glab[duplicated(glab)])
if (length(dups)) {
  log_line("  duplicated gene labels (", length(dups), "): ",
           paste(dups, collapse = ", "), logfile = LOG)
  log_line("  collapsing on the LINEAR scale; summing log2 values would be",
           logfile = LOG)
  log_line("  meaningless", logfile = LOG)
  hit <- intersect(dups, c(TARGETS, MEL_MK))
  if (length(hit))
    log_line("  [FLAG] a target or marker is duplicated: ",
             paste(hit, collapse = ", "), logfile = LOG)
  lin <- 2^G - 1
  lin[lin < 0] <- 0
  agg <- rowsum(lin, group = rownames(lin), reorder = FALSE)
  G <- log2(agg + 1)
  rm(lin, agg); gc(verbose = FALSE)
  log_line("  genes after collapsing: ", nrow(G), logfile = LOG)
}
for (g in TARGETS) if (!g %in% rownames(G)) stop(g, " missing after collapsing")

# ==============================================================================
# 2. Cell classes
# ==============================================================================
hdr("2. Cell classes")

CT <- c("0" = "Non-malignant, unassigned", "1" = "T cell", "2" = "B cell",
        "3" = "Macrophage", "4" = "Endothelial", "5" = "CAF", "6" = "NK")
cls <- ifelse(malig == 2, "Malignant",
       ifelse(malig == 0, "Unresolved", CT[as.character(ctype)]))
cls[is.na(cls)] <- "Unresolved"

CELL <- tibble::tibble(cell = cells, tumour = tumour, malignant = malig,
                       ctype = ctype, class = cls,
                       n_det = colSums(G > 0, na.rm = TRUE))
for (g in TARGETS) CELL[[g]] <- as.numeric(G[g, ])

tb <- sort(table(CELL$class), decreasing = TRUE)
for (i in seq_along(tb))
  log_line(sprintf("    %-28s %5d cells", names(tb)[i], as.integer(tb[i])),
           logfile = LOG)
log_line("  tumours: ", dplyr::n_distinct(CELL$tumour), logfile = LOG)
log_line("  Unresolved cells are reported but excluded from every test",
         logfile = LOG)

# ==============================================================================
# 3. Donor coverage BEFORE choosing a reference
# ==============================================================================
hdr("3. Donor coverage by class -- checked first, not after")

COV <- CELL |> dplyr::filter(class != "Unresolved") |>
  dplyr::count(class, tumour) |>
  dplyr::group_by(class) |>
  dplyr::summarise(cells = sum(n), tumours = dplyr::n_distinct(tumour),
                   tumours_ge = sum(n >= MIN_CELLS), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(cells))
for (i in seq_len(nrow(COV)))
  log_line(sprintf("    %-28s %5d cells | %2d tumours | %2d with >= %d cells",
                   COV$class[i], COV$cells[i], COV$tumours[i],
                   COV$tumours_ge[i], MIN_CELLS), logfile = LOG)

pair_n <- function(a, b) {
  x <- CELL |> dplyr::filter(class %in% c(a, b)) |>
    dplyr::count(tumour, class) |> dplyr::filter(n >= MIN_CELLS)
  length(intersect(x$tumour[x$class == a], x$tumour[x$class == b]))
}
log_line("", logfile = LOG)
for (r in c("Endothelial", "CAF", "T cell", "Macrophage"))
  log_line(sprintf("    tumours with >= %d cells in BOTH Malignant and %-12s %2d",
                   MIN_CELLS, r, pair_n("Malignant", r)), logfile = LOG)

# ==============================================================================
# 4. Marker gate
# ==============================================================================
hdr("4. Marker gate")

mean_in <- function(g, k) mean(as.numeric(G[g, CELL$class == k]), na.rm = TRUE)
gate <- c()
for (g in MEL_MK) {
  a <- mean_in(g, "Malignant"); b <- mean_in(g, "T cell")
  ok2 <- is.finite(a) && is.finite(b) && a > b
  gate <- c(gate, ok2)
  log_line(sprintf("    %-7s malignant %7.3f vs T cell %7.3f  %s", g, a, b,
                   if (ok2) "ok" else "FAIL"), logfile = LOG)
}
for (pr in list(c("PTPRC", "T cell"), c("PECAM1", "Endothelial"),
                c("COL1A1", "CAF"))) {
  if (!pr[1] %in% rownames(G)) next
  a <- mean_in(pr[1], pr[2]); b <- mean_in(pr[1], "Malignant")
  ok2 <- is.finite(a) && is.finite(b) && a > b
  gate <- c(gate, ok2)
  log_line(sprintf("    %-7s %-12s %7.3f vs malignant %7.3f  %s", pr[1], pr[2],
                   a, b, if (ok2) "ok" else "FAIL"), logfile = LOG)
}
pass <- mean(gate, na.rm = TRUE)
log_line(sprintf("  gate: %.0f%% passed", 100 * pass), logfile = LOG)
if (pass < 0.75) stop("Marker gate failed; see ", LOG)

# ==============================================================================
# 5. Depth, and the common band that replaces thinning
# ==============================================================================
hdr("5. Depth control without thinning")

log_line("  Values are log2(TPM/10+1), so counts are unavailable and binomial",
         logfile = LOG)
log_line("  thinning is impossible. Detected-gene count is the depth proxy,",
         logfile = LOG)
log_line("  and comparisons are restricted to a band covered by every class.",
         logfile = LOG)

DEP <- CELL |> dplyr::filter(class != "Unresolved") |>
  dplyr::group_by(class) |>
  dplyr::summarise(n = dplyr::n(), p10 = stats::quantile(n_det, 0.1),
                   med = stats::median(n_det), p90 = stats::quantile(n_det, 0.9),
                   .groups = "drop") |> dplyr::arrange(dplyr::desc(n))
for (i in seq_len(nrow(DEP)))
  log_line(sprintf("    %-28s n = %5d | genes detected p10 %5.0f med %5.0f p90 %5.0f",
                   DEP$class[i], DEP$n[i], DEP$p10[i], DEP$med[i], DEP$p90[i]),
           logfile = LOG)

focus <- c("Malignant", "Endothelial", "CAF", "T cell")
sub <- DEP |> dplyr::filter(class %in% focus)
band_lo <- max(sub$p10); band_hi <- min(sub$p90)
log_line(sprintf("  common band: %.0f-%.0f detected genes", band_lo, band_hi),
         logfile = LOG)
if (!is.finite(band_lo) || !is.finite(band_hi) || band_lo >= band_hi) {
  log_line("  [FLAG] the classes do not overlap in depth; the band is empty.",
           logfile = LOG)
  log_line("  Falling back to the overall interquartile range and reporting",
           logfile = LOG)
  log_line("  the raw comparison only, with depth NOT controlled.", logfile = LOG)
  band_lo <- stats::quantile(CELL$n_det, 0.25)
  band_hi <- stats::quantile(CELL$n_det, 0.75)
}
CELL$in_band <- CELL$n_det >= band_lo & CELL$n_det <= band_hi
for (k in focus)
  log_line(sprintf("    %-28s %4d of %4d cells in band", k,
                   sum(CELL$in_band & CELL$class == k),
                   sum(CELL$class == k)), logfile = LOG)

# ==============================================================================
# 6. Detection by class
# ==============================================================================
hdr("6. Detection by class")

det_by <- function(d, lbl) {
  dplyr::bind_rows(lapply(TARGETS, function(g) {
    d |> dplyr::filter(class != "Unresolved") |>
      dplyr::group_by(class) |>
      dplyr::summarise(gene = g, scale = lbl, n_cells = dplyr::n(),
                       n_det = sum(.data[[g]] > 0, na.rm = TRUE),
                       pct_det = 100 * mean(.data[[g]] > 0, na.rm = TRUE),
                       mean_log = mean(.data[[g]], na.rm = TRUE),
                       .groups = "drop")
  }))
}
ALL <- dplyr::bind_rows(det_by(CELL, "all cells"),
                        det_by(CELL[CELL$in_band, ], "depth band"))
for (lbl in c("all cells", "depth band")) {
  for (g in TARGETS) {
    s <- ALL |> dplyr::filter(gene == g, scale == lbl) |>
      dplyr::arrange(dplyr::desc(n_cells))
    log_line("\n  ", g, " (", lbl, "):", logfile = LOG)
    log_line(sprintf("    %-28s %7s %7s %8s %10s", "class", "cells", "det",
                     "% det", "mean log"), logfile = LOG)
    for (i in seq_len(nrow(s)))
      log_line(sprintf("    %-28s %7d %7d %7.2f%% %10.4f", s$class[i],
                       s$n_cells[i], s$n_det[i], s$pct_det[i], s$mean_log[i]),
               logfile = LOG)
  }
}
readr::write_csv(ALL, file.path(RES, "tirosh_detection_by_class.csv"))

# ==============================================================================
# 7. The pre-registered prediction
# ==============================================================================
hdr("7. Pre-registered prediction: endothelial and CAF > malignant")

getp <- function(g, k, lbl = "depth band") {
  v <- ALL$pct_det[ALL$gene == g & ALL$class == k & ALL$scale == lbl]
  if (length(v)) v[1] else NA_real_
}
SUM <- list()
for (g in TARGETS) {
  m <- getp(g, "Malignant")
  e <- getp(g, "Endothelial"); c2 <- getp(g, "CAF")
  tc <- getp(g, "T cell")
  pooled <- {
    d <- CELL[CELL$in_band & CELL$class %in% c("Endothelial", "CAF"), ]
    if (nrow(d)) 100 * mean(d[[g]] > 0, na.rm = TRUE) else NA_real_
  }
  log_line(sprintf("\n  %s within the depth band:", g), logfile = LOG)
  log_line(sprintf("    malignant %.2f%% | endothelial %.2f%% | CAF %.2f%% | endo+CAF %.2f%% | T cell %.2f%%",
                   m, e, c2, pooled, tc), logfile = LOG)
  ratio <- if (is.finite(m) && m > 0) pooled / m else
    if (isTRUE(pooled > 0)) Inf else NA_real_
  log_line(sprintf("    (endo+CAF) / malignant = %s",
                   if (is.finite(ratio)) sprintf("%.2fx", ratio) else
                     if (is.infinite(ratio)) "malignant rate is 0" else "n/a"),
           logfile = LOG)
  log_line(sprintf("    immune negative control, T cell = %.2f%%", tc),
           logfile = LOG)
  n_caf <- sum(CELL$in_band & CELL$class == "CAF")
  n_end <- sum(CELL$in_band & CELL$class == "Endothelial")
  log_line(sprintf("    cells behind those: endothelial %d, CAF %d -- a CAF rate",
                   n_end, n_caf), logfile = LOG)
  log_line("    from so few cells carries a wide interval and the all-cells",
           logfile = LOG)
  log_line("    figure is the more informative one for that class",
           logfile = LOG)
  SUM[[g]] <- tibble::tibble(gene = g, malignant = m, endo = e, caf = c2,
                             endo_caf = pooled, tcell = tc, ratio = ratio)
}

DON <- CELL |> dplyr::filter(class != "Unresolved") |>
  dplyr::group_by(tumour, class) |>
  dplyr::summarise(n = dplyr::n(),
                   pct9  = 100 * mean(RASSF9 > 0, na.rm = TRUE),
                   pct10 = 100 * mean(RASSF10 > 0, na.rm = TRUE),
                   .groups = "drop")
readr::write_csv(DON, file.path(RES, "tirosh_donor_level.csv"))

paired <- function(ref, col, lbl) {
  a <- as.data.frame(DON)[DON$class == "Malignant" & DON$n >= MIN_CELLS,
                          c("tumour", "n", col)]
  names(a) <- c("tumour", "n_mal", "mal")
  b <- as.data.frame(DON)[DON$class == ref & DON$n >= MIN_CELLS,
                          c("tumour", "n", col)]
  names(b) <- c("tumour", "n_ref", "ref")
  p <- dplyr::inner_join(a, b, by = "tumour")
  log_line("\n  ", lbl, " vs Malignant, reference = ", ref, logfile = LOG)
  if (!nrow(p)) {
    log_line("    no tumour has >= ", MIN_CELLS, " cells in both; no test",
             logfile = LOG); return(NULL)
  }
  for (i in seq_len(nrow(p)))
    log_line(sprintf("    tumour %-5s n mal %4d  n ref %4d | mal %6.2f%% | ref %6.2f%%",
                     p$tumour[i], p$n_mal[i], p$n_ref[i], p$mal[i], p$ref[i]),
             logfile = LOG)
  higher <- sum(p$ref > p$mal)
  log_line(sprintf("    reference higher in %d of %d tumours", higher, nrow(p)),
           logfile = LOG)
  if (nrow(p) >= 4) {
    wt <- suppressWarnings(stats::wilcox.test(p$ref, p$mal, paired = TRUE))
    log_line("    paired Wilcoxon p = ", fmt_p(wt$p.value), logfile = LOG)
    tibble::tibble(gene = lbl, reference = ref, tumours = nrow(p),
                   higher = higher, p = wt$p.value)
  } else {
    log_line("    fewer than 4 tumours; no paired test attempted", logfile = LOG)
    tibble::tibble(gene = lbl, reference = ref, tumours = nrow(p),
                   higher = higher, p = NA_real_)
  }
}
PD <- dplyr::bind_rows(
  paired("Endothelial", "pct9",  "RASSF9"),
  paired("CAF",         "pct9",  "RASSF9"),
  paired("T cell",      "pct9",  "RASSF9"),
  paired("Endothelial", "pct10", "RASSF10"),
  paired("CAF",         "pct10", "RASSF10"),
  paired("T cell",      "pct10", "RASSF10"))

# ==============================================================================
# 8. Cross-dataset, clearly subordinate
# ==============================================================================
hdr("8. Tirosh malignant vs Belote normal melanocytes")

log_line("  CONFOUNDED BY DESIGN: different studies, labs and processing.",
         logfile = LOG)
log_line("  Shared Smart-seq2 chemistry removes the largest source of",
         logfile = LOG)
log_line("  difference but not batch, and Belote is counts while Tirosh is",
         logfile = LOG)
log_line("  log TPM, so only DETECTION can be compared, never level. Depth is",
         logfile = LOG)
log_line("  matched on detected-gene count using the authors' own n_genes for",
         logfile = LOG)
log_line("  Belote. Suggestive only.", logfile = LOG)

f_bcell <- file.path(RES, "belote_cell_level.csv")
f_bmeta <- file.path(B_DIR, "GSE151091_Metadata.csv.gz")
XD <- NULL
if (file.exists(f_bcell) && file.exists(f_bmeta)) {
  BC <- readr::read_csv(f_bcell, show_col_types = FALSE, progress = FALSE)
  BM <- data.table::fread(f_bmeta, showProgress = FALSE)
  if (all(c("cell_id", "n_genes") %in% names(BM))) {
    BC <- dplyr::left_join(BC, as.data.frame(BM)[, c("cell_id", "n_genes")],
                           by = "cell_id")
    BMel <- BC |> dplyr::filter(group == "Melanocyte", !is.na(n_genes))
    TMal <- CELL |> dplyr::filter(class == "Malignant")
    lo <- max(stats::quantile(BMel$n_genes, 0.1),
              stats::quantile(TMal$n_det, 0.1))
    hi <- min(stats::quantile(BMel$n_genes, 0.9),
              stats::quantile(TMal$n_det, 0.9))
    log_line(sprintf("  Belote melanocytes: n_genes p10 %.0f med %.0f p90 %.0f",
                     stats::quantile(BMel$n_genes, 0.1),
                     stats::median(BMel$n_genes),
                     stats::quantile(BMel$n_genes, 0.9)), logfile = LOG)
    log_line(sprintf("  Tirosh malignant:   n_genes p10 %.0f med %.0f p90 %.0f",
                     stats::quantile(TMal$n_det, 0.1),
                     stats::median(TMal$n_det),
                     stats::quantile(TMal$n_det, 0.9)), logfile = LOG)
    if (is.finite(lo) && is.finite(hi) && lo < hi) {
      log_line(sprintf("  common band: %.0f-%.0f detected genes", lo, hi),
               logfile = LOG)
      Bb <- BMel |> dplyr::filter(n_genes >= lo, n_genes <= hi)
      Tb <- TMal |> dplyr::filter(n_det >= lo, n_det <= hi)
      rows <- list()
      for (g in TARGETS) {
        bk <- sum(Bb[[paste0(g, "_raw")]] > 0, na.rm = TRUE); bn <- nrow(Bb)
        tk <- sum(Tb[[g]] > 0, na.rm = TRUE); tn <- nrow(Tb)
        if (!bn || !tn) next
        bt <- stats::binom.test(bk, bn); tt <- stats::binom.test(tk, tn)
        ft <- try(stats::fisher.test(matrix(c(tk, tn - tk, bk, bn - bk),
                                            nrow = 2)), silent = TRUE)
        log_line(sprintf("\n  %s: Belote normal melanocytes %d/%d = %.3f%% [%.3f, %.3f]",
                         g, bk, bn, 100 * bk / bn, 100 * bt$conf.int[1],
                         100 * bt$conf.int[2]), logfile = LOG)
        log_line(sprintf("     Tirosh malignant melanoma  %d/%d = %.3f%% [%.3f, %.3f]",
                         tk, tn, 100 * tk / tn, 100 * tt$conf.int[1],
                         100 * tt$conf.int[2]), logfile = LOG)
        if (!inherits(ft, "try-error"))
          log_line("     Fisher p = ", fmt_p(ft$p.value),
                   " (batch-confounded, pseudo-replicated)", logfile = LOG)
        rows[[length(rows) + 1L]] <- tibble::tibble(
          gene = g, belote_det = bk, belote_n = bn, belote_pct = 100 * bk / bn,
          tirosh_det = tk, tirosh_n = tn, tirosh_pct = 100 * tk / tn,
          band_lo = lo, band_hi = hi,
          fisher_p = if (inherits(ft, "try-error")) NA_real_ else ft$p.value)
      }
      XD <- dplyr::bind_rows(rows)
      if (!is.null(XD) && nrow(XD))
        readr::write_csv(XD, file.path(RES, "tirosh_vs_belote.csv"))
    } else {
      log_line("  [FLAG] the two datasets do not overlap in detected-gene",
               logfile = LOG)
      log_line("  count, so no honest depth-matched comparison is possible.",
               logfile = LOG)
    }
  } else {
    log_line("  Belote metadata lacks cell_id/n_genes; skipped", logfile = LOG)
  }
} else {
  log_line("  belote_cell_level.csv or the Belote metadata is missing;",
           logfile = LOG)
  log_line("  run 37c first if the cross-dataset panel is wanted", logfile = LOG)
}

# ==============================================================================
# 9. Verdict
# ==============================================================================
hdr("9. Verdict against the pre-set criterion")

log_line("  SUPPORTED     = endo+CAF >= ", RATIO_MIN,
         "x malignant in the depth band AND", logfile = LOG)
log_line("                  the donor-level test agrees where available",
         logfile = LOG)
log_line("  SUGGESTIVE    = ratio holds but donor coverage cannot support a test",
         logfile = LOG)
log_line("  NOT SUPPORTED = the ratio does not hold", logfile = LOG)
log_line("", logfile = LOG)

for (g in TARGETS) {
  s <- SUM[[g]]; if (is.null(s)) next
  pd <- PD |> dplyr::filter(gene == g, reference %in% c("Endothelial", "CAF"),
                            is.finite(p))
  agree <- if (nrow(pd)) all(pd$p < 0.05) else NA
  verdict <- if (isTRUE(s$ratio >= RATIO_MIN) && isTRUE(agree)) "SUPPORTED"
  else if (isTRUE(s$ratio >= RATIO_MIN)) "SUGGESTIVE (no donor-level test available)"
  else "NOT SUPPORTED"
  log_line("  ", g, ": ", verdict, logfile = LOG)
  log_line(sprintf("    malignant %.2f%% | endo+CAF %.2f%% | ratio %s | immune %.2f%%",
                   s$malignant, s$endo_caf,
                   if (is.finite(s$ratio)) sprintf("%.2fx", s$ratio) else "n/a",
                   s$tcell), logfile = LOG)
}
log_line("", logfile = LOG)
log_line("  If the immune classes are NOT low, the positive-control logic is",
         logfile = LOG)
log_line("  unreliable and none of this should be trusted.", logfile = LOG)

# ==============================================================================
# 10. Figure
# ==============================================================================
hdr("10. Figure")

reset_legends("# Tirosh GSE72056 figure legends\n")

pd <- ALL |> dplyr::filter(gene == "RASSF9", scale == "depth band",
                           class %in% c("Malignant", "Endothelial", "CAF",
                                        "T cell", "B cell", "Macrophage", "NK"))
pa <- ggplot(pd, aes(x = reorder(class, pct_det), y = pct_det,
                     fill = class == "Malignant")) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c(`TRUE` = COL_MAL, `FALSE` = COL_REF),
                    guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "Cells with RASSF9 detected (%)") +
  theme_journal() +
  theme(plot.margin = margin(4, 12, 4, 4, unit = "pt"))

plots <- list(pa)
if (!is.null(XD) && nrow(XD)) {
  # Built directly rather than with transmute(): the source is one row and the
  # target is two, which transmute() cannot do -- it recycles to size 1.
  x1 <- XD[XD$gene == "RASSF9", ]
  xd <- tibble::tibble(
    set = factor(c("Normal melanocytes\n(Belote)",
                   "Malignant melanoma\n(Tirosh)"),
                 levels = c("Normal melanocytes\n(Belote)",
                            "Malignant melanoma\n(Tirosh)")),
    pct = c(x1$belote_pct[1], x1$tirosh_pct[1]),
    lab = c(sprintf("%d/%d", x1$belote_det[1], x1$belote_n[1]),
            sprintf("%d/%d", x1$tirosh_det[1], x1$tirosh_n[1])))
  pb <- ggplot(xd, aes(x = set, y = pct, fill = set)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = lab), vjust = -0.6,
              size = (FIG$base_size - 3) / .pt, colour = COL_INK) +
    scale_fill_manual(values = c(COL_REF, COL_MAL), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "RASSF9 detected (%)") +
    theme_journal() +
    theme(plot.margin = margin(4, 12, 4, 4, unit = "pt"))
  plots <- list(pa, pb)
}
fig <- patchwork::wrap_plots(plots, nrow = 1) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = FIG$base_size + 2, face = "bold",
                                colour = COL_INK))
save_fig_to(fig, "geo_rassf9_tirosh", FIG_GEO,
            width = FIG$width_double, height = 86)
log_line("  wrote: geo_rassf9_tirosh", logfile = LOG)

add_legend("geo_rassf9_tirosh",
  "Figure. RASSF9 detection in malignant melanoma cells and their stromal neighbours.",
  paste0("(A) Percentage of cells with RASSF9 detected, by author-assigned class, ",
         "within a common detected-gene band of ", round(band_lo), "-",
         round(band_hi), " genes. Because this matrix is distributed as ",
         "log2(TPM/10+1) rather than counts, binomial thinning is not possible ",
         "and detected-gene count is used as the depth proxy -- a weaker control. ",
         "Endothelial cells and CAFs serve as a positive control, predicted from ",
         "an independent Visium deconvolution; immune classes serve as a negative ",
         "control. ",
         if (!is.null(XD) && nrow(XD))
           paste0("(B) Detection in normal melanocytes (Belote, GSE151091) against ",
                  "malignant melanoma cells (Tirosh, GSE72056), restricted to a ",
                  "shared detected-gene range. This contrast is confounded by study ",
                  "and processing and is suggestive only. ") else "",
         "Data from Tirosh et al., Science 2016 (GSE72056), ",
         dplyr::n_distinct(CELL$tumour), " tumours."))

report_warnings(logfile = LOG)
write_session_info("38b_tirosh_detection")
log_line("\nFull log: ", LOG, logfile = NULL)
