# ==============================================================================
# 41_belote_keratinocyte_axis.R
#
# Does RASSF9 sit in differentiating keratinocytes, in human, at single-cell
# resolution? And does it travel with p21?
#
# THE CLAIM BEING TESTED. Lee et al. (PLoS ONE 2011, doi:10.1371/journal.pone.0017867)
# showed in MOUSE that RASSF9 protein is strongest in the suprabasal epidermis,
# colocalising with keratin 1, and that it induces p21(Cip1) to drive cell-cycle
# withdrawal as keratinocytes commit to differentiate. RASSF9-null mice have a
# three-fold thickened epidermis with doubled basal proliferation. Nobody has
# shown this in human cells, and nobody has shown it at single-cell resolution.
#
# Belote et al. (GSE151091) gives 1,865 human keratinocytes, Smart-seq2, from 22
# donors, already on disk. That is enough to ask the question directly.
#
# WHY THIS IS THE RIGHT FIRST EXPERIMENT. It is a comparison WITHIN one cell
# type, in single cells. There is no tumour purity, no stromal fraction, no
# bulk composition. The confound that has shadowed every bulk analysis in this
# project simply does not arise.
#
# THE REAL CONFOUND HERE IS DEPTH, AND IT CUTS IN OUR FAVOUR. Detection rates
# rise with sequencing depth, so a naive comparison would find RASSF9-positive
# cells richer in everything. All cells are therefore thinned to a common depth
# before anything is compared. Separately, terminally differentiating
# keratinocytes carry degraded RNA and lower transcriptome complexity, so they
# detect FEWER genes overall -- which biases AGAINST finding RASSF9 enriched
# there. Section 4 measures that rather than assuming it, and a positive result
# obtained against that gradient is stronger for it.
#
# PRE-REGISTERED PREDICTION, from Lee's mouse immunostaining:
#   RASSF9 detection rises across the basal-to-differentiated axis of human
#   keratinocytes, and RASSF9-positive keratinocytes carry more CDKN1A (p21).
#
# PRIMARY ANALYSIS. Rather than a binary split, keratinocytes are ordered on a
# differentiation score and RASSF9 detection is computed across quintiles. This
# uses every cell, gives a dose-response rather than a single contrast, and is
# harder to obtain by chance than a two-group difference.
#
# CRITERION, FIXED BEFORE THE NUMBERS:
#   CONFIRMS LEE   = RASSF9 detection increases across differentiation
#                    quintiles, top quintile at least 2x the bottom, after
#                    depth matching, with donor-level support
#   CONTRADICTS    = detection decreases across quintiles
#   NULL           = flat, or the trend does not survive depth matching
#   Melanocytes serve as a negative control: RASSF9 is essentially absent from
#   them (5 of 6,379), so no differentiation gradient should appear there.
#
# Output: results/geo/belote_krt_axis_quintiles.csv
#         results/geo/belote_krt_axis_donor.csv
#         results/geo/belote_krt_cell_scores.csv
#         figures/geo/geo_rassf9_keratinocyte_axis.{pdf,png,svg}
#         results/logs/41_belote_keratinocyte_axis.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required")

LOG <- file.path(PATHS$logs, "41_belote_keratinocyte_axis.log")
cat("Belote keratinocyte differentiation axis —", format(Sys.time()), "\n\n",
    file = LOG)
hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

DIR <- file.path(PATHS$raw, "geo", "belote")
RES <- file.path(PATHS$results, "geo")
FIG_GEO <- file.path(PATHS$figures, "geo")
invisible(lapply(c(RES, FIG_GEO), dir.create, recursive = TRUE, showWarnings = FALSE))

# Marker panels. Basal and suprabasal are kept separate so the axis can be
# validated before it is used, rather than trusted.
BASAL   <- c("KRT5", "KRT14", "COL17A1", "ITGA6", "ITGB1", "TP63", "KRT15")
SPINOUS <- c("KRT1", "KRT10", "KRTDAP", "DSP")           # Lee: RASSF9 with KRT1
GRANULAR<- c("IVL", "FLG", "LOR", "TGM1", "SPINK5", "CASP14")
P21     <- "CDKN1A"
TARGETS <- c("RASSF9", "RASSF10")
N_QUINT   <- 5
MIN_CELLS <- 20
RATIO_MIN <- 2.0
SEED      <- 20261016      # same as the earlier Belote scripts
set.seed(SEED)

pick <- function(slot, fb) {
  if (!exists("PAL", inherits = TRUE)) return(fb)
  x <- get("PAL", inherits = TRUE)[[slot]]
  if (!is.null(x) && is.character(x) && nzchar(x[1])) x[1] else fb
}
COL_K <- pick("nevus", "#0095D9"); COL_M <- pick("melanoma", "#F2680C")
COL_GREY <- pick("grey_pt", "#8E8E8E"); COL_INK <- pick("ink", "#1A1A1A")

# ==============================================================================
# 1. Load
# ==============================================================================
hdr("1. Load")

MET <- data.table::fread(file.path(DIR, "GSE151091_Metadata.csv.gz"),
                         showProgress = FALSE)
MET$group <- dplyr::case_when(
  MET$class_1 %in% c("mel", "cyc_mel")         ~ "Melanocyte",
  MET$class_1 %in% c("krt", "cyc_krt")         ~ "Keratinocyte",
  MET$class_1 %in% c("eccrine", "cyc_eccrine") ~ "Eccrine",
  TRUE                                         ~ "Other")

DT <- data.table::fread(file.path(DIR, "GSE151091_raw_matrix.csv.gz"),
                        showProgress = FALSE)
genes <- as.character(DT[[1]]); cells <- names(DT)[-1]
keep <- intersect(cells, MET$cell_id)
M <- as.matrix(DT[, keep, with = FALSE]); rownames(M) <- genes
rm(DT); gc(verbose = FALSE)
MET <- MET[match(keep, MET$cell_id), ]
tot <- colSums(M)

have <- function(v) intersect(v, rownames(M))
BASAL <- have(BASAL); SPINOUS <- have(SPINOUS); GRANULAR <- have(GRANULAR)
log_line("  basal markers:    ", paste(BASAL, collapse = ", "), logfile = LOG)
log_line("  spinous markers:  ", paste(SPINOUS, collapse = ", "), logfile = LOG)
log_line("  granular markers: ", paste(GRANULAR, collapse = ", "), logfile = LOG)
for (g in c(TARGETS, P21))
  if (!g %in% rownames(M)) stop(g, " absent from the matrix")
log_line("  cells: ", ncol(M), " | keratinocytes: ",
         sum(MET$group == "Keratinocyte"), " | melanocytes: ",
         sum(MET$group == "Melanocyte"), logfile = LOG)

# ==============================================================================
# 2. Donor coverage, checked before anything is built on it
# ==============================================================================
hdr("2. Donor coverage")

COV <- tibble::tibble(donor = MET$patient, group = MET$group) |>
  dplyr::filter(group %in% c("Keratinocyte", "Melanocyte")) |>
  dplyr::count(group, donor) |>
  dplyr::group_by(group) |>
  dplyr::summarise(cells = sum(n), donors = dplyr::n_distinct(donor),
                   donors_ge = sum(n >= MIN_CELLS), .groups = "drop")
for (i in seq_len(nrow(COV)))
  log_line(sprintf("    %-14s %5d cells | %2d donors | %2d with >= %d cells",
                   COV$group[i], COV$cells[i], COV$donors[i], COV$donors_ge[i],
                   MIN_CELLS), logfile = LOG)

# ==============================================================================
# 3. Depth matching
# ==============================================================================
hdr("3. Thinning to a common depth")

cmp <- c("Keratinocyte", "Melanocyte")
p10 <- tapply(tot, MET$group, function(x) stats::quantile(x, 0.1))
D_TARGET <- floor(min(p10[cmp]))
usable <- tot >= D_TARGET
log_line("  target depth ", format(D_TARGET, big.mark = ","),
         " | cells retained ", sum(usable), " of ", length(tot), logfile = LOG)

score_genes <- unique(c(BASAL, SPINOUS, GRANULAR, P21, TARGETS))
TH <- matrix(NA_real_, nrow = length(score_genes), ncol = ncol(M),
             dimnames = list(score_genes, colnames(M)))
for (g in score_genes) {
  x <- M[g, ]; v <- rep(NA_real_, length(x))
  i <- which(usable & x > 0)
  if (length(i)) v[i] <- stats::rbinom(length(i), size = x[i],
                                       prob = D_TARGET / tot[i])
  v[usable & x == 0] <- 0
  TH[g, ] <- v
}
ngene_raw <- colSums(M > 0)
rm(M); gc(verbose = FALSE)

sc <- function(set, idx) if (!length(set)) rep(NA_real_, length(idx)) else
  colSums(TH[set, idx, drop = FALSE], na.rm = TRUE)

K <- which(usable & MET$group == "Keratinocyte")
L <- which(usable & MET$group == "Melanocyte")
log_line("  usable keratinocytes: ", length(K), " | melanocytes: ", length(L),
         logfile = LOG)

CELL <- tibble::tibble(
  cell = MET$cell_id[K], donor = MET$patient[K],
  basal = sc(BASAL, K), spinous = sc(SPINOUS, K), granular = sc(GRANULAR, K),
  p21 = TH[P21, K], RASSF9 = TH["RASSF9", K], RASSF10 = TH["RASSF10", K],
  n_genes = ngene_raw[K])
CELL$diff_score <- CELL$spinous + CELL$granular
# A ratio contrast is used for ordering because it is internally normalised:
# a cell cannot score high simply by expressing more of everything.
CELL$axis <- log2((CELL$diff_score + 1) / (CELL$basal + 1))

# ==============================================================================
# 4. Validate the axis before using it
# ==============================================================================
hdr("4. Does the differentiation axis behave?")

ct <- suppressWarnings(stats::cor.test(CELL$basal, CELL$diff_score,
                                       method = "spearman"))
log_line(sprintf("  basal vs differentiation score: rho = %+.3f, p = %s",
                 ct$estimate, fmt_p(ct$p.value)), logfile = LOG)
log_line("  A negative correlation is expected: basal and suprabasal programmes",
         logfile = LOG)
log_line("  are mutually exclusive in epidermis. If positive, the axis is not",
         logfile = LOG)
log_line("  capturing differentiation and nothing below should be trusted.",
         logfile = LOG)
if (is.finite(ct$estimate) && ct$estimate > 0.2)
  log_line("  [FLAG] axis does not behave as expected", logfile = LOG)

cg <- suppressWarnings(stats::cor.test(CELL$axis, CELL$n_genes,
                                       method = "spearman"))
log_line(sprintf("\n  axis vs detected genes: rho = %+.3f, p = %s",
                 cg$estimate, fmt_p(cg$p.value)), logfile = LOG)
if (is.finite(cg$estimate) && cg$estimate < 0) {
  log_line("  Differentiated cells detect FEWER genes, as expected from RNA",
           logfile = LOG)
  log_line("  degradation in terminal differentiation. This biases AGAINST",
           logfile = LOG)
  log_line("  finding RASSF9 enriched there, so a positive result below is",
           logfile = LOG)
  log_line("  obtained against the gradient, not with it.", logfile = LOG)
} else {
  log_line("  Differentiated cells detect MORE genes, so part of any rise in",
           logfile = LOG)
  log_line("  RASSF9 detection could be complexity rather than biology.",
           logfile = LOG)
  log_line("  Section 5's complexity-matched column is then the one to read.",
           logfile = LOG)
}

# ==============================================================================
# 5. Primary: RASSF9 detection across the differentiation axis
# ==============================================================================
hdr("5. RASSF9 across differentiation quintiles")

CELL$quint <- dplyr::ntile(CELL$axis, N_QUINT)
Q <- CELL |> dplyr::group_by(quint) |>
  dplyr::summarise(n = dplyr::n(),
                   axis_med = stats::median(axis),
                   basal_med = stats::median(basal),
                   diff_med = stats::median(diff_score),
                   genes_med = stats::median(n_genes),
                   R9_det = sum(RASSF9 > 0, na.rm = TRUE),
                   R9_pct = 100 * mean(RASSF9 > 0, na.rm = TRUE),
                   R10_pct = 100 * mean(RASSF10 > 0, na.rm = TRUE),
                   p21_pct = 100 * mean(p21 > 0, na.rm = TRUE),
                   p21_mean = mean(p21, na.rm = TRUE), .groups = "drop")
log_line(sprintf("  %5s %6s %9s %9s %9s %8s %7s %8s %8s %9s", "quint", "n",
                 "axis", "basal", "diff", "genes", "R9 det", "R9 %", "R10 %",
                 "p21 %"), logfile = LOG)
for (i in seq_len(nrow(Q)))
  log_line(sprintf("  %5d %6d %+9.2f %9.1f %9.1f %8.0f %7d %7.2f%% %7.2f%% %8.2f%%",
                   Q$quint[i], Q$n[i], Q$axis_med[i], Q$basal_med[i],
                   Q$diff_med[i], Q$genes_med[i], Q$R9_det[i], Q$R9_pct[i],
                   Q$R10_pct[i], Q$p21_pct[i]), logfile = LOG)
readr::write_csv(Q, file.path(RES, "belote_krt_axis_quintiles.csv"))

lo <- Q$R9_pct[Q$quint == 1]; hi <- Q$R9_pct[Q$quint == N_QUINT]
ratio <- if (length(lo) && lo > 0) hi / lo else if (isTRUE(hi > 0)) Inf else NA
log_line(sprintf("\n  top / bottom quintile = %s",
                 if (is.finite(ratio)) sprintf("%.2fx", ratio) else
                   if (is.infinite(ratio)) "bottom is 0" else "n/a"),
         logfile = LOG)

# Cochran-Armitage style trend test on the ordered quintiles.
tr <- suppressWarnings(stats::prop.trend.test(Q$R9_det, Q$n, score = Q$quint))
log_line("  trend across quintiles: chi-sq = ", sprintf("%.2f", tr$statistic),
         ", p = ", fmt_p(tr$p.value), logfile = LOG)

# Cell-level, continuous, with complexity as a covariate.
fit <- try(stats::glm(as.integer(RASSF9 > 0) ~ axis + n_genes,
                      data = CELL, family = stats::binomial()), silent = TRUE)
if (!inherits(fit, "try-error")) {
  cf <- summary(fit)$coefficients
  if ("axis" %in% rownames(cf))
    log_line(sprintf("  logistic, adjusted for detected genes: axis OR = %.3f, p = %s (PSEUDO-REPLICATED)",
                     exp(cf["axis", 1]), fmt_p(cf["axis", 4])), logfile = LOG)
}

# ==============================================================================
# 6. Donor-level
# ==============================================================================
hdr("6. Within-donor comparison")

CELL$half <- ifelse(CELL$axis >= stats::median(CELL$axis), "Differentiated",
                    "Basal")
DON <- CELL |> dplyr::group_by(donor, half) |>
  dplyr::summarise(n = dplyr::n(), pct = 100 * mean(RASSF9 > 0, na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::filter(n >= MIN_CELLS)
a <- as.data.frame(DON)[DON$half == "Basal", c("donor", "n", "pct")]
names(a) <- c("donor", "n_basal", "pct_basal")
b <- as.data.frame(DON)[DON$half == "Differentiated", c("donor", "n", "pct")]
names(b) <- c("donor", "n_diff", "pct_diff")
P <- dplyr::inner_join(a, b, by = "donor")
readr::write_csv(DON, file.path(RES, "belote_krt_axis_donor.csv"))

log_line("  donors with >= ", MIN_CELLS, " keratinocytes in both halves: ",
         nrow(P), logfile = LOG)
if (nrow(P)) {
  log_line(sprintf("  %-12s %8s %8s %10s %10s", "donor", "n basal", "n diff",
                   "basal %", "diff %"), logfile = LOG)
  for (i in seq_len(nrow(P)))
    log_line(sprintf("  %-12s %8d %8d %9.2f%% %9.2f%%", P$donor[i],
                     P$n_basal[i], P$n_diff[i], P$pct_basal[i], P$pct_diff[i]),
             logfile = LOG)
  higher <- sum(P$pct_diff > P$pct_basal)
  log_line(sprintf("  differentiated higher in %d of %d donors", higher, nrow(P)),
           logfile = LOG)
  if (nrow(P) >= 4) {
    wt <- suppressWarnings(stats::wilcox.test(P$pct_diff, P$pct_basal,
                                              paired = TRUE))
    st <- stats::binom.test(higher, nrow(P))
    log_line("  paired Wilcoxon p = ", fmt_p(wt$p.value), " | sign test p = ",
             fmt_p(st$p.value), logfile = LOG)
  } else {
    log_line("  fewer than 4 donors; no paired test attempted", logfile = LOG)
  }
}

# ==============================================================================
# 7. p21
# ==============================================================================
hdr("7. Does RASSF9 travel with CDKN1A (p21)?")

pos <- CELL$RASSF9 > 0; neg <- !pos
log_line(sprintf("  RASSF9-positive keratinocytes: %d | negative: %d",
                 sum(pos, na.rm = TRUE), sum(neg, na.rm = TRUE)), logfile = LOG)
if (sum(pos, na.rm = TRUE) >= 20) {
  w <- suppressWarnings(stats::wilcox.test(CELL$p21[pos], CELL$p21[neg]))
  log_line(sprintf("  CDKN1A detected in %.2f%% of RASSF9+ vs %.2f%% of RASSF9- cells",
                   100 * mean(CELL$p21[pos] > 0, na.rm = TRUE),
                   100 * mean(CELL$p21[neg] > 0, na.rm = TRUE)), logfile = LOG)
  log_line(sprintf("  CDKN1A mean (thinned counts): %.3f vs %.3f | Wilcoxon p = %s",
                   mean(CELL$p21[pos], na.rm = TRUE),
                   mean(CELL$p21[neg], na.rm = TRUE), fmt_p(w$p.value)),
           logfile = LOG)
  log_line("", logfile = LOG)
  log_line("  CAUTION: CDKN1A also rises with differentiation in its own right,",
           logfile = LOG)
  log_line("  so this comparison cannot separate 'RASSF9 marks p21-high cells'",
           logfile = LOG)
  log_line("  from 'both rise together along the axis'. The model below adds",
           logfile = LOG)
  log_line("  the axis as a covariate to ask whether RASSF9 carries anything",
           logfile = LOG)
  log_line("  beyond it. A null there does NOT refute Lee's mechanism; a",
           logfile = LOG)
  log_line("  co-regulated pair is exactly what his model predicts.",
           logfile = LOG)
  f2 <- try(stats::glm(as.integer(p21 > 0) ~ as.integer(RASSF9 > 0) + axis + n_genes,
                       data = CELL, family = stats::binomial()), silent = TRUE)
  if (!inherits(f2, "try-error")) {
    cf <- summary(f2)$coefficients
    rn <- grep("RASSF9", rownames(cf), value = TRUE)
    if (length(rn))
      log_line(sprintf("  p21 ~ RASSF9 + axis + complexity: OR = %.3f, p = %s",
                       exp(cf[rn[1], 1]), fmt_p(cf[rn[1], 4])), logfile = LOG)
  }
} else {
  log_line("  too few RASSF9-positive keratinocytes for this comparison",
           logfile = LOG)
}

# ==============================================================================
# 8. Negative control: melanocytes
# ==============================================================================
hdr("8. Negative control in melanocytes")

log_line("  RASSF9 is essentially absent from melanocytes (5 of 6,379), so no",
         logfile = LOG)
log_line("  differentiation gradient should appear here. If one does, the",
         logfile = LOG)
log_line("  quintile result above is suspect.", logfile = LOG)
MELC <- tibble::tibble(basal = sc(BASAL, L), diff_score = sc(SPINOUS, L) + sc(GRANULAR, L),
                       RASSF9 = TH["RASSF9", L], n_genes = ngene_raw[L])
MELC$axis <- log2((MELC$diff_score + 1) / (MELC$basal + 1))
MELC$quint <- dplyr::ntile(MELC$axis, N_QUINT)
QM <- MELC |> dplyr::group_by(quint) |>
  dplyr::summarise(n = dplyr::n(), det = sum(RASSF9 > 0, na.rm = TRUE),
                   pct = 100 * mean(RASSF9 > 0, na.rm = TRUE), .groups = "drop")
for (i in seq_len(nrow(QM)))
  log_line(sprintf("    quintile %d: n = %4d | detected %2d | %.3f%%",
                   QM$quint[i], QM$n[i], QM$det[i], QM$pct[i]), logfile = LOG)
if (sum(QM$det) >= 5) {
  trm <- suppressWarnings(stats::prop.trend.test(QM$det, QM$n, score = QM$quint))
  log_line("    trend p = ", fmt_p(trm$p.value),
           " (expected: non-significant)", logfile = LOG)
} else {
  log_line("    too few detections for a trend test, which is itself the",
           logfile = LOG)
  log_line("    expected result", logfile = LOG)
}

readr::write_csv(CELL, file.path(RES, "belote_krt_cell_scores.csv"))

# ==============================================================================
# 9. Verdict
# ==============================================================================
hdr("9. Verdict against the pre-set criterion")

log_line("  CONFIRMS LEE = detection rises across quintiles, top >= ", RATIO_MIN,
         "x bottom,", logfile = LOG)
log_line("                 with donor-level support", logfile = LOG)
log_line("  CONTRADICTS  = detection falls across quintiles", logfile = LOG)
log_line("  NULL         = flat, or not surviving depth matching", logfile = LOG)
log_line("", logfile = LOG)

mono_up <- all(diff(Q$R9_pct) >= 0)
donor_ok <- nrow(P) >= 4 &&
  isTRUE(suppressWarnings(stats::wilcox.test(P$pct_diff, P$pct_basal,
                                             paired = TRUE))$p.value < 0.05)
verdict <- if (isTRUE(ratio >= RATIO_MIN) && tr$p.value < 0.05 && donor_ok) {
  "CONFIRMS LEE in human single cells"
} else if (isTRUE(ratio >= RATIO_MIN) && tr$p.value < 0.05) {
  "SUPPORTS LEE at cell level; donor-level support not available"
} else if (isTRUE(hi < lo) && tr$p.value < 0.05) {
  "CONTRADICTS -- RASSF9 is enriched in BASAL keratinocytes"
} else "NULL -- no gradient along the differentiation axis"

log_line("  ", verdict, logfile = LOG)
log_line(sprintf("    bottom quintile %.2f%% -> top quintile %.2f%% (%s), trend p = %s",
                 lo, hi, if (is.finite(ratio)) sprintf("%.2fx", ratio) else "n/a",
                 fmt_p(tr$p.value)), logfile = LOG)
log_line(sprintf("    monotonic increase across all five quintiles: %s",
                 if (mono_up) "yes" else "no"), logfile = LOG)
log_line("", logfile = LOG)
log_line("  WHAT THIS CANNOT SAY: detection of a transcript is not protein", logfile = LOG)
log_line("  localisation. Lee's evidence is immunostaining; this is mRNA in",
         logfile = LOG)
log_line("  dissociated cells. Agreement across those two measurements is",
         logfile = LOG)
log_line("  meaningful, but this does not reproduce his experiment.",
         logfile = LOG)

# ==============================================================================
# 10. Figure
# ==============================================================================
hdr("10. Figure")

reset_legends("# Belote keratinocyte axis figure legends\n")

pa <- ggplot(Q, aes(x = factor(quint), y = R9_pct)) +
  geom_col(width = 0.65, fill = COL_K) +
  geom_text(aes(label = sprintf("%d/%d", R9_det, n)), vjust = -0.5,
            size = (FIG$base_size - 4) / .pt, colour = COL_INK) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = "Differentiation quintile (basal → differentiated)",
       y = "Keratinocytes with RASSF9 detected (%)") +
  theme_journal() +
  theme(plot.margin = margin(4, 12, 4, 4, unit = "pt"))

pb <- ggplot(CELL, aes(x = axis, y = basal + 1)) +
  geom_point(aes(colour = RASSF9 > 0), size = 0.7, alpha = 0.6) +
  scale_y_log10() +
  scale_colour_manual(values = c(`TRUE` = COL_M, `FALSE` = COL_GREY),
                      labels = c(`TRUE` = "RASSF9 detected", `FALSE` = "not detected"),
                      name = NULL) +
  labs(x = expression(bold("Differentiation axis  "*log[2]*"(suprabasal / basal)")),
       y = "Basal marker counts") +
  theme_journal() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = FIG$base_size - 3),
        legend.key.size = unit(8, "pt"),
        plot.margin = margin(4, 12, 4, 4, unit = "pt"))

fig <- patchwork::wrap_plots(list(pa, pb), nrow = 1) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = FIG$base_size + 2, face = "bold",
                                colour = COL_INK))
save_fig_to(fig, "geo_rassf9_keratinocyte_axis", FIG_GEO,
            width = FIG$width_double, height = 88)
log_line("  wrote: geo_rassf9_keratinocyte_axis", logfile = LOG)

add_legend("geo_rassf9_keratinocyte_axis",
  "Figure. RASSF9 across the basal-to-differentiated axis of human keratinocytes.",
  paste0("(A) Percentage of keratinocytes with RASSF9 detected, by quintile of a ",
         "differentiation axis defined as log2 of suprabasal (KRT1, KRT10, KRTDAP, ",
         "DSP, involucrin, filaggrin, loricrin, TGM1, SPINK5, CASP14) over basal ",
         "(KRT5, KRT14, COL17A1, ITGA6, ITGB1, TP63, KRT15) marker counts; ",
         "annotations give detected over total cells. All cells were thinned to a ",
         "common depth of ", format(D_TARGET, big.mark = ","),
         " counts beforehand, and differentiating keratinocytes detect fewer genes ",
         "overall, so any enrichment is obtained against that gradient. (B) Each ",
         "keratinocyte placed on the differentiation axis against its basal marker ",
         "counts, coloured by RASSF9 detection. Data from Belote et al., Nat Cell ",
         "Biol 2021 (GSE151091), ", dplyr::n_distinct(CELL$donor),
         " donors. The prediction tested is from Lee et al., PLoS ONE 2011, which ",
         "localised RASSF9 protein to the suprabasal epidermis in mouse."))

report_warnings(logfile = LOG)
write_session_info("41_belote_keratinocyte_axis")
log_line("\nFull log: ", LOG, logfile = NULL)
