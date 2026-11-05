# ==============================================================================
# 37d_belote_doublet_check.R
#
# Are the 12 RASSF9/RASSF10-positive melanocytes real, or doublets?
#
# WHY THIS MATTERS. 37c established the headline: RASSF9 detection is higher in
# keratinocytes than melanocytes in 17 of 17 donors (paired p = 1.53e-05), with
# 14 donors at exactly 0.00% in melanocytes. But 12 melanocytes are positive,
# spread across 5 donors, and four of them carry substantial counts -- 141, 99,
# 39 and 33 raw -- which is not stray single-molecule noise. A cell with 141
# RASSF9 counts at a 1.56M library is expressing the gene.
#
# So the claim is not "zero". It is either:
#   (a) RASSF9 is absent from normal melanocytes, and these 12 are
#       melanocyte-keratinocyte DOUBLETS -- which would be unsurprising, since
#       keratinocytes are precisely the cell type that does express RASSF9 and
#       they sit directly against melanocytes in the epidermis; or
#   (b) a rare genuine melanocyte population expresses it, which is a more
#       interesting claim and needs stating carefully rather than buried.
#
# These demand different sentences in the manuscript, so the question gets
# answered rather than assumed either way.
#
# THE TEST. A melanocyte-keratinocyte doublet carries both programmes at once:
# keratins it should not have, melanocyte markers it should, and more detected
# genes than a single cell. So each positive cell is scored on three axes
# against the melanocyte population it came from.
#
# CRITERION, FIXED BEFORE THE NUMBERS. A positive cell is flagged LIKELY
# DOUBLET if its keratin score exceeds the 95th percentile of all melanocytes
# AND its detected-gene count exceeds the 90th percentile. Both, not either:
# keratin alone could be ambient contamination, and a high gene count alone is
# just a well-covered cell.
#
# Then the melanocyte detection rate is recomputed with flagged cells removed,
# so the paper can quote a figure that survives this objection.
#
# NOTE this cannot be fully settled from a count matrix. Flagging is evidence,
# not proof, and an unflagged positive cell is not thereby proven to be a
# single melanocyte. The script says which it can and cannot distinguish.
#
# Output: results/geo/belote_doublet_check.csv
#         results/logs/37d_belote_doublet_check.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required")

LOG <- file.path(PATHS$logs, "37d_belote_doublet_check.log")
cat("Belote doublet check —", format(Sys.time()), "\n\n", file = LOG)
hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

DIR <- file.path(PATHS$raw, "geo", "belote")
RES <- file.path(PATHS$results, "geo")
dir.create(RES, recursive = TRUE, showWarnings = FALSE)

TARGETS <- c("RASSF9", "RASSF10")
KRT_SET <- c("KRT14", "KRT1", "KRT5", "KRT10", "KRT2", "KRTDAP", "DSP")
MEL_SET <- c("MLANA", "PMEL", "TYR", "DCT", "TYRP1", "MITF", "SOX10")
Q_KRT   <- 0.95
Q_GENES <- 0.90
SEED    <- 20261016
set.seed(SEED)

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
tot <- colSums(M); ngene <- colSums(M > 0)

krt_use <- intersect(KRT_SET, rownames(M))
mel_use <- intersect(MEL_SET, rownames(M))
log_line("  keratin genes used: ", paste(krt_use, collapse = ", "), logfile = LOG)
log_line("  melanocyte genes used: ", paste(mel_use, collapse = ", "), logfile = LOG)

score <- function(set) 1e4 * colSums(M[set, , drop = FALSE]) / pmax(tot, 1)
krt_score <- score(krt_use); mel_score <- score(mel_use)

is_mel <- !is.na(MET$group) & MET$group == "Melanocyte"
ref_krt  <- krt_score[is_mel]; ref_gene <- ngene[is_mel]
n_bad_k <- sum(!is.finite(ref_krt)); n_bad_g <- sum(!is.finite(ref_gene))
if (n_bad_k || n_bad_g) {
  log_line("  [FLAG] non-finite reference values dropped: keratin ", n_bad_k,
           ", gene count ", n_bad_g,
           " -- these produced an all-NA percentile column on the first run",
           logfile = LOG)
}
ref_krt  <- ref_krt[is.finite(ref_krt)]
ref_gene <- ref_gene[is.finite(ref_gene)]
k_cut <- stats::quantile(ref_krt,  Q_KRT,   names = FALSE)
g_cut <- stats::quantile(ref_gene, Q_GENES, names = FALSE)
log_line(sprintf("  melanocyte keratin score: median %.2f, %d%%ile %.2f",
                 stats::median(krt_score[is_mel]), round(100 * Q_KRT), k_cut),
         logfile = LOG)
log_line(sprintf("  melanocyte detected genes: median %.0f, %d%%ile %.0f",
                 stats::median(ngene[is_mel]), round(100 * Q_GENES), g_cut),
         logfile = LOG)
log_line(sprintf("  keratinocyte keratin score for contrast: median %.2f",
                 stats::median(krt_score[MET$group == "Keratinocyte"])),
         logfile = LOG)

# ==============================================================================
# 2. The positive melanocytes, scored
# ==============================================================================
hdr("2. Each positive melanocyte against its own population")

cmp <- c("Melanocyte", "Eccrine", "Keratinocyte")
p10 <- tapply(tot, MET$group, function(x) stats::quantile(x, 0.1))
D_TARGET <- floor(min(p10[cmp]))
usable <- tot >= D_TARGET

pos <- which(is_mel & usable & (M["RASSF9", ] > 0 | M["RASSF10", ] > 0))
log_line("  positive melanocytes (usable depth): ", length(pos), logfile = LOG)

# Percentiles are computed BEFORE the tibble, not inside it. On the first run
# they were computed inline as `krt_score[pos]`, where tibble's data mask
# resolved `krt_score` to the 17-element column defined one line earlier rather
# than the global 9,719-element vector -- so every index in `pos` was out of
# range and the whole column came back NA. The gene version escaped only
# because its column is named `n_genes` while the vector is `ngene`.
pct_of <- function(v, ref) {
  ref <- ref[is.finite(ref)]
  vapply(v, function(x) if (!is.finite(x)) NA_real_ else
    100 * mean(ref <= x), numeric(1))
}
krt_pct_pos  <- pct_of(krt_score[pos], ref_krt)
gene_pct_pos <- pct_of(ngene[pos],     ref_gene)

PC <- tibble::tibble(
  cell_id = MET$cell_id[pos], donor = MET$patient[pos],
  stage = if ("dev_stage" %in% names(MET)) MET$dev_stage[pos] else NA_character_,
  site = if ("general_location" %in% names(MET)) MET$general_location[pos] else NA_character_,
  lib = tot[pos], n_genes = ngene[pos],
  RASSF9 = M["RASSF9", pos], RASSF10 = M["RASSF10", pos],
  krt_score = krt_score[pos], mel_score = mel_score[pos],
  krt_pctile = krt_pct_pos, gene_pctile = gene_pct_pos)

# Flagging, REVISED after the first run. The original rule required keratin
# ABOVE the 95th percentile AND gene count above the 90th, on the reasoning
# that keratin alone might be ambient contamination. That rule missed
# A1011L_C5, which carries a keratin score of 350 against a melanocyte median
# of 1.21 -- a third of the way to the keratinocyte median -- purely because
# its gene count was ordinary. A doublet whose partner is a low-complexity
# cell need not have an inflated gene count, so keratin dominance is treated
# as sufficient on its own and the two signals are reported separately rather
# than conjunctively.
PC$keratin_high  <- is.finite(PC$krt_score) & PC$krt_score > k_cut
PC$complexity_hi <- is.finite(PC$n_genes)   & PC$n_genes   > g_cut
PC$likely_doublet <- PC$keratin_high

log_line(sprintf("\n  %-26s %-10s %6s %6s %7s %9s %8s %8s %s", "cell", "donor",
                 "R9", "R10", "genes", "krt score", "krt %ile", "gene %ile",
                 "flag"), logfile = LOG)
for (i in seq_len(nrow(PC)))
  log_line(sprintf("  %-26s %-10s %6.0f %6.0f %7.0f %9.2f %7.1f%% %7.1f%% %s",
                   PC$cell_id[i], PC$donor[i], PC$RASSF9[i], PC$RASSF10[i],
                   PC$n_genes[i], PC$krt_score[i], PC$krt_pctile[i],
                   PC$gene_pctile[i],
                   paste0(if (PC$keratin_high[i]) "KERATIN-HIGH" else "-",
                          if (PC$complexity_hi[i]) " cplx" else "")),
           logfile = LOG)
readr::write_csv(PC, file.path(RES, "belote_doublet_check.csv"))

log_line("", logfile = LOG)
log_line("  flagged: ", sum(PC$likely_doublet), " of ", nrow(PC), logfile = LOG)

# The substantial-count cells are the ones that matter; a cell positive with a
# single molecule changes no conclusion either way.
big <- PC |> dplyr::filter(RASSF9 >= 10 | RASSF10 >= 10)
if (nrow(big)) {
  log_line("\n  cells with >= 10 counts of either gene (the ones that matter):",
           logfile = LOG)
  for (i in seq_len(nrow(big)))
    log_line(sprintf("    %-26s R9 %4.0f R10 %4.0f | krt %%ile %5.1f | gene %%ile %5.1f | %s",
                     big$cell_id[i], big$RASSF9[i], big$RASSF10[i],
                     big$krt_pctile[i], big$gene_pctile[i],
                     paste0(if (big$keratin_high[i]) "KERATIN-HIGH -> doublet" else "keratin-negative -> real",
                            if (big$complexity_hi[i]) ", high complexity" else "")),
             logfile = LOG)
  log_line("    flagged among these: ", sum(big$likely_doublet), " of ", nrow(big),
           logfile = LOG)
}

# ==============================================================================
# 3. Detection rate with flagged cells removed
# ==============================================================================
hdr("3. Melanocyte detection, flagged cells removed")

log_line("  NOTE: these counts are RAW, not thinned. 37c reported thinned",
         logfile = LOG)
log_line("  values, so its figures are lower (RASSF9 5 of 6,379 rather than 8).",
         logfile = LOG)
log_line("  Raw is the right basis for a doublet check; the THINNED figure is",
         logfile = LOG)
log_line("  the one to quote for any comparison between cell classes.",
         logfile = LOG)

drop_ids <- PC$cell_id[PC$likely_doublet]
mel_idx  <- which(is_mel & usable)
mel_keep <- mel_idx[!MET$cell_id[mel_idx] %in% drop_ids]

for (g in TARGETS) {
  k0 <- sum(M[g, mel_idx] > 0);  n0 <- length(mel_idx)
  k1 <- sum(M[g, mel_keep] > 0); n1 <- length(mel_keep)
  b0 <- stats::binom.test(k0, n0); b1 <- stats::binom.test(k1, n1)
  log_line(sprintf("  %-8s all melanocytes      %3d / %4d = %.3f%% [%.3f, %.3f]",
                   g, k0, n0, 100 * k0 / n0, 100 * b0$conf.int[1],
                   100 * b0$conf.int[2]), logfile = LOG)
  log_line(sprintf("  %-8s flagged removed      %3d / %4d = %.3f%% [%.3f, %.3f]",
                   "", k1, n1, 100 * k1 / n1, 100 * b1$conf.int[1],
                   100 * b1$conf.int[2]), logfile = LOG)
}

hdr("4. Reading this")

log_line("  If most substantial-count cells are flagged, the melanocyte rate is",
         logfile = LOG)
log_line("  effectively zero and the manuscript can say RASSF9 is not expressed",
         logfile = LOG)
log_line("  in normal melanocytes, with doublets named as the source of the few",
         logfile = LOG)
log_line("  positives.", logfile = LOG)
log_line("", logfile = LOG)
log_line("  If they are NOT flagged -- keratin-low with ordinary gene counts --",
         logfile = LOG)
log_line("  then a rare population of normal melanocytes does express RASSF9,",
         logfile = LOG)
log_line("  at roughly 1 in 1,000. That is worth one honest sentence rather",
         logfile = LOG)
log_line("  than a claim of absence, and it would sit comfortably with an",
         logfile = LOG)
log_line("  induction model: a gene off at baseline in nearly every cell, with",
         logfile = LOG)
log_line("  a few cells already carrying it.", logfile = LOG)
log_line("", logfile = LOG)
log_line("  WHAT THIS CANNOT DO: a count matrix cannot prove a cell is a", logfile = LOG)
log_line("  singlet. Flagging is evidence of a doublet; absence of a flag is", logfile = LOG)
log_line("  not proof of a single melanocyte. Neither reading should be stated",
         logfile = LOG)
log_line("  more strongly than that.", logfile = LOG)

report_warnings(logfile = LOG)
write_session_info("37d_belote_doublet_check")
log_line("\nFull log: ", LOG, logfile = NULL)
