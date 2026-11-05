# ==============================================================================
# 35_methylation_full_cohort.R
#
# The complete TCGA-SKCM cohort, and every probe at the locus.
#
# TWO THINGS AT ONCE, because both need the raw objects rather than the saved
# cohort files:
#
# 1. THE TOTAL COHORT. Script 34's "all cutaneous" tier came out at 411, which
#    is every tumour with mutational signature data, UV-positive or not. But
#    methylation exists for 470 samples, and the 59 missing ones are absent not
#    because they lack RNA-seq but because they were never in the cohort files:
#    those were built by joining to the SigProfiler signature table, so any
#    tumour without signature data dropped out before anything was measured.
#    Rebuilding from skcm_star_counts.rds and skcm_hm450.rds directly removes
#    that requirement entirely.
#
# 2. EVERY PROBE AT THE LOCUS. The saved probe table holds 23 probes, of which
#    RASSF9 has 7 and all 7 are promoter. That left RASSF9 with no external
#    specificity control -- a gap flagged in script 34 but not closeable from
#    the saved files. The HM450 object carries all 485,577 probes, so the
#    probes near RASSF9 can be looked up directly rather than inferred from
#    what a previous script happened to retain. If non-promoter RASSF9 probes
#    exist, this finds them; if there genuinely are none, that is now a
#    verified fact rather than an assumption.
#
# WHAT IS PRE-SPECIFIED, before the numbers:
#   - The promoter-mean correlation in the full cohort is expected near -0.37,
#     between the 411-sample value (-0.374) and the no-UV-signature subset
#     (-0.409). A materially different value means the 59 added samples differ
#     systematically, so section 5 characterises them rather than leaving that
#     as a surprise.
#   - The specificity claim is that probes near the TSS correlate negatively
#     while probes elsewhere at the locus do not behave the same way. With the
#     full probe set this becomes a real test for RASSF9 rather than an
#     argument from the island-shore pattern alone.
#
# RUNTIME: the methylation object is 1.5 GB and takes a couple of minutes to
# load. Everything after that is fast.
#
# Output: results/tcga/methylation_full_tiers.csv
#         results/tcga/methylation_full_probes.csv
#         results/tcga/methylation_full_addedsamples.csv
#         figures/tcga/tcga_<gene>_methylation_full.{pdf,png,svg}
#         results/logs/35_methylation_full_cohort.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

LOG <- file.path(PATHS$logs, "35_methylation_full_cohort.log")
cat("TCGA-SKCM methylation, complete cohort —", format(Sys.time()), "\n\n",
    file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

PROC <- PATHS$processed
RES  <- file.path(PATHS$results, "tcga")
FIG_TCGA <- file.path(PATHS$figures, "tcga")
invisible(lapply(c(RES, FIG_TCGA), dir.create, recursive = TRUE,
                 showWarnings = FALSE))

GENES   <- c("RASSF9", "RASSF10")
# Window for collecting probes. 10 kb around the TSS was too narrow to settle
# the specificity question: RASSF9 spans roughly 90 kb, so probes in its distal
# gene body fell outside it, and "no non-promoter probes within 10 kb" is a
# weaker claim than "none at the locus". 150 kb covers both gene bodies with
# room to spare, so the non-promoter probe count is now a fact about the array
# rather than about an arbitrary radius.
WINDOW  <- 150000
PROM_LO <- -1500      # promoter window, as the original analysis defined it
PROM_HI <- 500

pick <- function(slot, fallback) {
  if (!exists("PAL", inherits = TRUE)) return(fallback)
  x <- get("PAL", inherits = TRUE)[[slot]]
  if (!is.null(x) && is.character(x) && nzchar(x[1])) x[1] else fallback
}
COL_A    <- pick("nevus",    "#0095D9")
COL_B    <- pick("melanoma", "#F2680C")
COL_GREY <- pick("grey_pt",  "#8E8E8E")
COL_INK  <- pick("ink",      "#1A1A1A")

RNA_F  <- file.path(PATHS$raw, "skcm_star_counts.rds")
METH_F <- file.path(PATHS$raw, "skcm_hm450.rds")

# ==============================================================================
# 1. Expression
# ==============================================================================
hdr("1. Expression")

if (!file.exists(RNA_F)) stop("Missing ", RNA_F)
log_line("  loading ", basename(RNA_F), sprintf(" (%.0f MB)",
         file.size(RNA_F) / 1024^2), logfile = LOG)
SE <- readRDS(RNA_F)
log_line("  class: ", paste(class(SE), collapse = "/"), logfile = LOG)
log_line("  dim: ", nrow(SE), " x ", ncol(SE), logfile = LOG)
log_line("  assays: ", paste(assayNames(SE), collapse = ", "), logfile = LOG)

a_tpm <- intersect(c("tpm_unstrand", "tpm_unstranded", "TPM"), assayNames(SE))
if (!length(a_tpm)) {
  stop("No TPM assay found; assays present: ",
       paste(assayNames(SE), collapse = ", "))
}
a_tpm <- a_tpm[1]
log_line("  using assay: ", a_tpm, logfile = LOG)

rd <- as.data.frame(rowData(SE))
gcol <- intersect(c("gene_name", "external_gene_name", "symbol"), names(rd))
if (!length(gcol)) {
  stop("No gene-symbol column in rowData; columns: ",
       paste(names(rd), collapse = ", "))
}
gcol <- gcol[1]
log_line("  gene symbols from rowData$", gcol, logfile = LOG)

# Gene coordinates come from this same object, so probe and gene positions
# share one genome build -- the point the original methods made.
GR <- rowRanges(SE)
tss_of <- function(g) {
  i <- which(rd[[gcol]] == g)
  if (!length(i)) return(NULL)
  gr <- GR[i[1]]
  list(chr = as.character(GenomicRanges::seqnames(gr)),
       strand = as.character(BiocGenerics::strand(gr)),
       tss = if (as.character(BiocGenerics::strand(gr)) == "-") {
         BiocGenerics::end(gr)
       } else {
         BiocGenerics::start(gr)
       })
}

EXPR <- list()
for (g in GENES) {
  i <- which(rd[[gcol]] == g)
  if (!length(i)) {
    log_line("  [FLAG] ", g, " not found in the expression object",
             logfile = LOG); next
  }
  v <- SummarizedExperiment::assay(SE, a_tpm)[i[1], ]
  EXPR[[g]] <- log2(as.numeric(v) + 1)
  tg <- tss_of(g)
  log_line(sprintf("  %-8s %s:%d (%s strand), median log2(TPM+1) = %.3f",
                   g, tg$chr, tg$tss, tg$strand,
                   stats::median(EXPR[[g]], na.rm = TRUE)), logfile = LOG)
}

cd <- as.data.frame(colData(SE))
stcol <- intersect(c("sample_type", "shortLetterCode", "definition"), names(cd))
RNA <- tibble::tibble(
  rna_barcode = colnames(SE),
  patient = substr(colnames(SE), 1, 12),
  sample_type = if (length(stcol)) as.character(cd[[stcol[1]]]) else NA_character_)
for (g in names(EXPR)) RNA[[paste0(g, "_log2tpm")]] <- EXPR[[g]]
log_line("  expression samples: ", nrow(RNA), " | patients: ",
         dplyr::n_distinct(RNA$patient), logfile = LOG)
rm(SE); gc(verbose = FALSE)

# ==============================================================================
# 2. Methylation, and every probe at the locus
# ==============================================================================
hdr("2. Methylation and the full probe set")

if (!file.exists(METH_F)) stop("Missing ", METH_F)
log_line("  loading ", basename(METH_F), sprintf(" (%.1f GB)",
         file.size(METH_F) / 1024^3), " -- this takes a minute or two",
         logfile = LOG)
ME <- readRDS(METH_F)
log_line("  class: ", paste(class(ME), collapse = "/"), logfile = LOG)
log_line("  dim: ", nrow(ME), " x ", ncol(ME), logfile = LOG)

MGR <- try(rowRanges(ME), silent = TRUE)
if (inherits(MGR, "try-error") || !length(MGR)) {
  stop("No rowRanges on the methylation object; probe positions unavailable.")
}

probe_rows <- list()
for (g in GENES) {
  i <- which(rd[[gcol]] == g)
  if (!length(i)) next
  gr <- GR[i[1]]
  chr <- as.character(GenomicRanges::seqnames(gr))
  str <- as.character(BiocGenerics::strand(gr))
  tss <- if (str == "-") BiocGenerics::end(gr) else BiocGenerics::start(gr)

  sel <- which(as.character(GenomicRanges::seqnames(MGR)) == chr &
               abs(BiocGenerics::start(MGR) - tss) <= WINDOW)
  if (!length(sel)) {
    log_line("  [FLAG] no probes within ", WINDOW, " bp of ", g, logfile = LOG)
    next
  }
  pos <- BiocGenerics::start(MGR[sel])
  # Strand-aware distance: upstream is negative on both strands.
  d <- if (str == "-") tss - pos else pos - tss
  probe_rows[[g]] <- tibble::tibble(
    gene = g, probe = names(MGR)[sel], chr = chr, pos = pos, dist_tss = d,
    promoter = d >= PROM_LO & d <= PROM_HI,
    region = dplyr::case_when(
      d >= PROM_LO & d <= PROM_HI ~ "Promoter",
      d < PROM_LO ~ "Upstream",
      d > PROM_HI & d <= BiocGenerics::width(gr) ~ "Gene body",
      TRUE ~ "Downstream")) |>
    dplyr::arrange(dist_tss)
  log_line(sprintf("  %-8s %d probes within %d bp of the TSS, %d in the promoter window",
                   g, length(sel), WINDOW, sum(probe_rows[[g]]$promoter)),
           logfile = LOG)
}
PRF <- dplyr::bind_rows(probe_rows)

# THE VERIFICATION. Script 34 reported that RASSF9 has no non-promoter probes,
# which was true of the retained subset but untested against the array itself.
log_line("", logfile = LOG)
for (g in GENES) {
  s <- PRF[PRF$gene == g, ]
  if (!nrow(s)) next
  np <- sum(!s$promoter)
  log_line("  ", g, ": ", np, " non-promoter probes on the array within ",
           WINDOW, " bp", logfile = LOG)
  if (np == 0) {
    log_line("    -> confirmed: no external specificity control exists for ",
             g, logfile = LOG)
  } else {
    log_line("    -> an external specificity control exists for ", g,
             ": these probes", logfile = LOG)
    log_line("       should NOT reproduce the promoter correlation if the",
             logfile = LOG)
    log_line("       signal is promoter-specific. Section 6 tests that.",
             logfile = LOG)
  }
}

keep <- intersect(PRF$probe, rownames(ME))
log_line("", logfile = LOG)
log_line("  probes retrievable from the beta matrix: ", length(keep), " of ",
         nrow(PRF), logfile = LOG)
B <- SummarizedExperiment::assay(ME)[keep, , drop = FALSE]
BE <- tibble::tibble(meth_barcode = colnames(ME),
                     patient = substr(colnames(ME), 1, 12))
BE <- dplyr::bind_cols(BE, tibble::as_tibble(t(as.matrix(B))))
log_line("  methylation samples: ", nrow(BE), " | patients: ",
         dplyr::n_distinct(BE$patient), logfile = LOG)
rm(ME, B); gc(verbose = FALSE)

PRF <- PRF[PRF$probe %in% keep, ]
readr::write_csv(PRF, file.path(RES, "methylation_full_probes.csv"))

# ==============================================================================
# 3. The complete matched cohort
# ==============================================================================
hdr("3. Matching")

RNA1 <- RNA |> dplyr::distinct(patient, .keep_all = TRUE)
BE1  <- BE  |> dplyr::distinct(patient, .keep_all = TRUE)
ALL  <- dplyr::inner_join(RNA1, BE1, by = "patient")
log_line("  patients with BOTH expression and methylation: ", nrow(ALL),
         logfile = LOG)

f_ana  <- file.path(PROC, "analysis_cohort.csv")
f_full <- file.path(PROC, "full_cohort.csv")
AN <- if (file.exists(f_ana))  readr::read_csv(f_ana,  show_col_types = FALSE,
                                               progress = FALSE) else NULL
FU <- if (file.exists(f_full)) readr::read_csv(f_full, show_col_types = FALSE,
                                               progress = FALSE) else NULL
if (!is.null(AN)) log_line("  T1 reference (analysis cohort): ", nrow(AN),
                           logfile = LOG)
if (!is.null(FU)) log_line("  T2 reference (full cohort): ", nrow(FU),
                           logfile = LOG)

ADDED <- if (!is.null(FU)) ALL[!ALL$patient %in% FU$patient, ] else ALL[0, ]
log_line("  patients in the complete set but not in the cohort files: ",
         nrow(ADDED), logfile = LOG)

# Purity, where the cohort files carry it; the added samples will not have it.
if (!is.null(FU) && "CPE" %in% names(FU)) {
  ALL <- dplyr::left_join(ALL, FU |> dplyr::select(patient, CPE),
                          by = "patient")
}

# ==============================================================================
# 4. The association, by tier
# ==============================================================================
hdr("4. Promoter methylation and expression")

promoter_mean <- function(d, g) {
  p <- intersect(PRF$probe[PRF$gene == g & PRF$promoter], names(d))
  if (!length(p)) return(rep(NA_real_, nrow(d)))
  rowMeans(as.matrix(d[, p, drop = FALSE]), na.rm = TRUE)
}

run <- function(d, g, label) {
  ycol <- paste0(g, "_log2tpm")
  if (!ycol %in% names(d)) return(NULL)
  d$.beta <- promoter_mean(d, g)
  d$.y <- suppressWarnings(as.numeric(d[[ycol]]))
  s <- d[is.finite(d$.beta) & is.finite(d$.y), ]
  if (nrow(s) < 20) {
    log_line("  ", g, " / ", label, ": only ", nrow(s), "; skipped",
             logfile = LOG)
    return(NULL)
  }
  ct <- suppressWarnings(stats::cor.test(s$.beta, s$.y, method = "spearman"))
  # Covariate selection. A covariate is only usable if it varies AND is
  # populated: CPE purity came from the cohort files, so the samples those
  # files never contained have none, and lm would silently drop them as
  # incomplete cases. That produced an "n = 469" adjusted estimate that was
  # really the 411-sample model under a wrong label. A covariate that costs
  # more than 5% of the rows is therefore dropped, with the loss reported.
  cov <- intersect(c("CPE", "sample_type"), names(s))
  cov <- cov[vapply(cov, function(cn)
    sum(!is.na(s[[cn]])) > 20 &&
      dplyr::n_distinct(s[[cn]][!is.na(s[[cn]])]) > 1, logical(1))]
  dropped <- character(0)
  for (cn in cov) {
    lost <- sum(is.na(s[[cn]]))
    if (lost > 0.05 * nrow(s)) {
      dropped <- c(dropped, sprintf("%s (missing in %d)", cn, lost))
      cov <- setdiff(cov, cn)
    }
  }
  fml <- stats::as.formula(paste(".y ~ .beta",
                                 if (length(cov)) paste("+", paste(cov, collapse = " + ")) else ""))
  fit <- try(stats::lm(fml, data = s), silent = TRUE)
  b <- lo <- hi <- pv <- NA_real_; n_mod <- NA_integer_
  if (!inherits(fit, "try-error")) {
    n_mod <- stats::nobs(fit)
    cf <- summary(fit)$coefficients
    if (".beta" %in% rownames(cf)) {
      b <- cf[".beta", 1]; pv <- cf[".beta", 4]
      ci <- try(stats::confint(fit, ".beta"), silent = TRUE)
      if (!inherits(ci, "try-error")) { lo <- ci[1]; hi <- ci[2] }
    }
  }
  if (length(dropped)) {
    log_line("           dropped covariate: ", paste(dropped, collapse = "; "),
             " -- would have silently removed those rows", logfile = LOG)
  }
  if (!is.na(n_mod) && n_mod < nrow(s)) {
    log_line(sprintf("           [FLAG] adjusted model used %d of %d rows",
                     n_mod, nrow(s)), logfile = LOG)
  }
  log_line(sprintf("  %-8s %-34s n = %3d | rho = %+.3f, p = %-10s | adj beta = %+.3f [%.3f, %.3f], p = %-10s (model n = %s)",
                   g, label, nrow(s), ct$estimate, fmt_p(ct$p.value), b, lo, hi,
                   fmt_p(pv), ifelse(is.na(n_mod), "NA", as.character(n_mod))),
           logfile = LOG)
  tibble::tibble(gene = g, tier = label, n = nrow(s), n_model = n_mod,
                 rho = unname(ct$estimate),
                 rho_p = ct$p.value, adj_beta = b, adj_lo = lo, adj_hi = hi,
                 adj_p = pv, covariates = paste(cov, collapse = " + "),
                 covariates_dropped = paste(dropped, collapse = "; "))
}

tiers <- list()
for (g in GENES) {
  sets <- list(list(ALL, "T4 COMPLETE (no filters)"))
  if (!is.null(FU)) sets <- c(sets,
    list(list(ALL[ALL$patient %in% FU$patient, ], "T2 all cutaneous")))
  if (!is.null(AN)) sets <- c(sets,
    list(list(ALL[ALL$patient %in% AN$patient, ], "T1 UV-positive (original)")))
  if (nrow(ADDED)) sets <- c(sets, list(list(ADDED, "Added samples only")))
  for (ss in sets) {
    r <- run(ss[[1]], g, ss[[2]])
    if (!is.null(r)) tiers[[length(tiers) + 1L]] <- r
  }
}
TI <- dplyr::bind_rows(tiers)
if (nrow(TI)) readr::write_csv(TI, file.path(RES, "methylation_full_tiers.csv"))

log_line("", logfile = LOG)
log_line("  PRE-SPECIFIED EXPECTATION: the complete-cohort rho for RASSF9 lands",
         logfile = LOG)
log_line("  near -0.37, between the 411-sample -0.374 and the no-UV-signature",
         logfile = LOG)
log_line("  subset's -0.409. A materially different value means the added",
         logfile = LOG)
log_line("  samples differ systematically -- see the next section.",
         logfile = LOG)

# ==============================================================================
# 5. Do the added samples differ?
# ==============================================================================
hdr("5. Characterising the added samples")

if (nrow(ADDED)) {
  base <- if (!is.null(FU)) ALL[ALL$patient %in% FU$patient, ] else ALL
  cmp <- list()
  for (g in GENES) {
    yc <- paste0(g, "_log2tpm")
    if (!yc %in% names(ALL)) next
    cmp[[length(cmp) + 1L]] <- tibble::tibble(
      measure = paste0(g, " log2(TPM+1), median"),
      in_cohort = stats::median(base[[yc]], na.rm = TRUE),
      added = stats::median(ADDED[[yc]], na.rm = TRUE))
    bb <- promoter_mean(base, g); ba <- promoter_mean(ADDED, g)
    cmp[[length(cmp) + 1L]] <- tibble::tibble(
      measure = paste0(g, " promoter beta, median"),
      in_cohort = stats::median(bb, na.rm = TRUE),
      added = stats::median(ba, na.rm = TRUE))
  }
  CMP <- dplyr::bind_rows(cmp)
  for (i in seq_len(nrow(CMP))) {
    log_line(sprintf("  %-34s in cohort %8.3f | added %8.3f", CMP$measure[i],
                     CMP$in_cohort[i], CMP$added[i]), logfile = LOG)
  }
  if ("sample_type" %in% names(ALL)) {
    log_line("", logfile = LOG)
    log_line("  sample type, in cohort: ",
             paste(names(table(base$sample_type)), as.integer(table(base$sample_type)),
                   sep = "=", collapse = ", "), logfile = LOG)
    log_line("  sample type, added:     ",
             paste(names(table(ADDED$sample_type)), as.integer(table(ADDED$sample_type)),
                   sep = "=", collapse = ", "), logfile = LOG)
  }
  readr::write_csv(ADDED |> dplyr::select(dplyr::any_of(
    c("patient", "sample_type", paste0(GENES, "_log2tpm")))),
    file.path(RES, "methylation_full_addedsamples.csv"))
  log_line("", logfile = LOG)
  log_line("  These samples were absent from the cohort files because those",
           logfile = LOG)
  log_line("  were built by joining to the mutational signature table, not",
           logfile = LOG)
  log_line("  because anything was wrong with them.", logfile = LOG)
}

# ==============================================================================
# 6. Specificity, with the full probe set
# ==============================================================================
hdr("6. Specificity across the whole locus")

pr_rows <- list()
for (g in GENES) {
  yc <- paste0(g, "_log2tpm")
  if (!yc %in% names(ALL)) next
  ps <- PRF[PRF$gene == g, ]
  y <- suppressWarnings(as.numeric(ALL[[yc]]))
  for (i in seq_len(nrow(ps))) {
    if (!ps$probe[i] %in% names(ALL)) next
    v <- suppressWarnings(as.numeric(ALL[[ps$probe[i]]]))
    ok <- is.finite(v) & is.finite(y)
    if (sum(ok) < 20) next
    ct <- suppressWarnings(stats::cor.test(v[ok], y[ok], method = "spearman"))
    pr_rows[[length(pr_rows) + 1L]] <- tibble::tibble(
      gene = g, probe = ps$probe[i], dist_tss = ps$dist_tss[i],
      region = ps$region[i], promoter = ps$promoter[i], n = sum(ok),
      rho = unname(ct$estimate), p = ct$p.value,
      median_beta = stats::median(v[ok]))
  }
}
PB <- dplyr::bind_rows(pr_rows)
if (nrow(PB)) {
  PB <- PB |> dplyr::group_by(gene) |>
    dplyr::mutate(p_bh = stats::p.adjust(p, method = "BH")) |>
    dplyr::ungroup()
  readr::write_csv(PB, file.path(RES, "methylation_full_probes.csv"))
  for (g in GENES) {
    sub <- PB |> dplyr::filter(gene == g) |> dplyr::arrange(dist_tss)
    if (!nrow(sub)) next
    log_line("", logfile = LOG)
    log_line("  ", g, ", every probe within ", WINDOW, " bp:", logfile = LOG)
    log_line(sprintf("    %-12s %9s %-11s %9s %8s %10s", "probe", "TSS dist",
                     "region", "med beta", "rho", "BH p"), logfile = LOG)
    for (i in seq_len(nrow(sub))) {
      log_line(sprintf("    %-12s %+9d %-11s %9.3f %+8.3f %10s", sub$probe[i],
                       sub$dist_tss[i], sub$region[i], sub$median_beta[i],
                       sub$rho[i], fmt_p(sub$p_bh[i])), logfile = LOG)
    }
    pr <- sub |> dplyr::filter(promoter); nq <- sub |> dplyr::filter(!promoter)
    if (nrow(pr) && nrow(nq)) {
      log_line(sprintf("    promoter:     mean rho %+.3f (%d/%d significant)",
                       mean(pr$rho), sum(pr$p_bh < 0.05), nrow(pr)),
               logfile = LOG)
      log_line(sprintf("    non-promoter: mean rho %+.3f (%d/%d significant)",
                       mean(nq$rho), sum(nq$p_bh < 0.05), nrow(nq)),
               logfile = LOG)
      if (sign(mean(pr$rho)) != sign(mean(nq$rho))) {
        log_line("    -> OPPOSITE SIGNS: promoter-specific, the expected",
                 logfile = LOG)
        log_line("       promoter-versus-gene-body contrast", logfile = LOG)
      } else if (abs(mean(nq$rho)) > 0.7 * abs(mean(pr$rho))) {
        log_line("    -> SIMILAR: methylation anywhere at the locus tracks",
                 logfile = LOG)
        log_line("       expression about as well, so the silencing reading",
                 logfile = LOG)
        log_line("       weakens considerably", logfile = LOG)
      } else {
        log_line("    -> promoter probes correlate substantially more strongly",
                 logfile = LOG)
      }
    }
  }
}

# ==============================================================================
# 7. Figure
# ==============================================================================
hdr("7. Figure")

reset_legends("# TCGA methylation (complete cohort) figure legends\n")

for (g in GENES) {
  sub <- TI |> dplyr::filter(gene == g)
  pb  <- PB |> dplyr::filter(gene == g)
  if (!nrow(sub)) next

  pa <- ggplot(sub, aes(x = rho, y = factor(tier, levels = rev(tier)))) +
    geom_vline(xintercept = 0, linetype = "22", colour = COL_GREY,
               linewidth = 0.4) +
    geom_point(aes(size = n), colour = COL_B) +
    scale_size_continuous(range = c(1.8, 3.6), guide = "none") +
    labs(x = "Spearman rho", y = NULL) +
    theme_journal() +
    theme(panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey92"),
          plot.margin = margin(4, 10, 4, 4, unit = "pt"))

  plots <- list(pa)
  if (nrow(pb)) {
    pbp <- ggplot(pb, aes(x = dist_tss, y = rho)) +
      annotate("rect", xmin = PROM_LO, xmax = PROM_HI, ymin = -Inf, ymax = Inf,
               fill = COL_GREY, alpha = 0.12) +
      geom_hline(yintercept = 0, linetype = "22", colour = COL_GREY,
                 linewidth = 0.4) +
      geom_point(aes(colour = promoter, shape = p_bh < 0.05), size = 1.8) +
      scale_colour_manual(values = c(`TRUE` = COL_B, `FALSE` = COL_A),
                          labels = c(`TRUE` = "Promoter", `FALSE` = "Elsewhere"),
                          name = NULL) +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                         labels = c(`TRUE` = "BH p < 0.05",
                                    `FALSE` = "not significant"), name = NULL) +
      labs(x = "Distance from TSS (bp)", y = "Spearman rho") +
      theme_journal() +
      theme(legend.position = "bottom",
            legend.text = element_text(size = FIG$base_size - 3),
            legend.key.size = unit(8, "pt"),
            plot.margin = margin(4, 10, 4, 4, unit = "pt"))
    plots <- list(pa, pbp)
  }
  fig <- patchwork::wrap_plots(plots, nrow = 1) +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = FIG$base_size + 2, face = "bold",
                                  colour = COL_INK))
  stem <- paste0("tcga_", tolower(g), "_methylation_full")
  save_fig_to(fig, stem, FIG_TCGA, width = FIG$width_double, height = 88)
  log_line("  wrote: ", stem, logfile = LOG)

  t4 <- sub |> dplyr::filter(grepl("^T4", tier))
  add_legend(stem,
    paste0("Figure. ", g,
           " promoter methylation and expression in the complete TCGA-SKCM cohort."),
    paste0("(A) Spearman correlation between mean promoter methylation and ", g,
           " expression across cohort definitions, rebuilt from the raw ",
           "expression and methylation objects so that no mutational-signature ",
           "requirement applies; point size is sample number. (B) Per-probe ",
           "correlation against distance from the transcription start site for ",
           "every HM450 probe within ", format(WINDOW, big.mark = ","),
           " bp, with the promoter window shaded; filled points are significant ",
           "after Benjamini-Hochberg correction. ",
           if (nrow(t4)) paste0("In the complete cohort (n = ", t4$n[1],
                                "), rho = ", sprintf("%.3f", t4$rho[1]),
                                ", p = ", fmt_p(t4$rho_p[1]), ".") else ""))
}

report_warnings(logfile = LOG)
write_session_info("35_methylation_full_cohort")
log_line("\nFull log: ", LOG, logfile = NULL)
