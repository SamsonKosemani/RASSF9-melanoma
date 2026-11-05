# ==============================================================================
# 36_methylation_specificity.R
#
# Is the RASSF9 promoter result promoter-specific, or a sample-level artefact?
#
# WHY THIS SCRIPT EXISTS. Script 35 widened the probe window to 150 kb and found
# five non-promoter RASSF9 probes that correlate with expression in the SAME
# direction as the promoter (mean rho -0.144 against the promoter's -0.243, with
# one probe 112 kb upstream reaching rho = -0.214, p = 6.8e-06). RASSF10 inverts
# cleanly; RASSF9 does not. That is the pattern pre-registered in script 34 as
# the one that would UNDERCUT the silencing reading:
#
#   "every probe on the locus correlates equally, which would point at a
#    sample-level property (purity, cell composition) rather than promoter
#    silencing"
#
# It is not "equally" -- the promoter is about 1.7x stronger and its five
# TSS-proximal probes (-0.294 to -0.485) beat any distal probe -- but the
# margin rests on a threshold I chose, and the honest position is that this is
# unresolved rather than passed.
#
# THE SPECIFIC ALTERNATIVE. Our own Visium result found only 17.1% of RASSF9
# molecules in melanocytic spots; the gene is largely stromal, endothelial and
# eccrine. In a bulk tumour, higher purity therefore means less RASSF9 by
# composition alone. Purity also raises tumour methylation genome-wide. The two
# together produce a negative expression-methylation correlation at ANY probe,
# with no silencing in the mechanism. Distal probes behaving like the promoter
# is exactly that confound's signature.
#
# TWO TESTS.
#   (1) ADJUSTMENT. A per-sample global methylation index from 20,000 random
#       probes, plus CPE purity where available. Does the promoter effect
#       survive?
#   (2) EMPIRICAL NULL. Correlate RASSF9 expression against 2,000 random probes
#       genome-wide. Where do its own promoter probes fall in that
#       distribution? Extreme tail = locus-specific. Bulk = sample-level.
#
# CRITERION, FIXED BEFORE THE NUMBERS:
#   PROMOTER-SPECIFIC  = adjusted effect keeps >60% of its unadjusted magnitude
#                        at p < 0.05, AND the TSS-proximal probes fall in the
#                        extreme 5% of the null
#   GLOBAL CONFOUND    = adjusted effect loses >60% of its magnitude, OR the
#                        probes sit inside the null's bulk (outside the 10% tail)
#   PARTIAL            = anything in between; reported as unresolved, not as
#                        support
#
# A note on what "survives" means. The global index is itself partly driven by
# purity, so adjusting for it is conservative: it removes real promoter signal
# along with the confound. An effect that survives this is therefore strong
# evidence; one that does not is NOT proof of artefact, only failure to exclude
# it. That asymmetry is stated here so it cannot be read either way later.
#
# Output: results/tcga/methylation_specificity.csv
#         results/tcga/methylation_null_distribution.csv
#         figures/tcga/tcga_rassf9_specificity.{pdf,png,svg}
#         results/logs/36_methylation_specificity.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

LOG <- file.path(PATHS$logs, "36_methylation_specificity.log")
cat("RASSF9 methylation specificity —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

PROC <- PATHS$processed
RES  <- file.path(PATHS$results, "tcga")
FIG_TCGA <- file.path(PATHS$figures, "tcga")
invisible(lapply(c(RES, FIG_TCGA), dir.create, recursive = TRUE,
                 showWarnings = FALSE))

GENES    <- c("RASSF9", "RASSF10")
WINDOW   <- 150000
PROM_LO  <- -1500
PROM_HI  <- 500
N_GLOBAL <- 20000     # random probes for the global methylation index
N_NULL   <- 2000      # random probes for the empirical null
SEED     <- 20261016
RETAIN   <- 0.60      # pre-set: fraction of effect that must survive
TAIL     <- 0.05      # pre-set: null tail defining "locus-specific"

set.seed(SEED)
log_line("  seed: ", SEED, " | global index probes: ", N_GLOBAL,
         " | null probes: ", N_NULL, logfile = LOG)

pick <- function(slot, fallback) {
  if (!exists("PAL", inherits = TRUE)) return(fallback)
  x <- get("PAL", inherits = TRUE)[[slot]]
  if (!is.null(x) && is.character(x) && nzchar(x[1])) x[1] else fallback
}
COL_A <- pick("nevus", "#0095D9"); COL_B <- pick("melanoma", "#F2680C")
COL_GREY <- pick("grey_pt", "#8E8E8E"); COL_INK <- pick("ink", "#1A1A1A")

RNA_F  <- file.path(PATHS$raw, "skcm_star_counts.rds")
METH_F <- file.path(PATHS$raw, "skcm_hm450.rds")

# ==============================================================================
# 1. Expression and gene coordinates
# ==============================================================================
hdr("1. Expression")

SE <- readRDS(RNA_F)
a_tpm <- intersect(c("tpm_unstrand", "tpm_unstranded", "TPM"), assayNames(SE))[1]
rd <- as.data.frame(rowData(SE))
gcol <- intersect(c("gene_name", "external_gene_name", "symbol"), names(rd))[1]
GR <- rowRanges(SE)

loci <- list()
EXPR <- list()
for (g in GENES) {
  i <- which(rd[[gcol]] == g)
  if (!length(i)) { log_line("  [FLAG] ", g, " absent", logfile = LOG); next }
  EXPR[[g]] <- log2(as.numeric(
    SummarizedExperiment::assay(SE, a_tpm)[i[1], ]) + 1)
  gr <- GR[i[1]]
  loci[[g]] <- list(
    chr = as.character(GenomicRanges::seqnames(gr)),
    strand = as.character(BiocGenerics::strand(gr)),
    tss = if (as.character(BiocGenerics::strand(gr)) == "-")
      BiocGenerics::end(gr) else BiocGenerics::start(gr))
}

cd <- as.data.frame(colData(SE))
stcol <- intersect(c("sample_type", "shortLetterCode", "definition"), names(cd))
RNA <- tibble::tibble(
  patient = substr(colnames(SE), 1, 12),
  sample_type = if (length(stcol)) as.character(cd[[stcol[1]]]) else NA_character_)
for (g in names(EXPR)) RNA[[paste0(g, "_log2tpm")]] <- EXPR[[g]]
RNA <- RNA |> dplyr::distinct(patient, .keep_all = TRUE)
log_line("  expression patients: ", nrow(RNA), logfile = LOG)
rm(SE); gc(verbose = FALSE)

# ==============================================================================
# 2. Methylation: locus probes, a global index, and a null set
# ==============================================================================
hdr("2. Methylation")

ME <- readRDS(METH_F)
MGR <- rowRanges(ME)
pn <- rownames(ME)
log_line("  probes on array: ", length(pn), " | samples: ", ncol(ME),
         logfile = LOG)

PRF <- list()
for (g in names(loci)) {
  L <- loci[[g]]
  sel <- which(as.character(GenomicRanges::seqnames(MGR)) == L$chr &
               abs(BiocGenerics::start(MGR) - L$tss) <= WINDOW)
  pos <- BiocGenerics::start(MGR[sel])
  d <- if (L$strand == "-") L$tss - pos else pos - L$tss
  PRF[[g]] <- tibble::tibble(
    gene = g, probe = pn[sel], dist_tss = d,
    promoter = d >= PROM_LO & d <= PROM_HI,
    tss_proximal = d >= -400 & d <= PROM_HI)
  log_line("  ", g, ": ", length(sel), " probes, ",
           sum(PRF[[g]]$promoter), " promoter, ",
           sum(PRF[[g]]$tss_proximal), " TSS-proximal", logfile = LOG)
}
PRF <- dplyr::bind_rows(PRF)
locus_probes <- unique(PRF$probe)

# The global index and the null set must exclude the loci under test, or the
# comparison is partly against itself.
excl <- unique(unlist(lapply(names(loci), function(g) {
  L <- loci[[g]]
  pn[as.character(GenomicRanges::seqnames(MGR)) == L$chr &
     abs(BiocGenerics::start(MGR) - L$tss) <= 1e6]
})))
pool <- setdiff(pn, excl)
log_line("  pool for index/null (loci +/-1 Mb removed): ", length(pool),
         logfile = LOG)

g_probes <- sample(pool, min(N_GLOBAL, length(pool)))
n_probes <- sample(setdiff(pool, g_probes), min(N_NULL, length(pool) - length(g_probes)))

A <- SummarizedExperiment::assay(ME)
GB <- colMeans(as.matrix(A[g_probes, , drop = FALSE]), na.rm = TRUE)
BL <- as.matrix(A[locus_probes, , drop = FALSE])
BN <- as.matrix(A[n_probes, , drop = FALSE])
meth_pat <- substr(colnames(ME), 1, 12)
rm(ME, A); gc(verbose = FALSE)

log_line("  global methylation index: median ", sprintf("%.3f", stats::median(GB, na.rm = TRUE)),
         ", range ", sprintf("%.3f-%.3f", min(GB, na.rm = TRUE), max(GB, na.rm = TRUE)),
         logfile = LOG)

keep1 <- !duplicated(meth_pat)
MET <- tibble::tibble(patient = meth_pat[keep1], global_beta = GB[keep1])
MET <- dplyr::bind_cols(MET, tibble::as_tibble(t(BL[, keep1, drop = FALSE])))
NUL <- t(BN[, keep1, drop = FALSE]); rownames(NUL) <- meth_pat[keep1]

D <- dplyr::inner_join(RNA, MET, by = "patient")
f_full <- file.path(PROC, "full_cohort.csv")
if (file.exists(f_full)) {
  FU <- readr::read_csv(f_full, show_col_types = FALSE, progress = FALSE)
  if ("CPE" %in% names(FU))
    D <- dplyr::left_join(D, FU |> dplyr::select(patient, CPE), by = "patient")
}
NUL <- NUL[match(D$patient, rownames(NUL)), , drop = FALSE]
log_line("  matched patients: ", nrow(D), " | with CPE: ",
         if ("CPE" %in% names(D)) sum(!is.na(D$CPE)) else 0, logfile = LOG)

# ==============================================================================
# 3. Does purity behave as the confound requires?
# ==============================================================================
hdr("3. The confound's own premises")

log_line("  For the purity explanation to work, three things must hold:",
         logfile = LOG)
log_line("    (a) purity correlates NEGATIVELY with RASSF9 expression",
         logfile = LOG)
log_line("    (b) purity correlates POSITIVELY with global methylation",
         logfile = LOG)
log_line("    (c) global methylation correlates negatively with expression",
         logfile = LOG)
log_line("  If any fails, the confound cannot generate the observed pattern.",
         logfile = LOG)
log_line("", logfile = LOG)

prem <- list()
show_cor <- function(x, y, lab) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 20) { log_line("    ", lab, ": too few"); return(NULL) }
  ct <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman"))
  log_line(sprintf("    %-46s n = %3d  rho = %+.3f  p = %s", lab, sum(ok),
                   ct$estimate, fmt_p(ct$p.value)), logfile = LOG)
  tibble::tibble(comparison = lab, n = sum(ok), rho = unname(ct$estimate),
                 p = ct$p.value)
}
if ("CPE" %in% names(D)) {
  for (g in names(loci)) {
    prem[[length(prem) + 1L]] <- show_cor(D$CPE, D[[paste0(g, "_log2tpm")]],
      paste0("(a) CPE purity vs ", g, " expression"))
  }
  prem[[length(prem) + 1L]] <- show_cor(D$CPE, D$global_beta,
    "(b) CPE purity vs global methylation")
}
for (g in names(loci)) {
  prem[[length(prem) + 1L]] <- show_cor(D$global_beta, D[[paste0(g, "_log2tpm")]],
    paste0("(c) global methylation vs ", g, " expression"))
}
PREM <- dplyr::bind_rows(prem)

# ==============================================================================
# 4. Nested adjustment
# ==============================================================================
hdr("4. Does the promoter effect survive adjustment?")

mean_beta <- function(probes) {
  probes <- intersect(probes, names(D))
  if (!length(probes)) return(rep(NA_real_, nrow(D)))
  rowMeans(as.matrix(D[, probes, drop = FALSE]), na.rm = TRUE)
}

fit_one <- function(y, beta, extra, lab, gene, set) {
  dd <- data.frame(y = y, beta = beta)
  for (nm in names(extra)) dd[[nm]] <- extra[[nm]]
  dd <- dd[stats::complete.cases(dd), , drop = FALSE]
  if (nrow(dd) < 30 || dplyr::n_distinct(dd$beta) < 5) return(NULL)
  keep <- vapply(setdiff(names(dd), c("y", "beta")),
                 function(cn) dplyr::n_distinct(dd[[cn]]) > 1, logical(1))
  use <- c("beta", names(keep)[keep])
  fit <- try(stats::lm(stats::as.formula(paste("y ~", paste(use, collapse = " + "))),
                       data = dd), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  cf <- summary(fit)$coefficients
  if (!"beta" %in% rownames(cf)) return(NULL)
  ci <- try(stats::confint(fit, "beta"), silent = TRUE)
  tibble::tibble(gene = gene, probe_set = set, model = lab, n = nrow(dd),
                 beta = cf["beta", 1],
                 lo = if (inherits(ci, "try-error")) NA_real_ else ci[1],
                 hi = if (inherits(ci, "try-error")) NA_real_ else ci[2],
                 p = cf["beta", 4])
}

rows <- list()
for (g in names(loci)) {
  y <- D[[paste0(g, "_log2tpm")]]
  sets <- list(
    promoter = PRF$probe[PRF$gene == g & PRF$promoter],
    tss_proximal = PRF$probe[PRF$gene == g & PRF$tss_proximal],
    distal = PRF$probe[PRF$gene == g & !PRF$promoter])
  for (sn in names(sets)) {
    b <- mean_beta(sets[[sn]])
    if (all(is.na(b))) next
    mods <- list(
      list("M0 unadjusted", list()),
      list("M1 + sample type", list(sample_type = D$sample_type)),
      list("M2 + global methylation", list(global_beta = D$global_beta)),
      list("M3 + global + sample type",
           list(global_beta = D$global_beta, sample_type = D$sample_type)))
    if ("CPE" %in% names(D)) mods <- c(mods, list(
      list("M4 + global + type + purity",
           list(global_beta = D$global_beta, sample_type = D$sample_type,
                CPE = D$CPE))))
    for (m in mods) {
      r <- fit_one(y, b, m[[2]], m[[1]], g, sn)
      if (!is.null(r)) rows[[length(rows) + 1L]] <- r
    }
  }
}
AD <- dplyr::bind_rows(rows)

for (g in names(loci)) {
  for (sn in c("promoter", "tss_proximal", "distal")) {
    s <- AD |> dplyr::filter(gene == g, probe_set == sn)
    if (!nrow(s)) next
    b0 <- s$beta[s$model == "M0 unadjusted"]
    log_line("", logfile = LOG)
    log_line("  ", g, " / ", sn, ":", logfile = LOG)
    for (i in seq_len(nrow(s))) {
      # A ratio against a near-zero baseline is noise, not a percentage.
      ret <- if (length(b0) && is.finite(b0) && abs(b0) >= 0.2)
        sprintf("%3.0f%% retained", 100 * s$beta[i] / b0) else
        "  n/a (baseline ~ 0)"
      log_line(sprintf("    %-30s n = %3d  beta = %+.3f [%+.3f, %+.3f]  p = %-10s %s",
                       s$model[i], s$n[i], s$beta[i], s$lo[i], s$hi[i],
                       fmt_p(s$p[i]), ret), logfile = LOG)
    }
  }
}
if (nrow(AD)) readr::write_csv(AD, file.path(RES, "methylation_specificity.csv"))

# ==============================================================================
# 5. Empirical null
# ==============================================================================
hdr("5. Empirical null: 2,000 random probes genome-wide")

log_line("  A locus-specific effect puts RASSF9's own promoter probes in the",
         logfile = LOG)
log_line("  extreme tail of this distribution. A sample-level property puts",
         logfile = LOG)
log_line("  them in the bulk, because then any probe would do as well.",
         logfile = LOG)

null_rows <- list()
for (g in names(loci)) {
  y <- D[[paste0(g, "_log2tpm")]]
  nr <- apply(NUL, 2, function(v) {
    ok <- is.finite(v) & is.finite(y)
    if (sum(ok) < 50 || stats::sd(v[ok]) == 0) return(NA_real_)
    suppressWarnings(stats::cor(v[ok], y[ok], method = "spearman"))
  })
  nr <- nr[is.finite(nr)]
  if (!length(nr)) next

  own <- PRF |> dplyr::filter(gene == g) |>
    dplyr::rowwise() |>
    dplyr::mutate(rho = {
      v <- if (probe %in% names(D)) suppressWarnings(as.numeric(D[[probe]])) else NA_real_
      ok <- is.finite(v) & is.finite(y)
      if (sum(ok) < 50) NA_real_ else
        suppressWarnings(stats::cor(v[ok], y[ok], method = "spearman"))
    }) |> dplyr::ungroup() |> dplyr::filter(is.finite(rho))

  q_lo <- stats::quantile(nr, TAIL, names = FALSE)
  q_lo10 <- stats::quantile(nr, 0.10, names = FALSE)
  log_line("", logfile = LOG)
  log_line(sprintf("  %s null: n = %d probes, median rho %+.3f, %d%%ile %+.3f, 10%%ile %+.3f",
                   g, length(nr), stats::median(nr), round(TAIL * 100), q_lo,
                   q_lo10), logfile = LOG)
  log_line(sprintf("    %-13s %9s %-13s %8s %10s", "probe", "TSS dist", "set",
                   "rho", "null pctile"), logfile = LOG)
  for (i in seq_len(nrow(own))) {
    pct <- 100 * mean(nr <= own$rho[i])
    setl <- if (own$tss_proximal[i]) "TSS-proximal" else
      if (own$promoter[i]) "promoter" else "distal"
    log_line(sprintf("    %-13s %+9d %-13s %+8.3f %9.2f%%", own$probe[i],
                     own$dist_tss[i], setl, own$rho[i], pct), logfile = LOG)
    null_rows[[length(null_rows) + 1L]] <- tibble::tibble(
      gene = g, probe = own$probe[i], dist_tss = own$dist_tss[i], set = setl,
      rho = own$rho[i], null_pctile = pct, null_n = length(nr),
      null_median = stats::median(nr), null_p05 = q_lo)
  }
  tp <- own |> dplyr::filter(tss_proximal)
  if (nrow(tp)) {
    inside <- mean(tp$rho <= q_lo)
    log_line(sprintf("    -> %d of %d TSS-proximal probes fall in the extreme %d%% tail",
                     sum(tp$rho <= q_lo), nrow(tp), round(TAIL * 100)),
             logfile = LOG)
  }
  ds <- own |> dplyr::filter(!promoter)
  if (nrow(ds)) {
    log_line(sprintf("    -> %d of %d distal probes also fall in that tail",
                     sum(ds$rho <= q_lo), nrow(ds)), logfile = LOG)
    log_line("       (distal probes in the tail argue the locus matters;",
             logfile = LOG)
    log_line("        distal probes in the bulk argue only the promoter does)",
             logfile = LOG)
  }
  if (g == "RASSF9") assign(".null_rassf9", nr)
}
ND <- dplyr::bind_rows(null_rows)
if (nrow(ND)) readr::write_csv(ND, file.path(RES, "methylation_null_distribution.csv"))

# ==============================================================================
# 6. Verdict
# ==============================================================================
hdr("6. Verdict against the pre-set criterion")

log_line("  PROMOTER-SPECIFIC = adjusted effect keeps >", round(RETAIN * 100),
         "% at p < 0.05 AND", logfile = LOG)
log_line("                      TSS-proximal probes in the extreme ",
         round(TAIL * 100), "% of the null", logfile = LOG)
log_line("  GLOBAL CONFOUND   = keeps <", round(RETAIN * 100),
         "%, OR probes outside the 10% tail", logfile = LOG)
log_line("  PARTIAL           = in between; unresolved, NOT support",
         logfile = LOG)
log_line("", logfile = LOG)

for (g in names(loci)) {
  s <- AD |> dplyr::filter(gene == g, probe_set == "promoter")
  if (!nrow(s)) next
  b0 <- s$beta[s$model == "M0 unadjusted"]
  mfull <- s |> dplyr::filter(grepl("^M3|^M4", model)) |> dplyr::slice_tail(n = 1)
  if (!nrow(mfull) || !length(b0) || !is.finite(b0)) next
  if (abs(b0) < 0.2) {
    log_line("  ", g, ": unadjusted slope ~ 0; retention undefined, no verdict",
             logfile = LOG)
    next
  }
  ret <- mfull$beta / b0
  tp <- ND |> dplyr::filter(gene == g, set == "TSS-proximal")
  in_tail <- if (nrow(tp)) all(tp$null_pctile <= TAIL * 100) else NA
  in_bulk <- if (nrow(tp)) any(tp$null_pctile > 10) else NA

  verdict <- if (isTRUE(ret >= RETAIN && mfull$p < 0.05 && isTRUE(in_tail))) {
    "PROMOTER-SPECIFIC"
  } else if (isTRUE(ret < RETAIN) || isTRUE(in_bulk)) {
    "GLOBAL CONFOUND NOT EXCLUDED"
  } else "PARTIAL -- unresolved"

  log_line("  ", g, ": ", verdict, logfile = LOG)
  log_line(sprintf("    unadjusted beta %+.3f -> fully adjusted %+.3f (%.0f%% retained, p = %s)",
                   b0, mfull$beta, 100 * ret, fmt_p(mfull$p)), logfile = LOG)
  if (nrow(tp)) {
    log_line(sprintf("    TSS-proximal null percentiles: %s",
                     paste(sprintf("%.2f%%", tp$null_pctile), collapse = ", ")),
             logfile = LOG)
  }
  log_line("", logfile = LOG)
}

log_line("  ASYMMETRY, stated so it cannot be read both ways later: the global", logfile = LOG)
log_line("  index is itself partly driven by purity, so adjusting for it also", logfile = LOG)
log_line("  removes real promoter signal. Surviving is therefore strong", logfile = LOG)
log_line("  evidence; not surviving is failure to EXCLUDE the confound, not", logfile = LOG)
log_line("  proof of one.", logfile = LOG)

# ==============================================================================
# 7. Figure
# ==============================================================================
hdr("7. Figure")

reset_legends("# TCGA methylation specificity figure legends\n")

if (exists(".null_rassf9")) {
  nr <- get(".null_rassf9")
  own <- ND |> dplyr::filter(gene == "RASSF9")
  pa <- ggplot(data.frame(rho = nr), aes(x = rho)) +
    geom_histogram(bins = 60, fill = COL_GREY, alpha = 0.45,
                   colour = NA) +
    geom_vline(data = own, aes(xintercept = rho, colour = set),
               linewidth = 0.5) +
    scale_colour_manual(values = c(`TSS-proximal` = COL_B,
                                   promoter = "#8C2D04", distal = COL_A),
                        name = NULL) +
    labs(x = "Spearman rho with RASSF9 expression",
         y = "Random probes") +
    theme_journal() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = FIG$base_size - 3),
          legend.key.size = unit(8, "pt"),
          plot.margin = margin(4, 12, 4, 4, unit = "pt"))

  sub <- AD |> dplyr::filter(gene == "RASSF9", probe_set == "promoter")
  pb <- ggplot(sub, aes(x = beta, y = factor(model, levels = rev(model)))) +
    geom_vline(xintercept = 0, linetype = "22", colour = COL_GREY,
               linewidth = 0.4) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0,
                  linewidth = 0.4, colour = COL_INK) +
    geom_point(size = 1.9, colour = COL_B) +
    labs(x = "Adjusted slope, promoter methylation", y = NULL) +
    theme_journal() +
    theme(plot.margin = margin(4, 12, 4, 4, unit = "pt"))

  fig <- patchwork::wrap_plots(list(pa, pb), nrow = 1) +
    patchwork::plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = FIG$base_size + 2, face = "bold",
                                  colour = COL_INK))
  save_fig_to(fig, "tcga_rassf9_specificity", FIG_TCGA,
              width = FIG$width_double, height = 90)
  log_line("  wrote: tcga_rassf9_specificity", logfile = LOG)

  add_legend("tcga_rassf9_specificity",
    "Figure. RASSF9 promoter methylation is tested against a genome-wide null and against global methylation.",
    paste0("(A) Distribution of Spearman correlations between RASSF9 expression and ",
           format(length(nr), big.mark = ","), " randomly chosen HM450 probes drawn from ",
           "outside the RASSF9 and RASSF10 loci (both excluded with a 1 Mb margin); ",
           "vertical lines mark RASSF9's own probes, separated into TSS-proximal, ",
           "remaining promoter, and distal. A locus-specific association places the ",
           "TSS-proximal probes in the extreme tail of this distribution. ",
           "(B) Slope of RASSF9 expression on mean promoter methylation under ",
           "progressive adjustment for sample type, a global methylation index built ",
           "from ", format(N_GLOBAL, big.mark = ","), " random probes, and consensus ",
           "purity; bars are 95% confidence intervals. Adjusting for the global index ",
           "is conservative, since that index is itself partly driven by purity."))
}

report_warnings(logfile = LOG)
write_session_info("36_methylation_specificity")
log_line("\nFull log: ", LOG, logfile = NULL)
