# ==============================================================================
# 39a_burden_inspect.R
#
# INSPECT ONLY. No analysis, no conclusions, no p values.
#
# Gustavo has asked for an analysis of mutational burden against RASSF9 and
# RASSF10. Before any of that is written, this reports what burden data we
# actually hold and under what column names, because guessing a column name has
# already cost this project two wrong joins today.
#
# It answers:
#   1. which files carry mutational-signature or burden data
#   2. every column in them, with type, range, and missingness
#   3. which columns plausibly represent total burden (TMB), UV burden (SBS7),
#      and the UV fraction
#   4. whether purity (CPE) is present, since purity is the central confound
#   5. how many samples would survive a join to expression
#
# WHY PURITY MATTERS HERE, stated before any numbers exist. Three platforms now
# agree that RASSF9 in bulk tumour tissue is largely NOT melanocyte-derived, and
# TCGA purity correlates with RASSF9 expression at rho = -0.204 (p = 3.05e-05).
# Mutational burden estimates are themselves purity-dependent, because variant
# allele frequencies fall as normal tissue contaminates the sample. A negative
# burden-expression correlation is therefore exactly what composition alone
# would produce, with no repression in the mechanism. Any association found in
# 39b must be tested against that, the same way the methylation association was.
#
# Output: results/logs/39a_burden_inspect.log
# ==============================================================================

source(here::here("R", "00_setup.R"))

LOG <- file.path(PATHS$logs, "39a_burden_inspect.log")
cat("Mutational burden data inventory —", format(Sys.time()), "\n\n", file = LOG)
hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

PROC <- PATHS$processed

hdr("1. Where burden data might live")

# Candidate files, in the order they are most likely to carry burden columns.
cands <- c(file.path(PROC, "full_cohort.csv"),
           file.path(PROC, "analysis_cohort.csv"))
sig <- list.files(PATHS$raw, pattern = "sigProfiler|SBS|signature", recursive = TRUE,
                  full.names = TRUE, ignore.case = TRUE)
if (length(sig)) {
  log_line("  signature files found under data/raw:", logfile = LOG)
  for (f in sig)
    log_line(sprintf("    %-62s %8.2f MB", basename(f), file.size(f) / 1024^2),
             logfile = LOG)
  cands <- c(cands, sig)
} else {
  log_line("  [FLAG] no file matching sigProfiler/SBS/signature under data/raw.",
           logfile = LOG)
  log_line("  If the SigProfiler table lives elsewhere (it was originally in a",
           logfile = LOG)
  log_line("  Downloads folder), 39b will need its path.", logfile = LOG)
}

cands <- cands[file.exists(cands)]
if (!length(cands)) stop("None of the candidate files exist.")

hdr("2. Columns in each file")

BURDEN_PAT <- "tmb|burden|mutation|sbs|signature|snv|nmut|count"
PURITY_PAT <- "cpe|purity|absolute|estimate|tumou?r_?frac"
ID_PAT     <- "patient|barcode|sample|case|id$"

peek_cols <- function(f) {
  d <- try(suppressWarnings(readr::read_csv(f, n_max = 5000,
            show_col_types = FALSE, progress = FALSE)), silent = TRUE)
  if (inherits(d, "try-error")) {
    d <- try(suppressWarnings(data.table::fread(f, nrows = 5000)), silent = TRUE)
    if (inherits(d, "try-error")) return(NULL)
    d <- tibble::as_tibble(d)
  }
  d
}

for (f in cands) {
  d <- peek_cols(f)
  log_line("\n  ", basename(f), logfile = LOG)
  if (is.null(d)) { log_line("    could not be parsed", logfile = LOG); next }
  log_line("    ", nrow(d), " rows (first 5,000) x ", ncol(d), " columns",
           logfile = LOG)

  hits_b <- grep(BURDEN_PAT, names(d), value = TRUE, ignore.case = TRUE)
  hits_p <- grep(PURITY_PAT, names(d), value = TRUE, ignore.case = TRUE)
  hits_i <- grep(ID_PAT,     names(d), value = TRUE, ignore.case = TRUE)
  log_line("    identifier-like columns: ",
           if (length(hits_i)) paste(utils::head(hits_i, 8), collapse = ", ") else "none",
           logfile = LOG)
  log_line("    purity-like columns:     ",
           if (length(hits_p)) paste(hits_p, collapse = ", ") else "none",
           logfile = LOG)
  log_line("    burden-like columns (", length(hits_b), "): ",
           if (length(hits_b)) paste(utils::head(hits_b, 25), collapse = ", ") else "none",
           logfile = LOG)
  if (length(hits_b) > 25) log_line("      ... and ", length(hits_b) - 25, " more",
                                    logfile = LOG)

  # Full column list only when the file is narrow enough to be readable.
  if (ncol(d) <= 60) {
    log_line("    all columns: ", paste(names(d), collapse = ", "), logfile = LOG)
  }

  show <- unique(c(hits_b, hits_p))
  show <- show[seq_len(min(length(show), 18))]
  if (length(show)) {
    log_line(sprintf("    %-26s %-10s %10s %10s %10s %8s", "column", "type",
                     "min", "median", "max", "NA"), logfile = LOG)
    for (cn in show) {
      v <- d[[cn]]
      if (is.numeric(v)) {
        log_line(sprintf("    %-26s %-10s %10.3g %10.3g %10.3g %8d", cn,
                         "numeric", suppressWarnings(min(v, na.rm = TRUE)),
                         suppressWarnings(stats::median(v, na.rm = TRUE)),
                         suppressWarnings(max(v, na.rm = TRUE)), sum(is.na(v))),
                 logfile = LOG)
      } else {
        u <- unique(as.character(v))
        log_line(sprintf("    %-26s %-10s  %d distinct, e.g. %s", cn,
                         class(v)[1], length(u),
                         paste(utils::head(u, 3), collapse = " / ")),
                 logfile = LOG)
      }
    }
  }
}

hdr("3. What 39b needs")

log_line("  1. the file and column holding TOTAL burden (TMB), and whether it",
         logfile = LOG)
log_line("     is a raw count or already per-megabase", logfile = LOG)
log_line("  2. the columns holding SBS7a/b/c/d, so UV burden and the UV",
         logfile = LOG)
log_line("     FRACTION can be derived rather than assumed", logfile = LOG)
log_line("  3. whether a NON-UV burden (TMB minus SBS7) can be computed, which",
         logfile = LOG)
log_line("     is the cleanest way to separate a UV-specific effect from a",
         logfile = LOG)
log_line("     general mutation-load effect", logfile = LOG)
log_line("  4. the purity column, without which the central confound cannot be",
         logfile = LOG)
log_line("     tested", logfile = LOG)
log_line("  5. the identifier to join on, and how many samples survive",
         logfile = LOG)
log_line("", logfile = LOG)
log_line("  PRE-SPECIFICATION, recorded now, before any association is computed:",
         logfile = LOG)
log_line("    PRIMARY   Spearman correlation of RASSF9 expression with TOTAL",
         logfile = LOG)
log_line("              mutational burden, in the largest cohort with both.",
         logfile = LOG)
log_line("    FAMILY    total burden, UV (SBS7) burden, non-UV burden, UV",
         logfile = LOG)
log_line("              fraction -- four definitions x two genes = 8 tests,",
         logfile = LOG)
log_line("              Benjamini-Hochberg corrected across the family. Every",
         logfile = LOG)
log_line("              one is reported whatever it shows.", logfile = LOG)
log_line("    CONFOUND  the association must survive adjustment for purity and",
         logfile = LOG)
log_line("              sample type, and the premises of the purity explanation",
         logfile = LOG)
log_line("              are tested individually as they were for methylation.",
         logfile = LOG)
log_line("    ROBUSTNESS Spearman alongside Pearson on log burden, since the",
         logfile = LOG)
log_line("              original analysis found rho = -0.130 (p = 0.021) while",
         logfile = LOG)
log_line("              r = -0.033 (p = 0.564) -- a gap that usually means a",
         logfile = LOG)
log_line("              few extreme-burden tumours drive the rank result.",
         logfile = LOG)
log_line("    NOT CLAIMABLE Cross-sectional bulk data cannot establish that",
         logfile = LOG)
log_line("              burden causes repression, nor the direction of any",
         logfile = LOG)
log_line("              arrow. 39b will report association only.", logfile = LOG)

report_warnings(logfile = LOG)
log_line("\nFull log: ", LOG, logfile = NULL)
