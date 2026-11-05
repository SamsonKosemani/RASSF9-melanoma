# ==============================================================================
# 12_geo_inventory.R
#
# Purpose: find out what the four requested GEO datasets actually contain before
# any analysis is written. Inventory only -- no hypothesis is tested here.
#
# For each series this reports: title, platform, number of samples, the sample
# grouping fields, whether a processed expression matrix is available through
# GEOquery, what supplementary files exist, and whether RASSF9 and RASSF10 are
# represented on the platform at all.
#
# That last point is not a formality. GSE4587 is on Affymetrix HG-U133 Plus 2.0,
# designed in 2003; RASSF10 was characterised later and may have no probe set.
# If a gene is not on the array, no amount of analysis will produce an answer,
# and that is a result worth knowing before writing code around it.
#
# Output: data/raw/geo/              cached GEO downloads
#         results/logs/12_geo_inventory.log
#         results/geo_inventory.csv  one row per dataset
# ==============================================================================

source(here::here("R", "00_setup.R"))

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))

LOG <- file.path(PATHS$logs, "12_geo_inventory.log")
cat("GEO inventory —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

GEO_DIR <- file.path(PATHS$raw, "geo")
dir.create(GEO_DIR, recursive = TRUE, showWarnings = FALSE)

GENES <- c("RASSF9", "RASSF10")

SERIES <- c(
  "GSE112509",   # nevi (23) vs primary melanoma (57), RNA-seq
  "GSE4587",     # melanoma progression series, Affymetrix U133 Plus 2.0
  "GSE109245",   # melanocyte vs metastatic melanoma lines (unverified)
  "GSE230574"    # scRNA-seq, melanocytes / nevus cells / melanoma lines
)

options(timeout = 1200)   # GEO can be slow

inv <- list()

for (gse in SERIES) {
  hdr(gse)

  obj <- try(getGEO(gse, destdir = GEO_DIR, GSEMatrix = TRUE,
                    getGPL = TRUE), silent = TRUE)

  if (inherits(obj, "try-error")) {
    log_line("  [FAILED] could not retrieve ", gse, logfile = LOG)
    log_line("  message: ", as.character(obj), logfile = LOG)
    log_line("  -> check the accession is correct, or whether the series has",
             logfile = LOG)
    log_line("     no series matrix (common for RNA-seq, where only",
             logfile = LOG)
    log_line("     supplementary count files are deposited).", logfile = LOG)

    # Even without a matrix, the supplementary file listing is informative.
    sup <- try(getGEOSuppFiles(gse, baseDir = GEO_DIR, fetch_files = FALSE),
               silent = TRUE)
    if (!inherits(sup, "try-error") && !is.null(sup) && nrow(sup)) {
      log_line("  supplementary files available:", logfile = LOG)
      for (i in seq_len(nrow(sup))) {
        log_line("    ", rownames(sup)[i], logfile = LOG)
      }
    } else {
      log_line("  no supplementary file listing returned either", logfile = LOG)
    }

    inv[[gse]] <- tibble::tibble(
      series = gse, status = "failed", platform = NA_character_,
      n_samples = NA_integer_, expr_matrix = FALSE,
      rassf9 = NA, rassf10 = NA, note = "retrieval failed"
    )
    next
  }

  # getGEO returns a list, one ExpressionSet per platform.
  for (j in seq_along(obj)) {
    es <- obj[[j]]
    pf <- annotation(es)
    n  <- ncol(es)

    log_line("  title:    ", as.character(experimentData(es)@title), logfile = LOG)
    log_line("  platform: ", pf, logfile = LOG)
    log_line("  samples:  ", n, logfile = LOG)
    log_line("  features: ", nrow(es), logfile = LOG)

    pd <- pData(es)

    # --- which column describes the groups? -----------------------------------
    cand <- names(pd)[grepl("title|characteristic|source|group|tissue|cell|type|stage|diagnosis",
                            names(pd), ignore.case = TRUE)]
    cand <- cand[!grepl("supplementary|relation|contact|protocol|data_processing|extract",
                        cand, ignore.case = TRUE)]

    log_line("\n  candidate grouping columns:", logfile = LOG)
    for (cc in head(cand, 8)) {
      vals <- unique(as.character(pd[[cc]]))
      if (length(vals) > 1 && length(vals) <= n) {
        log_line(sprintf("    %-32s %d distinct value(s)", cc, length(vals)),
                 logfile = LOG)
        for (v in head(vals, 8)) log_line("        ", substr(v, 1, 90),
                                          logfile = LOG)
        if (length(vals) > 8) log_line("        ... and ", length(vals) - 8,
                                       " more", logfile = LOG)
      }
    }

    # --- sample titles, which usually encode the group ------------------------
    if ("title" %in% names(pd)) {
      log_line("\n  sample titles (first 20):", logfile = LOG)
      for (t in head(as.character(pd$title), 20)) {
        log_line("    ", t, logfile = LOG)
      }
      if (n > 20) log_line("    ... and ", n - 20, " more", logfile = LOG)
    }

    # --- is there usable expression data? -------------------------------------
    ex <- try(exprs(es), silent = TRUE)
    has_expr <- !inherits(ex, "try-error") && nrow(ex) > 0 && ncol(ex) > 0
    log_line("\n  processed expression matrix via GEOquery: ",
             ifelse(has_expr, "YES", "NO (supplementary files needed)"),
             logfile = LOG)
    if (has_expr) {
      log_line(sprintf("    range %.2f to %.2f, %.1f%% missing",
                       min(ex, na.rm = TRUE), max(ex, na.rm = TRUE),
                       100 * mean(is.na(ex))), logfile = LOG)
    }

    # --- are our genes on this platform? --------------------------------------
    fd <- fData(es)
    found <- setNames(rep(FALSE, length(GENES)), GENES)
    if (ncol(fd) > 0) {
      sym_cols <- names(fd)[grepl("symbol|gene|assign", names(fd),
                                  ignore.case = TRUE)]
      log_line("\n  feature annotation columns: ",
               paste(head(names(fd), 12), collapse = ", "), logfile = LOG)

      for (g in GENES) {
        hits <- character(0)
        for (sc in sym_cols) {
          v <- as.character(fd[[sc]])
          # whole-token match: RASSF9 must not match RASSF9P1
          idx <- which(vapply(strsplit(v, "[ ///;,|]+"),
                              function(p) g %in% trimws(p), logical(1)))
          if (length(idx)) hits <- union(hits, rownames(fd)[idx])
        }
        found[g] <- length(hits) > 0
        if (length(hits)) {
          log_line(sprintf("    [ ok ] %-8s %d probe(s)/feature(s): %s", g,
                           length(hits),
                           paste(head(hits, 6), collapse = ", ")), logfile = LOG)
        } else {
          log_line("    [FLAG] ", g, " NOT found in the platform annotation",
                   logfile = LOG)
        }
      }
    } else {
      log_line("\n  no feature annotation returned; gene presence must be",
               logfile = LOG)
      log_line("  checked in the supplementary files instead.", logfile = LOG)
    }

    # --- supplementary files --------------------------------------------------
    sup <- try(getGEOSuppFiles(gse, baseDir = GEO_DIR, fetch_files = FALSE),
               silent = TRUE)
    if (!inherits(sup, "try-error") && !is.null(sup) && nrow(sup)) {
      log_line("\n  supplementary files:", logfile = LOG)
      for (i in seq_len(nrow(sup))) log_line("    ", rownames(sup)[i],
                                             logfile = LOG)
    }

    inv[[paste0(gse, "_", j)]] <- tibble::tibble(
      series = gse, status = "ok", platform = pf, n_samples = n,
      expr_matrix = has_expr,
      rassf9 = unname(found["RASSF9"]), rassf10 = unname(found["RASSF10"]),
      note = ""
    )
  }
}

# ==============================================================================
# Summary
# ==============================================================================
hdr("SUMMARY")

inv_tbl <- dplyr::bind_rows(inv)
readr::write_csv(inv_tbl, file.path(PATHS$results, "geo_inventory.csv"))

cat(capture.output(print(as.data.frame(inv_tbl), row.names = FALSE)),
    sep = "\n", file = LOG, append = TRUE)
print(as.data.frame(inv_tbl), row.names = FALSE)

log_line("\n  Read the per-dataset sections above before the analysis is",
         logfile = LOG)
log_line("  written. Three things decide what is feasible:", logfile = LOG)
log_line("    1. whether both genes are on the platform at all", logfile = LOG)
log_line("    2. how many samples per group (GSE4587 and GSE230574 are small)",
         logfile = LOG)
log_line("    3. whether expression comes from GEOquery or supplementary files",
         logfile = LOG)

report_warnings(LOG)
write_session_info("12_geo_inventory")

log_line("\nFull log: ", LOG, logfile = LOG)
