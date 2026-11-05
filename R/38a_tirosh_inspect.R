# ==============================================================================
# 38a_tirosh_inspect.R
#
# INSPECT ONLY. No analysis, no conclusions.
#
# Target: Tirosh et al., Science 2016 (GSE72056) -- Smart-seq2 melanoma, ~4,645
# cells from ~19 patients. The reason to use this one rather than any other
# melanoma single-cell dataset is that it shares Belote's chemistry, so a
# normal-melanocyte versus malignant-melanocyte contrast across the two is the
# one cross-dataset comparison where platform difference does not sink it.
#
# TWO THINGS ARE EXPECTED TO BE AWKWARD, and this script is written to find out
# rather than assume:
#
# 1. THE VALUE SCALE. This matrix is, by reputation, distributed as
#    log2(TPM/10 + 1) rather than raw counts. If so, binomial thinning is
#    IMPOSSIBLE -- you cannot downsample counts that are not counts -- and the
#    depth control that made the Belote result trustworthy is unavailable in
#    that form. The substitute would be to match cells on detected-gene count,
#    the best remaining depth proxy, and compare only within a common range.
#    That is weaker and would have to be labelled as such in the manuscript.
#    So the scale decides the whole design and is read off the data here.
#
# 2. EMBEDDED ANNOTATION ROWS. The first few rows are believed to carry the
#    tumour of origin, a malignant/non-malignant flag and a cell-type code,
#    sitting above the gene rows in the same file. Treated as genes they would
#    silently corrupt everything. This script detects non-gene-looking row
#    labels and prints their distinct values, so the encoding is learned from
#    the file rather than recalled.
#
# It scans the first field of EVERY row, not a capped sample -- the lesson from
# 37a, where a 4,000-row cap would have reported RASSF9 absent purely for
# sitting further down a 26,576-gene matrix.
#
# USAGE: download the series supplementary file(s) into DATA_DIR, then run.
#        Do not unzip; .gz is read directly.
#
# Output: results/logs/38a_tirosh_inspect.log
# ==============================================================================

source(here::here("R", "00_setup.R"))

LOG <- file.path(PATHS$logs, "38a_tirosh_inspect.log")
cat("Tirosh GSE72056 inspection —", format(Sys.time()), "\n\n", file = LOG)
hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

# ---- EDIT IF YOU PUT THE FILES SOMEWHERE ELSE --------------------------------
DATA_DIR <- file.path(PATHS$raw, "geo", "tirosh")
# ------------------------------------------------------------------------------

TARGETS <- c("RASSF9", "RASSF10")
MARKERS <- c("MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10", "TYRP1",
             "PTPRC", "COL1A1", "PECAM1", "CD3D", "KRT14")

hdr("1. Files")

if (!dir.exists(DATA_DIR)) {
  log_line("  DATA_DIR does not exist: ", DATA_DIR, logfile = LOG)
  log_line("  Create it, download the series supplementary file(s) there,",
           logfile = LOG)
  log_line("  and rerun. Do not unzip.", logfile = LOG)
  log_line("\nFull log: ", LOG, logfile = NULL)
  quit(save = "no", status = 0)
}
ff <- list.files(DATA_DIR, recursive = TRUE, full.names = TRUE)
if (!length(ff)) { log_line("  directory is empty", logfile = LOG)
  quit(save = "no", status = 0) }
for (f in ff)
  log_line(sprintf("  %-62s %10.2f MB", basename(f), file.size(f) / 1024^2),
           logfile = LOG)

open_any <- function(f) if (grepl("\\.gz$", f, ignore.case = TRUE))
  gzfile(f, "rt") else file(f, "rt")
peek <- function(f, n = 4) {
  con <- open_any(f); on.exit(close(con), add = TRUE)
  readLines(con, n = n, warn = FALSE)
}
guess_delim <- function(line) {
  cnt <- c("\t" = lengths(regmatches(line, gregexpr("\t", line))),
           "," = lengths(regmatches(line, gregexpr(",", line))),
           ";" = lengths(regmatches(line, gregexpr(";", line))))
  names(cnt)[which.max(cnt)]
}

hdr("2. Head of each text file")

TXT <- ff[grepl("\\.(csv|tsv|txt|tab)(\\.gz)?$", ff, ignore.case = TRUE)]
if (!length(TXT)) {
  log_line("  no delimited text files found. If the download is a .tar,",
           logfile = LOG)
  log_line("  extract it into its own folder and rerun.", logfile = LOG)
}
info <- list()
for (f in TXT) {
  h <- try(peek(f, 4), silent = TRUE)
  if (inherits(h, "try-error") || !length(h)) next
  d <- guess_delim(h[1])
  nf <- length(strsplit(h[1], d, fixed = TRUE)[[1]])
  log_line("\n  ", basename(f), logfile = LOG)
  log_line("    delimiter: ", if (d == "\t") "TAB" else d,
           " | header fields: ", nf, logfile = LOG)
  for (i in seq_along(h))
    log_line("    [", i, "] ", substr(h[i], 1, 200),
             if (nchar(h[i]) > 200) " ..." else "", logfile = LOG)
  info[[f]] <- list(delim = d, nfields = nf)
}

hdr("3. Row labels: genes, and anything that is not a gene")

log_line("  Every row's first field is read. Rows whose label does not look",
         logfile = LOG)
log_line("  like a gene symbol are candidate ANNOTATION rows -- their distinct",
         logfile = LOG)
log_line("  values are printed so the encoding is learned from the file.",
         logfile = LOG)

gene_like <- function(x) grepl("^(ENSG[0-9]+|[A-Za-z][A-Za-z0-9\\.\\-_]{0,19})$", x)

for (f in names(info)) {
  d <- info[[f]]$delim
  lab <- try({
    con <- open_any(f); on.exit(try(close(con), silent = TRUE), add = TRUE)
    acc <- character(0)
    repeat {
      L <- readLines(con, n = 20000, warn = FALSE)
      if (!length(L)) break
      acc <- c(acc, vapply(strsplit(L, d, fixed = TRUE), function(p)
        if (length(p)) gsub('^"|"$', "", p[1]) else NA_character_,
        character(1)))
      if (length(acc) > 2e6) break
    }
    acc
  }, silent = TRUE)
  if (inherits(lab, "try-error")) next

  body <- lab[-1]
  log_line("\n  ", basename(f), ": ", length(body), " data rows (entire file)",
           logfile = LOG)
  log_line("    first 12 row labels: ",
           paste(utils::head(body, 12), collapse = " | "), logfile = LOG)

  # Annotation rows are usually at the top and have long or spaced labels.
  odd <- which(!gene_like(body) | nchar(body) > 20)
  odd_top <- odd[odd <= 20]
  if (length(odd_top)) {
    log_line("    CANDIDATE ANNOTATION ROWS in the first 20 rows: ",
             paste(odd_top, collapse = ", "), logfile = LOG)
    con <- open_any(f); on.exit(try(close(con), silent = TRUE), add = TRUE)
    L <- readLines(con, n = max(odd_top) + 1, warn = FALSE); close(con)
    for (r in odd_top) {
      p <- strsplit(L[r + 1], d, fixed = TRUE)[[1]]
      vals <- gsub('^"|"$', "", p[-1])
      u <- sort(unique(vals))
      log_line("      row ", r, ' label "', body[r], '"', logfile = LOG)
      log_line("        distinct values (", length(u), "): ",
               paste(utils::head(u, 30), collapse = ", "),
               if (length(u) > 30) " ..." else "", logfile = LOG)
      if (length(u) <= 30) {
        tb <- sort(table(vals), decreasing = TRUE)
        log_line("        counts: ", paste(names(tb), as.integer(tb), sep = "=",
                                           collapse = ", "), logfile = LOG)
      }
    }
    log_line("    -> 38b must SKIP these rows. Treated as genes they would",
             logfile = LOG)
    log_line("       corrupt every downstream number silently.", logfile = LOG)
  } else {
    log_line("    no annotation rows detected in the first 20", logfile = LOG)
  }

  for (g in TARGETS) {
    hit <- grep(paste0("^", g, "$"), body, value = TRUE)
    log_line("    ", g, ": ", if (length(hit)) paste(hit, collapse = ", ") else
             "NOT FOUND as an exact symbol", logfile = LOG)
    if (!length(hit)) {
      loose <- grep(g, body, value = TRUE, ignore.case = TRUE)
      if (length(loose))
        log_line("      loose matches: ", paste(utils::head(loose, 8),
                 collapse = ", "), logfile = LOG)
    }
  }
  present <- MARKERS[MARKERS %in% body]
  log_line("    markers present: ",
           if (length(present)) paste(present, collapse = ", ") else "NONE",
           logfile = LOG)
  dup <- sum(duplicated(body[gene_like(body)]))
  if (dup) log_line("    [FLAG] ", dup, " duplicated gene labels -- 38b must",
                    " decide how to collapse them", logfile = LOG)
}

hdr("4. Value scale -- this decides whether thinning is possible")

for (f in names(info)) {
  d <- info[[f]]$delim
  v <- try({
    con <- open_any(f); on.exit(try(close(con), silent = TRUE), add = TRUE)
    L <- readLines(con, n = 400, warn = FALSE)
    if (length(L) < 25) NULL else {
      p <- strsplit(L[25:length(L)], d, fixed = TRUE)   # skip annotation rows
      suppressWarnings(as.numeric(unlist(lapply(p, function(x)
        utils::head(x[-1], 60)))))
    }
  }, silent = TRUE)
  if (inherits(v, "try-error") || is.null(v)) next
  v <- v[is.finite(v)]
  if (!length(v)) next
  int_frac <- mean(abs(v - round(v)) < 1e-9)
  log_line("\n  ", basename(f), logfile = LOG)
  log_line(sprintf("    min %.4f | median %.4f | max %.3f | zeros %.1f%% | integer-valued %.0f%%",
                   min(v), stats::median(v), max(v), 100 * mean(v == 0),
                   100 * int_frac), logfile = LOG)
  verdict <- if (int_frac > 0.98 && max(v) > 100) {
    "RAW COUNTS -> binomial thinning possible, same design as 37b"
  } else if (max(v) < 25 && int_frac < 0.95) {
    "LOG-TRANSFORMED (likely log2(TPM/10+1)) -> THINNING IMPOSSIBLE"
  } else if (max(v) > 1000 && int_frac < 0.95) {
    "LINEAR TPM/CPM -> thinning impossible; counts not recoverable"
  } else "UNCLEAR -- inspect manually before designing 38b"
  log_line("    scale: ", verdict, logfile = LOG)
  if (grepl("IMPOSSIBLE|not recoverable", verdict)) {
    log_line("    CONSEQUENCE: depth cannot be equalised by downsampling. The",
             logfile = LOG)
    log_line("    substitute is to match cells on DETECTED-GENE COUNT between",
             logfile = LOG)
    log_line("    datasets and compare only within a common range. That is a", logfile = LOG)
    log_line("    weaker control than thinning and must be stated as such.",
             logfile = LOG)
    log_line("    The raw counts exist in SRA if the comparison needs to be",
             logfile = LOG)
    log_line("    airtight, at considerably more work.", logfile = LOG)
  }
}

hdr("5. What 38b needs from this")

log_line("  1. which file holds the matrix, its delimiter, and how many",
         logfile = LOG)
log_line("     annotation rows to skip", logfile = LOG)
log_line("  2. the encoding of the malignant flag and the cell-type code,",
         logfile = LOG)
log_line("     taken from the printed distinct values rather than recalled",
         logfile = LOG)
log_line("  3. the exact spelling of RASSF9 and RASSF10, or confirmation that",
         logfile = LOG)
log_line("     they are absent -- which would be a fact about this matrix's",
         logfile = LOG)
log_line("     gene filtering, NOT about melanoma biology", logfile = LOG)
log_line("  4. the value scale, which decides whether depth can be equalised",
         logfile = LOG)
log_line("     at all and therefore how strong the comparison can be", logfile = LOG)
log_line("  5. whether gene labels are duplicated, and the patient identifier",
         logfile = LOG)
log_line("     so donors remain the unit of analysis", logfile = LOG)
log_line("", logfile = LOG)
log_line("  REMINDER OF THE COMPARISON'S LIMIT: Belote melanocytes and Tirosh",
         logfile = LOG)
log_line("  malignant cells come from different studies, different labs and",
         logfile = LOG)
log_line("  different processing. Shared chemistry removes the largest source",
         logfile = LOG)
log_line("  of difference but not batch. Any contrast is suggestive, not a",
         logfile = LOG)
log_line("  controlled experiment, and should be reported that way.",
         logfile = LOG)

report_warnings(logfile = LOG)
log_line("\nFull log: ", LOG, logfile = NULL)
