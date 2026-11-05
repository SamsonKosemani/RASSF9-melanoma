# ==============================================================================
# 37a_belote_inspect.R
#
# INSPECT ONLY. No analysis, no conclusions.
#
# This reports what the downloaded Belote et al. (Nat Cell Biol 2021) Smart-seq2
# files actually contain, so that 37b can be written against their real shape
# instead of a guess. The Slieker metadata join matched 0 of 39,481 spots on the
# first attempt because the file's structure was assumed rather than read; this
# script exists so that does not happen twice.
#
# It answers, and nothing more:
#   - which files are present, their size and delimiter
#   - which dimension holds genes and which holds cells
#   - whether RASSF9 and RASSF10 are present as features AT ALL, and under what
#     exact spelling (symbol, Ensembl ID, versioned ID)
#   - what the value scale looks like (raw integer counts, TPM, log-transformed)
#   - which metadata columns exist, and which one carries the cell-type label
#   - how many cells carry a melanocyte label, by whatever vocabulary the
#     authors used
#
# Run it, send the output, and 37b follows.
#
# USAGE: edit DATA_DIR below to wherever the files were unzipped.
#
# Output: results/logs/37a_belote_inspect.log
# ==============================================================================

source(here::here("R", "00_setup.R"))

LOG <- file.path(PATHS$logs, "37a_belote_inspect.log")
cat("Belote Smart-seq2 inspection —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

# ---- EDIT THIS ---------------------------------------------------------------
DATA_DIR <- file.path(PATHS$raw, "geo", "belote")
# ------------------------------------------------------------------------------

TARGETS <- c("RASSF9", "RASSF10")
MEL_WORDS <- c("melanocyt", "melan", "mc", "pigment")
MARKERS <- c("MLANA", "PMEL", "TYR", "DCT", "MITF", "SOX10", "TYRP1",
             "KRT14", "KRT1", "PTPRC", "COL1A1")

hdr("1. Files")

if (!dir.exists(DATA_DIR)) {
  log_line("  DATA_DIR does not exist: ", DATA_DIR, logfile = LOG)
  log_line("  Create it, put the downloaded files there, and edit DATA_DIR.",
           logfile = LOG)
  log_line("\nFull log: ", LOG, logfile = NULL)
  quit(save = "no", status = 0)
}

ff <- list.files(DATA_DIR, recursive = TRUE, full.names = TRUE)
if (!length(ff)) {
  log_line("  directory is empty: ", DATA_DIR, logfile = LOG)
  quit(save = "no", status = 0)
}
for (f in ff) {
  log_line(sprintf("  %-58s %10.2f MB", basename(f), file.size(f) / 1024^2),
           logfile = LOG)
}

# A gz file is read directly; nothing is unpacked, so the download stays pristine.
peek <- function(f, n = 3) {
  con <- if (grepl("\\.gz$", f, ignore.case = TRUE)) gzfile(f, "rt") else file(f, "rt")
  on.exit(close(con), add = TRUE)
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
info <- list()
for (f in TXT) {
  h <- try(peek(f, 3), silent = TRUE)
  if (inherits(h, "try-error") || !length(h)) {
    log_line("\n  ", basename(f), ": could not read as text", logfile = LOG)
    next
  }
  d <- guess_delim(h[1])
  nf <- length(strsplit(h[1], d, fixed = TRUE)[[1]])
  log_line("\n  ", basename(f), logfile = LOG)
  log_line("    delimiter: ", if (d == "\t") "TAB" else d, " | header fields: ",
           nf, logfile = LOG)
  for (i in seq_along(h)) {
    s <- substr(h[i], 1, 220)
    log_line("    [", i, "] ", s, if (nchar(h[i]) > 220) " ..." else "",
             logfile = LOG)
  }
  info[[f]] <- list(delim = d, nfields = nf)
}

hdr("3. Which dimension holds genes?")

log_line("  Strategy: read the first column of each file and test whether its",
         logfile = LOG)
log_line("  values look like gene identifiers. Whichever axis carries the",
         logfile = LOG)
log_line("  markers is the gene axis -- this is not assumed.", logfile = LOG)

looks_gene <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(0)
  mean(grepl("^(ENSG[0-9]+|[A-Z][A-Z0-9orf\\.\\-]{1,14})$", x))
}

for (f in names(info)) {
  d <- info[[f]]$delim
  # Read the first field of EVERY line, streamed in chunks. An earlier version
  # capped this at 4,000 rows, which would have reported RASSF9 absent purely
  # because it sits past row 4,000 of a ~25,000-gene matrix -- a false negative
  # on the one question this whole script exists to answer.
  first <- try({
    con <- if (grepl("\\.gz$", f, ignore.case = TRUE)) gzfile(f, "rt") else file(f, "rt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    acc <- character(0)
    repeat {
      L <- readLines(con, n = 20000, warn = FALSE)
      if (!length(L)) break
      acc <- c(acc, vapply(strsplit(L, d, fixed = TRUE), function(p)
        if (length(p)) gsub('^"|"$', "", p[1]) else NA_character_, character(1)))
      if (length(acc) > 2e6) break      # runaway guard only
    }
    acc
  }, silent = TRUE)
  if (inherits(first, "try-error")) next

  hdr_fields <- gsub('^"|"$', "", strsplit(peek(f, 1)[1], d, fixed = TRUE)[[1]])
  r1 <- looks_gene(first[-1])
  r2 <- looks_gene(hdr_fields[-1])
  log_line("\n  ", basename(f), logfile = LOG)
  log_line(sprintf("    rows look like genes: %.0f%% | columns look like genes: %.0f%%",
                   100 * r1, 100 * r2), logfile = LOG)
  log_line("    rows read: ", length(first) - 1, " (entire file)",
           logfile = LOG)

  # Exact spelling matters: a versioned Ensembl ID will not match a plain one.
  for (g in TARGETS) {
    hits_r <- grep(paste0("^", g, "$|^", g, "\\."), first, value = TRUE)
    hits_c <- grep(paste0("^", g, "$|^", g, "\\."), hdr_fields, value = TRUE)
    loose_r <- grep(g, first, value = TRUE, ignore.case = TRUE)
    log_line("    ", g, ": exact in rows = ",
             if (length(hits_r)) paste(hits_r, collapse = ", ") else "no",
             " | exact in cols = ",
             if (length(hits_c)) paste(hits_c, collapse = ", ") else "no",
             logfile = LOG)
    if (!length(hits_r) && length(loose_r)) {
      log_line("      loose row matches (could be paralogues): ",
               paste(utils::head(loose_r, 6), collapse = ", "), logfile = LOG)
    }
  }
  present <- MARKERS[MARKERS %in% first]
  log_line("    markers found in rows: ",
           if (length(present)) paste(present, collapse = ", ") else "none",
           logfile = LOG)
  log_line("    first 8 row names: ",
           paste(utils::head(first[-1], 8), collapse = ", "), logfile = LOG)
  log_line("    first 6 column names: ",
           paste(utils::head(hdr_fields[-1], 6), collapse = ", "), logfile = LOG)
}

hdr("4. Value scale")

log_line("  Raw counts, TPM and log-transformed values need different handling,",
         logfile = LOG)
log_line("  so the scale is read off the data rather than assumed from the",
         logfile = LOG)
log_line("  filename.", logfile = LOG)

for (f in names(info)) {
  d <- info[[f]]$delim
  v <- try({
    con <- if (grepl("\\.gz$", f, ignore.case = TRUE)) gzfile(f, "rt") else file(f, "rt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    L <- readLines(con, n = 200, warn = FALSE)
    if (length(L) < 2) NULL else {
      p <- strsplit(L[-1], d, fixed = TRUE)
      suppressWarnings(as.numeric(unlist(lapply(p, function(x)
        utils::head(x[-1], 40)))))
    }
  }, silent = TRUE)
  if (inherits(v, "try-error") || is.null(v)) next
  v <- v[is.finite(v)]
  if (!length(v)) next
  int_frac <- mean(abs(v - round(v)) < 1e-9)
  log_line("\n  ", basename(f), logfile = LOG)
  log_line(sprintf("    min %.4f | median %.4f | max %.2f | zeros %.1f%% | integer-valued %.0f%%",
                   min(v), stats::median(v), max(v), 100 * mean(v == 0),
                   100 * int_frac), logfile = LOG)
  guess <- if (int_frac > 0.98 && max(v) > 100) "raw counts" else
    if (max(v) <= 25 && int_frac < 0.9) "log-transformed" else
    if (max(v) > 1000) "TPM/CPM (linear)" else "unclear"
  log_line("    scale looks like: ", guess, logfile = LOG)
}

hdr("5. Metadata and cell-type labels")

for (f in names(info)) {
  if (info[[f]]$nfields > 80) next     # wide = matrix, not metadata
  M <- try(suppressWarnings(readr::read_delim(
    f, delim = info[[f]]$delim, show_col_types = FALSE, progress = FALSE,
    n_max = 200000)), silent = TRUE)
  if (inherits(M, "try-error") || !ncol(M)) next
  log_line("\n  ", basename(f), ": ", nrow(M), " rows x ", ncol(M), " cols",
           logfile = LOG)
  log_line("    columns: ", paste(names(M), collapse = ", "), logfile = LOG)
  for (cn in names(M)) {
    x <- M[[cn]]
    if (!is.character(x) && !is.factor(x)) next
    u <- unique(as.character(x)); u <- u[!is.na(u)]
    if (length(u) > 40 || !length(u)) next
    looks_ct <- any(vapply(MEL_WORDS, function(w)
      any(grepl(w, u, ignore.case = TRUE)), logical(1)))
    log_line("    ", cn, if (looks_ct) "  <-- CANDIDATE cell-type column" else "",
             logfile = LOG)
    tb <- sort(table(as.character(x)), decreasing = TRUE)
    log_line("      ", paste(names(tb), as.integer(tb), sep = "=",
                             collapse = ", "), logfile = LOG)
  }
}

hdr("6. What 37b needs from this")

log_line("  1. the matrix file name, its delimiter, and which axis is genes",
         logfile = LOG)
log_line("  2. the EXACT feature spelling for RASSF9 (symbol or Ensembl,",
         logfile = LOG)
log_line("     versioned or not) -- or confirmation that it is absent, which",
         logfile = LOG)
log_line("     would itself be a finding about the reference annotation rather",
         logfile = LOG)
log_line("     than about melanocyte biology", logfile = LOG)
log_line("  3. the value scale, so counts are not treated as TPM", logfile = LOG)
log_line("  4. the cell-type column and the exact melanocyte label string",
         logfile = LOG)
log_line("  5. whether donor, age/stage and anatomical site are present, since",
         logfile = LOG)
log_line("     those decide whether any stratified test is possible at all",
         logfile = LOG)
log_line("", logfile = LOG)
log_line("  NOTE: if RASSF9 is absent as a feature, that is NOT evidence it is",
         logfile = LOG)
log_line("  unexpressed. It would mean the processed matrix was filtered, and", logfile = LOG)
log_line("  the test would have to move to the raw data. 37b will distinguish",
         logfile = LOG)
log_line("  these two cases rather than reporting either as absence.",
         logfile = LOG)

report_warnings(logfile = LOG)
log_line("\nFull log: ", LOG, logfile = NULL)
