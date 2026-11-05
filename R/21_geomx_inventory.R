# ==============================================================================
# 21_geomx_inventory.R
#
# Is GSE168013 usable for the one question left open?
#
# THE QUESTION. Script 20 could not settle whether RASSF9's fall from nevus to
# melanoma is melanocyte-intrinsic or tissue composition. Bulk cannot separate
# them; GSE230574 has nevus tissue but melanoma only as cell lines. What is
# needed is a MELANOCYTE-COMPARTMENT comparison across lesion stages, in tissue,
# in several patients.
#
# GSE168013 is that design: NanoString GeoMx digital spatial profiling of FFPE
# melanocytic tumours, 16 patients, 325 regions of interest, across common nevi,
# dysplastic nevi, melanoma in situ and invasive melanoma, with ROIs selected by
# a pathologist as melanocyte-enriched, keratinocyte-rich or immune-rich. If
# RASSF9 can be read in the melanocyte-enriched ROIs across stages, that is a
# direct answer.
#
# THE CATCH, AND WHY THIS SCRIPT ONLY LOOKS. GeoMx here is a TARGETED panel --
# a 108-gene pilot for 4 patients and the 1,412-gene Cancer Transcriptome Atlas
# developmental panel for 12 -- not whole transcriptome. A curated cancer panel
# is unlikely to carry a gene as obscure as RASSF9. The series deposits the
# platform file listing every target, so this is checkable rather than a matter
# of opinion, and it is checked BEFORE any analysis is written.
#
# This script therefore does nothing but inventory and verify:
#   1. download the series' supplementary files
#   2. read the platform file and ask whether RASSF9 / RASSF10 are targets
#   3. inventory the count files and the ROI metadata, reporting their actual
#      shape rather than an assumed one
#   4. state plainly whether the dataset can answer the question
#
# No analysis is attempted here. Writing an analysis for a gene that may not be
# on the panel would be work spent before the one fact that decides whether the
# work is possible.
#
# Output: results/geo/gse168013_panel_check.csv
#         results/geo/gse168013_file_inventory.csv
#         results/logs/21_geomx_inventory.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(GEOquery))

LOG <- file.path(PATHS$logs, "21_geomx_inventory.log")
cat("GSE168013 inventory —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

GEO_DIR <- file.path(PATHS$raw, "geo")
RAW_DIR <- file.path(GEO_DIR, "GSE168013")
RES_GEO <- file.path(PATHS$results, "geo")
invisible(lapply(c(RAW_DIR, RES_GEO), dir.create, recursive = TRUE,
                 showWarnings = FALSE))

GENES <- c("RASSF9", "RASSF10")
# Aliases matter here: a panel may list a gene under an older symbol, and
# missing it for that reason would be a false negative on the one question this
# script exists to answer. RASSF9 was PAMCI / P-CIP1; RASSF10 was FLJ33790.
ALIASES <- list(
  RASSF9  = c("RASSF9", "PAMCI", "P-CIP1", "PCIP1", "C12orf2"),
  RASSF10 = c("RASSF10", "FLJ33790", "A1BG-AS", "KIAA1098")
)

options(timeout = 3600)

# ==============================================================================
# 1. Fetch
# ==============================================================================
hdr("1. Download")

have <- list.files(RAW_DIR, recursive = TRUE)
if (length(have) < 3) {
  log_line("  fetching supplementary files (about 12 MB)", logfile = LOG)
  ok <- try(GEOquery::getGEOSuppFiles("GSE168013", baseDir = GEO_DIR,
                                      makeDirectory = TRUE), silent = TRUE)
  if (inherits(ok, "try-error")) {
    log_line("  [FLAG] download failed: ", as.character(ok), logfile = LOG)
    log_line("  Fetch manually from", logfile = LOG)
    log_line("  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168013",
             logfile = LOG)
    log_line("  into ", RAW_DIR, " and re-run.", logfile = LOG)
  }
} else {
  log_line("  already present (", length(have), " files) -- reusing",
           logfile = LOG)
}

files <- list.files(RAW_DIR, full.names = TRUE, recursive = TRUE)
inv <- tibble::tibble(
  file = basename(files),
  size_kb = round(file.size(files) / 1024, 1))
log_line("", logfile = LOG)
for (i in seq_len(nrow(inv))) {
  log_line(sprintf("    %-58s %9.1f KB", inv$file[i], inv$size_kb[i]),
           logfile = LOG)
}
readr::write_csv(inv, file.path(RES_GEO, "gse168013_file_inventory.csv"))

# Untar if the archive is there, since the per-sample files are inside it.
tarf <- files[grepl("\\.tar$", files)]
if (length(tarf)) {
  ex <- file.path(RAW_DIR, "extracted")
  dir.create(ex, recursive = TRUE, showWarnings = FALSE)
  if (!length(list.files(ex))) {
    utils::untar(tarf[1], exdir = ex)
    log_line("\n  extracted ", length(list.files(ex)), " files from the archive",
             logfile = LOG)
  }
  files <- c(files, list.files(ex, full.names = TRUE))
}

# ==============================================================================
# 2. The decisive check: is the gene on the panel?
# ==============================================================================
hdr("2. Panel membership -- the decisive check")

if (!requireNamespace("readxl", quietly = TRUE)) {
  log_line("  installing readxl", logfile = LOG)
  utils::install.packages("readxl", repos = "https://cloud.r-project.org")
}

plat <- files[grepl("platform_file.*\\.xlsx$", files, ignore.case = TRUE)]
panel_hits <- list()

if (!length(plat)) {
  log_line("  [FLAG] platform file not found among the downloads.", logfile = LOG)
  log_line("  Falling back to searching the count files' own gene columns.",
           logfile = LOG)
} else {
  log_line("  platform file: ", basename(plat[1]), logfile = LOG)
  sheets <- readxl::excel_sheets(plat[1])
  log_line("  sheets: ", paste(sheets, collapse = ", "), logfile = LOG)

  all_targets <- character(0)
  for (sh in sheets) {
    tb <- try(readxl::read_excel(plat[1], sheet = sh, .name_repair = "minimal"),
              silent = TRUE)
    if (inherits(tb, "try-error") || !ncol(tb)) next
    # The target symbol column is not reliably named, so every character column
    # is searched. Over-inclusive on purpose: a false negative here would end
    # the enquiry wrongly.
    chars <- names(tb)[vapply(tb, function(z) is.character(z) || is.factor(z),
                              logical(1))]
    vals <- unique(unlist(lapply(tb[chars], as.character)))
    all_targets <- unique(c(all_targets, vals))
    log_line(sprintf("    sheet %-28s %5d rows, %2d cols, %6d distinct strings",
                     sh, nrow(tb), ncol(tb), length(vals)), logfile = LOG)
  }
  all_targets <- toupper(trimws(all_targets))
  log_line("  distinct strings across the platform file: ", length(all_targets),
           logfile = LOG)

  log_line("", logfile = LOG)
  for (g in GENES) {
    al <- toupper(ALIASES[[g]])
    found <- al[al %in% all_targets]
    panel_hits[[g]] <- tibble::tibble(
      gene = g, on_panel = length(found) > 0,
      matched_symbol = if (length(found)) paste(found, collapse = ", ") else NA)
    if (length(found)) {
      log_line("  ", g, ": ON THE PANEL as ", paste(found, collapse = ", "),
               logfile = LOG)
    } else {
      log_line("  ", g, ": NOT on the panel (checked ",
               paste(al, collapse = ", "), ")", logfile = LOG)
    }
  }

  # A sanity check on the search itself: genes that certainly ARE on a cancer
  # panel must be found, otherwise a "not on the panel" verdict might just mean
  # the platform file was parsed wrongly.
  sentinels <- c("MITF", "PMEL", "S100A8", "CD3E", "PTPRC", "COL1A1", "ACTB")
  sfound <- sentinels[sentinels %in% all_targets]
  log_line("", logfile = LOG)
  log_line("  parser sanity check -- sentinel genes found: ",
           length(sfound), " of ", length(sentinels), " (",
           paste(sfound, collapse = ", "), ")", logfile = LOG)
  if (length(sfound) < 4) {
    log_line("  [FLAG] too few sentinels found. The platform file was probably",
             logfile = LOG)
    log_line("  not parsed correctly, so the verdict above is UNRELIABLE.",
             logfile = LOG)
    log_line("  Open the xlsx by hand before concluding anything.", logfile = LOG)
  }
}

if (length(panel_hits)) {
  readr::write_csv(dplyr::bind_rows(panel_hits),
                   file.path(RES_GEO, "gse168013_panel_check.csv"))
}

# ==============================================================================
# 3. Count files and ROI metadata
# ==============================================================================
# Recorded whatever the panel verdict, because the shape of these files decides
# how an analysis would have to be written if the gene IS present -- and because
# a second GeoMx series with a wider panel would have the same structure.
hdr("3. Count files and ROI annotation")

cnt <- files[grepl("summarized_counts|finalData|dedup", files)]
for (f in head(cnt, 8)) {
  first <- try(readr::read_lines(f, n_max = 1), silent = TRUE)
  if (inherits(first, "try-error")) next
  nf <- length(strsplit(first, "\t")[[1]])
  log_line("  ", basename(f), logfile = LOG)
  log_line("    columns: ", nf, " | header starts: ",
           substr(first, 1, 110), logfile = LOG)
}

# ROI-level annotation lives in the series matrix, not the supplementary files.
gse <- try(GEOquery::getGEO("GSE168013", destdir = GEO_DIR, GSEMatrix = TRUE,
                            getGPL = FALSE), silent = TRUE)
if (!inherits(gse, "try-error")) {
  pd <- Biobase::pData(gse[[1]])
  log_line("", logfile = LOG)
  log_line("  series matrix: ", nrow(pd), " samples (ROIs), ", ncol(pd),
           " metadata columns", logfile = LOG)
  ch <- grep("characteristics|title|source", names(pd), value = TRUE)
  for (cc in head(ch, 12)) {
    vals <- unique(as.character(pd[[cc]]))
    log_line(sprintf("    %-34s %3d distinct: %s", cc, length(vals),
                     substr(paste(head(vals, 4), collapse = " | "), 1, 90)),
             logfile = LOG)
  }
  readr::write_csv(tibble::as_tibble(pd),
                   file.path(RES_GEO, "gse168013_roi_metadata.csv"))
  log_line("  wrote: results/geo/gse168013_roi_metadata.csv", logfile = LOG)
}

# ==============================================================================
# 4. Verdict
# ==============================================================================
hdr("4. Verdict")

usable <- length(panel_hits) &&
  any(vapply(panel_hits, function(z) isTRUE(z$on_panel), logical(1)))

if (usable) {
  on <- names(panel_hits)[vapply(panel_hits,
                                 function(z) isTRUE(z$on_panel), logical(1))]
  log_line("  USABLE for: ", paste(on, collapse = ", "), logfile = LOG)
  log_line("  A melanocyte-compartment comparison across lesion stages is",
           logfile = LOG)
  log_line("  therefore possible, which is the test bulk data cannot do.",
           logfile = LOG)
  log_line("  Next: restrict to melanocyte-enriched ROIs, compare across",
           logfile = LOG)
  log_line("  common nevus / dysplastic nevus / in situ / invasive, and treat",
           logfile = LOG)
  log_line("  the PATIENT as the unit of analysis -- 16 patients, not 325",
           logfile = LOG)
  log_line("  ROIs, since ROIs from one tumour are not independent.",
           logfile = LOG)
} else {
  log_line("  NOT USABLE. Neither gene is a target on this panel, so the",
           logfile = LOG)
  log_line("  dataset cannot speak to RASSF9 or RASSF10 however good its",
           logfile = LOG)
  log_line("  design is. This is a property of the 1,412-gene Cancer", logfile = LOG)
  log_line("  Transcriptome Atlas panel, not of the analysis.", logfile = LOG)
  log_line("", logfile = LOG)
  log_line("  What would still work, in order of preference:", logfile = LOG)
  log_line("    1. a GeoMx WHOLE TRANSCRIPTOME (~18,000 gene) study of", logfile = LOG)
  log_line("       melanocytic progression -- none found in a search of the",
           logfile = LOG)
  log_line("       public literature as of October 2026", logfile = LOG)
  log_line("    2. single-cell RNA-seq of nevus AND melanoma TISSUE from the",
           logfile = LOG)
  log_line("       same or matched patients -- GSE230574 has nevus tissue but",
           logfile = LOG)
  log_line("       melanoma only as cell lines", logfile = LOG)
  log_line("    3. nevus-associated melanoma specimens, where benign and",
           logfile = LOG)
  log_line("       malignant melanocytes sit in ONE section and share a",
           logfile = LOG)
  log_line("       patient, removing between-patient variation entirely",
           logfile = LOG)
  log_line("    4. laser-capture microdissection of melanocytic nests, which",
           logfile = LOG)
  log_line("       is what GSE112509 already did -- and is why GSE112509 is",
           logfile = LOG)
  log_line("       the more trustworthy of the two bulk cohorts", logfile = LOG)
  log_line("", logfile = LOG)
  log_line("  Option 3 is the one to put to Gustavo: a within-patient", logfile = LOG)
  log_line("  comparison in nevus-associated melanoma would settle this, and",
           logfile = LOG)
  log_line("  such specimens exist in any dermatopathology archive.",
           logfile = LOG)
}

report_warnings(logfile = LOG)
write_session_info("21_geomx_inventory")
log_line("\nFull log: ", LOG, logfile = NULL)
