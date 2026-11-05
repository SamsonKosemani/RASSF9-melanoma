# ==============================================================================
# 23_visium_structure.R
#
# Download GSE298774 and establish exactly what the deposit contains, then load
# ONE single-patient slide end to end to prove the pipeline works.
#
# WHAT THE GEO RECORD SAYS. 11 samples (GSM9022963-GSM9022973), each a SLIDE
# carrying one to three patients' lesions, covering the 18 Visium IDs in the
# paper's Table S1. Supplement: GSE298774_RAW.tar, 1.3 GB, holding CSV, JPG,
# JSON, MTX, PNG and TSV -- the shape of standard spaceranger output. One file
# sits at series level, GSE298774_A2_features.tsv.gz, named by VISIUM ID rather
# than by GSM, which suggests the per-section files are already split by Visium
# ID. If that holds, patient assignment is done for us.
#
# TWO THINGS THIS SCRIPT SETTLES.
#   1. The file layout: are matrices per Visium ID (so per patient section) or
#      per slide? Is there any per-spot annotation, or only tissue positions?
#   2. Whether RASSF9 and RASSF10 are in the feature list at all. Visium is
#      whole-transcriptome so they should be, but it is one grep and it decides
#      whether anything else is worth doing.
#
# THEN IT LOADS ONE SLIDE. GSM9022964 (patient A7) carries a single patient, so
# it needs no section separation. Loading it proves the reader works and gives
# the spot count, depth distribution and marker detection we need before
# committing to all 11.
#
# NO CLASSIFICATION AND NO RASSF9 COMPARISON HERE. Spot labels have to be
# derived (their 15-topic assignment is not in the deposit) and a derived
# classifier must be validated against the paper's own published findings
# before any RASSF9 number is read off it. That is the next script. Reading a
# result off an unvalidated classifier is how this analysis would go wrong for
# a fifth time.
#
# Output: results/geo/gse298774_tar_contents.csv
#         results/geo/gse298774_section_inventory.csv
#         results/geo/gse298774_A7_spot_qc.csv
#         results/logs/23_visium_structure.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(GEOquery))

LOG <- file.path(PATHS$logs, "23_visium_structure.log")
cat("GSE298774 structure —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

GEO_DIR <- file.path(PATHS$raw, "geo")
RAW_DIR <- file.path(GEO_DIR, "GSE298774")
EXT_DIR <- file.path(RAW_DIR, "extracted")
RES_GEO <- file.path(PATHS$results, "geo")
invisible(lapply(c(RAW_DIR, EXT_DIR, RES_GEO), dir.create, recursive = TRUE,
                 showWarnings = FALSE))

GENES <- c("RASSF9", "RASSF10")

# Markers for the QC readout only. The classifier is built in the next script.
MARKERS <- list(
  melanocytic  = c("MITF", "PMEL", "MLANA", "DCT", "TYR", "TYRP1", "SOX10"),
  melanoma_up  = c("PRAME", "CDK2", "MKI67"),   # their reported discriminators
  fibroblast   = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"),
  keratinocyte = c("KRT14", "KRT5", "KRT1", "KRT10"),
  immune       = c("PTPRC", "CD3E", "CD68", "LYZ", "CD163"),
  endothelial  = c("PECAM1", "VWF", "CDH5"),
  housekeeping = c("ACTB", "GAPDH")
)

options(timeout = 7200)

# ==============================================================================
# 1. Download
# ==============================================================================
hdr("1. Download")

tarf <- file.path(RAW_DIR, "GSE298774_RAW.tar")
if (!file.exists(tarf) || file.size(tarf) < 1e9) {
  log_line("  fetching GSE298774_RAW.tar (1.3 GB) -- this takes a while",
           logfile = LOG)
  ok <- try(GEOquery::getGEOSuppFiles("GSE298774", baseDir = GEO_DIR,
                                      makeDirectory = TRUE), silent = TRUE)
  if (inherits(ok, "try-error")) {
    log_line("  [FLAG] download failed: ", as.character(ok), logfile = LOG)
    log_line("  Fetch GSE298774_RAW.tar by hand into ", RAW_DIR, logfile = LOG)
    log_line("  from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298774",
             logfile = LOG)
    stop("Download failed; see the log.")
  }
  tarf <- list.files(RAW_DIR, pattern = "RAW\\.tar$", full.names = TRUE)[1]
} else {
  log_line("  already downloaded: ", basename(tarf), " (",
           round(file.size(tarf) / 1024^3, 2), " GB)", logfile = LOG)
}

# ==============================================================================
# 2. What is inside the archive
# ==============================================================================
# Listed before extraction, so the layout is known before 1.3 GB is unpacked.
hdr("2. Archive contents")

contents <- utils::untar(tarf, list = TRUE)
log_line("  files in the archive: ", length(contents), logfile = LOG)

cinv <- tibble::tibble(path = contents, base = basename(contents)) |>
  dplyr::mutate(
    ext = tolower(sub(".*\\.", "", sub("\\.gz$", "", base))),
    kind = dplyr::case_when(
      grepl("matrix", base, ignore.case = TRUE)            ~ "matrix",
      grepl("barcode", base, ignore.case = TRUE)           ~ "barcodes",
      grepl("feature|genes", base, ignore.case = TRUE)     ~ "features",
      grepl("tissue_position|positions", base, ignore.case = TRUE) ~ "positions",
      grepl("scalefactor", base, ignore.case = TRUE)       ~ "scalefactors",
      grepl("\\.png$|\\.jpg$|\\.jpeg$|image", base, ignore.case = TRUE) ~ "image",
      grepl("annot|cluster|metad|label|topic|deconv", base,
            ignore.case = TRUE)                            ~ "ANNOTATION",
      TRUE                                                 ~ "other"))

log_line("", logfile = LOG)
tab <- cinv |> dplyr::count(kind, name = "n") |> dplyr::arrange(dplyr::desc(n))
for (i in seq_len(nrow(tab))) {
  log_line(sprintf("    %-14s %4d", tab$kind[i], tab$n[i]), logfile = LOG)
}
readr::write_csv(cinv, file.path(RES_GEO, "gse298774_tar_contents.csv"))

log_line("", logfile = LOG)
log_line("  first 20 names as deposited:", logfile = LOG)
for (b in utils::head(cinv$base, 20)) log_line("    ", b, logfile = LOG)

# -- the decisive question about annotation ------------------------------------
ann <- cinv |> dplyr::filter(kind == "ANNOTATION")
log_line("", logfile = LOG)
if (nrow(ann)) {
  log_line("  PER-SPOT ANNOTATION FOUND (", nrow(ann), " file(s)):", logfile = LOG)
  for (b in ann$base) log_line("    ", b, logfile = LOG)
  log_line("  If these carry the authors' nevus/tumour calls, we use THEM and",
           logfile = LOG)
  log_line("  skip deriving a classifier entirely -- far better.", logfile = LOG)
} else {
  log_line("  NO per-spot annotation in the archive. Only spaceranger output,",
           logfile = LOG)
  log_line("  so spot labels must be derived and then validated against the",
           logfile = LOG)
  log_line("  paper's published findings before use.", logfile = LOG)
  log_line("", logfile = LOG)
  log_line("  WORTH DOING IN PARALLEL: email the corresponding author",
           logfile = LOG)
  log_line("  (r.c.slieker@lumc.nl) for the per-spot topic assignments from",
           logfile = LOG)
  log_line("  Semla_Data_Deconvoluted.RData. Authors commonly share these, it",
           logfile = LOG)
  log_line("  costs one email, and it would remove the single largest source",
           logfile = LOG)
  log_line("  of error in this analysis.", logfile = LOG)
}

# -- how sections are keyed ----------------------------------------------------
# The series-level file was GSE298774_A2_features.tsv.gz, so the prefix before
# the file role is probably the Visium ID. Extracted here rather than assumed.
VIS_IDS <- c("B1","A2","D3","B4","A7","B8","B13","B14","B15","A16","A17",
             "D20","A21","A22","A23","A24","A34","B38")
cinv$visium_id <- NA_character_
for (v in VIS_IDS) {
  hit <- grepl(paste0("(^|[^A-Za-z0-9])", v, "([^0-9]|$)"), cinv$base)
  cinv$visium_id[hit & is.na(cinv$visium_id)] <- v
}
log_line("", logfile = LOG)
log_line("  files whose name carries a Visium ID: ",
         sum(!is.na(cinv$visium_id)), " of ", nrow(cinv), logfile = LOG)

sect <- cinv |> dplyr::filter(!is.na(visium_id)) |>
  dplyr::count(visium_id, kind) |>
  tidyr::pivot_wider(names_from = kind, values_from = n, values_fill = 0) |>
  dplyr::arrange(visium_id)
if (nrow(sect)) {
  log_line("", logfile = LOG)
  log_line("  per-section file counts:", logfile = LOG)
  cn <- setdiff(names(sect), "visium_id")
  log_line(sprintf("    %-10s %s", "visium",
                   paste(sprintf("%12s", cn), collapse = " ")), logfile = LOG)
  for (i in seq_len(nrow(sect))) {
    log_line(sprintf("    %-10s %s", sect$visium_id[i],
                     paste(sprintf("%12s", unlist(sect[i, cn])),
                           collapse = " ")), logfile = LOG)
  }
  readr::write_csv(sect, file.path(RES_GEO, "gse298774_section_inventory.csv"))
  log_line("", logfile = LOG)
  log_line("  -> matrices are keyed by VISIUM ID, so each section is already",
           logfile = LOG)
  log_line("     attributable to one patient and no geometric separation of",
           logfile = LOG)
  log_line("     multi-patient slides is needed.", logfile = LOG)
} else {
  log_line("", logfile = LOG)
  log_line("  [FLAG] no Visium IDs in the file names -- files are keyed by GSM,",
           logfile = LOG)
  log_line("  so a slide holding two or three patients will have to be split",
           logfile = LOG)
  log_line("  geometrically and matched to patients using the deposited", logfile = LOG)
  log_line("  images. That is a substantial extra step; report before doing it.",
           logfile = LOG)
}

# ==============================================================================
# 3. Extract
# ==============================================================================
hdr("3. Extract")

if (length(list.files(EXT_DIR, recursive = TRUE)) < length(contents) / 2) {
  utils::untar(tarf, exdir = EXT_DIR)
  log_line("  extracted into data/raw/geo/GSE298774/extracted/", logfile = LOG)
} else {
  log_line("  already extracted -- reusing", logfile = LOG)
}
files <- list.files(EXT_DIR, full.names = TRUE, recursive = TRUE)
log_line("  files on disk: ", length(files), logfile = LOG)

# ==============================================================================
# 4. Are the two genes in the feature list?
# ==============================================================================
hdr("4. Gene presence")

featf <- c(list.files(RAW_DIR, pattern = "features.*tsv", full.names = TRUE),
           files[grepl("features.*tsv|genes.*tsv", files, ignore.case = TRUE)])
featf <- featf[file.exists(featf)]

if (!length(featf)) {
  log_line("  [FLAG] no features/genes file found", logfile = LOG)
} else {
  f1 <- featf[1]
  log_line("  reading ", basename(f1), logfile = LOG)
  ft <- readr::read_tsv(f1, col_names = FALSE, show_col_types = FALSE,
                        progress = FALSE)
  log_line("  features: ", nrow(ft), " rows, ", ncol(ft), " columns",
           logfile = LOG)
  sym_col <- if (ncol(ft) >= 2) 2 else 1
  syms <- toupper(as.character(ft[[sym_col]]))

  log_line("", logfile = LOG)
  for (g in GENES) {
    log_line("  ", g, ": ", if (g %in% syms) "PRESENT" else "ABSENT",
             logfile = LOG)
  }
  log_line("", logfile = LOG)
  for (nm in names(MARKERS)) {
    present <- MARKERS[[nm]][MARKERS[[nm]] %in% syms]
    log_line(sprintf("    %-13s %d of %d present%s", nm, length(present),
                     length(MARKERS[[nm]]),
                     if (length(present) < length(MARKERS[[nm]]))
                       paste0(" (missing: ",
                              paste(setdiff(MARKERS[[nm]], present),
                                    collapse = ", "), ")") else ""),
             logfile = LOG)
  }
  if (!all(GENES %in% syms)) {
    log_line("", logfile = LOG)
    log_line("  [STOP] a gene of interest is not in the feature list. Nothing",
             logfile = LOG)
    log_line("  further is possible for it in this dataset.", logfile = LOG)
  }
}

# ==============================================================================
# 5. Load one single-patient section end to end
# ==============================================================================
# A7, A21, A34 and B4 are the sections whose slide carries a single patient.
# A7 is loaded here: it proves the reader works and gives the spot count and
# depth distribution before all 11 slides are committed to.
hdr("5. One section loaded (A7)")

PILOT <- "A7"

find_one <- function(vid, pattern) {
  hits <- files[grepl(paste0("(^|[^A-Za-z0-9])", vid, "([^0-9]|$)"),
                      basename(files)) &
                grepl(pattern, basename(files), ignore.case = TRUE)]
  if (!length(hits)) NA_character_ else hits[1]
}

mtx <- find_one(PILOT, "matrix.*mtx")
bcs <- find_one(PILOT, "barcode")
fts <- find_one(PILOT, "features|genes")
pos <- find_one(PILOT, "tissue_position|positions")

log_line("  matrix:    ", ifelse(is.na(mtx), "NOT FOUND", basename(mtx)),
         logfile = LOG)
log_line("  barcodes:  ", ifelse(is.na(bcs), "NOT FOUND", basename(bcs)),
         logfile = LOG)
log_line("  features:  ", ifelse(is.na(fts), "NOT FOUND", basename(fts)),
         logfile = LOG)
log_line("  positions: ", ifelse(is.na(pos), "NOT FOUND", basename(pos)),
         logfile = LOG)

if (any(is.na(c(mtx, bcs, fts)))) {
  log_line("", logfile = LOG)
  log_line("  [FLAG] the triplet for ", PILOT, " is incomplete. Read section 2",
           logfile = LOG)
  log_line("  above for the actual naming, and the pilot will need its file",
           logfile = LOG)
  log_line("  patterns adjusted before anything is loaded.", logfile = LOG)
} else {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    utils::install.packages("Matrix", repos = "https://cloud.r-project.org")
  }
  m <- Matrix::readMM(mtx)
  bc <- readr::read_tsv(bcs, col_names = FALSE, show_col_types = FALSE,
                        progress = FALSE)[[1]]
  ft <- readr::read_tsv(fts, col_names = FALSE, show_col_types = FALSE,
                        progress = FALSE)
  sym <- toupper(as.character(ft[[min(2, ncol(ft))]]))

  # Orientation is checked rather than assumed: MTX from spaceranger is
  # genes x spots, but a transposed deposit would silently ruin everything.
  if (nrow(m) == length(sym) && ncol(m) == length(bc)) {
    log_line("  matrix is genes x spots: ", nrow(m), " x ", ncol(m),
             logfile = LOG)
  } else if (ncol(m) == length(sym) && nrow(m) == length(bc)) {
    log_line("  matrix is spots x genes; transposing", logfile = LOG)
    m <- Matrix::t(m)
  } else {
    log_line("  [FLAG] dimensions do not match the barcode and feature files:",
             logfile = LOG)
    log_line("    matrix ", nrow(m), " x ", ncol(m), ", features ", length(sym),
             ", barcodes ", length(bc), logfile = LOG)
    stop("Dimension mismatch; see the log.")
  }

  lib <- Matrix::colSums(m)
  log_line("", logfile = LOG)
  log_line("  spots: ", length(lib), logfile = LOG)
  log_line(sprintf("  UMIs per spot: median %.0f, IQR %.0f-%.0f, range %.0f-%.0f",
                   median(lib), quantile(lib, .25), quantile(lib, .75),
                   min(lib), max(lib)), logfile = LOG)
  log_line("  genes detected per spot: median ",
           round(median(Matrix::colSums(m > 0))), logfile = LOG)

  # Detection of the two genes and of every marker, as a fraction of spots.
  frac_of <- function(g) {
    i <- which(sym == g)
    if (!length(i)) return(NA_real_)
    v <- if (length(i) > 1) Matrix::colSums(m[i, , drop = FALSE]) else m[i, ]
    mean(v > 0)
  }
  qc <- tibble::tibble(
    gene = c(GENES, unlist(MARKERS, use.names = FALSE)),
    panel = c(rep("gene of interest", length(GENES)),
              rep(names(MARKERS), lengths(MARKERS)))) |>
    dplyr::distinct(gene, .keep_all = TRUE) |>
    dplyr::mutate(frac_spots_positive = vapply(gene, frac_of, numeric(1)),
                  section = PILOT)

  log_line("", logfile = LOG)
  log_line("  percent of spots positive:", logfile = LOG)
  for (i in seq_len(nrow(qc))) {
    log_line(sprintf("    %-13s %-18s %6.1f%%", qc$panel[i], qc$gene[i],
                     100 * qc$frac_spots_positive[i]), logfile = LOG)
  }
  readr::write_csv(qc, file.path(RES_GEO, "gse298774_A7_spot_qc.csv"))
  log_line("", logfile = LOG)
  log_line("  wrote: results/geo/gse298774_A7_spot_qc.csv", logfile = LOG)

  r9 <- qc$frac_spots_positive[qc$gene == "RASSF9"]
  log_line("", logfile = LOG)
  if (!is.na(r9) && r9 >= 0.01) {
    log_line(sprintf("  RASSF9 is detected in %.1f%% of spots in this section.",
                     100 * r9), logfile = LOG)
    log_line("  That is enough to support a spot-level comparison once spots",
             logfile = LOG)
    log_line("  are classified.", logfile = LOG)
  } else {
    log_line(sprintf("  RASSF9 is detected in only %.2f%% of spots.",
                     100 * ifelse(is.na(r9), 0, r9)), logfile = LOG)
    log_line("  A spot-level comparison may be underpowered even across 18",
             logfile = LOG)
    log_line("  sections; pseudobulk per compartment per patient would then be",
             logfile = LOG)
    log_line("  the better design. Decide from this number, not from hope.",
             logfile = LOG)
  }
}

# ==============================================================================
# 6. Next
# ==============================================================================
hdr("6. Next")

log_line("  Three things decide the next script:", logfile = LOG)
log_line("    1. whether any per-spot annotation exists (section 2)", logfile = LOG)
log_line("    2. whether the files are keyed by Visium ID (section 2)", logfile = LOG)
log_line("    3. how often RASSF9 is detected per spot (section 5)", logfile = LOG)
log_line("", logfile = LOG)
log_line("  If detection is low, the analysis becomes pseudobulk: sum counts",
         logfile = LOG)
log_line("  within each compartment within each patient, giving 18 paired", logfile = LOG)
log_line("  observations and a paired test -- far less powerful per spot but",
         logfile = LOG)
log_line("  robust, and the pairing is what makes it work.", logfile = LOG)

report_warnings(logfile = LOG)
write_session_info("23_visium_structure")
log_line("\nFull log: ", LOG, logfile = NULL)
