# ==============================================================================
# 24_visium_qc.R
#
# Quality control and three specific diagnoses on GSE298774, across all
# sections. No classification and no RASSF9 comparison: this script exists to
# find out whether the matrices can be trusted at all.
#
# THREE THINGS TO SETTLE, all raised by the A7 pilot in script 23.
#
# 1. GAPDH READ AS 0.0% OF SPOTS while ACTB read as 39.9%. GAPDH is normally
#    about as abundant as ACTB in any tissue, so a clean zero for one
#    housekeeper and not the other is the signature of a ROW-ALIGNMENT problem
#    between the matrix and the features file -- the matrix row order not
#    matching the features file order would silently corrupt every gene-level
#    number taken from it, including both genes of interest. Diagnosed directly:
#    find GAPDH's row, print its row sum, check for duplicate symbols, and list
#    the most-expressed genes in the matrix. In skin the top of that list should
#    be mitochondrial, ribosomal, keratin and collagen genes. If it is not, the
#    alignment is wrong and nothing else in this dataset can be used yet.
#
# 2. OUT-OF-TISSUE SPOTS. 4,985 spots against a standard Visium grid of 4,992
#    means the matrix holds every spot on the slide, not only those over
#    tissue, which is why median depth looked like 218 UMIs and the minimum was
#    0. tissue_positions.csv carries an in_tissue flag. Everything is
#    recomputed on in-tissue spots only, and the two sets are reported side by
#    side so the effect of the filter is visible rather than assumed.
#
# 3. SECTION A2 has a matrix but NO barcodes and NO features file in the
#    archive, and only a features file at series level. Checked explicitly: a
#    matrix whose rows and columns cannot be named is not usable, and it is
#    better to drop one section knowingly than to mis-join it.
#
# Also recorded for the analysis design: per-section spot counts, depth, and
# detection of both genes and every marker panel, so the choice between a
# spot-level model and a per-compartment pseudobulk rests on numbers.
#
# NOTE: B38 has TWO sections (B38ex and B38_2), so there are 19 sections for 18
# patients. They are kept separate here and must be treated as repeated
# measures on one patient in any model.
#
# Output: results/geo/gse298774_qc_sections.csv
#         results/geo/gse298774_qc_markers.csv
#         results/geo/gse298774_qc_topgenes.csv
#         results/logs/24_visium_qc.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(Matrix))

LOG <- file.path(PATHS$logs, "24_visium_qc.log")
cat("GSE298774 QC —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

EXT_DIR <- file.path(PATHS$raw, "geo", "GSE298774", "extracted")
RES_GEO <- file.path(PATHS$results, "geo")
dir.create(RES_GEO, recursive = TRUE, showWarnings = FALSE)

GENES <- c("RASSF9", "RASSF10")
MARKERS <- list(
  melanocytic  = c("MITF", "PMEL", "MLANA", "DCT", "TYR", "TYRP1", "SOX10"),
  melanoma_up  = c("PRAME", "CDK2", "MKI67"),
  fibroblast   = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"),
  keratinocyte = c("KRT14", "KRT5", "KRT1", "KRT10"),
  immune       = c("PTPRC", "CD3E", "CD68", "LYZ", "CD163"),
  endothelial  = c("PECAM1", "VWF", "CDH5"),
  housekeeping = c("ACTB", "GAPDH", "B2M", "TMSB4X", "EEF1A1")
)

files <- list.files(EXT_DIR, full.names = TRUE)
if (!length(files)) stop("No extracted files in ", EXT_DIR,
                         " -- run script 23 first.")

# Section id = the token between the GSM prefix and the file role.
sec_of <- function(b) sub("^GSM\\d+_(.+?)_(aligned_fiducials|barcodes|features|matrix|scalefactors_json|tissue_hires_image|tissue_lowres_image|tissue_positions).*$", "\\1", b)

inv <- tibble::tibble(path = files, base = basename(files)) |>
  dplyr::mutate(section = sec_of(base),
                role = dplyr::case_when(
                  grepl("_matrix\\.mtx", base)        ~ "matrix",
                  grepl("_barcodes\\.tsv", base)      ~ "barcodes",
                  grepl("_features\\.tsv", base)      ~ "features",
                  grepl("_tissue_positions", base)    ~ "positions",
                  TRUE                                ~ "other")) |>
  dplyr::filter(role != "other")

SECTIONS <- sort(unique(inv$section))
log_line("  sections found: ", length(SECTIONS), " (",
         paste(SECTIONS, collapse = ", "), ")", logfile = LOG)

grab <- function(sec, role) {
  p <- inv$path[inv$section == sec & inv$role == role]
  if (!length(p)) NA_character_ else p[1]
}

# ==============================================================================
# 1. Which sections are complete?
# ==============================================================================
hdr("1. Section completeness")

comp <- tibble::tibble(section = SECTIONS) |>
  dplyr::rowwise() |>
  dplyr::mutate(matrix = !is.na(grab(section, "matrix")),
                barcodes = !is.na(grab(section, "barcodes")),
                features = !is.na(grab(section, "features")),
                positions = !is.na(grab(section, "positions"))) |>
  dplyr::ungroup() |>
  dplyr::mutate(usable = matrix & barcodes & features)

log_line(sprintf("  %-10s %8s %9s %9s %10s %8s", "section", "matrix",
                 "barcodes", "features", "positions", "usable"), logfile = LOG)
for (i in seq_len(nrow(comp))) {
  log_line(sprintf("  %-10s %8s %9s %9s %10s %8s", comp$section[i],
                   comp$matrix[i], comp$barcodes[i], comp$features[i],
                   comp$positions[i],
                   ifelse(comp$usable[i], "yes", "NO")), logfile = LOG)
}
bad <- comp$section[!comp$usable]
if (length(bad)) {
  log_line("", logfile = LOG)
  log_line("  [FLAG] not usable: ", paste(bad, collapse = ", "), logfile = LOG)
  log_line("  A matrix whose rows and columns cannot be named cannot be used.",
           logfile = LOG)
  log_line("  These sections are DROPPED, and that is stated in the report",
           logfile = LOG)
  log_line("  rather than passed over.", logfile = LOG)
}
USE <- comp$section[comp$usable]

# ==============================================================================
# 2. The GAPDH diagnosis, on one section
# ==============================================================================
hdr("2. Row-alignment diagnosis (A7)")

DIAG <- if ("A7" %in% USE) "A7" else USE[1]
log_line("  section: ", DIAG, logfile = LOG)

m <- Matrix::readMM(grab(DIAG, "matrix"))
bc <- readr::read_tsv(grab(DIAG, "barcodes"), col_names = FALSE,
                      show_col_types = FALSE, progress = FALSE)[[1]]
ft <- readr::read_tsv(grab(DIAG, "features"), col_names = FALSE,
                      show_col_types = FALSE, progress = FALSE)
log_line("  matrix ", nrow(m), " x ", ncol(m), " | features ", nrow(ft),
         " | barcodes ", length(bc), logfile = LOG)

if (ncol(m) == nrow(ft) && nrow(m) == length(bc)) {
  log_line("  orientation: spots x genes -> transposing", logfile = LOG)
  m <- Matrix::t(m)
}
stopifnot(nrow(m) == nrow(ft), ncol(m) == length(bc))

sym <- toupper(as.character(ft[[min(2, ncol(ft))]]))
rownames(m) <- make.unique(sym)

# -- duplicate symbols --------------------------------------------------------
dups <- sum(duplicated(sym))
log_line("  duplicate gene symbols in the features file: ", dups, logfile = LOG)

# -- the two housekeepers, row by row -----------------------------------------
log_line("", logfile = LOG)
for (g in c("ACTB", "GAPDH", "B2M", "TMSB4X", "EEF1A1")) {
  idx <- which(sym == g)
  if (!length(idx)) {
    log_line("  ", g, ": not in the features file", logfile = LOG); next
  }
  for (i in idx) {
    rs <- sum(m[i, ])
    log_line(sprintf("  %-8s row %6d  total counts %10.0f  spots positive %5.1f%%",
                     g, i, rs, 100 * mean(m[i, ] > 0)), logfile = LOG)
  }
}

# -- what the matrix says is most abundant ------------------------------------
# The decisive check. In skin, the most-expressed genes should be mitochondrial,
# ribosomal, keratin and collagen. If the top of this list looks arbitrary, the
# matrix rows are not the features file rows.
rs <- Matrix::rowSums(m)
top <- tibble::tibble(gene = rownames(m), total = as.numeric(rs)) |>
  dplyr::arrange(dplyr::desc(total)) |> dplyr::slice_head(n = 25) |>
  dplyr::mutate(section = DIAG)
log_line("", logfile = LOG)
log_line("  25 most-expressed genes in this section:", logfile = LOG)
for (i in seq_len(nrow(top))) {
  log_line(sprintf("    %2d. %-14s %12.0f", i, top$gene[i], top$total[i]),
           logfile = LOG)
}
readr::write_csv(top, file.path(RES_GEO, "gse298774_qc_topgenes.csv"))

plausible <- sum(grepl("^MT-|^RP[LS]|^KRT|^COL|^B2M$|^ACTB$|^EEF|^TMSB|^MALAT1$|^FTL$|^FTH1$",
                       top$gene))
log_line("", logfile = LOG)
log_line("  of the top 25, ", plausible,
         " are mitochondrial, ribosomal, keratin, collagen or housekeeping",
         logfile = LOG)
if (plausible >= 10) {
  log_line("  ROW ALIGNMENT LOOKS CORRECT. The matrix's abundant genes are the",
           logfile = LOG)
  log_line("  ones skin should express, so GAPDH reading low is a property of",
           logfile = LOG)
  log_line("  this assay rather than a join error. Visium 3' capture on FFPE",
           logfile = LOG)
  log_line("  does vary by gene and probe set, and GAPDH is not a reliable",
           logfile = LOG)
  log_line("  reference here -- which is itself worth knowing, since it means",
           logfile = LOG)
  log_line("  ACTB and B2M should carry the quantification control.",
           logfile = LOG)
} else {
  log_line("  [STOP] ROW ALIGNMENT IS SUSPECT. The most-expressed genes are",
           logfile = LOG)
  log_line("  not what skin should express, so the matrix rows probably do not",
           logfile = LOG)
  log_line("  correspond to the features file. Do not read any gene-level",
           logfile = LOG)
  log_line("  number from this dataset until that is resolved.", logfile = LOG)
}

# ==============================================================================
# 3. In-tissue filtering, and what it changes
# ==============================================================================
hdr("3. In-tissue spots")

read_pos <- function(p) {
  if (is.na(p)) return(NULL)
  # Newer spaceranger writes a header; older does not. Detected rather than
  # assumed, because guessing wrong silently shifts every column.
  first <- readr::read_lines(p, n_max = 1)
  has_hdr <- grepl("barcode", first, ignore.case = TRUE)
  tb <- readr::read_csv(p, col_names = has_hdr, show_col_types = FALSE,
                        progress = FALSE)
  if (!has_hdr) {
    names(tb)[seq_len(min(6, ncol(tb)))] <-
      c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")[seq_len(min(6, ncol(tb)))]
  } else {
    names(tb) <- tolower(gsub("[^a-z_]", "", tolower(names(tb))))
  }
  bcol <- grep("barcode", names(tb), value = TRUE)[1]
  icol <- grep("in_?tissue", names(tb), value = TRUE)[1]
  if (is.na(bcol) || is.na(icol)) return(NULL)
  tibble::tibble(barcode = as.character(tb[[bcol]]),
                 in_tissue = as.integer(tb[[icol]]))
}

secq <- list(); mrk <- list()
for (sec in USE) {
  mm <- Matrix::readMM(grab(sec, "matrix"))
  bb <- readr::read_tsv(grab(sec, "barcodes"), col_names = FALSE,
                        show_col_types = FALSE, progress = FALSE)[[1]]
  ff <- readr::read_tsv(grab(sec, "features"), col_names = FALSE,
                        show_col_types = FALSE, progress = FALSE)
  if (ncol(mm) == nrow(ff) && nrow(mm) == length(bb)) mm <- Matrix::t(mm)
  if (nrow(mm) != nrow(ff) || ncol(mm) != length(bb)) {
    log_line("  [FLAG] ", sec, ": dimensions do not match; skipped",
             logfile = LOG)
    next
  }
  sy <- toupper(as.character(ff[[min(2, ncol(ff))]]))

  pos <- read_pos(grab(sec, "positions"))
  keep <- rep(TRUE, length(bb))
  n_pos <- NA_integer_
  if (!is.null(pos)) {
    it <- pos$in_tissue[match(bb, pos$barcode)]
    keep <- !is.na(it) & it == 1L
    n_pos <- sum(keep)
  }

  lib_all <- Matrix::colSums(mm)
  lib_in  <- if (any(keep)) lib_all[keep] else numeric(0)

  secq[[sec]] <- tibble::tibble(
    section = sec,
    spots_total = length(lib_all),
    spots_in_tissue = ifelse(is.na(n_pos), NA_integer_, n_pos),
    median_umi_all = median(lib_all),
    median_umi_in = ifelse(length(lib_in), median(lib_in), NA_real_),
    median_genes_in = ifelse(length(lib_in),
                             median(Matrix::colSums(mm[, keep, drop = FALSE] > 0)),
                             NA_real_))

  sub <- if (any(keep)) mm[, keep, drop = FALSE] else mm
  frac <- function(g) {
    i <- which(sy == g)
    if (!length(i)) return(NA_real_)
    v <- if (length(i) > 1) Matrix::colSums(sub[i, , drop = FALSE]) else sub[i, ]
    mean(v > 0)
  }
  all_g <- c(GENES, unlist(MARKERS, use.names = FALSE))
  pan <- c(rep("gene of interest", length(GENES)),
           rep(names(MARKERS), lengths(MARKERS)))
  mrk[[sec]] <- tibble::tibble(section = sec, gene = all_g, panel = pan) |>
    dplyr::distinct(section, gene, .keep_all = TRUE) |>
    dplyr::mutate(frac_spots = vapply(gene, frac, numeric(1)))

  rm(mm, sub); gc(verbose = FALSE)
}
secq <- dplyr::bind_rows(secq); mrk <- dplyr::bind_rows(mrk)

log_line(sprintf("  %-10s %8s %10s %11s %10s %11s", "section", "spots",
                 "in tissue", "med UMI all", "med UMI in", "med genes in"),
         logfile = LOG)
for (i in seq_len(nrow(secq))) {
  log_line(sprintf("  %-10s %8d %10s %11.0f %10s %11s", secq$section[i],
                   secq$spots_total[i],
                   ifelse(is.na(secq$spots_in_tissue[i]), "?",
                          format(secq$spots_in_tissue[i])),
                   secq$median_umi_all[i],
                   ifelse(is.na(secq$median_umi_in[i]), "?",
                          format(round(secq$median_umi_in[i]))),
                   ifelse(is.na(secq$median_genes_in[i]), "?",
                          format(round(secq$median_genes_in[i])))),
           logfile = LOG)
}
readr::write_csv(secq, file.path(RES_GEO, "gse298774_qc_sections.csv"))

log_line("", logfile = LOG)
log_line("  total in-tissue spots across usable sections: ",
         sum(secq$spots_in_tissue, na.rm = TRUE), logfile = LOG)

# ==============================================================================
# 4. Detection of the genes and the markers, in-tissue
# ==============================================================================
hdr("4. Detection, in-tissue spots only")

summ <- mrk |>
  dplyr::group_by(panel, gene) |>
  dplyr::summarise(median_frac = median(frac_spots, na.rm = TRUE),
                   min_frac = min(frac_spots, na.rm = TRUE),
                   max_frac = max(frac_spots, na.rm = TRUE),
                   .groups = "drop") |>
  dplyr::arrange(factor(panel, levels = c("gene of interest", names(MARKERS))),
                 dplyr::desc(median_frac))

log_line(sprintf("  %-17s %-14s %9s %9s %9s", "panel", "gene", "median",
                 "min", "max"), logfile = LOG)
for (i in seq_len(nrow(summ))) {
  log_line(sprintf("  %-17s %-14s %8.1f%% %8.1f%% %8.1f%%", summ$panel[i],
                   summ$gene[i], 100 * summ$median_frac[i],
                   100 * summ$min_frac[i], 100 * summ$max_frac[i]),
           logfile = LOG)
}
readr::write_csv(mrk, file.path(RES_GEO, "gse298774_qc_markers.csv"))

r9 <- summ$median_frac[summ$gene == "RASSF9"]
log_line("", logfile = LOG)
log_line(sprintf("  RASSF9 median detection across sections: %.1f%% of in-tissue spots",
                 100 * r9), logfile = LOG)
log_line("", logfile = LOG)
if (r9 >= 0.05) {
  log_line("  DESIGN: spot-level modelling is viable. A zero-inflated negative",
           logfile = LOG)
  log_line("  binomial with patient as a random effect, as Kreuger et al. used,",
           logfile = LOG)
  log_line("  is the right model once spots are classified.", logfile = LOG)
} else {
  log_line("  DESIGN: detection per spot is low, so a spot-level test would",
           logfile = LOG)
  log_line("  lean heavily on the zero-inflation model. PSEUDOBULK is the",
           logfile = LOG)
  log_line("  safer primary: sum counts within each compartment within each",
           logfile = LOG)
  log_line("  patient, giving one paired observation per patient and a paired",
           logfile = LOG)
  log_line("  test across 18. Spot-level then becomes the sensitivity",
           logfile = LOG)
  log_line("  analysis, not the headline.", logfile = LOG)
}

report_warnings(logfile = LOG)
write_session_info("24_visium_qc")
log_line("\nFull log: ", LOG, logfile = NULL)
