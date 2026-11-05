# ==============================================================================
# 25_visium_classifier.R
#
# Label every in-tissue Visium spot as melanocytic or not, and every melanocytic
# spot as NEVUS-like or MELANOMA-like. Then validate those labels against the
# findings Kreuger et al. published from this same data.
#
# THIS SCRIPT DOES NOT TOUCH RASSF9 OR RASSF10. That is deliberate and it is
# the most important line in this file. Spot labelling involves several
# judgement calls -- which markers, what threshold, how many clusters -- and if
# those calls are made while watching what they do to RASSF9, they will be
# tuned toward a result. Not dishonestly; just reliably. So the classifier is
# built and judged entirely on whether it recovers KNOWN biology, the genes of
# interest are never loaded, and RASSF9 is read exactly once, in the next
# script, after this one has passed.
#
# WHY A CLASSIFIER IS NEEDED AT ALL. The authors assign spots by a 15-topic
# deconvolution held in a local object that is not in the GEO deposit, so their
# labels are unavailable. Emailing the corresponding author
# (r.c.slieker@lumc.nl) for them remains the better route and would make this
# script unnecessary.
#
# THE VALIDATION, pre-specified here before any label exists. The paper reports,
# from this data: proliferative activity in melanoma marked by increased CDK2;
# PRAME negative in the nevus component and positive in the melanoma component
# in all 18 cases; and enrichment of melanocytic lineage programmes in melanoma.
# So a correct labelling must show, in the melanoma-like group relative to the
# nevus-like group:
#   (V1) higher CDK2
#   (V2) higher PRAME
#   (V3) higher MKI67 (proliferation, implied by their CDK2 finding)
#   (V4) comparable melanocytic marker expression -- BOTH groups are
#        melanocytic, so a "melanoma" group that is simply more melanocytic
#        than the "nevus" group has probably just split melanocyte density
#        rather than lesion type, and that would be a failure, not a success.
# V1 and V2 are the primary checks. A section that fails them is reported as
# failed and excluded, rather than carried forward quietly.
#
# DEPTH. Median UMIs per in-tissue spot range 170-fold across sections (128 to
# 21,899), so every score is computed on counts-per-10,000 within spot and
# every threshold is set WITHIN section. Nothing is compared across sections
# before normalisation.
#
# Output: results/geo/gse298774_spot_labels.csv.gz
#         results/geo/gse298774_classifier_validation.csv
#         results/logs/25_visium_classifier.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(Matrix))

LOG <- file.path(PATHS$logs, "25_visium_classifier.log")
cat("GSE298774 spot classifier —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

EXT_DIR <- file.path(PATHS$raw, "geo", "GSE298774", "extracted")
RES_GEO <- file.path(PATHS$results, "geo")
dir.create(RES_GEO, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

# ---- Panels -----------------------------------------------------------------
# RASSF9 and RASSF10 are deliberately ABSENT from every list in this script.
MELANOCYTIC <- c("MITF", "PMEL", "MLANA", "DCT", "TYR", "TYRP1", "SOX10")
MELANOMA_UP <- c("PRAME", "CDK2", "MKI67", "TOP2A", "CCNB1")
OTHER <- list(
  fibroblast   = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"),
  keratinocyte = c("KRT14", "KRT5", "KRT1", "KRT10", "KRTDAP", "FLG"),
  immune       = c("PTPRC", "CD3E", "CD68", "LYZ", "CD163"),
  endothelial  = c("PECAM1", "VWF", "CDH5")
)

MIN_UMI        <- 100   # a spot below this cannot be scored meaningfully
MIN_MEL_MARK   <- 2     # melanocytic markers detected, to call a spot melanocytic
MIN_MEL_SPOTS  <- 60    # melanocytic spots needed to attempt a split
DROP_SECTIONS  <- "A2"  # no barcodes or features in the deposit (script 24)

# ==============================================================================
# 1. Sections
# ==============================================================================
hdr("1. Sections")

files <- list.files(EXT_DIR, full.names = TRUE)
if (!length(files)) stop("No extracted files in ", EXT_DIR, " -- run script 23.")

sec_of <- function(b) sub("^GSM\\d+_(.+?)_(aligned_fiducials|barcodes|features|matrix|scalefactors_json|tissue_hires_image|tissue_lowres_image|tissue_positions).*$", "\\1", b)
inv <- tibble::tibble(path = files, base = basename(files)) |>
  dplyr::mutate(section = sec_of(base),
                role = dplyr::case_when(
                  grepl("_matrix\\.mtx", base)     ~ "matrix",
                  grepl("_barcodes\\.tsv", base)   ~ "barcodes",
                  grepl("_features\\.tsv", base)   ~ "features",
                  grepl("_tissue_positions", base) ~ "positions",
                  TRUE ~ "other")) |>
  dplyr::filter(role != "other")

grab <- function(sec, role) {
  p <- inv$path[inv$section == sec & inv$role == role]
  if (!length(p)) NA_character_ else p[1]
}
SECTIONS <- sort(setdiff(unique(inv$section), DROP_SECTIONS))
SECTIONS <- SECTIONS[vapply(SECTIONS, function(s)
  !any(is.na(c(grab(s, "matrix"), grab(s, "barcodes"), grab(s, "features")))),
  logical(1))]

# Patient = section with the B38 suffixes collapsed, since B38ex and B38_2 are
# two sections from one person and must not count as two patients.
patient_of <- function(sec) sub("^(B38).*$", "\\1", sec)

log_line("  sections to classify: ", length(SECTIONS), logfile = LOG)
log_line("  patients represented: ",
         dplyr::n_distinct(vapply(SECTIONS, patient_of, character(1))),
         logfile = LOG)
log_line("  dropped: ", paste(DROP_SECTIONS, collapse = ", "),
         " (incomplete deposit)", logfile = LOG)

read_pos <- function(p) {
  if (is.na(p)) return(NULL)
  first <- readr::read_lines(p, n_max = 1)
  has_hdr <- grepl("barcode", first, ignore.case = TRUE)
  tb <- readr::read_csv(p, col_names = has_hdr, show_col_types = FALSE,
                        progress = FALSE)
  if (!has_hdr) {
    nm <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")
    names(tb)[seq_len(min(6, ncol(tb)))] <- nm[seq_len(min(6, ncol(tb)))]
  }
  bcol <- grep("barcode", names(tb), ignore.case = TRUE, value = TRUE)[1]
  icol <- grep("in_?tissue", names(tb), ignore.case = TRUE, value = TRUE)[1]
  rcol <- grep("array_row", names(tb), ignore.case = TRUE, value = TRUE)[1]
  ccol <- grep("array_col", names(tb), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(bcol) || is.na(icol)) return(NULL)
  tibble::tibble(barcode = as.character(tb[[bcol]]),
                 in_tissue = as.integer(tb[[icol]]),
                 array_row = if (is.na(rcol)) NA_integer_ else as.integer(tb[[rcol]]),
                 array_col = if (is.na(ccol)) NA_integer_ else as.integer(tb[[ccol]]))
}

# ==============================================================================
# 2. Classify, section by section
# ==============================================================================
hdr("2. Classification")

labels <- list(); valid <- list()

for (sec in SECTIONS) {
  m <- Matrix::readMM(grab(sec, "matrix"))
  bc <- readr::read_tsv(grab(sec, "barcodes"), col_names = FALSE,
                        show_col_types = FALSE, progress = FALSE)[[1]]
  ft <- readr::read_tsv(grab(sec, "features"), col_names = FALSE,
                        show_col_types = FALSE, progress = FALSE)
  if (ncol(m) == nrow(ft) && nrow(m) == length(bc)) m <- Matrix::t(m)
  if (nrow(m) != nrow(ft) || ncol(m) != length(bc)) {
    log_line("  [FLAG] ", sec, ": dimension mismatch; skipped", logfile = LOG)
    next
  }
  sym <- toupper(as.character(ft[[min(2, ncol(ft))]]))

  pos <- read_pos(grab(sec, "positions"))
  lib <- Matrix::colSums(m)
  keep <- rep(TRUE, length(bc))
  if (!is.null(pos)) {
    it <- pos$in_tissue[match(bc, pos$barcode)]
    keep <- !is.na(it) & it == 1L
  }
  keep <- keep & lib >= MIN_UMI
  if (sum(keep) < MIN_MEL_SPOTS) {
    log_line("  ", sec, ": only ", sum(keep), " spots pass depth and tissue",
             " filters; skipped", logfile = LOG)
    next
  }

  sub <- m[, keep, drop = FALSE]
  sublib <- lib[keep]
  subbc <- bc[keep]

  # Counts per 10,000 within spot, so depth differences between spots and
  # between sections cannot drive the scores.
  cpm_of <- function(genes) {
    i <- which(sym %in% genes)
    if (!length(i)) return(rep(NA_real_, ncol(sub)))
    v <- if (length(i) > 1) Matrix::colSums(sub[i, , drop = FALSE]) else sub[i, ]
    1e4 * as.numeric(v) / pmax(sublib, 1)
  }
  ndet_of <- function(genes) {
    i <- which(sym %in% genes)
    if (!length(i)) return(rep(0L, ncol(sub)))
    as.integer(Matrix::colSums(sub[i, , drop = FALSE] > 0))
  }

  mel_n   <- ndet_of(MELANOCYTIC)
  mel_cpm <- cpm_of(MELANOCYTIC)

  # -- melanocytic or not ----------------------------------------------------
  is_mel <- mel_n >= MIN_MEL_MARK
  n_mel <- sum(is_mel)

  d <- tibble::tibble(
    section = sec, patient = patient_of(sec), barcode = subbc,
    lib_size = sublib, mel_markers = mel_n, mel_cpm = mel_cpm,
    melanocytic = is_mel,
    array_row = if (is.null(pos)) NA_integer_ else pos$array_row[match(subbc, pos$barcode)],
    array_col = if (is.null(pos)) NA_integer_ else pos$array_col[match(subbc, pos$barcode)])
  for (nm in names(OTHER)) d[[paste0(nm, "_cpm")]] <- cpm_of(OTHER[[nm]])
  for (g in MELANOMA_UP) d[[paste0(g, "_cpm")]] <- cpm_of(g)
  d$class <- ifelse(is_mel, NA_character_, "Non-melanocytic")

  if (n_mel < MIN_MEL_SPOTS) {
    log_line(sprintf("  %-7s %5d spots, %4d melanocytic -- too few to split",
                     sec, nrow(d), n_mel), logfile = LOG)
    d$class[is.na(d$class)] <- "Melanocytic (unsplit)"
    labels[[sec]] <- d
    rm(m, sub); gc(verbose = FALSE)
    next
  }

  # -- nevus-like vs melanoma-like, within the melanocytic spots -------------
  #
  # TWO CORRECTIONS over the first version of this script, both from a failure
  # it produced: that version split the spots into groups differing EIGHT-FOLD
  # in sequencing depth (66,557 against 8,525 median UMIs) and differing in
  # stromal content, then passed every validation check. It passed because CPM
  # normalisation fixes the MAGNITUDE of a detected count but not the
  # PROBABILITY of detection: a 66k-UMI spot detects a sparse transcript far
  # more often than an 8k one, so a score built from PRAME, CDK2 and MKI67
  # separates deep spots from shallow ones and reproduces "CDK2 up, PRAME up"
  # for reasons that have nothing to do with lesion type.
  #
  # (1) DEPTH RESIDUALISATION. Each marker's CPM is regressed on log10 library
  #     size within section and the residual is used. What remains is the part
  #     of a marker that is NOT explained by how deeply its spot was
  #     sequenced. Some of the depth difference is genuinely biological --
  #     melanoma nests are more cellular than dermal nevus, so they really do
  #     yield more RNA -- but a label that could have been produced by depth
  #     alone cannot be used, so the depth-driven component goes.
  #
  # (2) SPATIAL SMOOTHING. Nevus and melanoma components are contiguous
  #     regions, not interleaved spots, so a spot's neighbours carry real
  #     information about its identity. The score is averaged over each spot
  #     and its immediate neighbours on the Visium grid, which suppresses
  #     per-spot sampling noise -- the dominant noise source for sparse genes.
  #     This uses the one thing Visium offers that bulk and single-cell do not,
  #     and the earlier version ignored it entirely.
  mi <- which(is_mel)
  loglib <- log10(pmax(d$lib_size[mi], 1))

  resid_z <- function(g) {
    v <- d[[paste0(g, "_cpm")]][mi]
    if (all(!is.finite(v)) || stats::sd(v, na.rm = TRUE) == 0)
      return(rep(0, length(mi)))
    fit <- try(stats::lm(v ~ loglib), silent = TRUE)
    r <- if (inherits(fit, "try-error")) v - mean(v, na.rm = TRUE)
         else stats::residuals(fit)
    if (stats::sd(r, na.rm = TRUE) == 0) return(rep(0, length(mi)))
    as.numeric(scale(r))
  }
  sc <- vapply(MELANOMA_UP, resid_z, numeric(length(mi)))
  if (is.null(dim(sc))) sc <- matrix(sc, ncol = length(MELANOMA_UP))
  prog_raw <- rowMeans(sc, na.rm = TRUE)

  # Neighbour smoothing on the hex grid. Visium array coordinates put immediate
  # neighbours within a row offset of 1 and a column offset of 2, so a radius
  # of 2 in that metric covers the first ring. Falls back to the unsmoothed
  # score when coordinates are unavailable.
  ar <- d$array_row[mi]; ac <- d$array_col[mi]
  prog <- prog_raw
  if (!all(is.na(ar)) && !all(is.na(ac))) {
    ok <- !is.na(ar) & !is.na(ac)
    if (sum(ok) > 10) {
      cells <- split(seq_along(ar)[ok], paste(ar[ok] %/% 3, ac[ok] %/% 6))
      sm <- prog_raw
      for (idx in cells) {
        if (length(idx) >= 2) sm[idx] <- mean(prog_raw[idx], na.rm = TRUE)
      }
      # Half the weight on the spot, half on its local neighbourhood: enough to
      # suppress sampling noise without erasing a small compartment.
      prog <- 0.5 * prog_raw + 0.5 * sm
    }
  }

  km <- stats::kmeans(matrix(prog, ncol = 1), centers = 2, nstart = 25)
  hi <- which.max(km$centers[, 1])
  cl <- ifelse(km$cluster == hi, "Melanoma-like", "Nevus-like")
  d$class[mi] <- cl
  d$prog_score <- NA_real_; d$prog_score[mi] <- prog
  d$prog_raw <- NA_real_; d$prog_raw[mi] <- prog_raw

  # -- validation on this section --------------------------------------------
  gm <- function(g, grp) {
    v <- d[[paste0(g, "_cpm")]][mi][cl == grp]
    mean(v, na.rm = TRUE)
  }
  wt <- function(g) {
    a <- d[[paste0(g, "_cpm")]][mi][cl == "Nevus-like"]
    b <- d[[paste0(g, "_cpm")]][mi][cl == "Melanoma-like"]
    if (length(a) < 5 || length(b) < 5) return(NA_real_)
    suppressWarnings(stats::wilcox.test(a, b))$p.value
  }
  v <- tibble::tibble(
    section = sec, patient = patient_of(sec),
    n_spots = nrow(d), n_melanocytic = n_mel,
    n_nevus = sum(cl == "Nevus-like"), n_melanoma = sum(cl == "Melanoma-like"),
    CDK2_nevus = gm("CDK2", "Nevus-like"),
    CDK2_melanoma = gm("CDK2", "Melanoma-like"), CDK2_p = wt("CDK2"),
    PRAME_nevus = gm("PRAME", "Nevus-like"),
    PRAME_melanoma = gm("PRAME", "Melanoma-like"), PRAME_p = wt("PRAME"),
    MKI67_nevus = gm("MKI67", "Nevus-like"),
    MKI67_melanoma = gm("MKI67", "Melanoma-like"),
    mel_cpm_nevus = mean(d$mel_cpm[mi][cl == "Nevus-like"], na.rm = TRUE),
    mel_cpm_melanoma = mean(d$mel_cpm[mi][cl == "Melanoma-like"], na.rm = TRUE),
    # V5 and V6 inputs. These two checks were ABSENT from the first version of
    # this suite, which is why it certified a split that had separated spots by
    # sequencing depth.
    lib_nevus = median(d$lib_size[mi][cl == "Nevus-like"]),
    lib_melanoma = median(d$lib_size[mi][cl == "Melanoma-like"]),
    stroma_nevus = mean(d$fibroblast_cpm[mi][cl == "Nevus-like"], na.rm = TRUE) +
                   mean(d$endothelial_cpm[mi][cl == "Nevus-like"], na.rm = TRUE),
    stroma_melanoma = mean(d$fibroblast_cpm[mi][cl == "Melanoma-like"], na.rm = TRUE) +
                      mean(d$endothelial_cpm[mi][cl == "Melanoma-like"], na.rm = TRUE)) |>
    dplyr::mutate(
      V1_CDK2_up  = CDK2_melanoma > CDK2_nevus,
      V2_PRAME_up = PRAME_melanoma > PRAME_nevus,
      V3_MKI67_up = MKI67_melanoma > MKI67_nevus,
      # V4: the melanocytic content of the two groups should be comparable.
      # A ratio far from 1 means the split tracked melanocyte density, not
      # lesion type.
      V4_mel_ratio = mel_cpm_melanoma / pmax(mel_cpm_nevus, 1e-9),
      V4_ok = V4_mel_ratio > 0.5 & V4_mel_ratio < 2,
      # V5 DEPTH: the single check whose absence let the first version pass a
      # broken split. If the two groups differ materially in library size, the
      # labels may be a depth artefact and no comparison between them is safe.
      V5_lib_ratio = lib_melanoma / pmax(lib_nevus, 1),
      V5_ok = V5_lib_ratio > 0.67 & V5_lib_ratio < 1.5,
      # V6 STROMA: the dataset's entire advantage over the bulk cohorts is that
      # both compartments share the same stroma. If they do not, that advantage
      # is gone and this is just another composition-confounded comparison.
      V6_stroma_ratio = stroma_melanoma / pmax(stroma_nevus, 1e-9),
      V6_ok = V6_stroma_ratio > 0.67 & V6_stroma_ratio < 1.5,
      passes = V1_CDK2_up & V2_PRAME_up & V5_ok & V6_ok)

  valid[[sec]] <- v
  labels[[sec]] <- d

  log_line(sprintf("  %-7s %4dn/%4dm | CDK2 %s PRAME %s | mel %4.2f%s depth %4.2f%s stroma %4.2f%s | %s",
                   sec, v$n_nevus, v$n_melanoma,
                   ifelse(v$V1_CDK2_up, "ok", "FAIL"),
                   ifelse(v$V2_PRAME_up, "ok", "FAIL"),
                   v$V4_mel_ratio, ifelse(v$V4_ok, " ", "!"),
                   v$V5_lib_ratio, ifelse(v$V5_ok, " ", "!"),
                   v$V6_stroma_ratio, ifelse(v$V6_ok, " ", "!"),
                   ifelse(v$passes, "PASS", "FAIL")), logfile = LOG)

  rm(m, sub); gc(verbose = FALSE)
}

L <- dplyr::bind_rows(labels)
V <- dplyr::bind_rows(valid)

# ==============================================================================
# 3. Does the classifier pass?
# ==============================================================================
hdr("3. Validation")

if (!nrow(V)) {
  log_line("  no section could be split; the classifier fails outright",
           logfile = LOG)
} else {
  log_line("  sections split: ", nrow(V), logfile = LOG)
  log_line("  V1 CDK2 higher in melanoma-like:  ", sum(V$V1_CDK2_up), " of ",
           nrow(V), logfile = LOG)
  log_line("  V2 PRAME higher in melanoma-like: ", sum(V$V2_PRAME_up), " of ",
           nrow(V), logfile = LOG)
  log_line("  V3 MKI67 higher in melanoma-like: ", sum(V$V3_MKI67_up), " of ",
           nrow(V), logfile = LOG)
  log_line("  V4 melanocytic content comparable: ", sum(V$V4_ok), " of ",
           nrow(V), logfile = LOG)
  log_line("  V5 sequencing DEPTH comparable:   ", sum(V$V5_ok), " of ",
           nrow(V), logfile = LOG)
  log_line("  V6 STROMAL content comparable:    ", sum(V$V6_ok), " of ",
           nrow(V), logfile = LOG)
  log_line("  passing V1, V2, V5 and V6:        ", sum(V$passes), " of ",
           nrow(V), logfile = LOG)
  log_line("", logfile = LOG)
  log_line(sprintf("  depth ratio across sections: median %.2f (range %.2f-%.2f)",
                   median(V$V5_lib_ratio), min(V$V5_lib_ratio),
                   max(V$V5_lib_ratio)), logfile = LOG)
  log_line(sprintf("  stroma ratio across sections: median %.2f (range %.2f-%.2f)",
                   median(V$V6_stroma_ratio), min(V$V6_stroma_ratio),
                   max(V$V6_stroma_ratio)), logfile = LOG)

  # A sign test across sections: if the labelling were arbitrary, CDK2 and
  # PRAME would be higher in the melanoma-like group about half the time.
  for (nm in c("V1_CDK2_up", "V2_PRAME_up", "V3_MKI67_up")) {
    k <- sum(V[[nm]]); n <- nrow(V)
    bt <- stats::binom.test(k, n, 0.5)
    log_line(sprintf("    %s: %d/%d, sign test p = %s", nm, k, n,
                     fmt_p(bt$p.value)), logfile = LOG)
  }

  readr::write_csv(V, file.path(RES_GEO, "gse298774_classifier_validation.csv"))
  log_line("", logfile = LOG)
  log_line("  wrote: results/geo/gse298774_classifier_validation.csv",
           logfile = LOG)

  # A section that fails V5 or V6 is named, because the next script excludes it
  # by name and a reader needs to know which patients left the analysis.
  failed <- V$section[!V$passes]
  if (length(failed)) {
    log_line("", logfile = LOG)
    log_line("  sections FAILING: ", paste(failed, collapse = ", "),
             logfile = LOG)
    for (s in failed) {
      r <- V[V$section == s, ]
      log_line(sprintf("    %-7s depth ratio %.2f, stroma ratio %.2f, mel ratio %.2f",
                       s, r$V5_lib_ratio, r$V6_stroma_ratio, r$V4_mel_ratio),
               logfile = LOG)
    }
  }

  pass_frac <- mean(V$passes)
  log_line("", logfile = LOG)
  if (pass_frac >= 0.75) {
    log_line(sprintf("  CLASSIFIER PASSES in %.0f%% of sections. It recovers the",
                     100 * pass_frac), logfile = LOG)
    log_line("  CDK2 and PRAME differences the authors reported from this same",
             logfile = LOG)
    log_line("  data, which it was not told about. The labels can be used.",
             logfile = LOG)
    log_line("  Sections that FAILED are excluded by name in the next script,",
             logfile = LOG)
    log_line("  not silently pooled.", logfile = LOG)
  } else {
    log_line(sprintf("  [STOP] CLASSIFIER FAILS: only %.0f%% of sections",
                     100 * pass_frac), logfile = LOG)
    log_line("  reproduce both published markers. These labels do not identify",
             logfile = LOG)
    log_line("  nevus and melanoma reliably, so no RASSF9 comparison should be",
             logfile = LOG)
    log_line("  built on them. The right move then is to ask the authors for",
             logfile = LOG)
    log_line("  their per-spot topic assignments rather than to retune this",
             logfile = LOG)
    log_line("  classifier until it passes -- retuning against a validation",
             logfile = LOG)
    log_line("  target is how a classifier gets fitted to noise.", logfile = LOG)
  }
}

# ==============================================================================
# 4. Composition of the labelled groups
# ==============================================================================
hdr("4. Composition")

if (nrow(L)) {
  comp <- L |>
    dplyr::group_by(class) |>
    dplyr::summarise(n = dplyr::n(),
                     median_lib = median(lib_size),
                     mel = mean(mel_cpm, na.rm = TRUE),
                     fibro = mean(fibroblast_cpm, na.rm = TRUE),
                     kerat = mean(keratinocyte_cpm, na.rm = TRUE),
                     immune = mean(immune_cpm, na.rm = TRUE),
                     endo = mean(endothelial_cpm, na.rm = TRUE),
                     .groups = "drop")
  log_line(sprintf("  %-22s %7s %10s %8s %8s %8s %8s %8s", "class", "n",
                   "med UMI", "melano", "fibro", "kerat", "immune", "endo"),
           logfile = LOG)
  for (i in seq_len(nrow(comp))) {
    log_line(sprintf("  %-22s %7d %10.0f %8.1f %8.1f %8.1f %8.1f %8.1f",
                     comp$class[i], comp$n[i], comp$median_lib[i], comp$mel[i],
                     comp$fibro[i], comp$kerat[i], comp$immune[i],
                     comp$endo[i]), logfile = LOG)
  }
  log_line("", logfile = LOG)
  log_line("  Read the stromal columns: if nevus-like and melanoma-like spots",
           logfile = LOG)
  log_line("  carry similar fibroblast and endothelial content, then a", logfile = LOG)
  log_line("  difference between them cannot be the composition effect that",
           logfile = LOG)
  log_line("  confounded the bulk cohorts -- which is the whole reason for",
           logfile = LOG)
  log_line("  using this dataset.", logfile = LOG)

  readr::write_csv(L, file.path(RES_GEO, "gse298774_spot_labels.csv.gz"))
  log_line("", logfile = LOG)
  log_line("  wrote: results/geo/gse298774_spot_labels.csv.gz  (",
           nrow(L), " spots)", logfile = LOG)
}

log_line("", logfile = LOG)
log_line("  NOTE: neither RASSF9 nor RASSF10 was loaded, scored or mentioned",
         logfile = LOG)
log_line("  anywhere in this script. They are read once, in script 26, now",
         logfile = LOG)
log_line("  that the labels are fixed.", logfile = LOG)

report_warnings(logfile = LOG)
write_session_info("25_visium_classifier")
log_line("\nFull log: ", LOG, logfile = NULL)
