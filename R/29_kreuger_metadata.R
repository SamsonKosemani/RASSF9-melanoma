# ==============================================================================
# 29_kreuger_metadata.R
#
# Parse and verify the per-spot metadata Roderick Slieker sent on 14 Oct 2026,
# and build the slim join table the RASSF9 readout will use.
#
# WHAT THIS FILE IS. 39,481 spots x 80 columns, tab-separated. It carries, for
# every spot in their Seurat object:
#   - the published STdeconvolve topic proportions, X1 to X15 -- their actual
#     theta, not a reproduction
#   - MaxTopics, the winning topic, and MaxTopicsName, their published name for
#     it (Nevus, Tumor 1, Tumor 2, Fibroblasts 1/2, Keratinocytes 1/2/3,
#     Immune cells 1/2, Endothelial cells, Eccrine gland, Eccrine gland duct,
#     Sebaceous gland, Hair follicle)
#   - the Kunz, Badal and Shain nevus/tumour signature scores
#   - three independent reference-based cell-type scorings
#   - per-spot QC (nCount_Spatial, nFeature_Spatial, percent.mito/ribo)
#
# This supersedes our own topic reproduction as the SOURCE OF LABELS. Scripts
# 26b and 27 become a validation exercise instead: we can now ask how closely
# our independent K=15 fit agrees with theirs, which is a real check on whether
# the pipeline was faithfully reproduced and was not available before.
#
# TWO THINGS THIS SCRIPT EXISTS TO SETTLE, both of which would silently corrupt
# the readout if guessed wrong:
#
#   1. COLUMN ALIGNMENT. The header has 80 names. If the file was written with
#      write.table(row.names = TRUE) the data rows carry 81 fields and every
#      column is shifted by one -- which would put topic proportions under the
#      wrong names and go completely unnoticed. So field counts are compared
#      directly, and the parse is only accepted if an independent sanity check
#      passes: X1..X15 must be numeric, non-negative, and sum to 1 per spot.
#
#   2. THE JOIN KEY. Visium barcodes are drawn from a fixed whitelist, so the
#      SAME barcode occurs in every section. A join on barcode alone would
#      scatter each spot's label across 18 sections. The key must be
#      section + barcode, and this script identifies which column holds the
#      section by matching its values against the section names we already have,
#      rather than assuming it is 'sample' or 'library_id'.
#
# AND ONE THING IT REPORTS RATHER THAN SOLVES: their 39,481 spots against our
# 61,734 in-tissue spots. The shortfall is stated per section, because whether
# it is uniform QC filtering or whole sections missing changes what the readout
# can claim.
#
# NEITHER RASSF9 NOR RASSF10 APPEARS IN THIS SCRIPT.
#
# Output: data/processed/gse298774_kreuger_labels.rds  (slim join table)
#         results/geo/gse298774_kreuger_coverage.csv
#         results/geo/gse298774_kreuger_label_counts.csv
#         results/logs/29_kreuger_metadata.log
# ==============================================================================

source(here::here("R", "00_setup.R"))

LOG <- file.path(PATHS$logs, "29_kreuger_metadata.log")
cat("GSE298774 Kreuger metadata —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

META_F <- file.path(PATHS$raw, "geo", "GSE298774", "kreuger_metadata",
                    "MetaData_Visium_Kreuger.txt")
PROC    <- PATHS$processed
RES_GEO <- file.path(PATHS$results, "geo")

TOPIC_COLS <- paste0("X", 1:15)

# ==============================================================================
# 1. Column alignment
# ==============================================================================
hdr("1. Column alignment")

if (!file.exists(META_F)) {
  stop("Metadata file not found at ", META_F)
}

hl <- readLines(META_F, n = 1L)
d2 <- readLines(META_F, n = 3L)[2:3]
n_head <- length(strsplit(hl, "\t", fixed = TRUE)[[1]])
n_data <- vapply(d2, function(x) length(strsplit(x, "\t", fixed = TRUE)[[1]]),
                 integer(1))

log_line("  header fields: ", n_head, logfile = LOG)
log_line("  data fields (rows 1-2): ", paste(n_data, collapse = ", "),
         logfile = LOG)

row_names_present <- all(n_data == n_head + 1L)
if (row_names_present) {
  log_line("  -> data rows carry one MORE field than the header.", logfile = LOG)
  log_line("     The file was written with row names. The first field is the",
           logfile = LOG)
  log_line("     spot barcode and every named column sits one to its right.",
           logfile = LOG)
} else if (all(n_data == n_head)) {
  log_line("  -> field counts match; no row-name column.", logfile = LOG)
} else {
  log_line("  [FLAG] field counts are inconsistent: header ", n_head,
           ", data ", paste(unique(n_data), collapse = "/"), logfile = LOG)
  log_line("  Ragged rows. Do not proceed on a guess -- inspect the file.",
           logfile = LOG)
  stop("Inconsistent field counts in the metadata file.")
}

M <- utils::read.delim(META_F, sep = "\t", header = TRUE,
                       row.names = if (row_names_present) 1 else NULL,
                       check.names = FALSE, quote = "",
                       stringsAsFactors = FALSE)
log_line("  parsed: ", nrow(M), " rows x ", ncol(M), " columns", logfile = LOG)

# ---- The independent check on the alignment ---------------------------------
# If the columns are shifted, X1..X15 are not topic proportions and will not sum
# to 1. This is a check the parse cannot pass by accident.
if (!all(TOPIC_COLS %in% names(M))) {
  stop("Topic columns X1..X15 not found after parsing; names are: ",
       paste(utils::head(names(M), 20), collapse = ", "), " ...")
}
TH <- as.matrix(M[, TOPIC_COLS])
mode(TH) <- "numeric"
rs <- rowSums(TH, na.rm = TRUE)
log_line("", logfile = LOG)
log_line(sprintf("  theta row sums: min %.4f, median %.4f, max %.4f",
                 min(rs), stats::median(rs), max(rs)), logfile = LOG)
log_line("  any negative topic proportion: ", any(TH < 0, na.rm = TRUE),
         logfile = LOG)
ok_theta <- abs(stats::median(rs) - 1) < 0.01 && !any(TH < 0, na.rm = TRUE)
if (!ok_theta) {
  log_line("  [STOP] X1..X15 do not behave like a composition. Either the",
           logfile = LOG)
  log_line("  columns are misaligned or these are not proportions. Nothing",
           logfile = LOG)
  log_line("  downstream can be trusted until this is understood.",
           logfile = LOG)
  stop("Topic columns failed the composition check; see ", LOG)
}
log_line("  -> columns are correctly aligned (theta sums to 1, no negatives)",
         logfile = LOG)

# The barcode: the row names if present, otherwise whichever column holds
# barcode-shaped strings.
if (row_names_present) {
  M$.barcode <- rownames(M)
} else {
  cand <- names(M)[vapply(M, function(x)
    is.character(x) && mean(grepl("^[ACGT]{14,18}-\\d+$", x)) > 0.9,
    logical(1))]
  if (!length(cand)) stop("No column holds Visium-shaped barcodes.")
  M$.barcode <- M[[cand[1]]]
  log_line("  barcode column: ", cand[1], logfile = LOG)
}
log_line("  barcode example: ", M$.barcode[1], logfile = LOG)
log_line("  barcodes unique across the whole file: ",
         !any(duplicated(M$.barcode)), logfile = LOG)

# ---- The suffix is a section index, not part of the barcode -----------------
# The first run exposed this. Seurat appends an incrementing index when objects
# are merged, so the SAME physical spot is "AACAATGTGCTCCGAG-1" in their first
# section and "AACAATGTGCTCCGAG-2" in their second, while ours are all "-1"
# straight from Space Ranger. Joining on the barcode as written matched only
# A16 -- their first object -- at 53% and every other section at zero, which
# looked like a coverage problem and was in fact a key problem.
#
# The 16-nucleotide core is the real spot identity and is unique within a
# section, so section + core is the key. Both sides are stripped the same way.
core16 <- function(x) sub("-\\d+$", "", x)
M$.core <- core16(M$.barcode)
log_line("", logfile = LOG)
log_line("  barcode suffixes present: ",
         paste(utils::head(sort(unique(sub("^.*-", "-", M$.barcode))), 25),
               collapse = " "), logfile = LOG)
log_line("  -> the suffix is a merge index, not spot identity. Stripping it.",
         logfile = LOG)
log_line("  16-mer core example: ", M$.core[1], logfile = LOG)
log_line("  core unique within each section: ",
         !any(duplicated(paste0(M[[if ("library_id" %in% names(M)) "library_id" else names(M)[1]]], "|", M$.core))),
         logfile = LOG)

# ==============================================================================
# 2. Which column is the section
# ==============================================================================
# Identified by matching against the section names we already have, not by
# assuming a column name.
hdr("2. The section column")

CORPUS_F <- file.path(PROC, "gse298774_corpus.rds")
if (!file.exists(CORPUS_F)) {
  stop("Need data/processed/gse298774_corpus.rds for our section names ",
       "-- run script 26 (or 26b) first.")
}
OURS <- readRDS(CORPUS_F)$meta          # spot, section, barcode
our_sections <- sort(unique(OURS$section))
log_line("  our sections (", length(our_sections), "): ",
         paste(our_sections, collapse = ", "), logfile = LOG)

# Their section names are not spelled exactly as ours: "B38.2" against our
# "B38_2". A dot for an underscore is the kind of difference that silently
# drops a whole patient, so names are normalised on both sides before matching
# rather than compared literally.
norm_sec <- function(x) gsub("[.\\-]", "_", as.character(x))
our_norm <- norm_sec(our_sections)

score_col <- function(x) {
  if (!is.character(x) && !is.factor(x)) return(0)
  v <- unique(norm_sec(x))
  if (length(v) > 40 || length(v) < 2) return(0)
  mean(v %in% our_norm)
}
sc <- vapply(M, score_col, numeric(1))
sc <- sort(sc[sc > 0], decreasing = TRUE)
log_line("", logfile = LOG)
log_line("  columns whose values look like our section names:", logfile = LOG)
for (nm in names(sc)) {
  log_line(sprintf("    %-20s %.0f%% of its values match", nm, 100 * sc[nm]),
           logfile = LOG)
}
if (!length(sc) || sc[1] < 0.5) {
  log_line("", logfile = LOG)
  log_line("  [STOP] no column matches our section names. Their section", logfile = LOG)
  log_line("  identifiers must be spelled differently from ours. Print the",
           logfile = LOG)
  log_line("  unique values of 'sample' and 'library_id' and map them by hand",
           logfile = LOG)
  log_line("  before going further -- a wrong mapping silently mislabels every",
           logfile = LOG)
  log_line("  spot.", logfile = LOG)
  for (nm in intersect(c("sample", "library_id", "orig.ident", "Batch",
                         "slide", "area"), names(M))) {
    v <- sort(unique(as.character(M[[nm]])))
    log_line("    ", nm, " (", length(v), "): ",
             paste(utils::head(v, 30), collapse = ", "), logfile = LOG)
  }
  stop("Section column could not be identified; see ", LOG)
}
SEC_COL <- names(sc)[1]
log_line("", logfile = LOG)
log_line("  using: ", SEC_COL, logfile = LOG)

# Map each of their (normalised) names back to OUR spelling, so everything
# downstream speaks one vocabulary.
back <- stats::setNames(our_sections, our_norm)
their_raw <- as.character(M[[SEC_COL]])
M$.section <- unname(back[norm_sec(their_raw)])

renamed <- unique(data.frame(theirs = their_raw, ours = M$.section))
renamed <- renamed[!is.na(renamed$ours) & renamed$theirs != renamed$ours, ]
if (nrow(renamed)) {
  log_line("  reconciled spellings:", logfile = LOG)
  for (i in seq_len(nrow(renamed))) {
    log_line("    their ", renamed$theirs[i], " -> our ", renamed$ours[i],
             logfile = LOG)
  }
}

extra <- sort(unique(their_raw[is.na(M$.section)]))
if (length(extra)) {
  log_line("", logfile = LOG)
  log_line("  sections in their file with no counterpart of ours: ",
           paste(extra, collapse = ", "), logfile = LOG)
  log_line("  (A2 is expected here: GEO deposited no barcodes or features for",
           logfile = LOG)
  log_line("   it, so script 26 could not load it. Their labels for A2 exist",
           logfile = LOG)
  log_line("   but there is no count matrix to read expression from, so A2",
           logfile = LOG)
  log_line("   cannot contribute whatever we do.)", logfile = LOG)
  log_line("  dropping ", sum(is.na(M$.section)), " spots from those sections",
           logfile = LOG)
  M <- M[!is.na(M$.section), , drop = FALSE]
  TH <- as.matrix(M[, TOPIC_COLS]); mode(TH) <- "numeric"
}
missing <- setdiff(our_sections, unique(M$.section))
if (length(missing)) {
  log_line("  sections of ours ABSENT from their metadata: ",
           paste(missing, collapse = ", "), logfile = LOG)
  log_line("  those sections cannot contribute to a labels-based comparison.",
           logfile = LOG)
}

M$.spot <- paste0(M$.section, "|", M$.core)
log_line("", logfile = LOG)
log_line("  spot keys unique: ", !any(duplicated(M$.spot)), logfile = LOG)
if (any(duplicated(M$.spot))) {
  log_line("  [FLAG] ", sum(duplicated(M$.spot)),
           " duplicated section+core keys -- the key is not unique and the",
           logfile = LOG)
  log_line("  join would multiply rows. Stopping.", logfile = LOG)
  stop("Non-unique spot key; see ", LOG)
}

# Our side gets the same treatment.
OURS$core <- core16(OURS$barcode)
OURS$key  <- paste0(OURS$section, "|", OURS$core)
if (any(duplicated(OURS$key))) {
  log_line("  [FLAG] our own section+core key is not unique (",
           sum(duplicated(OURS$key)), " duplicates)", logfile = LOG)
}

# ==============================================================================
# 3. Coverage: their spots against ours
# ==============================================================================
hdr("3. Coverage")

cov <- OURS |>
  dplyr::group_by(section) |>
  dplyr::summarise(ours_in_tissue = dplyr::n(), .groups = "drop") |>
  dplyr::left_join(
    M |> dplyr::count(.section, name = "theirs") |>
      dplyr::rename(section = .section), by = "section") |>
  dplyr::mutate(theirs = tidyr::replace_na(theirs, 0L),
                matched = NA_integer_)

matched_keys <- intersect(M$.spot, OURS$key)
mt <- table(sub("\\|.*$", "", matched_keys))
cov$matched <- as.integer(mt[match(cov$section, names(mt))])
cov$matched[is.na(cov$matched)] <- 0L
cov$pct_of_ours <- 100 * cov$matched / cov$ours_in_tissue

log_line(sprintf("  %-8s %12s %10s %10s %9s", "section", "ours in-tiss",
                 "theirs", "matched", "% of ours"), logfile = LOG)
for (i in seq_len(nrow(cov))) {
  log_line(sprintf("  %-8s %12d %10d %10d %8.1f%%", cov$section[i],
                   cov$ours_in_tissue[i], cov$theirs[i], cov$matched[i],
                   cov$pct_of_ours[i]), logfile = LOG)
}
log_line("", logfile = LOG)
log_line(sprintf("  TOTAL   %12d %10d %10d %8.1f%%", sum(cov$ours_in_tissue),
                 nrow(M), length(matched_keys),
                 100 * length(matched_keys) / sum(cov$ours_in_tissue)),
         logfile = LOG)
readr::write_csv(cov, file.path(RES_GEO, "gse298774_kreuger_coverage.csv"))

unmatched <- sum(!M$.spot %in% OURS$key)
log_line("", logfile = LOG)
log_line("  their spots with no match in ours: ", unmatched,
         sprintf(" (%.1f%%)", 100 * unmatched / nrow(M)), logfile = LOG)
if (unmatched / nrow(M) > 0.02) {
  log_line("  [FLAG] more than 2% of their spots do not match. Barcode",
           logfile = LOG)
  log_line("  suffixes may differ (-1 against -1_1, say) or the section", logfile = LOG)
  log_line("  mapping is wrong. Examples of theirs that did not match:",
           logfile = LOG)
  ex <- utils::head(M$.spot[!M$.spot %in% OURS$key], 5)
  for (e in ex) log_line("    ", e, logfile = LOG)
  log_line("  Examples of ours:", logfile = LOG)
  for (e in utils::head(OURS$key, 5)) log_line("    ", e, logfile = LOG)
}

# Uniform filtering looks different from missing sections, so it is worth
# saying which this is rather than only giving the total.
# The verdict is computed only over sections they actually cover, and only
# when there are enough of them to say anything. The first run drew a
# "roughly uniform" conclusion from ONE matched section, which was meaningless.
present <- cov$theirs > 0
log_line("", logfile = LOG)
log_line("  sections they cover: ", sum(present), " of ", nrow(cov),
         logfile = LOG)
if (sum(present) >= 3) {
  r <- range(cov$pct_of_ours[present])
  log_line(sprintf("  retention relative to our in-tissue spots: %.1f%% to %.1f%% (median %.1f%%)",
                   r[1], r[2], stats::median(cov$pct_of_ours[present])),
           logfile = LOG)
  if (r[2] - r[1] < 15) {
    log_line("  -> roughly uniform, consistent with per-spot QC filtering",
             logfile = LOG)
  } else {
    log_line("  -> UNEVEN across sections. This is not simple uniform QC.",
             logfile = LOG)
    log_line("     Sections retaining least:", logfile = LOG)
    o <- order(cov$pct_of_ours)
    for (i in utils::head(o[cov$pct_of_ours[o] > 0 | cov$theirs[o] > 0], 5)) {
      log_line(sprintf("       %-8s %5.1f%% retained (%d of %d)",
                       cov$section[i], cov$pct_of_ours[i], cov$matched[i],
                       cov$ours_in_tissue[i]), logfile = LOG)
    }
    log_line("     Script 24 found a 170-fold depth range across sections, so",
             logfile = LOG)
    log_line("     a depth-based filter would hit shallow sections hardest.",
             logfile = LOG)
    log_line("     This costs spots, and so patients, but it does NOT bias the",
             logfile = LOG)
    log_line("     paired comparison: both compartments of a patient are",
             logfile = LOG)
    log_line("     filtered by the same rule in the same section.", logfile = LOG)
  }
}

# If their spots are a QC-filtered subset, their spots should be DEEPER than the
# ones they dropped. Worth confirming, since it tells us what the filter was.
if ("nCount_Spatial" %in% names(M)) {
  q <- stats::quantile(as.numeric(M$nCount_Spatial), c(0, .25, .5, .75, 1),
                       na.rm = TRUE)
  log_line("", logfile = LOG)
  log_line(sprintf("  their nCount_Spatial: min %.0f, Q1 %.0f, median %.0f, Q3 %.0f, max %.0f",
                   q[1], q[2], q[3], q[4], q[5]), logfile = LOG)
}

# ==============================================================================
# 4. Their labels
# ==============================================================================
hdr("4. Labels")

LBL <- "MaxTopicsName"
if (!LBL %in% names(M)) stop("MaxTopicsName not found after parsing.")

tabl <- M |>
  dplyr::count(.data[[LBL]], name = "spots") |>
  dplyr::arrange(dplyr::desc(spots)) |>
  dplyr::rename(label = 1) |>
  dplyr::mutate(pct = 100 * spots / sum(spots))
log_line(sprintf("  %-24s %8s %7s", "label", "spots", "%"), logfile = LOG)
for (i in seq_len(nrow(tabl))) {
  log_line(sprintf("  %-24s %8d %6.1f%%", tabl$label[i], tabl$spots[i],
                   tabl$pct[i]), logfile = LOG)
}
readr::write_csv(tabl, file.path(RES_GEO, "gse298774_kreuger_label_counts.csv"))

# Which topic number carries which name, from their own file. Recorded because
# our topic numbering will differ and this is the key to comparing the two fits.
if ("MaxTopics" %in% names(M)) {
  map <- M |>
    dplyr::count(MaxTopics, .data[[LBL]], name = "n") |>
    dplyr::rename(label = 2) |>
    dplyr::group_by(MaxTopics) |>
    dplyr::slice_max(n, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::arrange(MaxTopics)
  log_line("", logfile = LOG)
  log_line("  their topic number -> name:", logfile = LOG)
  for (i in seq_len(nrow(map))) {
    log_line(sprintf("    %-5s %s", map$MaxTopics[i], map$label[i]),
             logfile = LOG)
  }
}

# ==============================================================================
# 5. The paired compartments, per patient
# ==============================================================================
# This is the number that decides the readout. Not the spot count.
hdr("5. Nevus and tumour spots per patient")

NEVUS_LAB <- "Nevus"
TUMOUR_LAB <- grep("^Tumor", unique(M[[LBL]]), value = TRUE)
log_line("  nevus label:   ", NEVUS_LAB, logfile = LOG)
log_line("  tumour labels: ", paste(TUMOUR_LAB, collapse = ", "), logfile = LOG)
log_line("  (Tumor 1 and Tumor 2 are pooled as 'tumour': both are melanoma by",
         logfile = LOG)
log_line("   their naming, and splitting them would halve an already small",
         logfile = LOG)
log_line("   per-patient count on a distinction we have no stake in.)",
         logfile = LOG)

# One patient can contribute more than one section (B38 has three in the
# deposit), so patient is the section name stripped of its suffix.
M$.patient <- sub("_.*$", "", sub("ex$", "", M$.section))

pp <- M |>
  dplyr::filter(.data[[LBL]] %in% c(NEVUS_LAB, TUMOUR_LAB)) |>
  dplyr::mutate(grp = ifelse(.data[[LBL]] == NEVUS_LAB, "nevus", "tumour")) |>
  dplyr::count(.patient, grp, name = "n") |>
  tidyr::pivot_wider(names_from = grp, values_from = n, values_fill = 0L) |>
  dplyr::rename(patient = .patient)
if (!"nevus" %in% names(pp)) pp$nevus <- 0L
if (!"tumour" %in% names(pp)) pp$tumour <- 0L
pp <- pp |>
  dplyr::mutate(usable_pair = nevus >= 10 & tumour >= 10) |>
  dplyr::arrange(dplyr::desc(usable_pair), dplyr::desc(pmin(nevus, tumour)))

log_line("", logfile = LOG)
log_line(sprintf("  %-10s %8s %8s %8s", "patient", "nevus", "tumour", "usable"),
         logfile = LOG)
for (i in seq_len(nrow(pp))) {
  log_line(sprintf("  %-10s %8d %8d %8s", pp$patient[i], pp$nevus[i],
                   pp$tumour[i], ifelse(pp$usable_pair[i], "yes", "no")),
           logfile = LOG)
}
n_use <- sum(pp$usable_pair)
log_line("", logfile = LOG)
log_line("  PATIENTS WITH >=10 SPOTS IN BOTH COMPARTMENTS: ", n_use, " of ",
         nrow(pp), logfile = LOG)
log_line("  That is the sample size of the paired test. Not the spot count.",
         logfile = LOG)
if (n_use < 8) {
  log_line("", logfile = LOG)
  log_line("  [FLAG] fewer than 8 usable patients. A paired test on this many",
           logfile = LOG)
  log_line("  can only detect a large and consistent effect; a null result",
           logfile = LOG)
  log_line("  would be uninformative rather than negative, and must be",
           logfile = LOG)
  log_line("  reported as such.", logfile = LOG)
}

# ==============================================================================
# 6. Save the slim join table
# ==============================================================================
hdr("6. Save")

slim <- tibble::tibble(
  spot      = M$.spot,
  section   = M$.section,
  patient   = M$.patient,
  barcode   = M$.barcode,
  label     = M[[LBL]],
  topic     = if ("MaxTopics" %in% names(M)) as.character(M$MaxTopics) else NA,
  nCount    = suppressWarnings(as.numeric(M$nCount_Spatial)),
  nFeature  = suppressWarnings(as.numeric(M$nFeature_Spatial)))
slim <- dplyr::bind_cols(slim, tibble::as_tibble(TH))
slim$purity <- apply(TH, 1, max)

log_line(sprintf("  purity (max topic proportion): median %.3f, Q1 %.3f, Q3 %.3f",
                 stats::median(slim$purity),
                 stats::quantile(slim$purity, .25),
                 stats::quantile(slim$purity, .75)), logfile = LOG)
log_line("  spots with purity > 0.50: ", sum(slim$purity > 0.50), logfile = LOG)
log_line("  spots with purity >= 0.75: ", sum(slim$purity >= 0.75),
         logfile = LOG)
log_line("  (purity thresholds now come from THEIR theta, so a 'pure nevus'",
         logfile = LOG)
log_line("   spot is pure by their fit, not by our reproduction of it)",
         logfile = LOG)

saveRDS(list(labels = slim, coverage = cov, per_patient = pp,
             label_counts = tabl,
             source = list(file = META_F,
                           section_col = SEC_COL,
                           row_names = row_names_present,
                           received = "Slieker, 14 Oct 2026, by email")),
        file.path(PROC, "gse298774_kreuger_labels.rds"))
log_line("", logfile = LOG)
log_line("  saved: data/processed/gse298774_kreuger_labels.rds", logfile = LOG)

# ==============================================================================
# 7. What changes downstream
# ==============================================================================
hdr("7. Next")

log_line("  Script 27 is now a VALIDATION, not the source of labels: our",
         logfile = LOG)
log_line("  independent K=15 fit can be cross-tabulated against their", logfile = LOG)
log_line("  MaxTopicsName. Agreement is evidence the reproduction was", logfile = LOG)
log_line("  faithful; disagreement is worth knowing and would have been",
         logfile = LOG)
log_line("  invisible without this file.", logfile = LOG)
log_line("", logfile = LOG)
log_line("  Script 28 uses THEIR labels. The readout becomes 'RASSF9 in the",
         logfile = LOG)
log_line("  nevus and melanoma compartments as defined by the authors of the",
         logfile = LOG)
log_line("  dataset', which is a materially stronger claim than one resting",
         logfile = LOG)
log_line("  on our own classifier, and it retires the limitation that the",
         logfile = LOG)
log_line("  labels were our reproduction.", logfile = LOG)

report_warnings(logfile = LOG)
write_session_info("29_kreuger_metadata")
log_line("\nFull log: ", LOG, logfile = NULL)
