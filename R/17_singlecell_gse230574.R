# ==============================================================================
# 17_singlecell_gse230574.R
#
# GSE230574: two BD Rhapsody single-cell runs (GSM7226899 = 19145X1,
# GSM7226900 = 19145X2). Each run has a molecule-count matrix
# (*_RSEC_MolsPerCell.csv.gz, cells x genes, HGNC symbols) and a sample-tag
# call file (*_Sample_Tag_Calls.csv.gz) that says which multiplexed sample each
# cell came from. The tags are what make this dataset usable: without them the
# pooled run is an unlabelled mixture of melanoma lines and melanocytes.
#
# WHAT THIS CAN AND CANNOT DO.
#
# CAN: say whether RASSF9 / RASSF10 are DETECTED in melanocytes versus melanoma
# cells, and in what fraction of cells. For a gene that sat at the detection
# floor in bulk cell lines (GSE109245), single-cell detection fraction is the
# more informative readout -- bulk averages a gene expressed in 5% of cells
# down into nothing.
#
# CANNOT: provide a p value for melanocyte-vs-melanoma. There are TWO runs and
# a handful of cultures. Cells from one culture are not independent biological
# replicates of "melanoma"; treating thousands of cells as n = thousands is
# pseudo-replication, and it manufactures p < 1e-50 out of a sample size of
# one. The script therefore reports per-tag descriptives and runs NO group
# test. That is a design limit of the dataset, not a gap to be patched.
#
# MEMORY. These matrices are cells x ~25,000 genes as plain CSV. Reading one
# whole is not safe on a machine with little free RAM, so the file is read in
# chunks: each chunk contributes its row sums (the per-cell library size, needed
# for normalisation) and its two gene columns, and is then discarded.
#
# Output: results/geo/table_geo_rassf9_singlecell_gse230574.csv  (and rassf10)
#         results/geo/gse230574_tag_summary.csv
#         figures/geo/geo_rassf9_singlecell_gse230574.{pdf,png} (and rassf10)
#         results/logs/17_singlecell.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(patchwork))

LOG <- file.path(PATHS$logs, "17_singlecell.log")
cat("GSE230574 single cell —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

GEO_DIR <- file.path(PATHS$raw, "geo")
RAW_DIR <- file.path(GEO_DIR, "GSE230574")
EXT_DIR <- file.path(RAW_DIR, "extracted")
RES_GEO <- file.path(PATHS$results, "geo")
FIG_GEO <- file.path(PATHS$figures, "geo")
invisible(lapply(c(EXT_DIR, RES_GEO, FIG_GEO), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

GENES <- c("RASSF9", "RASSF10")

# POSITIVE CONTROLS. A detection fraction of 0.0% has two possible causes: the
# gene is off, or the column lookup is wrong and an empty column is being read.
# Nothing in the output distinguishes them. So genes that MUST be present are
# quantified alongside the two of interest: two housekeepers, and four
# melanocytic-lineage markers that must be on in melanocytes, nevus cells and
# most melanoma lines. If ACTB is not near 100% of cells in every tag, the
# negative result is a bug and must not be reported.
CONTROLS  <- c("ACTB", "GAPDH", "MITF", "PMEL", "SOX10", "TYR")
ALL_GENES <- c(GENES, CONTROLS)

# A gene is called DETECTED in a tag group only if it clears BOTH bars: present
# in at least this fraction of cells, and at least this mean normalised level.
# One stray molecule in one cell out of 4,000 is ambient RNA, not expression.
MIN_FRAC   <- 0.01     # 1% of cells
MIN_CPM    <- 1.0      # mean counts per 10k molecules

# DEPTH FLOOR for the sensitivity pass. The two runs differ ~37-fold in median
# library size, and at a few hundred molecules per cell most genes read zero
# whether or not they are expressed. A zero from a shallow cell is uninformative
# and must not be pooled with a zero from a deep one, so detection is recomputed
# among cells at or above this depth and the two answers are compared.
DEEP_MIN   <- 2000     # molecules per cell
CHUNK_SIZE <- 25000    # cells read per chunk (64 GB machine)

# ==============================================================================
# 1. Unpack the archive
# ==============================================================================
hdr("1. Unpack")

tarf <- list.files(RAW_DIR, pattern = "\\.tar$", full.names = TRUE)
if (!length(tarf)) {
  stop("No GSE230574_RAW.tar found in ", RAW_DIR,
       "\n  Download it from GEO and place it there, then re-run.")
}
log_line("  archive: ", basename(tarf[1]), logfile = LOG)

already <- list.files(EXT_DIR, pattern = "\\.csv\\.gz$")
if (length(already) < 4) {
  utils::untar(tarf[1], exdir = EXT_DIR)
  log_line("  extracted into data/raw/geo/GSE230574/extracted/", logfile = LOG)
} else {
  log_line("  already extracted (", length(already), " files) -- reusing",
           logfile = LOG)
}

mols_files <- sort(list.files(EXT_DIR, pattern = "MolsPerCell\\.csv\\.gz$",
                              full.names = TRUE))
tag_files  <- sort(list.files(EXT_DIR, pattern = "Sample_Tag_Calls\\.csv\\.gz$",
                              full.names = TRUE))

log_line("  count matrices: ", length(mols_files), logfile = LOG)
log_line("  tag call files: ", length(tag_files), logfile = LOG)
if (!length(mols_files)) stop("No *_RSEC_MolsPerCell.csv.gz after extraction.")

# Run id = the GSM + array prefix shared by a matrix and its tag file.
run_id <- function(f) sub("_(RSEC_MolsPerCell|Sample_Tag_Calls)\\.csv\\.gz$", "",
                          basename(f))

# ==============================================================================
# 2. Series metadata -- what the sample tags are supposed to mean
# ==============================================================================
# The tag names inside the file are usually generic ("SampleTag01_hs"). The
# mapping from tag to culture lives in the GEO series description, so it is
# logged here for Samson to read rather than guessed at in code. Guessing which
# tag is the melanocyte would be the single worst error available in this
# script, so nothing downstream depends on a guess.
hdr("2. Series metadata (tag meanings)")

meta_ok <- requireNamespace("GEOquery", quietly = TRUE)
if (meta_ok) {
  gse <- try(GEOquery::getGEO("GSE230574", GSEMatrix = FALSE,
                              destdir = RAW_DIR), silent = TRUE)
  if (!inherits(gse, "try-error")) {
    for (nm in names(GEOquery::GSMList(gse))) {
      g <- GEOquery::GSMList(gse)[[nm]]
      log_line("  ", nm, ": ", GEOquery::Meta(g)$title, logfile = LOG)
      for (ch in GEOquery::Meta(g)$characteristics_ch1) {
        log_line("      ", ch, logfile = LOG)
      }
    }
    ov <- GEOquery::Meta(gse)$overall_design
    if (!is.null(ov)) {
      log_line("  overall design:", logfile = LOG)
      for (ln in strwrap(ov, width = 68)) log_line("      ", ln, logfile = LOG)
    }
  } else {
    log_line("  [FLAG] GEO metadata fetch failed -- read the series page by hand",
             logfile = LOG)
  }
} else {
  log_line("  GEOquery not installed; skipping metadata", logfile = LOG)
}

# ==============================================================================
# 3. Sample tag calls
# ==============================================================================
hdr("3. Sample tag calls")

read_tags <- function(f) {
  # BD files carry '#'-prefixed experiment metadata above the real header.
  raw <- readr::read_lines(f)
  start <- which(!startsWith(raw, "#"))[1]
  tb <- readr::read_csv(I(raw[start:length(raw)]), show_col_types = FALSE,
                        progress = FALSE)
  names(tb) <- make.names(names(tb))
  cell_col <- grep("^Cell_Index$|^Cell.Index$", names(tb), value = TRUE)[1]
  tag_col  <- grep("Sample_Name|Sample.Name", names(tb), value = TRUE)[1]
  if (is.na(tag_col)) {
    tag_col <- grep("Sample_Tag|Sample.Tag", names(tb), value = TRUE)[1]
  }
  tibble::tibble(run  = run_id(f),
                 cell = as.character(tb[[cell_col]]),
                 tag  = as.character(tb[[tag_col]]))
}

tags <- purrr::map_dfr(tag_files, read_tags)

tag_counts <- tags |>
  dplyr::count(run, tag, name = "n_cells") |>
  dplyr::arrange(run, dplyr::desc(n_cells))

log_line("  cells with a tag call: ", nrow(tags), logfile = LOG)
log_line("", logfile = LOG)
for (i in seq_len(nrow(tag_counts))) {
  log_line(sprintf("    %-28s %-22s %6d cells",
                   tag_counts$run[i], tag_counts$tag[i], tag_counts$n_cells[i]),
           logfile = LOG)
}

# Multiplet and Undetermined are BD's own quality calls, not samples. They are
# counted and then dropped: a multiplet is two cells in one well, so its counts
# belong to no single culture.
DROP_TAGS <- c("Multiplet", "Undetermined")
n_drop <- sum(tags$tag %in% DROP_TAGS)
log_line("", logfile = LOG)
log_line("  dropping ", n_drop, " cells called Multiplet/Undetermined (",
         sprintf("%.1f%%", 100 * n_drop / nrow(tags)), " of tagged cells)",
         logfile = LOG)
tags_use <- tags |> dplyr::filter(!tag %in% DROP_TAGS)
log_line("  cells carried forward: ", nrow(tags_use), logfile = LOG)

# ==============================================================================
# 4. Quantify the two genes, chunk by chunk
# ==============================================================================
hdr("4. Expression extraction")

extract_genes <- function(f) {
  rid <- run_id(f)
  log_line("  ", rid, logfile = LOG)

  # -- header only, to locate the gene columns before reading any data --------
  con <- gzfile(f, "r"); on.exit(try(close(con), silent = TRUE), add = TRUE)
  repeat {
    ln <- readLines(con, n = 1)
    if (!length(ln)) stop("Reached end of file without a header: ", basename(f))
    if (!startsWith(ln, "#")) break
  }
  hdr_names <- make.names(strsplit(ln, ",", fixed = TRUE)[[1]])
  close(con)

  cell_idx <- grep("^Cell_Index$|^Cell.Index$", hdr_names)[1]
  if (is.na(cell_idx)) cell_idx <- 1L
  gene_pos <- stats::setNames(match(ALL_GENES, hdr_names), ALL_GENES)

  log_line("    columns: ", length(hdr_names), logfile = LOG)
  log_line("    found:  ", paste(ALL_GENES[!is.na(gene_pos)], collapse = ", "),
           logfile = LOG)
  if (any(is.na(gene_pos))) {
    log_line("    ABSENT: ", paste(ALL_GENES[is.na(gene_pos)], collapse = ", "),
             logfile = LOG)
  }
  # The column index is logged so the lookup can be checked by eye against the
  # header. A silently wrong index is the one failure mode that would turn a
  # real signal into a clean-looking zero.
  for (g in ALL_GENES[!is.na(gene_pos)]) {
    log_line(sprintf("      %-8s -> column %d (%s)", g, gene_pos[[g]],
                     hdr_names[gene_pos[[g]]]), logfile = LOG)
  }

  # -- chunked pass: library size from every column, values from the named -----
  acc <- list()
  cb <- function(chunk, pos) {
    names(chunk) <- hdr_names[seq_len(ncol(chunk))]
    cells <- as.character(chunk[[cell_idx]])
    mat <- as.matrix(chunk[, -cell_idx, drop = FALSE])
    mat[is.na(mat)] <- 0
    out <- tibble::tibble(cell = cells, lib_size = rowSums(mat))
    for (g in ALL_GENES) {
      out[[g]] <- if (is.na(gene_pos[[g]])) NA_real_
                  else as.numeric(chunk[[gene_pos[[g]]]])
    }
    acc[[length(acc) + 1L]] <<- out
    invisible(NULL)
  }

  readr::read_csv_chunked(
    f, callback = readr::SideEffectChunkCallback$new(cb),
    chunk_size = CHUNK_SIZE, comment = "#", col_names = TRUE,
    col_types = readr::cols(.default = readr::col_double()),
    progress = FALSE)

  res <- dplyr::bind_rows(acc) |> dplyr::mutate(run = rid, .before = 1)
  log_line("    cells read: ", nrow(res),
           " | median library size: ", round(median(res$lib_size)),
           logfile = LOG)
  res
}

expr <- purrr::map_dfr(mols_files, extract_genes)

# Counts per 10,000 molecules. Standard single-cell depth normalisation: raw
# molecule counts are not comparable between cells sequenced at different
# depths, and depth varies several-fold within a run.
for (g in ALL_GENES) {
  expr[[paste0(g, "_cpm")]] <- 1e4 * expr[[g]] / pmax(expr$lib_size, 1)
}

d <- dplyr::inner_join(tags_use, expr, by = c("run", "cell"))
log_line("", logfile = LOG)
log_line("  cells with BOTH a tag call and counts: ", nrow(d), logfile = LOG)
if (!nrow(d)) {
  stop("No cells matched between the tag files and the count matrices.\n",
       "  Cell_Index values differ between the two file types -- inspect the\n",
       "  first few rows of each before going further.")
}

# ==============================================================================
# 5a. Sequencing depth per tag
# ==============================================================================
# The run-level median reported in section 4 is dragged down by Undetermined
# cells, most of which are empty droplets. Depth must be judged per tag, among
# the cells actually being compared.
hdr("5a. Sequencing depth per tag")

depth <- d |>
  dplyr::group_by(run, tag) |>
  dplyr::summarise(n_cells = dplyr::n(),
                   median_lib = median(lib_size),
                   q25_lib = quantile(lib_size, 0.25),
                   q75_lib = quantile(lib_size, 0.75),
                   n_deep = sum(lib_size >= DEEP_MIN),
                   frac_deep = n_deep / n_cells,
                   .groups = "drop") |>
  dplyr::arrange(run, dplyr::desc(median_lib))

log_line(sprintf("  %-26s %-20s %6s %10s %8s %7s",
                 "run", "tag", "cells", "med_lib", "n_deep", "%deep"),
         logfile = LOG)
for (i in seq_len(nrow(depth))) {
  log_line(sprintf("  %-26s %-20s %6d %10d %8d %6.1f%%",
                   depth$run[i], depth$tag[i], depth$n_cells[i],
                   round(depth$median_lib[i]), depth$n_deep[i],
                   100 * depth$frac_deep[i]), logfile = LOG)
}
readr::write_csv(depth, file.path(RES_GEO, "gse230574_depth_per_tag.csv"))
log_line("", logfile = LOG)
log_line("  A tag whose median library size is a few hundred molecules cannot",
         logfile = LOG)
log_line("  support a claim that any particular gene is absent: most genes read",
         logfile = LOG)
log_line("  zero at that depth. Such tags are reported but must not be used as",
         logfile = LOG)
log_line("  evidence of absence. Depth floor for the sensitivity pass: ",
         DEEP_MIN, " molecules.", logfile = LOG)

# ==============================================================================
# 5b. Positive controls -- does the pipeline detect genes that must be there?
# ==============================================================================
# This section decides whether the RASSF9/RASSF10 result below can be believed
# at all. ACTB and GAPDH must be near-ubiquitous in every tag. MITF, PMEL,
# SOX10 and TYR must be on in melanocytes and nevus cells, and in most melanoma
# lines (TYR is the one legitimately variable marker -- dedifferentiated lines
# lose it). If the housekeepers fail, STOP: the lookup or the parsing is wrong
# and every zero below is an artefact.
hdr("5b. Positive controls")

ctrl <- list()
for (g in CONTROLS) {
  if (all(is.na(d[[g]]))) {
    log_line("  ", g, ": not in the matrices", logfile = LOG); next
  }
  ctrl[[g]] <- d |>
    dplyr::group_by(run, tag) |>
    dplyr::summarise(n_cells = dplyr::n(),
                     frac_pos = mean(.data[[g]] > 0, na.rm = TRUE),
                     mean_cpm = mean(.data[[paste0(g, "_cpm")]], na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::mutate(gene = g, .before = 1)
}
ctrl_all <- dplyr::bind_rows(ctrl)

if (nrow(ctrl_all)) {
  wide <- ctrl_all |>
    dplyr::mutate(cell_txt = sprintf("%5.1f%%", 100 * frac_pos)) |>
    dplyr::select(run, tag, gene, cell_txt) |>
    tidyr::pivot_wider(names_from = gene, values_from = cell_txt)

  gcols <- setdiff(names(wide), c("run", "tag"))
  log_line("  percent of cells positive:", logfile = LOG)
  log_line(sprintf("  %-26s %-20s %s", "run", "tag",
                   paste(sprintf("%8s", gcols), collapse = " ")), logfile = LOG)
  for (i in seq_len(nrow(wide))) {
    log_line(sprintf("  %-26s %-20s %s", wide$run[i], wide$tag[i],
                     paste(sprintf("%8s", unlist(wide[i, gcols])),
                           collapse = " ")), logfile = LOG)
  }
  readr::write_csv(ctrl_all, file.path(RES_GEO, "gse230574_positive_controls.csv"))

  # The automated verdict. ACTB in deep tags is the single decisive check.
  deep_tags <- depth$tag[depth$median_lib >= DEEP_MIN]
  actb <- ctrl_all |> dplyr::filter(gene == "ACTB", tag %in% deep_tags)
  log_line("", logfile = LOG)
  if (!nrow(actb)) {
    log_line("  [FLAG] no tag reaches the depth floor -- controls cannot be",
             logfile = LOG)
    log_line("         evaluated, and NO absence claim is supportable.",
             logfile = LOG)
  } else if (min(actb$frac_pos) >= 0.80) {
    log_line("  CONTROLS PASS: ACTB is positive in >=80% of cells in every",
             logfile = LOG)
    log_line("  adequately sequenced tag (minimum ",
             sprintf("%.1f%%", 100 * min(actb$frac_pos)),
             "). The pipeline detects genes", logfile = LOG)
    log_line("  that are present, so a zero below is evidence of absence.",
             logfile = LOG)
  } else {
    log_line("  [FLAG] CONTROLS FAIL: ACTB reaches only ",
             sprintf("%.1f%%", 100 * min(actb$frac_pos)),
             " of cells in the", logfile = LOG)
    log_line("         worst adequately sequenced tag. Something is wrong with",
             logfile = LOG)
    log_line("         the column lookup or the parsing. DO NOT report the",
             logfile = LOG)
    log_line("         RASSF9/RASSF10 result until this is resolved.",
             logfile = LOG)
  }
} else {
  log_line("  [FLAG] no control genes found in the matrices at all -- the gene",
           logfile = LOG)
  log_line("         symbols in this file are not HGNC. Inspect the header.",
           logfile = LOG)
}

# ==============================================================================
# 5c. The two genes of interest
# ==============================================================================
hdr("5c. Per-tag detection and level")

summ <- list()
for (g in GENES) {
  if (all(is.na(d[[g]]))) {
    log_line("  ", g, ": ABSENT from the count matrices -- not quantifiable",
             logfile = LOG)
    next
  }
  s <- d |>
    dplyr::group_by(run, tag) |>
    dplyr::summarise(
      n_cells      = dplyr::n(),
      n_pos        = sum(.data[[g]] > 0, na.rm = TRUE),
      frac_pos     = n_pos / n_cells,
      mean_cpm     = mean(.data[[paste0(g, "_cpm")]], na.rm = TRUE),
      median_cpm_pos = ifelse(n_pos > 0,
                              median(.data[[paste0(g, "_cpm")]][.data[[g]] > 0]),
                              NA_real_),
      max_counts   = max(.data[[g]], na.rm = TRUE),
      median_lib   = median(lib_size),
      # Sensitivity pass: the same two measures among adequately sequenced
      # cells only. A gene can look absent in a shallow tag and present once
      # the shallow cells are removed; if the deep-cell figure matches the
      # all-cell figure, depth is not driving the result.
      n_deep       = sum(lib_size >= DEEP_MIN),
      n_pos_deep   = sum(lib_size >= DEEP_MIN & .data[[g]] > 0, na.rm = TRUE),
      frac_pos_deep = ifelse(n_deep > 0, n_pos_deep / n_deep, NA_real_),
      mean_cpm_deep = ifelse(n_deep > 0,
                             mean(.data[[paste0(g, "_cpm")]][lib_size >= DEEP_MIN],
                                  na.rm = TRUE), NA_real_),
      .groups = "drop") |>
    dplyr::mutate(
      gene = g,
      detected = frac_pos >= MIN_FRAC & mean_cpm >= MIN_CPM,
      # A tag too shallow to detect anything is marked as such rather than
      # being counted as evidence that the gene is off.
      depth_ok = median_lib >= DEEP_MIN,
      detected_deep = !is.na(frac_pos_deep) &
                      frac_pos_deep >= MIN_FRAC & mean_cpm_deep >= MIN_CPM,
      .before = 1) |>
    dplyr::arrange(run, dplyr::desc(frac_pos))

  summ[[g]] <- s

  log_line("", logfile = LOG)
  log_line("  ", g, logfile = LOG)
  log_line(sprintf("    %-26s %-20s %6s %8s %9s %5s %7s %9s %5s",
                   "run", "tag", "cells", "%pos", "mean_cpm", "det",
                   "%pos_dp", "cpm_deep", "dp_ok"), logfile = LOG)
  for (i in seq_len(nrow(s))) {
    log_line(sprintf("    %-26s %-20s %6d %7.1f%% %9.3f %5s %6.1f%% %9.3f %5s",
                     s$run[i], s$tag[i], s$n_cells[i], 100 * s$frac_pos[i],
                     s$mean_cpm[i], ifelse(s$detected[i], "yes", "no"),
                     100 * ifelse(is.na(s$frac_pos_deep[i]), 0, s$frac_pos_deep[i]),
                     ifelse(is.na(s$mean_cpm_deep[i]), 0, s$mean_cpm_deep[i]),
                     ifelse(s$depth_ok[i], "yes", "NO")), logfile = LOG)
  }

  n_ok <- sum(s$depth_ok)
  log_line("", logfile = LOG)
  log_line("    tags with adequate depth: ", n_ok, " of ", nrow(s), logfile = LOG)
  if (!any(s$detected) && !any(s$detected_deep)) {
    if (n_ok == 0) {
      log_line("    -> NOT DETECTED, but NO tag has adequate depth. This is",
               logfile = LOG)
      log_line("       uninformative: absence cannot be distinguished from",
               logfile = LOG)
      log_line("       undersequencing. Report it as such.", logfile = LOG)
    } else {
      log_line("    -> NOT DETECTED in any tag, including the ", n_ok,
               " adequately", logfile = LOG)
      log_line("       sequenced one(s) (floor: >=", 100 * MIN_FRAC,
               "% of cells and >= ", MIN_CPM, " mean CPM).", logfile = LOG)
      log_line("       In those tags this IS evidence of absence, provided the",
               logfile = LOG)
      log_line("       controls in 5b passed.", logfile = LOG)
    }
  } else if (any(s$detected_deep) && !any(s$detected)) {
    log_line("    -> DETECTED only after the depth filter. The all-cell result",
             logfile = LOG)
    log_line("       was a depth artefact; use the deep-cell columns.",
             logfile = LOG)
  }
}
summ_all <- dplyr::bind_rows(summ)

if (nrow(summ_all)) {
  readr::write_csv(summ_all, file.path(RES_GEO, "gse230574_tag_summary.csv"))
  log_line("", logfile = LOG)
  log_line("  wrote: results/geo/gse230574_tag_summary.csv", logfile = LOG)
}

log_line("", logfile = LOG)
log_line("  NO group test is run. Cells from one culture are not independent", logfile = LOG)
log_line("  replicates of their lesion type, and two runs cannot support a", logfile = LOG)
log_line("  melanocyte-vs-melanoma p value. Detection fraction per tag is the", logfile = LOG)
log_line("  honest readout, and its value is that a gene expressed in a small", logfile = LOG)
log_line("  subset of cells is invisible to bulk RNA-seq -- which is exactly", logfile = LOG)
log_line("  the scenario GSE109245 could not rule out.", logfile = LOG)

# ==============================================================================
# 6. Figures, one per gene
# ==============================================================================
hdr("6. Figures")

# Fifteen tag names will not fit under a vertical bar chart at any font size
# that stays legible, so the bars run HORIZONTALLY: the names sit on the y axis
# and read left to right at full size, with no rotation and no overlap. Rows
# are grouped into panels by biological class, which is also what colour
# carries -- colour showing "detected / not detected" would be wasted here,
# since nothing is detected and the bar length already says so.
make_fig <- function(g) {
  s <- summ[[g]]
  if (is.null(s) || !nrow(s)) return(NULL)

  s <- s |> dplyr::mutate(
    class = dplyr::case_when(
      grepl("^SampleTag", tag)  ~ "Nevus (primary)",
      grepl("^melanocyte", tag) ~ "Melanocyte (cultured)",
      TRUE                      ~ "Melanoma line"),
    # A dagger marks a tag too shallow to support a claim of absence.
    lab = paste0(tag, ifelse(depth_ok, "", " †")))

  s$class <- factor(s$class, levels = names(PAL_CLASS))
  s <- s |> dplyr::arrange(class, frac_pos)
  s$lab <- factor(s$lab, levels = s$lab)

  any_det <- any(s$detected & s$depth_ok)

  # Long format so the two measures become facet COLUMNS. Each column gets its
  # own x scale and its own floor line; each row block is sized to the number
  # of tags it holds, so no panel is padded with empty space.
  long <- dplyr::bind_rows(
    s |> dplyr::transmute(class, lab, depth_ok,
                          metric = "Cells positive (%)",
                          value = 100 * frac_pos),
    s |> dplyr::transmute(class, lab, depth_ok,
                          metric = "Mean counts per 10,000 molecules",
                          value = mean_cpm))
  # Column headers kept SHORT. A long facet label is what pushes a panel past
  # the page edge, and the units belong in the legend text anyway.
  long$metric <- factor(long$metric,
                        levels = c("Cells positive (%)",
                                   "Mean counts per 10,000 molecules"),
                        labels = c("Cells positive (%)",
                                   "Mean molecules per 10,000"))

  floors <- tibble::tibble(
    metric = factor(levels(long$metric), levels = levels(long$metric)),
    at = c(100 * MIN_FRAC, MIN_CPM))

  ggplot(long, aes(y = lab, x = value, fill = class)) +
    geom_vline(data = floors, aes(xintercept = at), linetype = "22",
               colour = PAL$grey_pt, linewidth = 0.4) +
    geom_col(width = 0.68, colour = PAL$ink, linewidth = 0.25) +
    # Shallow tags get a dotted outline so they read as provisional even in
    # greyscale, on top of the dagger on the label.
    geom_col(data = dplyr::filter(long, !depth_ok), width = 0.68,
             fill = NA, colour = PAL$ink, linewidth = 0.25, linetype = "11") +
    facet_grid(class ~ metric, scales = "free", space = "free_y",
               switch = "y") +
    scale_fill_manual(values = PAL_CLASS, name = NULL, drop = FALSE) +
    # Generous right-hand expansion: the panel can then never be the thing that
    # clips, whatever the longest bar turns out to be.
    scale_x_continuous(expand = expansion(mult = c(0, 0.22)),
                       n.breaks = 4) +
    coord_cartesian(clip = "off") +
    labs(y = NULL, x = NULL) +
    theme_journal(grid = FALSE) +
    theme(legend.position = "top",
          legend.justification = "left",
          panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey92"),
          panel.spacing.y = unit(5, "pt"),
          panel.spacing.x = unit(11, "pt"),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 0,
                                           size = FIG$base_size - 1),
          strip.text.x = element_text(size = FIG$base_size - 0.5),
          axis.text.y = element_text(size = FIG$base_size - 1, hjust = 1),
          axis.text.x = element_text(size = FIG$base_size - 1.5),
          plot.margin = margin(4, 10, 4, 4, unit = "pt"))
}

reset_legends("# GEO figure legends\n\nPaste into the manuscript. The figures\nthemselves carry no titles or captions, per journal convention.\n")

for (g in GENES) {
  fig <- make_fig(g)
  if (is.null(fig)) { log_line("  ", g, ": skipped", logfile = LOG); next }

  stem <- paste0("geo_", tolower(g), "_singlecell_gse230574")
  save_fig_to(fig, stem, FIG_GEO, width = FIG$width_double, height = 150)

  s <- summ[[g]]
  add_legend(
    stem,
    paste0("Figure. ", g, " across multiplexed samples in GSE230574."),
    paste0(
      "Single-cell molecule counts from two BD Rhapsody runs, demultiplexed ",
      "by sample tag and grouped by biological class (rows). Left column, the ",
      "percentage of cells in which ", g, " was detected; right column, the ",
      "mean number of ", g, " molecules per 10,000 total molecules. Dashed ",
      "lines mark the detection floor used throughout (",
      100 * MIN_FRAC, "% of cells and ", MIN_CPM,
      " molecules per 10,000); a gene must clear both to be called detected. ",
      "Daggers mark samples whose median library size is below ", DEEP_MIN,
      " molecules per cell, in which absence cannot be distinguished from ",
      "undersequencing, and these bars are outlined rather than filled. ",
      "Multiplet and Undetermined cells were excluded (",
      sprintf("%.1f%%", 100 * n_drop / nrow(tags)), " of tagged cells). ",
      "n = ", sum(s$n_cells), " cells across ", nrow(s), " samples. ",
      "Descriptive only: cells within one culture are not independent ",
      "biological replicates, so no group comparison is performed."))

  log_line("  wrote: ", stem, ".{pdf,png,svg}", logfile = LOG)

  readr::write_csv(
    d |> dplyr::select(run, tag, cell, lib_size,
                       dplyr::all_of(c(g, paste0(g, "_cpm")))),
    file.path(RES_GEO, paste0("table_geo_", tolower(g),
                              "_singlecell_gse230574.csv")))
}

report_warnings(logfile = LOG)
write_session_info("17_singlecell")
log_line("\nFull log: ", LOG, logfile = NULL)
