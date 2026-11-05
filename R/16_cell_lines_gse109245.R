# ==============================================================================
# 16_cell_lines_gse109245.R
#
# GSE109245: 5 metastatic melanoma cell lines (MRA2, MRA4, MRA5, MRA6, MRA9)
# versus 1 primary melanocyte culture. A cell-intrinsic comparison: no stroma,
# no immune cells, so any difference is melanocytic-lineage biology rather than
# tissue composition. That is its value, and it is the one thing the bulk tumour
# datasets cannot tell us.
#
# WHAT IT CANNOT DO. There is ONE melanocyte sample. With n = 1 in a group no
# statistical test is possible or meaningful, and none is run. The melanocyte
# value is a single reference point, not a group mean, and the figure shows it
# as such. Five melanoma lines give a spread; the melanocyte gives a line.
#
# TRANSCRIPT-LEVEL INPUT. This series deposits salmon output at TRANSCRIPT level
# (ENST ids), log2-transformed. Gene-level values therefore require: reverse the
# log, sum transcripts belonging to the gene, re-log. Summing log values
# directly would be wrong (it multiplies rather than adds).
#
# Transcript-to-gene mapping uses EnsDb.Hsapiens.v86 (offline, hg38). The salmon
# index release is unknown, so the script reports how many of the gene's
# transcripts were actually found -- a low match rate is a real caveat, not a
# detail to skip past.
#
# Output: results/geo/table_geo_rassf9_celllines_gse109245.csv  (and rassf10)
#         figures/geo/geo_rassf9_celllines_gse109245.{pdf,png} (and rassf10)
#         results/logs/16_celllines.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(patchwork))

LOG <- file.path(PATHS$logs, "16_celllines.log")
cat("GSE109245 cell lines —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

GEO_DIR <- file.path(PATHS$raw, "geo")
RES_GEO <- file.path(PATHS$results, "geo")
FIG_GEO <- file.path(PATHS$figures, "geo")
dir.create(RES_GEO, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_GEO, recursive = TRUE, showWarnings = FALSE)

ENS   <- c(RASSF9 = "ENSG00000198774", RASSF10 = "ENSG00000189431")
GENES <- names(ENS)

# ==============================================================================
# 1. Load the salmon matrix
# ==============================================================================
hdr("1. Input")

f <- list.files(file.path(GEO_DIR, "GSE109245"), pattern = "salmon.*csv",
                full.names = TRUE)
if (!length(f)) stop("GSE109245 salmon file not found in ",
                     file.path(GEO_DIR, "GSE109245"))

x <- as.data.frame(readr::read_csv(f[1], show_col_types = FALSE))
rownames(x) <- x[[1]]
m <- as.matrix(x[, -1, drop = FALSE])

log_line("  file: ", basename(f[1]), logfile = LOG)
log_line("  transcripts: ", nrow(m), " | samples: ", ncol(m), logfile = LOG)
log_line("  columns: ", paste(colnames(m), collapse = ", "), logfile = LOG)
log_line("  value range: ", round(min(m, na.rm = TRUE), 2), " to ",
         round(max(m, na.rm = TRUE), 2), " (log2 scale, per the file name)",
         logfile = LOG)

# ==============================================================================
# 2. Transcript -> gene
# ==============================================================================
hdr("2. Transcript-to-gene mapping")

# Two routes to a transcript-to-gene map. The local annotation package is
# preferred (offline, fixed release). If it is unavailable or its database is
# damaged, fall back to a live Ensembl BioMart query, which needs the network
# but returns the same mapping. Whichever is used is recorded, because the
# Ensembl release affects which transcript ids exist.
tx <- NULL
map_source <- NA_character_

if (requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
  attempt <- try({
    suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
    ensembldb::transcripts(EnsDb.Hsapiens.v86,
                           columns = c("tx_id", "gene_id"),
                           return.type = "data.frame")
  }, silent = TRUE)
  if (!inherits(attempt, "try-error")) {
    tx <- attempt
    map_source <- "EnsDb.Hsapiens.v86 (local, Ensembl 86)"
  } else {
    log_line("  [note] EnsDb.Hsapiens.v86 present but failed to load:",
             logfile = LOG)
    log_line("         ", trimws(as.character(attempt)), logfile = LOG)
  }
}

if (is.null(tx) && requireNamespace("biomaRt", quietly = TRUE)) {
  log_line("  falling back to Ensembl BioMart (needs internet)...",
           logfile = LOG)
  attempt <- try({
    mart <- biomaRt::useEnsembl(biomart = "genes",
                                dataset = "hsapiens_gene_ensembl")
    bm <- biomaRt::getBM(
      attributes = c("ensembl_transcript_id", "ensembl_gene_id"),
      filters = "ensembl_gene_id", values = unname(ENS), mart = mart)
    data.frame(tx_id = bm$ensembl_transcript_id,
               gene_id = bm$ensembl_gene_id, stringsAsFactors = FALSE)
  }, silent = TRUE)
  if (!inherits(attempt, "try-error") && nrow(attempt)) {
    tx <- attempt
    map_source <- "Ensembl BioMart (live, current release)"
  } else {
    log_line("  [note] BioMart query failed", logfile = LOG)
  }
}

if (is.null(tx)) {
  stop("No transcript-to-gene mapping available.\n",
       "Either repair the local annotation:\n",
       "  Rscript -e \"BiocManager::install('EnsDb.Hsapiens.v86', ask=FALSE, force=TRUE)\"\n",
       "or install biomaRt for the online route:\n",
       "  Rscript -e \"BiocManager::install('biomaRt', ask=FALSE)\"")
}

log_line("  annotation source: ", map_source, logfile = LOG)
log_line("  transcripts in map: ", nrow(tx), logfile = LOG)

row_tx <- sub("\\..*$", "", rownames(m))   # drop any version suffix

gene_values <- list()
for (g in GENES) {
  want <- tx$tx_id[tx$gene_id == ENS[[g]]]
  hit  <- which(row_tx %in% want)

  log_line(sprintf("  %-8s %d transcripts annotated, %d present in this matrix",
                   g, length(want), length(hit)), logfile = LOG)

  if (!length(hit)) {
    log_line("    [FLAG] no transcripts found -- gene cannot be quantified here",
             logfile = LOG)
    next
  }
  if (length(hit) < length(want) * 0.5) {
    log_line("    [FLAG] fewer than half the annotated transcripts matched;",
             logfile = LOG)
    log_line("           the salmon index release probably differs from v86.",
             logfile = LOG)
    log_line("           Gene-level values are still computed but are a",
             logfile = LOG)
    log_line("           partial sum -- state this if the result is used.",
             logfile = LOG)
  }

  # Reverse the log, sum on the linear scale, re-log. Summing logs would
  # multiply the transcripts together instead of adding them.
  lin <- 2^m[hit, , drop = FALSE] - 1
  lin[lin < 0] <- 0
  gene_values[[g]] <- log2(colSums(lin, na.rm = TRUE) + 1)

  log_line("    gene-level log2 values: ",
           paste(sprintf("%s=%.2f", names(gene_values[[g]]),
                         gene_values[[g]]), collapse = "  "), logfile = LOG)
}

if (!length(gene_values)) stop("Neither gene could be quantified in GSE109245.")

# ==============================================================================
# 3. Assemble
# ==============================================================================
hdr("3. Samples")

# "mcf" is the melanocyte culture; m2/m4/m5/m6/m9 are the metastatic lines.
d <- tibble::tibble(sample = colnames(m)) |>
  dplyr::mutate(
    group = ifelse(grepl("^mcf", sample, ignore.case = TRUE),
                   "Melanocyte", "Melanoma line"),
    label = ifelse(group == "Melanocyte", "Melanocyte",
                   toupper(sub("^m", "MRA", sample)))
  )
for (g in names(gene_values)) d[[g]] <- as.numeric(gene_values[[g]][d$sample])

cat(capture.output(print(as.data.frame(d), row.names = FALSE, digits = 3)),
    sep = "\n", file = LOG, append = TRUE)
print(as.data.frame(d), row.names = FALSE, digits = 3)

log_line("\n  n = 1 melanocyte sample: NO test is performed. The melanocyte is",
         logfile = LOG)
log_line("  a single reference value, and the five lines are shown",
         logfile = LOG)
log_line("  individually rather than summarised.", logfile = LOG)

# ==============================================================================
# 4. Descriptive comparison
# ==============================================================================
hdr("4. Descriptive summary")

# A DETECTION FLOOR must be checked before any direction is claimed. If a gene
# sits at or near zero in every sample, "all lines above the melanocyte" is an
# artefact of nothing being below zero -- not a biological direction. The
# threshold below (max < 1 on the log2 scale, i.e. under one summed TPM) marks
# the gene as undetected and suppresses the directional wording entirely.
DETECT_FLOOR <- 1.0

summ <- list()
for (g in names(gene_values)) {
  mc <- d[[g]][d$group == "Melanocyte"]
  ml <- d[[g]][d$group == "Melanoma line"]
  allv <- c(mc, ml)

  undetected <- max(allv, na.rm = TRUE) < DETECT_FLOOR
  n_below <- sum(ml < mc, na.rm = TRUE)
  n_equal <- sum(ml == mc, na.rm = TRUE)
  n_above <- sum(ml > mc, na.rm = TRUE)

  direction <- if (undetected) {
    "NOT DETECTED in any sample -- no direction can be read"
  } else if (n_below == length(ml)) {
    "all lines below melanocyte"
  } else if (n_above == length(ml)) {
    "all lines above melanocyte"
  } else {
    "mixed"
  }

  summ[[g]] <- tibble::tibble(
    gene = g, melanocyte = mc,
    lines_median = median(ml, na.rm = TRUE),
    lines_min = min(ml, na.rm = TRUE), lines_max = max(ml, na.rm = TRUE),
    n_lines = length(ml),
    n_below_melanocyte = n_below, n_equal_melanocyte = n_equal,
    n_above_melanocyte = n_above,
    detected = !undetected, direction = direction
  )

  log_line(sprintf("  %-8s melanocyte %.2f | lines median %.2f (range %.2f-%.2f)",
                   g, mc, median(ml), min(ml), max(ml)), logfile = LOG)
  log_line(sprintf("           %d below / %d equal / %d above the melanocyte value",
                   n_below, n_equal, n_above), logfile = LOG)
  log_line("           -> ", direction, logfile = LOG)

  if (undetected) {
    log_line(sprintf("           every value is under %.1f on the log2 scale:",
                     DETECT_FLOOR), logfile = LOG)
    log_line("           the gene is at the detection floor in these cultures,",
             logfile = LOG)
    log_line("           so this dataset is UNINFORMATIVE for it. Report that,",
             logfile = LOG)
    log_line("           not a direction.", logfile = LOG)
  }
}
summ <- dplyr::bind_rows(summ)

log_line("\n  Where a gene IS detected, 'all five lines on the same side of the",
         logfile = LOG)
log_line("  single melanocyte' is the strongest statement this design",
         logfile = LOG)
log_line("  supports -- suggestive, not a test, since one donor cannot",
         logfile = LOG)
log_line("  represent melanocytes. Where it is not detected, nothing can be",
         logfile = LOG)
log_line("  said at all.", logfile = LOG)

# ==============================================================================
# 5. Figures, one per gene
# ==============================================================================
hdr("5. Figures")

make_fig <- function(g) {
  if (!g %in% names(d)) return(NULL)
  pd <- d[!is.na(d[[g]]), ]
  if (!nrow(pd)) return(NULL)
  mc <- pd[[g]][pd$group == "Melanocyte"]

  # Carry the detection verdict onto the figure itself. A bar chart of values
  # that all sit within rounding distance of zero looks like data unless the
  # panel says otherwise, and a reader seeing the file out of context would
  # read a direction out of it that is not there.
  is_undetected <- !isTRUE(summ$detected[summ$gene == g])

  pd <- pd |> dplyr::arrange(group != "Melanocyte", .data[[g]])
  pd$label <- factor(pd$label, levels = pd$label)

  sub <- if (is_undetected) {
    paste0("At the detection floor in every sample (all values < ",
           DETECT_FLOOR, " on the log2 scale) — uninformative")
  } else {
    "Dashed line = the melanocyte reference value"
  }

  p <- ggplot(pd, aes(x = label, y = .data[[g]], fill = group)) +
    # The melanocyte reference drawn as a rule across the panel, so every line
    # can be read against it directly.
    geom_hline(yintercept = mc, linetype = "22", colour = PAL$nevus,
               linewidth = 0.5) +
    geom_col(width = 0.62, colour = PAL$ink, linewidth = 0.3) +
    scale_fill_manual(values = c("Melanocyte" = PAL$nevus,
                                 "Melanoma line" = PAL$melanoma),
                      name = NULL) +
    labs(x = NULL, y = bquote(bold(italic(.(g))~log[2]*"(summed TPM + 1)")),
         subtitle = sub) +
    theme_pub(grid = TRUE) +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 30, hjust = 1))

  if (is_undetected) {
    # Fix the axis to the detection floor so the bars are shown at their true
    # (negligible) height instead of being stretched to fill the panel. This is
    # the whole point: free scaling turns noise into a convincing-looking plot.
    p <- p +
      scale_y_continuous(limits = c(0, DETECT_FLOOR),
                         expand = expansion(mult = c(0, 0.02))) +
      annotate("text", x = (length(levels(pd$label)) + 1) / 2,
               y = DETECT_FLOOR * 0.82, label = "NOT DETECTED",
               size = FIG$base_size / .pt, fontface = "bold",
               colour = PAL$grey_pt) +
      geom_hline(yintercept = DETECT_FLOOR, linetype = "11",
                 colour = PAL$grey_pt, linewidth = 0.4)
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = c(0, 0.10)))
  }

  p
}

for (g in GENES) {
  p <- make_fig(g)
  if (is.null(p)) { log_line("  ", g, ": skipped", logfile = LOG); next }

  fig <- p + plot_annotation(
    title = bquote(bold(italic(.(g))~"in melanoma cell lines versus melanocytes (GSE109245)")),
    caption = paste("Five metastatic melanoma lines and ONE melanocyte culture.",
                    "Descriptive only: n = 1 in the reference group, so no test",
                    "is performed.")) &
    theme(plot.title = element_text(size = FIG$base_size + 1, face = "bold"),
          plot.caption = element_text(size = FIG$base_size - 1,
                                      colour = PAL$grey_pt, hjust = 0))

  for (fmt in c("pdf", "png")) {
    fp <- file.path(FIG_GEO, paste0("geo_", tolower(g),
                                    "_celllines_gse109245.", fmt))
    if (file.exists(fp)) file.remove(fp)
    args <- list(filename = fp, plot = fig, width = FIG$width_1_5,
                 height = 85, units = "mm")
    if (fmt == "pdf") args$device <- grDevices::cairo_pdf
    else { args$dpi <- 600; args$bg <- "white" }
    do.call(ggplot2::ggsave, args)
  }
  log_line("  wrote: geo_", tolower(g), "_celllines_gse109245.{pdf,png}",
           logfile = LOG)

  readr::write_csv(
    d |> dplyr::select(sample, label, group, dplyr::all_of(g)),
    file.path(RES_GEO, paste0("table_geo_", tolower(g),
                              "_celllines_gse109245.csv")))
}

readr::write_csv(summ, file.path(RES_GEO, "gse109245_summary.csv"))
log_line("  wrote: results/geo/gse109245_summary.csv", logfile = LOG)

# ==============================================================================
# 6. What is in the GSE230574 archive
# ==============================================================================
hdr("6. GSE230574 archive contents (for the next script)")

tarf <- list.files(file.path(GEO_DIR, "GSE230574"), pattern = "\\.tar$",
                   full.names = TRUE)
if (length(tarf)) {
  contents <- try(utils::untar(tarf[1], list = TRUE), silent = TRUE)
  if (!inherits(contents, "try-error")) {
    log_line("  ", length(contents), " files inside ", basename(tarf[1]),
             logfile = LOG)
    for (cn in head(contents, 30)) log_line("    ", cn, logfile = LOG)
    if (length(contents) > 30) log_line("    ... and ", length(contents) - 30,
                                        " more", logfile = LOG)
  } else {
    log_line("  [FLAG] could not list the archive", logfile = LOG)
  }
} else {
  log_line("  archive not found; run script 13 to fetch it", logfile = LOG)
}

report_warnings(LOG)
write_session_info("16_celllines")
log_line("\nFull log: ", LOG, logfile = LOG)
