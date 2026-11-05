# ==============================================================================
# 19_lineage_gse230574.R
#
# THE DECISIVE CONTROL for the whole GEO arm of this project.
#
# Script 18 established that RASSF9 is ~2.5 TPM in primary nevus cells and below
# 0.14 TPM in cultured melanocytes and below 0.02 TPM in melanoma lines. That
# result is only about melanocytic biology IF the RASSF9 signal in the nevus
# samples comes from melanocytes.
#
# There is strong reason to doubt it. In the nevus tags only 23-38% of cells are
# MITF-positive and 10-15% PMEL-positive, versus 98-100% in cultured
# melanocytes. The KIT-sorted nevus preparation evidently carries substantial
# non-melanocytic content. If RASSF9's molecules come from fibroblasts,
# keratinocytes, immune or endothelial cells, then:
#
#   - "RASSF9 higher in nevus than in cultured melanocytes" means only that
#     nevus preparations contain stroma and cultures do not;
#   - the replicated "RASSF9 down in melanoma versus nevus" result in
#     GSE112509 and GSE98394 would be explained by melanomas having higher
#     tumour-cell content than nevi, not by anything happening in the
#     melanocytic lineage;
#   - and that is exactly the direction of the TCGA purity correlation already
#     found in script 06 (RASSF9 vs CPE purity rho = -0.123, p = 0.027).
#
# Three findings would collapse into one artefact. So this is tested directly:
# every nevus cell is assigned a lineage from marker genes, and RASSF9
# detection is computed WITHIN each lineage. The question is single and sharp:
# are the RASSF9-positive nevus cells melanocytes or not?
#
# WHAT WOULD FALSIFY THE MELANOCYTIC INTERPRETATION: RASSF9 detection
# concentrated in non-melanocytic cells, with melanocytic nevus cells at or
# near the cultured-melanocyte level. The script states which way the result
# went rather than leaving it to be read off a table.
#
# The cultured run doubles as a positive control for the lineage assignment
# itself: those cells are known to be ~pure melanocytic, so the classifier must
# call them melanocytic. If it does not, the classifier is broken.
#
# Output: results/geo/gse230574_lineage_assignment.csv
#         results/geo/gse230574_rassf9_by_lineage.csv
#         figures/geo/geo_rassf9_lineage_gse230574.{pdf,png} (and rassf10)
#         results/logs/19_lineage.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(patchwork))

LOG <- file.path(PATHS$logs, "19_lineage.log")
cat("GSE230574 lineage attribution —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

GEO_DIR <- file.path(PATHS$raw, "geo")
EXT_DIR <- file.path(GEO_DIR, "GSE230574", "extracted")
RES_GEO <- file.path(PATHS$results, "geo")
FIG_GEO <- file.path(PATHS$figures, "geo")

GENES <- c("RASSF9", "RASSF10")

# ---- Lineage marker panels ---------------------------------------------------
# Deliberately several markers per lineage: at single-cell depth any one gene is
# missing from most cells by chance, so a panel is needed. Markers chosen to be
# unambiguous and non-overlapping between lineages.
MARKERS <- list(
  Melanocytic  = c("MITF", "PMEL", "MLANA", "DCT", "TYR", "TYRP1", "SOX10"),
  Fibroblast   = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB"),
  Keratinocyte = c("KRT14", "KRT5", "KRT1", "KRT10", "KRT6A"),
  Immune       = c("PTPRC", "CD3E", "CD68", "LYZ", "CD14", "AIF1"),
  Endothelial  = c("PECAM1", "VWF", "CDH5", "CLDN5")
)
HOUSEKEEPING <- c("ACTB", "GAPDH")
ALL_GENES <- unique(c(GENES, unlist(MARKERS), HOUSEKEEPING))

# Lineage calling needs enough molecules per cell for the markers to be seen at
# all. Below this the cell is "undetermined" by lack of information, which is
# reported rather than silently pooled.
LINEAGE_MIN_LIB <- 2000
CHUNK_SIZE      <- 25000

# ==============================================================================
# 1. Extract all marker genes per cell
# ==============================================================================
hdr("1. Extraction")

mols_files <- sort(list.files(EXT_DIR, pattern = "MolsPerCell\\.csv\\.gz$",
                              full.names = TRUE))
tag_files  <- sort(list.files(EXT_DIR, pattern = "Sample_Tag_Calls\\.csv\\.gz$",
                              full.names = TRUE))
if (!length(mols_files)) {
  stop("No extracted matrices in ", EXT_DIR,
       "\n  Run R/17_singlecell_gse230574.R first.")
}
run_id <- function(f) sub("_(RSEC_MolsPerCell|Sample_Tag_Calls)\\.csv\\.gz$", "",
                          basename(f))

read_tags <- function(f) {
  raw <- readr::read_lines(f)
  start <- which(!startsWith(raw, "#"))[1]
  tb <- readr::read_csv(I(raw[start:length(raw)]), show_col_types = FALSE,
                        progress = FALSE)
  names(tb) <- make.names(names(tb))
  cc <- grep("^Cell_Index$|^Cell.Index$", names(tb), value = TRUE)[1]
  tc <- grep("Sample_Name|Sample.Name", names(tb), value = TRUE)[1]
  if (is.na(tc)) tc <- grep("Sample_Tag|Sample.Tag", names(tb), value = TRUE)[1]
  tibble::tibble(run = run_id(f), cell = as.character(tb[[cc]]),
                 tag = as.character(tb[[tc]]))
}
tags <- purrr::map_dfr(tag_files, read_tags) |>
  dplyr::filter(!tag %in% c("Multiplet", "Undetermined"))

extract <- function(f) {
  rid <- run_id(f)
  con <- gzfile(f, "r")
  repeat {
    ln <- readLines(con, n = 1)
    if (!length(ln)) stop("No header in ", basename(f))
    if (!startsWith(ln, "#")) break
  }
  hdr_names <- make.names(strsplit(ln, ",", fixed = TRUE)[[1]])
  close(con)

  cell_idx <- grep("^Cell_Index$|^Cell.Index$", hdr_names)[1]
  if (is.na(cell_idx)) cell_idx <- 1L
  pos <- stats::setNames(match(ALL_GENES, hdr_names), ALL_GENES)
  found <- ALL_GENES[!is.na(pos)]
  miss  <- ALL_GENES[is.na(pos)]

  log_line("  ", rid, ": ", length(found), " of ", length(ALL_GENES),
           " genes found", logfile = LOG)
  if (length(miss)) {
    log_line("    absent: ", paste(miss, collapse = ", "), logfile = LOG)
  }

  acc <- list()
  cb <- function(chunk, p) {
    names(chunk) <- hdr_names[seq_len(ncol(chunk))]
    mat <- as.matrix(chunk[, -cell_idx, drop = FALSE]); mat[is.na(mat)] <- 0
    out <- tibble::tibble(run = rid,
                          cell = as.character(chunk[[cell_idx]]),
                          lib_size = rowSums(mat))
    for (g in ALL_GENES) {
      out[[g]] <- if (is.na(pos[[g]])) NA_real_ else as.numeric(chunk[[pos[[g]]]])
    }
    acc[[length(acc) + 1L]] <<- out
    invisible(NULL)
  }
  readr::read_csv_chunked(
    f, callback = readr::SideEffectChunkCallback$new(cb),
    chunk_size = CHUNK_SIZE, comment = "#", col_names = TRUE,
    col_types = readr::cols(.default = readr::col_double()), progress = FALSE)
  dplyr::bind_rows(acc)
}

expr <- purrr::map_dfr(mols_files, extract)
d <- dplyr::inner_join(tags, expr, by = c("run", "cell"))
log_line("  cells with tag and counts: ", nrow(d), logfile = LOG)

d <- d |> dplyr::mutate(
  class = dplyr::case_when(grepl("^SampleTag", tag)  ~ "Nevus (primary)",
                           grepl("^melanocyte", tag) ~ "Melanocyte (cultured)",
                           TRUE                      ~ "Melanoma line"))

# ==============================================================================
# 2. Lineage assignment
# ==============================================================================
hdr("2. Lineage assignment")

score_panel <- function(df, genes) {
  g <- intersect(genes, names(df))
  if (!length(g)) return(rep(0, nrow(df)))
  m <- as.matrix(df[, g, drop = FALSE]); m[is.na(m)] <- 0
  rowSums(m > 0)                       # number of panel markers detected
}

for (nm in names(MARKERS)) {
  d[[paste0("n_", nm)]] <- score_panel(d, MARKERS[[nm]])
}
score_cols <- paste0("n_", names(MARKERS))

# A cell is assigned to the lineage with the most markers detected. Ties
# between a melanocytic and a non-melanocytic panel are called "Ambiguous"
# rather than resolved in either direction -- resolving them towards
# melanocytic would bias the very question being asked.
sc <- as.matrix(d[, score_cols])
colnames(sc) <- names(MARKERS)
best_n   <- apply(sc, 1, max)
n_at_max <- apply(sc, 1, function(x) sum(x == max(x)))
best_nm  <- names(MARKERS)[apply(sc, 1, which.max)]

d$lineage <- dplyr::case_when(
  d$lib_size < LINEAGE_MIN_LIB ~ "Undetermined (shallow)",
  best_n == 0                  ~ "Undetermined (no markers)",
  n_at_max > 1                 ~ "Ambiguous",
  TRUE                         ~ best_nm)

tab <- d |> dplyr::count(class, lineage, name = "n_cells") |>
  dplyr::group_by(class) |>
  dplyr::mutate(pct = 100 * n_cells / sum(n_cells)) |>
  dplyr::ungroup() |> dplyr::arrange(class, dplyr::desc(n_cells))

log_line(sprintf("  %-24s %-26s %8s %8s", "class", "lineage", "cells", "pct"),
         logfile = LOG)
for (i in seq_len(nrow(tab))) {
  log_line(sprintf("  %-24s %-26s %8d %7.1f%%", tab$class[i], tab$lineage[i],
                   tab$n_cells[i], tab$pct[i]), logfile = LOG)
}
readr::write_csv(tab, file.path(RES_GEO, "gse230574_lineage_assignment.csv"))

# ---- Positive control for the classifier ------------------------------------
cult <- tab |> dplyr::filter(class != "Nevus (primary)", lineage == "Melanocytic")
cult_pct <- d |> dplyr::filter(class != "Nevus (primary)") |>
  dplyr::summarise(p = 100 * mean(lineage == "Melanocytic")) |> dplyr::pull(p)
log_line("", logfile = LOG)
if (cult_pct >= 80) {
  log_line(sprintf("  CLASSIFIER PASSES: %.1f%% of cultured cells called", cult_pct),
           logfile = LOG)
  log_line("  Melanocytic, as they must be. The assignment can be trusted.",
           logfile = LOG)
} else {
  log_line(sprintf("  [FLAG] CLASSIFIER FAILS: only %.1f%% of cultured cells", cult_pct),
           logfile = LOG)
  log_line("  called Melanocytic, though they are pure cultures. The marker",
           logfile = LOG)
  log_line("  panels or the depth threshold are wrong. Do not interpret the",
           logfile = LOG)
  log_line("  nevus result below until this is fixed.", logfile = LOG)
}

# ==============================================================================
# 3. The decisive question: where do the RASSF9 molecules come from?
# ==============================================================================
hdr("3. Gene detection within lineage")

by_lin <- list()
for (g in GENES) {
  r <- d |>
    dplyr::group_by(class, lineage) |>
    dplyr::summarise(n_cells = dplyr::n(),
                     n_pos = sum(.data[[g]] > 0, na.rm = TRUE),
                     frac_pos = n_pos / n_cells,
                     molecules = sum(.data[[g]], na.rm = TRUE),
                     total_mol = sum(lib_size),
                     tpm = 1e6 * molecules / sum(lib_size),
                     .groups = "drop") |>
    dplyr::mutate(gene = g, .before = 1) |>
    dplyr::arrange(class, dplyr::desc(frac_pos))
  by_lin[[g]] <- r

  log_line("", logfile = LOG)
  log_line("  ", g, logfile = LOG)
  log_line(sprintf("    %-24s %-26s %7s %7s %9s %8s",
                   "class", "lineage", "cells", "%pos", "molecules", "TPM"),
           logfile = LOG)
  for (i in seq_len(nrow(r))) {
    log_line(sprintf("    %-24s %-26s %7d %6.2f%% %9d %8.3f",
                     r$class[i], r$lineage[i], r$n_cells[i],
                     100 * r$frac_pos[i], r$molecules[i], r$tpm[i]),
             logfile = LOG)
  }
}
by_lin_all <- dplyr::bind_rows(by_lin)
readr::write_csv(by_lin_all, file.path(RES_GEO, "gse230574_rassf9_by_lineage.csv"))

# ---- Where the nevus molecules actually came from ---------------------------
hdr("4. Verdict")

for (g in GENES) {
  nv <- by_lin[[g]] |> dplyr::filter(class == "Nevus (primary)")
  tot <- sum(nv$molecules)
  log_line("", logfile = LOG)
  log_line("  ", g, ": ", tot, " molecules in the nevus samples, split by lineage:",
           logfile = LOG)
  if (tot == 0) { log_line("    none to attribute", logfile = LOG); next }
  nv2 <- nv |> dplyr::mutate(share = 100 * molecules / tot) |>
    dplyr::arrange(dplyr::desc(share))
  for (i in seq_len(nrow(nv2))) {
    log_line(sprintf("    %-26s %5.1f%% of molecules   (%6.3f TPM within lineage)",
                     nv2$lineage[i], nv2$share[i], nv2$tpm[i]), logfile = LOG)
  }

  mel_tpm <- nv$tpm[nv$lineage == "Melanocytic"]
  non_lin <- c("Fibroblast", "Keratinocyte", "Immune", "Endothelial")
  non <- nv |> dplyr::filter(lineage %in% non_lin)
  non_tpm <- if (nrow(non)) 1e6 * sum(non$molecules) / sum(non$total_mol) else NA_real_
  cult_tpm <- by_lin[[g]] |>
    dplyr::filter(class != "Nevus (primary)", lineage == "Melanocytic") |>
    dplyr::summarise(t = 1e6 * sum(molecules) / sum(total_mol)) |> dplyr::pull(t)

  log_line("", logfile = LOG)
  log_line(sprintf("    nevus melanocytic      %7.3f TPM",
                   ifelse(length(mel_tpm), mel_tpm, NA)), logfile = LOG)
  log_line(sprintf("    nevus non-melanocytic  %7.3f TPM", non_tpm), logfile = LOG)
  log_line(sprintf("    cultured melanocytic   %7.3f TPM", cult_tpm), logfile = LOG)

  if (!length(mel_tpm) || is.na(mel_tpm)) {
    log_line("    -> no melanocytic nevus cells to judge from", logfile = LOG)
  } else if (mel_tpm > 5 * max(cult_tpm, 1e-6) &&
             (is.na(non_tpm) || mel_tpm >= non_tpm)) {
    log_line("    -> SURVIVES. The signal is in melanocytic nevus cells, at a",
             logfile = LOG)
    log_line("       level well above cultured melanocytes and not lower than",
             logfile = LOG)
    log_line("       in non-melanocytic cells. Tissue composition does NOT",
             logfile = LOG)
    log_line("       explain it, and the melanocytic reading stands.",
             logfile = LOG)
  } else if (!is.na(non_tpm) && non_tpm > 3 * mel_tpm) {
    log_line("    -> FAILS. The molecules are concentrated in NON-melanocytic",
             logfile = LOG)
    log_line("       cells. This is a tissue-composition effect, and the",
             logfile = LOG)
    log_line("       nevus-versus-melanoma results in GSE112509 / GSE98394",
             logfile = LOG)
    log_line("       must be re-read in that light before any claim is made.",
             logfile = LOG)
  } else {
    log_line("    -> EQUIVOCAL. Melanocytic and non-melanocytic levels are",
             logfile = LOG)
    log_line("       comparable, so this dataset cannot separate the two",
             logfile = LOG)
    log_line("       explanations. Report it as unresolved.", logfile = LOG)
  }
}

# ==============================================================================
# 5. Figures
# ==============================================================================
hdr("5. Figures")

# PAL_LINEAGE comes from 00_setup.R, where it was selected by colour-vision
# simulation rather than by eye (worst-case dE 23.1 across protanopia,
# deuteranopia and tritanopia). Both panels run horizontally so the lineage
# names read at full size without rotation.
LIN_ORDER <- names(PAL_LINEAGE)

make_fig <- function(g) {
  r <- by_lin[[g]]
  if (is.null(r) || !nrow(r)) return(NULL)

  # Panel A -- composition. A stacked proportion bar per class answers the
  # question the reader has first: how melanocytic is each preparation?
  comp <- tab |>
    dplyr::mutate(lineage = factor(lineage, levels = rev(LIN_ORDER)),
                  class = factor(class, levels = rev(names(PAL_CLASS))))
  pa <- ggplot(comp, aes(y = class, x = pct, fill = lineage)) +
    geom_col(width = 0.6, colour = "white", linewidth = 0.25) +
    scale_fill_manual(values = PAL_LINEAGE, name = NULL, drop = FALSE,
                      breaks = LIN_ORDER) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.01)),
                       breaks = seq(0, 100, 25)) +
    coord_cartesian(clip = "off") +
    labs(y = NULL, x = "Cells (%)") +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
    theme_journal(grid = FALSE) +
    theme(legend.position = "top",
          legend.justification = "left",
          legend.text = element_text(size = FIG$base_size - 1.5),
          legend.key.size = unit(8, "pt"),
          axis.text.y = element_text(size = FIG$base_size - 0.5, hjust = 1),
          plot.margin = margin(4, 10, 4, 4, unit = "pt"))

  # Panel B -- the decisive panel. Undetermined cells are excluded because a
  # cell with no markers carries no information about lineage.
  rr <- r |>
    dplyr::filter(!grepl("^Undetermined", lineage)) |>
    dplyr::mutate(lineage = factor(lineage, levels = rev(LIN_ORDER)),
                  class = factor(class, levels = names(PAL_CLASS)))
  # Value labels sit inside the long bars and outside the short ones, so a
  # label can never run past the panel edge whatever the largest value is.
  rr <- rr |> dplyr::mutate(
    lbl = ifelse(tpm >= 0.005, sprintf("%.2f", tpm), "0"),
    inside = tpm >= 0.45 * max(tpm))
  pb <- ggplot(rr, aes(y = lineage, x = tpm, fill = lineage)) +
    geom_col(width = 0.68, colour = PAL$ink, linewidth = 0.25) +
    geom_text(data = dplyr::filter(rr, inside), aes(label = lbl),
              hjust = 1.15, colour = "white", fontface = "bold",
              size = (FIG$base_size - 2) / .pt) +
    geom_text(data = dplyr::filter(rr, !inside), aes(label = lbl),
              hjust = -0.18, colour = PAL$ink,
              size = (FIG$base_size - 2) / .pt) +
    facet_grid(class ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_manual(values = PAL_LINEAGE, guide = "none", drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.24)), n.breaks = 5) +
    coord_cartesian(clip = "off") +
    labs(y = NULL, x = bquote(bold(italic(.(g))~"(TPM within lineage)"))) +
    theme_journal(grid = FALSE) +
    theme(panel.grid.major.x = element_line(linewidth = 0.25, colour = "grey92"),
          panel.spacing.y = unit(5, "pt"),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 0,
                                           size = FIG$base_size - 1),
          axis.text.y = element_text(size = FIG$base_size - 1, hjust = 1),
          plot.margin = margin(4, 10, 4, 4, unit = "pt"))

  (pa / pb) + patchwork::plot_layout(heights = c(1, 2.3))
}

for (g in GENES) {
  p <- make_fig(g)
  if (is.null(p)) { log_line("  ", g, ": skipped", logfile = LOG); next }
  fig <- p + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = FIG$base_size + 3, face = "bold",
                                  colour = PAL$ink))

  stem <- paste0("geo_", tolower(g), "_lineage_gse230574")
  save_fig_to(fig, stem, FIG_GEO, width = FIG$width_double, height = 175)

  # The legend states the attribution outcome in numbers, because this figure
  # is the one that decides whether the gene's nevus signal is melanocytic.
  nv <- by_lin[[g]] |> dplyr::filter(class == "Nevus (primary)")
  tot <- sum(nv$molecules)
  mel <- nv$tpm[nv$lineage == "Melanocytic"]
  non <- nv |> dplyr::filter(lineage %in% c("Fibroblast", "Keratinocyte",
                                            "Immune", "Endothelial"))
  non_tpm <- if (nrow(non)) 1e6 * sum(non$molecules) / sum(non$total_mol) else NA
  mel_mol <- nv$molecules[nv$lineage == "Melanocytic"]
  mel_share <- if (tot > 0 && length(mel_mol)) 100 * mel_mol / tot else NA_real_
  cult_tpm_g <- by_lin[[g]] |>
    dplyr::filter(class != "Nevus (primary)", lineage == "Melanocytic") |>
    dplyr::summarise(t = 1e6 * sum(molecules) / sum(total_mol)) |> dplyr::pull(t)
  nv_assigned <- d |> dplyr::filter(class == "Nevus (primary)") |>
    dplyr::summarise(p = 100 * mean(!grepl("^Undetermined", lineage))) |>
    dplyr::pull(p)

  add_legend(
    stem,
    paste0("Figure. Lineage attribution of ", g, " in GSE230574."),
    paste0(
      "Each cell was assigned to the lineage for which the most panel markers ",
      "were detected, from five panels: melanocytic (MITF, PMEL, MLANA, DCT, ",
      "TYR, TYRP1, SOX10), fibroblast (COL1A1, COL1A2, COL3A1, DCN, LUM, ",
      "PDGFRB), keratinocyte (KRT14, KRT5, KRT1, KRT10, KRT6A), immune ",
      "(PTPRC, CD3E, CD68, LYZ, CD14, AIF1) and endothelial (PECAM1, VWF, ",
      "CDH5, CLDN5). Ties between panels were called ambiguous rather than ",
      "resolved, and cells below ", LINEAGE_MIN_LIB,
      " molecules were left undetermined for want of information. The ",
      "cultured samples are the positive control for the classifier, which ",
      "assigned ", sprintf("%.1f%%", cult_pct), " of them to the melanocytic ",
      "lineage as expected. (A) Lineage composition of each preparation. ",
      "(B) ", g, " abundance within each lineage, in molecules per million, ",
      "excluding undetermined cells. In the nevus samples ",
      sprintf("%.1f%%", ifelse(length(mel_share), mel_share, NA)),
      " of the ", tot, " ", g, " molecules were attributable to melanocytic ",
      "cells (", sprintf("%.3f", ifelse(length(mel), mel, NA)),
      " TPM within that lineage) against ", sprintf("%.3f", non_tpm),
      " TPM across the non-melanocytic lineages combined, with cultured ",
      "melanocytes at ", sprintf("%.3f", cult_tpm_g), " TPM. ",
      "Attribution rests on small molecule counts and on the ",
      sprintf("%.1f%%", nv_assigned), " of nevus cells deep enough to assign, ",
      "so the direction is interpretable but the magnitudes are not precise."))

  log_line("  wrote: ", stem, ".{pdf,png,svg}", logfile = LOG)
}

readr::write_csv(
  d |> dplyr::select(run, tag, class, cell, lib_size, lineage,
                     dplyr::all_of(score_cols), dplyr::all_of(GENES)),
  file.path(RES_GEO, "table_geo_lineage_percell_gse230574.csv"))
log_line("  wrote: results/geo/table_geo_lineage_percell_gse230574.csv",
         logfile = LOG)

report_warnings(logfile = LOG)
write_session_info("19_lineage")
log_line("\nFull log: ", LOG, logfile = NULL)
