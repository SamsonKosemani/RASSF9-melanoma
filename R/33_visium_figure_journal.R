# ==============================================================================
# 33_visium_figure_journal.R
#
# The publication figure: RASSF9 against PRAME in situ, drawn to journal
# standard.
#
# WHY NOT semla's MapFeatures. Script 32's panels were right scientifically and
# wrong as figures: titles baked into the image, colour bars scaled per panel so
# nothing is comparable, section numbers that mean nothing to a reader, and a
# white-to-red ramp on white paper where low values vanish. Those are not
# cosmetic complaints -- a per-panel scale makes a cross-panel comparison
# misleading, which is the one thing this figure exists to do.
#
# So the spots are drawn directly with ggplot2 from the Visium coordinates. The
# Seurat object is NOT needed: coordinates come from the GEO tissue_positions
# files, counts from the GEO matrices, and compartments from the authors' own
# labels. Every element is then under our control.
#
# WHAT THE FIGURE SHOWS, and why it is built this way:
#
#   PANEL A -- the compartments, from Kreuger et al.'s own spot assignments.
#   The reader sees WHERE the nevus and melanoma are before being shown what is
#   expressed in them. In script 32's version PRAME was doing double duty as
#   both the map and the finding, which is circular.
#
#   PANEL B -- RASSF9.  PANEL C -- PRAME.
#   Same section, same spots, same transformation. PRAME fills the melanoma
#   compartment; RASSF9 does not follow it. That contrast is the figure.
#
# DECISIONS THAT MATTER FOR HONESTY, not just for looks:
#
#   ONE SCALE PER GENE, SHARED ACROSS SECTIONS. In the supplementary figure a
#   gene's colour scale is computed once over every section drawn, so a spot of
#   a given colour means the same thing everywhere. Per-section scaling would
#   let a section with little signal look as intense as one with a lot.
#
#   SCALES ARE NOT SHARED BETWEEN GENES. RASSF9 and PRAME differ roughly
#   thirty-fold in abundance; forcing one scale would render RASSF9 invisible
#   and the figure would say more about dynamic range than about localisation.
#   The legend states this explicitly rather than leaving the reader to assume.
#
#   RAW COUNTS, normalised to counts per 10,000 within each spot, then log1p.
#   Script 32 plotted the object's SCT values, which are their preprocessing;
#   this matches the quantification in script 30 instead.
#
#   HIGH VALUES DRAWN LAST, so a few positive spots are not buried under the
#   thousands of zeros. Drawing order changes what a reader sees and is
#   therefore declared, not left to the order rows happen to arrive in.
#
#   A REAL SCALE BAR, computed from each section's own scalefactors_json: a
#   Visium spot is 55 um, so spot_diameter_fullres gives micrometres per pixel.
#
# JOURNAL MODE: no titles, subtitles or captions inside the image. Panel tags
# only. Legends are written to results/figure_legends_geo.md for the manuscript
# text.
#
# Output: figures/geo/fig_visium_main_<section>.{pdf,png,svg}
#         figures/geo/fig_visium_supp_<gene>.{pdf,png,svg}
#         results/logs/33_visium_figure_journal.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(Matrix)
  library(patchwork)
})

LOG <- file.path(PATHS$logs, "33_visium_figure_journal.log")
cat("GSE298774 journal figure —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

EXT_DIR <- file.path(PATHS$raw, "geo", "GSE298774", "extracted")
PROC    <- PATHS$processed
FIG_GEO <- file.path(PATHS$figures, "geo")
dir.create(FIG_GEO, recursive = TRUE, showWarnings = FALSE)

MAIN_SECTION <- "B4"                       # clearest melanoma mass of the four
SUPP_SECTIONS <- c("A16", "A21", "B1", "B4")
GENES_MAIN <- c("RASSF9", "PRAME")
GENES_SUPP <- c("RASSF9", "RASSF10")
SPOT_UM <- 55                              # Visium spot diameter
BAR_UM  <- 1000                            # 1 mm scale bar

# ---- Palette ----------------------------------------------------------------
# Compartments use the CVD-validated three-class set already in use elsewhere in
# this project, so the figure is consistent with the others.
COL_NEVUS <- "#0095D9"
COL_TUM   <- "#F2680C"
COL_OTHER <- "#C2C2C2"
COL_NONE  <- "#EDEDED"   # in the section but not in the authors' object
COL_INK   <- "#1A1A1A"

# Sequential ramp for expression. A single-hue ramp that rises in saturation AND
# falls in lightness stays readable under every form of colour-vision
# deficiency, because lightness alone carries the ordering. The pale end is a
# light warm grey rather than white, so zero spots still render the tissue
# outline on white paper -- the failure in script 32's panels.
RAMP_EXPR <- c("#EDEDEA", "#FCD9A8", "#F59B42", "#E05C00", "#8C2D04")

# ==============================================================================
# 1. Coordinates
# ==============================================================================
hdr("1. Coordinates")

files <- list.files(EXT_DIR, full.names = TRUE)
sec_of <- function(b) sub("^GSM\\d+_(.+?)_(aligned_fiducials|barcodes|features|matrix|scalefactors_json|tissue_hires_image|tissue_lowres_image|tissue_positions).*$", "\\1", b)
inv <- tibble::tibble(path = files, base = basename(files)) |>
  dplyr::mutate(section = sec_of(base),
                role = dplyr::case_when(
                  grepl("_matrix\\.mtx", base)      ~ "matrix",
                  grepl("_barcodes\\.tsv", base)    ~ "barcodes",
                  grepl("_features\\.tsv", base)    ~ "features",
                  grepl("_tissue_positions", base)  ~ "positions",
                  grepl("_scalefactors_json", base) ~ "scalefactors",
                  TRUE ~ "other")) |>
  dplyr::filter(role != "other")
g2 <- function(sec, role) {
  p <- inv$path[inv$section == sec & inv$role == role]
  if (!length(p)) NA_character_ else p[1]
}
core16 <- function(x) sub("-\\d+$", "", x)

read_pos <- function(p) {
  if (is.na(p)) return(NULL)
  first <- readr::read_lines(p, n_max = 1)
  has_hdr <- grepl("barcode", first, ignore.case = TRUE)
  tb <- readr::read_csv(p, col_names = has_hdr, show_col_types = FALSE,
                        progress = FALSE)
  if (!has_hdr) {
    nm <- c("barcode", "in_tissue", "array_row", "array_col",
            "pxl_row_in_fullres", "pxl_col_in_fullres")
    names(tb)[seq_len(min(6, ncol(tb)))] <- nm[seq_len(min(6, ncol(tb)))]
  }
  bc <- grep("barcode", names(tb), ignore.case = TRUE, value = TRUE)[1]
  it <- grep("in_?tissue", names(tb), ignore.case = TRUE, value = TRUE)[1]
  yr <- grep("pxl_row", names(tb), ignore.case = TRUE, value = TRUE)[1]
  xc <- grep("pxl_col", names(tb), ignore.case = TRUE, value = TRUE)[1]
  if (any(is.na(c(bc, it, yr, xc)))) return(NULL)
  tibble::tibble(barcode = as.character(tb[[bc]]),
                 in_tissue = as.integer(tb[[it]]),
                 y = as.numeric(tb[[yr]]), x = as.numeric(tb[[xc]]))
}

# Micrometres per pixel, from the section's own scale factors. A spot is 55 um
# across, so its diameter in full-resolution pixels fixes the conversion.
um_per_px <- function(sec) {
  p <- g2(sec, "scalefactors")
  if (is.na(p)) return(NA_real_)
  txt <- paste(readr::read_lines(p), collapse = " ")
  m <- regmatches(txt, regexpr('"spot_diameter_fullres"\\s*:\\s*[0-9.]+', txt))
  if (!length(m)) return(NA_real_)
  d <- as.numeric(sub('.*:\\s*', '', m))
  if (!is.finite(d) || d <= 0) return(NA_real_)
  SPOT_UM / d
}

SECTIONS <- union(MAIN_SECTION, SUPP_SECTIONS)
POS <- list()
for (sec in SECTIONS) {
  p <- read_pos(g2(sec, "positions"))
  if (is.null(p)) {
    log_line("  [FLAG] ", sec, ": no usable tissue_positions; skipped",
             logfile = LOG)
    next
  }
  p <- p |>
    dplyr::filter(in_tissue == 1L) |>
    dplyr::mutate(section = sec, spot = paste0(sec, "|", core16(barcode)))
  POS[[sec]] <- p
  log_line(sprintf("  %-6s %5d in-tissue spots | %.3f um/px",
                   sec, nrow(p), um_per_px(sec)), logfile = LOG)
}
if (!length(POS)) stop("No coordinates could be read.")
XY <- dplyr::bind_rows(POS)

# ==============================================================================
# 2. Counts and compartments
# ==============================================================================
hdr("2. Counts and compartments")

WANT <- unique(c(GENES_MAIN, GENES_SUPP))
vals <- list()
for (sec in names(POS)) {
  mp <- g2(sec, "matrix"); bp <- g2(sec, "barcodes"); fp <- g2(sec, "features")
  if (any(is.na(c(mp, bp, fp)))) next
  m  <- Matrix::readMM(mp)
  bc <- readr::read_tsv(bp, col_names = FALSE, show_col_types = FALSE,
                        progress = FALSE)[[1]]
  ft <- readr::read_tsv(fp, col_names = FALSE, show_col_types = FALSE,
                        progress = FALSE)
  if (ncol(m) == nrow(ft) && nrow(m) == length(bc)) m <- Matrix::t(m)
  if (nrow(m) != nrow(ft) || ncol(m) != length(bc)) {
    log_line("  [FLAG] ", sec, ": dimension mismatch; skipped", logfile = LOG)
    rm(m); next
  }
  sy  <- toupper(as.character(ft[[min(2, ncol(ft))]]))
  lib <- as.numeric(Matrix::colSums(m))
  out <- tibble::tibble(spot = paste0(sec, "|", core16(bc)), lib_size = lib)
  for (g in WANT) {
    i <- which(sy == g)
    v <- if (!length(i)) {
      NA_real_
    } else if (length(i) > 1) {
      as.numeric(Matrix::colSums(m[i, , drop = FALSE]))
    } else {
      as.numeric(m[i, ])
    }
    # Counts per 10,000 within spot, then log1p. Stated in the legend.
    out[[g]] <- log1p(1e4 * v / pmax(lib, 1))
  }
  vals[[sec]] <- out
  rm(m); gc(verbose = FALSE)
}
E <- dplyr::bind_rows(vals)

LAB <- readRDS(file.path(PROC, "gse298774_kreuger_labels.rds"))$labels |>
  dplyr::select(spot, label)

D <- XY |>
  dplyr::inner_join(E, by = "spot") |>
  dplyr::left_join(LAB, by = "spot") |>
  # FOUR categories, not three. A spot with no label is not the same as a spot
  # the authors assigned to keratinocytes or fibroblasts: their object excludes
  # 15-99% of in-tissue spots per section, most likely by depth. Drawing both in
  # one grey would tell the reader the tissue was fully classified when about a
  # fifth of it was never in their object at all.
  dplyr::mutate(compartment = dplyr::case_when(
    is.na(label) ~ "Not assigned",
    label == "Nevus" ~ "Nevus",
    grepl("^Tumor", label) ~ "Melanoma",
    TRUE ~ "Other tissue"),
    compartment = factor(compartment,
                         levels = c("Nevus", "Melanoma", "Other tissue",
                                    "Not assigned")))
log_line("  spots with coordinates and counts: ", nrow(D), logfile = LOG)
log_line("  unlabelled spots (drawn as other tissue): ",
         sum(is.na(D$label)), logfile = LOG)

# ---- Shared scale limits, computed ONCE per gene across all sections drawn ---
LIMS <- lapply(WANT, function(g) {
  v <- D[[g]][is.finite(D[[g]])]
  if (!length(v)) return(c(0, 1))
  c(0, as.numeric(stats::quantile(v, 0.995, na.rm = TRUE)))
})
names(LIMS) <- WANT
log_line("", logfile = LOG)
log_line("  colour scale upper limits (99.5th percentile, shared across sections):",
         logfile = LOG)
for (g in WANT) {
  log_line(sprintf("    %-8s 0 to %.2f   [max observed %.2f]", g, LIMS[[g]][2],
                   max(D[[g]], na.rm = TRUE)), logfile = LOG)
}
log_line("  the 99.5th percentile, not the maximum, so a handful of extreme",
         logfile = LOG)
log_line("  spots cannot compress the rest of the range into invisibility",
         logfile = LOG)

# ==============================================================================
# 3. Drawing
# ==============================================================================
hdr("3. Drawing")

# Spacing is set so that nothing can collide: a top margin that reserves room
# for the panel tag above the tissue, a bottom margin that reserves room for the
# scale bar and its label, and legend text small enough that three compartment
# labels fit one line across a third of a double-column figure.
theme_spatial <- function() {
  theme_void(base_family = FIG$base_family) +
    theme(legend.position = "bottom",
          legend.justification = "center",
          legend.title = element_text(size = FIG$base_size - 2,
                                      colour = COL_INK, hjust = 0.5),
          legend.text = element_text(size = FIG$base_size - 3.5,
                                     colour = COL_INK),
          legend.key.height = unit(3.0, "pt"),
          legend.key.width = unit(15, "pt"),
          legend.spacing.x = unit(2, "pt"),
          legend.spacing.y = unit(1, "pt"),
          legend.margin = margin(2, 2, 0, 2, unit = "pt"),
          legend.box.spacing = unit(3, "pt"),
          legend.box.margin = margin(0, 2, 0, 2, unit = "pt"),
          # A wider side margin keeps each panel's colourbar inside its own
          # column. At 6 pt the RASSF9 bar's upper tick and the PRAME bar's
          # lower tick met between panels B and C and read as one number.
          plot.margin = margin(12, 14, 4, 14, unit = "pt"))
}

# A scale bar in data units, from the section's own pixel size.
# The bar sits below the tissue with its label below the bar. Under
# scale_y_reverse a LARGER y draws LOWER, so both offsets are positive and the
# label never rides back over the bar -- the first version used vjust against a
# reversed axis, which is how text ends up on top of a line.
scale_bar <- function(d, sec) {
  upp <- um_per_px(sec)
  if (!is.finite(upp)) return(NULL)
  len <- BAR_UM / upp
  span <- diff(range(d$y))
  x0 <- max(d$x)
  x1 <- x0 - len
  y_bar <- max(d$y) + 0.055 * span
  y_txt <- max(d$y) + 0.115 * span
  list(annotate("segment", x = x1, xend = x0, y = y_bar, yend = y_bar,
                linewidth = 0.7, colour = COL_INK, lineend = "butt"),
       annotate("text", x = (x1 + x0) / 2, y = y_txt, label = "1 mm",
                vjust = 0.5, hjust = 0.5,
                size = (FIG$base_size - 3.5) / .pt,
                colour = COL_INK, family = FIG$base_family))
}

panel_compartment <- function(d, sec, pt) {
  # Other tissue first, so the compartments of interest sit on top.
  d <- d |> dplyr::arrange(compartment %in% c("Nevus", "Melanoma"))
  ggplot(d, aes(x, y, colour = compartment)) +
    geom_point(size = pt, shape = 16) +
    scale_colour_manual(values = c(Nevus = COL_NEVUS, Melanoma = COL_TUM,
                                   `Other tissue` = COL_OTHER,
                                   `Not assigned` = COL_NONE),
                        name = NULL, drop = FALSE) +
    # ONE ITEM PER ROW. Laid out horizontally this legend is wider than a
    # third of a 180 mm figure, so it ran sideways into the neighbouring
    # panel's colourbar title. Stacked, its width is bounded by the longest
    # single label ("Other tissue") and it cannot reach the next panel.
    guides(colour = guide_legend(override.aes = list(size = 1.5),
                                 ncol = 1, byrow = TRUE,
                                 keywidth = unit(6, "pt"),
                                 keyheight = unit(6, "pt"))) +
    scale_y_reverse(expand = expansion(mult = c(0.03, 0.17))) +
    scale_x_continuous(expand = expansion(mult = 0.04)) +
    coord_fixed(clip = "off") +
    scale_bar(d, sec) +
    theme_spatial()
}

panel_gene <- function(d, g, sec, pt) {
  # Highest values last, so positive spots are never hidden under zeros.
  d <- d[order(d[[g]], na.last = FALSE), ]
  ggplot(d, aes(x, y, colour = .data[[g]])) +
    geom_point(size = pt, shape = 16) +
    scale_colour_gradientn(colours = RAMP_EXPR, limits = LIMS[[g]],
                           oob = scales::squish,
                           name = bquote(italic(.(g))),
                           breaks = scales::pretty_breaks(3)) +
    guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                    ticks.colour = NA, frame.colour = NA)) +
    scale_y_reverse(expand = expansion(mult = c(0.03, 0.17))) +
    scale_x_continuous(expand = expansion(mult = 0.04)) +
    coord_fixed(clip = "off") +
    scale_bar(d, sec) +
    theme_spatial() +
    theme(legend.direction = "horizontal")
}

# Spot size has to suit the section: a dense 11 mm capture area needs smaller
# points than a sparse one, or the spots merge into a solid mass.
pt_for <- function(n) max(0.18, min(0.55, 900 / n))

reset_legends("# GSE298774 figure legends\n")

# ---- Main figure ------------------------------------------------------------
if (MAIN_SECTION %in% names(POS)) {
  d <- D |> dplyr::filter(section == MAIN_SECTION)
  pt <- pt_for(nrow(d))
  log_line("  main figure: ", MAIN_SECTION, " (", nrow(d), " spots, point size ",
           round(pt, 2), ")", logfile = LOG)

  pa <- panel_compartment(d, MAIN_SECTION, pt)
  pb <- panel_gene(d, GENES_MAIN[1], MAIN_SECTION, pt)
  pc <- panel_gene(d, GENES_MAIN[2], MAIN_SECTION, pt)

  fig <- (pa | pb | pc) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = FIG$base_size + 2, face = "bold",
                                  colour = COL_INK, hjust = 0, vjust = 1),
          plot.tag.position = c(0.005, 0.995))

  stem <- paste0("fig_visium_main_", MAIN_SECTION)
  save_fig_to(fig, stem, FIG_GEO, width = FIG$width_double, height = 94)
  log_line("  wrote: ", stem, logfile = LOG)

  n_nev <- sum(d$compartment == "Nevus")
  n_tum <- sum(d$compartment == "Melanoma")
  log_line("  panel A: ", n_nev, " nevus, ", n_tum, " melanoma, ",
           sum(d$compartment == "Other tissue"), " other identities, ",
           sum(d$compartment == "Not assigned"), " not in their object",
           logfile = LOG)
  add_legend(stem,
    paste0("Figure. RASSF9 is not enriched in the melanoma compartment of a ",
           "nevus-associated melanoma."),
    paste0(
      "Spatial transcriptomics of lesion ", MAIN_SECTION, " (GSE298774, ",
      "Kreuger et al.), ", nrow(d), " spots under tissue. ",
      "(A) Spot identities as assigned by the original authors: ", n_nev,
      " nevus spots and ", n_tum, " melanoma spots, with their other twelve ",
      "identities in mid grey. Spots in pale grey lie under the tissue but are ",
      "absent from the authors' object, which retains between 15% and 99% of ",
      "in-tissue spots per section, and are shown so that classified and ",
      "unclassified tissue are not confused. ",
      "(B) RASSF9. (C) PRAME, a melanoma marker the authors ",
      "report as positive in melanoma and negative in nevus in every case. ",
      "Counts are expressed per 10,000 within each spot and log-transformed. ",
      "Colour scales are shared across sections within a gene, but NOT between ",
      "genes: PRAME is roughly thirtyfold more abundant than RASSF9, and a ",
      "common scale would render RASSF9 invisible. Scale bar 1 mm. PRAME marks ",
      "the melanoma compartment densely while RASSF9 is distributed across the ",
      "section without enrichment there, consistent with the quantification ",
      "showing that only 17.1% of RASSF9 molecules lie in melanocytic spots. ",
      "Note that RASSF9 is detected in 11% of spots overall, so its scattered ",
      "appearance reflects sparsity as well as distribution; this panel ",
      "illustrates the quantitative result rather than establishing it."))
}

# ---- Supplementary: one gene across sections, one shared scale --------------
for (g in GENES_SUPP) {
  secs <- intersect(SUPP_SECTIONS, names(POS))
  if (!length(secs)) next
  ps <- lapply(secs, function(sec) {
    d <- D |> dplyr::filter(section == sec)
    panel_gene(d, g, sec, pt_for(nrow(d))) +
      labs(tag = sec) +
      theme(plot.tag = element_text(size = FIG$base_size - 1, face = "plain",
                                    colour = COL_INK, hjust = 0, vjust = 1),
            plot.tag.position = c(0.005, 0.995))
  })
  fig <- patchwork::wrap_plots(ps, nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  stem <- paste0("fig_visium_supp_", tolower(g))
  save_fig_to(fig, stem, FIG_GEO, width = FIG$width_double, height = 86)
  log_line("  wrote: ", stem, " (", paste(secs, collapse = ", "), ")",
           logfile = LOG)

  add_legend(stem,
    paste0("Supplementary Figure. ", g,
           " across four nevus-associated melanomas."),
    paste0(
      g, " in sections ", paste(secs, collapse = ", "),
      " (GSE298774). Counts per 10,000 within spot, log-transformed. A single ",
      "colour scale is used across all four sections, so a given colour means ",
      "the same expression level in every panel; the scale is capped at the ",
      "99.5th percentile so that a few extreme spots do not compress the ",
      "range. Scale bar 1 mm."))
}

# ==============================================================================
# 4. What was decided, and what is left to judgement
# ==============================================================================
hdr("4. Figure decisions")

log_line("  FIXED HERE, and stated in the legends:", logfile = LOG)
log_line("    - one colour scale per gene, shared across sections, capped at",
         logfile = LOG)
log_line("      the 99.5th percentile", logfile = LOG)
log_line("    - scales NOT shared between genes, because PRAME is ~30x more",
         logfile = LOG)
log_line("      abundant and a common scale would hide RASSF9 entirely",
         logfile = LOG)
log_line("    - raw counts per 10,000 within spot, log1p -- matching script 30,",
         logfile = LOG)
log_line("      not the object's SCT values", logfile = LOG)
log_line("    - high values drawn last, so positives are not buried by zeros",
         logfile = LOG)
log_line("    - scale bar from each section's own scalefactors_json",
         logfile = LOG)
log_line("    - no titles or captions in the image; legends in",
         logfile = LOG)
log_line("      results/figure_legends_geo.md", logfile = LOG)
log_line("", logfile = LOG)
log_line("  STILL A JUDGEMENT CALL: which section leads. ", MAIN_SECTION,
         " was chosen", logfile = LOG)
log_line("  for having the most clearly delimited melanoma mass. Change",
         logfile = LOG)
log_line("  MAIN_SECTION at the top and rerun to see another.", logfile = LOG)

report_warnings(logfile = LOG)
write_session_info("33_visium_figure_journal")
log_line("\nFull log: ", LOG, logfile = NULL)
