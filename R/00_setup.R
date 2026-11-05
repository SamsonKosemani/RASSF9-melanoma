# ==============================================================================
# 00_setup.R
# Shared configuration for the RASSF9 / TCGA-SKCM analysis.
# Sourced at the top of every other script. Produces no outputs of its own.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(janitor)
})

# ---- Reproducibility ---------------------------------------------------------
set.seed(123)
options(
  stringsAsFactors = FALSE,
  dplyr.summarise.inform = FALSE,
  readr.show_col_types = FALSE
)

# ---- Project paths -----------------------------------------------------------
PATHS <- list(
  raw        = here("data", "raw"),
  processed  = here("data", "processed"),
  figures    = here("figures"),
  results    = here("results"),
  logs       = here("results", "logs"),
  docs       = here("docs")
)

invisible(lapply(PATHS, dir.create, recursive = TRUE, showWarnings = FALSE))

# ---- Figure specification ----------------------------------------------------
# Cell Reports / Nature column widths. Multi-panel figures use the full double
# column so each panel keeps a legible size -- cramming four panels into a
# single column is what makes text look small.
FIG <- list(
  width_single = 90,
  width_1_5    = 140,
  width_double = 183,
  units        = "mm",
  base_size    = 9,      # 7 pt is the journal floor; 9 reads far better
  base_family  = "sans", # Arial on Windows maps through "sans" reliably
  dpi          = 600
)

# ---- Colour palette ----------------------------------------------------------
# Chosen to be distinguishable in greyscale and for common colour-vision
# deficiencies: a cool teal against a warm coral, well separated in both hue
# and lightness.
# Validated with the six-check palette validator (lightness band, chroma floor,
# CVD separation, normal-vision floor, contrast vs surface). The bright
# azure/orange pair passes ALL FIVE, with protanope dE 24.7 and contrast above
# 3:1 against the page -- brighter hues usually lose the contrast check, so
# these two were chosen as the most vivid pair that still clears it.
# For reference, the teal/coral pair used earlier FAILED the chroma floor and
# scored dE 7.1 under protanopia: two near-identical greys for a red-green
# colourblind reader, which is roughly 8% of male readers.
PAL <- list(
  uv_low   = "#0095D9",   # bright azure  (low / benign / reference group)
  uv_high  = "#F2680C",   # bright orange (high / malignant / contrast group)
  mid      = "#C9CBC8",   # muted grey for groups shown but not analysed
  hist     = "#0095D9",   # single-series distributions
  accent   = "#F2B705",   # bright gold, highlights only -- never a third series
  ink      = "#1A1A1A",   # text and axes
  grey_pt  = "#5A5A5A",
  tumor    = "#F2680C",
  normal   = "#0095D9",
  nevus    = "#0095D9",
  melanoma = "#F2680C"
)

PAL_UV <- c("UV-low" = PAL$uv_low, "Middle" = PAL$mid, "UV-high" = PAL$uv_high)
PAL_LESION <- c("Nevus" = PAL$nevus, "Melanoma" = PAL$melanoma)

# Three biological classes in the single-cell work. The violet was selected by
# running every candidate through protanope, deuteranope and tritanope
# simulation with azure and orange held fixed: #7E3FF2 scored the HIGHEST
# worst-case separation of any colour tested (dE 35.8) while also being the
# most saturated (C* 102), and it is darker than both others (L* 44.8 against
# 58.5 and 60.6) so the three separate in greyscale too.
PAL_CLASS <- c(
  "Nevus (primary)"       = "#0095D9",
  "Melanocyte (cultured)" = "#7E3FF2",
  "Melanoma line"         = "#F2680C"
)

# Five lineages for the marker-based cell assignment, chosen over four
# alternative sets by the same simulation: worst-case CVD dE 23.1, minimum
# chroma 40. For comparison, a hand-picked set of teals and blues scored 4.5 --
# two colours a deuteranope could not separate at all. Melanocytic takes the
# warm focal colour; greys carry the non-categories so they recede.
PAL_LINEAGE <- c(
  "Melanocytic"               = "#F2680C",
  "Fibroblast"                = "#0095D9",
  "Keratinocyte"              = "#7E3FF2",
  "Immune"                    = "#1F3E7C",
  "Endothelial"               = "#A33F00",
  "Ambiguous"                 = "#8E8E8E",
  "Undetermined (no markers)" = "#C9CBC8",
  "Undetermined (shallow)"    = "#E4E4E2"
)

# Sequential ramp for ORDERED categories (progression stages). Magnitude and
# order are carried by lightness, so it stays readable in greyscale and under
# any colour-vision deficiency -- a categorical palette would imply the stages
# are unrelated identities, which they are not.
pal_sequential <- function(n) {
  grDevices::colorRampPalette(c("#DFF1FB", "#9AD9F5", "#4CBCEC",
                                "#0095D9", "#0B5C94"))(n)
}

# ---- Publication theme -------------------------------------------------------
# Bolder than ggplot defaults: heavier axis lines, bold axis titles, generous
# tick labels. Designed to stay legible when a figure is reduced to column
# width in print.
theme_pub <- function(base_size = FIG$base_size,
                      base_family = FIG$base_family,
                      grid = FALSE) {
  th <- theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line         = element_line(linewidth = 0.45, colour = PAL$ink),
      axis.ticks        = element_line(linewidth = 0.45, colour = PAL$ink),
      axis.ticks.length = unit(2.2, "pt"),
      axis.text         = element_text(size = base_size, colour = PAL$ink),
      axis.title        = element_text(size = base_size + 0.5, colour = PAL$ink,
                                       face = "bold"),
      axis.title.x      = element_text(margin = margin(t = 4)),
      axis.title.y      = element_text(margin = margin(r = 4)),
      plot.title        = element_text(size = base_size + 1, face = "bold",
                                       hjust = 0, margin = margin(b = 3)),
      plot.subtitle     = element_text(size = base_size - 0.5,
                                       colour = PAL$grey_pt,
                                       margin = margin(b = 4)),
      plot.tag          = element_text(size = base_size + 3, face = "bold",
                                       colour = PAL$ink),
      legend.title      = element_text(size = base_size, face = "bold"),
      legend.text       = element_text(size = base_size),
      legend.key.size   = unit(9, "pt"),
      legend.background = element_blank(),
      legend.margin     = margin(0, 0, 0, 0),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, face = "bold",
                                       colour = PAL$ink),
      plot.margin       = margin(4, 5, 4, 4, unit = "pt")
    )
  if (grid) {
    th <- th + theme(panel.grid.major.y = element_line(linewidth = 0.25,
                                                       colour = "grey90"))
  }
  th
}

theme_set(theme_pub())
update_geom_defaults("point", list(size = 0.7))
update_geom_defaults("line",  list(linewidth = 0.45))

# ---- Figure housekeeping -----------------------------------------------------

# Remove previously generated versions of a figure so the folder never holds a
# stale PNG next to a fresh PDF.
purge_fig <- function(name, formats = c("pdf", "svg", "png", "tiff", "eps")) {
  files <- file.path(PATHS$figures, paste0(name, ".", formats))
  gone  <- files[file.exists(files)]
  if (length(gone)) {
    file.remove(gone)
    message("  removed stale figure files: ",
            paste(basename(gone), collapse = ", "))
  }
  invisible(gone)
}

# Wipe every figure matching a prefix (e.g. "figure1") before a rebuild.
purge_fig_prefix <- function(prefix) {
  files <- list.files(PATHS$figures, pattern = paste0("^", prefix),
                      full.names = TRUE)
  if (length(files)) {
    file.remove(files)
    message("  removed ", length(files), " stale file(s) matching '", prefix, "'")
  }
  invisible(files)
}

# ---- Journal mode ------------------------------------------------------------
# Nature / Cell / Science figures carry ONLY panel letters, axis labels with
# units, and legend keys. Titles, subtitles and explanatory captions do not
# belong in the image: every word of explanation goes in the figure legend in
# the manuscript text, which the journal typesets itself. A figure with prose
# baked in has to be rebuilt at submission, and reads as a slide rather than a
# figure.
#
# theme_journal() therefore blanks those elements outright, so a stray
# subtitle= left in a plot call cannot reach the page. Legends are written to a
# separate markdown file by add_legend(), ready to paste into the manuscript.
theme_journal <- function(base_size = FIG$base_size,
                          base_family = FIG$base_family,
                          grid = FALSE) {
  theme_pub(base_size = base_size, base_family = base_family, grid = grid) +
    theme(plot.title    = element_blank(),
          plot.subtitle = element_blank(),
          plot.caption  = element_blank())
}

LEGEND_FILE <- here("results", "figure_legends_geo.md")

reset_legends <- function(header = NULL) {
  dir.create(dirname(LEGEND_FILE), recursive = TRUE, showWarnings = FALSE)
  cat(if (is.null(header)) "# Figure legends\n" else header,
      "\n", sep = "", file = LEGEND_FILE)
  invisible(LEGEND_FILE)
}

# One legend per figure, in manuscript form: bold lead-in, then the text. Kept
# beside the figures rather than inside them.
add_legend <- function(file_stem, lead, text) {
  if (!file.exists(LEGEND_FILE)) reset_legends()
  cat("\n**", lead, "** ", text, "\n\n`", file_stem,
      ".pdf` / `.png` / `.svg`\n", sep = "", file = LEGEND_FILE, append = TRUE)
  invisible(NULL)
}

# Save into an arbitrary directory (the GEO figures live in figures/geo) with
# the same three formats and no stale leftovers.
save_fig_to <- function(plot, name, dir, width = FIG$width_double,
                        height = 120, formats = c("pdf", "svg", "png")) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (fmt in c("pdf", "svg", "png", "tiff", "eps")) {
    f <- file.path(dir, paste0(name, ".", fmt))
    if (file.exists(f)) file.remove(f)
  }
  for (fmt in formats) {
    file <- file.path(dir, paste0(name, ".", fmt))
    args <- list(filename = file, plot = plot,
                 width = width, height = height, units = FIG$units)
    if (fmt == "pdf") {
      args$device <- grDevices::cairo_pdf
    } else if (fmt == "svg") {
      if (requireNamespace("svglite", quietly = TRUE)) {
        args$device <- svglite::svglite
      } else next
    } else if (fmt == "png") {
      args$dpi <- FIG$dpi
      args$bg  <- "white"
    }
    do.call(ggplot2::ggsave, args)
  }
  invisible(plot)
}

# ---- Save helper -------------------------------------------------------------
# Writes PDF (submission), SVG (Illustrator editing) and PNG (drafts, slides)
# at identical dimensions. Old versions are deleted first.
save_fig <- function(plot, name, width = FIG$width_single, height = 60,
                     formats = c("pdf", "svg", "png")) {
  purge_fig(name)
  for (fmt in formats) {
    file <- file.path(PATHS$figures, paste0(name, ".", fmt))
    args <- list(filename = file, plot = plot,
                 width = width, height = height, units = FIG$units)
    if (fmt == "pdf") {
      args$device <- grDevices::cairo_pdf
    } else if (fmt == "svg") {
      if (requireNamespace("svglite", quietly = TRUE)) {
        args$device <- svglite::svglite
      } else {
        warning("svglite not installed; skipping SVG for ", name); next
      }
    } else if (fmt == "png") {
      args$dpi <- FIG$dpi
      args$bg  <- "white"
    }
    do.call(ggplot2::ggsave, args)
  }
  message("  saved figure: ", name, " (", width, " x ", height, " mm) [",
          paste(formats, collapse = ", "), "]")
  invisible(plot)
}

# ---- Small utilities ---------------------------------------------------------

tcga_patient <- function(x) substr(x, 1, 12)

# GDC capitalisation is used throughout so values match colData directly.
tcga_sample_type <- function(x) {
  code <- stringr::str_match(x, "^TCGA-[^-]+-[^-]+-(\\d{2})")[, 2]
  dplyr::case_when(
    code == "01" ~ "Primary Tumor",
    code == "06" ~ "Metastatic",
    code == "11" ~ "Solid Tissue Normal",
    is.na(code)  ~ NA_character_,
    TRUE         ~ paste0("Other (", code, ")")
  )
}

fmt_median_iqr <- function(x, digits = 2) {
  q <- stats::quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.*f [%.*f, %.*f]", digits, q[2], digits, q[1], digits, q[3])
}

fmt_p <- function(p) {
  ifelse(p < 2.2e-16, "< 2.2e-16",
         ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
                formatC(p, format = "f", digits = 3)))
}

# p value formatted as a plotmath expression for printing inside a panel.
fmt_p_plot <- function(p) {
  if (p < 0.0001) return("italic(P) < 0.0001")
  sprintf("italic(P) == %s", formatC(p, format = "f",
                                     digits = ifelse(p < 0.01, 4, 3)))
}

log_line <- function(..., logfile = NULL) {
  msg <- paste0(...)
  message(msg)
  if (!is.null(logfile)) cat(msg, "\n", file = logfile, append = TRUE)
  invisible(msg)
}

write_session_info <- function(tag) {
  f <- file.path(PATHS$logs, paste0("sessionInfo_", tag, ".txt"))
  capture.output(sessionInfo(), file = f)
  invisible(f)
}

# Summarise accumulated warnings into the log, so "50 or more warnings" never
# goes uninvestigated.
report_warnings <- function(logfile = NULL, n = 8) {
  w <- warnings()
  if (length(w) == 0) {
    log_line("  no warnings raised", logfile = logfile)
    return(invisible(NULL))
  }
  msgs <- unique(names(w))
  log_line("  ", length(w), " warning(s), ", length(msgs), " unique:",
           logfile = logfile)
  for (m in head(msgs, n)) log_line("    - ", m, logfile = logfile)
  invisible(msgs)
}

message("00_setup.R loaded | R ", getRversion(),
        " | project root: ", here::here())
