# ==============================================================================
# 13_gse4587_and_fetch.R
#
# Part A: analyse GSE4587 (melanoma progression series, Affymetrix U133 Plus 2.0)
# Part B: download the supplementary files for the four RNA-seq series and
#         report their structure, so the next script can parse them correctly.
#         GSE98394 is included alongside the series Gustavo listed; the reason
#         is documented in Part B.
#
# GSE4587 NOTES
#   18 samples spanning normal skin, cultured melanocytes and keratinocytes,
#   benign and atypical nevi, melanoma in situ, VGP and MGP melanoma, lymph node
#   metastases, and a cultured MGP line. That is 1-3 samples per group, so this
#   is DESCRIPTIVE ONLY. No statistical test is run, and none should be: with
#   n = 2 per group any p value would be meaningless and inviting one would be
#   worse than reporting the values themselves. Gustavo asked for direction of
#   change and individual values, which is the right call.
#
#   RASSF9 has two probe sets (210335_at, 241470_x_at) and RASSF10 one
#   (238755_at). The two RASSF9 probes are reported SEPARATELY rather than
#   averaged: 241470_x_at carries the "_x_" suffix, meaning Affymetrix flagged
#   it as cross-hybridising with other transcripts, so it is less trustworthy
#   than 210335_at. If the two disagree, that disagreement is the result.
#
# Output: results/geo/gse4587_values.csv
#         figures/geo/gse4587_progression.{pdf,png}
#         data/raw/geo/<series>/   supplementary files for the RNA-seq sets
#         results/logs/13_gse4587.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(patchwork)
})
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))

LOG <- file.path(PATHS$logs, "13_gse4587.log")
cat("GSE4587 analysis and supplementary fetch —", format(Sys.time()), "\n\n",
    file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

GEO_DIR <- file.path(PATHS$raw, "geo")
RES_GEO <- file.path(PATHS$results, "geo")
FIG_GEO <- file.path(PATHS$figures, "geo")
dir.create(RES_GEO, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_GEO, recursive = TRUE, showWarnings = FALSE)

options(timeout = 1800)

# ==============================================================================
# PART A: GSE4587
# ==============================================================================
hdr("A1. Load GSE4587")

es <- getGEO("GSE4587", destdir = GEO_DIR, GSEMatrix = TRUE, getGPL = TRUE)[[1]]
ex <- exprs(es)
fd <- fData(es)
pd <- pData(es)

log_line("Samples: ", ncol(ex), " | probes: ", nrow(ex), logfile = LOG)
log_line("Value range: ", round(min(ex, na.rm = TRUE), 2), " to ",
         round(max(ex, na.rm = TRUE), 2), logfile = LOG)

# Values are on a linear intensity scale; log2 for display and comparison.
if (max(ex, na.rm = TRUE) > 100) {
  log_line("Linear intensities detected -> log2(x + 1) applied", logfile = LOG)
  ex <- log2(ex + 1)
} else {
  log_line("Values already appear log-scaled; left as provided", logfile = LOG)
}

# ==============================================================================
# A2. Probes and groups
# ==============================================================================
hdr("A2. Probes and groups")

PROBES <- c(RASSF9_210335_at   = "210335_at",
            RASSF9_241470_x_at = "241470_x_at",
            RASSF10_238755_at  = "238755_at")

for (nm in names(PROBES)) {
  p <- PROBES[[nm]]
  if (p %in% rownames(ex)) {
    log_line(sprintf("  [ ok ] %-20s present", nm), logfile = LOG)
  } else {
    log_line(sprintf("  [FLAG] %-20s NOT in matrix", nm), logfile = LOG)
  }
}
PROBES <- PROBES[PROBES %in% rownames(ex)]

grp_raw <- as.character(pd$source_name_ch1)
log_line("\n  groups present:", logfile = LOG)
for (g in unique(grp_raw)) {
  log_line(sprintf("    %-32s n = %d", g, sum(grp_raw == g)), logfile = LOG)
}

# Order along the progression axis. Keratinocytes are a lineage control, not a
# progression stage, and are kept at the end rather than dropped silently.
LEVELS <- c("normal skin", "melanocytes", "keratinocytes", "benign nevi",
            "atypical nevi", "melanoma in situ", "VGP melanoma",
            "MGP melanoma", "lymph node metastasis",
            "metastatic melanoma culture")

norm_group <- function(x) {
  x <- tolower(trimws(x))
  dplyr::case_when(
    grepl("normal skin", x)                 ~ "normal skin",
    grepl("melanocyte", x)                  ~ "melanocytes",
    grepl("keratinocyte", x)                ~ "keratinocytes",
    grepl("benign nev", x)                  ~ "benign nevi",
    grepl("atypical nev|dysplastic", x)     ~ "atypical nevi",
    grepl("in situ", x)                     ~ "melanoma in situ",
    grepl("vgp", x)                         ~ "VGP melanoma",
    grepl("mgp melanoma$|^mgp", x) & !grepl("culture", x) ~ "MGP melanoma",
    grepl("lymph node", x)                  ~ "lymph node metastasis",
    grepl("culture", x)                     ~ "metastatic melanoma culture",
    TRUE                                    ~ x
  )
}

dat <- tibble::tibble(
  sample = colnames(ex),
  title  = as.character(pd$title),
  group  = factor(norm_group(grp_raw), levels = LEVELS)
)

unmapped <- unique(grp_raw[is.na(dat$group)])
if (length(unmapped)) {
  log_line("\n  [FLAG] unmapped group labels: ",
           paste(unmapped, collapse = " | "), logfile = LOG)
}

for (nm in names(PROBES)) dat[[nm]] <- ex[PROBES[[nm]], ]

# ==============================================================================
# A3. Report values
# ==============================================================================
hdr("A3. Values by stage (log2 intensity)")

long <- dat |>
  tidyr::pivot_longer(dplyr::all_of(names(PROBES)),
                      names_to = "probe", values_to = "value") |>
  dplyr::mutate(gene = ifelse(grepl("^RASSF10", probe), "RASSF10", "RASSF9"))

for (p in names(PROBES)) {
  log_line("\n  ", p, ":", logfile = LOG)
  sub <- long |> dplyr::filter(probe == p) |> dplyr::arrange(group)
  for (g in levels(dat$group)) {
    v <- sub$value[sub$group == g]
    if (!length(v)) next
    log_line(sprintf("    %-30s n=%d  %s", g, length(v),
                     paste(sprintf("%.2f", v), collapse = ", ")), logfile = LOG)
  }
}

log_line("\n  No statistical test is reported: with 1-3 samples per group a p",
         logfile = LOG)
log_line("  value would be uninterpretable. Direction and individual values",
         logfile = LOG)
log_line("  are what this dataset can support.", logfile = LOG)

readr::write_csv(long, file.path(RES_GEO, "gse4587_values.csv"))
log_line("\n  Wrote: results/geo/gse4587_values.csv", logfile = LOG)

# ==============================================================================
# A4. Figure
# ==============================================================================
hdr("A4. Figure")

panel <- function(p) {
  d <- long |> dplyr::filter(probe == p, !is.na(group))
  ggplot(d, aes(x = group, y = value)) +
    geom_point(size = 1.8, shape = 21, fill = PAL$hist, colour = PAL$ink,
               stroke = 0.4,
               position = position_jitter(width = 0.08, height = 0)) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.45, linewidth = 0.3,
                 colour = PAL$uv_high) +
    labs(x = NULL, y = expression(bold(log[2]~"intensity")), subtitle = p) +
    theme_pub(grid = TRUE) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
}

fig <- wrap_plots(lapply(names(PROBES), panel), ncol = 1) +
  plot_annotation(
    title = "GSE4587: melanocytic progression series",
    caption = paste("Descriptive only; 1-3 samples per group.",
                    "Red bar = group mean."),
    tag_levels = "A") &
  theme(plot.tag = element_text(size = FIG$base_size + 3, face = "bold"),
        plot.caption = element_text(size = FIG$base_size - 1,
                                    colour = PAL$grey_pt, hjust = 0))

for (fmt in c("pdf", "png")) {
  f <- file.path(FIG_GEO, paste0("gse4587_progression.", fmt))
  if (file.exists(f)) file.remove(f)
  args <- list(filename = f, plot = fig,
               width = 150, height = 190, units = "mm")
  if (fmt == "pdf") {
    args$device <- grDevices::cairo_pdf
  } else {
    args$dpi <- 600          # ggplot2 4.x rejects dpi = NA, so set it only here
    args$bg  <- "white"
  }
  do.call(ggplot2::ggsave, args)
}
log_line("  Wrote: figures/geo/gse4587_progression.{pdf,png}", logfile = LOG)

# ==============================================================================
# PART B: fetch supplementary files for the RNA-seq series
# ==============================================================================
hdr("B. Supplementary files for the RNA-seq series")

# GSE98394 is ADDED to the four series Gustavo listed. Rationale, for the record:
# GSE230574 comprises only two 10x runs, and one of them mixes primary
# melanocytes with four melanoma lines in a single sample, so no pseudobulk
# comparison with replication is possible there -- it can illustrate detection
# but cannot test a change along progression. GSE98394 (27 benign melanocytic
# lesions vs 51 primary melanomas, Illumina RNA-seq) asks the SAME question as
# GSE112509 in an INDEPENDENT cohort, which is what "reproducible change along
# transformation" actually requires. GSE230574 is retained as illustration.

RNASEQ_SERIES <- c("GSE112509", "GSE98394", "GSE109245", "GSE230574")

for (gse in RNASEQ_SERIES) {
  log_line("\n  ", gse, ":", logfile = LOG)

  # Series-matrix metadata first: confirms the group labels before any parsing.
  meta <- try(getGEO(gse, destdir = GEO_DIR, GSEMatrix = TRUE, getGPL = FALSE),
              silent = TRUE)
  if (!inherits(meta, "try-error") && length(meta)) {
    pdm <- pData(meta[[1]])
    log_line("    samples: ", nrow(pdm), logfile = LOG)
    gcols <- names(pdm)[grepl("characteristics|source_name|:ch1$",
                              names(pdm))]
    for (cc in head(gcols, 6)) {
      vals <- unique(as.character(pdm[[cc]]))
      if (length(vals) > 1 && length(vals) <= nrow(pdm)) {
        log_line(sprintf("    %-28s %d value(s): %s", cc, length(vals),
                         paste(substr(head(vals, 4), 1, 46), collapse = " | ")),
                 logfile = LOG)
      }
    }
  }

  out <- try(getGEOSuppFiles(gse, baseDir = GEO_DIR, makeDirectory = TRUE),
             silent = TRUE)
  if (inherits(out, "try-error")) {
    log_line("    [FAILED] ", as.character(out), logfile = LOG)
    next
  }
  files <- list.files(file.path(GEO_DIR, gse), full.names = TRUE)
  for (f in files) {
    log_line(sprintf("    %-58s %.1f MB", basename(f),
                     file.info(f)$size / 1024^2), logfile = LOG)

    # Peek inside text-like files so the next script knows the layout.
    if (grepl("\\.(txt|tsv|csv)(\\.gz)?$", f, ignore.case = TRUE)) {
      con <- try(if (grepl("\\.gz$", f)) gzfile(f, "rt") else file(f, "rt"),
                 silent = TRUE)
      if (!inherits(con, "try-error")) {
        head_lines <- try(readLines(con, n = 3), silent = TRUE)
        close(con)
        if (!inherits(head_lines, "try-error")) {
          for (hl in head_lines) {
            log_line("        | ", substr(hl, 1, 160), logfile = LOG)
          }
        }
      }
    }
  }
}

hdr("DONE")
log_line("GSE4587 is analysed. The supplementary files above are what the next",
         logfile = LOG)
log_line("script will parse for GSE112509, GSE98394, GSE109245 and GSE230574 --",
         logfile = LOG)
log_line("the first lines printed show whether they carry gene symbols or",
         logfile = LOG)
log_line("Ensembl IDs, and whether they are counts, TPM or something else.",
         logfile = LOG)
log_line("\nGSE98394 was added to the list: see the note in Part B for why.",
         logfile = LOG)

report_warnings(LOG)
write_session_info("13_gse4587")
log_line("\nFull log: ", LOG, logfile = LOG)
