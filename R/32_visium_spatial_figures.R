# ==============================================================================
# 32_visium_spatial_figures.R
#
# The spatial panels, from the authors' own Seurat/semla object.
#
# WHY. Script 30 established that only 17.1% of RASSF9 molecules sit in
# melanocytic spots, with the highest expression in endothelium, eccrine duct
# and fibroblasts. That is the finding this figure SHOWS rather than asserts:
# PRAME filling the tumour mass while RASSF9 scatters across the whole section
# with no tumour enrichment. Roderick Slieker made exactly this plot himself
# and suggested the method, so it is the dataset author's own reading of his
# own data, reproduced.
#
# THE OBJECT'S STRUCTURE IS UNKNOWN UNTIL WE OPEN IT. A 1.64 GB RData can hold
# anything, under any object name, with any assay layout. So this script
# INVENTORIES FIRST and prints what it found before trying to use it. Every
# later section is wrapped so that a failure there still leaves the inventory
# in the log -- which is the thing needed to fix it.
#
# THREE JOBS, in descending order of value:
#   1. The spatial panels (section 5) -- the figure for the paper.
#   2. An independent check of the barcode join (section 4). Script 29 matched
#      39,264 of 39,264 spots and script 30's marker gate passed, so the join is
#      already well supported; reading RASSF9 straight from their object and
#      correlating it against our GEO-derived values would make it airtight.
#   3. A record of their preprocessing (section 3) -- which assays exist, and
#      whether raw counts are retained alongside SCT. We fitted LDA on raw
#      integers where they used SCT counts, a deviation recorded at the time;
#      this says how far apart the two actually are.
#
# NOTHING HERE CHANGES THE RESULT. Script 30 is the readout and stands on the
# GEO matrices. This is corroboration and illustration.
#
# Output: results/geo/gse298774_object_inventory.txt
#         results/geo/gse298774_join_crosscheck.csv
#         figures/geo/geo_spatial_<section>.{pdf,png}
#         results/logs/32_visium_spatial_figures.log
# ==============================================================================

source(here::here("R", "00_setup.R"))

LOG <- file.path(PATHS$logs, "32_visium_spatial_figures.log")
cat("GSE298774 spatial figures —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

RES_GEO <- file.path(PATHS$results, "geo")
FIG_GEO <- file.path(PATHS$figures, "geo")
PROC    <- PATHS$processed
invisible(lapply(c(RES_GEO, FIG_GEO), dir.create, recursive = TRUE,
                 showWarnings = FALSE))

OBJ_F <- file.path(PATHS$raw, "geo", "GSE298774", "Semla_Data_Deconvoluted.RData")
GENES <- c("RASSF9", "PRAME")        # the contrast that makes the point
EXTRA <- c("RASSF10", "MLANA")       # a SECOND figure, not crammed into the first

# Layout. The first run put four genes in a 10 x 3.2 inch canvas and the panels
# collided. Size is now derived from the panel count instead of fixed: each
# panel gets its own space and the canvas grows to fit.
PANEL_W <- 4.2      # inches per panel
PANEL_H <- 4.6      # inches per row
PT_SIZE <- 0.9      # spot size; semla's default is tuned for single sections

# Sections to draw: the patients with the largest paired compartments from
# script 29, so the nevus and tumour regions are both clearly visible.
WANT_SECTIONS <- c("A16", "B1", "B4", "A21")
N_FALLBACK    <- 4                   # if names cannot be matched, draw the first N

# ==============================================================================
# 1. Packages
# ==============================================================================
hdr("1. Packages")

need <- function(p, gh = NULL) {
  if (requireNamespace(p, quietly = TRUE)) return(TRUE)
  log_line("  installing: ", p, logfile = LOG)
  if (is.null(gh)) {
    try(utils::install.packages(p, repos = "https://cloud.r-project.org"),
        silent = TRUE)
  } else {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      utils::install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    try(remotes::install_github(gh, upgrade = "never"), silent = TRUE)
  }
  requireNamespace(p, quietly = TRUE)
}

has_seurat <- need("Seurat")
log_line("  Seurat: ", has_seurat, logfile = LOG)
if (has_seurat) {
  log_line("    version ", as.character(utils::packageVersion("Seurat")),
           logfile = LOG)
}

# semla's canonical repository is ludvigla/semla; Slieker's email pointed at
# spatial-research/semla. Both are tried.
has_semla <- requireNamespace("semla", quietly = TRUE)
if (!has_semla) has_semla <- need("semla", "spatial-research/semla")
if (!has_semla) has_semla <- need("semla", "ludvigla/semla")
log_line("  semla: ", has_semla, logfile = LOG)
if (has_semla) {
  log_line("    version ", as.character(utils::packageVersion("semla")),
           logfile = LOG)
}
if (!has_semla) {
  log_line("  semla could not be installed. The inventory and the join",
           logfile = LOG)
  log_line("  cross-check below still run; only the spatial panels are lost.",
           logfile = LOG)
  log_line("  Report whatever the install printed -- the message matters more",
           logfile = LOG)
  log_line("  than the failure.", logfile = LOG)
}

# ==============================================================================
# 2. Load, and see what is actually in there
# ==============================================================================
hdr("2. The file")

if (!file.exists(OBJ_F)) {
  stop("Object not found at ", OBJ_F)
}
log_line("  file: ", basename(OBJ_F), sprintf("  (%.2f GB)",
         file.size(OBJ_F) / 1024^3), logfile = LOG)
log_line("  loading -- a 1.6 GB RData takes a minute or two", logfile = LOG)

ENV <- new.env(parent = emptyenv())
t0 <- Sys.time()
loaded <- load(OBJ_F, envir = ENV)
log_line("  loaded in ", round(as.numeric(difftime(Sys.time(), t0,
                                                   units = "secs"))), "s",
         logfile = LOG)
log_line("  objects in the file: ", paste(loaded, collapse = ", "),
         logfile = LOG)
for (nm in loaded) {
  o <- get(nm, envir = ENV)
  log_line(sprintf("    %-20s %-28s %8.1f MB", nm,
                   paste(class(o), collapse = "/"),
                   as.numeric(utils::object.size(o)) / 1024^2), logfile = LOG)
}

# The Seurat object is whichever one carries assays; chosen by inspection
# rather than by assuming a name.
is_seurat <- vapply(loaded, function(nm)
  inherits(get(nm, envir = ENV), "Seurat"), logical(1))
if (!any(is_seurat)) {
  log_line("", logfile = LOG)
  log_line("  [STOP] no Seurat object among the loaded objects. The classes",
           logfile = LOG)
  log_line("  listed above are what the file actually holds -- the rest of",
           logfile = LOG)
  log_line("  this script assumes a Seurat object and cannot proceed.",
           logfile = LOG)
  report_warnings(logfile = LOG)
  stop("No Seurat object in ", basename(OBJ_F))
}
SE_NAME <- loaded[is_seurat][1]
se <- get(SE_NAME, envir = ENV)
log_line("  using object: ", SE_NAME, logfile = LOG)

# ==============================================================================
# 3. Inventory
# ==============================================================================
hdr("3. Inventory")

INV <- file.path(RES_GEO, "gse298774_object_inventory.txt")
sink(INV); print(se); cat("\n\n--- meta.data columns ---\n")
print(utils::str(se@meta.data, max.level = 1, list.len = 200)); sink()
log_line("  full print written to: results/geo/", basename(INV), logfile = LOG)

log_line("  spots: ", ncol(se), " | features: ", nrow(se), logfile = LOG)
log_line("  assays: ", paste(names(se@assays), collapse = ", "), logfile = LOG)
log_line("  default assay: ", Seurat::DefaultAssay(se), logfile = LOG)

for (a in names(se@assays)) {
  lay <- try(SeuratObject::Layers(se[[a]]), silent = TRUE)
  if (inherits(lay, "try-error")) lay <- slotNames(se@assays[[a]])
  log_line("    ", a, ": ", paste(lay, collapse = ", "), logfile = LOG)
}
log_line("  -> whether RAW COUNTS are retained alongside SCT matters for the",
         logfile = LOG)
log_line("     deviation recorded in script 26 (we fitted LDA on raw integers,",
         logfile = LOG)
log_line("     they used SCT counts)", logfile = LOG)

md <- se@meta.data
log_line("", logfile = LOG)
log_line("  meta.data: ", ncol(md), " columns", logfile = LOG)
key_cols <- intersect(c("library_id", "sample", "MaxTopicsName", "MaxTopics",
                        "nCount_Spatial", "nFeature_Spatial"), names(md))
log_line("  key columns present: ", paste(key_cols, collapse = ", "),
         logfile = LOG)

SEC_COL <- if ("library_id" %in% names(md)) "library_id" else key_cols[1]
secs <- sort(unique(as.character(md[[SEC_COL]])))
log_line("  sections (", length(secs), "): ", paste(secs, collapse = ", "),
         logfile = LOG)

if ("MaxTopicsName" %in% names(md)) {
  tb <- sort(table(md$MaxTopicsName), decreasing = TRUE)
  log_line("", logfile = LOG)
  log_line("  spot identities in the object:", logfile = LOG)
  for (i in seq_along(tb)) {
    log_line(sprintf("    %-24s %6d", names(tb)[i], as.integer(tb[i])),
             logfile = LOG)
  }
  log_line("  -> compare with script 29's counts from the text file; they",
           logfile = LOG)
  log_line("     should agree, since both come from this same fit", logfile = LOG)
}

# ==============================================================================
# 4. Independent check of the join
# ==============================================================================
hdr("4. Cross-check: their object against our GEO-derived values")

cross <- try({
  norm_sec <- function(x) gsub("[.\\-]", "_", as.character(x))
  core16   <- function(x) sub("-\\d+$", "", x)

  cn <- colnames(se)
  key_their <- paste0(norm_sec(md[[SEC_COL]]), "|", core16(cn))

  cand <- intersect(c("Spatial", "RNA"), names(se@assays))
  a <- if (length(cand)) cand[1] else names(se@assays)[1]
  log_line("  reading from assay: ", a, logfile = LOG)

  present <- intersect(GENES, rownames(se[[a]]))
  if (!length(present)) stop("neither gene is in assay ", a)

  m <- Seurat::GetAssayData(se, assay = a, layer = "counts")
  if (!nrow(m)) m <- Seurat::GetAssayData(se, assay = a, layer = "data")

  theirs <- tibble::tibble(key = key_their)
  for (g in present) theirs[[g]] <- as.numeric(m[g, ])

  # Script 29 already built a correct section + 16-mer key in `spot`; use it
  # rather than rebuilding one. The first version of this section rebuilt it
  # from `barcode`, which still carries the merge suffix, so it compared
  # "A16|AACAATGTGCTCCGAG" against "A16|AACAATGTGCTCCGAG-1" and matched zero of
  # 39,481 spots.
  LAB <- readRDS(file.path(PROC, "gse298774_kreuger_labels.rds"))$labels
  ours <- LAB |>
    dplyr::select(key = spot, label, patient)
  log_line("  our key example:   ", ours$key[1], logfile = LOG)
  log_line("  their key example: ", key_their[1], logfile = LOG)

  j <- dplyr::inner_join(theirs, ours, by = "key")
  log_line("  spots matched between their object and our label table: ",
           nrow(j), " of ", nrow(theirs), logfile = LOG)
  if (nrow(j) == 0) {
    stop("zero spots matched -- the two keys above are not in the same format")
  }
  readr::write_csv(j, file.path(RES_GEO, "gse298774_join_crosscheck.csv"))

  for (g in present) {
    byl <- j |>
      dplyr::group_by(label) |>
      dplyr::summarise(mean_count = mean(.data[[g]], na.rm = TRUE),
                       .groups = "drop") |>
      dplyr::arrange(dplyr::desc(mean_count))
    log_line("", logfile = LOG)
    log_line("  ", g, " mean raw count per spot, by their label:", logfile = LOG)
    for (i in seq_len(nrow(byl))) {
      log_line(sprintf("    %-24s %8.3f", byl$label[i], byl$mean_count[i]),
               logfile = LOG)
    }
  }
  log_line("", logfile = LOG)
  log_line("  If this ordering matches script 30 section 4 -- RASSF9 highest in",
           logfile = LOG)
  log_line("  endothelial and adnexal spots, lowest in Tumor 1 -- then the",
           logfile = LOG)
  log_line("  barcode join is confirmed from a completely independent route.",
           logfile = LOG)
  TRUE
}, silent = TRUE)
if (inherits(cross, "try-error")) {
  log_line("  cross-check could not be run: ", as.character(cross), logfile = LOG)
  log_line("  (the inventory above says why -- assay names or barcode format)",
           logfile = LOG)
}

# ==============================================================================
# 5. The spatial panels
# ==============================================================================
hdr("5. Spatial panels")

if (!has_semla) {
  log_line("  semla unavailable; panels skipped", logfile = LOG)
} else {
  suppressPackageStartupMessages(library(semla))

  # MapFeatures indexes sections by position, so position must be mapped back
  # to the section name rather than assumed.
  sec_order <- unique(as.character(md[[SEC_COL]]))
  log_line("  section index -> name:", logfile = LOG)
  for (i in seq_along(sec_order)) {
    log_line(sprintf("    %2d  %s", i, sec_order[i]), logfile = LOG)
  }

  idx <- which(sec_order %in% WANT_SECTIONS)
  if (!length(idx)) {
    idx <- seq_len(min(N_FALLBACK, length(sec_order)))
    log_line("  requested sections not found by name; drawing the first ",
             length(idx), logfile = LOG)
  }

  # TWO figures per section rather than one crowded one: the RASSF9/PRAME
  # contrast is the argument and deserves its own panel pair; RASSF10 and MLANA
  # are context and go in a second file.
  sets <- list(main = intersect(GENES, rownames(se)),
               extra = intersect(EXTRA, rownames(se)))
  sets <- sets[vapply(sets, length, integer(1)) > 0]
  log_line("  figure sets:", logfile = LOG)
  for (nm in names(sets)) {
    log_line("    ", nm, ": ", paste(sets[[nm]], collapse = ", "), logfile = LOG)
  }

  # semla's argument list varies by version, so the richer call is tried first
  # and a minimal one is the fallback. A plotting argument that does not exist
  # should not cost us the figure.
  map_try <- function(feats, i) {
    p <- try(semla::MapFeatures(se, features = feats, section = i,
                                pt_size = PT_SIZE, ncol = length(feats),
                                override_plot_dims = TRUE), silent = TRUE)
    if (!inherits(p, "try-error")) return(p)
    p <- try(semla::MapFeatures(se, features = feats, section = i,
                                pt_size = PT_SIZE), silent = TRUE)
    if (!inherits(p, "try-error")) return(p)
    try(semla::MapFeatures(se, features = feats, section = i), silent = TRUE)
  }

  for (i in idx) {
    nm <- sec_order[i]
    for (setnm in names(sets)) {
      feats <- sets[[setnm]]
      p <- map_try(feats, i)
      if (inherits(p, "try-error")) {
        log_line("  [FLAG] ", nm, " / ", setnm, " failed: ", as.character(p),
                 logfile = LOG)
        next
      }
      # Breathing room between panels and a legend that does not crowd them.
      p <- try(p & ggplot2::theme(
                 plot.margin = ggplot2::margin(6, 10, 6, 10, unit = "pt"),
                 legend.key.width = ggplot2::unit(10, "pt"),
                 legend.text = ggplot2::element_text(size = 7),
                 legend.title = ggplot2::element_text(size = 8),
                 plot.title = ggplot2::element_text(size = 10)),
               silent = TRUE)
      if (inherits(p, "try-error")) p <- map_try(feats, i)

      stem <- paste0("geo_spatial_", nm,
                     if (setnm == "extra") "_extra" else "")
      w <- PANEL_W * length(feats)
      ok <- try({
        ggplot2::ggsave(file.path(FIG_GEO, paste0(stem, ".png")), p,
                        width = w, height = PANEL_H, dpi = 400, bg = "white",
                        limitsize = FALSE)
        ggplot2::ggsave(file.path(FIG_GEO, paste0(stem, ".pdf")), p,
                        width = w, height = PANEL_H, bg = "white",
                        limitsize = FALSE)
        TRUE
      }, silent = TRUE)
      if (inherits(ok, "try-error")) {
        log_line("  [FLAG] ", stem, ": saving failed", logfile = LOG)
      } else {
        log_line(sprintf("  wrote: %-28s %.1f x %.1f in, %d panel(s)",
                         paste0(stem, ".{png,pdf}"), w, PANEL_H, length(feats)),
                 logfile = LOG)
      }
    }
  }

  log_line("", logfile = LOG)
  log_line("  These come out with semla's own titles and colour bars. They are",
           logfile = LOG)
  log_line("  the SCIENTIFIC draft, not the final figure: the journal-mode",
           logfile = LOG)
  log_line("  pass (no titles in the image, legend in the text) happens once",
           logfile = LOG)
  log_line("  the panel selection is settled, so that cosmetics are not redone",
           logfile = LOG)
  log_line("  for sections that get dropped.", logfile = LOG)
}

# ==============================================================================
# 6. What to look for
# ==============================================================================
hdr("6. Reading the panels")

log_line("  THE POINT OF THE FIGURE: PRAME marks the melanoma compartment and",
         logfile = LOG)
log_line("  should fill it densely. RASSF9 should NOT follow it -- if script",
         logfile = LOG)
log_line("  30's lineage result is right, RASSF9 appears scattered across the",
         logfile = LOG)
log_line("  whole section, including dermis and adnexal structures, with no",
         logfile = LOG)
log_line("  enrichment in the tumour.", logfile = LOG)
log_line("", logfile = LOG)
log_line("  WHAT WOULD CONTRADICT IT: RASSF9 tracking PRAME, or RASSF9", logfile = LOG)
log_line("  concentrated in the nevus region specifically. Either would mean",
         logfile = LOG)
log_line("  the pseudobulk summary is hiding spatial structure, and would need",
         logfile = LOG)
log_line("  explaining before the figure is used.", logfile = LOG)
log_line("", logfile = LOG)
log_line("  CAUTION ON READING SPARSE GENES SPATIALLY: RASSF9 is detected in",
         logfile = LOG)
log_line("  11% of spots. A scattered appearance is partly sparsity and not",
         logfile = LOG)
log_line("  only biology, so the figure ILLUSTRATES the quantified result in",
         logfile = LOG)
log_line("  script 30 section 4 -- it does not independently establish it, and",
         logfile = LOG)
log_line("  the legend must say so.", logfile = LOG)

report_warnings(logfile = LOG)
write_session_info("32_visium_spatial_figures")
log_line("\nFull log: ", LOG, logfile = NULL)
