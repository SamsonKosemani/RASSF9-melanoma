# ==============================================================================
# 22_visium_inventory.R
#
# GSE298774 (Kreuger et al., Cancer Res 2026;86(17):4400-4413): transcriptome-
# wide Visium spatial transcriptomics of 18 NEVUS-ASSOCIATED MELANOMAS, each
# section containing both a dermal nevus component and a melanoma component
# from the same patient.
#
# WHY THIS MATTERS. Scripts 14 and 20 left one question open: is RASSF9's fall
# from nevus to melanoma melanocyte-intrinsic, or tissue composition? Bulk
# cannot separate them. GSE230574 has nevus tissue but melanoma only as cell
# lines. GSE168013 has the right design but a targeted panel without RASSF9.
#
# This dataset is the design that settles it, and it is transcriptome-wide, so
# RASSF9 and RASSF10 should both be measurable:
#   - benign and malignant melanocytes in the SAME section, so they share the
#     stroma -- a difference between them cannot be a difference in stromal
#     content;
#   - PAIRED within patient, so between-patient variation drops out;
#   - one block, one fixation, one run per patient, so no batch confound;
#   - tissue rather than culture, so the near-zero baseline seen in cell lines
#     does not apply.
#
# THIS SCRIPT ONLY LOOKS. Visium deposits can run to many gigabytes across 18
# samples, and the analysis depends on what per-spot annotation the deposit
# carries. So: list the files and their sizes WITHOUT downloading them, pull
# only the small series matrix for sample metadata, and report. Nothing is
# analysed and nothing large is fetched until we have seen the inventory.
#
# NOTE ON THEIR METHOD, for the analysis that follows. Kreuger et al. assign
# each SPOT to a cell type by a 15-topic deconvolution (topic X15 = nevus,
# X3 = tumour 1, X4 = tumour 2) and then test with a zero-inflated negative
# binomial random-effects model (R package twosigma) with PATIENT as the random
# effect, adjusting for age, sex and Breslow thickness. That topic object lives
# in a local .RData file in their repository's working directory and is NOT part
# of the GEO deposit, so unless the deposit carries per-spot annotations we must
# derive our own spot classification. Their paper supplies a discriminator for
# that: PRAME is negative in the nevus component and positive in the melanoma
# component in all 18 cases.
#
# Output: results/geo/gse298774_file_inventory.csv
#         results/geo/gse298774_sample_metadata.csv
#         results/geo/gse298774_clinical_tableS1.csv
#         results/logs/22_visium_inventory.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

LOG <- file.path(PATHS$logs, "22_visium_inventory.log")
cat("GSE298774 inventory —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

GEO_DIR <- file.path(PATHS$raw, "geo")
RES_GEO <- file.path(PATHS$results, "geo")
invisible(lapply(c(GEO_DIR, RES_GEO), dir.create, recursive = TRUE,
                 showWarnings = FALSE))

options(timeout = 3600)
ACC <- "GSE298774"

# ==============================================================================
# 1. What is in the deposit, without downloading it
# ==============================================================================
# fetch_files = FALSE returns the file list and sizes only. This is the whole
# point of running an inventory first: an 18-sample Visium series could be tens
# of gigabytes, and we do not yet know which parts we need.
hdr("1. File inventory (no download)")

inv <- try(GEOquery::getGEOSuppFiles(ACC, baseDir = GEO_DIR,
                                     makeDirectory = TRUE,
                                     fetch_files = FALSE), silent = TRUE)

if (inherits(inv, "try-error")) {
  log_line("  [FLAG] could not reach GEO: ", as.character(inv), logfile = LOG)
  log_line("  Check the network, then re-run. If GEO is reachable in a browser",
           logfile = LOG)
  log_line("  but not from R, the series page is at", logfile = LOG)
  log_line("  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", ACC,
           logfile = LOG)
} else {
  inv <- as.data.frame(inv)
  log_line("  files in the series-level supplement: ", nrow(inv), logfile = LOG)
  log_line("", logfile = LOG)
  cols <- intersect(c("fname", "size", "url"), names(inv))
  if (!length(cols)) {
    log_line("  columns returned: ", paste(names(inv), collapse = ", "),
             logfile = LOG)
    print(utils::head(inv, 40))
    utils::capture.output(print(utils::head(inv, 40)), file = LOG, append = TRUE)
  } else {
    show <- inv[, cols, drop = FALSE]
    if ("size" %in% cols) {
      show$size_MB <- round(suppressWarnings(as.numeric(show$size)) / 1024^2, 1)
    }
    for (i in seq_len(nrow(show))) {
      log_line(sprintf("    %-62s %10s MB", show$fname[i],
                       ifelse(is.null(show$size_MB), "?",
                              format(show$size_MB[i]))), logfile = LOG)
    }
    readr::write_csv(inv, file.path(RES_GEO, "gse298774_file_inventory.csv"))
    log_line("", logfile = LOG)
    log_line("  wrote: results/geo/gse298774_file_inventory.csv", logfile = LOG)

    tot <- sum(suppressWarnings(as.numeric(inv$size)), na.rm = TRUE) / 1024^3
    log_line(sprintf("  total series-level supplement: %.2f GB", tot),
             logfile = LOG)
  }

  # What to look for in that list, flagged so it is not missed:
  #   *_filtered_feature_bc_matrix.h5  or barcodes/features/matrix triplets
  #       -> the per-spot expression we need
  #   *tissue_positions*.csv           -> spot coordinates
  #   *annotation* / *metadata* / *clusters* / *deconvolution*
  #       -> per-spot labels; if present, we can use THEIR nevus/tumour calls
  #          instead of deriving our own, which would be far better
  pats <- c(matrix = "matrix|\\.h5$", positions = "tissue_positions|positions",
            annotation = "annot|metad|cluster|deconv|topic|label",
            image = "image|hires|lowres|scalefactors")
  log_line("", logfile = LOG)
  log_line("  what the names suggest is present:", logfile = LOG)
  for (nm in names(pats)) {
    hit <- grep(pats[[nm]], inv$fname, ignore.case = TRUE, value = TRUE)
    log_line(sprintf("    %-11s %d file(s)%s", nm, length(hit),
                     if (length(hit)) paste0(": ",
                       paste(utils::head(basename(hit), 3), collapse = ", "))
                     else ""), logfile = LOG)
  }
  if (!length(grep(pats[["annotation"]], inv$fname, ignore.case = TRUE))) {
    log_line("", logfile = LOG)
    log_line("  [NOTE] no obvious per-spot annotation file. Their nevus/tumour",
             logfile = LOG)
    log_line("  topic assignment is not in the deposit, so spot classification",
             logfile = LOG)
    log_line("  will have to be derived -- see the note at the top of this",
             logfile = LOG)
    log_line("  script. Check the per-SAMPLE supplements below too, since",
             logfile = LOG)
    log_line("  annotations are often deposited per GSM rather than per series.",
             logfile = LOG)
  }
}

# ==============================================================================
# 2. Sample-level metadata from the series matrix
# ==============================================================================
# Small, and it tells us the Visium IDs, which must match the paper's
# Supplementary Table 1, plus any per-sample characteristics and the per-GSM
# supplementary file names.
hdr("2. Sample metadata")

gse <- try(GEOquery::getGEO(ACC, destdir = GEO_DIR, GSEMatrix = TRUE,
                            getGPL = FALSE), silent = TRUE)
if (inherits(gse, "try-error")) {
  log_line("  [FLAG] series matrix not available: ", as.character(gse),
           logfile = LOG)
} else {
  pd <- Biobase::pData(if (is.list(gse)) gse[[1]] else gse)
  log_line("  samples: ", nrow(pd), " | metadata columns: ", ncol(pd),
           logfile = LOG)
  log_line("", logfile = LOG)
  keep <- grep("title|geo_accession|characteristics|source|supplementary",
               names(pd), value = TRUE, ignore.case = TRUE)
  for (cc in keep) {
    vals <- unique(as.character(pd[[cc]]))
    log_line(sprintf("    %-36s %3d distinct: %s", cc, length(vals),
                     substr(paste(utils::head(vals, 3), collapse = " | "),
                            1, 95)), logfile = LOG)
  }
  readr::write_csv(tibble::as_tibble(pd),
                   file.path(RES_GEO, "gse298774_sample_metadata.csv"))
  log_line("", logfile = LOG)
  log_line("  wrote: results/geo/gse298774_sample_metadata.csv", logfile = LOG)

  log_line("", logfile = LOG)
  log_line("  sample titles:", logfile = LOG)
  for (t in as.character(pd$title)) log_line("    ", t, logfile = LOG)
}

# ==============================================================================
# 3. The paper's clinical table, transcribed
# ==============================================================================
# Supplementary Table S1 of Kreuger et al., transcribed from the published
# table. These are the covariates their random-effects model used (age, sex,
# Breslow thickness), and they let the Visium IDs be matched to patients. Kept
# here so the analysis does not depend on reaching the journal again.
#
# Patients 5, 6, 9-12, 18, 19, 25-33 and 35-37 are absent from the published
# table; the 18 rows below are the lesions profiled by Visium. Every case is a
# WHO class I low-CSD melanoma and every one includes BOTH a melanoma in situ
# component AND a dermal nevus component.
hdr("3. Clinical table (Supplementary Table S1, transcribed)")

clin <- tibble::tribble(
  ~patient, ~visium_id, ~sex,      ~age, ~site,            ~breslow_mm, ~stage,
  1L,  "B1",  "Male",   48L, "mid back",       3.1, "T3aN0M0, Stage IIA",
  2L,  "A2",  "Female", 49L, "shoulder",       0.4, "T1aN0M0, Stage IA",
  3L,  "D3",  "Female", 41L, "face, jaw",      0.7, "T1aN0M0, Stage IA",
  4L,  "B4",  "Male",   72L, "upper back",     2.0, "T2aN0M0, Stage IB",
  7L,  "A7",  "Female", 31L, "lower back",     0.9, "T1bN0M0, Stage IA",
  8L,  "B8",  "Male",   62L, "mid back",       0.6, "T1aN0M0, Stage IA",
  13L, "B13", "Female", 38L, "face, temporal", 3.0, "T3bN3cM0, Stage IIIC",
  14L, "B14", "Female", 47L, "upper arm",      0.9, "T1bN0M0, Stage IA",
  15L, "B15", "Female", 53L, "back",           1.1, "T2aN0M0, Stage IB",
  16L, "A16", "Female", 72L, "shoulder",       1.2, "T2aN0M0, Stage IB",
  17L, "A17", "Female", 36L, "face, cheek",    1.1, "T2aN0M0, Stage IB",
  20L, "D20", "Male",   43L, "abdomen",        1.0, "T1bN0M0, Stage IA",
  21L, "A21", "Female", 40L, "lower back",     1.3, "T2aN0M0, Stage IB",
  22L, "A22", "Male",   48L, "lower back",     1.5, "T2aN0M0, Stage IB",
  23L, "A23", "Male",   37L, "shoulder",       0.9, "T1bN0M0, Stage IA",
  24L, "A24", "Male",   71L, "upper arm",      2.4, "T3aN0M0, Stage IIA",
  34L, "A34", "Female", 17L, "face, temporal", 0.9, "T1bN0M0, Stage IA",
  38L, "B38", "Male",   47L, "face, jaw",      0.7, "T1aN0M0, Stage IA"
)

# Anatomical site collapsed the way they did it, for use as a covariate.
clin <- clin |>
  dplyr::mutate(site_grouped = dplyr::case_when(
    grepl("back", site)  ~ "back",
    grepl("face", site)  ~ "face",
    TRUE                 ~ site))

log_line("  lesions in the published table: ", nrow(clin), logfile = LOG)
log_line("  sex: ", sum(clin$sex == "Female"), " female, ",
         sum(clin$sex == "Male"), " male", logfile = LOG)
log_line("  age: median ", median(clin$age), " (range ", min(clin$age), "-",
         max(clin$age), ")", logfile = LOG)
log_line("  Breslow: median ", median(clin$breslow_mm), " mm (range ",
         min(clin$breslow_mm), "-", max(clin$breslow_mm), ")", logfile = LOG)
log_line("  sites: ", paste(names(table(clin$site_grouped)), " n=",
                            as.integer(table(clin$site_grouped)),
                            sep = "", collapse = ", "), logfile = LOG)
log_line("  stages: ", paste(names(table(sub(".*Stage ", "Stage ", clin$stage))),
                             collapse = ", "), logfile = LOG)
log_line("", logfile = LOG)
log_line("  NOTE: one patient is 17 years old. If a per-patient covariate",
         logfile = LOG)
log_line("  model is fitted, age enters as a covariate as they did it; no",
         logfile = LOG)
log_line("  individual is identifiable from these aggregate fields.",
         logfile = LOG)

readr::write_csv(clin, file.path(RES_GEO, "gse298774_clinical_tableS1.csv"))
log_line("  wrote: results/geo/gse298774_clinical_tableS1.csv", logfile = LOG)

# ==============================================================================
# 4. What to decide from this
# ==============================================================================
hdr("4. Next step")

log_line("  Read section 1 and answer two questions:", logfile = LOG)
log_line("", logfile = LOG)
log_line("    (a) Is there per-spot ANNOTATION in the deposit (nevus vs", logfile = LOG)
log_line("        melanoma, or cluster/topic labels)? If yes, we use their", logfile = LOG)
log_line("        classification and the analysis is a paired comparison of", logfile = LOG)
log_line("        RASSF9 between annotated compartments within each patient.", logfile = LOG)
log_line("        If no, we derive spot classes from melanocytic markers plus", logfile = LOG)
log_line("        PRAME, which their paper reports as negative in the nevus", logfile = LOG)
log_line("        component and positive in the melanoma component in all 18", logfile = LOG)
log_line("        cases -- a weaker assignment, and it must be validated", logfile = LOG)
log_line("        before any RASSF9 number is read off it.", logfile = LOG)
log_line("", logfile = LOG)
log_line("    (b) How large is the download, and are the matrices per sample", logfile = LOG)
log_line("        or one combined file? That decides whether we pull all 18", logfile = LOG)
log_line("        samples or start with two or three.", logfile = LOG)
log_line("", logfile = LOG)
log_line("  NOTHING has been downloaded by this script beyond the series", logfile = LOG)
log_line("  matrix. The next script fetches only what these answers call for.",
         logfile = LOG)

report_warnings(logfile = LOG)
write_session_info("22_visium_inventory")
log_line("\nFull log: ", LOG, logfile = NULL)
