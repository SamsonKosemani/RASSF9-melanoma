# ==============================================================================
# 02_download_expression.R
#
# Purpose: download TCGA-SKCM RNA-seq (STAR - Counts) from the GDC and save it
# as a SummarizedExperiment object. Download only -- no filtering, no analysis.
#
# This is slow (expect 15-45 min on a first run, several GB of traffic) but it
# only needs to happen once. The script detects an existing download and skips
# straight to reporting, so re-running it is cheap and safe.
#
# Output: data/raw/GDCdata/              (raw GDC files, left untouched)
#         data/raw/skcm_star_counts.rds  (SummarizedExperiment)
#         results/logs/02_download.log
# ==============================================================================

source(here::here("R", "00_setup.R"))

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

LOG <- file.path(PATHS$logs, "02_download.log")
cat("TCGA-SKCM expression download —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

SE_FILE  <- file.path(PATHS$raw, "skcm_star_counts.rds")
GDC_DIR  <- file.path(PATHS$raw, "GDCdata")

# ==============================================================================
# 1. Skip if already downloaded
# ==============================================================================
hdr("1. Checking for existing download")

if (file.exists(SE_FILE)) {
  log_line("Found existing file: ", SE_FILE, logfile = LOG)
  log_line("Size: ", round(file.info(SE_FILE)$size / 1024^2, 1), " MB",
           logfile = LOG)
  log_line("Loading it instead of re-downloading.", logfile = LOG)
  log_line("(Delete this file and re-run if you want a fresh download.)",
           logfile = LOG)
  se <- readRDS(SE_FILE)
} else {

  # ============================================================================
  # 2. Build the query
  # ============================================================================
  hdr("2. Querying the GDC")

  log_line("Project:       TCGA-SKCM", logfile = LOG)
  log_line("Data category: Transcriptome Profiling", logfile = LOG)
  log_line("Data type:     Gene Expression Quantification", logfile = LOG)
  log_line("Workflow:      STAR - Counts", logfile = LOG)

  query <- GDCquery(
    project      = "TCGA-SKCM",
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  n_files <- nrow(getResults(query))
  log_line("Files matched: ", n_files, logfile = LOG)
  if (n_files == 0) {
    stop("GDC query returned no files. Check internet access and that the GDC ",
         "is reachable (https://portal.gdc.cancer.gov).")
  }

  # ============================================================================
  # 3. Download
  # ============================================================================
  hdr("3. Downloading")
  log_line("This is the slow step. Progress is printed by TCGAbiolinks below.",
           logfile = LOG)

  GDCdownload(query, directory = GDC_DIR, method = "api", files.per.chunk = 25)

  # ============================================================================
  # 4. Assemble into a SummarizedExperiment
  # ============================================================================
  hdr("4. Preparing SummarizedExperiment")

  se <- GDCprepare(query, directory = GDC_DIR)
  saveRDS(se, SE_FILE)
  log_line("Saved: ", SE_FILE, logfile = LOG)
}

# ==============================================================================
# 5. Describe what we got
# ==============================================================================
hdr("5. Contents")

log_line("Genes:   ", nrow(se), logfile = LOG)
log_line("Samples: ", ncol(se), logfile = LOG)
log_line("Assays:  ", paste(assayNames(se), collapse = ", "), logfile = LOG)

if (!"tpm_unstrand" %in% assayNames(se)) {
  log_line("WARNING: 'tpm_unstrand' assay not present. Available assays are ",
           "listed above; script 03 will need adjusting.", logfile = LOG)
} else {
  log_line("[ ok ] tpm_unstrand assay present (used for expression).",
           logfile = LOG)
}

# --- sample types -------------------------------------------------------------
cd <- as.data.frame(colData(se))
if ("sample_type" %in% names(cd)) {
  st <- as.data.frame(table(cd$sample_type))
  names(st) <- c("sample_type", "n")
  log_line("\nSample types:", logfile = LOG)
  cat(capture.output(print(st, row.names = FALSE)),
      sep = "\n", file = LOG, append = TRUE)
  print(st)
}

# --- are our genes of interest present? --------------------------------------
rd <- as.data.frame(rowData(se))
gene_col <- intersect(c("gene_name", "external_gene_name", "symbol"), names(rd))[1]

if (is.na(gene_col)) {
  log_line("WARNING: could not find a gene-symbol column in rowData. ",
           "Columns are: ", paste(names(rd), collapse = ", "), logfile = LOG)
} else {
  log_line("\nGene symbol column: ", gene_col, logfile = LOG)
  for (g in c("RASSF9", "RASSF10")) {
    idx <- which(rd[[gene_col]] == g)
    if (length(idx) == 0) {
      log_line("  [FLAG] ", g, " NOT FOUND", logfile = LOG)
    } else {
      log_line("  [ ok ] ", g, " found (", length(idx), " row",
               ifelse(length(idx) > 1, "s", ""), "): ",
               paste(rownames(se)[idx], collapse = ", "), logfile = LOG)
    }
  }
}

# --- duplicate patients -------------------------------------------------------
if ("barcode" %in% names(cd)) {
  pat <- tcga_patient(cd$barcode)
  n_dup_pat <- sum(duplicated(pat))
  log_line("\nPatients with more than one RNA-seq sample: ", n_dup_pat,
           logfile = LOG)
  if (n_dup_pat > 0) {
    log_line("  These need a documented rule in script 03 (e.g. keep the ",
             "primary tumour, or the first barcode).", logfile = LOG)
  }
}

write_session_info("02_download")

hdr("DONE")
log_line("Next: script 03 builds the analysis cohort by matching this ",
         "expression data to data/processed/skcm_uv_proxies.csv", logfile = LOG)
log_line("Full log: ", LOG, logfile = LOG)
