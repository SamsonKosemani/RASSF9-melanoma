# ==============================================================================
# 03_build_cohort.R
#
# Purpose: build the analysis cohort by matching RASSF9/RASSF10 expression to
# UV mutational burden, resolving duplicate patients by a documented rule, and
# attaching tumour purity and clinical covariates.
#
# Implements decisions D4, D5, D7, D9, D10 from docs/decisions-log.md.
# No hypothesis is tested here. Expression values are assembled but NOT
# correlated with UV burden -- that happens in 04 and 05, after this cohort is
# fixed and inspected.
#
# Input:  data/raw/skcm_star_counts.rds
#         data/processed/skcm_uv_proxies.csv
# Output: data/processed/analysis_cohort.csv   (primary, SBS7 > 0)
#         data/processed/full_cohort.csv       (incl. SBS7 == 0, for sensitivity)
#         results/logs/03_cohort.log
# ==============================================================================

source(here::here("R", "00_setup.R"))

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(TCGAbiolinks)
})

# SummarizedExperiment attaches S4Vectors/BiocGenerics, which mask several
# dplyr verbs (count, rename, first, filter among them). Re-attaching dplyr
# here puts it back at the front of the search path. Without this, calls like
# count() dispatch to the Bioconductor generic and fail on data frames.
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))

LOG <- file.path(PATHS$logs, "03_cohort.log")
cat("Cohort construction —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

GENES <- c(RASSF9 = "RASSF9", RASSF10 = "RASSF10")

# ==============================================================================
# 1. Load inputs
# ==============================================================================
hdr("1. Inputs")

se <- readRDS(file.path(PATHS$raw, "skcm_star_counts.rds"))
uv <- readr::read_csv(file.path(PATHS$processed, "skcm_uv_proxies.csv"))

log_line("Expression object: ", nrow(se), " genes x ", ncol(se), " samples",
         logfile = LOG)
log_line("UV proxy table:    ", nrow(uv), " melanoma samples", logfile = LOG)

# ==============================================================================
# 2. Extract expression for the genes of interest
# ==============================================================================
hdr("2. Expression extraction")

rd  <- as.data.frame(rowData(se))
tpm <- assay(se, "tpm_unstrand")

expr_list <- list()
for (g in GENES) {
  idx <- which(rd$gene_name == g)
  if (length(idx) == 0) {
    log_line("  [FLAG] ", g, " not found -- skipped", logfile = LOG)
    next
  }
  if (length(idx) > 1) {
    log_line("  [FLAG] ", g, " maps to ", length(idx), " rows; summing TPM ",
             "across them", logfile = LOG)
    v <- colSums(tpm[idx, , drop = FALSE])
  } else {
    v <- tpm[idx, ]
  }
  expr_list[[g]] <- log2(v + 1)
  log_line("  [ ok ] ", g, " (", paste(rownames(se)[idx], collapse = ","),
           ") extracted; log2(TPM+1) median = ",
           round(median(log2(v + 1)), 3), logfile = LOG)
}

cd <- colData(se)

# colData from the GDC can contain list / CharacterList columns, which silently
# poison downstream dplyr verbs. Pull the two fields we need explicitly and
# force them to plain character vectors.
safe_chr <- function(x, fallback = NA_character_) {
  if (is.null(x)) return(fallback)
  if (is.list(x)) x <- vapply(x, function(e) if (length(e)) as.character(e)[1] else NA_character_,
                              character(1))
  as.character(x)
}

bc <- safe_chr(if ("barcode" %in% names(cd)) cd$barcode else colnames(se))
stv <- safe_chr(if ("sample_type" %in% names(cd)) cd$sample_type else NULL)

log_line("  barcode column class: ",
         paste(class(if ("barcode" %in% names(cd)) cd$barcode else colnames(se)),
               collapse = "/"),
         " -> coerced to character", logfile = LOG)
log_line("  example barcode: ", bc[1], logfile = LOG)

if (all(is.na(bc)) || any(nchar(bc) < 16, na.rm = TRUE)) {
  stop("Barcodes look malformed. First few: ", paste(head(bc, 3), collapse = " | "))
}

expr <- tibble(
  barcode     = bc,
  sample_type = stv
)

for (g in names(expr_list)) {
  expr[[paste0(g, "_log2tpm")]] <- as.numeric(unname(expr_list[[g]]))
}

expr <- expr |>
  mutate(
    patient   = tcga_patient(barcode),
    sample_id = substr(barcode, 1, 16)   # patient + sample type + vial
  )

stopifnot(is.character(expr$patient), is.character(expr$sample_id))
log_line("  expression table assembled: ", nrow(expr), " rows x ",
         ncol(expr), " columns", logfile = LOG)

# ==============================================================================
# 3. Resolve duplicate patients  (decision D9)
# ==============================================================================
hdr("3. Duplicate resolution (D9)")

# NOTE on D9 as originally written: an EXACT full-barcode match between the two
# datasets is impossible. The signature file uses DNA aliquots (…-11D-…) while
# expression uses RNA aliquots (…-11R-…), so the analyte and plate fields always
# differ. The correct level at which to match "the same tumour material" is the
# 16-character sample barcode (TCGA-XX-XXXX-NNV), which identifies participant,
# sample type and vial. That is what is used here. The intent of D9 is unchanged;
# the implementation is at sample rather than aliquot resolution.

uv_sample_ids <- unique(substr(uv$sample_full, 1, 16))

pat_tab <- table(expr$patient)
dup_patients <- names(pat_tab)[pat_tab > 1]
log_line("Patients with >1 expression sample: ", length(dup_patients), logfile = LOG)

if (length(dup_patients) > 0) {
  log_line("\nBefore resolution:", logfile = LOG)
  shown <- expr |>
    dplyr::filter(patient %in% dup_patients) |>
    dplyr::mutate(matches_sig_sample = sample_id %in% uv_sample_ids) |>
    dplyr::arrange(patient, barcode) |>
    dplyr::select(patient, barcode, sample_type, matches_sig_sample)
  cat(capture.output(print(as.data.frame(shown), row.names = FALSE)),
      sep = "\n", file = LOG, append = TRUE)
  print(as.data.frame(shown))
}

# Priority: (1) sample barcode matches the signature file, (2) primary tumour,
# (3) metastatic, (4) anything else. Ties broken by barcode for determinism.
expr_dedup <- expr |>
  dplyr::mutate(
    .match_sig = as.integer(sample_id %in% uv_sample_ids),
    .type_rank = dplyr::case_when(
      sample_type == "Primary Tumor"         ~ 3L,
      sample_type == "Metastatic"            ~ 2L,
      sample_type == "Additional Metastatic" ~ 1L,
      TRUE                                   ~ 0L
    )
  ) |>
  dplyr::arrange(patient, dplyr::desc(.match_sig), dplyr::desc(.type_rank), barcode) |>
  dplyr::group_by(patient) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(-.match_sig, -.type_rank)

log_line("\nExpression samples: ", nrow(expr), " -> ", nrow(expr_dedup),
         " after one-per-patient resolution", logfile = LOG)

# Drop the normal sample and any non-tumour material from the tumour analysis.
n_normal <- sum(expr_dedup$sample_type == "Solid Tissue Normal")
if (n_normal > 0) {
  log_line("Removing ", n_normal, " 'Solid Tissue Normal' sample(s) from the ",
           "tumour cohort (too few for a tumour-vs-normal test; see D8).",
           logfile = LOG)
  expr_dedup <- expr_dedup |> dplyr::filter(sample_type != "Solid Tissue Normal")
}

# ==============================================================================
# 4. Resolve duplicates on the UV side too
# ==============================================================================
hdr("4. UV table, one row per patient")

uv_pat_tab <- table(uv$patient)
uv_dup <- names(uv_pat_tab)[uv_pat_tab > 1]
log_line("Patients with >1 signature sample: ", length(uv_dup), logfile = LOG)

uv_dedup <- uv |>
  dplyr::arrange(patient, dplyr::desc(accuracy), sample_full) |>
  dplyr::group_by(patient) |>
  dplyr::slice(1) |>
  dplyr::ungroup()

if (length(uv_dup) > 0) {
  log_line("Resolved by keeping the highest reconstruction accuracy.",
           logfile = LOG)
}
log_line("UV rows: ", nrow(uv), " -> ", nrow(uv_dedup), logfile = LOG)

# ==============================================================================
# 5. Merge
# ==============================================================================
hdr("5. Merge")

cohort <- dplyr::inner_join(
  uv_dedup |> dplyr::select(patient, sig_barcode = sample_full,
                            sig_sample_type = sample_type,
                            accuracy, SBS7_total, TMB_sig, SBS7_fraction,
                            log10_SBS7, log10_TMB, has_uv_signature),
  expr_dedup |> dplyr::select(patient, rna_barcode = barcode, sample_type,
                              dplyr::ends_with("_log2tpm")),
  by = "patient"
)

log_line("UV-only patients:         ", nrow(uv_dedup), logfile = LOG)
log_line("Expression-only patients: ", nrow(expr_dedup), logfile = LOG)
log_line("Matched (merged) cohort:  ", nrow(cohort), logfile = LOG)

# How often did the two datasets refer to the same tumour sample?
cohort <- cohort |>
  dplyr::mutate(same_sample = substr(sig_barcode, 1, 16) == substr(rna_barcode, 1, 16))
log_line("Pairs where DNA and RNA came from the same sample vial: ",
         sum(cohort$same_sample), " / ", nrow(cohort),
         " (", round(100 * mean(cohort$same_sample), 1), "%)", logfile = LOG)

# Do the two sources agree on sample type?
disagree <- sum(cohort$sig_sample_type != cohort$sample_type, na.rm = TRUE)
log_line("Sample-type disagreement between DNA and RNA records: ", disagree,
         logfile = LOG)

# ==============================================================================
# 6. Tumour purity  (decision D10)
# ==============================================================================
hdr("6. Tumour purity")

purity_ok <- FALSE
try({
  data("Tumor.purity", package = "TCGAbiolinks", envir = environment())
  tp <- as.data.frame(Tumor.purity)

  # Some releases store these as comma-decimal character columns.
  num_cols <- intersect(c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE"), names(tp))
  tp[num_cols] <- lapply(tp[num_cols], function(x)
    suppressWarnings(as.numeric(gsub(",", ".", as.character(x)))))

  tp <- tp |>
    dplyr::mutate(patient = tcga_patient(as.character(Sample.ID))) |>
    dplyr::group_by(patient) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(num_cols),
                                   ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  cohort <- dplyr::left_join(cohort, tp, by = "patient")
  purity_ok <- TRUE

  for (cl in num_cols) {
    log_line(sprintf("  %-9s available for %d / %d samples (median %.3f)",
                     cl, sum(!is.na(cohort[[cl]])), nrow(cohort),
                     median(cohort[[cl]], na.rm = TRUE)), logfile = LOG)
  }
  log_line("\nPrimary purity measure: ABSOLUTE (per D10). CPE is a consensus ",
           "estimate and is retained as a fallback where ABSOLUTE is missing.",
           logfile = LOG)
}, silent = TRUE)

if (!purity_ok) {
  log_line("  [FLAG] Could not load Tumor.purity from TCGAbiolinks. ",
           "Purity will need to be sourced separately before the adjusted ",
           "model in script 06.", logfile = LOG)
}

# ==============================================================================
# 7. Apply exclusions and describe  (decisions D4, D5)
# ==============================================================================
hdr("7. Cohort definition")

full_cohort <- cohort                                          # D5 sensitivity set
primary_cohort <- cohort |> dplyr::filter(has_uv_signature)    # D5 primary set

log_line("Full cohort (incl. SBS7 = 0):  n = ", nrow(full_cohort), logfile = LOG)
log_line("Primary cohort (SBS7 > 0):     n = ", nrow(primary_cohort), logfile = LOG)
log_line("Excluded for SBS7 = 0:         n = ",
         nrow(full_cohort) - nrow(primary_cohort), logfile = LOG)
log_line("Low signature accuracy (<0.80) within primary cohort: ",
         sum(primary_cohort$accuracy < 0.80, na.rm = TRUE),
         "  (D6 sensitivity analysis)", logfile = LOG)

st <- as.data.frame(table(primary_cohort$sample_type))
names(st) <- c("sample_type", "n")
st$pct <- round(100 * st$n / sum(st$n), 1)
log_line("\nSample types in primary cohort:", logfile = LOG)
cat(capture.output(print(st, row.names = FALSE)),
    sep = "\n", file = LOG, append = TRUE)
print(st)

log_line("\nUV burden in primary cohort:", logfile = LOG)
log_line("  SBS7_total  median [IQR] = ",
         fmt_median_iqr(primary_cohort$SBS7_total, 0),
         "  range = [", min(primary_cohort$SBS7_total), ", ",
         max(primary_cohort$SBS7_total), "]", logfile = LOG)

tert <- quantile(primary_cohort$SBS7_total, c(1/3, 2/3), na.rm = TRUE)
log_line("  Tertile cutoffs: lower = ", round(tert[1]),
         ", upper = ", round(tert[2]),
         "   <-- these are the cutoffs script 04 will use", logfile = LOG)

# Expression descriptives. Reported so the cohort can be checked, NOT broken
# down by UV group -- that comparison belongs in script 04.
for (g in names(GENES)) {
  col <- paste0(g, "_log2tpm")
  if (col %in% names(primary_cohort)) {
    x <- primary_cohort[[col]]
    log_line(sprintf("\n  %s log2(TPM+1): median [IQR] = %s, range = [%.2f, %.2f], %%zero = %.1f",
                     g, fmt_median_iqr(x, 2), min(x), max(x),
                     100 * mean(x == 0)), logfile = LOG)
  }
}

# ==============================================================================
# 8. Write out
# ==============================================================================
hdr("8. Output")

readr::write_csv(primary_cohort, file.path(PATHS$processed, "analysis_cohort.csv"))
readr::write_csv(full_cohort,    file.path(PATHS$processed, "full_cohort.csv"))

log_line("Wrote: data/processed/analysis_cohort.csv (n = ",
         nrow(primary_cohort), ")", logfile = LOG)
log_line("Wrote: data/processed/full_cohort.csv (n = ",
         nrow(full_cohort), ")", logfile = LOG)

write_session_info("03_cohort")

hdr("DONE")
log_line("Add the final cohort numbers above to docs/decisions-log.md ",
         "under 'Still to document'.", logfile = LOG)
log_line("Full log: ", LOG, logfile = LOG)
