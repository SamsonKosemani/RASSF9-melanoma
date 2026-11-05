# ==============================================================================
# 01_check_sbs.R
#
# Purpose: validate the SigProfiler SBS signature file BEFORE anything depends
# on it, then derive the UV-exposure proxies for TCGA-SKCM.
#
# This script makes no assumptions about exact column names. It detects them,
# reports what it found, and stops loudly if something is wrong. Read the
# console output — it is the point of this script, not a side effect.
#
# Input:  data/raw/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv
# Output: data/processed/skcm_uv_proxies.csv
#         results/logs/01_sbs_qc.log
#         figures/qc_01_sbs7_distribution.png
# ==============================================================================

source(here::here("R", "00_setup.R"))

LOG <- file.path(PATHS$logs, "01_sbs_qc.log")
cat("SBS signature QC —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ", strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

# ==============================================================================
# 1. Locate the file
# ==============================================================================
hdr("1. Input file")

# The file has shipped under slightly different names across releases
# ("signature" vs "signatures", with or without "WES"). Match on a pattern
# rather than an exact name, and search the usual download locations.
sbs_pattern <- "sigProfiler.*SBS.*signature.*samples.*\\.csv$"

find_sbs <- function(dirs, pattern) {
  hits <- unlist(lapply(dirs, function(d) {
    if (!dir.exists(d)) return(character(0))
    list.files(d, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  }))
  hits[file.exists(hits)]
}

search_dirs <- c(
  PATHS$raw,
  file.path(path.expand("~"), "Downloads"),
  file.path(path.expand("~"), "Desktop"),
  file.path(path.expand("~"), "Documents")
)

hits <- find_sbs(search_dirs, sbs_pattern)

if (length(hits) == 0) {
  stop("Could not find the SigProfiler SBS signatures file.\n",
       "Searched (pattern '", sbs_pattern, "'):\n  ",
       paste(search_dirs, collapse = "\n  "), "\n",
       "Copy the CSV into data/raw/ and re-run.\n",
       "Source: https://www.synapse.org/#!Synapse:syn11726620 ",
       "(PCAWG/TCGA SigProfiler SBS signatures in samples)")
}

if (length(hits) > 1) {
  log_line("NOTE: more than one candidate file found:", logfile = LOG)
  for (h in hits) log_line("    ", h, logfile = LOG)
  log_line("Using the first. If that is wrong, move the correct file into ",
           "data/raw/ and delete the others from the search path.", logfile = LOG)
}

sbs_file <- hits[1]

# Keep the project self-contained: copy into data/raw/ if it came from elsewhere.
if (normalizePath(dirname(sbs_file)) != normalizePath(PATHS$raw)) {
  dest <- file.path(PATHS$raw, basename(sbs_file))
  log_line("Found outside the project at: ", sbs_file, logfile = LOG)
  log_line("Copying into data/raw/ so the analysis is self-contained.", logfile = LOG)
  file.copy(sbs_file, dest, overwrite = TRUE)
  sbs_file <- dest
}

log_line("File:  ", sbs_file, logfile = LOG)
log_line("Size:  ", round(file.info(sbs_file)$size / 1024^2, 2), " MB", logfile = LOG)
log_line("Mtime: ", format(file.info(sbs_file)$mtime), logfile = LOG)
log_line("MD5:   ", tools::md5sum(sbs_file), logfile = LOG)

# ==============================================================================
# 2. Read and inspect structure
# ==============================================================================
hdr("2. File structure")

sbs_raw <- readr::read_csv(sbs_file)

log_line("Dimensions: ", nrow(sbs_raw), " rows x ", ncol(sbs_raw), " columns",
         logfile = LOG)
log_line("First 6 column names: ",
         paste(head(names(sbs_raw), 6), collapse = " | "), logfile = LOG)
log_line("Last 3 column names:  ",
         paste(tail(names(sbs_raw), 3), collapse = " | "), logfile = LOG)

# --- Detect the key columns rather than hard-coding them ----------------------
nm <- names(sbs_raw)

col_cancer <- nm[str_detect(tolower(nm), "cancer|tumou?r.?type|project|cohort")][1]
col_sample <- nm[str_detect(tolower(nm), "sample")][1]
col_acc    <- nm[str_detect(tolower(nm), "accuracy")][1]
cols_sbs   <- nm[str_detect(nm, "^SBS")]

if (is.na(col_cancer)) stop("Could not identify a cancer-type column. Column names: ",
                            paste(nm, collapse = ", "))
if (is.na(col_sample)) stop("Could not identify a sample-name column. Column names: ",
                            paste(nm, collapse = ", "))
if (length(cols_sbs) == 0) stop("No columns matching '^SBS' found. This may not be ",
                                "the signatures-in-samples file.")

log_line("Detected cancer-type column: ", col_cancer, logfile = LOG)
log_line("Detected sample column:      ", col_sample, logfile = LOG)
log_line("Detected accuracy column:    ",
         ifelse(is.na(col_acc), "(none)", col_acc), logfile = LOG)
log_line("Detected ", length(cols_sbs), " SBS signature columns: ",
         paste(head(cols_sbs, 5), collapse = ", "), " ... ",
         paste(tail(cols_sbs, 3), collapse = ", "), logfile = LOG)

# --- The UV signatures must be present ----------------------------------------
sbs7 <- intersect(c("SBS7a", "SBS7b", "SBS7c", "SBS7d"), cols_sbs)
if (length(sbs7) == 0) {
  stop("No SBS7 columns found. Cannot derive a UV proxy from this file.")
}
if (length(sbs7) < 4) {
  log_line("WARNING: only ", length(sbs7), " of 4 SBS7 subtypes present (",
           paste(sbs7, collapse = ", "), "). UV burden will be the sum of ",
           "those available.", logfile = LOG)
} else {
  log_line("All four UV signatures present: ", paste(sbs7, collapse = ", "),
           logfile = LOG)
}

# ==============================================================================
# 3. Cancer types and SKCM selection
# ==============================================================================
hdr("3. Cancer types")

type_counts <- sbs_raw |>
  count(.data[[col_cancer]], sort = TRUE, name = "n_samples")

log_line("Distinct cancer types: ", nrow(type_counts), logfile = LOG)
cat(capture.output(print(as.data.frame(type_counts), row.names = FALSE)),
    sep = "\n", file = LOG, append = TRUE)
print(head(type_counts, 10))

# Melanoma label differs between releases ("Skin-Melanoma", "SKCM", ...).
#
# CRITICAL: uveal (eye) melanoma must NOT be included. It is a biologically
# distinct disease -- driven by GNAQ/GNA11 rather than BRAF/NRAS, arising in
# tissue shielded from UV, and carrying essentially no SBS7. Including it would
# add a block of low-SBS7 tumours with unrelated expression biology at one end
# of the UV axis and could manufacture a spurious association. We therefore
# match melanoma labels, then explicitly drop ocular ones.
mel_all <- type_counts[[col_cancer]][
  str_detect(tolower(type_counts[[col_cancer]]), "melanoma|skcm")
]

ocular <- mel_all[str_detect(tolower(mel_all), "eye|uvea|uveal|ocular|conjunctiv|uvm")]
mel_label <- setdiff(mel_all, ocular)

if (length(mel_label) == 0) {
  stop("No cutaneous melanoma cancer type found. Available types printed above.")
}

if (length(ocular) > 0) {
  n_ocular <- sum(sbs_raw[[col_cancer]] %in% ocular)
  log_line("EXCLUDED (not UV-driven): ", paste(ocular, collapse = ", "),
           " -- ", n_ocular, " samples", logfile = LOG)
}

if (length(mel_label) > 1) {
  log_line("NOTE: more than one cutaneous melanoma label found (",
           paste(mel_label, collapse = ", "), "). All will be included. ",
           "Check that this is intended.", logfile = LOG)
}
log_line("INCLUDED cutaneous melanoma label(s): ",
         paste(mel_label, collapse = ", "), logfile = LOG)

skcm <- sbs_raw |> filter(.data[[col_cancer]] %in% mel_label)
log_line("Melanoma samples in file: ", nrow(skcm), logfile = LOG)

# ==============================================================================
# 4. Integrity checks
# ==============================================================================
hdr("4. Integrity checks")

problems <- character(0)
flag <- function(cond, msg) {
  if (isTRUE(cond)) {
    problems <<- c(problems, msg)
    log_line("  [FLAG] ", msg, logfile = LOG)
  } else {
    log_line("  [ ok ] ", msg, logfile = LOG)
  }
}

# -- duplicates --
n_dup <- sum(duplicated(skcm[[col_sample]]))
flag(n_dup > 0, paste0("duplicate sample names: ", n_dup))

# -- missing values in SBS columns --
n_na <- sum(is.na(as.matrix(skcm[, cols_sbs])))
flag(n_na > 0, paste0("missing values across SBS columns: ", n_na))

# -- negative values (signature exposures must be >= 0) --
n_neg <- sum(as.matrix(skcm[, cols_sbs]) < 0, na.rm = TRUE)
flag(n_neg > 0, paste0("negative signature exposures: ", n_neg))

# -- are values counts or proportions? --
row_sums <- rowSums(skcm[, cols_sbs], na.rm = TRUE)
looks_like_proportions <- stats::median(row_sums, na.rm = TRUE) < 1.5
log_line("  Row sums: median = ", round(stats::median(row_sums), 1),
         ", range = ", round(min(row_sums), 1), " to ", round(max(row_sums), 1),
         logfile = LOG)
if (looks_like_proportions) {
  log_line("  [FLAG] Row sums ~1 — these look like PROPORTIONS, not counts. ",
           "The analysis assumes absolute mutation counts. Check the file ",
           "version before proceeding.", logfile = LOG)
  problems <- c(problems, "values appear to be proportions, not counts")
} else {
  log_line("  [ ok ] Values look like absolute mutation counts.", logfile = LOG)
}

# -- barcode format --
bc_ok <- str_detect(skcm[[col_sample]], "^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}")
flag(any(!bc_ok), paste0("sample names not matching TCGA barcode format: ",
                         sum(!bc_ok)))
log_line("  Example sample names: ",
         paste(head(skcm[[col_sample]], 3), collapse = " | "), logfile = LOG)
log_line("  Barcode field count (first row): ",
         length(strsplit(skcm[[col_sample]][1], "-")[[1]]),
         " — 3 means patient-level only, no sample-type code.", logfile = LOG)

# -- accuracy of the signature fit --
if (!is.na(col_acc)) {
  acc <- skcm[[col_acc]]
  log_line("  Reconstruction accuracy: median = ", round(median(acc, na.rm = TRUE), 3),
           ", min = ", round(min(acc, na.rm = TRUE), 3), logfile = LOG)
  flag(sum(acc < 0.80, na.rm = TRUE) > 0,
       paste0("samples with reconstruction accuracy < 0.80: ",
              sum(acc < 0.80, na.rm = TRUE),
              " (consider excluding in a sensitivity analysis)"))
}

# ==============================================================================
# 5. Derive UV proxies
# ==============================================================================
hdr("5. UV proxy derivation")

uv <- skcm |>
  transmute(
    sample_full   = .data[[col_sample]],
    patient       = tcga_patient(.data[[col_sample]]),
    sample_type   = tcga_sample_type(.data[[col_sample]]),
    accuracy      = if (!is.na(col_acc)) .data[[col_acc]] else NA_real_,
    SBS7_total    = rowSums(across(all_of(sbs7)), na.rm = TRUE),
    TMB_sig       = rowSums(across(all_of(cols_sbs)), na.rm = TRUE)
  ) |>
  mutate(
    SBS7_fraction = ifelse(TMB_sig > 0, SBS7_total / TMB_sig, NA_real_),
    log10_SBS7    = log10(SBS7_total + 1),
    log10_TMB     = log10(TMB_sig + 1)
  )

log_line("Samples with SBS7_total == 0: ", sum(uv$SBS7_total == 0), logfile = LOG)
log_line("Samples with SBS7_total  > 0: ", sum(uv$SBS7_total > 0), logfile = LOG)

# NOTE: whether to drop SBS7 == 0 samples is a real analytical decision.
# Dropping them restricts the cohort to UV-attributable melanomas (acral and
# mucosal cases are largely UV-independent). Keeping them widens the dynamic
# range of the proxy. We retain a flag here and decide in 02, so the choice is
# explicit and reversible rather than baked into the input file.
uv <- uv |> mutate(has_uv_signature = SBS7_total > 0)

desc <- function(x, label) {
  log_line(sprintf("  %-14s n=%d  mean=%.3g  sd=%.3g  median=%.3g  IQR=[%.3g, %.3g]  range=[%.3g, %.3g]",
                   label, sum(!is.na(x)), mean(x, na.rm = TRUE), sd(x, na.rm = TRUE),
                   median(x, na.rm = TRUE),
                   quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE),
                   min(x, na.rm = TRUE), max(x, na.rm = TRUE)),
           logfile = LOG)
}

log_line("\nDescriptive statistics (all melanoma samples):", logfile = LOG)
desc(uv$SBS7_total,    "SBS7_total")
desc(uv$TMB_sig,       "TMB_sig")
desc(uv$SBS7_fraction, "SBS7_fraction")

log_line("\nDescriptive statistics (SBS7 > 0 only):", logfile = LOG)
uvp <- uv |> filter(has_uv_signature)
desc(uvp$SBS7_total,    "SBS7_total")
desc(uvp$SBS7_fraction, "SBS7_fraction")

# Tertile cutoffs, reported now so they are documented before use in 04.
tert <- quantile(uvp$SBS7_total, c(1/3, 2/3), na.rm = TRUE)
log_line("\nTertile cutoffs on SBS7_total (SBS7>0 samples): lower = ",
         round(tert[1]), ", upper = ", round(tert[2]), logfile = LOG)

# ==============================================================================
# 6. QC figure
# ==============================================================================
hdr("6. QC figure")

p <- ggplot(uvp, aes(x = log10_SBS7)) +
  geom_histogram(bins = 40, fill = PAL$neutral, colour = NA) +
  geom_vline(xintercept = log10(tert + 1), linetype = "dashed",
             colour = PAL$uv_high, linewidth = 0.3) +
  labs(x = expression(log[10]~"(SBS7 mutations + 1)"),
       y = "Tumours",
       title = "UV mutational burden, TCGA melanoma",
       subtitle = paste0("n = ", nrow(uvp),
                         "; dashed lines = tertile cutoffs")) +
  theme_pub()

save_fig(p, "qc_01_sbs7_distribution", width = FIG$width_single, height = 55,
         formats = "png")

# ==============================================================================
# 7. Write output
# ==============================================================================
hdr("7. Output")

out <- file.path(PATHS$processed, "skcm_uv_proxies.csv")
readr::write_csv(uv, out)
log_line("Wrote: ", out, " (", nrow(uv), " rows)", logfile = LOG)

write_session_info("01_check_sbs")

hdr("SUMMARY")
if (length(problems) == 0) {
  log_line("No integrity problems detected. Safe to proceed to 02.", logfile = LOG)
} else {
  log_line("Review these before proceeding to 02:", logfile = LOG)
  for (p_ in problems) log_line("  - ", p_, logfile = LOG)
}
log_line("\nFull log: ", LOG, logfile = LOG)
