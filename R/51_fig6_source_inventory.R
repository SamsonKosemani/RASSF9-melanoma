## =============================================================================
## 51_fig6_source_inventory.R
##
## Inventory only. Draws nothing, computes nothing, changes nothing.
##
## PURPOSE
##   Figure 6 is being rebuilt to the structure agreed on 3 Nov 2026:
##     A  RASSF9 expression, nevi vs primary melanoma, GSE112509 and GSE98394
##     B  Cohen's d with 95% CI for the same two cohorts
##     C  RASSF9 detection in malignant melanoma cells, GSE72056
##     D  Promoter methylation vs expression, TCGA-SKCM, n = 469
##     E  Promoter-core methylation, 56 melanoma lines vs 940 other lineages
##     F  Promoter-core methylation vs expression, 741 cancer cell lines
##
##   Panels E and F can be drawn today from 46_line_level_data.csv. Panels A
##   to D need per-sample data produced by earlier scripts in this project,
##   and the plotting code cannot be written without knowing which files hold
##   it and what the columns are called.
##
##   This script finds those files, prints their columns and a few rows, and
##   says plainly what is missing and which script produces it. Writing the
##   figure code against guessed column names is how the last two scripts
##   broke, so this step exists to make that impossible.
##
## OUTPUT
##   results/logs/51_fig6_source_inventory.log   (send this back)
## =============================================================================

if (file.exists("R/00_setup.R")) source("R/00_setup.R")
suppressPackageStartupMessages(library(data.table))

ROOT    <- if (exists("PROJECT_ROOT")) PROJECT_ROOT else getwd()
LOG_DIR <- file.path(ROOT, "results", "logs")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_PATH <- file.path(LOG_DIR, "51_fig6_source_inventory.log")
.logcon  <- file(LOG_PATH, open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = .logcon) }
on.exit({ say(""); say("Full log: ", LOG_PATH); close(.logcon) }, add = TRUE)

say("")
say("================================================================")
say(" Figure 6 source inventory")
say(" project root: ", ROOT)
say("================================================================")

## ---- helpers ----------------------------------------------------------------
## Search the whole results/ and data/processed trees rather than assuming a
## layout. Returns every match so that duplicates are visible rather than
## silently resolved to the first one.
search_all <- function(pattern) {
  roots <- c(file.path(ROOT, "results"), file.path(ROOT, "data", "processed"),
             file.path(ROOT, "figures"))
  roots <- roots[dir.exists(roots)]
  if (!length(roots)) return(character(0))
  unlist(lapply(roots, function(r)
    list.files(r, pattern = pattern, recursive = TRUE, full.names = TRUE,
               ignore.case = TRUE)), use.names = FALSE)
}

preview <- function(path, n = 3) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    say("      (RDS; not loaded here)")
    return(invisible(NULL))
  }
  x <- try(fread(path, nrows = n, showProgress = FALSE), silent = TRUE)
  if (inherits(x, "try-error")) { say("      could not be read as a table"); return(invisible(NULL)) }
  full <- try(fread(path, select = 1L, showProgress = FALSE), silent = TRUE)
  nrw  <- if (inherits(full, "try-error")) NA_integer_ else nrow(full)
  say("      rows: ", ifelse(is.na(nrw), "?", nrw), " | columns: ", ncol(x))
  say("      names: ", paste(names(x), collapse = ", "))
  for (i in seq_len(nrow(x))) {
    vals <- vapply(x[i], function(v) {
      v <- as.character(v)[1]
      if (is.na(v)) "NA" else if (nchar(v) > 18) paste0(substr(v, 1, 15), "...") else v
    }, character(1))
    say("      row ", i, ": ", paste(vals, collapse = " | "))
  }
  invisible(NULL)
}

report <- function(panel, what, patterns, produced_by, needed_cols) {
  say("")
  say("----------------------------------------------------------------")
  say(" PANEL ", panel, " -- ", what)
  say("----------------------------------------------------------------")
  say("  needs: ", needed_cols)
  hits <- unique(unlist(lapply(patterns, search_all), use.names = FALSE))
  if (!length(hits)) {
    say("  STATUS: NOT FOUND")
    say("  produced by: ", produced_by)
    return(invisible(FALSE))
  }
  say("  STATUS: ", length(hits), " candidate file(s)")
  for (h in hits) {
    say("    ", sub(paste0("^", ROOT, "/?"), "", h),
        sprintf("  (%.1f KB)", file.size(h) / 1024))
    preview(h)
  }
  invisible(TRUE)
}

## ---- Panels A and B : bulk nevus vs melanoma --------------------------------
report("A/B", "GSE112509 and GSE98394, per-sample RASSF9 expression",
       c("gse112509.*\\.csv$", "gse98394.*\\.csv$",
         "nevus.*melanoma.*\\.csv$", "bulk.*cohort.*\\.csv$",
         "bulk_vs_visium.*\\.csv$"),
       "the GEO cohort scripts (14 / 20 era) and script 31",
       "one row per sample: sample id, group (nevus / melanoma), cohort, log2 expression")

report("B", "Cohen's d and 95% CI per cohort",
       c("cohen.*\\.csv$", "effect.*size.*\\.csv$", "bulk_vs_visium.*\\.csv$"),
       "the same scripts; the artwork quotes d = -0.89 [-1.39, -0.38] and -1.00 [-1.49, -0.51]",
       "cohort, d, lo, hi, p, n_nevus, n_melanoma")

## ---- Panel C : single-cell ---------------------------------------------------
report("C", "GSE72056 depth-matched detection by cell type",
       c("tirosh.*\\.csv$", "38[a-z]?_.*\\.csv$", "gse72056.*\\.csv$",
         "singlecell.*\\.csv$"),
       "script 38b",
       "cell type, n cells, n RASSF9-positive, percent  (malignant 1/494 = 0.20%; endothelial 26.83%)")

say("")
say("  NOTE for panel C: the figure should carry the endothelial comparator")
say("  alongside the malignant count. 0.20% on its own invites the question")
say("  of whether detection was sensitive enough; 26.83% in the same")
say("  depth-matched analysis answers it. If no table exists, the numbers can")
say("  be read from results/logs/38b_tirosh_detection.log.")

## ---- Panel D : TCGA per-patient ----------------------------------------------
report("D", "TCGA-SKCM per-patient promoter methylation and RASSF9 expression",
       c("methylation_full_addedsamples\\.csv$", "tcga.*per.*sample.*\\.csv$",
         "skcm.*patient.*\\.csv$", "methylation.*sample.*\\.csv$"),
       "script 35 -- which wrote summary tables but may not have saved the per-patient values",
       "patient, RASSF9 log2(TPM+1), mean promoter beta, sample type, cohort flag")

say("")
say("  NOTE for panel D: if no per-patient table exists, it has to be")
say("  regenerated from skcm_hm450.rds and skcm_star_counts.rds, which is a")
say("  couple of minutes of loading. Say so and I will write that script")
say("  rather than have you hunt for a file that was never written.")

## ---- Panels E and F : cell lines ---------------------------------------------
report("E/F", "cell-line methylation and expression",
       c("46_line_level_data\\.csv$"),
       "script 46 -- this one is known to exist",
       "model_id, name, lineage, disease, subtype, RASSF9, RASSF9_core/shore/distal, RASSF9_cn")

## ---- Supplementary spatial figure --------------------------------------------
say("")
say("----------------------------------------------------------------")
say(" SUPPLEMENTARY -- GSE298774 spatial panels")
say("----------------------------------------------------------------")
say("  Moving to supplementary per the decision of 3 Nov. Existing artwork")
say("  and outputs that may be reusable:")
for (p in c("visium", "gse298774", "bulk_vs_visium")) {
  hits <- search_all(paste0(p, ".*\\.(csv|png|pdf|rds)$"))
  if (!length(hits)) { say("    ", p, ": nothing found"); next }
  say("    ", p, ": ", length(hits), " file(s)")
  for (h in head(hits, 12))
    say("      ", sub(paste0("^", ROOT, "/?"), "", h),
        sprintf("  (%.1f KB)", file.size(h) / 1024))
  if (length(hits) > 12) say("      ... and ", length(hits) - 12, " more")
}

## ---- Environment -------------------------------------------------------------
say("")
say("----------------------------------------------------------------")
say(" Environment")
say("----------------------------------------------------------------")
say("  R: ", R.version.string)
for (p in c("ggplot2", "patchwork", "data.table", "scales", "ggrepel")) {
  v <- try(as.character(utils::packageVersion(p)), silent = TRUE)
  say(sprintf("  %-12s %s", p, ifelse(inherits(v, "try-error"), "NOT INSTALLED", v)))
}

say("")
say("================================================================")
say(" Send this log back. With it I can write the figure script once,")
say(" against the real column names, instead of guessing and iterating.")
say("================================================================")
