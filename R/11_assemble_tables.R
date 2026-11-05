# ==============================================================================
# 11_assemble_tables.R
#
# Purpose: turn the raw analysis CSVs into a manuscript-numbered set of tables,
# ready to attach to an email or drop into a submission.
#
# Two numbering schemes are supported. Set SCHEME at the top:
#
#   "collapsed"  (default) Tables 1-3 of the one-to-one scheme are merged into a
#                single "RASSF9 and UV mutational burden" table, because they are
#                three views of one analysis. Gives 3 main tables:
#                  Table 1  RASSF9 and UV mutational burden (groups + correlations + models)
#                  Table 2  RASSF9 promoter methylation
#                  Table 3  RASSF9 and markers of aggressive disease
#
#   "one_to_one" One file, one table, no merging. Gives 5 main tables:
#                  Table 1  Group comparison
#                  Table 2  Correlations
#                  Table 3  Adjusted models
#                  Table 4  Promoter methylation
#                  Table 5  Aggressiveness
#
# Supplementary tables (RASSF10, then the exploratory mediation) are numbered to
# follow whichever scheme is chosen.
#
# The original CSVs are COPIED, never renamed or moved: the scripts that produce
# them must stay reproducible, and a renamed input would break a re-run.
#
# Input:  results/*.csv  (produced by scripts 04-10)
# Output: results/manuscript_tables/   numbered CSVs + TABLE_LEGENDS.txt
#                                      + MANIFEST.csv mapping source to table
#         results/logs/11_assemble.log
# ==============================================================================

source(here::here("R", "00_setup.R"))

SCHEME <- "collapsed"      # "collapsed" or "one_to_one"

LOG <- file.path(PATHS$logs, "11_assemble.log")
cat("Table assembly —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

# Top-level folder so it is easy to find and open: E:\RASSF9 Project\FINAL_TABLES
OUT <- here::here("FINAL_TABLES")
if (dir.exists(OUT)) unlink(OUT, recursive = TRUE)   # rebuild cleanly each run
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

log_line("Scheme: ", SCHEME, logfile = LOG)
log_line("Output: ", OUT, logfile = LOG)

# ==============================================================================
# 1. Load the analysis outputs
# ==============================================================================
hdr("1. Source files")

SRC <- c(
  groups_9    = "table_main_04_rassf9_groups.csv",
  corr_9      = "table_main_05_rassf9_correlation.csv",
  adj_9       = "table_main_06_rassf9_adjusted.csv",
  meth_9      = "table_main_08_rassf9_methylation.csv",
  aggr_9      = "table_main_09_rassf9_aggressiveness.csv",
  groups_10   = "table_supp_04_rassf10_groups.csv",
  corr_10     = "table_supp_05_rassf10_correlation.csv",
  adj_10      = "table_supp_06_rassf10_adjusted.csv",
  meth_10     = "table_supp_08_rassf10_methylation.csv",
  aggr_10     = "table_supp_09_rassf10_aggressiveness.csv",
  mediation   = "table_exploratory_10_mediation.csv"
)

dat <- list()
for (k in names(SRC)) {
  f <- file.path(PATHS$results, SRC[[k]])
  if (!file.exists(f)) {
    log_line("  [MISSING] ", SRC[[k]], " -- run the script that produces it",
             logfile = LOG)
    next
  }
  dat[[k]] <- readr::read_csv(f)
  log_line(sprintf("  [ ok ] %-44s %d rows", SRC[[k]], nrow(dat[[k]])),
           logfile = LOG)
}

if (length(dat) < length(SRC)) {
  log_line("\n  Some sources are missing; the tables below cover what exists.",
           logfile = LOG)
}

# ==============================================================================
# 2. Build the merged UV table
# ==============================================================================
hdr("2. Merged UV table")

# The three analyses report different quantities, so they are harmonised into
# one long format: what was tested, on what, with what estimate and p value.
# Nothing is recomputed here -- every number is carried through unchanged.

num <- function(x) suppressWarnings(as.numeric(x))

merge_uv <- function(groups, corr, adj, gene) {
  out <- list()

  if (!is.null(groups)) {
    out[[1]] <- tibble::tibble(
      section   = "Group comparison (UV-high vs UV-low)",
      analysis  = groups$analysis,
      measure   = "Cohen's d",
      n         = groups$n_low + groups$n_high,
      detail    = sprintf("median %.3f vs %.3f (n = %d vs %d)",
                          groups$median_low, groups$median_high,
                          groups$n_low, groups$n_high),
      estimate  = num(groups$cohens_d),
      ci_low    = num(groups$d_lower),
      ci_high   = num(groups$d_upper),
      statistic = sprintf("W = %.0f", num(groups$W)),
      p_value   = num(groups$p_value)
    )
  }

  if (!is.null(corr)) {
    out[[2]] <- tibble::tibble(
      section   = "Correlation with burden measures",
      analysis  = corr$predictor,
      measure   = corr$method,
      n         = corr$n,
      detail    = corr$variable,
      estimate  = num(corr$estimate),
      ci_low    = num(corr$conf_low),
      ci_high   = num(corr$conf_high),
      statistic = NA_character_,
      p_value   = num(corr$p_value)
    )
  }

  if (!is.null(adj)) {
    out[[3]] <- tibble::tibble(
      section   = "Linear models (UV coefficient)",
      analysis  = adj$model,
      measure   = "beta per log10 SBS7",
      n         = adj$n,
      detail    = ifelse(is.na(adj$purity_used), "no covariates",
                         paste0("purity = ", adj$purity_used,
                                ", + sample type")),
      estimate  = num(adj$uv_beta),
      ci_low    = num(adj$uv_ci_low),
      ci_high   = num(adj$uv_ci_high),
      statistic = ifelse(is.na(adj$adj_r2), NA_character_,
                         sprintf("adj R2 = %.3f", num(adj$adj_r2))),
      p_value   = num(adj$uv_p)
    )
  }

  res <- dplyr::bind_rows(out)
  if (nrow(res)) res <- dplyr::mutate(res, gene = gene, .before = 1)
  res
}

uv_9  <- merge_uv(dat$groups_9,  dat$corr_9,  dat$adj_9,  "RASSF9")
uv_10 <- merge_uv(dat$groups_10, dat$corr_10, dat$adj_10, "RASSF10")

log_line("  RASSF9 merged table:  ", nrow(uv_9), " rows", logfile = LOG)
log_line("  RASSF10 merged table: ", nrow(uv_10), " rows", logfile = LOG)

# ==============================================================================
# 3. Assign table numbers
# ==============================================================================
hdr("3. Table numbering")

# Each entry: number, title, and either a data frame or a source key to copy.
if (SCHEME == "collapsed") {
  # NOTE ON NUMBERING. The manuscript already contains Table 1 (the qRT-PCR
  # primer list), so the TCGA tables start at Table 2. Supplementary tables
  # start at S1: the existing Supplementary 1 and 2 are FIGURES, not tables,
  # so there is no collision there.
  plan <- list(
    list(num = "Table_2",  title = "RASSF9 expression and UV mutational burden",
         data = uv_9,  src = "table_main_04/05/06 merged"),
    list(num = "Table_3",  title = "RASSF9 promoter methylation and expression",
         key = "meth_9"),
    list(num = "Table_S1", title = "RASSF9 expression and markers of aggressive disease",
         key = "aggr_9"),
    list(num = "Table_S2", title = "RASSF10 expression and UV mutational burden",
         data = uv_10, src = "table_supp_04/05/06 merged"),
    list(num = "Table_S3", title = "RASSF10 promoter methylation and expression",
         key = "meth_10"),
    list(num = "Table_S4", title = "RASSF10 expression and markers of aggressive disease",
         key = "aggr_10"),
    list(num = "Table_S5", title = "Exploratory mediation analysis",
         key = "mediation")
  )
} else if (SCHEME == "one_to_one") {
  plan <- list(
    list(num = "Table_1",  title = "RASSF9 expression in UV-high versus UV-low melanomas",
         key = "groups_9"),
    list(num = "Table_2",  title = "Correlation of RASSF9 expression with burden measures",
         key = "corr_9"),
    list(num = "Table_3",  title = "Adjusted models of RASSF9 expression on UV burden",
         key = "adj_9"),
    list(num = "Table_4",  title = "RASSF9 promoter methylation and expression",
         key = "meth_9"),
    list(num = "Table_5",  title = "RASSF9 expression and markers of aggressive disease",
         key = "aggr_9"),
    list(num = "Table_S1", title = "RASSF10 expression in UV-high versus UV-low melanomas",
         key = "groups_10"),
    list(num = "Table_S2", title = "Correlation of RASSF10 expression with burden measures",
         key = "corr_10"),
    list(num = "Table_S3", title = "Adjusted models of RASSF10 expression on UV burden",
         key = "adj_10"),
    list(num = "Table_S4", title = "RASSF10 promoter methylation and expression",
         key = "meth_10"),
    list(num = "Table_S5", title = "RASSF10 expression and markers of aggressive disease",
         key = "aggr_10"),
    list(num = "Table_S6", title = "Exploratory mediation analysis",
         key = "mediation")
  )
} else {
  stop("SCHEME must be \"collapsed\" or \"one_to_one\"; got: ", SCHEME)
}

# ==============================================================================
# 4. Write the numbered tables
# ==============================================================================
hdr("4. Writing tables")

slug <- function(x) {
  x <- gsub("[^A-Za-z0-9 ]", "", x)
  gsub(" +", "_", trimws(x))
}

manifest <- list()
for (p in plan) {
  d   <- if (!is.null(p$data)) p$data else dat[[p$key]]
  src <- if (!is.null(p$src))  p$src  else SRC[[p$key]]

  if (is.null(d) || !nrow(d)) {
    log_line("  [skip] ", p$num, " -- source not available", logfile = LOG)
    next
  }

  fname <- paste0(p$num, "_", slug(p$title), ".csv")
  readr::write_csv(d, file.path(OUT, fname))

  manifest[[length(manifest) + 1]] <- tibble::tibble(
    table = gsub("_", " ", p$num),
    title = p$title,
    file  = fname,
    rows  = nrow(d),
    source = src
  )
  log_line(sprintf("  %-9s %-58s %d rows", gsub("_", " ", p$num), fname,
                   nrow(d)), logfile = LOG)
}

man <- dplyr::bind_rows(manifest)
readr::write_csv(man, file.path(OUT, "MANIFEST.csv"))

# ==============================================================================
# 5. Legends file
# ==============================================================================
hdr("5. Legends")

legend_for <- function(tab, title) {
  base <- paste0(tab, ". ", title, ".")
  body <- dplyr::case_when(
    grepl("UV mutational burden", title) & SCHEME == "collapsed" ~
      paste("Group comparison, correlation and linear-model results for the",
            "association between expression and UV mutational burden.",
            "UV-low and UV-high groups are the outer tertiles of total SBS7",
            "burden (at or below 324 and at or above 714 mutations;",
            "n = 107 and 108). Correlations are given for total SBS7 burden,",
            "SBS7 as a fraction of total signature burden, and total",
            "mutational burden. Linear models are unadjusted and adjusted for",
            "tumour purity and sample type, with four purity measures and a",
            "rank-inverse-normal sensitivity analysis. n = 320."),
    grepl("UV-high versus UV-low", title) ~
      paste("Wilcoxon rank-sum comparison between the outer tertiles of total",
            "SBS7 burden (UV-low at or below 324 mutations, n = 107; UV-high",
            "at or above 714, n = 108), with Cohen's d and 95% confidence",
            "intervals, for the primary cohort and two pre-specified",
            "sensitivity analyses."),
    grepl("Correlation of", title) ~
      paste("Spearman and Pearson correlations with total SBS7 burden, SBS7",
            "fraction and total mutational burden. Pearson coefficients are",
            "given on raw and log10-transformed burden scales. n = 320."),
    grepl("Adjusted models", title) ~
      paste("Linear regression coefficients for log10(SBS7 + 1) predicting",
            "log2(TPM + 1) expression, unadjusted and adjusted for tumour",
            "purity and sample type, across four purity measures and a",
            "rank-inverse-normal transformed outcome. n = 320 (310 for LUMP)."),
    grepl("promoter methylation", title) ~
      paste("Spearman correlations between probe beta value and expression",
            "for each promoter probe, with position relative to the",
            "transcription start site, median beta and Benjamini-Hochberg",
            "adjusted p values; and for the promoter mean, unadjusted and",
            "covariate-adjusted linear models. n = 320."),
    grepl("aggressive disease", title) ~
      paste("Tests against tumour site, AJCC stage, Breslow thickness and",
            "overall survival, with Benjamini-Hochberg correction across the",
            "four measures. Sample sizes differ by measure because clinical",
            "annotation is incomplete in TCGA-SKCM (stage n = 86, Breslow",
            "n = 20, survival n = 309 with 125 deaths)."),
    grepl("mediation", title) ~
      paste("Path estimates for an indirect effect of UV mutational burden on",
            "expression through promoter methylation, with percentile",
            "bootstrap confidence intervals from 5,000 resamples. All paths",
            "adjusted for tumour purity and sample type; variables",
            "standardised. Not pre-specified; reported as exploratory."),
    TRUE ~ ""
  )
  paste0(base, " ", body)
}

leg_file <- file.path(OUT, "TABLE_LEGENDS.txt")
cat("TABLE LEGENDS\n",
    "Generated ", format(Sys.Date()), " | scheme: ", SCHEME, "\n\n",
    sep = "", file = leg_file)

for (i in seq_len(nrow(man))) {
  cat(legend_for(man$table[i], man$title[i]), "\n\n",
      sep = "", file = leg_file, append = TRUE)
}
log_line("  Wrote: TABLE_LEGENDS.txt (", nrow(man), " legends)", logfile = LOG)

# --- README, so the folder explains itself to anyone who opens it ------------
readme <- file.path(OUT, "README.txt")
cat(
  "TCGA-SKCM TABLES FOR THE RASSF9 MANUSCRIPT\n",
  "Generated ", format(Sys.Date()), "\n\n",
  "NUMBERING\n",
  "These tables continue the numbering of the existing manuscript, which\n",
  "already contains Table 1 (the qRT-PCR primer list). The TCGA tables\n",
  "therefore begin at Table 2. Supplementary tables begin at S1, because the\n",
  "manuscript's existing Supplementary 1 and 2 are FIGURES, not tables.\n\n",
  "CONTENTS\n",
  paste0("  ", man$table, "  ", man$title, "\n", collapse = ""),
  "\nTABLE_LEGENDS.txt holds a draft legend for each table.\n",
  "MANIFEST.csv records which analysis output each table was built from.\n\n",
  "SOURCE\n",
  "Built by R/11_assemble_tables.R from the CSVs in results/. Those originals\n",
  "are untouched and remain the reproducible outputs of scripts 04 to 10.\n",
  sep = "", file = readme)
log_line("  Wrote: README.txt", logfile = LOG)

# ==============================================================================
# 6. Summary
# ==============================================================================
hdr("6. Summary")

cat(capture.output(print(as.data.frame(man[, c("table", "title", "file")]),
                         row.names = FALSE)),
    sep = "\n", file = LOG, append = TRUE)
print(as.data.frame(man[, c("table", "file", "rows")]), row.names = FALSE)

log_line("\n  ", nrow(man), " tables written to results/manuscript_tables/",
         logfile = LOG)
log_line("  Original analysis CSVs are untouched in results/.", logfile = LOG)
log_line("  To switch numbering, set SCHEME at the top of this script and",
         logfile = LOG)
log_line("  re-run; the output folder is rebuilt from scratch each time.",
         logfile = LOG)

report_warnings(LOG)
write_session_info("11_assemble")

log_line("\nFull log: ", LOG, logfile = LOG)
