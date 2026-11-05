## =============================================================================
## 42_ccle_expression_gate.R
##
## STEP 1 of Gustavo's cell-line question (email 23 Oct 2026).
##
## PURPOSE
##   Decide whether a methylation-expression correlation in melanoma cell lines
##   is testable AT ALL, before any methylation data is downloaded.
##
## WHY THIS IS THE FIRST STEP
##   Our own single-cell work (Tirosh, script 38b) found RASSF9 in 0.20% of
##   malignant melanoma cells against 26.83% of endothelial cells. A melanoma
##   cell line is pure malignant cells with no stroma. If the bulk-tumour
##   RASSF9 signal in TCGA is substantially non-malignant in origin, then
##   melanoma lines should express little or no RASSF9, and a Spearman
##   correlation against methylation would be computed on noise and ties.
##
##   Gustavo asked for exactly this check (his point 4) and asked that we
##   report inadequacy rather than force a conclusion. So the gate is
##   pre-specified below and the script states its own verdict.
##
## INPUTS  (download into data/raw/ccle/ from depmap.org, release 25Q2)
##   OmicsExpressionProteinCodingGenesTPMLogp1.csv   ~500 MB, log2(TPM+1)
##   Model.csv                                        model annotation
##
## NOTE ON MEMORY
##   The expression matrix is never loaded whole. The header is read alone,
##   the handful of wanted columns are located by pattern, and only those are
##   read back. Peak memory is a few MB.
##
## OUTPUTS
##   results/ccle/42_rassf9_melanoma_lines.csv    per-line values, RASSF9
##   results/ccle/42_rassf10_melanoma_lines.csv   per-line values, RASSF10
##   results/ccle/42_gate_summary.csv             the numbers behind the verdict
##   figures/ccle/42_rassf9_detection.png         RASSF9 only
##   figures/ccle/42_rassf10_detection.png        RASSF10 only
##   results/logs/42_ccle_expression_gate.log
## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------
if (file.exists("R/00_setup.R")) source("R/00_setup.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT     <- if (exists("PROJECT_ROOT")) PROJECT_ROOT else getwd()
RAW_DIR  <- file.path(ROOT, "data", "raw", "ccle")
RES_DIR  <- file.path(ROOT, "results", "ccle")
FIG_DIR  <- file.path(ROOT, "figures", "ccle")
LOG_DIR  <- file.path(ROOT, "results", "logs")
for (d in c(RES_DIR, FIG_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_PATH <- file.path(LOG_DIR, "42_ccle_expression_gate.log")
.logcon  <- file(LOG_PATH, open = "wt")
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = .logcon)
}
WARN <- character(0)
warn <- function(...) {
  msg <- paste0(...)
  WARN <<- c(WARN, msg)
  say("  [WARN] ", msg)
}
on.exit({
  if (length(WARN)) {
    say("")
    say("WARNINGS RAISED (", length(WARN), "):")
    for (w in WARN) say("  - ", w)
  } else say("  no warnings raised")
  say("")
  say("Full log: ", LOG_PATH)
  close(.logcon)
}, add = TRUE)

set.seed(1234)

GENES <- c("RASSF9", "RASSF10")

## Technical reference genes. These are NOT a biological comparison -- they
## exist only to prove the file parsed correctly and that a gene genuinely
## expressed in these lines looks expressed. Gustavo asked us not to add
## normal-melanocyte or skin comparisons, and these are not that.
REF_HIGH <- c("ACTB", "GAPDH")   # should be high everywhere
REF_MEL  <- c("MITF", "SOX10")   # melanocytic lineage, expected on in most lines
REF_OFF  <- c("MYOD1")           # skeletal muscle TF, expected near zero

## ---- Pre-specified gate -----------------------------------------------------
## Thresholds fixed before seeing any numbers.
##
##   PROCEED  : >= 70% of melanoma lines at TPM >= 1, AND SD of log2(TPM+1)
##              >= 0.50, AND < 40% of values tied at the file's floor
##   MARGINAL : 40-70% detected, or SD 0.30-0.50
##   STOP     : < 40% detected, or SD < 0.30
##
## Rationale for the SD cut: Spearman on a variable with SD < 0.3 log2 units
## means most lines sit within ~20% of one another, so the ranks being
## correlated are measurement noise rather than biology. The tie check matters
## for the same reason -- a Spearman dominated by a block of identical zeros
## reports the tie structure, not an association.
GATE_DET_HI  <- 0.70
GATE_DET_LO  <- 0.40
GATE_SD_HI   <- 0.50
GATE_SD_LO   <- 0.30
GATE_TIE_MAX <- 0.40
TPM_CUT      <- 1        # TPM >= 1  <=>  log2(TPM+1) >= 1

say("")
say("== 1. Inputs ===================================================")

## DepMap renames files between releases. Accept any of the known names for
## the protein-coding log2(TPM+1) matrix and report which one was found, so the
## release actually used is recorded rather than assumed.
EXPR_CANDIDATES <- c(
  "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",          # 26Q3 name
  "OmicsExpressionProteinCodingGenesTPMLogp1.csv",               # 25Q2 and earlier
  "OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv"   # stranded variant
)

found  <- EXPR_CANDIDATES[file.exists(file.path(RAW_DIR, EXPR_CANDIDATES))]
f_expr <- if (length(found)) file.path(RAW_DIR, found[1]) else NA_character_

f_model <- file.path(RAW_DIR, "Model.csv")

if (is.na(f_expr) || !file.exists(f_model)) {
  say("  MISSING INPUT")
  say("")
  say("  Expected in: ", RAW_DIR)
  say("    Model.csv                                            ",
      ifelse(file.exists(f_model), "FOUND", "missing"))
  for (cand in EXPR_CANDIDATES)
    say(sprintf("    %-52s %s", cand,
                ifelse(file.exists(file.path(RAW_DIR, cand)), "FOUND", "missing")))
  say("")
  say("  Download from https://depmap.org/portal/ -> Data -> All Data Downloads")
  stop("required input missing", call. = FALSE)
}

if (length(found) > 1)
  say("  NOTE: ", length(found), " expression files present; using ", found[1])
if (grepl("Stranded", basename(f_expr)))
  warn("using the STRANDED expression matrix -- note this in any write-up, ",
       "as it is not the DepMap default")

for (f in c(f_expr, f_model))
  say(sprintf("  %-52s %8.1f MB", basename(f), file.size(f) / 1024^2))

## ---- 2. Model annotation ----------------------------------------------------
say("")
say("== 2. Melanoma lines in the annotation =========================")

MD <- fread(f_model, showProgress = FALSE)

## Column names have changed across DepMap releases, so locate them by pattern
## and print what was found rather than hard-coding and silently mismatching.
pick <- function(dt, patterns, what) {
  for (p in patterns) {
    hit <- grep(p, names(dt), ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  warn("could not locate a column for ", what,
       " -- available: ", paste(head(names(dt), 40), collapse = ", "))
  NA_character_
}

col_id      <- pick(MD, c("^ModelID$", "^DepMap_ID$", "ModelID"),            "model ID")
col_disease <- pick(MD, c("OncotreePrimaryDisease", "primary_disease"),      "primary disease")
col_lineage <- pick(MD, c("OncotreeLineage", "^lineage$"),                   "lineage")
col_name    <- pick(MD, c("StrippedCellLineName", "CellLineName", "^cell_line_name$"), "cell line name")
col_subtype <- pick(MD, c("OncotreeSubtype", "^Subtype$"),                   "subtype")

say("  model ID column   : ", col_id)
say("  disease column    : ", col_disease)
say("  lineage column    : ", col_lineage)
say("  name column       : ", col_name)
say("  subtype column    : ", col_subtype)

if (is.na(col_id) || is.na(col_disease)) stop("cannot identify models without ID and disease columns", call. = FALSE)

MD[[col_disease]] <- as.character(MD[[col_disease]])
is_mel <- grepl("melanoma", MD[[col_disease]], ignore.case = TRUE)

say("")
say("  total models in annotation: ", nrow(MD))
say("  models with 'melanoma' in primary disease: ", sum(is_mel))
say("")
say("  breakdown of those by primary disease:")
tb <- sort(table(MD[[col_disease]][is_mel]), decreasing = TRUE)
for (i in seq_along(tb)) say(sprintf("    %-52s %4d", names(tb)[i], tb[i]))

if (!is.na(col_subtype)) {
  say("")
  say("  by subtype (acral / mucosal / uveal are biologically distinct):")
  tb2 <- sort(table(as.character(MD[[col_subtype]][is_mel])), decreasing = TRUE)
  for (i in seq_along(tb2)) say(sprintf("    %-52s %4d", names(tb2)[i], tb2[i]))
}

## Uveal melanoma is a different disease with a different mutational and
## epigenetic background. Keep it identifiable but exclude it from the primary
## set, and report both so the choice is visible rather than buried.
uveal <- rep(FALSE, nrow(MD))
if (!is.na(col_subtype)) uveal <- grepl("uveal", as.character(MD[[col_subtype]]), ignore.case = TRUE)
uveal <- uveal | grepl("uveal", MD[[col_disease]], ignore.case = TRUE)

MEL <- data.table(
  model_id = MD[[col_id]][is_mel],
  name     = if (!is.na(col_name))    as.character(MD[[col_name]][is_mel])    else NA_character_,
  disease  = MD[[col_disease]][is_mel],
  subtype  = if (!is.na(col_subtype)) as.character(MD[[col_subtype]][is_mel]) else NA_character_,
  lineage  = if (!is.na(col_lineage)) as.character(MD[[col_lineage]][is_mel]) else NA_character_,
  uveal    = uveal[is_mel]
)

say("")
say("  cutaneous/other melanoma (primary set): ", sum(!MEL$uveal))
say("  uveal melanoma (reported separately)  : ", sum(MEL$uveal))

## ---- 3. Locate the wanted gene columns --------------------------------------
say("")
say("== 3. Locating gene columns in the expression matrix ============")

## Column names in this file are "SYMBOL (ENTREZID)", e.g. "RASSF9 (9182)".
hdr <- names(fread(f_expr, nrows = 0, showProgress = FALSE))
say("  columns in file: ", length(hdr))
say("  first four columns: ", paste(sprintf("'%s'", head(hdr, 4)), collapse = ", "))

## Layout differs by release:
##   26Q3  : ProfileID, is_default_entry, ModelID, then one column per gene.
##           A model can appear more than once (one row per profile).
##   <=25Q2: a single unnamed model-ID column, then one column per gene.
## Detect rather than assume -- getting this wrong silently merges on the
## wrong key and reports zero melanoma lines.
id_col   <- if ("ModelID" %in% hdr) "ModelID" else hdr[1]
flag_col <- if ("is_default_entry" %in% hdr) "is_default_entry" else NA_character_

say("  model ID column   : '", id_col, "'")
say("  default-row flag  : ",
    ifelse(is.na(flag_col), "absent (one row per model assumed)", paste0("'", flag_col, "'")))
if (id_col == hdr[1] && !("ModelID" %in% hdr))
  say("  (older layout: using the first column as the model ID)")

want    <- c(GENES, REF_HIGH, REF_MEL, REF_OFF)
col_map <- setNames(rep(NA_character_, length(want)), want)

for (g in want) {
  ## anchor on the symbol followed by space-paren, so RASSF9 does not match
  ## RASSF9P1 or similar
  hit <- grep(paste0("^", g, " \\("), hdr, value = TRUE)
  if (length(hit) == 1) {
    col_map[g] <- hit
  } else if (length(hit) > 1) {
    col_map[g] <- hit[1]
    warn(g, " matched ", length(hit), " columns (", paste(hit, collapse = ", "),
         ") -- using the first")
  } else {
    ## fall back to a bare symbol column name in case the format changed
    hit2 <- grep(paste0("^", g, "$"), hdr, value = TRUE)
    if (length(hit2) == 1) col_map[g] <- hit2
    else warn(g, " not found in the expression matrix")
  }
}

for (g in want) say(sprintf("  %-9s -> %s", g, ifelse(is.na(col_map[g]), "NOT FOUND", col_map[g])))

if (is.na(col_map["RASSF9"])) stop("RASSF9 column not found -- cannot run the gate", call. = FALSE)

## ---- 4. Read only those columns ---------------------------------------------
say("")
say("== 4. Reading the selected columns =============================")

sel <- unique(c(id_col, if (!is.na(flag_col)) flag_col, unname(col_map[!is.na(col_map)])))
EX  <- fread(f_expr, select = sel, showProgress = FALSE)
setnames(EX, id_col, "model_id")
for (g in want) if (!is.na(col_map[g])) setnames(EX, col_map[g], g)

say("  rows read: ", nrow(EX))

## One row per sequencing profile in 26Q3. Keep the row DepMap designates as
## representing the model, then drop any remaining duplicate so each line
## contributes exactly one value. Both steps are counted in the log.
if (!is.na(flag_col)) {
  fv <- tolower(trimws(as.character(EX[[flag_col]])))
  keep <- fv %in% c("yes", "true", "t", "1")
  say("  rows flagged as the default entry for their model: ", sum(keep))
  if (sum(keep) == 0) {
    warn("no rows matched the default-entry flag (values seen: ",
         paste(unique(head(fv, 20)), collapse = ", "),
         ") -- keeping all rows and de-duplicating by model instead")
  } else {
    EX <- EX[keep]
  }
  EX[[flag_col]] <- NULL
}

n_before <- nrow(EX)
EX <- EX[!duplicated(EX$model_id)]
if (n_before != nrow(EX))
  say("  dropped ", n_before - nrow(EX), " further duplicate model rows (kept the first)")

say("  distinct models retained: ", nrow(EX))
say("  units: log2(TPM + 1), per the DepMap file description")

## Sanity check the units before trusting anything downstream. log2(TPM+1)
## should be non-negative and ACTB should be high.
if (!is.na(col_map["ACTB"])) {
  a <- EX$ACTB[is.finite(EX$ACTB)]
  say(sprintf("  ACTB across all models: median %.2f  range %.2f to %.2f",
              median(a), min(a), max(a)))
  if (median(a) < 5) warn("ACTB median is ", round(median(a), 2),
                          " on the log2(TPM+1) scale -- unexpectedly low; check the file/units")
}
if (any(EX$RASSF9 < 0, na.rm = TRUE)) warn("negative RASSF9 values present -- the file may not be log2(TPM+1)")

## ---- 5. Restrict to melanoma ------------------------------------------------
say("")
say("== 5. Matching annotation to expression ========================")

M <- merge(MEL, EX, by = "model_id", all.x = TRUE)
say("  melanoma models in annotation      : ", nrow(MEL))
say("  of those, present in expression    : ", sum(!is.na(M$RASSF9)))
say("  of those, missing from expression  : ", sum(is.na(M$RASSF9)))

M <- M[!is.na(M$RASSF9)]
PRI <- M[!M$uveal]      # primary set: non-uveal melanoma
UVE <- M[ M$uveal]

say("")
say("  PRIMARY SET (non-uveal melanoma) n = ", nrow(PRI))
say("  uveal set                        n = ", nrow(UVE))
if (nrow(PRI) < 20) warn("only ", nrow(PRI), " non-uveal melanoma lines with expression -- ",
                         "a Spearman correlation on this many lines is poorly powered")

## ---- 6. Detection and variance ----------------------------------------------
say("")
say("== 6. Detection and dynamic range ==============================")

## The file's floor is the smallest value present; with log2(TPM+1) this is
## normally exactly 0. Values tied at the floor are the ones that destroy a
## rank correlation, so they are counted explicitly.
describe <- function(v, label, n_all) {
  v  <- v[is.finite(v)]
  fl <- min(c(v, 0))
  data.table(
    gene        = label,
    n           = length(v),
    n_detected  = sum(v >= log2(TPM_CUT + 1)),
    pct_det     = 100 * mean(v >= log2(TPM_CUT + 1)),
    pct_zero    = 100 * mean(v <= fl + 1e-9),
    median      = median(v),
    iqr         = IQR(v),
    sd          = sd(v),
    min         = min(v),
    max         = max(v),
    n_distinct  = length(unique(round(v, 6))),
    pct_tied    = 100 * (1 - length(unique(round(v, 6))) / max(length(v), 1))
  )
}

SUM <- rbindlist(lapply(want[!is.na(col_map[want])],
                        function(g) describe(PRI[[g]], g, nrow(PRI))))

say("  non-uveal melanoma lines, n = ", nrow(PRI))
say("")
say(sprintf("  %-9s %5s %7s %7s %7s %7s %7s %7s %7s %8s",
            "gene", "n", "%TPM>1", "%floor", "median", "IQR", "SD", "min", "max", "%tied"))
for (i in seq_len(nrow(SUM))) {
  s <- SUM[i]
  say(sprintf("  %-9s %5d %6.1f%% %6.1f%% %7.2f %7.2f %7.2f %7.2f %7.2f %7.1f%%",
              s$gene, s$n, s$pct_det, s$pct_zero, s$median, s$iqr, s$sd, s$min, s$max, s$pct_tied))
}

say("")
say("  Read ACTB/GAPDH as a file-integrity check, MITF/SOX10 as evidence the")
say("  lines are melanocytic, MYOD1 as a near-zero reference. They are not a")
say("  biological comparison for RASSF9.")

if (nrow(UVE) >= 5) {
  say("")
  say("  uveal melanoma, n = ", nrow(UVE), " (reported, not pooled):")
  for (g in GENES) if (!is.na(col_map[g])) {
    v <- UVE[[g]][is.finite(UVE[[g]])]
    say(sprintf("    %-9s %%TPM>1 %5.1f%%  median %5.2f  SD %5.2f",
                g, 100 * mean(v >= log2(TPM_CUT + 1)), median(v), sd(v)))
  }
}

## ---- 7. Verdict -------------------------------------------------------------
say("")
say("== 7. Verdict against the pre-set gate =========================")
say("  PROCEED  : >= ", round(100 * GATE_DET_HI), "% of lines at TPM >= ", TPM_CUT,
    ", SD >= ", GATE_SD_HI, ", tied < ", round(100 * GATE_TIE_MAX), "%")
say("  MARGINAL : ", round(100 * GATE_DET_LO), "-", round(100 * GATE_DET_HI),
    "% detected, or SD ", GATE_SD_LO, "-", GATE_SD_HI)
say("  STOP     : < ", round(100 * GATE_DET_LO), "% detected, or SD < ", GATE_SD_LO)
say("")

verdicts <- character(0)
for (g in GENES) {
  if (is.na(col_map[g])) next
  s   <- SUM[SUM$gene == g]
  det <- s$pct_det / 100
  v <- if (det >= GATE_DET_HI && s$sd >= GATE_SD_HI && s$pct_tied / 100 < GATE_TIE_MAX) {
    "PROCEED"
  } else if (det < GATE_DET_LO || s$sd < GATE_SD_LO) {
    "STOP"
  } else "MARGINAL"
  verdicts[g] <- v

  say("  ", g, ": ", v)
  say(sprintf("    detected %.1f%% of %d lines | SD %.2f | tied %.1f%% | median %.2f (log2 TPM+1)",
              s$pct_det, s$n, s$sd, s$pct_tied, s$median))
  if (v == "STOP") {
    say("    -> a methylation-expression correlation in these lines would be")
    say("       computed largely on the detection floor. Report as inadequate.")
  } else if (v == "MARGINAL") {
    say("    -> testable but underpowered; any result needs the detection")
    say("       fraction reported alongside it, and a sensitivity analysis")
    say("       restricted to detected lines.")
  } else {
    say("    -> proceed to the methylation coverage check (step 2).")
  }
  say("")
}

say("  WHAT THIS CANNOT SAY: this is a feasibility check on expression only.")
say("  A PROCEED verdict does not mean the methylation data covers the RASSF9")
say("  promoter, which is a separate question answered in step 2. A STOP")
say("  verdict is informative in its own right: it is consistent with the")
say("  bulk-tumour RASSF9 signal in TCGA being substantially non-malignant in")
say("  origin, which is what our Tirosh analysis already indicated.")

## ---- 8. Write results -------------------------------------------------------
say("")
say("== 8. Output ===================================================")

keep <- c("model_id", "name", "disease", "subtype", "uveal")
fwrite(PRI[, c(keep, "RASSF9"), with = FALSE],
       file.path(RES_DIR, "42_rassf9_melanoma_lines.csv"))
if (!is.na(col_map["RASSF10"]))
  fwrite(PRI[, c(keep, "RASSF10"), with = FALSE],
         file.path(RES_DIR, "42_rassf10_melanoma_lines.csv"))

SUM$verdict <- verdicts[SUM$gene]
SUM$n_lines_primary <- nrow(PRI)
SUM$n_lines_uveal   <- nrow(UVE)
fwrite(SUM, file.path(RES_DIR, "42_gate_summary.csv"))

say("  wrote 42_rassf9_melanoma_lines.csv")
if (!is.na(col_map["RASSF10"])) say("  wrote 42_rassf10_melanoma_lines.csv")
say("  wrote 42_gate_summary.csv")

## ---- 9. Figures -------------------------------------------------------------
## RASSF9 and RASSF10 in separate files, per the project convention.
say("")
say("== 9. Figures ==================================================")

plot_gene <- function(g, outfile) {
  if (is.na(col_map[g])) return(invisible(NULL))
  d <- data.table(v = PRI[[g]][is.finite(PRI[[g]])])
  s <- SUM[SUM$gene == g]
  cut_y <- log2(TPM_CUT + 1)

  p <- ggplot(d, aes(x = v)) +
    geom_histogram(bins = 30, fill = "grey35", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = cut_y, linetype = "dashed", linewidth = 0.4, colour = "firebrick") +
    annotate("text", x = cut_y, y = Inf, label = paste0("  TPM = ", TPM_CUT),
             hjust = 0, vjust = 1.6, size = 3.1, colour = "firebrick") +
    labs(
      title    = paste0(g, " expression across melanoma cell lines"),
      ## name the source file rather than a release label, so the figure
      ## cannot carry a stale version string
      subtitle = sprintf("%s | n = %d non-uveal melanoma lines | %.1f%% at TPM ≥ %d | SD %.2f | verdict %s",
                         basename(f_expr), s$n, s$pct_det, TPM_CUT, s$sd, verdicts[g]),
      x        = expression(log[2]*"(TPM + 1)"),
      y        = "cell lines"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", size = 11),
      plot.subtitle    = element_text(size = 8.5, colour = "grey25"),
      plot.margin      = margin(10, 14, 8, 10)
    )

  ggsave(outfile, p, width = 6.4, height = 4.0, dpi = 400)
  say("  wrote: ", basename(outfile))
}

plot_gene("RASSF9",  file.path(FIG_DIR, "42_rassf9_detection.png"))
plot_gene("RASSF10", file.path(FIG_DIR, "42_rassf10_detection.png"))
