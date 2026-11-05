## =============================================================================
## 43_ccle_gate_robustness.R
##
## Hardens the STOP verdict from script 42 against the two weaknesses in it.
##
## WEAKNESS 1 -- arbitrary profile choice
##   The 26Q3 expression matrix has one row per sequencing profile and, contrary
##   to its own file description, carries no is_default_entry column. Script 42
##   kept the first row per model and dropped 184 duplicates arbitrarily.
##
##   Fix: collapse profiles per model by MAXIMUM rather than by first. The max
##   is by construction >= whichever profile DepMap would have designated, so
##   this is the most generous possible reading of RASSF9 expression and makes
##   OmicsProfiles.csv unnecessary. If the gate still says STOP under the max,
##   no choice of profile could have changed it.
##
## WEAKNESS 2 -- pooled subtypes
##   Script 42 separated uveal only. Acral (6) and vulvovaginal mucosal (3)
##   lines stayed pooled with cutaneous.
##
##   Fix: report every subtype separately. These arise on non-UV-exposed sites
##   and are biologically distinct, and unlike the TCGA metadata the cell-line
##   panel labels them explicitly.
##
## This script changes no conclusion on its own. It either confirms the STOP
## under the most favourable assumptions available, or it reveals that the
## profile collapse mattered -- in which case script 42's result needs redoing
## with OmicsProfiles.csv rather than defending.
##
## INPUTS   data/raw/ccle/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
##          data/raw/ccle/Model.csv
##          (both already present from script 42 -- no new download)
##
## OUTPUTS  results/ccle/43_profile_collapse_comparison.csv
##          results/ccle/43_subtype_breakdown_rassf9.csv
##          results/ccle/43_subtype_breakdown_rassf10.csv
##          figures/ccle/43_rassf9_by_subtype.png
##          figures/ccle/43_rassf10_by_subtype.png
##          results/logs/43_ccle_gate_robustness.log
## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------
if (file.exists("R/00_setup.R")) source("R/00_setup.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT    <- if (exists("PROJECT_ROOT")) PROJECT_ROOT else getwd()
RAW_DIR <- file.path(ROOT, "data", "raw", "ccle")
RES_DIR <- file.path(ROOT, "results", "ccle")
FIG_DIR <- file.path(ROOT, "figures", "ccle")
LOG_DIR <- file.path(ROOT, "results", "logs")
for (d in c(RES_DIR, FIG_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_PATH <- file.path(LOG_DIR, "43_ccle_gate_robustness.log")
.logcon  <- file(LOG_PATH, open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = .logcon) }
WARN <- character(0)
warn <- function(...) { m <- paste0(...); WARN <<- c(WARN, m); say("  [WARN] ", m) }
on.exit({
  if (length(WARN)) { say(""); say("WARNINGS RAISED (", length(WARN), "):"); for (w in WARN) say("  - ", w) }
  else say("  no warnings raised")
  say(""); say("Full log: ", LOG_PATH); close(.logcon)
}, add = TRUE)

set.seed(1234)

GENES   <- c("RASSF9", "RASSF10")
TPM_CUT <- 1                 # TPM >= 1  <=>  log2(TPM+1) >= 1

## Same pre-set gate as script 42. Not re-derived, so the comparison is honest.
GATE_DET_HI <- 0.70; GATE_DET_LO <- 0.40
GATE_SD_HI  <- 0.50; GATE_SD_LO  <- 0.30
GATE_TIE_MAX <- 0.40

verdict_of <- function(pct_det, sd_v, pct_tied) {
  if (pct_det / 100 >= GATE_DET_HI && sd_v >= GATE_SD_HI && pct_tied / 100 < GATE_TIE_MAX) "PROCEED"
  else if (pct_det / 100 < GATE_DET_LO || sd_v < GATE_SD_LO) "STOP"
  else "MARGINAL"
}

## ---- 1. Inputs --------------------------------------------------------------
say("")
say("== 1. Inputs ===================================================")

EXPR_CANDIDATES <- c(
  "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
  "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
  "OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv"
)
found   <- EXPR_CANDIDATES[file.exists(file.path(RAW_DIR, EXPR_CANDIDATES))]
f_expr  <- if (length(found)) file.path(RAW_DIR, found[1]) else NA_character_
f_model <- file.path(RAW_DIR, "Model.csv")

if (is.na(f_expr) || !file.exists(f_model))
  stop("inputs from script 42 not found in ", RAW_DIR, call. = FALSE)

say(sprintf("  %-52s %8.1f MB", basename(f_expr),  file.size(f_expr)  / 1024^2))
say(sprintf("  %-52s %8.1f MB", basename(f_model), file.size(f_model) / 1024^2))

## ---- 2. Gene columns --------------------------------------------------------
say("")
say("== 2. Gene and ID columns ======================================")

hdr <- names(fread(f_expr, nrows = 0, showProgress = FALSE))
id_col <- if ("ModelID" %in% hdr) "ModelID" else hdr[1]
seq_col <- if ("SequencingID" %in% hdr) "SequencingID" else NA_character_

col_map <- setNames(rep(NA_character_, length(GENES)), GENES)
for (g in GENES) {
  hit <- grep(paste0("^", g, " \\("), hdr, value = TRUE)
  if (length(hit) >= 1) col_map[g] <- hit[1] else warn(g, " not found")
}
say("  model ID column: '", id_col, "'")
say("  profile column : ", ifelse(is.na(seq_col), "none", paste0("'", seq_col, "'")))
for (g in GENES) say(sprintf("  %-8s -> %s", g, ifelse(is.na(col_map[g]), "NOT FOUND", col_map[g])))
if (is.na(col_map["RASSF9"])) stop("RASSF9 column not found", call. = FALSE)

sel <- unique(c(id_col, if (!is.na(seq_col)) seq_col, unname(col_map[!is.na(col_map)])))
EX  <- fread(f_expr, select = sel, showProgress = FALSE)
setnames(EX, id_col, "model_id")
for (g in GENES) if (!is.na(col_map[g])) setnames(EX, col_map[g], g)

say("  rows (profiles) read: ", nrow(EX))
say("  distinct models      : ", length(unique(EX$model_id)))

## ---- 3. How much does the profile choice matter? ----------------------------
say("")
say("== 3. Profile collapse: first vs max vs mean ===================")

dup <- EX[, .N, by = model_id]
say("  models with a single profile : ", sum(dup$N == 1))
say("  models with >1 profile       : ", sum(dup$N > 1))
say("  largest number of profiles for one model: ", max(dup$N))

collapse <- function(dt, g) {
  dt[, .(first = .SD[[1]][1],
         max   = max(.SD[[1]], na.rm = TRUE),
         mean  = mean(.SD[[1]], na.rm = TRUE),
         n_prof = .N),
     by = model_id, .SDcols = g]
}

CMP <- list()
for (g in GENES) {
  if (is.na(col_map[g])) next
  CMP[[g]] <- collapse(EX, g)
  setnames(CMP[[g]], c("first", "max", "mean"), paste0(g, "_", c("first", "max", "mean")))
}

## ---- 4. Melanoma annotation, with every subtype separated -------------------
say("")
say("== 4. Melanoma lines by subtype ================================")

MD <- fread(f_model, showProgress = FALSE)
need <- c("ModelID", "StrippedCellLineName", "OncotreePrimaryDisease", "OncotreeSubtype")
miss <- setdiff(need, names(MD))
if (length(miss)) warn("Model.csv is missing: ", paste(miss, collapse = ", "))

MEL <- MD[grepl("melanoma", OncotreePrimaryDisease, ignore.case = TRUE),
          .(model_id = ModelID,
            name     = StrippedCellLineName,
            disease  = OncotreePrimaryDisease,
            subtype  = OncotreeSubtype)]

## Classify explicitly. Anything unmatched is kept visible rather than folded
## into cutaneous -- the error that produced the overstated site verdict in
## script 40.
cls <- function(sub, dis) {
  s <- tolower(paste(sub, dis))
  fcase(
    grepl("uveal|ocular", s),        "Uveal",
    grepl("acral", s),               "Acral",
    grepl("mucosal|vulva|vagina", s),"Mucosal",
    grepl("cutaneous", s),           "Cutaneous",
    default = "Melanoma NOS"
  )
}
MEL[, group := cls(subtype, disease)]

say("  melanoma models in annotation: ", nrow(MEL))
tb <- MEL[, .N, by = group][order(-N)]
for (i in seq_len(nrow(tb))) say(sprintf("    %-16s %4d", tb$group[i], tb$N[i]))
say("")
say("  'Melanoma NOS' is skin melanoma of unstated subtype. It is NOT verified")
say("  to exclude acral or mucosal cases, so it is reported as its own group")
say("  rather than counted as cutaneous.")

## ---- 5. Gate under the max ---------------------------------------------------
say("")
say("== 5. Gate under each collapse rule ============================")
say("  The MAX column is the decisive one: it is >= whichever profile DepMap")
say("  would have designated, so it is the most generous possible reading.")
say("")

describe <- function(v) {
  v  <- v[is.finite(v)]
  fl <- min(c(v, 0))
  list(n = length(v),
       pct_det  = 100 * mean(v >= log2(TPM_CUT + 1)),
       pct_zero = 100 * mean(v <= fl + 1e-9),
       median   = median(v), sd = sd(v),
       max      = max(v),
       pct_tied = 100 * (1 - length(unique(round(v, 6))) / max(length(v), 1)))
}

ROWS <- list()
for (g in GENES) {
  if (is.na(col_map[g])) next
  M <- merge(MEL, CMP[[g]], by = "model_id")
  PRI <- M[group != "Uveal"]          # same primary set definition as script 42

  say("  ", g, "  (non-uveal melanoma lines with expression, n = ", nrow(PRI), ")")
  say(sprintf("    %-7s %5s %8s %8s %8s %7s %7s %8s   %s",
              "rule", "n", "%TPM>1", "%floor", "median", "SD", "max", "%tied", "verdict"))
  for (rule in c("first", "max", "mean")) {
    d <- describe(PRI[[paste0(g, "_", rule)]])
    v <- verdict_of(d$pct_det, d$sd, d$pct_tied)
    say(sprintf("    %-7s %5d %7.1f%% %7.1f%% %8.3f %7.3f %7.2f %7.1f%%   %s",
                rule, d$n, d$pct_det, d$pct_zero, d$median, d$sd, d$max, d$pct_tied, v))
    ROWS[[length(ROWS) + 1]] <- data.table(gene = g, rule = rule, n = d$n,
                                           pct_det = d$pct_det, pct_floor = d$pct_zero,
                                           median = d$median, sd = d$sd, max = d$max,
                                           pct_tied = d$pct_tied, verdict = v)
  }

  ## Did the collapse rule change any individual line materially?
  dd <- PRI[[paste0(g, "_max")]] - PRI[[paste0(g, "_first")]]
  n_moved <- sum(dd > 0.5, na.rm = TRUE)
  n_cross <- sum(PRI[[paste0(g, "_first")]] < 1 & PRI[[paste0(g, "_max")]] >= 1, na.rm = TRUE)
  say(sprintf("    lines whose value rose > 0.5 log2 units under max: %d", n_moved))
  say(sprintf("    lines crossing the TPM >= %d threshold under max  : %d", TPM_CUT, n_cross))
  if (n_cross > 0)
    warn(g, ": ", n_cross, " line(s) cross the detection threshold under the max rule -- ",
         "quote the max-rule numbers, not script 42's, in anything sent out")
  say("")
}

RES <- rbindlist(ROWS)
fwrite(RES, file.path(RES_DIR, "43_profile_collapse_comparison.csv"))

## ---- 6. By subtype -----------------------------------------------------------
say("")
say("== 6. Detection by subtype (max rule) ==========================")

for (g in GENES) {
  if (is.na(col_map[g])) next
  M <- merge(MEL, CMP[[g]], by = "model_id")
  vcol <- paste0(g, "_max")

  BY <- M[, .(n = .N,
              n_det   = sum(get(vcol) >= log2(TPM_CUT + 1), na.rm = TRUE),
              pct_det = 100 * mean(get(vcol) >= log2(TPM_CUT + 1), na.rm = TRUE),
              median  = median(get(vcol), na.rm = TRUE),
              max     = max(get(vcol), na.rm = TRUE)),
          by = group][order(-n)]

  say("  ", g, ":")
  say(sprintf("    %-16s %5s %6s %8s %8s %7s", "group", "n", "det", "%det", "median", "max"))
  for (i in seq_len(nrow(BY)))
    say(sprintf("    %-16s %5d %6d %7.1f%% %8.3f %7.2f",
                BY$group[i], BY$n[i], BY$n_det[i], BY$pct_det[i], BY$median[i], BY$max[i]))

  small <- BY[n < 10]
  if (nrow(small))
    say("    NOTE: ", paste(small$group, collapse = ", "),
        " have fewer than 10 lines -- percentages there are indicative only.")
  say("")

  fwrite(BY, file.path(RES_DIR, sprintf("43_subtype_breakdown_%s.csv", tolower(g))))
}

## ---- 7. Figures --------------------------------------------------------------
say("")
say("== 7. Figures ==================================================")

for (g in GENES) {
  if (is.na(col_map[g])) next
  M <- merge(MEL, CMP[[g]], by = "model_id")
  vcol <- paste0(g, "_max")
  D <- data.table(group = M$group, v = M[[vcol]])[is.finite(v)]

  ## order groups by size so the sparse subtypes do not read as equals
  ord <- D[, .N, by = group][order(-N)]
  D[, group := factor(group, levels = ord$group)]
  lab <- setNames(sprintf("%s\n(n = %d)", ord$group, ord$N), ord$group)

  p <- ggplot(D, aes(x = group, y = v)) +
    geom_hline(yintercept = log2(TPM_CUT + 1), linetype = "dashed",
               linewidth = 0.4, colour = "firebrick") +
    geom_jitter(width = 0.18, height = 0, size = 1.3, alpha = 0.65, colour = "grey25") +
    stat_summary(fun = median, geom = "crossbar", width = 0.45,
                 linewidth = 0.3, colour = "steelblue4", fatten = 0) +
    annotate("text", x = 0.55, y = log2(TPM_CUT + 1), label = "TPM = 1",
             hjust = 0, vjust = -0.6, size = 3, colour = "firebrick") +
    scale_x_discrete(labels = lab) +
    labs(
      title    = paste0(g, " expression in melanoma cell lines, by subtype"),
      subtitle = paste0(basename(f_expr),
                        " | max across sequencing profiles per model | blue bar = median"),
      x = NULL, y = expression(log[2]*"(TPM + 1)")
    ) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8, colour = "grey25"),
          plot.margin   = margin(10, 14, 8, 10))

  out <- file.path(FIG_DIR, sprintf("43_%s_by_subtype.png", tolower(g)))
  ggsave(out, p, width = 6.8, height = 4.2, dpi = 400)
  say("  wrote: ", basename(out))
}

## ---- 8. What this settles ----------------------------------------------------
say("")
say("== 8. What this settles ========================================")
say("  If the max rule returns STOP for RASSF9, then no choice of sequencing")
say("  profile could have produced a testable correlation, and the arbitrary")
say("  duplicate handling in script 42 is no longer a caveat worth stating.")
say("")
say("  If the max rule returns MARGINAL or PROCEED, script 42's verdict was an")
say("  artefact of keeping the first profile. In that case do NOT defend it --")
say("  download OmicsProfiles.csv, use IsDefaultEntryForModel, and redo the")
say("  gate properly before anything goes to Gustavo.")
say("")
say("  WHAT THIS CANNOT SAY: nothing here addresses whether the RASSF9")
say("  promoter is methylated in these lines. Low expression with uniform")
say("  hypermethylation and low expression with an unmethylated promoter look")
say("  identical in this script and mean opposite things. That is step 2.")
