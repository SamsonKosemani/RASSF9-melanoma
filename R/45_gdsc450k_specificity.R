## =============================================================================
## 45_gdsc450k_specificity.R
##
## Tests the observation from 44b. No new download; reads the extracted matrix
## from data/processed/44b_locus_betas.rds.
##
## WHAT 44b FOUND (medians, untested)
##   RASSF9 promoter core:  melanoma lines 0.588 | TCGA tumours 0.374 |
##                          972 other cancer lines 0.161
##   shore:                 melanoma 0.464 | other lines 0.573
##   distal:                melanoma 0.453 | other lines 0.461
##   i.e. elevation confined to the promoter core, at all 5 core probes.
##
## WHAT THIS SCRIPT ASKS
##   1. Is the core elevation statistically supported, or is it 56 lines of
##      noise against a large comparison group?
##   2. Is it confined to the core? Shore and distal are the pre-specified
##      negative controls: if they move with the core, the result is global
##      hypermethylation from culture, not locus-specific.
##   3. Is it melanoma-SPECIFIC, or shared with related lineages? A single
##      pooled "other cancer" baseline cannot tell neural-crest or skin-wide
##      from melanoma-specific. So the core score is resolved by lineage and
##      melanoma's rank among them is reported.
##
## THE KEY STATISTIC
##   core minus distal, per line. Each line is its own reference, so a line
##   that is globally hypermethylated contributes nothing to this difference.
##   It is the cleanest single number for "is the promoter core specifically
##   methylated in this line", and it is immune to the culture-drift objection
##   Gustavo would otherwise be right to raise.
##
## PRE-SPECIFIED CRITERION  (fixed before any test is run)
##   MELANOMA-SPECIFIC  core higher in melanoma, BH p < 0.05
##                      AND core-minus-distal higher in melanoma, p < 0.05
##                      AND melanoma among the top 3 lineages by core score
##   LINEAGE-SHARED     core elevated but melanoma not top-3, or other
##                      lineages at a similar level
##   NOT SUPPORTED      core difference not significant
##
##   Multiple testing: BH applied within gene, across that gene's probe family
##   (12 probes for RASSF9, 36 for RASSF10). The two genes are separate
##   families and separate tables, per project convention.
##
## INPUTS
##   data/processed/44b_locus_betas.rds
##   results/ccle/44a_sample_map.csv
##   data/raw/ccle/Model.csv
##
## OUTPUTS
##   results/ccle/45_probe_tests_rassf9.csv
##   results/ccle/45_probe_tests_rassf10.csv
##   results/ccle/45_lineage_core_rassf9.csv
##   results/ccle/45_lineage_core_rassf10.csv
##   results/ccle/45_set_tests.csv
##   figures/ccle/45_rassf9_lineage.png
##   figures/ccle/45_rassf10_lineage.png
##   results/logs/45_gdsc450k_specificity.log
## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------
if (file.exists("R/00_setup.R")) source("R/00_setup.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT     <- if (exists("PROJECT_ROOT")) PROJECT_ROOT else getwd()
CCLE_DIR <- file.path(ROOT, "data", "raw", "ccle")
RES_DIR  <- file.path(ROOT, "results", "ccle")
FIG_DIR  <- file.path(ROOT, "figures", "ccle")
PROC_DIR <- file.path(ROOT, "data", "processed")
LOG_DIR  <- file.path(ROOT, "results", "logs")
for (d in c(RES_DIR, FIG_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_PATH <- file.path(LOG_DIR, "45_gdsc450k_specificity.log")
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
MIN_LINEAGE <- 10     # lineages with fewer lines are reported but not ranked

## ---- 1. Load ----------------------------------------------------------------
say("")
say("== 1. Inputs ===================================================")

f_rds <- file.path(PROC_DIR, "44b_locus_betas.rds")
if (!file.exists(f_rds)) stop("missing ", f_rds, " -- run 44b first", call. = FALSE)
OBJ    <- readRDS(f_rds)
B      <- OBJ$beta
PROBES <- as.data.table(OBJ$probes)

f_model <- file.path(CCLE_DIR, "Model.csv")
if (!file.exists(f_model)) stop("Model.csv not found in ", CCLE_DIR, call. = FALSE)
MD <- fread(f_model, showProgress = FALSE)

say("  beta matrix: ", nrow(B), " probes x ", ncol(B), " cell lines")
say("  probes: RASSF9 ", sum(PROBES$gene == "RASSF9"),
    " | RASSF10 ", sum(PROBES$gene == "RASSF10"))

## ---- 2. Lineage for every line ----------------------------------------------
## 44b only had lineage for the melanoma lines. For the specificity question we
## need it for all 1,028, so every GDSC sample name is matched to Model.csv.
say("")
say("== 2. Assigning a lineage to all ", ncol(B), " lines ================")

strip <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

ANN <- MD[, .(match_key = strip(StrippedCellLineName),
              lineage   = OncotreeLineage,
              disease   = OncotreePrimaryDisease,
              subtype   = OncotreeSubtype)]
ANN <- ANN[!duplicated(match_key)]

LINES <- data.table(geo_name = colnames(B), match_key = strip(colnames(B)))
LINES <- merge(LINES, ANN, by = "match_key", all.x = TRUE)

## Restore the matrix column order with an explicit match rather than a keyed
## subset: data.table drops the key after some subsetting operations, and a
## silently misordered annotation would mislabel every line's lineage.
LINES <- LINES[match(colnames(B), LINES$geo_name)]
stopifnot(identical(LINES$geo_name, colnames(B)))

## Named lookups, used instead of keyed indexing for the same reason.
lin_of <- setNames(LINES$lineage, LINES$geo_name)

n_lin <- sum(!is.na(LINES$lineage))
say("  lines matched to a DepMap lineage: ", n_lin, " of ", ncol(B),
    sprintf(" (%.1f%%)", 100 * n_lin / ncol(B)))
if (n_lin < 0.5 * ncol(B))
  warn("fewer than half the lines got a lineage -- the lineage analysis is ",
       "weakly powered and the pooled comparison should be preferred")

LINES[, is_mel := !is.na(disease) & grepl("melanoma", disease, ignore.case = TRUE)]
LINES[, is_uveal := !is.na(subtype) &
        grepl("uveal|ocular", paste(subtype, disease), ignore.case = TRUE)]

grp_mel   <- LINES$geo_name[LINES$is_mel & !LINES$is_uveal]
grp_other <- LINES$geo_name[!LINES$is_mel & !is.na(LINES$lineage)]

say("  melanoma (non-uveal)            : ", length(grp_mel))
say("  other lines WITH a lineage      : ", length(grp_other))
say("  unmatched (excluded from tests) : ", ncol(B) - length(grp_mel) - length(grp_other))
say("")
say("  NOTE: 44b used all 972 non-melanoma columns as the baseline. Here the")
say("  baseline is restricted to lines with a known lineage, so the pooled")
say("  numbers may differ slightly from 44b. That is intentional -- the")
say("  lineage test requires the annotation.")

if (length(grp_mel) < 20) warn("only ", length(grp_mel), " melanoma lines -- underpowered")

## ---- 3. Per-line set scores -------------------------------------------------
say("")
say("== 3. Per-line scores ==========================================")

set_score <- function(g, st, cols) {
  pr <- PROBES[gene == g & set == st & probe %in% rownames(B), probe]
  if (!length(pr)) return(rep(NA_real_, length(cols)))
  colMeans(B[pr, cols, drop = FALSE], na.rm = TRUE)
}

SC <- list()
for (g in c("RASSF9", "RASSF10")) {
  cols <- c(grp_mel, grp_other)
  SC[[g]] <- data.table(
    geo_name = cols,
    group    = c(rep("Melanoma", length(grp_mel)), rep("Other", length(grp_other))),
    lineage  = unname(lin_of[cols]),
    core     = set_score(g, "tss_proximal", cols),
    shore    = set_score(g, "shore",        cols),
    distal   = set_score(g, "distal",       cols)
  )
  SC[[g]][, core_minus_distal := core - distal]
  say("  ", g, ": scored ", nrow(SC[[g]]), " lines")
}

## ---- 4. Set-level tests -----------------------------------------------------
say("")
say("== 4. Melanoma vs other lineages, by probe set =================")
say("  core is the test; shore and distal are pre-specified negative controls")
say("  core_minus_distal normalises each line against itself")
say("")

set_rows <- list()
for (g in c("RASSF9", "RASSF10")) {
  D <- SC[[g]]
  say("  ", g, ":")
  say(sprintf("    %-18s %8s %8s %9s %10s   %s",
              "measure", "melMed", "othMed", "shift", "p", "n mel / n oth"))
  for (m in c("core", "shore", "distal", "core_minus_distal")) {
    a <- D[group == "Melanoma"][[m]]; a <- a[is.finite(a)]
    b <- D[group == "Other"   ][[m]]; b <- b[is.finite(b)]
    if (length(a) < 5 || length(b) < 5) { say("    ", m, ": too few values"); next }
    wt <- suppressWarnings(stats::wilcox.test(a, b, conf.int = TRUE))
    say(sprintf("    %-18s %8.3f %8.3f %+9.3f %10.2e   %d / %d",
                m, median(a), median(b), unname(wt$estimate), wt$p.value,
                length(a), length(b)))
    set_rows[[length(set_rows) + 1]] <- data.table(
      gene = g, measure = m, n_mel = length(a), n_oth = length(b),
      med_mel = median(a), med_oth = median(b),
      shift = unname(wt$estimate), lo = wt$conf.int[1], hi = wt$conf.int[2],
      p = wt$p.value)
  }
  say("")
}
SETT <- rbindlist(set_rows)
fwrite(SETT, file.path(RES_DIR, "45_set_tests.csv"))

say("  'shift' is the Hodges-Lehmann estimate of the median difference")
say("  (melanoma minus other), with a 95% CI in the output CSV.")

## ---- 5. Per-probe tests, BH within gene -------------------------------------
say("")
say("== 5. Per-probe tests (BH within each gene's probe family) =====")

for (g in c("RASSF9", "RASSF10")) {
  pp <- PROBES[gene == g & probe %in% rownames(B)][order(dist_tss)]
  rows <- rbindlist(lapply(seq_len(nrow(pp)), function(i) {
    pr <- pp$probe[i]
    a <- B[pr, grp_mel];   a <- a[is.finite(a)]
    b <- B[pr, grp_other]; b <- b[is.finite(b)]
    if (length(a) < 5 || length(b) < 5)
      return(data.table(gene = g, probe = pr, dist_tss = pp$dist_tss[i],
                        set = pp$set[i], n_mel = length(a), n_oth = length(b),
                        med_mel = NA_real_, med_oth = NA_real_,
                        shift = NA_real_, p = NA_real_))
    wt <- suppressWarnings(stats::wilcox.test(a, b, conf.int = TRUE))
    data.table(gene = g, probe = pr, dist_tss = pp$dist_tss[i], set = pp$set[i],
               n_mel = length(a), n_oth = length(b),
               med_mel = median(a), med_oth = median(b),
               shift = unname(wt$estimate), p = wt$p.value)
  }))
  rows[, p_bh := p.adjust(p, method = "BH")]

  say("")
  say("  ", g, "  (family of ", nrow(rows), " probes):")
  say(sprintf("    %-12s %8s %-13s %8s %8s %+8s %10s %10s",
              "probe", "TSSdist", "set", "melMed", "othMed", "shift", "p", "BH p"))
  for (i in seq_len(nrow(rows))) {
    r <- rows[i]
    say(sprintf("    %-12s %+8d %-13s %8.3f %8.3f %+8.3f %10.2e %10.2e",
                r$probe, r$dist_tss, r$set, r$med_mel, r$med_oth, r$shift, r$p, r$p_bh))
  }
  for (st in c("tss_proximal", "shore", "distal")) {
    ss <- rows[set == st & is.finite(p_bh)]
    if (nrow(ss))
      say(sprintf("    %-13s %d of %d probes at BH p < 0.05 | mean shift %+.3f",
                  st, sum(ss$p_bh < 0.05), nrow(ss), mean(ss$shift)))
  }
  fwrite(rows, file.path(RES_DIR, sprintf("45_probe_tests_%s.csv", tolower(g))))
}

## ---- 6. Lineage-resolved core score -----------------------------------------
say("")
say("== 6. Core score by lineage ====================================")
say("  A pooled baseline cannot distinguish melanoma-specific from")
say("  neural-crest-wide or skin-wide. If melanoma does not rank near the top")
say("  here, the pooled result is a lineage effect rather than a melanoma one.")

ranks <- list()
for (g in c("RASSF9", "RASSF10")) {
  D <- SC[[g]][is.finite(core)]
  D[, lin := ifelse(group == "Melanoma", "Melanoma (non-uveal)", lineage)]
  L <- D[, .(n = .N, median_core = median(core), iqr = IQR(core)), by = lin][order(-median_core)]

  L_rank <- L[n >= MIN_LINEAGE]
  L_rank[, rank := seq_len(.N)]
  mel_rank <- L_rank[lin == "Melanoma (non-uveal)", rank]
  n_rank   <- nrow(L_rank)

  say("")
  say("  ", g, " promoter-core score by lineage (n >= ", MIN_LINEAGE, " ranked):")
  say(sprintf("    %-34s %5s %12s %8s", "lineage", "n", "medianCore", "IQR"))
  for (i in seq_len(nrow(L))) {
    l <- L[i]
    mark <- if (l$lin == "Melanoma (non-uveal)") "  <<<" else
            if (l$n < MIN_LINEAGE) "   (not ranked)" else ""
    say(sprintf("    %-34s %5d %12.3f %8.3f%s", l$lin, l$n, l$median_core, l$iqr, mark))
  }
  if (length(mel_rank))
    say("    melanoma rank: ", mel_rank, " of ", n_rank, " lineages with n >= ", MIN_LINEAGE)
  else
    say("    melanoma not rankable")

  ranks[[g]] <- list(rank = if (length(mel_rank)) mel_rank else NA_integer_, n = n_rank)
  fwrite(L, file.path(RES_DIR, sprintf("45_lineage_core_%s.csv", tolower(g))))
}

## ---- 7. Verdict -------------------------------------------------------------
say("")
say("== 7. Verdict against the pre-set criterion ====================")

for (g in c("RASSF9", "RASSF10")) {
  pt   <- fread(file.path(RES_DIR, sprintf("45_probe_tests_%s.csv", tolower(g))),
                showProgress = FALSE)
  core_sig <- pt[set == "tss_proximal" & is.finite(p_bh) & p_bh < 0.05 & shift > 0, .N]
  core_n   <- pt[set == "tss_proximal" & is.finite(p_bh), .N]
  cmd      <- SETT[gene == g & measure == "core_minus_distal"]
  r        <- ranks[[g]]

  cond_core <- core_sig > 0
  cond_cmd  <- nrow(cmd) == 1 && is.finite(cmd$p) && cmd$p < 0.05 && cmd$shift > 0
  cond_rank <- is.finite(r$rank) && r$rank <= 3

  verdict <- if (!cond_core) "NOT SUPPORTED" else
             if (cond_core && cond_cmd && cond_rank) "MELANOMA-SPECIFIC" else "LINEAGE-SHARED"

  say("")
  say("  ", g, ": ", verdict)
  say(sprintf("    core probes raised at BH p < 0.05 : %d of %d", core_sig, core_n))
  if (nrow(cmd) == 1)
    say(sprintf("    core-minus-distal shift           : %+.3f, p = %.2e", cmd$shift, cmd$p))
  say(sprintf("    melanoma rank among lineages      : %s of %d",
              ifelse(is.finite(r$rank), as.character(r$rank), "n/a"), r$n))

  sh <- pt[set == "shore"  & is.finite(p_bh) & p_bh < 0.05, .N]
  ds <- pt[set == "distal" & is.finite(p_bh) & p_bh < 0.05, .N]
  say(sprintf("    control sets also significant     : shore %d, distal %d", sh, ds))
  if (ds > 0.5 * pt[set == "distal", .N])
    warn(g, ": more than half the distal probes also differ -- part of this may be ",
         "global methylation difference rather than promoter-specific")
}

## ---- 8. Figures --------------------------------------------------------------
say("")
say("== 8. Figures ==================================================")

for (g in c("RASSF9", "RASSF10")) {
  D <- SC[[g]][is.finite(core)]
  D[, lin := ifelse(group == "Melanoma", "Melanoma (non-uveal)", lineage)]
  keep <- D[, .N, by = lin][N >= MIN_LINEAGE, lin]
  D <- D[lin %in% keep]
  ord <- D[, .(m = median(core)), by = lin][order(m)]
  D[, lin := factor(lin, levels = ord$lin)]
  D[, hl := lin == "Melanoma (non-uveal)"]

  p <- ggplot(D, aes(x = lin, y = core)) +
    geom_boxplot(aes(fill = hl), outlier.size = 0.5, linewidth = 0.3, width = 0.65) +
    scale_fill_manual(values = c(`FALSE` = "grey88", `TRUE` = "firebrick"), guide = "none") +
    coord_flip() +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title    = paste0(g, " promoter-core methylation across cancer lineages"),
      subtitle = sprintf(paste0("GSE68379 (HM450) | mean beta over %d TSS-proximal probes | ",
                                "lineages with n ≥ %d | melanoma in red"),
                         PROBES[gene == g & set == "tss_proximal" & probe %in% rownames(B), .N],
                         MIN_LINEAGE),
      x = NULL, y = "promoter-core beta (per line)"
    ) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 7.6, colour = "grey25"),
          plot.margin   = margin(10, 14, 8, 10))

  out <- file.path(FIG_DIR, sprintf("45_%s_lineage.png", tolower(g)))
  ggsave(out, p, width = 7.2, height = 5.4, dpi = 400)
  say("  wrote: ", basename(out))
}

## ---- 9. What this still cannot say -------------------------------------------
say("")
say("== 9. What this still cannot say ===============================")
say("  1. Nothing here is a methylation-expression relationship. It is a")
say("     comparison of methylation levels between groups of cell lines.")
say("     The expression side remains at the floor in melanoma (0 of 39 lines")
say("     with both measurements reached TPM >= 1).")
say("  2. Cell lines are not tumours. The core-minus-distal statistic removes")
say("     per-line global differences, which is the main culture objection, but")
say("     it does not make a cultured line a melanoma in a patient.")
say("  3. Lineage is an observational grouping. Melanoma lines differ from")
say("     other lineages in many ways besides being melanoma.")
say("  4. A melanoma-specific methylated-and-silent promoter is consistent")
say("     with methylation-associated silencing. It does not establish that")
say("     methylation CAUSED the silencing; that needs demethylation in vitro.")
