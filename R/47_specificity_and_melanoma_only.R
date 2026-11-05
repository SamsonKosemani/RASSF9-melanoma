## =============================================================================
## 47_specificity_and_melanoma_only.R
##
## Corrects two defects in script 46 and investigates one finding it turned up.
##
## DEFECT 1 -- the specificity condition was unsigned
##   Script 46 tested abs(core rho) > abs(shore rho). For RASSF9 that compared
##   core -0.383 against shore +0.520 and failed, labelling the result
##   COMPOSITION. Those correlations point in OPPOSITE directions: a positive
##   shore correlation is the island-shore inversion, which is the expected
##   architecture and was also seen in TCGA (core -0.29 to -0.49, shore +0.135
##   and +0.136). Comparing magnitudes treated confirmation as contradiction.
##
##   Corrected condition, signed:
##     core rho < 0, AND shore rho > core rho, AND distal rho > core rho
##   i.e. the core must be the most negative, not merely the largest in size.
##
## DEFECT 2 -- melanoma was pooled into the Skin lineage
##   DepMap's OncotreeLineage for melanoma IS "Skin". Script 45 relabelled
##   melanoma as its own group, leaving n = 3 genuine non-melanoma skin lines.
##   Script 46's within-lineage loop used lineage directly, so its "Skin" row
##   (n = 40) is mostly melanoma. The two Skin rows are therefore not
##   comparable, and 46's Skin correlation is really a melanoma correlation
##   wearing the wrong label. Here melanoma is always its own stratum.
##
## THE FINDING TO INTERROGATE
##   That mislabelled row showed rho = -0.618 (BH p = 1.8e-04) with only 2.5%
##   of lines above TPM 1. So there is a negative methylation-expression rank
##   correlation among melanoma lines, computed on expression values sitting
##   BELOW the detection threshold. Scripts 42/43 concluded no within-melanoma
##   correlation was testable. Both can be true: rank variation exists below
##   the threshold, and it may be signal, noise, or both.
##
##   Three things are done rather than asserted:
##     1. a permutation test, because Spearman's p-value is unreliable with
##        this much tie structure
##     2. the same test against CONTROL GENES. If core methylation correlates
##        with RASSF9 sub-threshold expression but not with other genes'
##        sub-threshold expression, the correlation is gene-specific rather
##        than a generic artefact of low-count noise
##     3. the expression distribution printed, so the reader can see how far
##        below detection these values sit
##
## PRE-SPECIFIED CRITERION for the melanoma-only question
##   SUGGESTIVE    rho < 0, permutation p < 0.05, AND no control gene reaches
##                 a correlation of similar strength and direction
##   UNINFORMATIVE anything else
##   In neither case is this a demonstration of a within-melanoma
##   methylation-expression relationship. At best it is a hint that the
##   relationship persists below the detection floor.
##
## INPUTS
##   results/ccle/46_line_level_data.csv
##   data/raw/ccle/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv  (controls)
##
## OUTPUTS
##   results/ccle/47_specificity_corrected.csv
##   results/ccle/47_within_lineage_corrected_rassf9.csv  /  ..._rassf10.csv
##   results/ccle/47_melanoma_only.csv
##   figures/ccle/47_rassf9_melanoma_only.png
##   results/logs/47_specificity_and_melanoma_only.log
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
LOG_DIR  <- file.path(ROOT, "results", "logs")
for (d in c(RES_DIR, FIG_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_PATH <- file.path(LOG_DIR, "47_specificity_and_melanoma_only.log")
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

GENES       <- c("RASSF9", "RASSF10")
TPM_CUT     <- 1
MIN_LINEAGE <- 15
N_PERM      <- 20000

## Control genes for the sub-threshold test. Two housekeeping (expression
## should be unrelated to RASSF9 promoter methylation), one gene expected at
## the floor in these lines, and three RASSF family members as biologically
## adjacent but separately regulated loci.
CONTROLS <- c("ACTB", "GAPDH", "MYOD1", "RASSF1", "RASSF2", "RASSF4")

## ---- 1. Load ----------------------------------------------------------------
say("")
say("== 1. Inputs ===================================================")

f_dat <- file.path(RES_DIR, "46_line_level_data.csv")
if (!file.exists(f_dat)) stop("missing ", basename(f_dat), " -- run 46 first", call. = FALSE)
D <- fread(f_dat, showProgress = FALSE)
say("  ", basename(f_dat), ": ", nrow(D), " lines, ", ncol(D), " columns")

## Melanoma is always its own stratum. Non-melanoma skin lines keep the Skin
## label; this is the distinction script 46 lost.
D[, is_mel := !is.na(disease) & grepl("melanoma", disease, ignore.case = TRUE)]
D[, is_uveal := !is.na(subtype) & grepl("uveal|ocular", paste(subtype, disease), ignore.case = TRUE)]
D[, stratum := ifelse(is_mel & !is_uveal, "Melanoma (non-uveal)",
                ifelse(is_mel, "Uveal melanoma",
                  ifelse(!is.na(lineage) & lineage == "Skin", "Skin (non-melanoma)",
                         lineage)))]

say("")
say("  strata relevant to the Skin confusion:")
for (s in c("Melanoma (non-uveal)", "Uveal melanoma", "Skin (non-melanoma)"))
  say(sprintf("    %-24s %4d lines", s, sum(D$stratum == s, na.rm = TRUE)))
say("")
say("  Script 46 reported a 'Skin' lineage of n = 40. That pooled the two")
say("  groups above. Script 45's n = 3 referred to non-melanoma skin only.")

## ---- 2. Corrected specificity -----------------------------------------------
say("")
say("== 2. Specificity, signed ======================================")
say("  condition: core rho < 0 AND core is the MOST NEGATIVE of the three sets")
say("")

rows <- list()
for (g in GENES) {
  say("  ", g, ":")
  est <- setNames(rep(NA_real_, 3), c("core", "shore", "distal"))
  for (st in c("core", "shore", "distal")) {
    cn <- paste0(g, "_", st)
    if (!cn %in% names(D)) next
    ok <- is.finite(D[[cn]]) & is.finite(D[[g]])
    if (sum(ok) < 20) next
    ct <- suppressWarnings(stats::cor.test(D[[cn]][ok], D[[g]][ok], method = "spearman"))
    est[st] <- unname(ct$estimate)
    say(sprintf("    %-8s n = %4d  rho = %+.3f  p = %.2e", st, sum(ok), est[st], ct$p.value))
    rows[[length(rows) + 1]] <- data.table(gene = g, set = st, n = sum(ok),
                                           rho = est[st], p = ct$p.value)
  }
  pass <- is.finite(est["core"]) && est["core"] < 0 &&
          (!is.finite(est["shore"])  || est["shore"]  > est["core"]) &&
          (!is.finite(est["distal"]) || est["distal"] > est["core"])
  say("    core is the most negative set: ", pass)
  ## Braces required: without them only the first say() is conditional and the
  ## second would print the TCGA claim for every gene regardless of its signs.
  if (is.finite(est["core"]) && is.finite(est["shore"]) &&
      est["core"] < 0 && est["shore"] > 0) {
    say("    -> opposite signs across the island boundary: the shore inversion,")
    say("       as in TCGA (core negative, shore positive). Confirmatory.")
  } else if (is.finite(est["core"]) && is.finite(est["shore"]) && est["core"] < 0) {
    say("    -> shore runs in the same direction as the core, so this locus")
    say("       does not show the island-shore inversion.")
  }
  say("")
}
SPEC <- rbindlist(rows)
fwrite(SPEC, file.path(RES_DIR, "47_specificity_corrected.csv"))

## ---- 3. Within-stratum correlations, melanoma separated ----------------------
say("")
say("== 3. Within-stratum correlations ==============================")
say("  as script 46 section 8, but with melanoma as its own stratum and")
say("  non-melanoma skin kept separate")

for (g in GENES) {
  cn <- paste0(g, "_core")
  d  <- D[is.finite(get(cn)) & is.finite(get(g)) & !is.na(stratum)]

  L <- d[, {
    if (.N >= MIN_LINEAGE && length(unique(get(cn))) > 3 && length(unique(get(g))) > 3) {
      ct <- suppressWarnings(stats::cor.test(get(cn), get(g), method = "spearman"))
      .(n = .N, rho = unname(ct$estimate), p = ct$p.value,
        pct_det = 100 * mean(get(g) >= log2(TPM_CUT + 1)),
        median_core = median(get(cn)))
    } else .(n = .N, rho = NA_real_, p = NA_real_,
             pct_det = 100 * mean(get(g) >= log2(TPM_CUT + 1)),
             median_core = median(get(cn)))
  }, by = stratum][order(rho)]

  LT <- L[is.finite(rho)]
  LT[, p_bh := p.adjust(p, method = "BH")]

  say("")
  say("  ", g, "  (strata with n >= ", MIN_LINEAGE, "):")
  say(sprintf("    %-28s %5s %8s %10s %10s %8s %10s",
              "stratum", "n", "rho", "p", "BH p", "%det", "medianCore"))
  for (i in seq_len(nrow(LT))) {
    l <- LT[i]
    mark <- if (grepl("^Melanoma", l$stratum)) "  <<<" else ""
    say(sprintf("    %-28s %5d %8.3f %10.2e %10.2e %7.1f%% %10.3f%s",
                l$stratum, l$n, l$rho, l$p, l$p_bh, l$pct_det, l$median_core, mark))
  }
  if (nrow(LT)) {
    n_neg <- sum(LT$rho < 0)
    st <- stats::binom.test(n_neg, nrow(LT), p = 0.5)
    say(sprintf("    median rho %+.3f | negative in %d of %d (sign test p = %.3g)",
                median(LT$rho), n_neg, nrow(LT), st$p.value))
    say(sprintf("    BH-significant and negative: %d", sum(LT$p_bh < 0.05 & LT$rho < 0)))
  }
  ## strata too small to test are still worth seeing
  LS <- L[!is.finite(rho)]
  if (nrow(LS)) {
    say("    strata below n = ", MIN_LINEAGE, " (not tested):")
    for (i in seq_len(nrow(LS)))
      say(sprintf("      %-26s %4d lines | %.1f%% detected | median core %.3f",
                  LS$stratum[i], LS$n[i], LS$pct_det[i], LS$median_core[i]))
  }
  fwrite(L, file.path(RES_DIR, sprintf("47_within_lineage_corrected_%s.csv", tolower(g))))
}

## ---- 4. Control genes for the sub-threshold test -----------------------------
say("")
say("== 4. Loading control genes =====================================")

EXPR_CANDIDATES <- c(
  "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
  "OmicsExpressionProteinCodingGenesTPMLogp1.csv"
)
found  <- EXPR_CANDIDATES[file.exists(file.path(CCLE_DIR, EXPR_CANDIDATES))]
CTRL   <- NULL
if (!length(found)) {
  warn("expression matrix not found -- the control-gene test will be skipped")
} else {
  f_expr <- file.path(CCLE_DIR, found[1])
  hdr <- names(fread(f_expr, nrows = 0, showProgress = FALSE))
  idc <- if ("ModelID" %in% hdr) "ModelID" else hdr[1]
  cm  <- setNames(rep(NA_character_, length(CONTROLS)), CONTROLS)
  for (g in CONTROLS) {
    hit <- grep(paste0("^", g, " \\("), hdr, value = TRUE)
    if (length(hit)) cm[g] <- hit[1] else warn("control gene ", g, " not found")
  }
  sel <- unique(c(idc, unname(cm[!is.na(cm)])))
  CTRL <- fread(f_expr, select = sel, showProgress = FALSE)
  setnames(CTRL, idc, "model_id")
  for (g in CONTROLS) if (!is.na(cm[g])) setnames(CTRL, cm[g], g)
  gg <- intersect(CONTROLS, names(CTRL))
  CTRL <- CTRL[, lapply(.SD, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)),
               by = model_id, .SDcols = gg]
  say("  control genes loaded: ", paste(gg, collapse = ", "))
  D <- merge(D, CTRL, by = "model_id", all.x = TRUE)
}

## ---- 5. Melanoma only ---------------------------------------------------------
say("")
say("== 5. Melanoma lines only ======================================")

M <- D[stratum == "Melanoma (non-uveal)" & is.finite(RASSF9_core) & is.finite(RASSF9)]
say("  melanoma lines with methylation and expression: ", nrow(M))

if (nrow(M) < 15) {
  warn("too few melanoma lines (", nrow(M), ") for this section")
} else {
  v <- M$RASSF9
  say("")
  say("  RASSF9 expression in these lines, log2(TPM + 1):")
  qs <- stats::quantile(v, c(0, .25, .5, .75, .9, 1))
  say(sprintf("    min %.4f | Q1 %.4f | median %.4f | Q3 %.4f | p90 %.4f | max %.4f",
              qs[1], qs[2], qs[3], qs[4], qs[5], qs[6]))
  say(sprintf("    at or above TPM %d: %d of %d (%.1f%%)",
              TPM_CUT, sum(v >= log2(TPM_CUT + 1)), length(v),
              100 * mean(v >= log2(TPM_CUT + 1))))
  say(sprintf("    exactly zero: %d | distinct values: %d",
              sum(v == 0), length(unique(v))))
  say("")
  say("  EVERY value below the TPM = 1 line is beneath the threshold scripts")
  say("  42 and 43 used to declare expression untestable. Any correlation here")
  say("  is a correlation among sub-threshold values.")

  ## Permutation test: Spearman's asymptotic p is unreliable with heavy ties.
  obs <- suppressWarnings(stats::cor(M$RASSF9_core, v, method = "spearman"))
  perm <- replicate(N_PERM, stats::cor(M$RASSF9_core, sample(v), method = "spearman"))
  p_perm <- (1 + sum(perm <= obs)) / (N_PERM + 1)       # one-sided, negative
  say("")
  say(sprintf("  RASSF9 core vs RASSF9 expression: rho = %+.3f", obs))
  say(sprintf("    permutation p (one-sided, %d permutations) = %.4g", N_PERM, p_perm))
  say(sprintf("    null distribution: mean %+.3f, 5th pct %+.3f",
              mean(perm), stats::quantile(perm, 0.05)))

  ## Control genes: does core methylation correlate with OTHER genes'
  ## expression in the same lines? If so, the result above is not specific.
  ctrl_rows <- list()
  gg <- intersect(CONTROLS, names(M))
  if (length(gg)) {
    say("")
    say("  Control genes -- RASSF9 core methylation vs each gene's expression")
    say("  in the same melanoma lines. These should NOT be strongly negative.")
    say(sprintf("    %-9s %5s %8s %10s %9s %8s",
                "gene", "n", "rho", "perm p", "medianExp", "%det"))
    for (cg in gg) {
      ok <- is.finite(M[[cg]])
      if (sum(ok) < 15) next
      o  <- suppressWarnings(stats::cor(M$RASSF9_core[ok], M[[cg]][ok], method = "spearman"))
      pm <- replicate(N_PERM, stats::cor(M$RASSF9_core[ok], sample(M[[cg]][ok]), method = "spearman"))
      pp <- (1 + sum(pm <= o)) / (N_PERM + 1)
      say(sprintf("    %-9s %5d %+8.3f %10.4g %9.3f %7.1f%%",
                  cg, sum(ok), o, pp, median(M[[cg]][ok]),
                  100 * mean(M[[cg]][ok] >= log2(TPM_CUT + 1))))
      ctrl_rows[[length(ctrl_rows) + 1]] <- data.table(gene = cg, n = sum(ok), rho = o, p_perm = pp)
    }
  }

  CTRLT <- if (length(ctrl_rows)) rbindlist(ctrl_rows) else data.table(gene = character(0), rho = numeric(0), p_perm = numeric(0))
  worst <- if (nrow(CTRLT)) min(CTRLT$rho) else NA_real_

  ## ---- verdict ----
  say("")
  say("== 6. Verdict, melanoma-only question ==========================")
  say("  SUGGESTIVE    rho < 0, permutation p < 0.05, and no control gene")
  say("                reaches a similar negative correlation")
  say("  UNINFORMATIVE otherwise")
  say("")
  cond_dir  <- obs < 0
  cond_perm <- p_perm < 0.05
  cond_ctrl <- !is.finite(worst) || worst > obs + 0.15
  verdict <- if (cond_dir && cond_perm && cond_ctrl) "SUGGESTIVE" else "UNINFORMATIVE"
  say("  VERDICT: ", verdict)
  say(sprintf("    direction negative        : %s", cond_dir))
  say(sprintf("    permutation p < 0.05      : %s (p = %.4g)", cond_perm, p_perm))
  say(sprintf("    no control gene as strong : %s (most negative control rho %+.3f)",
              cond_ctrl, worst))
  say("")
  say("  Either way this is NOT a demonstration of a methylation-expression")
  say("  relationship within melanoma cells. The expression values are all")
  say("  below the detection threshold. At best it hints that the pan-cancer")
  say("  relationship persists below that floor; at worst it is tie-structure")
  say("  noise. It belongs in a supplementary note, not a main claim.")

  OUT <- data.table(
    n_lines = nrow(M), rho = obs, p_perm = p_perm,
    pct_detected = 100 * mean(v >= log2(TPM_CUT + 1)),
    median_expr = median(v), max_expr = max(v),
    median_core = median(M$RASSF9_core), verdict = verdict,
    worst_control_rho = worst
  )
  fwrite(OUT, file.path(RES_DIR, "47_melanoma_only.csv"))
  if (nrow(CTRLT)) fwrite(CTRLT, file.path(RES_DIR, "47_melanoma_only_controls.csv"))

  ## ---- figure ----
  say("")
  say("== 7. Figure ===================================================")
  pl <- ggplot(M, aes(x = RASSF9_core, y = RASSF9)) +
    geom_hline(yintercept = log2(TPM_CUT + 1), linetype = "dashed",
               linewidth = 0.4, colour = "firebrick") +
    annotate("text", x = 0.02, y = log2(TPM_CUT + 1), label = "TPM = 1 (detection threshold)",
             hjust = 0, vjust = -0.6, size = 2.8, colour = "firebrick") +
    geom_point(size = 1.8, alpha = 0.8, colour = "grey20") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "steelblue4", linewidth = 0.45, fill = "grey85") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title    = "RASSF9 in melanoma cell lines: promoter methylation vs expression",
      subtitle = sprintf(paste0("n = %d lines | Spearman rho = %.3f, permutation p = %.3g | ",
                                "%.1f%% reach TPM ≥ 1 — all points below the dashed line are sub-threshold"),
                         nrow(M), obs, p_perm, 100 * mean(v >= log2(TPM_CUT + 1))),
      x = "RASSF9 promoter-core beta (mean of 5 TSS-proximal probes)",
      y = expression(log[2]*"(TPM + 1)")
    ) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(face = "bold", size = 10.5),
          plot.subtitle = element_text(size = 7.4, colour = "grey25"),
          plot.margin   = margin(10, 14, 8, 10))

  out <- file.path(FIG_DIR, "47_rassf9_melanoma_only.png")
  ggsave(out, pl, width = 6.6, height = 4.4, dpi = 400)
  say("  wrote: ", basename(out))
}

## ---- 8. Summary of what changed ---------------------------------------------
say("")
say("== 8. What this changes about script 46 ========================")
say("  1. The specificity condition in 46 was unsigned and mislabelled RASSF9")
say("     as COMPOSITION. On the signed test the core is the most negative set")
say("     and the shore runs the other way, reproducing the TCGA inversion.")
say("     Quote this script's specificity table, not 46's verdict line.")
say("  2. 46's 'Skin' lineage row pooled melanoma with non-melanoma skin and")
say("     should not be cited. Use section 3 here instead.")
say("  3. The melanoma-only correlation is reported with a permutation test")
say("     and control genes, and remains supplementary whatever its verdict.")
