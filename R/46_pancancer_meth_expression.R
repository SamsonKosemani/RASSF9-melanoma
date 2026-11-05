## =============================================================================
## 46_pancancer_meth_expression.R
##
## The methylation-expression test Gustavo asked for, run where it is actually
## computable.
##
## WHY NOT IN MELANOMA ALONE
##   Scripts 42/43: RASSF9 is at the detection floor in melanoma lines (0 of 39
##   lines with both measurements reach TPM >= 1). No variance, no correlation.
##   Scripts 44b/45: the RASSF9 promoter core is methylated in melanoma lines
##   (0.588) and unmethylated in most other lineages (0.183), 5 of 5 core probes
##   BH-significant with 0 of 2 shore and 0 of 5 distal -- so melanoma sits at
##   the methylated-and-silent extreme of a range that other lineages span.
##
##   That range is what makes the test possible. Across the full panel both
##   methylation and expression vary, so the question becomes:
##
##     Does RASSF9 promoter-core methylation track RASSF9 expression in cancer
##     cells, with melanoma as the extreme case?
##
##   This is not "within melanoma cells", and it must not be presented as if it
##   were. It is the strongest version of the question these data can support.
##
## THE CONFOUND THAT MATTERS MOST
##   Lineages differ in both methylation and expression, so a pooled
##   correlation could be pure composition -- a between-lineage effect with no
##   within-lineage relationship at all (Simpson's paradox). The analysis
##   therefore does three things, in increasing strictness:
##     1. pooled Spearman
##     2. lineage as a covariate in a linear model
##     3. Spearman computed SEPARATELY within each lineage, then summarised
##   Only (3) can show a relationship that is not composition.
##
## ZERO INFLATION
##   Expression is log2(TPM+1) and many lines sit exactly at 0, which attenuates
##   Spearman and inflates ties. So detection (TPM >= 1) is also modelled by
##   logistic regression, which handles the floor properly. Both are reported.
##
## PRE-SPECIFIED CRITERION  (fixed before any number is seen)
##   SUPPORTED     pooled core rho < 0 at p < 0.05
##                 AND |core rho| > |shore rho| and > |distal rho|
##                 AND lineage-adjusted beta keeps > 50% of its unadjusted value
##                     at p < 0.05
##                 AND within-lineage: median rho < 0 and a majority of
##                     lineages (n >= 15) negative
##   COMPOSITION   pooled correlation significant but fails the within-lineage
##                 or adjustment conditions -- i.e. a between-lineage effect
##   NOT SUPPORTED pooled core rho not significant, or positive
##
##   Copy number enters as a further covariate, per Gustavo's request.
##
## INPUTS
##   data/processed/44b_locus_betas.rds
##   data/raw/ccle/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv
##   data/raw/ccle/OmicsCNGeneWGS.csv
##   data/raw/ccle/Model.csv
##
## OUTPUTS
##   results/ccle/46_pooled_tests.csv
##   results/ccle/46_models_rassf9.csv   /  46_models_rassf10.csv
##   results/ccle/46_within_lineage_rassf9.csv  /  ..._rassf10.csv
##   results/ccle/46_line_level_data.csv        (underlying data for Gustavo)
##   figures/ccle/46_rassf9_scatter.png  /  46_rassf10_scatter.png
##   results/logs/46_pancancer_meth_expression.log
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

LOG_PATH <- file.path(LOG_DIR, "46_pancancer_meth_expression.log")
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
MIN_LINEAGE <- 15     # minimum lines for a within-lineage correlation

## ---- 1. Methylation ----------------------------------------------------------
say("")
say("== 1. Methylation scores =======================================")

f_rds <- file.path(PROC_DIR, "44b_locus_betas.rds")
if (!file.exists(f_rds)) stop("missing ", f_rds, " -- run 44b first", call. = FALSE)
OBJ    <- readRDS(f_rds)
B      <- OBJ$beta
PROBES <- as.data.table(OBJ$probes)

set_score <- function(g, st) {
  pr <- PROBES[gene == g & set == st & probe %in% rownames(B), probe]
  if (!length(pr)) return(rep(NA_real_, ncol(B)))
  colMeans(B[pr, , drop = FALSE], na.rm = TRUE)
}

MET <- data.table(geo_name = colnames(B))
for (g in GENES) {
  MET[[paste0(g, "_core")]]   <- set_score(g, "tss_proximal")
  MET[[paste0(g, "_shore")]]  <- set_score(g, "shore")
  MET[[paste0(g, "_distal")]] <- set_score(g, "distal")
}
say("  methylation for ", nrow(MET), " cell lines")
for (g in GENES)
  say(sprintf("    %-8s core: median %.3f | range %.3f - %.3f",
              g, median(MET[[paste0(g, "_core")]], na.rm = TRUE),
              min(MET[[paste0(g, "_core")]], na.rm = TRUE),
              max(MET[[paste0(g, "_core")]], na.rm = TRUE)))

## ---- 2. Annotation ----------------------------------------------------------
say("")
say("== 2. Annotation and name mapping ==============================")

f_model <- file.path(CCLE_DIR, "Model.csv")
if (!file.exists(f_model)) stop("Model.csv not found", call. = FALSE)
MD <- fread(f_model, showProgress = FALSE)

strip <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

ANN <- MD[, .(match_key = strip(StrippedCellLineName), model_id = ModelID,
              name = StrippedCellLineName, lineage = OncotreeLineage,
              disease = OncotreePrimaryDisease, subtype = OncotreeSubtype)]
ANN <- ANN[!duplicated(match_key)]

MET[, match_key := strip(geo_name)]
DT <- merge(MET, ANN, by = "match_key", all.x = TRUE)
DT[, is_mel := !is.na(disease) & grepl("melanoma", disease, ignore.case = TRUE)]
DT[, is_uveal := !is.na(subtype) & grepl("uveal|ocular", paste(subtype, disease), ignore.case = TRUE)]
DT[, group := ifelse(is_mel & !is_uveal, "Melanoma", ifelse(is_mel, "Uveal", "Other"))]

say("  lines with a DepMap model ID: ", sum(!is.na(DT$model_id)), " of ", nrow(DT))
say("  melanoma (non-uveal)        : ", sum(DT$group == "Melanoma"))

## ---- 3. Expression ----------------------------------------------------------
say("")
say("== 3. Expression ===============================================")

EXPR_CANDIDATES <- c(
  "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
  "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
  "OmicsExpressionTPMLogp1HumanProteinCodingGenesStranded.csv"
)
found  <- EXPR_CANDIDATES[file.exists(file.path(CCLE_DIR, EXPR_CANDIDATES))]
if (!length(found)) stop("no expression matrix found in ", CCLE_DIR, call. = FALSE)
f_expr <- file.path(CCLE_DIR, found[1])
say("  ", basename(f_expr))

## Locate gene columns and the model ID column, as in scripts 42/43.
grab_genes <- function(path, genes, label) {
  hdr <- names(fread(path, nrows = 0, showProgress = FALSE))
  idc <- if ("ModelID" %in% hdr) "ModelID" else hdr[1]
  cm  <- setNames(rep(NA_character_, length(genes)), genes)
  for (g in genes) {
    hit <- grep(paste0("^", g, " \\("), hdr, value = TRUE)
    if (!length(hit)) hit <- grep(paste0("^", g, "$"), hdr, value = TRUE)
    if (length(hit)) cm[g] <- hit[1] else warn(g, " not found in ", label)
  }
  sel <- unique(c(idc, unname(cm[!is.na(cm)])))
  X <- fread(path, select = sel, showProgress = FALSE)
  setnames(X, idc, "model_id")
  for (g in genes) if (!is.na(cm[g])) setnames(X, cm[g], g)
  say("    ", label, ": ", nrow(X), " rows | ID column '", idc, "'")
  for (g in genes) say(sprintf("      %-8s -> %s", g, ifelse(is.na(cm[g]), "NOT FOUND", cm[g])))
  X
}

EX <- grab_genes(f_expr, GENES, "expression")

## One row per sequencing profile with no default flag in 26Q3 (established in
## script 43). Collapse by maximum: script 43 showed first/max/mean give
## identical melanoma results, and the max is the most generous reading.
n0 <- nrow(EX)
EX <- EX[, lapply(.SD, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)),
         by = model_id, .SDcols = intersect(GENES, names(EX))]
say("    collapsed ", n0, " profile rows to ", nrow(EX), " models (max per model)")

## ---- 4. Copy number ----------------------------------------------------------
say("")
say("== 4. Copy number ==============================================")

f_cn <- file.path(CCLE_DIR, "OmicsCNGeneWGS.csv")
HAS_CN <- file.exists(f_cn)
if (!HAS_CN) {
  warn("OmicsCNGeneWGS.csv not found -- copy-number adjustment will be skipped")
} else {
  say("  ", basename(f_cn), sprintf("  %.1f MB", file.size(f_cn) / 1024^2))
  CN <- grab_genes(f_cn, GENES, "copy number")
  setnames(CN, GENES[GENES %in% names(CN)], paste0(GENES[GENES %in% names(CN)], "_cn"))
  cnc <- grep("_cn$", names(CN), value = TRUE)
  n0 <- nrow(CN)
  CN <- CN[, lapply(.SD, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)),
           by = model_id, .SDcols = cnc]
  if (n0 != nrow(CN)) say("    collapsed ", n0, " rows to ", nrow(CN), " models")
  for (cc in cnc) {
    v <- CN[[cc]]; v <- v[is.finite(v)]
    say(sprintf("    %-12s linear scale: median %.3f | range %.3f - %.3f",
                cc, median(v), min(v), max(v)))
  }
  say("    values are linear relative copy ratio; log2 is used in the models")
}

## ---- 5. Assemble -------------------------------------------------------------
say("")
say("== 5. Matched dataset ==========================================")

D <- merge(DT[!is.na(model_id)], EX, by = "model_id", all.x = FALSE)
if (HAS_CN) D <- merge(D, CN, by = "model_id", all.x = TRUE)

say("  lines with methylation AND expression: ", nrow(D))
if (HAS_CN) say("  of those, with copy number          : ", sum(is.finite(D[[paste0(GENES[1], "_cn")]])))

for (g in GENES) {
  if (!g %in% names(D)) next
  v <- D[[g]]
  say(sprintf("  %-8s expression: %.1f%% at TPM >= %d | median %.3f | SD %.3f",
              g, 100 * mean(v >= log2(TPM_CUT + 1), na.rm = TRUE), TPM_CUT,
              median(v, na.rm = TRUE), sd(v, na.rm = TRUE)))
}
say("")
say("  Expression has real variance across the panel, unlike within melanoma")
say("  alone -- that is what makes this analysis possible.")

if (nrow(D) < 200) warn("only ", nrow(D), " matched lines -- weaker than expected")

keep_cols <- c("model_id", "name", "geo_name", "lineage", "disease", "subtype", "group",
               unlist(lapply(GENES, function(g)
                 intersect(c(g, paste0(g, c("_core","_shore","_distal","_cn"))), names(D)))))
fwrite(D[, intersect(keep_cols, names(D)), with = FALSE],
       file.path(RES_DIR, "46_line_level_data.csv"))
say("  wrote 46_line_level_data.csv (the underlying data Gustavo asked for)")

## ---- 6. Pooled correlations --------------------------------------------------
say("")
say("== 6. Pooled Spearman, core vs shore vs distal =================")
say("  shore and distal are pre-specified controls: the core should be the")
say("  strongest and most negative if the association is promoter-specific.")
say("")

pooled <- list()
for (g in GENES) {
  if (!g %in% names(D)) next
  say("  ", g, ":")
  say(sprintf("    %-8s %6s %8s %11s", "set", "n", "rho", "p"))
  for (st in c("core", "shore", "distal")) {
    cn <- paste0(g, "_", st)
    if (!cn %in% names(D)) next
    ok <- is.finite(D[[cn]]) & is.finite(D[[g]])
    if (sum(ok) < 20) next
    ct <- suppressWarnings(stats::cor.test(D[[cn]][ok], D[[g]][ok], method = "spearman"))
    say(sprintf("    %-8s %6d %8.3f %11.2e", st, sum(ok), unname(ct$estimate), ct$p.value))
    pooled[[length(pooled) + 1]] <- data.table(gene = g, set = st, n = sum(ok),
                                               rho = unname(ct$estimate), p = ct$p.value)
  }

  ## Detection model: handles the floor properly where Spearman is weakened by
  ## a block of tied zeros.
  cn <- paste0(g, "_core")
  ok <- is.finite(D[[cn]]) & is.finite(D[[g]])
  det <- as.integer(D[[g]][ok] >= log2(TPM_CUT + 1))
  if (length(unique(det)) == 2) {
    fit <- stats::glm(det ~ D[[cn]][ok], family = binomial())
    s <- summary(fit)$coefficients[2, ]
    say(sprintf("    detection ~ core: OR per 0.1 beta = %.3f, p = %.2e (n = %d, %.1f%% detected)",
                exp(s[1] * 0.1), s[4], sum(ok), 100 * mean(det)))
  } else say("    detection model not estimable (no variation in detection)")
  say("")
}
POOL <- rbindlist(pooled)
fwrite(POOL, file.path(RES_DIR, "46_pooled_tests.csv"))

## ---- 7. Nested models --------------------------------------------------------
say("")
say("== 7. Nested adjustment ========================================")
say("  Reported as in scripts 36 and 40: beta, CI, and percentage of the")
say("  unadjusted effect retained. Lineage is the critical covariate -- if the")
say("  effect vanishes with lineage in the model, it is composition.")

for (g in GENES) {
  if (!g %in% names(D)) next
  cn <- paste0(g, "_core")
  cncn <- paste0(g, "_cn")

  d <- D[is.finite(get(cn)) & is.finite(get(g))]
  d[, lin := factor(ifelse(is.na(lineage), "Unknown", lineage))]
  ## lineages with a single line cannot contribute to a factor adjustment
  keep <- d[, .N, by = lin][N >= 3, lin]
  d <- d[lin %in% keep]
  d[, lin := droplevels(lin)]

  mods <- list(
    "M0 unadjusted"      = as.formula(paste(g, "~", cn)),
    "M1 + lineage"       = as.formula(paste(g, "~", cn, "+ lin"))
  )
  if (HAS_CN && cncn %in% names(d)) {
    mods[["M2 + log2 CN"]]           <- as.formula(paste(g, "~", cn, "+ log2(pmax(", cncn, ", 0.01))"))
    mods[["M3 + lineage + log2 CN"]] <- as.formula(paste(g, "~", cn, "+ lin + log2(pmax(", cncn, ", 0.01))"))
  }

  say("")
  say("  ", g, "  (n = ", nrow(d), " lines, ", nlevels(d$lin), " lineages):")
  rows <- list(); b0 <- NA_real_
  for (nm in names(mods)) {
    fit <- try(stats::lm(mods[[nm]], data = d), silent = TRUE)
    if (inherits(fit, "try-error")) { say("    ", nm, ": model failed"); next }
    co <- summary(fit)$coefficients
    if (!(cn %in% rownames(co))) { say("    ", nm, ": core term dropped"); next }
    est <- co[cn, 1]; se <- co[cn, 2]; pv <- co[cn, 4]
    ci  <- est + c(-1.96, 1.96) * se
    if (is.na(b0)) b0 <- est
    ret <- if (is.finite(b0) && abs(b0) >= 0.05)
      sprintf("%3.0f%% retained", 100 * est / b0) else "  n/a (baseline ~ 0)"
    say(sprintf("    %-24s n = %4d  beta = %+.3f [%+.3f, %+.3f]  p = %.2e   %s",
                nm, stats::nobs(fit), est, ci[1], ci[2], pv, ret))
    rows[[length(rows) + 1]] <- data.table(gene = g, model = nm, n = stats::nobs(fit),
                                           beta = est, lo = ci[1], hi = ci[2], p = pv,
                                           pct_retained = 100 * est / b0)
  }
  if (length(rows))
    fwrite(rbindlist(rows), file.path(RES_DIR, sprintf("46_models_%s.csv", tolower(g))))

  ## Copy number premises, stated rather than assumed
  if (HAS_CN && cncn %in% names(d)) {
    o1 <- is.finite(d[[cncn]]) & is.finite(d[[g]])
    o2 <- is.finite(d[[cncn]]) & is.finite(d[[cn]])
    if (sum(o1) > 20) {
      c1 <- suppressWarnings(stats::cor.test(d[[cncn]][o1], d[[g]][o1], method = "spearman"))
      say(sprintf("    premise: CN vs expression   rho = %+.3f, p = %.3g", unname(c1$estimate), c1$p.value))
    }
    if (sum(o2) > 20) {
      c2 <- suppressWarnings(stats::cor.test(d[[cncn]][o2], d[[cn]][o2], method = "spearman"))
      say(sprintf("    premise: CN vs core beta    rho = %+.3f, p = %.3g", unname(c2$estimate), c2$p.value))
    }
    say("    A covariate unrelated to both cannot be confounding the result.")
  }
}

## ---- 8. Within-lineage -------------------------------------------------------
say("")
say("== 8. Within-lineage correlations ==============================")
say("  This is the strict test. A pooled correlation can be produced entirely")
say("  by lineages differing in both variables, with no relationship inside any")
say("  of them. Only a consistent within-lineage association rules that out.")

for (g in GENES) {
  if (!g %in% names(D)) next
  cn <- paste0(g, "_core")
  d  <- D[is.finite(get(cn)) & is.finite(get(g)) & !is.na(lineage)]

  L <- d[, {
    if (.N >= MIN_LINEAGE && length(unique(get(cn))) > 3 && length(unique(get(g))) > 3) {
      ct <- suppressWarnings(stats::cor.test(get(cn), get(g), method = "spearman"))
      .(n = .N, rho = unname(ct$estimate), p = ct$p.value,
        pct_det = 100 * mean(get(g) >= log2(TPM_CUT + 1)),
        median_core = median(get(cn)))
    } else .(n = .N, rho = NA_real_, p = NA_real_,
             pct_det = 100 * mean(get(g) >= log2(TPM_CUT + 1)),
             median_core = median(get(cn)))
  }, by = lineage][order(rho)]

  LT <- L[is.finite(rho)]
  LT[, p_bh := p.adjust(p, method = "BH")]

  say("")
  say("  ", g, "  (lineages with n >= ", MIN_LINEAGE, "):")
  say(sprintf("    %-28s %5s %8s %10s %10s %8s %10s",
              "lineage", "n", "rho", "p", "BH p", "%det", "medianCore"))
  for (i in seq_len(nrow(LT))) {
    l <- LT[i]
    say(sprintf("    %-28s %5d %8.3f %10.2e %10.2e %7.1f%% %10.3f",
                l$lineage, l$n, l$rho, l$p, l$p_bh, l$pct_det, l$median_core))
  }
  if (nrow(LT)) {
    n_neg <- sum(LT$rho < 0); n_tot <- nrow(LT)
    st <- stats::binom.test(n_neg, n_tot, p = 0.5)
    say(sprintf("    median rho across lineages: %+.3f", median(LT$rho)))
    say(sprintf("    negative in %d of %d lineages (sign test p = %.3g)", n_neg, n_tot, st$p.value))
    say(sprintf("    BH-significant and negative: %d", sum(LT$p_bh < 0.05 & LT$rho < 0)))
  } else warn(g, ": no lineage had enough lines for a within-lineage correlation")

  fwrite(L, file.path(RES_DIR, sprintf("46_within_lineage_%s.csv", tolower(g))))
}

## ---- 9. Verdict --------------------------------------------------------------
say("")
say("== 9. Verdict against the pre-set criterion ====================")

for (g in GENES) {
  if (!g %in% names(D)) next
  pc <- POOL[gene == g & set == "core"]
  ps <- POOL[gene == g & set == "shore"]
  pd <- POOL[gene == g & set == "distal"]
  mf <- file.path(RES_DIR, sprintf("46_models_%s.csv", tolower(g)))
  MM <- if (file.exists(mf)) fread(mf, showProgress = FALSE) else data.table()
  wf <- file.path(RES_DIR, sprintf("46_within_lineage_%s.csv", tolower(g)))
  WW <- if (file.exists(wf)) fread(wf, showProgress = FALSE)[is.finite(rho)] else data.table()

  cond_pool  <- nrow(pc) == 1 && pc$rho < 0 && pc$p < 0.05
  cond_spec  <- nrow(pc) == 1 &&
                (!nrow(ps) || abs(pc$rho) > abs(ps$rho)) &&
                (!nrow(pd) || abs(pc$rho) > abs(pd$rho))
  m1 <- MM[grepl("lineage", model) & !grepl("CN", model)]
  cond_adj   <- nrow(m1) == 1 && is.finite(m1$pct_retained) &&
                m1$pct_retained > 50 && m1$p < 0.05
  cond_within<- nrow(WW) > 0 && median(WW$rho) < 0 && sum(WW$rho < 0) > nrow(WW) / 2

  verdict <- if (!cond_pool) "NOT SUPPORTED" else
             if (cond_pool && cond_spec && cond_adj && cond_within) "SUPPORTED" else "COMPOSITION"

  say("")
  say("  ", g, ": ", verdict)
  if (nrow(pc) == 1) say(sprintf("    pooled core rho        %+.3f (p = %.2e, n = %d)", pc$rho, pc$p, pc$n))
  if (nrow(ps) == 1) say(sprintf("    shore rho              %+.3f", ps$rho))
  if (nrow(pd) == 1) say(sprintf("    distal rho             %+.3f", pd$rho))
  if (nrow(m1) == 1) say(sprintf("    lineage-adjusted       %.0f%% retained, p = %.2e", m1$pct_retained, m1$p))
  if (nrow(WW) > 0)  say(sprintf("    within-lineage         median rho %+.3f, negative in %d of %d",
                                 median(WW$rho), sum(WW$rho < 0), nrow(WW)))
  say("    conditions met: pooled ", cond_pool, " | specificity ", cond_spec,
      " | adjustment ", cond_adj, " | within-lineage ", cond_within)
}

## ---- 10. Figures -------------------------------------------------------------
say("")
say("== 10. Figures =================================================")

for (g in GENES) {
  if (!g %in% names(D)) next
  cn <- paste0(g, "_core")
  d  <- D[is.finite(get(cn)) & is.finite(get(g))]
  d[, grp := factor(ifelse(group == "Melanoma", "Melanoma", "Other cancer lineages"),
                    levels = c("Other cancer lineages", "Melanoma"))]
  pc <- POOL[gene == g & set == "core"]

  p <- ggplot(d, aes(x = get(cn), y = get(g))) +
    geom_point(aes(colour = grp, size = grp, alpha = grp)) +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                colour = "grey20", linewidth = 0.45, fill = "grey80") +
    geom_hline(yintercept = log2(TPM_CUT + 1), linetype = "dashed",
               linewidth = 0.3, colour = "grey45") +
    annotate("text", x = 0.02, y = log2(TPM_CUT + 1), label = "TPM = 1",
             hjust = 0, vjust = -0.6, size = 2.7, colour = "grey35") +
    scale_colour_manual(values = c("Other cancer lineages" = "grey55",
                                   "Melanoma" = "firebrick"), name = NULL) +
    scale_size_manual(values = c("Other cancer lineages" = 0.9, "Melanoma" = 1.7), guide = "none") +
    scale_alpha_manual(values = c("Other cancer lineages" = 0.45, "Melanoma" = 0.9), guide = "none") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title    = paste0(g, ": promoter-core methylation vs expression across cancer cell lines"),
      subtitle = if (nrow(pc) == 1)
        sprintf("GSE68379 HM450 + DepMap expression | n = %d lines | Spearman rho = %.3f, p = %.2e",
                pc$n, pc$rho, pc$p) else
        sprintf("GSE68379 HM450 + DepMap expression | n = %d lines", nrow(d)),
      x = "promoter-core beta (mean of TSS-proximal probes)",
      y = expression(log[2]*"(TPM + 1)")
    ) +
    theme_bw(base_size = 9) +
    theme(legend.position = "top",
          legend.margin = margin(0, 0, -4, 0),
          panel.grid.minor = element_blank(),
          plot.title    = element_text(face = "bold", size = 10.5),
          plot.subtitle = element_text(size = 7.6, colour = "grey25"),
          plot.margin   = margin(10, 14, 8, 10))

  out <- file.path(FIG_DIR, sprintf("46_%s_scatter.png", tolower(g)))
  ggsave(out, p, width = 6.8, height = 4.6, dpi = 400)
  say("  wrote: ", basename(out))
}

## ---- 11. Limits --------------------------------------------------------------
say("")
say("== 11. Limits to carry into the write-up =======================")
say("  1. This is a PAN-CANCER analysis. It does not establish a")
say("     methylation-expression relationship WITHIN melanoma cells, which")
say("     remains untestable because expression there has no variance.")
say("  2. A SUPPORTED verdict means methylation tracks expression across")
say("     cancer cells, with melanoma at the methylated-and-silent end. It")
say("     does not show methylation CAUSED the silencing. Demethylation in")
say("     vitro is the experiment that would.")
say("  3. A COMPOSITION verdict means the pooled correlation reflects lineages")
say("     differing in both variables. That is still informative -- it is how")
say("     a lineage-restricted epigenetic state looks -- but it is a weaker")
say("     claim and must be labelled as such.")
say("  4. Methylation (Sanger) and expression and copy number (Broad) come")
say("     from different laboratories and different cell stocks.")
say("  5. Cell lines are not tumours.")
