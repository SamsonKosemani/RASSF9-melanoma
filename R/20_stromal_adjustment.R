# ==============================================================================
# 20_stromal_adjustment.R
#
# Does the nevus-versus-melanoma difference in RASSF9 survive adjustment for
# tissue composition?
#
# WHY THIS EXISTS. Script 14 found RASSF9 significantly lower in primary
# melanoma than in nevi, in two independent cohorts, with large effects
# (GSE112509 d = -0.886, p = 8.8e-04; GSE98394 d = -1.000, p = 1.8e-05). Script
# 19 then showed, at single-cell resolution, that RASSF9 in nevus tissue is
# concentrated in FIBROBLASTS and ENDOTHELIUM (11.0 and 24.9 TPM) rather than in
# melanocytes (0.57 TPM). Nevi are stroma-rich; melanomas are tumour-cell dense.
# So the bulk difference may measure composition rather than anything happening
# inside the melanocytic lineage -- which would also explain the negative
# RASSF9-purity correlation in TCGA (script 06, rho = -0.123, p = 0.027).
#
# THE LOGIC. Confounding by composition requires two things to be true, and both
# are tested separately rather than assumed:
#   (1) stromal content must DIFFER between nevi and melanomas, and
#   (2) RASSF9 must TRACK stromal content within a group.
# If either fails, composition cannot be the explanation. If both hold, the
# question becomes whether the group effect survives adjustment.
#
# PRE-SPECIFIED before any model was fitted:
#   - Primary: linear model expr ~ group + fibroblast + endothelial + immune
#     score, per cohort per gene. The group coefficient before and after is the
#     result.
#   - Secondary: add the melanocytic score, which is the closest available
#     proxy for how much melanocytic material each sample contains.
#   - Sensitivity: a stroma-matched subset, trimmed to the stromal-score range
#     the two groups share, which needs no model at all.
#   - Both genes, both cohorts, reported whatever happens.
#
# WHAT ADJUSTMENT CANNOT DO. A marker-score covariate is a proxy, not a
# deconvolution. It cannot recover melanocyte-intrinsic expression from bulk
# tissue, and over-adjustment is possible because stromal content and lesion
# type are deeply entangled -- in a sense being stroma-poor is part of what
# being a melanoma IS. A surviving effect is therefore not proof of melanocytic
# biology, and a vanishing effect is the more interpretable outcome. This
# asymmetry is stated up front so the result is not over-read in either
# direction.
#
# Output: results/geo/stromal_scores.csv
#         results/geo/stromal_adjustment_models.csv
#         figures/geo/geo_<gene>_stromal_adjustment.{pdf,png,svg}
#         results/logs/20_stromal.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(patchwork)
  library(EnsDb.Hsapiens.v86)
  library(ensembldb)
})
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))

LOG <- file.path(PATHS$logs, "20_stromal.log")
cat("Stromal adjustment —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

GEO_DIR <- file.path(PATHS$raw, "geo")
RES_GEO <- file.path(PATHS$results, "geo")
FIG_GEO <- file.path(PATHS$figures, "geo")
invisible(lapply(c(RES_GEO, FIG_GEO), dir.create, recursive = TRUE,
                 showWarnings = FALSE))

ENS <- c(RASSF9 = "ENSG00000198774", RASSF10 = "ENSG00000189431")
options(timeout = 1800)

# Same panels as script 19, so the bulk scores and the single-cell lineage
# assignment are built from identical marker definitions.
MARKERS <- list(
  Melanocytic  = c("MITF", "PMEL", "MLANA", "DCT", "TYR", "TYRP1", "SOX10"),
  Fibroblast   = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB"),
  Keratinocyte = c("KRT14", "KRT5", "KRT1", "KRT10", "KRT6A"),
  Immune       = c("PTPRC", "CD3E", "CD68", "LYZ", "CD14", "AIF1"),
  Endothelial  = c("PECAM1", "VWF", "CDH5", "CLDN5")
)
ADJ_PRIMARY <- c("Fibroblast", "Endothelial", "Immune")

norm_key <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

cohens_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  s <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  (mean(x) - mean(y)) / s
}

# ==============================================================================
# 1. Symbol -> Ensembl, offline
# ==============================================================================
hdr("1. Marker gene annotation")

sym <- unique(unlist(MARKERS))
ann <- ensembldb::genes(EnsDb.Hsapiens.v86,
                        filter = AnnotationFilterList(
                          GeneNameFilter(sym)),
                        columns = c("gene_id", "gene_name"),
                        return.type = "data.frame")
sym2ens <- ann |> dplyr::filter(gene_name %in% sym) |>
  dplyr::select(gene_name, gene_id) |> dplyr::distinct()

log_line("  markers requested: ", length(sym), logfile = LOG)
log_line("  mapped to Ensembl: ", dplyr::n_distinct(sym2ens$gene_name),
         logfile = LOG)
miss <- setdiff(sym, sym2ens$gene_name)
if (length(miss)) log_line("  [FLAG] unmapped: ", paste(miss, collapse = ", "),
                           logfile = LOG)

ens_of <- function(panel) {
  sym2ens$gene_id[sym2ens$gene_name %in% MARKERS[[panel]]]
}

# ==============================================================================
# 2. Load both cohorts (loaders identical to script 14)
# ==============================================================================
hdr("2. Cohorts")

read_matrix <- function(path) {
  x <- as.data.frame(readr::read_tsv(path, show_col_types = FALSE))
  rownames(x) <- x[[1]]
  as.matrix(x[, -1, drop = FALSE])
}

# ---- GSE112509 --------------------------------------------------------------
f1 <- list.files(file.path(GEO_DIR, "GSE112509"), pattern = "normalized_counts",
                 full.names = TRUE)
if (!length(f1)) stop("GSE112509 counts not found; run script 13 first.")
m1 <- as.data.frame(readr::read_tsv(f1[1], show_col_types = FALSE))
rownames(m1) <- m1[[1]]; m1 <- as.matrix(m1[, -1, drop = FALSE])
L1 <- log2(m1 + 1)
log_line("  GSE112509: ", nrow(L1), " genes x ", ncol(L1), " samples",
         logfile = LOG)

gse1 <- getGEO("GSE112509", destdir = GEO_DIR, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
pd1 <- pData(gse1)
tis1 <- grep("tissue", names(pd1), value = TRUE, ignore.case = TRUE)
tis1 <- tis1[grepl("nevi|melanoma",
                   tolower(sapply(pd1[tis1], function(z) paste(z, collapse = ""))))][1]
map1 <- tibble::tibble(
  key = norm_key(pd1$title),
  group = ifelse(grepl("nevi|nevus|naev", tolower(as.character(pd1[[tis1]]))),
                 "Nevus", "Melanoma"))
g1 <- tibble::tibble(sample = colnames(L1), key = norm_key(colnames(L1))) |>
  dplyr::left_join(map1, by = "key") |>
  dplyr::mutate(group = dplyr::coalesce(
    group,
    dplyr::case_when(grepl("_N$", sample) ~ "Nevus",
                     grepl("_M$", sample) ~ "Melanoma",
                     TRUE ~ NA_character_))) |>
  dplyr::filter(!is.na(group)) |>
  dplyr::mutate(dataset = "GSE112509")

# ---- GSE98394 ---------------------------------------------------------------
# The expression file covers ~12,800 genes and omits some markers, so the
# counts file (full gene set, converted to log2 CPM here) is used for the
# scores. Using ONE matrix for every marker keeps the scores internally
# comparable, which matters more here than matching script 14's per-gene source
# choice -- and the gene of interest is still taken the same way script 14 took
# it, so the headline values do not move.
f2e <- list.files(file.path(GEO_DIR, "GSE98394"), pattern = "expression",
                  full.names = TRUE)
f2c <- list.files(file.path(GEO_DIR, "GSE98394"), pattern = "counts",
                  full.names = TRUE)
m2e <- if (length(f2e)) read_matrix(f2e[1]) else NULL
m2c <- if (length(f2c)) read_matrix(f2c[1]) else NULL
if (is.null(m2c)) stop("GSE98394 counts file not found; run script 13 first.")

lib2 <- colSums(m2c, na.rm = TRUE)
L2 <- log2(t(t(m2c) / lib2) * 1e6 + 1)
log_line("  GSE98394:  ", nrow(L2), " genes x ", ncol(L2),
         " samples (log2 CPM from the counts file)", logfile = LOG)

gse2 <- getGEO("GSE98394", destdir = GEO_DIR, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
pd2 <- pData(gse2)
tis2 <- grep("tissue", names(pd2), value = TRUE, ignore.case = TRUE)
tis2 <- tis2[grepl("nevus|nevi|melanoma",
                   tolower(sapply(pd2[tis2], function(z) paste(z, collapse = ""))))][1]
hdrs <- colnames(L2)
score_col <- sapply(names(pd2), function(cc) {
  v <- as.character(pd2[[cc]])
  sum(vapply(hdrs, function(h) any(grepl(h, v, fixed = TRUE)), logical(1)))
})
best <- names(pd2)[which.max(score_col)]
log_line("  GSE98394 id column: ", best, " (matches ", max(score_col), " of ",
         length(hdrs), ")", logfile = LOG)
grp2 <- ifelse(grepl("nevus|nevi|naev", tolower(as.character(pd2[[tis2]]))),
               "Nevus", "Melanoma")
vals2 <- as.character(pd2[[best]])
g2 <- tibble::tibble(sample = hdrs) |>
  dplyr::mutate(group = vapply(sample, function(h) {
    i <- which(grepl(h, vals2, fixed = TRUE))
    if (length(i) == 1) grp2[i] else NA_character_
  }, character(1))) |>
  dplyr::filter(!is.na(group)) |>
  dplyr::mutate(dataset = "GSE98394")

COHORTS <- list(GSE112509 = list(mat = L1, meta = g1),
                GSE98394  = list(mat = L2, meta = g2))

# ==============================================================================
# 3. Composition scores
# ==============================================================================
# Each score is the mean z-scored log expression of its panel. z-scoring is done
# WITHIN cohort, across all samples of that cohort, so a score is a relative
# statement about samples in the same dataset and nothing is compared across
# platforms.
hdr("3. Composition scores")

panel_score <- function(mat, ids) {
  ids_nov <- sub("\\..*$", "", rownames(mat))
  idx <- which(ids_nov %in% ids)
  if (!length(idx)) return(NULL)
  sub <- mat[idx, , drop = FALSE]
  z <- t(scale(t(sub)))                    # per gene across samples
  z <- z[rowSums(is.na(z)) == 0, , drop = FALSE]
  if (!nrow(z)) return(NULL)
  list(score = colMeans(z), n_genes = nrow(z))
}

scores <- list()
for (cn in names(COHORTS)) {
  mat <- COHORTS[[cn]]$mat
  meta <- COHORTS[[cn]]$meta
  sc <- meta
  log_line("\n  ", cn, logfile = LOG)
  for (p in names(MARKERS)) {
    ps <- panel_score(mat, ens_of(p))
    if (is.null(ps)) {
      log_line("    [FLAG] ", p, ": no panel genes found", logfile = LOG)
      sc[[p]] <- NA_real_; next
    }
    sc[[p]] <- as.numeric(ps$score[sc$sample])
    log_line(sprintf("    %-13s %d of %d panel genes", p, ps$n_genes,
                     length(MARKERS[[p]])), logfile = LOG)
  }
  for (g in names(ENS)) {
    idx <- which(sub("\\..*$", "", rownames(mat)) == ENS[[g]])
    sc[[g]] <- if (length(idx) == 0) NA_real_
               else if (length(idx) > 1) colSums(mat[idx, sc$sample, drop = FALSE])
               else as.numeric(mat[idx, sc$sample])
  }
  scores[[cn]] <- sc
}
S <- dplyr::bind_rows(scores)
readr::write_csv(S, file.path(RES_GEO, "stromal_scores.csv"))
log_line("\n  wrote: results/geo/stromal_scores.csv", logfile = LOG)

# ==============================================================================
# 4. Condition (1): does stromal content differ between the groups?
# ==============================================================================
hdr("4. Condition 1 -- composition differs by group?")

cond1 <- list()
for (cn in names(scores)) {
  s <- scores[[cn]]
  log_line("\n  ", cn, logfile = LOG)
  log_line(sprintf("    %-13s %9s %9s %8s %10s", "panel", "nevus", "melanoma",
                   "d", "p"), logfile = LOG)
  for (p in names(MARKERS)) {
    if (all(is.na(s[[p]]))) next
    a <- s[[p]][s$group == "Nevus"]; b <- s[[p]][s$group == "Melanoma"]
    wt <- suppressWarnings(stats::wilcox.test(a, b))
    d  <- cohens_d(b, a)                    # melanoma minus nevus
    cond1[[length(cond1) + 1L]] <- tibble::tibble(
      dataset = cn, panel = p, mean_nevus = mean(a), mean_melanoma = mean(b),
      d_melanoma_vs_nevus = d, p = wt$p.value)
    log_line(sprintf("    %-13s %9.3f %9.3f %8.2f %10s", p, mean(a), mean(b),
                     d, fmt_p(wt$p.value)), logfile = LOG)
  }
}
cond1 <- dplyr::bind_rows(cond1)

strom_diff <- cond1 |> dplyr::filter(panel %in% ADJ_PRIMARY, p < 0.05)
log_line("", logfile = LOG)
if (nrow(strom_diff)) {
  log_line("  CONDITION 1 HOLDS: stromal content differs by group in ",
           nrow(strom_diff), " of ", length(ADJ_PRIMARY) * length(scores),
           " panel-cohort tests.", logfile = LOG)
  log_line("  Composition is therefore a live confounder and adjustment is",
           logfile = LOG)
  log_line("  warranted.", logfile = LOG)
} else {
  log_line("  CONDITION 1 FAILS: stromal content does NOT differ between nevi",
           logfile = LOG)
  log_line("  and melanomas in these cohorts. Composition cannot explain the",
           logfile = LOG)
  log_line("  RASSF9 difference, and the single-cell attribution does not",
           logfile = LOG)
  log_line("  transfer to this comparison.", logfile = LOG)
}

# ==============================================================================
# 5. Condition (2): does the gene track stromal content WITHIN group?
# ==============================================================================
# Within group, so the correlation cannot be produced by the group difference
# itself. This is the bulk-level echo of what script 19 saw per cell.
hdr("5. Condition 2 -- gene tracks stroma within group?")

cond2 <- list()
for (cn in names(scores)) {
  s <- scores[[cn]]
  for (g in names(ENS)) {
    if (all(is.na(s[[g]]))) next
    for (grp in c("Nevus", "Melanoma")) {
      ss <- s[s$group == grp, ]
      fib <- ss$Fibroblast
      ct <- suppressWarnings(stats::cor.test(ss[[g]], fib, method = "spearman"))
      cond2[[length(cond2) + 1L]] <- tibble::tibble(
        dataset = cn, gene = g, group = grp, n = nrow(ss),
        rho = unname(ct$estimate), p = ct$p.value)
    }
  }
}
cond2 <- dplyr::bind_rows(cond2)
log_line(sprintf("  %-11s %-9s %-9s %5s %8s %10s", "dataset", "gene", "group",
                 "n", "rho", "p"), logfile = LOG)
for (i in seq_len(nrow(cond2))) {
  log_line(sprintf("  %-11s %-9s %-9s %5d %8.3f %10s", cond2$dataset[i],
                   cond2$gene[i], cond2$group[i], cond2$n[i], cond2$rho[i],
                   fmt_p(cond2$p[i])), logfile = LOG)
}
log_line("", logfile = LOG)
log_line("  Spearman against the fibroblast score, computed WITHIN group so the",
         logfile = LOG)
log_line("  group difference itself cannot generate it.", logfile = LOG)

# ==============================================================================
# 6. The adjustment
# ==============================================================================
hdr("6. Group effect before and after adjustment")

fit_one <- function(s, g, covars, label) {
  dd <- s |> dplyr::select(dplyr::all_of(c(g, "group", covars))) |> stats::na.omit()
  if (!nrow(dd)) return(NULL)
  dd$group <- factor(dd$group, levels = c("Nevus", "Melanoma"))
  f <- stats::as.formula(paste0("`", g, "` ~ group",
                                if (length(covars))
                                  paste0(" + ", paste(sprintf("`%s`", covars),
                                                      collapse = " + ")) else ""))
  m <- stats::lm(f, data = dd)
  cf <- summary(m)$coefficients
  row <- grep("^groupMelanoma$", rownames(cf))
  ci <- stats::confint(m)[row, ]
  tibble::tibble(
    model = label, n = nrow(dd),
    beta = cf[row, 1], se = cf[row, 2], p = cf[row, 4],
    ci_lo = ci[1], ci_hi = ci[2],
    adj_r2 = summary(m)$adj.r.squared)
}

models <- list()
for (cn in names(scores)) {
  s <- scores[[cn]]
  for (g in names(ENS)) {
    if (all(is.na(s[[g]]))) next
    fits <- dplyr::bind_rows(
      fit_one(s, g, character(0), "unadjusted"),
      fit_one(s, g, ADJ_PRIMARY, "adjusted (stroma)"),
      fit_one(s, g, c(ADJ_PRIMARY, "Melanocytic"),
              "adjusted (stroma + melanocytic)"))
    if (is.null(fits) || !nrow(fits)) next
    fits <- fits |> dplyr::mutate(dataset = cn, gene = g, .before = 1)
    models[[length(models) + 1L]] <- fits

    log_line("\n  ", cn, " / ", g, logfile = LOG)
    log_line(sprintf("    %-32s %5s %9s %18s %10s %8s", "model", "n", "beta",
                     "95% CI", "p", "adjR2"), logfile = LOG)
    for (i in seq_len(nrow(fits))) {
      log_line(sprintf("    %-32s %5d %9.3f  [%7.3f, %7.3f] %10s %8.3f",
                       fits$model[i], fits$n[i], fits$beta[i], fits$ci_lo[i],
                       fits$ci_hi[i], fmt_p(fits$p[i]), fits$adj_r2[i]),
               logfile = LOG)
    }
    u <- fits |> dplyr::filter(model == "unadjusted")
    a <- fits |> dplyr::filter(model == "adjusted (stroma)")
    if (nrow(u) && nrow(a)) {
      shrink <- 100 * (1 - a$beta / u$beta)
      log_line(sprintf("    the group coefficient shrinks by %.0f%% on adjustment",
                       shrink), logfile = LOG)
      if (a$p >= 0.05 && u$p < 0.05) {
        log_line("    -> the group effect does NOT survive. Composition is a",
                 logfile = LOG)
        log_line("       sufficient explanation for this cohort.", logfile = LOG)
      } else if (a$p < 0.05 && sign(a$beta) == sign(u$beta)) {
        log_line("    -> the group effect SURVIVES adjustment, same direction.",
                 logfile = LOG)
      } else {
        log_line("    -> equivocal; read the interval, not the p value.",
                 logfile = LOG)
      }
    }
  }
}
models <- dplyr::bind_rows(models)
readr::write_csv(models, file.path(RES_GEO, "stromal_adjustment_models.csv"))
readr::write_csv(cond1, file.path(RES_GEO, "stromal_condition1_composition.csv"))
readr::write_csv(cond2, file.path(RES_GEO, "stromal_condition2_within_group.csv"))
log_line("\n  wrote: results/geo/stromal_adjustment_models.csv", logfile = LOG)

# ==============================================================================
# 7. Model-free sensitivity: caliper-matched subset
# ==============================================================================
# An earlier version of this section trimmed to the min-max overlap of the
# fibroblast score. That is a near-useless criterion: in GSE112509 it retained
# 78 of 80 samples, so the "matched" test was the full cohort under another
# name, and it duly contradicted the adjusted model. Replaced with 1:1 greedy
# nearest-neighbour matching, without replacement, on a composite stromal score,
# inside a caliper.
#
# Crucially this section now reports whether matching ACHIEVED balance, as the
# standardised mean difference before and after. A matched analysis whose SMD is
# still large has not controlled for anything, and saying so is the difference
# between a sensitivity analysis and a decoration.
hdr("7. Caliper-matched subset (model-free)")

CALIPER <- 0.20        # of the SD of the linear propensity, the usual convention
SMD_OK  <- 0.10        # balance achieved below this, per panel

smd <- function(a, b) {
  s <- sqrt((var(a) + var(b)) / 2)
  if (!is.finite(s) || s == 0) return(NA_real_)
  (mean(b) - mean(a)) / s
}

# Matching on the MEAN of the three stromal panels was the previous approach and
# it was unsound: in GSE112509 the immune score is higher in melanoma (d = +0.60)
# while the endothelial score is higher in nevus (d = -0.52), so the two cancel
# and a composite can read as balanced while both of its components are
# imbalanced. Balance has to be achieved and demonstrated PANEL BY PANEL.
#
# So: a propensity score from all three panels jointly (logistic regression of
# group on the panels), greedy 1:1 nearest-neighbour matching without
# replacement on the linear propensity inside a 0.20-SD caliper, and then the
# standardised mean difference reported for EVERY panel separately. This is the
# standard construction and, more to the point, it cannot report success while
# any single panel is still out of balance.
match_propensity <- function(s, panels) {
  nv_idx <- which(s$group == "Nevus")
  ml_idx <- which(s$group == "Melanoma")
  if (length(nv_idx) < 3 || length(ml_idx) < 3) return(NULL)

  dd <- s[, panels, drop = FALSE]
  dd$.trt <- as.integer(s$group == "Melanoma")
  f <- stats::as.formula(paste(".trt ~",
                               paste(sprintf("`%s`", panels), collapse = " + ")))
  fit <- suppressWarnings(stats::glm(f, data = dd, family = binomial()))
  ps <- stats::predict(fit, type = "link")        # linear propensity (logit)
  if (any(!is.finite(ps))) return(NULL)
  thr <- CALIPER * stats::sd(ps)

  avail <- ml_idx
  pairs <- list()
  # Nevi taken in order of propensity: deterministic and reproducible.
  for (i in nv_idx[order(ps[nv_idx])]) {
    if (!length(avail)) break
    dist <- abs(ps[avail] - ps[i])
    j <- which.min(dist)
    if (dist[j] > thr) next           # no partner close enough; nevus dropped
    pairs[[length(pairs) + 1L]] <- c(i, avail[j])
    avail <- avail[-j]
  }
  if (!length(pairs)) return(NULL)
  list(idx = unlist(pairs), caliper = thr, n_pairs = length(pairs),
       ps_sd = stats::sd(ps))
}

# Per-panel balance table, which is what decides whether a matched analysis
# controls for anything.
balance_table <- function(s, panels) {
  vapply(panels, function(p) {
    smd(s[[p]][s$group == "Nevus"], s[[p]][s$group == "Melanoma"])
  }, numeric(1))
}

matched <- list()
balance <- list()
for (cn in names(scores)) {
  s <- scores[[cn]] |>
    dplyr::filter(!is.na(Fibroblast), !is.na(Immune), !is.na(Endothelial))
  if (!nrow(s)) next

  pre <- balance_table(s, ADJ_PRIMARY)
  log_line("\n  ", cn, logfile = LOG)
  log_line("    per-panel SMD before matching:", logfile = LOG)
  for (p in ADJ_PRIMARY) {
    log_line(sprintf("      %-13s %+.3f", p, pre[[p]]), logfile = LOG)
  }

  mm <- match_propensity(s, ADJ_PRIMARY)
  if (is.null(mm)) {
    log_line("    no pairs within the caliper -- the groups do not overlap on",
             logfile = LOG)
    log_line("    stromal content, so NO model-free test is possible here and",
             logfile = LOG)
    log_line("    the adjusted model in section 6 is extrapolating across",
             logfile = LOG)
    log_line("    ranges the data do not share. Report it that way.",
             logfile = LOG)
    balance[[cn]] <- tibble::tibble(dataset = cn, panel = ADJ_PRIMARY,
                                    smd_pre = unname(pre[ADJ_PRIMARY]),
                                    smd_post = NA_real_, n_pairs = 0L,
                                    balanced = FALSE)
    next
  }
  sub <- s[mm$idx, ]
  post <- balance_table(sub, ADJ_PRIMARY)
  worst <- max(abs(post), na.rm = TRUE)
  ok <- worst < SMD_OK

  log_line(sprintf("    caliper %.3f on the linear propensity -> %d pairs (%d of %d samples)",
                   mm$caliper, mm$n_pairs, nrow(sub), nrow(s)), logfile = LOG)
  log_line("    per-panel SMD after matching:", logfile = LOG)
  for (p in ADJ_PRIMARY) {
    log_line(sprintf("      %-13s %+.3f %s", p, post[[p]],
                     if (abs(post[[p]]) < SMD_OK) "" else "  <- imbalanced"),
             logfile = LOG)
  }
  log_line(sprintf("    worst panel |SMD| = %.3f  [%s]", worst,
                   if (ok) "BALANCED" else "STILL IMBALANCED"), logfile = LOG)

  balance[[cn]] <- tibble::tibble(dataset = cn, panel = ADJ_PRIMARY,
                                  smd_pre = unname(pre[ADJ_PRIMARY]),
                                  smd_post = unname(post[ADJ_PRIMARY]),
                                  n_pairs = mm$n_pairs, balanced = ok)

  if (!ok) {
    log_line("    [FLAG] at least one panel is still imbalanced, so the test",
             logfile = LOG)
    log_line("           below controls for less than it appears to. Treat it",
             logfile = LOG)
    log_line("           as suggestive and rely on section 6.", logfile = LOG)
  }
  if (mm$n_pairs < 15) {
    log_line("    [FLAG] only ", mm$n_pairs, " pairs; underpowered, and the",
             logfile = LOG)
    log_line("           retained samples are a selected tail of the larger",
             logfile = LOG)
    log_line("           group rather than a representative subset.",
             logfile = LOG)
  }

  for (g in names(ENS)) {
    if (all(is.na(sub[[g]]))) {
      log_line("    ", g, ": not available in this cohort", logfile = LOG)
      next
    }
    a <- sub[[g]][sub$group == "Nevus"]; b <- sub[[g]][sub$group == "Melanoma"]
    wt <- suppressWarnings(stats::wilcox.test(a, b))
    d <- cohens_d(b, a)
    matched[[length(matched) + 1L]] <- tibble::tibble(
      dataset = cn, gene = g, n_nevus = length(a), n_melanoma = length(b),
      worst_smd_post = worst, balanced = ok, d = d, p = wt$p.value)
    log_line(sprintf("    %-9s d = %6.3f, p = %s", g, d, fmt_p(wt$p.value)),
             logfile = LOG)
  }
}
matched <- dplyr::bind_rows(matched)
balance <- dplyr::bind_rows(balance)
if (nrow(matched)) {
  readr::write_csv(matched, file.path(RES_GEO, "stromal_matched_subset.csv"))
  readr::write_csv(balance, file.path(RES_GEO, "stromal_matching_balance.csv"))
  log_line("\n  wrote: results/geo/stromal_matched_subset.csv and _balance.csv",
           logfile = LOG)
}

# ==============================================================================
# 8. Figures
# ==============================================================================
hdr("8. Figures")

make_fig <- function(g) {
  sg <- S |> dplyr::filter(!is.na(.data[[g]]), !is.na(Fibroblast)) |>
    dplyr::mutate(group = factor(group, levels = c("Nevus", "Melanoma")))
  if (!nrow(sg)) return(NULL)

  pal <- c(Nevus = PAL$nevus, Melanoma = PAL$melanoma)

  # A -- condition 1: does stromal content differ by group?
  pa <- ggplot(sg, aes(x = group, y = Fibroblast, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, colour = PAL$ink,
                 linewidth = 0.3, alpha = 0.85) +
    geom_jitter(width = 0.14, size = 0.6, colour = PAL$ink, alpha = 0.6) +
    facet_wrap(~ dataset, nrow = 1) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(x = NULL, y = "Fibroblast score (z)") +
    theme_journal(grid = TRUE)

  # B -- condition 2: does the gene track stroma within group?
  pb <- ggplot(sg, aes(x = Fibroblast, y = .data[[g]], colour = group)) +
    geom_point(size = 0.8, alpha = 0.8) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, formula = y ~ x) +
    facet_wrap(~ dataset, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = pal, name = NULL) +
    labs(x = "Fibroblast score (z)",
         y = bquote(bold(.(g)~"expression (log"[2]*")"))) +
    theme_journal(grid = TRUE) +
    theme(legend.position = "top", legend.justification = "left")

  # C -- the adjustment, as a forest of group coefficients
  mg <- models |> dplyr::filter(gene == g) |>
    dplyr::mutate(model = factor(model, levels = rev(c(
      "unadjusted", "adjusted (stroma)",
      "adjusted (stroma + melanocytic)"))))
  pc <- ggplot(mg, aes(y = model, x = beta, colour = dataset)) +
    geom_vline(xintercept = 0, linetype = "22", colour = PAL$grey_pt,
               linewidth = 0.4) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0,
                  linewidth = 0.5, orientation = "y",
                  position = position_dodge(width = 0.5)) +
    geom_point(size = 1.6, position = position_dodge(width = 0.5)) +
    scale_colour_manual(values = c(GSE112509 = PAL$uv_low,
                                   GSE98394 = "#7E3FF2"), name = NULL) +
    labs(y = NULL, x = "Melanoma vs nevus coefficient (95% CI)") +
    theme_journal(grid = FALSE) +
    theme(panel.grid.major.x = element_line(linewidth = 0.25,
                                            colour = "grey92"),
          legend.position = "top", legend.justification = "left",
          axis.text.y = element_text(size = FIG$base_size - 1, hjust = 1),
          plot.margin = margin(4, 10, 4, 4, unit = "pt"))

  (pa / pb / pc) + patchwork::plot_layout(heights = c(1, 1.2, 1))
}

for (g in names(ENS)) {
  p <- make_fig(g)
  if (is.null(p)) { log_line("  ", g, ": skipped", logfile = LOG); next }
  fig <- p + plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = FIG$base_size + 3, face = "bold",
                                  colour = PAL$ink))

  stem <- paste0("geo_", tolower(g), "_stromal_adjustment")
  save_fig_to(fig, stem, FIG_GEO, width = FIG$width_double, height = 190)

  mg <- models |> dplyr::filter(gene == g)
  txt <- paste(sprintf("%s: unadjusted beta %.3f (p %s), adjusted for stroma %.3f (p %s)",
                       mg$dataset[mg$model == "unadjusted"],
                       mg$beta[mg$model == "unadjusted"],
                       fmt_p(mg$p[mg$model == "unadjusted"]),
                       mg$beta[mg$model == "adjusted (stroma)"],
                       fmt_p(mg$p[mg$model == "adjusted (stroma)"])),
               collapse = "; ")

  add_legend(
    stem,
    paste0("Figure. Tissue composition and the nevus-versus-melanoma ",
           "difference in ", g, "."),
    paste0(
      "Composition scores are the mean z-scored log expression of marker ",
      "panels, computed within each cohort: fibroblast (COL1A1, COL1A2, ",
      "COL3A1, DCN, LUM, PDGFRB), endothelial (PECAM1, VWF, CDH5, CLDN5), ",
      "immune (PTPRC, CD3E, CD68, LYZ, CD14, AIF1) and melanocytic (MITF, ",
      "PMEL, MLANA, DCT, TYR, TYRP1, SOX10) — the same panels used for the ",
      "single-cell lineage assignment. (A) Fibroblast score by lesion type in ",
      "each cohort, testing whether stromal content differs between the groups ",
      "at all. (B) ", g, " expression against the fibroblast score, with lines ",
      "fitted within group, so that any association cannot be produced by the ",
      "group difference itself. (C) The melanoma-versus-nevus coefficient with ",
      "95% confidence interval, unadjusted and after adjustment for stromal and ",
      "then melanocytic content: ", txt, ". A marker-score covariate is a ",
      "proxy for composition rather than a deconvolution, and stromal content ",
      "is itself partly constitutive of lesion type, so a coefficient that ",
      "survives adjustment is not by itself evidence of melanocyte-intrinsic ",
      "regulation."))

  log_line("  wrote: ", stem, ".{pdf,png,svg}", logfile = LOG)
}

report_warnings(logfile = LOG)
write_session_info("20_stromal")
log_line("\nFull log: ", LOG, logfile = NULL)
