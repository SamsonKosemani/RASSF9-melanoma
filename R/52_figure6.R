## =============================================================================
## 52_figure6.R
##
## Builds Figure 6 to the structure agreed with GPA-M on 3 November 2026:
##   A  RASSF9 expression, nevi vs primary melanoma, GSE112509 and GSE98394
##   B  Cohen's d with 95% CI for the same two cohorts
##   C  RASSF9 detection by cell type, GSE72056, depth-matched
##   D  Promoter methylation vs expression, TCGA-SKCM, n = 469
##   E  Promoter-core methylation, 56 melanoma lines vs 940 other lineages
##   F  Promoter-core methylation vs expression, 741 cancer cell lines
##
## The GSE298774 spatial panels move to supplementary and are not built here.
##
## PANEL D NEEDS DATA THAT WAS NEVER SAVED
##   Script 35 wrote summary tables but no per-patient table, so the scatter
##   has to be rebuilt from skcm_hm450.rds and skcm_star_counts.rds. That is
##   done once and cached to results/tcga/fig6_tcga_perpatient.csv; later runs
##   read the cache.
##
##   Reconstructing it means re-deriving choices script 35 made, in particular
##   how samples were collapsed to patients (473 expression samples map to 469
##   patients). Rather than assume, the script computes both plausible rules
##   and keeps whichever reproduces the published rho = -0.379. If neither
##   does, it stops. A panel that silently disagreed with the Results text
##   would be worse than no panel.
##
## PRODUCTION
##   Drawn at 180 mm, the final printed width, so the point sizes written here
##   are the ones that print. Place at 100%; do not scale on the slide.
##
## OUTPUTS
##   figures/tcga/Figure6.pdf    vector, for submission
##   figures/tcga/Figure6.png    600 dpi
##   results/tcga/fig6_tcga_perpatient.csv   (cache, panel D)
##   results/logs/52_figure6.log
## =============================================================================

if (file.exists("R/00_setup.R")) source("R/00_setup.R")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

ROOT    <- if (exists("PROJECT_ROOT")) PROJECT_ROOT else getwd()
FIG_DIR <- file.path(ROOT, "figures", "tcga")
LOG_DIR <- file.path(ROOT, "results", "logs")
for (d in c(FIG_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_PATH <- file.path(LOG_DIR, "52_figure6.log")
.logcon  <- file(LOG_PATH, open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = .logcon) }
WARN <- character(0)
warn <- function(...) { m <- paste0(...); WARN <<- c(WARN, m); say("  [WARN] ", m) }
on.exit({
  if (length(WARN)) { say(""); say("WARNINGS (", length(WARN), "):"); for (w in WARN) say("  - ", w) }
  else say("  no warnings raised")
  say(""); say("Full log: ", LOG_PATH); close(.logcon)
}, add = TRUE)

set.seed(1234)

WIDTH_MM <- 180; HEIGHT_MM <- 228
BASE_PT <- 7; AXIS_PT <- 6; LAB_PT <- 5.5; TAG_PT <- 9

## ---- palette -----------------------------------------------------------------
## One hue, one meaning, in every panel:
##   red    melanoma / malignant cells
##   blue   benign or non-melanoma comparison (nevus, other lineages)
##   gold   vascular and stromal compartments
##   purple immune compartments
##   grey   unassigned
## Checked for colour-vision deficiency rather than chosen by eye: worst
## adjacent separation dE 20.4 (protan), 28.1 (normal vision). The gold falls
## just under 3:1 contrast against white, which is acceptable only because
## every bar carrying it is directly labelled.
COL_MEL    <- "#C0392B"   # melanoma / malignant
COL_NEVUS  <- "#1F6FB2"   # nevus / other lineages
COL_STROMA <- "#C98A2E"   # endothelial, CAF
COL_IMMUNE <- "#7D4F9E"   # B, T, macrophage, NK
COL_LO     <- "#9A9A9A"   # unassigned / background points
COL_PT     <- COL_MEL     # TCGA melanoma tumours

theme_fig <- function() {
  theme_classic(base_size = BASE_PT) +
    theme(axis.text  = element_text(size = AXIS_PT, colour = "black"),
          axis.title = element_text(size = BASE_PT),
          axis.line  = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          strip.background = element_blank(),
          strip.text = element_text(size = AXIS_PT, face = "bold"),
          plot.margin = margin(2, 3, 1, 1, "mm"),
          legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = AXIS_PT),
          legend.background = element_blank())
}

find1 <- function(pat) {
  roots <- c(file.path(ROOT, "results"), file.path(ROOT, "data", "processed"))
  roots <- roots[dir.exists(roots)]
  h <- unlist(lapply(roots, function(r)
    list.files(r, pattern = pat, recursive = TRUE, full.names = TRUE)), use.names = FALSE)
  if (length(h)) h[1] else NA_character_
}

say(""); say("== 1. Inputs ===================================================")
f_vals  <- find1("^nevus_vs_melanoma_values\\.csv$")
f_tests <- find1("^nevus_vs_melanoma_tests\\.csv$")
f_tir   <- find1("^tirosh_detection_by_class\\.csv$")
f_lines <- find1("^46_line_level_data\\.csv$")
f_cache <- file.path(ROOT, "results", "tcga", "fig6_tcga_perpatient.csv")
for (nm in c(f_vals, f_tests, f_tir, f_lines))
  say(sprintf("  %-44s %s", basename(nm), ifelse(is.na(nm), "NOT FOUND", "ok")))
if (any(is.na(c(f_vals, f_tests, f_tir, f_lines))))
  stop("required input missing -- see script 51's log", call. = FALSE)

## ---- 2. Panel D data: cache or rebuild --------------------------------------
say(""); say("== 2. Panel D source ===========================================")

PROMOTER <- c("cg18827107","cg19622623","cg27035734","cg25456477",
              "cg00310523","cg02569458","cg06740227")
TARGET_RHO <- -0.3789      # published value this reconstruction must reproduce

if (file.exists(f_cache)) {
  TD <- fread(f_cache, showProgress = FALSE)
  say("  using cache: ", basename(f_cache), " (", nrow(TD), " patients)")
} else {
  say("  no cache; rebuilding from the RDS objects (this takes a few minutes)")
  f_meth <- find1("^skcm_hm450\\.rds$"); f_expr <- find1("^skcm_star_counts\\.rds$")
  if (is.na(f_meth) || is.na(f_expr)) {
    f_meth <- list.files(ROOT, "^skcm_hm450\\.rds$", recursive = TRUE, full.names = TRUE)[1]
    f_expr <- list.files(ROOT, "^skcm_star_counts\\.rds$", recursive = TRUE, full.names = TRUE)[1]
  }
  if (is.na(f_meth) || is.na(f_expr) || !file.exists(f_meth) || !file.exists(f_expr))
    stop("skcm_hm450.rds / skcm_star_counts.rds not found; panel D cannot be rebuilt",
         call. = FALSE)
  say("    ", f_expr); say("    ", f_meth)

  suppressPackageStartupMessages(library(SummarizedExperiment))

  EX <- readRDS(f_expr)
  gi <- which(rowData(EX)$gene_name == "RASSF9")
  if (!length(gi)) stop("RASSF9 not found in rowData(gene_name)", call. = FALSE)
  if (length(gi) > 1) warn("RASSF9 matched ", length(gi), " rows; using the first")
  e <- log2(assay(EX, "tpm_unstrand")[gi[1], ] + 1)
  EXD <- data.table(sample = names(e), patient = substr(names(e), 1, 12),
                    expr = as.numeric(e),
                    stype = as.character(colData(EX)$sample_type))
  say("    expression samples: ", nrow(EXD), " | patients: ", uniqueN(EXD$patient))
  rm(EX); gc(verbose = FALSE)

  ME <- readRDS(f_meth)
  keep <- rownames(ME) %in% PROMOTER
  say("    promoter probes present in the beta matrix: ", sum(keep), " of ", length(PROMOTER))
  b <- colMeans(assay(ME)[keep, , drop = FALSE], na.rm = TRUE)
  MED <- data.table(sample = names(b), patient = substr(names(b), 1, 12),
                    beta = as.numeric(b))
  say("    methylation samples: ", nrow(MED), " | patients: ", uniqueN(MED$patient))
  rm(ME); gc(verbose = FALSE)

  ## Two candidate collapsing rules; keep whichever reproduces the published rho.
  mk <- function(fun, label) {
    a <- EXD[, .(expr = fun(expr), stype = stype[1]), by = patient]
    m <- MED[, .(beta = fun(beta)), by = patient]
    d <- merge(a, m, by = "patient")
    d <- d[is.finite(expr) & is.finite(beta)]
    r <- suppressWarnings(cor(d$beta, d$expr, method = "spearman"))
    say(sprintf("    rule '%s': n = %d, rho = %+.4f", label, nrow(d), r))
    list(d = d, rho = r, label = label)
  }
  cands <- list(mk(function(v) v[1], "first sample per patient"),
                mk(function(v) mean(v, na.rm = TRUE), "mean across samples"))
  dev <- vapply(cands, function(k) abs(k$rho - TARGET_RHO), numeric(1))
  best <- cands[[which.min(dev)]]
  say(sprintf("    closest to the published %.4f: '%s' (rho %+.4f, n %d)",
              TARGET_RHO, best$label, best$rho, nrow(best$d)))
  if (min(dev) > 0.005)
    stop("neither collapsing rule reproduces the published rho of ", TARGET_RHO,
         " (closest ", round(best$rho, 4), "). Panel D would disagree with the ",
         "Results text; resolve against script 35 before continuing.", call. = FALSE)
  TD <- best$d
  dir.create(dirname(f_cache), recursive = TRUE, showWarnings = FALSE)
  fwrite(TD, f_cache)
  say("    cached to ", basename(f_cache))
}

rhoD <- suppressWarnings(cor(TD$beta, TD$expr, method = "spearman"))
say(sprintf("  panel D: n = %d, Spearman rho = %+.4f", nrow(TD), rhoD))
if (abs(rhoD - TARGET_RHO) > 0.005)
  stop("cached panel D data gives rho = ", round(rhoD, 4), ", not ", TARGET_RHO,
       "; delete the cache and rerun", call. = FALSE)

## ---- 3. Panel A --------------------------------------------------------------
say(""); say("== 3. Panel A ==================================================")
AV <- fread(f_vals, showProgress = FALSE)
AV <- AV[, .(dataset, sample, group, RASSF9)]
AV[, group := factor(gsub("^(.)", "\\U\\1", tolower(group), perl = TRUE),
                     levels = c("Nevus", "Melanoma"))]
say("  ", nrow(AV), " samples across ", uniqueN(AV$dataset), " cohorts")
print_counts <- AV[, .N, by = .(dataset, group)][order(dataset, group)]
for (i in seq_len(nrow(print_counts)))
  say(sprintf("    %-12s %-9s n = %d", print_counts$dataset[i],
              as.character(print_counts$group[i]), print_counts$N[i]))

## Journal notation for p-values: "8.8 × 10⁻⁴", not R's "8.8e-04". Panel E
## already used that form, so without this the two panels disagree in style.
sci_p <- function(p) {
  sup <- c("⁰","¹","²","³","⁴",
           "⁵","⁶","⁷","⁸","⁹")
  vapply(p, function(x) {
    if (!is.finite(x) || x <= 0) return(NA_character_)
    e <- floor(log10(x)); m <- x / 10^e
    d <- strsplit(as.character(abs(e)), "")[[1]]
    sprintf("P = %.1f × 10%s%s", m,
            if (e < 0) "⁻" else "",
            paste(sup[as.integer(d) + 1], collapse = ""))
  }, character(1))
}

TS <- fread(f_tests, showProgress = FALSE)[gene == "RASSF9"]
lab <- TS[, .(dataset, p = p_value)]
lab[, ptxt := sci_p(p)]
AVlab <- merge(unique(AV[, .(dataset)]), lab, by = "dataset")
AVlab <- merge(AVlab, AV[, .(ytop = max(RASSF9)), by = dataset], by = "dataset")

pA <- ggplot(AV, aes(group, RASSF9, fill = group)) +
  geom_violin(linewidth = 0.25, width = 0.85, alpha = 0.45, colour = NA) +
  geom_boxplot(width = 0.17, outlier.shape = NA, linewidth = 0.3, fill = "white") +
  geom_jitter(width = 0.1, height = 0, size = 0.5, alpha = 0.6, colour = "grey25") +
  geom_text(data = AVlab, aes(x = 1.5, y = ytop * 1.1, label = ptxt),
            inherit.aes = FALSE, size = LAB_PT / .pt) +
  facet_wrap(~ dataset, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(Nevus = COL_NEVUS, Melanoma = COL_MEL), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.16))) +
  labs(x = NULL, y = expression(italic("RASSF9")*" expression, log"[2])) +
  theme_fig()

## ---- 4. Panel B --------------------------------------------------------------
say(""); say("== 4. Panel B ==================================================")
B <- TS[, .(dataset, d = cohens_d, lo = d_lower, hi = d_upper,
            n_nevus, n_melanoma)][order(-d)]
for (i in seq_len(nrow(B)))
  say(sprintf("    %-12s d = %+.3f [%+.3f, %+.3f]  (%d vs %d)",
              B$dataset[i], B$d[i], B$lo[i], B$hi[i], B$n_nevus[i], B$n_melanoma[i]))
B[, lab := factor(dataset, levels = dataset)]

pB <- ggplot(B, aes(d, lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0, linewidth = 0.45, colour = COL_MEL) +
  geom_point(size = 1.8, colour = COL_MEL) +
  geom_text(aes(label = sprintf("%.2f [%.2f, %.2f]", d, lo, hi)),
            vjust = -1.3, size = LAB_PT / .pt, colour = "grey25") +
  scale_x_continuous(limits = c(min(B$lo) - 0.25, 0.30)) +
  scale_y_discrete(expand = expansion(add = 0.7)) +
  labs(x = "Cohen's d (melanoma − nevus)", y = NULL) +
  theme_fig()

## ---- 5. Panel C --------------------------------------------------------------
say(""); say("== 5. Panel C ==================================================")
TI <- fread(f_tir, showProgress = FALSE)[gene == "RASSF9"]
say("  scales present: ", paste(unique(TI$scale), collapse = " | "))

## Prefer the depth-matched rows; "all cells" is the unmatched comparison.
dm <- setdiff(unique(TI$scale), "all cells")
if (!length(dm)) {
  warn("no depth-matched scale found; using 'all cells'")
  TIu <- TI[scale == "all cells"]
} else {
  if (length(dm) > 1) warn("several non-'all cells' scales: ", paste(dm, collapse = ", "))
  TIu <- TI[scale == dm[1]]
  say("  using scale: '", dm[1], "'")
}
TIu <- TIu[order(-pct_det)]
for (i in seq_len(nrow(TIu)))
  say(sprintf("    %-14s %5d cells  %4d positive  %6.2f%%",
              TIu$class[i], TIu$n_cells[i], TIu$n_det[i], TIu$pct_det[i]))

## Verify against the manuscript before drawing.
## Match "Malignant" exactly: a loose pattern also catches
## "Non-malignant, unassigned", which both defeats the check and highlights
## the wrong bar.
is_malignant <- function(x) grepl("^\\s*malignant\\s*$", x, ignore.case = TRUE)
mal <- TIu[is_malignant(class)]
end <- TIu[grepl("endothel", class, ignore.case = TRUE)]
if (nrow(mal) == 1) {
  say(sprintf("  malignant: %d of %d = %.2f%%", mal$n_det, mal$n_cells, mal$pct_det))
  if (abs(mal$pct_det - 0.20) > 0.05)
    warn("malignant detection is ", round(mal$pct_det, 2),
         "%, not the 0.20% quoted in the manuscript")
} else warn("could not identify a single malignant class")
if (nrow(end) == 1)
  say(sprintf("  endothelial: %d of %d = %.2f%%", end$n_det, end$n_cells, end$pct_det))

## Colour by compartment rather than highlighting one bar. This is the panel's
## actual finding: RASSF9 sits in vascular and stromal cells, barely in immune
## cells and essentially not in malignant ones. A single highlighted bar shows
## only where it is absent.
TIu[, compartment := fcase(
  is_malignant(class),                                        "Malignant",
  grepl("endothel|CAF|fibro", class, ignore.case = TRUE),     "Vascular / stromal",
  grepl("B cell|T cell|NK|macroph|monocy|mast", class, ignore.case = TRUE), "Immune",
  default                                                      = "Unassigned")]
say("  compartments: ", paste(sprintf("%s=%d", names(table(TIu$compartment)),
                                      table(TIu$compartment)), collapse = ", "))
if (sum(TIu$compartment == "Malignant") != 1)
  warn("expected exactly one malignant class, got ",
       sum(TIu$compartment == "Malignant"))
un <- TIu[compartment == "Unassigned", class]
if (length(un)) say("  classed as Unassigned: ", paste(un, collapse = ", "))

TIu[, compartment := factor(compartment,
      levels = c("Vascular / stromal", "Immune", "Malignant", "Unassigned"))]
TIu[, class := factor(class, levels = rev(class))]

pC <- ggplot(TIu, aes(pct_det, class, fill = compartment)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%% (%d/%d)", pct_det, n_det, n_cells)),
            hjust = -0.08, size = LAB_PT / .pt, colour = "grey25") +
  scale_fill_manual(values = c(`Vascular / stromal` = COL_STROMA,
                               Immune               = COL_IMMUNE,
                               Malignant            = COL_MEL,
                               Unassigned           = COL_LO), name = NULL) +
  scale_x_continuous(limits = c(0, max(TIu$pct_det) * 1.45),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression("Cells with detectable "*italic("RASSF9")*" (%)"), y = NULL) +
  theme_fig() +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.35),
        legend.justification = c(1, 0.5))

## ---- 6. Panel D --------------------------------------------------------------
say(""); say("== 6. Panel D ==================================================")
say(sprintf("  n = %d | rho = %+.3f", nrow(TD), rhoD))

pD <- ggplot(TD, aes(beta, expr)) +
  geom_point(size = 0.55, alpha = 0.45, colour = COL_PT) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "grey20", linewidth = 0.45, fill = "grey80") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5,
           size = LAB_PT / .pt,
           label = sprintf("rho == %.3f", rhoD), parse = TRUE) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 3.1,
           size = LAB_PT / .pt, label = sprintf("n = %d", nrow(TD))) +
  labs(x = "Mean promoter methylation (beta)",
       y = expression(italic("RASSF9")*" log"[2]*"(TPM + 1)")) +
  theme_fig()

## ---- 7. Panels E and F -------------------------------------------------------
say(""); say("== 7. Panels E and F ===========================================")
LL <- fread(f_lines, showProgress = FALSE)
LL[, is_mel := !is.na(disease) & grepl("melanoma", disease, ignore.case = TRUE)]
LL[, is_uv  := !is.na(subtype) & grepl("uveal|ocular", paste(subtype, disease), ignore.case = TRUE)]
LL[, grp := fifelse(is_mel & !is_uv, "Melanoma", "Other lineages")]

## PANEL E USES A DIFFERENT COHORT FROM PANEL F.
## 46_line_level_data.csv holds only lines with BOTH methylation and
## expression (741, of which 39 are melanoma). The methylation comparison in
## the manuscript is over every line with methylation and a lineage
## annotation: 56 melanoma against 940 others, medians 0.588 and 0.183.
## Building panel E from the matched file gives 39 vs 702 and medians that do
## not match the text or the quoted P value, so it is rebuilt here from the
## extracted beta matrix, exactly as script 45 did.
f_rds <- find1("^44b_locus_betas\\.rds$")
if (is.na(f_rds)) stop("44b_locus_betas.rds not found; panel E cannot be built",
                       call. = FALSE)
OBJ    <- readRDS(f_rds)
BM     <- OBJ$beta
PRB    <- as.data.table(OBJ$probes)
core9  <- PRB[gene == "RASSF9" & set == "tss_proximal" & probe %in% rownames(BM), probe]
say("  panel E rebuilt from ", basename(f_rds), " | core probes: ", length(core9))

f_model <- find1("^Model\\.csv$")
if (is.na(f_model)) f_model <- list.files(ROOT, "^Model\\.csv$", recursive = TRUE,
                                          full.names = TRUE)[1]
if (is.na(f_model) || !file.exists(f_model))
  stop("Model.csv not found; panel E cannot be built", call. = FALSE)
MD <- fread(f_model, showProgress = FALSE)
strip <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
ANN <- MD[, .(mk = strip(StrippedCellLineName), lineage = OncotreeLineage,
              disease = OncotreePrimaryDisease, subtype = OncotreeSubtype)]
ANN <- ANN[!duplicated(mk)]

E <- data.table(geo_name = colnames(BM), mk = strip(colnames(BM)),
                core = colMeans(BM[core9, , drop = FALSE], na.rm = TRUE))
E <- merge(E, ANN, by = "mk", all.x = TRUE)
E[, is_mel := !is.na(disease) & grepl("melanoma", disease, ignore.case = TRUE)]
E[, is_uv  := !is.na(subtype) &
     grepl("uveal|ocular", paste(subtype, disease), ignore.case = TRUE)]
E <- E[is.finite(core) & ((is_mel & !is_uv) | (!is_mel & !is.na(lineage)))]
E[, grp := factor(fifelse(is_mel, "Melanoma", "Other lineages"),
                  levels = c("Melanoma", "Other lineages"))]
setnames(E, "core", "RASSF9_core")

n_mel <- sum(E$grp == "Melanoma"); n_oth <- sum(E$grp == "Other lineages")
med_mel <- median(E[grp == "Melanoma", RASSF9_core])
med_oth <- median(E[grp == "Other lineages", RASSF9_core])
say("  panel E: ", n_mel, " melanoma vs ", n_oth, " other lines")
say(sprintf("    medians: %.3f vs %.3f", med_mel, med_oth))

## The quoted P value belongs to this exact comparison, so the panel must be
## plotting it and not some other subset.
if (n_mel != 56 || n_oth != 940)
  stop("panel E cohort is ", n_mel, " vs ", n_oth, ", expected 56 vs 940 -- ",
       "the quoted P = 8.7e-17 would not correspond to what is plotted",
       call. = FALSE)
if (abs(med_mel - 0.588) > 0.002 || abs(med_oth - 0.183) > 0.002)
  warn(sprintf("panel E medians are %.3f / %.3f, manuscript says 0.588 / 0.183",
               med_mel, med_oth))

pE <- ggplot(E, aes(grp, RASSF9_core, fill = grp)) +
  geom_violin(linewidth = 0.25, width = 0.85, alpha = 0.4, colour = NA) +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.3, fill = "white") +
  geom_jitter(width = 0.11, height = 0, size = 0.35, alpha = 0.35, colour = "grey25") +
  annotate("text", x = 1.5, y = 1.02, size = LAB_PT / .pt,
           label = "P = 8.7 × 10⁻¹⁷") +
  scale_fill_manual(values = c(Melanoma = COL_MEL, `Other lineages` = COL_NEVUS), guide = "none") +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = expression(italic("RASSF9")*" promoter-core methylation (beta)")) +
  theme_fig()

FF <- LL[is.finite(RASSF9_core) & is.finite(RASSF9)]
FF[, grp := factor(grp, levels = c("Other lineages", "Melanoma"))]
rhoF <- suppressWarnings(cor(FF$RASSF9_core, FF$RASSF9, method = "spearman"))
say("  panel F: ", nrow(FF), " lines | rho = ", sprintf("%+.3f", rhoF))

pF <- ggplot(FF, aes(RASSF9_core, RASSF9)) +
  geom_point(aes(colour = grp, size = grp, alpha = grp)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "grey20", linewidth = 0.45, fill = "grey80") +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  ## sits clear of the line rather than on it
  annotate("text", x = 0.99, y = 1, label = "TPM = 1", hjust = 1, vjust = -0.9,
           size = LAB_PT / .pt, colour = "grey35") +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.5, size = LAB_PT / .pt,
           label = sprintf("rho == %.3f", rhoF), parse = TRUE) +
  ## Background context, not a compared group: grey keeps the 702 other lines
  ## recessive so the 39 melanoma points are legible. Panel E, where the two
  ## groups ARE compared, keeps both hues.
  scale_colour_manual(values = c(`Other lineages` = COL_LO, Melanoma = COL_MEL), name = NULL) +
  scale_size_manual(values = c(`Other lineages` = 0.4, Melanoma = 1.3), guide = "none") +
  scale_alpha_manual(values = c(`Other lineages` = 0.30, Melanoma = 0.95), guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = expression(italic("RASSF9")*" promoter-core methylation (beta)"),
       y = expression(italic("RASSF9")*" log"[2]*"(TPM + 1)")) +
  theme_fig() +
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.93),
        legend.justification = c(1, 1))

## ---- 8. Assemble -------------------------------------------------------------
say(""); say("== 8. Assembling ===============================================")
say("  ", WIDTH_MM, " x ", HEIGHT_MM, " mm at final size")

FIG <- (pA | pB) / (pC | pD) / (pE | pF) +
  plot_layout(widths = c(1.25, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = TAG_PT, face = "bold", hjust = 0, vjust = 1))

## Windows locks a file that is open in a viewer, so a re-run while the PNG is
## being looked at fails at the very end, after all the work. Save to a
## temporary name and move it into place; if the move is blocked, keep the
## temporary file and say where it is rather than losing the run.
save_fig <- function(path, ...) {
  tmp <- file.path(dirname(path),
                   paste0(".tmp_", basename(path)))
  ok <- tryCatch({ ggsave(tmp, FIG, ...); TRUE },
                 error = function(e) { warn("could not write ", basename(tmp),
                                            ": ", conditionMessage(e)); FALSE })
  if (!ok) return(invisible(FALSE))
  if (file.exists(path) && !file.remove(path)) {
    warn(basename(path), " is locked (open in a viewer?). The new figure is at ",
         basename(tmp), " -- close the viewer and rename it, or re-run.")
    return(invisible(FALSE))
  }
  if (!file.rename(tmp, path)) {
    warn("could not rename ", basename(tmp), " to ", basename(path))
    return(invisible(FALSE))
  }
  say("  wrote ", basename(path))
  invisible(TRUE)
}

save_fig(file.path(FIG_DIR, "Figure6.pdf"),
         width = WIDTH_MM, height = HEIGHT_MM, units = "mm", device = cairo_pdf)
save_fig(file.path(FIG_DIR, "Figure6.png"),
         width = WIDTH_MM, height = HEIGHT_MM, units = "mm", dpi = 600)

say("")
say("  Place at 100%. The type sizes written here are the printed sizes;")
say("  scaling the figure on a slide reduces them below the journal minimum.")
say("")
say("  WHAT STILL NEEDS A HUMAN: panel C's cell-class labels come straight")
say("  from the Tirosh annotation. Check they read the way the manuscript")
say("  describes them before submission.")
