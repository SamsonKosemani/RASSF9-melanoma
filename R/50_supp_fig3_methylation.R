## =============================================================================
## 50_supp_fig3_methylation.R
##
## Builds Supplementary Figure 3: "Robustness and positional specificity of
## RASSF9 promoter methylation" (TCGA-SKCM).
##
##   A  Correlation with expression by CpG position relative to the TSS
##   B  Covariate adjustment ladder, M0-M4
##   C  Each probe against the genome-wide null distribution
##
## WHY THESE THREE
##   The manuscript claims the methylation-expression association is (i) real,
##   (ii) robust to covariates and (iii) specific to the promoter core. Panel A
##   shows the positional structure including the shore inversion, panel B the
##   robustness, panel C the specificity against an empirical null. Nothing in
##   the figure is new analysis; every value is read from the existing result
##   files and the script verifies them against the manuscript before plotting.
##
## PRODUCTION NOTES
##   The existing deck is 275 mm wide with 5.8 pt labels. Reproduced at a
##   journal's 180 mm column that becomes about 3.8 pt, below any publisher
##   minimum. This figure is therefore drawn AT FINAL SIZE (180 mm) so the
##   point sizes written here are the point sizes that print. Nothing is
##   scaled afterwards.
##
##   Outputs a vector PDF for submission and a 600 dpi PNG for circulation.
##
## INPUTS  (located automatically; see section 1)
##   methylation_full_probes.csv        probe, dist_tss, region, rho, median_beta
##   methylation_specificity.csv        gene, probe_set, model, n, beta, lo, hi, p
##   methylation_null_distribution.csv  gene, probe, set, rho, null_pctile
##
## OUTPUTS
##   figures/tcga/SupplementaryFigure3.pdf   vector, for submission
##   figures/tcga/SupplementaryFigure3.png   600 dpi
##   results/logs/50_supp_fig3_methylation.log
## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------
if (file.exists("R/00_setup.R")) source("R/00_setup.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT    <- if (exists("PROJECT_ROOT")) PROJECT_ROOT else getwd()
FIG_DIR <- file.path(ROOT, "figures", "tcga")
LOG_DIR <- file.path(ROOT, "results", "logs")
for (d in c(FIG_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_PATH <- file.path(LOG_DIR, "50_supp_fig3_methylation.log")
.logcon  <- file(LOG_PATH, open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = .logcon) }
WARN <- character(0)
warn <- function(...) { m <- paste0(...); WARN <<- c(WARN, m); say("  [WARN] ", m) }
on.exit({
  if (length(WARN)) { say(""); say("WARNINGS RAISED (", length(WARN), "):"); for (w in WARN) say("  - ", w) }
  else say("  no warnings raised")
  say(""); say("Full log: ", LOG_PATH); close(.logcon)
}, add = TRUE)

## Final-size geometry. Change WIDTH_MM only if the journal column differs.
WIDTH_MM  <- 180
HEIGHT_MM <- 165
BASE_PT   <- 7      # body text
AXIS_PT   <- 6      # tick labels; nothing smaller than this anywhere
LAB_PT    <- 5.5    # probe identifiers
TAG_PT    <- 9      # panel letters

COL_CORE   <- "#2C6E8E"
COL_SHORE  <- "#D1702F"
COL_DISTAL <- "#9A9A9A"

## ggplot2 3.5 deprecated the numeric form of legend.position in favour of
## legend.position.inside, and 4.0 no longer accepts it. This returns the
## right pair of theme arguments for whichever version is installed, so the
## script runs on both without a warning or an error.
inside_legend <- function(x, y, just) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    list(legend.position = "inside",
         legend.position.inside = c(x, y),
         legend.justification = just)
  } else {
    list(legend.position = c(x, y), legend.justification = just)
  }
}

## ---- 1. Locate the inputs ---------------------------------------------------
say("")
say("== 1. Inputs ===================================================")

find_one <- function(fname) {
  cand <- c(
    file.path(ROOT, "results", "tcga", fname),
    file.path(ROOT, "results", fname),
    file.path(ROOT, "results", "geo", fname),
    file.path(ROOT, fname)
  )
  hit <- cand[file.exists(cand)]
  if (length(hit)) return(hit[1])
  ## last resort: search the results tree
  deep <- list.files(file.path(ROOT, "results"), pattern = paste0("^", fname, "$"),
                     recursive = TRUE, full.names = TRUE)
  if (length(deep)) return(deep[1])
  NA_character_
}

f_probes <- find_one("methylation_full_probes.csv")
f_models <- find_one("methylation_specificity.csv")
f_null   <- find_one("methylation_null_distribution.csv")

found <- c("methylation_full_probes.csv"       = f_probes,
           "methylation_specificity.csv"       = f_models,
           "methylation_null_distribution.csv" = f_null)
for (nm in names(found))
  say(sprintf("  %-36s %s", nm, ifelse(is.na(found[[nm]]), "NOT FOUND", found[[nm]])))
if (any(is.na(c(f_probes, f_models, f_null)))) {
  say("")
  say("  Place the missing file(s) under results/tcga/ (or anywhere under")
  say("  results/) and re-run. They are the outputs of scripts 35 and 36.")
  stop("required input missing", call. = FALSE)
}

PR <- fread(f_probes, showProgress = FALSE)
MD <- fread(f_models, showProgress = FALSE)
NU <- fread(f_null,   showProgress = FALSE)

## ---- 2. Verify against the manuscript ---------------------------------------
## The figure must not quietly disagree with the text. These are the values the
## manuscript states; a mismatch stops the script rather than producing a
## figure that contradicts the paper.
say("")
say("== 2. Verifying values against the manuscript ==================")

expect <- data.table(
  probe = c("cg18827107","cg19622623","cg27035734","cg25456477",
            "cg00310523","cg02569458","cg06740227"),
  rho   = c( 0.135, 0.136, -0.368, -0.485, -0.429, -0.395, -0.294),
  dist  = c(-770, -638, -370, -180, 11, 94, 383))

P9 <- PR[gene == "RASSF9"]
bad <- 0L
for (i in seq_len(nrow(expect))) {
  row <- P9[probe == expect$probe[i]]
  if (!nrow(row)) { warn("probe ", expect$probe[i], " absent from the probe file"); bad <- bad + 1L; next }
  dr <- abs(round(row$rho[1], 3) - expect$rho[i])
  dd <- abs(row$dist_tss[1] - expect$dist[i])
  ok <- dr < 0.0015 && dd == 0
  say(sprintf("  %-12s dist %+5d  rho %+.3f   %s",
              expect$probe[i], row$dist_tss[1], row$rho[1], ifelse(ok, "ok", "*** MISMATCH")))
  if (!ok) bad <- bad + 1L
}
if (bad) stop(bad, " probe value(s) disagree with the manuscript -- resolve before plotting", call. = FALSE)
say("  all 7 promoter probes match the manuscript")

## ---- 3. Panel A : correlation by CpG position --------------------------------
say("")
say("== 3. Panel A ==================================================")

A <- P9[probe %in% expect$probe][order(dist_tss)]
A[, cls := ifelse(abs(dist_tss) <= 400, "Promoter core", "Promoter shore")]
A[, cls := factor(cls, levels = c("Promoter shore", "Promoter core"))]
say("  core probes: ", sum(A$cls == "Promoter core"),
    " | shore probes: ", sum(A$cls == "Promoter shore"))

## Probe identifiers are drawn vertically along the top. The closest pair is
## 83 bp apart, which on a ~1,250 bp axis across 180 mm is about 11 mm of
## separation -- far more than a rotated 5.5 pt label needs, so no repulsion
## package is required and nothing can overlap.
LAB_Y <- 0.26

## Value labels sit below negative points and above positive ones. Two core
## probes are only 83 bp apart (+11 and +94), so their labels would nearly
## touch at a single offset; any probe within 150 bp of its neighbour is
## pushed to a second row. Computed, not eyeballed, so it holds if the probe
## set ever changes.
A <- A[order(dist_tss)]
A[, gap_prev := c(Inf, diff(dist_tss))]
A[, row2 := gap_prev < 150]
A[, vj := fifelse(rho > 0, -1.1, fifelse(row2, 3.3, 1.9))]

pA <- ggplot(A, aes(x = dist_tss, y = rho, colour = cls)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey30") +
  geom_vline(xintercept = 0, linewidth = 0.25, linetype = "dotted", colour = "grey55") +
  geom_segment(aes(xend = dist_tss, yend = 0), linewidth = 0.45) +
  geom_point(size = 1.5) +
  geom_text(aes(y = LAB_Y, label = probe), angle = 90, hjust = 0,
            size = LAB_PT / .pt, colour = "grey30", show.legend = FALSE) +
  geom_text(aes(label = sprintf("%+.3f", rho), vjust = vj),
            size = LAB_PT / .pt, show.legend = FALSE) +
  scale_colour_manual(values = c("Promoter shore" = COL_SHORE,
                                 "Promoter core"  = COL_CORE), name = NULL) +
  scale_x_continuous(breaks = seq(-800, 400, 200), expand = expansion(mult = 0.06)) +
  ## Explicit labels. seq(-0.6, 0.6, 0.2) returns 1.110223e-16 rather than an
  ## exact zero, and ggplot's default labeller responds by rendering the whole
  ## axis in scientific notation. Supplying the labels avoids that entirely.
  scale_y_continuous(limits = c(-0.62, 0.62),
                     breaks = seq(-0.6, 0.6, 0.2),
                     labels = c("-0.6", "-0.4", "-0.2", "0", "0.2", "0.4", "0.6")) +
  labs(x = "CpG distance to RASSF9 transcription start site (bp)",
       y = expression("Spearman "*rho*" (methylation vs expression)")) +
  theme_classic(base_size = BASE_PT) +
  do.call(theme, c(inside_legend(0.015, 0.06, c(0, 0)), list(
    legend.key.size = unit(3, "mm"), legend.text = element_text(size = AXIS_PT),
    legend.background = element_blank(),
    axis.text  = element_text(size = AXIS_PT, colour = "black"),
    axis.title = element_text(size = BASE_PT),
    axis.line  = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(2, 3, 1, 1, "mm"))))

## ---- 4. Panel B : covariate adjustment --------------------------------------
say("")
say("== 4. Panel B ==================================================")

B <- MD[gene == "RASSF9" & probe_set == "promoter"]
if (!nrow(B)) stop("no RASSF9 promoter rows in the specificity file", call. = FALSE)

## Short labels; the full model definitions go in the legend text, not on the
## axis, so nothing is truncated.
short <- c("M0 unadjusted"                = "M0  unadjusted",
           "M1 + sample type"             = "M1  + sample type",
           "M2 + global methylation"      = "M2  + global methylation",
           "M3 + global + sample type"    = "M3  + global + sample type",
           "M4 + global + type + purity"  = "M4  + global + type + purity")
B[, lab := ifelse(model %in% names(short), short[model], model)]
B[, lab := factor(lab, levels = rev(unique(lab)))]

say(sprintf("  %-30s %5s %8s  [%s]", "model", "n", "beta", "95% CI"))
for (i in seq_len(nrow(B)))
  say(sprintf("  %-30s %5d %+8.3f  [%+.3f, %+.3f]", B$model[i], B$n[i], B$beta[i], B$lo[i], B$hi[i]))

## The x range must include zero: a forest plot whose null line is cropped out
## cannot show that every interval excludes it, which is the point of the
## panel. The earlier version clipped the line away and ggplot warned about it.
x_lo  <- min(B$lo) - 0.08 * diff(range(c(B$lo, 0)))
x_hi  <- 0 + 0.17 * diff(range(c(B$lo, 0)))   # room for the n labels
x_lab <- 0 + 0.03 * diff(range(c(B$lo, 0)))

pB <- ggplot(B, aes(x = beta, y = lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  ## geom_errorbarh is deprecated in ggplot2 4.0; geom_errorbar with
  ## orientation = "y" is the supported form and behaves identically here.
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0, linewidth = 0.45, colour = COL_CORE) +
  geom_point(size = 1.5, colour = COL_CORE) +
  geom_text(x = x_lab, aes(label = paste0("n = ", n)),
            hjust = 0, size = LAB_PT / .pt, colour = "grey30") +
  scale_x_continuous(limits = c(x_lo, x_hi), breaks = scales::pretty_breaks(5)) +
  labs(x = expression("Regression coefficient "*beta*" [95% CI]"), y = NULL) +
  theme_classic(base_size = BASE_PT) +
  theme(
    axis.text.y  = element_text(size = AXIS_PT, colour = "black", hjust = 0),
    axis.text.x  = element_text(size = AXIS_PT, colour = "black"),
    axis.title.x = element_text(size = BASE_PT),
    axis.line    = element_line(linewidth = 0.3),
    axis.ticks   = element_line(linewidth = 0.3),
    plot.margin  = margin(2, 2, 1, 1, "mm"))

## ---- 5. Panel C : against the empirical null ---------------------------------
say("")
say("== 5. Panel C ==================================================")

C <- NU[gene == "RASSF9"]
if (!nrow(C)) stop("no RASSF9 rows in the null-distribution file", call. = FALSE)

## The `set` labels in this file are "TSS-proximal", "promoter" and "distal".
## Match them case-insensitively and ignoring punctuation, because an exact
## string test here silently relabels every core probe as distal and inverts
## the specificity claim the panel exists to make.
set_key <- tolower(gsub("[^A-Za-z]", "", C$set))
C[, cls := fcase(
  grepl("tssproximal", set_key), "Promoter core",
  set_key == "promoter",         "Promoter shore",
  default                        = "Distal")]
C[, cls := factor(cls, levels = c("Promoter core", "Promoter shore", "Distal"))]

## Verify the classification before anything is drawn. Panel A checks its
## values against the manuscript; this panel had no such check, which is how
## the mismatch above went unnoticed.
n_core  <- sum(C$cls == "Promoter core")
n_shore <- sum(C$cls == "Promoter shore")
say("  classified: ", n_core, " core, ", n_shore, " shore, ",
    sum(C$cls == "Distal"), " distal")
say("  distinct set labels in the file: ", paste(sort(unique(C$set)), collapse = ", "))
if (n_core != 5 || n_shore != 2)
  stop("expected 5 core and 2 shore probes, got ", n_core, " and ", n_shore,
       " -- the set labels in ", basename(f_null), " do not match what this ",
       "script expects; fix the matching rather than the data", call. = FALSE)

## The core probes must be the ones in the extreme lower tail. If they are not,
## something upstream is wrong and the panel would misrepresent the result.
if (max(C[cls == "Promoter core", null_pctile]) > 5)
  warn("a promoter-core probe sits above the 5th percentile of the null -- ",
       "check before using this panel")
C <- C[order(null_pctile)]
C[, probe := factor(probe, levels = probe)]

n_null <- unique(C$null_n)
say("  probes plotted: ", nrow(C), " | null reference: ", n_null[1], " probes")
say(sprintf("  %-12s %-15s %8s %10s", "probe", "class", "rho", "percentile"))
for (i in seq_len(nrow(C)))
  say(sprintf("  %-12s %-15s %+8.3f %9.2f%%",
              as.character(C$probe[i]), as.character(C$cls[i]), C$rho[i], C$null_pctile[i]))

pC <- ggplot(C, aes(x = null_pctile, y = probe, colour = cls)) +
  geom_vline(xintercept = 5, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_segment(aes(xend = 0, yend = probe), linewidth = 0.45) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = c("Promoter core"  = COL_CORE,
                                 "Promoter shore" = COL_SHORE,
                                 "Distal"         = COL_DISTAL), name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0.01, 0.04))) +
  labs(x = paste0("Percentile of the null distribution (", n_null[1], " probes)"),
       y = NULL) +
  theme_classic(base_size = BASE_PT) +
  do.call(theme, c(inside_legend(0.98, 0.04, c(1, 0)), list(
    legend.key.size = unit(3, "mm"), legend.text = element_text(size = AXIS_PT),
    legend.background = element_blank(),
    axis.text.y  = element_text(size = AXIS_PT, colour = "black"),
    axis.text.x  = element_text(size = AXIS_PT, colour = "black"),
    axis.title.x = element_text(size = BASE_PT),
    axis.line    = element_line(linewidth = 0.3),
    axis.ticks   = element_line(linewidth = 0.3),
    plot.margin  = margin(2, 2, 1, 1, "mm"))))

## ---- 6. Assemble -------------------------------------------------------------
say("")
say("== 6. Assembling at final size =================================")
say("  ", WIDTH_MM, " x ", HEIGHT_MM, " mm; type sizes as written, no rescaling")

have_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (!have_patchwork) {
  say("  patchwork not installed -- attempting install")
  try(install.packages("patchwork", repos = "https://cloud.r-project.org", quiet = TRUE),
      silent = TRUE)
  have_patchwork <- requireNamespace("patchwork", quietly = TRUE)
}
if (!have_patchwork) stop("patchwork is required to assemble the panels", call. = FALSE)

library(patchwork)
FIG <- (pA / (pB | pC)) +
  plot_layout(heights = c(0.88, 1.12)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = TAG_PT, face = "bold", hjust = 0, vjust = 1))

out_pdf <- file.path(FIG_DIR, "SupplementaryFigure3.pdf")
out_png <- file.path(FIG_DIR, "SupplementaryFigure3.png")

ggsave(out_pdf, FIG, width = WIDTH_MM, height = HEIGHT_MM, units = "mm",
       device = cairo_pdf)
ggsave(out_png, FIG, width = WIDTH_MM, height = HEIGHT_MM, units = "mm", dpi = 600)

say("  wrote ", basename(out_pdf), " (vector, for submission)")
say("  wrote ", basename(out_png), " (600 dpi)")

## ---- 7. The legend, written from the plotted values --------------------------
say("")
say("== 7. Draft legend =============================================")
say("")
m0 <- B[grepl("^M0", model)]; m4 <- B[grepl("^M4", model)]
core_pct <- C[cls == "Promoter core", max(null_pctile)]
shore_pct <- C[cls == "Promoter shore", min(null_pctile)]

say("Supplementary Figure 3. Robustness and positional specificity of the")
say("RASSF9 promoter methylation-expression association in TCGA-SKCM.")
say("(A) Spearman correlation between methylation beta and RASSF9 expression")
say("for each of the seven probes in the promoter window, plotted against")
say("strand-aware distance to the transcription start site. RASSF9 lies on the")
say("minus strand of chromosome 12, so negative values are upstream. The five")
say("TSS-proximal probes that constitute the promoter core correlate")
say("negatively, whereas the two upstream shore probes correlate positively,")
say("the island-shore inversion. (B) Regression coefficient for mean")
say("promoter-region methylation under successive adjustment: M0 unadjusted;")
say("M1 plus sample type; M2 plus global methylation; M3 plus global")
say("methylation and sample type; M4 plus tumour purity. Points are point")
say("estimates, bars 95% confidence intervals, and n is given for each model;")
say("purity estimates are available for a subset, hence the smaller n for M4.")
say(sprintf("The coefficient does not attenuate under adjustment (%+.3f at M0 to %+.3f at M4).",
            m0$beta[1], m4$beta[1]))
say("(C) Each probe at the locus placed against an empirical null built from")
say(sprintf("%d probes outside the RASSF9 and RASSF10 loci; the dashed line marks", n_null[1]))
say("the fifth percentile. All five promoter-core probes fall in the extreme")
say(sprintf("lower tail (below the %.1f percentile), the distal probes lie within the", core_pct))
say(sprintf("body of the null, and the two shore probes sit near the %.0f percentile", shore_pct))
say("in the opposite tail. Correlations are two-sided Spearman; n = 469 tumours")
say("with paired RNA-seq and HM450 methylation.")
say("")
say("  NOTE: the figure is drawn at ", WIDTH_MM, " mm, the final printed width, so")
say("  the smallest type is ", LAB_PT, " pt as printed. Do not rescale it on the slide;")
say("  place it at 100%.")
