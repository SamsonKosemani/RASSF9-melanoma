## =============================================================================
## 53_supp_visium_figure.R
##
## Supplementary figure for the GSE298774 spatial analysis, moved out of
## Figure 6 by the decision of 3 November 2026.
##
##   A  RASSF9 by tissue compartment (the lineage result)
##   B  Paired nevus vs melanoma compartments within the same lesion
##   C  Bulk and paired designs on one fold-change scale
##
## WHY REDRAWN RATHER THAN COMPOSITED
##   figures/geo/ already holds geo_rassf9_visium_lineage, _paired and
##   _bulk_vs_visium as finished images. Pasting them together would carry
##   three different canvas sizes and three different type sizes into one
##   figure: the lineage chart alone was drawn with roughly 14 pt labels, which
##   land near 4 pt once it is shrunk to a panel. Embedding also replaces
##   vector text with pixels. All three panels are redrawn here from the source
##   CSVs at one consistent size instead.
##
##   fig_visium_main_B4 (the spatial map) is NOT included. It is already a
##   complete three-panel figure with its own scale bars, and reducing it to a
##   corner of this one would make the spots and the 1 mm bars unreadable. It
##   should stay a supplementary figure in its own right.
##
## COLUMN NAMES ARE DETECTED, NOT ASSUMED
##   I have not seen inside the two GSE298774 result files. Every column is
##   located by pattern, the detected mapping is printed, and the script stops
##   with the available names if anything cannot be found. Guessing column
##   names has broken three scripts in this project already.
##
## OUTPUTS
##   figures/geo/SupplementaryFigure_Visium.pdf   vector
##   figures/geo/SupplementaryFigure_Visium.png   600 dpi
##   results/logs/53_supp_visium_figure.log
## =============================================================================

if (file.exists("R/00_setup.R")) source("R/00_setup.R")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

ROOT    <- if (exists("PROJECT_ROOT")) PROJECT_ROOT else getwd()
FIG_DIR <- file.path(ROOT, "figures", "geo")
LOG_DIR <- file.path(ROOT, "results", "logs")
for (d in c(FIG_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_PATH <- file.path(LOG_DIR, "53_supp_visium_figure.log")
.logcon  <- file(LOG_PATH, open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = .logcon) }
WARN <- character(0)
warn <- function(...) { m <- paste0(...); WARN <<- c(WARN, m); say("  [WARN] ", m) }
on.exit({
  if (length(WARN)) { say(""); say("WARNINGS (", length(WARN), "):"); for (w in WARN) say("  - ", w) }
  else say("  no warnings raised")
  say(""); say("Full log: ", LOG_PATH); close(.logcon)
}, add = TRUE)

WIDTH_MM <- 180; HEIGHT_MM <- 150
BASE_PT <- 7; AXIS_PT <- 6; LAB_PT <- 5.5; TAG_PT <- 9

## Same scheme as Figure 6: red melanoma, blue benign/comparison, grey other.
COL_MEL <- "#C0392B"; COL_NEV <- "#1F6FB2"; COL_LO <- "#9A9A9A"
## Panel C needs its own two colours. Red already means "melanocytic" in A and
## "melanoma" in B, and blue means "nevus" in B, so reusing either for a STUDY
## DESIGN in C would give one colour three meanings inside one figure.
COL_INK <- "#2A2A2A"; COL_PAIR <- "#7D4F9E"

theme_fig <- function() {
  theme_classic(base_size = BASE_PT) +
    theme(axis.text  = element_text(size = AXIS_PT, colour = "black"),
          axis.title = element_text(size = BASE_PT),
          axis.line  = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          plot.margin = margin(2, 3, 1, 1, "mm"),
          legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = AXIS_PT),
          legend.background = element_blank())
}

find1 <- function(pat) {
  h <- list.files(file.path(ROOT, "results"), pattern = pat,
                  recursive = TRUE, full.names = TRUE)
  if (length(h)) h[1] else NA_character_
}

## Locate a column by trying patterns in order; report, and stop clearly.
pick <- function(dt, patterns, what, file) {
  for (p in patterns) {
    hit <- grep(p, names(dt), ignore.case = TRUE, value = TRUE)
    if (length(hit)) return(hit[1])
  }
  say("")
  say("  CANNOT FIND a column for ", what, " in ", basename(file))
  say("  columns present: ", paste(names(dt), collapse = ", "))
  stop("column for '", what, "' not found in ", basename(file),
       " -- send me the column list and I will map it explicitly", call. = FALSE)
}

say(""); say("== 1. Inputs ===================================================")
f_lin  <- find1("^gse298774_authorlabels_lineage\\.csv$")
f_pat  <- find1("^gse298774_authorlabels_perpatient\\.csv$")
f_test <- find1("^gse298774_authorlabels_tests\\.csv$")
f_cmp  <- find1("^bulk_vs_visium_comparison\\.csv$")
for (f in c(f_lin, f_pat, f_test, f_cmp))
  say(sprintf("  %-46s %s", ifelse(is.na(f), "(missing)", basename(f)),
              ifelse(is.na(f), "NOT FOUND", "ok")))
if (any(is.na(c(f_lin, f_pat, f_cmp))))
  stop("required input missing", call. = FALSE)

## ---- 2. Panel A : compartment lineage ---------------------------------------
say(""); say("== 2. Panel A : RASSF9 by compartment ==========================")
LIN <- fread(f_lin, showProgress = FALSE)
say("  columns: ", paste(names(LIN), collapse = ", "))

if ("gene" %in% names(LIN)) {
  LIN <- LIN[grepl("^RASSF9$", gene, ignore.case = TRUE)]
  say("  filtered to RASSF9: ", nrow(LIN), " rows")
}
c_name <- pick(LIN, c("^compartment$", "^label$", "^class$", "^topic", "^name$",
                      "MaxTopicsName", "celltype", "cell_type"), "compartment name", f_lin)
c_val  <- pick(LIN, c("^cpm$", "rassf9.*cpm", "cpm.*rassf9", "^value$",
                      "^mean_cpm$", "^expr"), "RASSF9 CPM", f_lin)
say("  using: name = '", c_name, "', value = '", c_val, "'")

A <- LIN[, .(compartment = as.character(get(c_name)), cpm = as.numeric(get(c_val)))]
A <- A[is.finite(cpm)][order(cpm)]
A[, melanocytic := grepl("nevus|tumou?r", compartment, ignore.case = TRUE)]
A[, compartment := factor(compartment, levels = compartment)]
say("  compartments: ", nrow(A), " | melanocytic: ", sum(A$melanocytic))
for (i in rev(seq_len(nrow(A))))
  say(sprintf("    %-24s %7.2f CPM%s", as.character(A$compartment[i]), A$cpm[i],
              ifelse(A$melanocytic[i], "   <- melanocytic", "")))

## Two different quantities, often confused:
##   mean of rates   sum of CPMs, melanocytic / all -- treats a tiny compartment
##                   and a huge one as equals, so it is NOT a share of molecules
##   share of counts sum of (CPM x total UMI) -- the actual fraction of RASSF9
##                   molecules sitting in melanocytic spots
pct_rate <- 100 * sum(A[melanocytic == TRUE, cpm]) / sum(A$cpm)
say(sprintf("  share of summed CPM in melanocytic compartments: %.1f%%", pct_rate))

u_col <- grep("^total_umi$|umi", names(LIN), ignore.case = TRUE, value = TRUE)
pct_mol <- NA_real_
if (length(u_col)) {
  M <- LIN[, .(compartment = as.character(get(c_name)),
               cpm = as.numeric(get(c_val)),
               umi = as.numeric(get(u_col[1])))]
  M <- M[is.finite(cpm) & is.finite(umi)]
  M[, mol := cpm * umi / 1e6]
  M[, melanocytic := grepl("nevus|tumou?r", compartment, ignore.case = TRUE)]
  pct_mol <- 100 * sum(M[melanocytic == TRUE, mol]) / sum(M$mol)
  say(sprintf("  share of RASSF9 MOLECULES in melanocytic compartments: %.1f%%", pct_mol))
  say(sprintf("    (estimated as CPM x %s / 1e6, summed)", u_col[1]))
} else {
  warn("no UMI column found, so the molecule share cannot be computed")
}

say("")
say("  The manuscript quotes 17.1%. Of the two definitions:")
say(sprintf("    summed CPM      %.1f%%  (sum of rates -- not a share of anything)", pct_rate))
say(sprintf("    molecule share  %s  (the quantity the wording implies)",
            ifelse(is.finite(pct_mol), sprintf("%.1f%%", pct_mol), "not computable")))
if (is.finite(pct_mol) && abs(pct_mol - 17.1) <= 0.2) {
  say("  The molecule share MATCHES the published value. The manuscript")
  say("  sentence is correct as written and panel A can quote it.")
} else if (is.finite(pct_mol) && abs(pct_mol - 17.1) > 1) {
  warn(sprintf("molecule share is %.1f%%, manuscript says 17.1%% -- unresolved; ",
               pct_mol),
       "trace the published value to the script that produced it before ",
       "this figure or that sentence is submitted")
}

pA <- ggplot(A, aes(cpm, compartment, fill = melanocytic)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c(`FALSE` = COL_LO, `TRUE` = COL_MEL),
                    labels = c("Other tissue", "Melanocytic"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression(italic("RASSF9")*" (counts per million)"), y = NULL) +
  theme_fig() +
  theme(legend.position = "inside", legend.position.inside = c(0.97, 0.12),
        legend.justification = c(1, 0))

## ---- 3. Panel B : paired within-lesion --------------------------------------
say(""); say("== 3. Panel B : paired comparison ==============================")
PP <- fread(f_pat, showProgress = FALSE)
say("  columns: ", paste(names(PP), collapse = ", "))

p_id <- pick(PP, c("^patient$", "^section$", "^sample$", "^donor$", "library"),
             "patient identifier", f_pat)

## Either long (one row per patient per compartment) or wide (one row per
## patient with a nevus and a melanoma column). Handle both.
long_grp <- grep("^(grp|compartment|group|class|label)$", names(PP),
                 ignore.case = TRUE, value = TRUE)
if (length(long_grp)) {
  ## CPM first: "RASSF9_n" (raw counts) sorts before "RASSF9_n_cpm", so a
  ## loose pattern would silently plot counts instead of a rate.
  p_val <- pick(PP, c("rassf9.*cpm", "^cpm$", "^rassf9$", "^value$"),
                "RASSF9 CPM", f_pat)
  s_col <- grep("^n_spots$", names(PP), ignore.case = TRUE, value = TRUE)
  B <- PP[, .(patient = as.character(get(p_id)),
              grp = tolower(as.character(get(long_grp[1]))),
              cpm = as.numeric(get(p_val)),
              n_spots = if (length(s_col)) as.numeric(get(s_col[1])) else NA_real_)]
  B <- B[grepl("nevus|tumou?r|melanoma", grp)]
  B[, grp := fifelse(grepl("nevus", grp), "Nevus", "Melanoma")]
  say("  long format via '", long_grp[1], "'")
} else {
  nv <- pick(PP, c("nevus"), "nevus value", f_pat)
  ml <- pick(PP, c("melanoma", "tumou?r"), "melanoma value", f_pat)
  B <- rbind(
    PP[, .(patient = as.character(get(p_id)), grp = "Nevus",    cpm = as.numeric(get(nv)))],
    PP[, .(patient = as.character(get(p_id)), grp = "Melanoma", cpm = as.numeric(get(ml)))])
  say("  wide format via '", nv, "' and '", ml, "'")
}
B <- B[is.finite(cpm)]
B[, grp := factor(grp, levels = c("Nevus", "Melanoma"))]
say("  patients in the file: ", uniqueN(B$patient))

## Script 30 applies a pre-specified minimum of 10 spots in EACH compartment
## (its MIN_SPOTS, set at line 92, enforced in paired_frame at line 235). The
## per-patient CSV is written at line 388, BEFORE that filter, so the file
## carries 17 patients while the published test carries 16: A24 has 9 nevus
## spots against 1,860 melanoma spots. Applying the same filter here makes
## this panel the cohort that was actually tested, rather than a 17th patient
## the manuscript never included.
MIN_SPOTS <- 10L
if ("n_spots" %in% names(B) && all(is.finite(B$n_spots))) {
  sp   <- dcast(B, patient ~ grp, value.var = "n_spots")
  keep <- sp[is.finite(Nevus) & is.finite(Melanoma) &
               Nevus >= MIN_SPOTS & Melanoma >= MIN_SPOTS, patient]
  dropped <- setdiff(unique(B$patient), keep)
  if (length(dropped)) {
    for (px in dropped)
      say(sprintf("  excluded by MIN_SPOTS = %d: %-6s (%.0f nevus, %.0f melanoma spots)",
                  MIN_SPOTS, px, sp[patient == px, Nevus],
                  sp[patient == px, Melanoma]))
  } else {
    say("  MIN_SPOTS = ", MIN_SPOTS, " excludes no patient")
  }
  B <- B[patient %in% keep]
} else {
  warn("no usable n_spots column, so MIN_SPOTS = ", MIN_SPOTS,
       " could not be applied; this cohort may not be the published one")
}

n_pat <- uniqueN(B$patient)
say("  patients after the spot filter: ", n_pat)
if (n_pat != 16)
  warn("expected 16 paired patients after MIN_SPOTS, found ", n_pat)

wide <- dcast(B, patient ~ grp, value.var = "cpm")
wide <- wide[is.finite(Nevus) & is.finite(Melanoma)]
wt <- suppressWarnings(stats::wilcox.test(wide$Melanoma, wide$Nevus, paired = TRUE))
say(sprintf("  medians: nevus %.2f, melanoma %.2f | lower in melanoma %d of %d",
            median(wide$Nevus), median(wide$Melanoma),
            sum(wide$Melanoma < wide$Nevus), nrow(wide)))
say(sprintf("  recomputed:  V = %g, p = %.4f", unname(wt$statistic), wt$p.value))

## The panel quotes script 30's own stored test rather than a recomputation,
## so the figure cannot drift from the manuscript. The recomputed value is
## logged above as a cross-check and a mismatch is reported, not overridden.
P_PUB <- NA_real_; N_PUB <- NA_integer_
if (!is.na(f_test)) {
  TT <- fread(f_test, showProgress = FALSE)
  rw <- TT[grepl("^RASSF9$", gene, ignore.case = TRUE) &
             grepl("^all labelled spots$", set, ignore.case = TRUE)]
  if (nrow(rw) == 1L) {
    P_PUB <- as.numeric(rw$wilcoxon_p[1])
    N_PUB <- as.integer(rw$n_patients[1])
    say(sprintf("  published:   V = %g, p = %.4f, n = %d",
                as.numeric(rw$wilcoxon_V[1]), P_PUB, N_PUB))
  }
}
if (is.na(P_PUB)) {
  warn("the published paired test could not be read from the tests table; ",
       "the panel falls back to the recomputed value")
  P_PUB <- wt$p.value; N_PUB <- nrow(wide)
} else {
  if (N_PUB != nrow(wide))
    warn("cohort is ", nrow(wide), " after MIN_SPOTS but the published n is ",
         N_PUB, " -- the filter does not reproduce the tested cohort")
  if (abs(P_PUB - wt$p.value) > 0.005)
    warn(sprintf("recomputed p %.4f differs from the published %.4f on the same %d patients. ",
                 wt$p.value, P_PUB, nrow(wide)),
         "The cause is how A34's zero difference is ranked: dropping it ",
         "leaves 15 ranks (V = 49), ranking it leaves 16 (V = 54). Both are ",
         "null; the methods should state which was used.")
}

pB <- ggplot(B, aes(grp, cpm)) +
  geom_line(aes(group = patient), colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = grp), size = 1.1) +
  scale_colour_manual(values = c(Nevus = COL_NEV, Melanoma = COL_MEL), guide = "none") +
  annotate("text", x = 1.5, y = max(B$cpm) * 1.08, size = LAB_PT / .pt,
           label = sprintf("n = %d paired lesions\nWilcoxon P = %.3f", N_PUB, P_PUB)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(x = NULL, y = expression(italic("RASSF9")*" (counts per million)")) +
  theme_fig()

## ---- 4. Panel C : bulk against paired ---------------------------------------
say(""); say("== 4. Panel C : designs on one scale ===========================")
CM <- fread(f_cmp, showProgress = FALSE)
CM <- CM[grepl("^RASSF9$", gene, ignore.case = TRUE)]
say("  rows for RASSF9: ", nrow(CM))
CM[, fc   := 2^log2fc]
CM[, flo  := 2^ci_lo]
CM[, fhi  := 2^ci_hi]
CM[, paired := grepl("paired", kind, ignore.case = TRUE) |
               grepl("298774", source)]
## "n = 16/16" on the paired row invites the reading that 32 independent
## samples were compared. They are the SAME 16 patients measured twice, which
## is the entire reason this row sits beside the bulk cohorts, so it says so.
CM[, lab := ifelse(paired,
                   sprintf("%s\npaired, n = %d pairs", source, n_nevus),
                   sprintf("%s\nbulk, n = %d/%d", source, n_nevus, n_melanoma))]
for (i in seq_len(nrow(CM)))
  say(sprintf("    %-12s %-8s fold %.2f [%.2f, %.2f]", CM$source[i],
              ifelse(CM$paired[i], "paired", "bulk"), CM$fc[i], CM$flo[i], CM$fhi[i]))
CM[, lab := factor(lab, levels = rev(lab))]

pC <- ggplot(CM, aes(fc, lab, colour = paired)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_errorbar(aes(xmin = flo, xmax = fhi), orientation = "y",
                width = 0, linewidth = 0.45) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = c(`FALSE` = COL_INK, `TRUE` = COL_PAIR), guide = "none") +
  ## Explicit limits with headroom past the widest interval (0.31 to 2.40).
  ## On the default 5% expansion the GSE298774 upper bound lands within a
  ## hair of the panel border, where a reader cannot tell a wide interval
  ## from a clipped one.
  scale_x_continuous(trans = "log2", breaks = c(0.25, 0.5, 1, 2),
                     labels = c("0.25", "0.5", "1", "2"),
                     limits = c(0.25, 3.1)) +
  ## Three categorical rows inherit the default half-row padding, which on a
  ## tall panel reads as three points adrift in white space.
  scale_y_discrete(expand = expansion(add = 0.62)) +
  labs(x = "Melanoma / nevus (log scale)", y = NULL) +
  theme_fig()

## ---- 5. Assemble -------------------------------------------------------------
say(""); say("== 5. Assembling ===============================================")
## Panel B carries 16 paired lines and needs the vertical room; panel C carries
## three rows and does not. Left at the default they share the height equally,
## which stretches C into mostly empty panel.
FIG <- (pA | ((pB / pC) + plot_layout(heights = c(1.5, 1)))) +
  plot_layout(widths = c(1.15, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = TAG_PT, face = "bold", hjust = 0, vjust = 1))

save_fig <- function(path, ...) {
  tmp <- file.path(dirname(path), paste0(".tmp_", basename(path)))
  ok <- tryCatch({ ggsave(tmp, FIG, ...); TRUE },
                 error = function(e) { warn("could not write ", basename(tmp), ": ",
                                            conditionMessage(e)); FALSE })
  if (!ok) return(invisible(FALSE))
  if (file.exists(path) && !file.remove(path)) {
    warn(basename(path), " is locked (open in a viewer?); new figure left as ", basename(tmp))
    return(invisible(FALSE))
  }
  if (!file.rename(tmp, path)) { warn("rename failed for ", basename(tmp)); return(invisible(FALSE)) }
  say("  wrote ", basename(path)); invisible(TRUE)
}

save_fig(file.path(FIG_DIR, "SupplementaryFigure_Visium.pdf"),
         width = WIDTH_MM, height = HEIGHT_MM, units = "mm", device = cairo_pdf)
save_fig(file.path(FIG_DIR, "SupplementaryFigure_Visium.png"),
         width = WIDTH_MM, height = HEIGHT_MM, units = "mm", dpi = 600)

## ---- 6. Layout report --------------------------------------------------------
## Everything below is TEXT, so the figure can be checked by pasting a log
## instead of uploading an image. It reports the things one would otherwise
## look for: whether anything is clipped, whether labels fit, whether colours
## carry more than one meaning, and -- via a low-resolution ASCII rendering of
## the actual PNG -- what the layout looks like.
say(""); say("== 6. Layout report ============================================")

report <- function() {

  ## -- 6.1 type sizes ---------------------------------------------------------
  say("")
  say("  6.1 Type sizes (pt). Journals require >= 5.5 pt at final size.")
  tp <- c(base = BASE_PT, axis = AXIS_PT, annotation = LAB_PT, tag = TAG_PT)
  for (n in names(tp))
    say(sprintf("    %-12s %4.1f  %s", n, tp[[n]],
                ifelse(tp[[n]] >= 5.5, "ok", "TOO SMALL")))
  if (any(tp < 5.5)) warn("a type size is below the 5.5 pt floor")

  ## -- 6.2 axis range against the data ----------------------------------------
  ## A point or interval end sitting within a whisker of the panel border is
  ## indistinguishable from one that has been clipped off.
  prange <- function(p, which = "x") {
    b <- try(ggplot2::ggplot_build(p), silent = TRUE)
    if (inherits(b, "try-error")) return(c(NA_real_, NA_real_))
    pp <- b$layout$panel_params[[1]]
    r  <- if (which == "x") pp$x.range else pp$y.range
    if (is.null(r)) r <- if (which == "x") pp$x$continuous_range else pp$y$continuous_range
    if (is.null(r)) return(c(NA_real_, NA_real_))
    as.numeric(r)
  }
  ## Margins must be measured in the scale's OWN units -- on panel C the axis
  ## is log2, where a linear margin calculation understates the low end badly
  ## enough to raise a false alarm.
  ## `ends` says which end is allowed to sit on the border. A bar chart is
  ## anchored at zero BY DESIGN, so flagging its baseline as clipped is a
  ## false alarm -- and a report that cries wolf gets ignored, which costs
  ## more than the check is worth.
  chk <- function(panel, axis, lo_d, hi_d, rng, fwd = identity, back = identity,
                  ends = c("both", "hi", "lo")) {
    ends <- match.arg(ends)
    if (any(is.na(rng))) { say(sprintf("    %-2s %s  range unavailable", panel, axis)); return(invisible()) }
    lo_a <- back(rng[1]); hi_a <- back(rng[2])
    span <- rng[2] - rng[1]
    m_lo <- 100 * (fwd(lo_d) - rng[1]) / span
    m_hi <- 100 * (rng[2] - fwd(hi_d)) / span
    mrg  <- switch(ends, both = min(m_lo, m_hi), hi = m_hi, lo = m_lo)
    tight <- mrg < 1.5
    say(sprintf("    %-2s %s  data %8.2f to %8.2f | axis %8.2f to %8.2f | margin %4.1f%% / %4.1f%%  %s%s",
                panel, axis, lo_d, hi_d, lo_a, hi_a, m_lo, m_hi,
                ifelse(tight, "TIGHT", "ok"),
                ifelse(ends == "both", "", sprintf(" (%s end only; the other is anchored)", ends))))
    if (tight)
      warn(sprintf("panel %s %s-axis leaves only %.1f%% margin; an end point may read as clipped",
                   panel, axis, mrg))
  }
  say("")
  say("  6.2 Axis range against the plotted data")
  ## Panel A's bars start at zero, so only the far end can be clipped.
  chk("A", "x", 0, max(A$cpm), prange(pA, "x"), ends = "hi")
  chk("B", "y", min(B$cpm), max(B$cpm), prange(pB, "y"))
  chk("C", "x", min(CM$flo), max(CM$fhi), prange(pC, "x"),
      fwd = log2, back = function(v) 2^v)

  ## -- 6.3 do the labels fit --------------------------------------------------
  ## Panel widths follow from the layout weights; the usable label column is
  ## what is left after the panel itself.
  say("")
  say("  6.3 Longest labels, rendered at final size")
  grDevices::pdf(NULL, width = WIDTH_MM / 25.4, height = HEIGHT_MM / 25.4)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  tmm <- function(s, pt)
    vapply(s, function(x) grid::convertWidth(
      grid::grobWidth(grid::textGrob(x, gp = grid::gpar(fontsize = pt))),
      "mm", valueOnly = TRUE), numeric(1))
  w_left  <- WIDTH_MM * 1.15 / 2.15
  w_right <- WIDTH_MM - w_left
  la <- as.character(A$compartment); wa <- tmm(la, AXIS_PT)
  say(sprintf("    panel A  widest '%s' = %.1f mm, of %.0f mm column  %s",
              la[which.max(wa)], max(wa), w_left,
              ifelse(max(wa) < 0.42 * w_left, "ok", "CROWDED")))
  lc <- unlist(strsplit(as.character(CM$lab), "\n")); wc <- tmm(lc, AXIS_PT)
  say(sprintf("    panel C  widest '%s' = %.1f mm, of %.0f mm column  %s",
              lc[which.max(wc)], max(wc), w_right,
              ifelse(max(wc) < 0.42 * w_right, "ok", "CROWDED")))
  if (max(wa) >= 0.42 * w_left || max(wc) >= 0.42 * w_right)
    warn("a label column is taking more than 42% of its panel's width")

  ## -- 6.4 overlap and crowding ------------------------------------------------
  say("")
  say("  6.4 Overlap and crowding")
  rb <- prange(pB, "y"); spanB <- if (any(is.na(rb))) NA else rb[2] - rb[1]
  ann_y <- max(B$cpm) * 1.08
  say(sprintf("    panel B  annotation sits at %.2f, highest point %.2f, axis top %.2f",
              ann_y, max(B$cpm), if (is.na(spanB)) NA else rb[2]))
  if (!is.na(spanB) && ann_y <= max(B$cpm))
    warn("panel B annotation overlaps the data")
  for (g in levels(B$grp)) {
    v <- sort(B[grp == g, cpm])
    close <- sum(diff(v) < 0.012 * ifelse(is.na(spanB), max(B$cpm), spanB))
    say(sprintf("    panel B  %-9s %2d points, %d adjacent pairs closer than 1.2%% of the axis",
                g, length(v), close))
  }
  say(sprintf("    panel A  %d bars; shortest %.2f CPM is %.1f%% of the longest",
              nrow(A), min(A$cpm), 100 * min(A$cpm) / max(A$cpm)))

  ## Vertical room per categorical row. A forest plot of three rows given the
  ## same height as a sixteen-line paired panel reads as three points adrift
  ## in white space -- a defect invisible to every check above, because the
  ## data, the margins and the labels are all perfectly correct.
  say("")
  say("  6.4b Vertical room per row (panel heights from the layout weights)")
  h_top <- HEIGHT_MM * 1.5 / 2.5      # panel B
  h_bot <- HEIGHT_MM - h_top          # panel C
  rowmm <- function(panel, h_gross, n_rows, usable = 0.78) {
    per <- h_gross * usable / n_rows
    verdict <- if (per > 18) "AIRY" else if (per < 3.5) "CRAMPED" else "ok"
    say(sprintf("    %-2s %5.1f mm gross, %2d rows -> %4.1f mm per row  %s",
                panel, h_gross, n_rows, per, verdict))
    if (verdict != "ok")
      warn(sprintf("panel %s gives %.1f mm per row (%s); adjust plot_layout heights",
                   panel, per, tolower(verdict)))
    invisible()
  }
  rowmm("A", HEIGHT_MM, nrow(A), usable = 0.88)
  rowmm("C", h_bot, nrow(CM))
  say(sprintf("    B %5.1f mm gross for %d paired lines", h_top, nrow(wide)))

  ## -- 6.5 colour roles --------------------------------------------------------
  ## One colour carrying two meanings inside one figure is the commonest way a
  ## multi-panel figure misleads, so the roles are listed rather than assumed.
  say("")
  say("  6.5 Colour roles, by panel")
  say(sprintf("    A fill    %s grey   'Other tissue'", COL_LO))
  say(sprintf("    A fill    %s red    'Melanocytic'", COL_MEL))
  say(sprintf("    B point   %s blue   'Nevus'", COL_NEV))
  say(sprintf("    B point   %s red    'Melanoma'", COL_MEL))
  say(sprintf("    C point   %s black  'bulk design'", COL_INK))
  say(sprintf("    C point   %s purple 'paired design'", COL_PAIR))
  say("    -> red means melanocytic tissue in A and the melanoma group in B,")
  say("       which agree. No colour appears with two unrelated meanings.")

  ## -- 6.6 ASCII rendering of the actual PNG -----------------------------------
  ## The figure is re-saved small, read back, and block-reduced to characters.
  ## Each character takes the DARKEST pixel in its block, so hairlines and
  ## small text survive instead of averaging away into white.
  say("")
  say("  6.6 Low-resolution view of the rendered figure")
  if (!requireNamespace("png", quietly = TRUE))
    try(utils::install.packages("png", repos = "https://cloud.r-project.org"),
        silent = TRUE)
  if (!requireNamespace("png", quietly = TRUE)) {
    say("    the 'png' package is unavailable, so the view is skipped")
    return(invisible())
  }
  tf <- tempfile(fileext = ".png")
  ok <- tryCatch({ ggsave(tf, FIG, width = WIDTH_MM, height = HEIGHT_MM,
                          units = "mm", dpi = 100); TRUE },
                 error = function(e) FALSE)
  if (!ok || !file.exists(tf)) { say("    preview could not be rendered"); return(invisible()) }
  im <- try(png::readPNG(tf), silent = TRUE); unlink(tf)
  if (inherits(im, "try-error")) { say("    preview could not be read"); return(invisible()) }
  ## readPNG returns grey (2D), RGB or RGBA (3D). An RGBA render composited
  ## onto nothing has meaningful colour only where alpha > 0, so the alpha
  ## channel is used to put transparent pixels on white rather than black.
  d <- dim(im)
  g <- if (length(d) == 2L) {
         im
       } else if (length(d) == 3L && d[3] >= 3L) {
         lum <- 0.299 * im[, , 1] + 0.587 * im[, , 2] + 0.114 * im[, , 3]
         if (d[3] >= 4L) lum * im[, , 4] + (1 - im[, , 4]) else lum
       } else if (length(d) == 3L) {
         if (d[3] >= 2L) im[, , 1] * im[, , 2] + (1 - im[, , 2]) else im[, , 1]
       } else {
         NULL
       }
  if (is.null(g) || !is.matrix(g <- as.matrix(g))) {
    say("    the render had an unexpected shape, so the view is skipped")
    return(invisible())
  }
  nr <- nrow(g); nc <- ncol(g)
  cols <- 104L
  rows <- max(12L, round(cols * (HEIGHT_MM / WIDTH_MM) / 2.05))
  ## Explicit integer block bounds. The previous split()/findInterval version
  ## could hand the loop an empty or malformed index and died on subsetting.
  rb <- round(seq(0, nr, length.out = rows + 1L))
  cb <- round(seq(0, nc, length.out = cols + 1L))
  blk <- matrix(1, rows, cols)
  for (i in seq_len(rows)) {
    r1 <- rb[i] + 1L; r2 <- rb[i + 1L]
    if (r2 < r1) next
    sub <- g[r1:r2, , drop = FALSE]
    for (j in seq_len(cols)) {
      c1 <- cb[j] + 1L; c2 <- cb[j + 1L]
      if (c2 >= c1) blk[i, j] <- min(sub[, c1:c2, drop = FALSE])
    }
  }
  ramp <- c("@", "%", "#", "*", "+", "=", "-", ":", ".", " ")
  ## Clamp by assignment, not pmin/pmax. Those take their attributes from the
  ## FIRST argument, so pmin(length(ramp), <matrix>) silently returns a plain
  ## vector and the row loop below then fails on idx[i, ]. Subset-assignment
  ## leaves the dim attribute alone.
  idx <- 1L + floor(blk * length(ramp))
  idx[idx < 1L] <- 1L
  idx[idx > length(ramp)] <- length(ramp)
  if (!is.matrix(idx)) dim(idx) <- dim(blk)
  say(sprintf("    %d x %d characters, darkest-pixel reduction of the %d x %d px render",
              cols, rows, nc, nr))
  say("    |", strrep("-", cols), "|")
  for (i in seq_len(rows)) say("    |", paste0(ramp[idx[i, ]], collapse = ""), "|")
  say("    |", strrep("-", cols), "|")
  ink <- 100 * mean(blk < 0.92)
  say(sprintf("    ink coverage %.1f%% of the canvas", ink))
  invisible()
}
tryCatch(report(), error = function(e)
  say("  [layout report failed: ", conditionMessage(e), " -- the figure itself is fine]"))

say("")
say("  Drawn at ", WIDTH_MM, " mm, the final width. Place at 100%.")
say("")
say("  fig_visium_main_B4 is deliberately NOT part of this figure. Keep it as")
say("  its own supplementary item -- shrinking a spot map with 1 mm scale bars")
say("  into a quarter panel makes both unreadable.")
