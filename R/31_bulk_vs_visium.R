# ==============================================================================
# 31_bulk_vs_visium.R
#
# Does the paired within-lesion result EXCLUDE the effect the bulk cohorts
# implied, or is it merely too small to adjudicate?
#
# WHY THIS IS NEEDED. Script 30 returned a null: RASSF9 does not differ between
# the nevus and melanoma compartments of the same lesion, fold change 1.16x
# [0.56, 2.40], n = 16 paired patients. A null is only informative if it can
# rule something out. Saying "the bulk fall is compositional" is an inference
# until the bulk effect is placed on the same scale and shown to lie outside
# that interval. This script turns the claim into a number.
#
# THE CRITERION, fixed before the numbers are computed:
#
#   The bulk and Visium estimates are both LOG2 RATIOS of melanoma to nevus, so
#   they are comparable as ratios even though the platforms are not comparable
#   in absolute units -- a ratio is dimensionless and cancels the scaling.
#
#   EXCLUDED: the bulk point estimate lies outside the Visium 95% interval AND
#   the two intervals do not overlap. Then the paired design rules out an
#   effect of the magnitude bulk implied.
#
#   CONSISTENT: the intervals overlap. Then the Visium data are underpowered to
#   adjudicate, and that is what gets reported -- not "no effect".
#
#   This is written down here, before the result, because deciding afterwards
#   what counts as exclusion is how a null gets over-sold.
#
# TWO ESTIMATORS ARE NOT THE SAME THING, and the script says so rather than
# hiding it. The bulk difference is UNPAIRED, between different people. The
# Visium difference is PAIRED, within one person. They answer the same question
# on the same scale but with different precision and different confounding, and
# that asymmetry is the entire point of the comparison rather than a flaw in it.
#
# GSE112509 DESERVES PARTICULAR WEIGHT. Both its arms are laser-microdissected,
# so it is the bulk cohort LEAST exposed to the composition confound and the
# closest bulk analogue to the Visium design. If its effect is excluded by the
# paired interval, that is the strongest form of the statement.
#
# ALSO REPORTED: the stroma-matched bulk estimates from script 20, which are
# already null (RASSF9 d = -0.395, p = 0.339 in GSE112509). Three independent
# routes to the same conclusion is worth stating plainly.
#
# Output: results/geo/bulk_vs_visium_comparison.csv
#         figures/geo/geo_rassf9_bulk_vs_visium.{pdf,png,svg}
#         figures/geo/geo_rassf10_bulk_vs_visium.{pdf,png,svg}
#         results/logs/31_bulk_vs_visium.log
# ==============================================================================

source(here::here("R", "00_setup.R"))

LOG <- file.path(PATHS$logs, "31_bulk_vs_visium.log")
cat("Bulk versus Visium, on one scale —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

RES_GEO <- file.path(PATHS$results, "geo")
FIG_GEO <- file.path(PATHS$figures, "geo")
GENES   <- c("RASSF9", "RASSF10")

pick <- function(slot, fallback) {
  if (!exists("PAL", inherits = TRUE)) return(fallback)
  x <- get("PAL", inherits = TRUE)[[slot]]
  if (!is.null(x) && is.character(x) && nzchar(x[1])) x[1] else fallback
}
COL_BULK <- pick("nevus",    "#0095D9")
COL_VIS  <- pick("melanoma", "#F2680C")
COL_GREY <- pick("grey_pt",  "#8E8E8E")
COL_INK  <- pick("ink",      "#1A1A1A")

# ==============================================================================
# 1. The bulk cohorts
# ==============================================================================
hdr("1. Bulk cohorts")

VF <- file.path(RES_GEO, "nevus_vs_melanoma_values.csv")
if (!file.exists(VF)) stop("Missing ", VF)
V <- readr::read_csv(VF, show_col_types = FALSE, progress = FALSE)
log_line("  samples: ", nrow(V), " across ",
         paste(sort(unique(V$dataset)), collapse = ", "), logfile = LOG)

# These are log2 values, so a difference of means IS a log2 ratio. Checked
# rather than assumed: log2 expression sits in single or low double digits,
# while linear CPM or TPM would not.
rng <- range(c(V$RASSF9, V$RASSF10), na.rm = TRUE)
log_line(sprintf("  value range across both genes: %.2f to %.2f", rng[1], rng[2]),
         logfile = LOG)
if (rng[2] > 100) {
  log_line("  [STOP] these do not look like log2 values. A difference of means",
           logfile = LOG)
  log_line("  would not be a log2 ratio and the comparison below would be",
           logfile = LOG)
  log_line("  meaningless.", logfile = LOG)
  stop("Values are not on a log2 scale; see ", LOG)
}
log_line("  -> consistent with log2 expression; differences are log2 ratios",
         logfile = LOG)

bulk <- list()
for (ds in sort(unique(V$dataset))) {
  for (g in GENES) {
    s <- V |> dplyr::filter(dataset == ds, !is.na(.data[[g]]))
    nv <- s[[g]][s$group == "Nevus"]
    mm <- s[[g]][s$group != "Nevus"]
    if (length(nv) < 3 || length(mm) < 3) next
    tt <- stats::t.test(mm, nv)                    # melanoma minus nevus
    wt <- suppressWarnings(stats::wilcox.test(mm, nv, conf.int = TRUE))
    sp <- sqrt(((length(nv) - 1) * stats::var(nv) +
                (length(mm) - 1) * stats::var(mm)) /
               (length(nv) + length(mm) - 2))
    bulk[[length(bulk) + 1L]] <- tibble::tibble(
      gene = g, source = ds, kind = "bulk, unadjusted",
      n_nevus = length(nv), n_melanoma = length(mm),
      log2fc = unname(diff(rev(tt$estimate))) * -1,
      ci_lo = tt$conf.int[1], ci_hi = tt$conf.int[2],
      p = tt$p.value,
      hl_shift = unname(wt$estimate),
      cohens_d = (mean(mm) - mean(nv)) / sp)
  }
}
B <- dplyr::bind_rows(bulk)
# t.test(mm, nv)$estimate is c(mean of x, mean of y); the difference we want is
# mean(mm) - mean(nv), which is what conf.int is built around. Recompute it
# directly rather than relying on the order of that vector.
B$log2fc <- vapply(seq_len(nrow(B)), function(i) {
  s <- V |> dplyr::filter(dataset == B$source[i])
  mean(s[[B$gene[i]]][s$group != "Nevus"], na.rm = TRUE) -
    mean(s[[B$gene[i]]][s$group == "Nevus"], na.rm = TRUE)
}, numeric(1))

for (i in seq_len(nrow(B))) {
  log_line(sprintf("  %-8s %-10s log2FC %+.3f [%.3f, %.3f] = %.2fx [%.2f, %.2f], p = %s, d = %+.2f",
                   B$gene[i], B$source[i], B$log2fc[i], B$ci_lo[i], B$ci_hi[i],
                   2^B$log2fc[i], 2^B$ci_lo[i], 2^B$ci_hi[i], fmt_p(B$p[i]),
                   B$cohens_d[i]), logfile = LOG)
}

# ==============================================================================
# 2. The stroma-matched bulk estimates
# ==============================================================================
hdr("2. Bulk after matching on stromal composition")

MF <- file.path(RES_GEO, "stromal_matched_subset.csv")
if (file.exists(MF)) {
  M <- readr::read_csv(MF, show_col_types = FALSE, progress = FALSE)
  for (i in seq_len(nrow(M))) {
    log_line(sprintf("  %-8s %-10s d = %+.3f, p = %s, n = %d vs %d, balanced = %s",
                     M$gene[i], M$dataset[i], M$d[i], fmt_p(M$p[i]),
                     M$n_nevus[i], M$n_melanoma[i], M$balanced[i]),
             logfile = LOG)
  }
  log_line("", logfile = LOG)
  log_line("  Reported as Cohen's d, the scale script 20 used. These are not",
           logfile = LOG)
  log_line("  converted to fold changes: the matched subset's sample identities",
           logfile = LOG)
  log_line("  were not saved, so a log2 ratio cannot be recomputed from this",
           logfile = LOG)
  log_line("  summary without re-running the matching. They enter the argument",
           logfile = LOG)
  log_line("  as a THIRD independent route, not as a point on the figure.",
           logfile = LOG)
} else {
  log_line("  stromal_matched_subset.csv not found; skipped", logfile = LOG)
}

# ==============================================================================
# 3. The Visium paired estimate
# ==============================================================================
hdr("3. Visium, paired within lesion")

TF <- file.path(RES_GEO, "gse298774_authorlabels_tests.csv")
if (!file.exists(TF)) stop("Missing ", TF, " -- run R/30 first.")
TT <- readr::read_csv(TF, show_col_types = FALSE, progress = FALSE)
VIS <- TT |>
  dplyr::filter(set == "all labelled spots") |>
  dplyr::transmute(gene, source = "GSE298774", kind = "Visium, paired",
                   n_nevus = n_patients, n_melanoma = n_patients,
                   log2fc, ci_lo, ci_hi, p = t_p,
                   hl_shift = NA_real_, cohens_d = NA_real_)
for (i in seq_len(nrow(VIS))) {
  log_line(sprintf("  %-8s log2FC %+.3f [%.3f, %.3f] = %.2fx [%.2f, %.2f], %d paired patients",
                   VIS$gene[i], VIS$log2fc[i], VIS$ci_lo[i], VIS$ci_hi[i],
                   2^VIS$log2fc[i], 2^VIS$ci_lo[i], 2^VIS$ci_hi[i],
                   VIS$n_nevus[i]), logfile = LOG)
}

ALL <- dplyr::bind_rows(B, VIS)
readr::write_csv(ALL, file.path(RES_GEO, "bulk_vs_visium_comparison.csv"))

# ==============================================================================
# 4. The pre-specified comparison
# ==============================================================================
hdr("4. Does the paired interval exclude the bulk effect?")

log_line("  CRITERION (fixed before the numbers):", logfile = LOG)
log_line("    EXCLUDED   = bulk point estimate outside the Visium interval AND",
         logfile = LOG)
log_line("                 the two intervals do not overlap", logfile = LOG)
log_line("    CONSISTENT = the intervals overlap; the Visium data are",
         logfile = LOG)
log_line("                 underpowered to adjudicate, and that is the claim",
         logfile = LOG)

verdicts <- list()
for (g in GENES) {
  v <- VIS |> dplyr::filter(gene == g)
  if (!nrow(v)) next
  log_line("", logfile = LOG)
  log_line("  ", g, "  --  Visium paired interval ",
           sprintf("%.2fx to %.2fx", 2^v$ci_lo[1], 2^v$ci_hi[1]), logfile = LOG)
  for (i in which(B$gene == g)) {
    pt_in  <- B$log2fc[i] >= v$ci_lo[1] && B$log2fc[i] <= v$ci_hi[1]
    overlap <- B$ci_lo[i] <= v$ci_hi[1] && B$ci_hi[i] >= v$ci_lo[1]
    verdict <- if (!pt_in && !overlap) "EXCLUDED" else "CONSISTENT"
    log_line(sprintf("    %-10s bulk %.2fx [%.2f, %.2f] -> point inside Visium CI: %-5s | intervals overlap: %-5s | %s",
                     B$source[i], 2^B$log2fc[i], 2^B$ci_lo[i], 2^B$ci_hi[i],
                     pt_in, overlap, verdict), logfile = LOG)
    verdicts[[length(verdicts) + 1L]] <- tibble::tibble(
      gene = g, bulk = B$source[i], point_inside = pt_in,
      intervals_overlap = overlap, verdict = verdict)
  }
}
VV <- dplyr::bind_rows(verdicts)

log_line("", logfile = LOG)
log_line("  WHAT EACH VERDICT LICENSES:", logfile = LOG)
log_line("    EXCLUDED   -> the paired within-lesion design rules out a",
         logfile = LOG)
log_line("                  melanocyte-intrinsic effect of the size the bulk",
         logfile = LOG)
log_line("                  cohort reported. The bulk difference must come",
         logfile = LOG)
log_line("                  substantially from something other than the",
         logfile = LOG)
log_line("                  melanocytes, which is the composition explanation.",
         logfile = LOG)
log_line("    CONSISTENT -> the paired design cannot distinguish the two",
         logfile = LOG)
log_line("                  explanations. Report the interval and say so.",
         logfile = LOG)
log_line("                  Do NOT report it as evidence of no effect.",
         logfile = LOG)

if (nrow(VV)) {
  lcm <- VV |> dplyr::filter(bulk == "GSE112509")
  if (nrow(lcm)) {
    log_line("", logfile = LOG)
    log_line("  GSE112509 CARRIES THE MOST WEIGHT: both its arms are laser-",
             logfile = LOG)
    log_line("  microdissected, so it is the bulk cohort least exposed to the",
             logfile = LOG)
    log_line("  composition confound and the closest analogue to the Visium",
             logfile = LOG)
    log_line("  design. Its verdicts:", logfile = LOG)
    for (i in seq_len(nrow(lcm))) {
      log_line("    ", lcm$gene[i], ": ", lcm$verdict[i], logfile = LOG)
    }
  }
}

# ==============================================================================
# 5. Figure: one forest plot per gene
# ==============================================================================
hdr("5. Figure")

reset_legends("# Bulk versus Visium figure legends\n")

for (g in GENES) {
  f <- ALL |>
    dplyr::filter(gene == g, !is.na(log2fc)) |>
    dplyr::mutate(
      # SHORT labels. The first version spelled out "bulk, unadjusted
      # (n = 23 vs 57)", which consumed nearly half the canvas, pushed the
      # plotting panel to the right and left the centred axis title running off
      # the edge. Shortening the title alone could not fix that, because the
      # cause was the width of these labels, not the length of the title.
      label = paste0(source, "\n",
                     ifelse(kind == "Visium, paired", "paired, n = ",
                            "bulk, n = "),
                     ifelse(kind == "Visium, paired",
                            as.character(n_nevus),
                            paste0(n_nevus, "/", n_melanoma))),
      fc = 2^log2fc, lo = 2^ci_lo, hi = 2^ci_hi,
      is_vis = kind == "Visium, paired")
  if (!nrow(f)) next
  f$label <- factor(f$label, levels = rev(f$label))

  p <- ggplot(f, aes(x = fc, y = label, colour = is_vis)) +
    geom_vline(xintercept = 1, linetype = "22", colour = COL_GREY,
               linewidth = 0.4) +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                  width = 0.18, linewidth = 0.5) +
    geom_point(size = 2.4) +
    scale_colour_manual(values = c(`TRUE` = COL_VIS, `FALSE` = COL_BULK),
                        guide = "none") +
    scale_x_continuous(trans = "log2",
                       breaks = c(0.125, 0.25, 0.5, 1, 2, 4),
                       labels = c("0.125", "0.25", "0.5", "1", "2", "4")) +
    coord_cartesian(clip = "off") +
    # The gene is already named in the figure legend, so the axis does not have
    # to carry it. The first version's label ran off the right edge of a
    # single-column canvas and "(log scale)" was clipped.
    labs(x = bquote(bold("Melanoma / nevus (log scale)")),
         y = NULL) +
    theme_journal(grid = FALSE) +
    theme(panel.grid.major.x = element_line(linewidth = 0.25,
                                            colour = "grey92"),
          axis.text.y = element_text(size = FIG$base_size - 1, lineheight = 1.1),
          plot.margin = margin(4, 12, 4, 4, unit = "pt"))

  stem <- paste0("geo_", tolower(g), "_bulk_vs_visium")
  save_fig_to(p, stem, FIG_GEO, width = FIG$width_double, height = 72)

  vv <- VV |> dplyr::filter(gene == g)
  add_legend(stem,
    paste0("Figure. ", g,
           " between nevus and melanoma, by study design."),
    paste0("Fold change of ", g, " in melanoma relative to nevus, on a log ",
           "scale, with 95% confidence intervals. The bulk cohorts compare ",
           "different people and cannot separate a melanocytic change from a ",
           "change in tissue composition. GSE298774 compares the nevus and ",
           "melanoma compartments of the SAME lesion, paired within patient, ",
           "so both share one stroma and one patient. ",
           if (nrow(vv)) paste0("Bulk estimates: ",
                                paste(paste0(vv$bulk, " ", tolower(vv$verdict)),
                                      collapse = "; "),
                                " by the paired interval. ") else "",
           "GSE112509 is laser-microdissected in both arms and is therefore ",
           "the bulk cohort least exposed to the composition confound."))
  log_line("  wrote: ", stem, logfile = LOG)
}

# ==============================================================================
# 6. The three routes
# ==============================================================================
hdr("6. Summary of evidence")

log_line("  Three independent routes bear on whether the nevus-to-melanoma",
         logfile = LOG)
log_line("  fall in these genes is melanocyte-intrinsic:", logfile = LOG)
log_line("    1. BULK, UNADJUSTED -- large and significant in both cohorts, but",
         logfile = LOG)
log_line("       cannot separate melanocytes from composition", logfile = LOG)
log_line("    2. BULK, STROMA-MATCHED (script 20) -- DISCORDANT between cohorts,",
         logfile = LOG)
log_line("       and this is not glossed: in GSE112509, the laser-", logfile = LOG)
log_line("       microdissected cohort, the RASSF9 difference collapses",
         logfile = LOG)
log_line("       (d = -0.395, p = 0.339); in GSE98394 it does not",
         logfile = LOG)
log_line("       (d = -2.750, p = 1.65e-04). But GSE98394's rests on 9 matched",
         logfile = LOG)
log_line("       pairs from a cohort whose pre-matching SMDs were 1.5 to 2.3,",
         logfile = LOG)
log_line("       so the matching extrapolates rather than balances, and",
         logfile = LOG)
log_line("       NEITHER matched analysis reached balance. The better-balanced",
         logfile = LOG)
log_line("       cohort is the null one.", logfile = LOG)
log_line("    3. VISIUM, PAIRED WITHIN LESION (script 30) -- null, in a design",
         logfile = LOG)
log_line("       where composition cannot be the explanation, but with an",
         logfile = LOG)
log_line("       interval too wide to exclude route 1's magnitude", logfile = LOG)
log_line("", logfile = LOG)
log_line("  WHERE THAT LEAVES THE QUESTION. The strongest evidence that the",
         logfile = LOG)
log_line("  bulk fall is compositional is NOT the paired null, which is",
         logfile = LOG)
log_line("  underpowered. It is the lineage result: only 17.1% of RASSF9",
         logfile = LOG)
log_line("  molecules sit in melanocytic spots at all (script 30 section 4),",
         logfile = LOG)
log_line("  replicating GSE230574 in a second technology. A gene that is",
         logfile = LOG)
log_line("  mostly stromal in skin will move with stromal content in bulk",
         logfile = LOG)
log_line("  tissue whether or not its melanocytic expression changes. That is",
         logfile = LOG)
log_line("  a large, directly measured effect rather than an inference from a",
         logfile = LOG)
log_line("  non-significant comparison, and it should carry the argument.",
         logfile = LOG)

report_warnings(logfile = LOG)
write_session_info("31_bulk_vs_visium")
log_line("\nFull log: ", LOG, logfile = NULL)
