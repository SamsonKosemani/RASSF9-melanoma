# ==============================================================================
# 37c_belote_donor_level.R
#
# Converting the 37b headline from pseudo-replicated to donor-supported.
#
# 37b found RASSF9 in 0.08% of melanocytes against 9.16% of keratinocytes and
# 6.32% of eccrine cells, depth-matched, in the same experiment. The pattern is
# large and the depth confound is excluded (1.9-fold range, and keratinocytes
# are the DEEPER group, so thinning worked against the result). But the only
# test behind it was a cell-level Fisher, and 6,379 melanocytes from 22 donors
# are not 6,379 independent observations. 37b declined the donor-level test
# because only three donors had >= 20 cells in both melanocyte and eccrine
# classes -- eccrine glands are dense in acral skin, so those cells are
# probably concentrated in a few donors, and a donor effect could masquerade as
# a cell-type effect.
#
# THE FIX, which should have been the design in the first place: keratinocytes,
# not eccrine, are the right reference. There are 1,790 of them, they come from
# the same epidermis as the melanocytes, and they showed the HIGHER RASSF9 rate
# of the two references anyway. Most donors should contribute both, which makes
# a within-donor paired test possible.
#
# This script asks four things and concludes nothing beyond them:
#   1. How many donors actually contribute each cell class? (the question 37b
#      should have asked before choosing a reference)
#   2. Within each donor that has both, is RASSF9 detected more often in
#      keratinocytes than melanocytes? Paired across donors.
#   3. Is melanocyte detection uniformly ~0, or concentrated in one stage or
#      site? A gene absent everywhere differs from one present in fetal cells
#      only.
#   4. Who are the 5 positive melanocytes? If all in one donor, noise.
#
# CRITERION, FIXED BEFORE THE NUMBERS:
#   DONOR-SUPPORTED = paired test p < 0.05 across donors AND keratinocyte rate
#                     higher in >= 80% of informative donors
#   NOT SUPPORTED   = otherwise; the cell-type claim then rests on cell-level
#                     counts only and must be reported with that caveat
#
# Thinning uses the same seed and target depth as 37b, so the numbers are
# directly comparable rather than a second random draw.
#
# Output: results/geo/belote_donor_level.csv
#         results/geo/belote_melanocyte_strata.csv
#         results/geo/belote_positive_melanocytes.csv
#         figures/geo/geo_rassf9_belote_donor.{pdf,png,svg}
#         results/logs/37c_belote_donor_level.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required")

LOG <- file.path(PATHS$logs, "37c_belote_donor_level.log")
cat("Belote donor-level detection —", format(Sys.time()), "\n\n", file = LOG)
hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

DIR <- file.path(PATHS$raw, "geo", "belote")
RES <- file.path(PATHS$results, "geo")
FIG_GEO <- file.path(PATHS$figures, "geo")
invisible(lapply(c(RES, FIG_GEO), dir.create, recursive = TRUE, showWarnings = FALSE))

TARGETS   <- c("RASSF9", "RASSF10")
MIN_CELLS <- 20
PROP_MIN  <- 0.80
SEED      <- 20261016          # same as 37b
set.seed(SEED)

pick <- function(slot, fb) {
  if (!exists("PAL", inherits = TRUE)) return(fb)
  x <- get("PAL", inherits = TRUE)[[slot]]
  if (!is.null(x) && is.character(x) && nzchar(x[1])) x[1] else fb
}
COL_MEL <- pick("melanoma", "#F2680C"); COL_KRT <- pick("nevus", "#0095D9")
COL_GREY <- pick("grey_pt", "#8E8E8E"); COL_INK <- pick("ink", "#1A1A1A")

# ==============================================================================
# 1. Load and rebuild the matched-depth values
# ==============================================================================
hdr("1. Load")

MET <- data.table::fread(file.path(DIR, "GSE151091_Metadata.csv.gz"),
                         showProgress = FALSE)
MET$group <- dplyr::case_when(
  MET$class_1 %in% c("mel", "cyc_mel")             ~ "Melanocyte",
  MET$class_1 %in% c("krt", "cyc_krt")             ~ "Keratinocyte",
  MET$class_1 %in% c("eccrine", "cyc_eccrine")     ~ "Eccrine",
  MET$class_1 %in% c("dendritic", "cyc_dendritic") ~ "Dendritic",
  MET$class_1 == "T-cell"                          ~ "T cell",
  MET$class_1 == "mast"                            ~ "Mast",
  TRUE                                             ~ as.character(MET$class_1))

DT <- data.table::fread(file.path(DIR, "GSE151091_raw_matrix.csv.gz"),
                        showProgress = FALSE)
genes <- as.character(DT[[1]]); cells <- names(DT)[-1]
keep <- intersect(cells, MET$cell_id)
M <- as.matrix(DT[, keep, with = FALSE]); rownames(M) <- genes
rm(DT); gc(verbose = FALSE)
MET <- MET[match(keep, MET$cell_id), ]
tot <- colSums(M)
log_line("  cells: ", ncol(M), " | donors: ", dplyr::n_distinct(MET$patient),
         logfile = LOG)

cmp <- c("Melanocyte", "Eccrine", "Keratinocyte")
p10 <- tapply(tot, MET$group, function(x) stats::quantile(x, 0.1))
D_TARGET <- floor(min(p10[cmp]))
usable <- tot >= D_TARGET
log_line("  target depth ", format(D_TARGET, big.mark = ","), " | cells retained ",
         sum(usable), " of ", length(tot), logfile = LOG)

thin_gene <- function(g) {
  x <- M[g, ]; v <- rep(NA_real_, length(x))
  i <- which(usable & x > 0)
  if (length(i)) v[i] <- stats::rbinom(length(i), size = x[i],
                                       prob = D_TARGET / tot[i])
  v[usable & x == 0] <- 0
  v
}
TH <- vapply(TARGETS, thin_gene, numeric(ncol(M)))
colnames(TH) <- TARGETS

CELL <- tibble::tibble(
  cell_id = MET$cell_id, donor = MET$patient, group = MET$group,
  dev_stage = if ("dev_stage" %in% names(MET)) MET$dev_stage else NA_character_,
  site = if ("general_location" %in% names(MET)) MET$general_location else NA_character_,
  subtype = if ("subtypes" %in% names(MET)) MET$subtypes else NA_character_,
  total = tot, usable = usable)
for (g in TARGETS) {
  CELL[[paste0(g, "_raw")]]  <- M[g, ]
  CELL[[paste0(g, "_thin")]] <- TH[, g]
}
readr::write_csv(CELL, file.path(RES, "belote_cell_level.csv"))

# ==============================================================================
# 2. Which donors contribute which classes
# ==============================================================================
hdr("2. Donor coverage by class -- what 37b should have checked first")

COV <- CELL |> dplyr::filter(usable) |>
  dplyr::count(group, donor) |>
  dplyr::group_by(group) |>
  dplyr::summarise(donors = dplyr::n_distinct(donor),
                   donors_ge_min = sum(n >= MIN_CELLS),
                   cells = sum(n), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(cells))
for (i in seq_len(nrow(COV)))
  log_line(sprintf("    %-14s %5d cells | %2d donors | %2d donors with >= %d cells",
                   COV$group[i], COV$cells[i], COV$donors[i],
                   COV$donors_ge_min[i], MIN_CELLS), logfile = LOG)

pair_n <- function(a, b) {
  x <- CELL |> dplyr::filter(usable, group %in% c(a, b)) |>
    dplyr::count(donor, group) |> dplyr::filter(n >= MIN_CELLS)
  length(intersect(x$donor[x$group == a], x$donor[x$group == b]))
}
log_line("", logfile = LOG)
log_line("  donors with >= ", MIN_CELLS, " cells in BOTH classes:", logfile = LOG)
log_line("    melanocyte + keratinocyte: ", pair_n("Melanocyte", "Keratinocyte"),
         logfile = LOG)
log_line("    melanocyte + eccrine:      ", pair_n("Melanocyte", "Eccrine"),
         logfile = LOG)
log_line("  -> whichever is larger is the reference the test can actually use",
         logfile = LOG)

# ==============================================================================
# 3. The paired donor-level test
# ==============================================================================
hdr("3. Within-donor comparison")

DON <- CELL |> dplyr::filter(usable) |>
  dplyr::group_by(donor, group) |>
  dplyr::summarise(n = dplyr::n(),
                   det9  = sum(RASSF9_thin > 0, na.rm = TRUE),
                   det10 = sum(RASSF10_thin > 0, na.rm = TRUE),
                   pct9  = 100 * mean(RASSF9_thin > 0, na.rm = TRUE),
                   pct10 = 100 * mean(RASSF10_thin > 0, na.rm = TRUE),
                   .groups = "drop")
readr::write_csv(DON, file.path(RES, "belote_donor_level.csv"))

run_paired <- function(ref, col, lab) {
  a <- as.data.frame(DON)[DON$group == "Melanocyte" & DON$n >= MIN_CELLS,
                          c("donor", "n", col)]
  names(a) <- c("donor", "n_mel", "mel")
  b <- as.data.frame(DON)[DON$group == ref & DON$n >= MIN_CELLS,
                          c("donor", "n", col)]
  names(b) <- c("donor", "n_ref", "ref")
  p <- dplyr::inner_join(a, b, by = "donor")
  log_line("\n  ", lab, " vs Melanocyte, reference = ", ref, logfile = LOG)
  if (!nrow(p)) { log_line("    no donors contribute both", logfile = LOG); return(NULL) }
  log_line(sprintf("    %-12s %8s %8s %10s %10s", "donor", "n mel", "n ref",
                   "mel %", "ref %"), logfile = LOG)
  for (i in seq_len(nrow(p)))
    log_line(sprintf("    %-12s %8d %8d %9.2f%% %9.2f%%", p$donor[i], p$n_mel[i],
                     p$n_ref[i], p$mel[i], p$ref[i]), logfile = LOG)
  higher <- sum(p$ref > p$mel); prop <- higher / nrow(p)
  log_line(sprintf("    reference higher in %d of %d donors (%.0f%%)", higher,
                   nrow(p), 100 * prop), logfile = LOG)
  pv <- NA_real_
  if (nrow(p) >= 4) {
    wt <- suppressWarnings(stats::wilcox.test(p$ref, p$mel, paired = TRUE))
    pv <- wt$p.value
    log_line("    paired Wilcoxon p = ", fmt_p(pv), logfile = LOG)
    st <- stats::binom.test(higher, nrow(p))
    log_line("    sign test p = ", fmt_p(st$p.value), logfile = LOG)
  } else {
    log_line("    fewer than 4 donors; no paired test attempted", logfile = LOG)
  }
  verdict <- if (is.finite(pv) && pv < 0.05 && prop >= PROP_MIN)
    "DONOR-SUPPORTED" else "NOT SUPPORTED at the donor level"
  log_line("    -> ", verdict, logfile = LOG)
  tibble::tibble(gene = lab, reference = ref, donors = nrow(p), higher = higher,
                 prop = prop, p_wilcox = pv, verdict = verdict)
}

VERD <- dplyr::bind_rows(
  run_paired("Keratinocyte", "pct9",  "RASSF9"),
  run_paired("Eccrine",      "pct9",  "RASSF9"),
  run_paired("Keratinocyte", "pct10", "RASSF10"),
  run_paired("Eccrine",      "pct10", "RASSF10"))

# ==============================================================================
# 4. Is melanocyte detection uniformly zero?
# ==============================================================================
hdr("4. Melanocyte detection by stage and site")

log_line("  A gene absent from every melanocyte population is a different",
         logfile = LOG)
log_line("  claim from one present in a single stage or site.", logfile = LOG)

strata <- function(col, lbl) {
  s <- CELL |> dplyr::filter(usable, group == "Melanocyte") |>
    dplyr::group_by(stratum = .data[[col]]) |>
    dplyr::summarise(n = dplyr::n(),
                     det9 = sum(RASSF9_thin > 0, na.rm = TRUE),
                     pct9 = 100 * mean(RASSF9_thin > 0, na.rm = TRUE),
                     det10 = sum(RASSF10_thin > 0, na.rm = TRUE),
                     pct10 = 100 * mean(RASSF10_thin > 0, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::mutate(variable = lbl) |> dplyr::arrange(dplyr::desc(n))
  log_line("\n  by ", lbl, ":", logfile = LOG)
  log_line(sprintf("    %-18s %7s %8s %9s %8s %9s", "stratum", "cells",
                   "R9 det", "R9 %", "R10 det", "R10 %"), logfile = LOG)
  for (i in seq_len(nrow(s)))
    log_line(sprintf("    %-18s %7d %8d %8.2f%% %8d %8.2f%%",
                     as.character(s$stratum[i]), s$n[i], s$det9[i], s$pct9[i],
                     s$det10[i], s$pct10[i]), logfile = LOG)
  s
}
STR <- dplyr::bind_rows(
  if ("dev_stage" %in% names(CELL)) strata("dev_stage", "developmental stage"),
  if ("site" %in% names(CELL)) strata("site", "anatomical site"),
  if ("subtype" %in% names(CELL)) strata("subtype", "melanocyte subtype"))
if (nrow(STR)) readr::write_csv(STR, file.path(RES, "belote_melanocyte_strata.csv"))

# ==============================================================================
# 5. The positive melanocytes
# ==============================================================================
hdr("5. The few positive melanocytes")

POS <- CELL |> dplyr::filter(usable, group == "Melanocyte",
                             RASSF9_thin > 0 | RASSF10_thin > 0) |>
  dplyr::select(cell_id, donor, dev_stage, site, subtype, total,
                RASSF9_raw, RASSF9_thin, RASSF10_raw, RASSF10_thin)
if (nrow(POS)) {
  readr::write_csv(POS, file.path(RES, "belote_positive_melanocytes.csv"))
  log_line("  ", nrow(POS), " melanocytes positive for either gene:", logfile = LOG)
  log_line(sprintf("    %-26s %-12s %-6s %-10s %10s %6s %6s", "cell", "donor",
                   "stage", "site", "lib", "R9", "R10"), logfile = LOG)
  for (i in seq_len(nrow(POS)))
    log_line(sprintf("    %-26s %-12s %-6s %-10s %10s %6.0f %6.0f",
                     POS$cell_id[i], POS$donor[i], as.character(POS$dev_stage[i]),
                     as.character(POS$site[i]),
                     format(round(POS$total[i]), big.mark = ","),
                     POS$RASSF9_raw[i], POS$RASSF10_raw[i]), logfile = LOG)
  nd <- dplyr::n_distinct(POS$donor)
  log_line("  spread across ", nd, " donor(s)", logfile = LOG)
  if (nd == 1) {
    log_line("  -> all in ONE donor: consistent with noise or a donor artefact,",
             logfile = LOG)
    log_line("     not a rare melanocyte population", logfile = LOG)
  } else {
    log_line("  -> spread across donors, so not a single-sample artefact;",
             logfile = LOG)
    log_line("     still far too few to describe a population", logfile = LOG)
  }
} else {
  log_line("  none", logfile = LOG)
}

# ==============================================================================
# 6. Summary
# ==============================================================================
hdr("6. Summary")

if (nrow(VERD)) {
  for (i in seq_len(nrow(VERD)))
    log_line(sprintf("  %-8s vs %-13s %2d donors | ref higher in %2d (%.0f%%) | p = %-10s %s",
                     VERD$gene[i], VERD$reference[i], VERD$donors[i],
                     VERD$higher[i], 100 * VERD$prop[i], fmt_p(VERD$p_wilcox[i]),
                     VERD$verdict[i]), logfile = LOG)
  readr::write_csv(VERD, file.path(RES, "belote_donor_verdict.csv"))
}
log_line("", logfile = LOG)
log_line("  NOTE: section 8 of 37b -- the abundance-matched calibration -- was",
         logfile = LOG)
log_line("  circular. It matched genes on MELANOCYTE expression, which for", logfile = LOG)
log_line("  RASSF9 is zero, so the matched genes were also near zero by", logfile = LOG)
log_line("  construction. It returned 1.60x and should not be cited. The", logfile = LOG)
log_line("  keratinocyte and eccrine references are the non-circular ones.",
         logfile = LOG)

# ==============================================================================
# 7. Figure
# ==============================================================================
hdr("7. Figure")

reset_legends("# Belote donor-level figure legends\n")

pdat <- DON |> dplyr::filter(n >= MIN_CELLS,
                             group %in% c("Melanocyte", "Keratinocyte", "Eccrine"))
pa <- ggplot(pdat, aes(x = group, y = pct9)) +
  geom_line(aes(group = donor), colour = COL_GREY, linewidth = 0.3, alpha = 0.6) +
  geom_point(aes(size = n, colour = group), alpha = 0.85) +
  scale_colour_manual(values = c(Melanocyte = COL_MEL, Keratinocyte = COL_KRT,
                                 Eccrine = "#8C2D04"), guide = "none") +
  scale_size_continuous(range = c(1.2, 3.4), guide = "none") +
  labs(x = NULL, y = "Cells with RASSF9 detected (%)") +
  theme_journal() +
  theme(plot.margin = margin(4, 12, 4, 4, unit = "pt"))

pb <- ggplot(pdat, aes(x = group, y = pct10)) +
  geom_line(aes(group = donor), colour = COL_GREY, linewidth = 0.3, alpha = 0.6) +
  geom_point(aes(size = n, colour = group), alpha = 0.85) +
  scale_colour_manual(values = c(Melanocyte = COL_MEL, Keratinocyte = COL_KRT,
                                 Eccrine = "#8C2D04"), guide = "none") +
  scale_size_continuous(range = c(1.2, 3.4), guide = "none") +
  labs(x = NULL, y = "Cells with RASSF10 detected (%)") +
  theme_journal() +
  theme(plot.margin = margin(4, 12, 4, 4, unit = "pt"))

fig <- patchwork::wrap_plots(list(pa, pb), nrow = 1) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = FIG$base_size + 2, face = "bold",
                                colour = COL_INK))
save_fig_to(fig, "geo_rassf9_belote_donor", FIG_GEO,
            width = FIG$width_double, height = 86)
log_line("  wrote: geo_rassf9_belote_donor", logfile = LOG)

add_legend("geo_rassf9_belote_donor",
  "Figure. Per-donor detection of RASSF9 and RASSF10 by cell class in normal human skin.",
  paste0("(A) RASSF9 and (B) RASSF10 detection rates computed separately within each ",
         "donor, for every donor contributing at least ", MIN_CELLS,
         " cells of a given class, after binomial thinning to a common depth of ",
         format(D_TARGET, big.mark = ","), " counts per cell. Lines join classes ",
         "within a donor, so the comparison is within-subject rather than pooled ",
         "across cells; point size is the number of cells behind each estimate. ",
         "Pooling cells across donors would treat ", sum(pdat$n),
         " cells as independent observations when they derive from ",
         dplyr::n_distinct(pdat$donor), " individuals. Data from Belote et al., ",
         "Nat Cell Biol 2021 (GSE151091)."))

report_warnings(logfile = LOG)
write_session_info("37c_belote_donor_level")
log_line("\nFull log: ", LOG, logfile = NULL)
