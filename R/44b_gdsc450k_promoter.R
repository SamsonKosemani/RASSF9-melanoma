## =============================================================================
## 44b_gdsc450k_promoter.R
##
## STEP 2b of Gustavo's cell-line question.
##
## THE QUESTION, AND WHY IT IS NOT A CORRELATION
##   Scripts 42 and 43 established that RASSF9 sits at the detection floor in
##   melanoma cell lines (4 of 113 at TPM >= 1; 0 of 15 explicitly cutaneous
##   lines). With no expression variance there is nothing to correlate, so this
##   script does NOT attempt a methylation-expression Spearman. Forcing one
##   would produce a number driven entirely by tie structure.
##
##   What survives zero expression variance is a one-sample question:
##
##     Is the RASSF9 promoter core methylated in melanoma cells?
##
##     HIGH -> silencing consistent with the TCGA direction. The correlation is
##             untestable in lines only because the gene is uniformly off.
##     LOW  -> methylation is not what holds RASSF9 off in melanoma cells, and
##             the TCGA association is more likely driven by the non-malignant
##             compartments where the gene is actually expressed (our Tirosh
##             analysis: 26.8% of endothelial cells vs 0.20% of malignant).
##
##   These are different conclusions for the paper, and the second is the more
##   interesting one.
##
## PRE-SPECIFIED CRITERION  (fixed before any number is seen)
##   Per line, the core score = mean beta across the 5 RASSF9 TSS-proximal
##   probes. Verdict on the MEDIAN core score across the 56 melanoma lines:
##     HYPERMETHYLATED        >= 0.60
##     INTERMEDIATE           0.30 - 0.60
##     UNMETHYLATED           <= 0.30
##   Reference point: the TCGA complete-cohort core median is 0.385, i.e.
##   already intermediate in bulk tumour. The comparison that matters is lines
##   versus TCGA per probe, not lines against an absolute threshold alone, so
##   both are reported.
##
## CONTROLS BUILT IN
##   - shore probes (cg19622623, cg18827107) carry the island-shore inversion
##     in TCGA (positive rho) and should behave differently from the core
##   - distal probes as an out-of-promoter reference
##   - detection p-values used to mask unreliable betas rather than ignored
##   - all 1,028 non-melanoma lines retained as a pan-cancer baseline, to show
##     whether any methylation level found is distinctive to melanoma or
##     generic across cancer cell lines
##
## INPUTS
##   data/raw/geo/gdsc450k/GSE68379_Matrix.processed.txt.gz
##   results/ccle/44a_probe_table.csv
##   results/ccle/44a_parse_plan.csv
##   results/ccle/44a_sample_map.csv
##   results/ccle/42_rassf9_melanoma_lines.csv   (optional, for the few
##   results/ccle/42_rassf10_melanoma_lines.csv   expressing lines)
##
## OUTPUTS
##   data/processed/44b_locus_betas.rds            extracted matrix, reusable
##   results/ccle/44b_probe_summary_rassf9.csv
##   results/ccle/44b_probe_summary_rassf10.csv
##   results/ccle/44b_line_scores_rassf9.csv
##   results/ccle/44b_line_scores_rassf10.csv
##   figures/ccle/44b_rassf9_promoter.png
##   figures/ccle/44b_rassf10_promoter.png
##   results/logs/44b_gdsc450k_promoter.log
##
## RUNTIME
##   One streaming pass over ~485,000 rows of a 3.8 Gb gzip, roughly 15 GB
##   decompressed. Expect 10-30 minutes. Only the probe ID is parsed for most
##   rows; full field splitting happens for the 48 wanted rows only.
## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------
if (file.exists("R/00_setup.R")) source("R/00_setup.R")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT     <- if (exists("PROJECT_ROOT")) PROJECT_ROOT else getwd()
GEO_DIR  <- file.path(ROOT, "data", "raw", "geo", "gdsc450k")
RES_DIR  <- file.path(ROOT, "results", "ccle")
FIG_DIR  <- file.path(ROOT, "figures", "ccle")
PROC_DIR <- file.path(ROOT, "data", "processed")
LOG_DIR  <- file.path(ROOT, "results", "logs")
for (d in c(RES_DIR, FIG_DIR, PROC_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_PATH <- file.path(LOG_DIR, "44b_gdsc450k_promoter.log")
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

## Pre-specified thresholds
HYPER_CUT  <- 0.60
UNMETH_CUT <- 0.30
PVAL_CUT   <- 0.01     # mask betas with detection p above this
CHUNK      <- 4000     # lines per read; each line is ~31 KB

## ---- 1. Load the plan from 44a ----------------------------------------------
say("")
say("== 1. Parse plan and probe table ===============================")

f_plan  <- file.path(RES_DIR, "44a_parse_plan.csv")
f_probe <- file.path(RES_DIR, "44a_probe_table.csv")
f_map   <- file.path(RES_DIR, "44a_sample_map.csv")
for (f in c(f_plan, f_probe, f_map))
  if (!file.exists(f)) stop("missing ", basename(f), " -- run 44a first", call. = FALSE)

PLAN   <- fread(f_plan,  showProgress = FALSE)
PROBES <- fread(f_probe, showProgress = FALSE)
MAP    <- fread(f_map,   showProgress = FALSE)

say("  file            : ", PLAN$file)
say("  field offset    : ", PLAN$field_offset)
say("  probe in field  : ", PLAN$probe_data_field)
say("  beta data fields: ", PLAN$beta_data_first, "-", PLAN$beta_data_last)
say("  pval data fields: ", PLAN$pval_data_first, "-", PLAN$pval_data_last)
say("  probes wanted   : ", nrow(PROBES))

if (PLAN$layout != "paired_beta_pval")
  warn("parse plan layout is '", PLAN$layout, "', not paired_beta_pval -- verify before trusting output")

f_met <- file.path(GEO_DIR, PLAN$file)
if (!file.exists(f_met)) stop("matrix not found: ", f_met, call. = FALSE)

## Sample names, in beta-column order, taken from the header the same way 44a did
con <- gzfile(f_met, "rt"); hdr_line <- readLines(con, n = 1, warn = FALSE); close(con)
hdr  <- gsub('^"|"$', "", strsplit(hdr_line, "\t", fixed = TRUE)[[1]])
beta_hdr_idx <- which(grepl("AVG.?Beta", hdr, ignore.case = TRUE))
pval_hdr_idx <- which(grepl("Detection", hdr, ignore.case = TRUE))

clean <- function(x) {
  x <- sub("^CellLine_", "", x)
  x <- sub("_AVG[._ ]?Beta$",       "", x, ignore.case = TRUE)
  x <- sub("_Detection[._ ]?PVal$", "", x, ignore.case = TRUE)
  trimws(x)
}
beta_names <- clean(hdr[beta_hdr_idx])
pval_names <- clean(hdr[pval_hdr_idx])

beta_dat <- beta_hdr_idx + PLAN$field_offset
pval_dat <- pval_hdr_idx + PLAN$field_offset

say("  beta columns    : ", length(beta_dat))
if (!identical(beta_names, pval_names))
  warn("beta and p-value column orders differ -- masking would misalign; ",
       "p-values will be matched by name instead of position")

## ---- 2. Stream the matrix ----------------------------------------------------
say("")
say("== 2. Extracting ", nrow(PROBES), " probe rows ===========================")
say("  one pass over the gzip; only the probe ID is parsed for most rows")

want     <- PROBES$probe
want_set <- new.env(hash = TRUE, parent = emptyenv())
for (w in want) assign(w, TRUE, envir = want_set)

con <- gzfile(f_met, "rt")
invisible(readLines(con, n = 1, warn = FALSE))   # discard header

keep_lines <- character(0)
keep_probe <- character(0)
n_rows <- 0L; n_chunk <- 0L
t0 <- Sys.time()

repeat {
  L <- readLines(con, n = CHUNK, warn = FALSE)
  if (!length(L)) break
  n_rows  <- n_rows + length(L)
  n_chunk <- n_chunk + 1L

  ## The probe ID lives in field 2, within the first ~40 characters (row index
  ## up to 6 digits + tab + a ~10-character cg ID). Truncating first means the
  ## regex never scans the 31 KB of values on each line.
  h  <- substr(L, 1L, 40L)
  pid <- sub("^[^\t]*\t([^\t]+)\t.*$", "\\1", h)

  hits <- which(vapply(pid, function(p) exists(p, envir = want_set, inherits = FALSE),
                       logical(1), USE.NAMES = FALSE))
  if (length(hits)) {
    keep_lines <- c(keep_lines, L[hits])
    keep_probe <- c(keep_probe, pid[hits])
  }

  if (n_chunk %% 25 == 0)
    say(sprintf("    %7d rows | %2d/%d probes found | %.1f min",
                n_rows, length(keep_probe), length(want),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  if (length(keep_probe) == length(want)) {
    say("    all probes found after ", n_rows, " rows -- stopping early")
    break
  }
}
close(con)

say(sprintf("  rows scanned: %d | probes found: %d of %d | %.1f min",
            n_rows, length(keep_probe), length(want),
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

missing <- setdiff(want, keep_probe)
if (length(missing)) {
  warn(length(missing), " probe(s) not found in the matrix: ",
       paste(missing, collapse = ", "))
  say("    these are dropped from the analysis and must be named in the methods")
}
if (!length(keep_probe)) stop("no probes extracted -- the parse plan is wrong", call. = FALSE)

## ---- 3. Build the beta matrix ------------------------------------------------
say("")
say("== 3. Parsing the extracted rows ===============================")

parse_row <- function(line) {
  f <- strsplit(line, "\t", fixed = TRUE)[[1]]
  f <- gsub('^"|"$', "", f)
  list(probe = f[PLAN$probe_data_field],
       beta  = suppressWarnings(as.numeric(f[beta_dat])),
       pval  = suppressWarnings(as.numeric(f[pval_dat])),
       nfield = length(f))
}

PR <- lapply(keep_lines, parse_row)
nf <- vapply(PR, function(x) x$nfield, integer(1))
if (length(unique(nf)) > 1)
  warn("extracted rows have differing field counts (", paste(unique(nf), collapse = ", "),
       ") -- the file may be ragged")
if (unique(nf)[1] != PLAN$n_data_fields)
  warn("field count ", unique(nf)[1], " differs from the ", PLAN$n_data_fields,
       " recorded by 44a")

B <- do.call(rbind, lapply(PR, function(x) x$beta))
P <- do.call(rbind, lapply(PR, function(x) x$pval))
rownames(B) <- rownames(P) <- vapply(PR, function(x) x$probe, character(1))
colnames(B) <- beta_names
colnames(P) <- pval_names

say("  matrix: ", nrow(B), " probes x ", ncol(B), " cell lines")
say(sprintf("  beta range across all extracted values: %.4f to %.4f",
            min(B, na.rm = TRUE), max(B, na.rm = TRUE)))
if (max(B, na.rm = TRUE) > 1.01) warn("beta values exceed 1 -- column mapping is wrong, stop here")

## Mask unreliable betas by detection p-value
if (!identical(colnames(B), colnames(P))) P <- P[, colnames(B), drop = FALSE]
n_bad <- sum(P > PVAL_CUT, na.rm = TRUE)
say(sprintf("  betas masked for detection p > %.2f: %d of %d (%.3f%%)",
            PVAL_CUT, n_bad, length(B), 100 * n_bad / length(B)))
B[P > PVAL_CUT] <- NA_real_

saveRDS(list(beta = B, pval = P, probes = PROBES, plan = PLAN),
        file.path(PROC_DIR, "44b_locus_betas.rds"))
say("  wrote data/processed/44b_locus_betas.rds (re-usable; no second 3.8 GB pass needed)")

## ---- 4. Split melanoma vs the rest -------------------------------------------
say("")
say("== 4. Cell-line groups =========================================")

MEL <- MAP[!is.na(field) & uveal == FALSE]
mel_names <- intersect(unique(MEL$geo_name), colnames(B))
oth_names <- setdiff(colnames(B), unique(MAP$geo_name[!is.na(MAP$field)]))

say("  melanoma lines (non-uveal) with methylation: ", length(mel_names))
say("  other cancer lines (pan-cancer baseline)   : ", length(oth_names))
if (length(mel_names) < 20)
  warn("only ", length(mel_names), " melanoma lines -- descriptive only")

sub_of <- setNames(MEL$subtype, MEL$geo_name)
tb <- table(sub_of[mel_names])
say("")
say("  by subtype:")
for (i in seq_along(tb)) say(sprintf("    %-40s %3d", names(tb)[i], tb[i]))

## ---- 5. Per-probe comparison against TCGA ------------------------------------
say("")
say("== 5. Per-probe: melanoma cell lines vs TCGA tumours ===========")

probe_summary <- function(g) {
  pp <- PROBES[gene == g][probe %in% rownames(B)]
  pp <- pp[order(dist_tss)]
  rbindlist(lapply(seq_len(nrow(pp)), function(i) {
    pr <- pp$probe[i]
    v  <- B[pr, mel_names]; v <- v[is.finite(v)]
    o  <- B[pr, oth_names]; o <- o[is.finite(o)]
    data.table(gene = g, probe = pr, dist_tss = pp$dist_tss[i], set = pp$set[i],
               n_lines = length(v),
               line_median = if (length(v)) median(v) else NA_real_,
               line_iqr    = if (length(v)) IQR(v)    else NA_real_,
               line_min    = if (length(v)) min(v)    else NA_real_,
               line_max    = if (length(v)) max(v)    else NA_real_,
               other_median = if (length(o)) median(o) else NA_real_,
               tcga_median = pp$tcga_beta[i],
               tcga_rho    = pp$tcga_rho[i],
               delta       = if (length(v)) median(v) - pp$tcga_beta[i] else NA_real_)
  }))
}

for (g in c("RASSF9", "RASSF10")) {
  S <- probe_summary(g)
  say("")
  say("  ", g, ":")
  say(sprintf("    %-12s %8s %-13s %5s %8s %7s %9s %9s %8s",
              "probe", "TSSdist", "set", "n", "lineMed", "IQR", "panCancer", "TCGAmed", "delta"))
  for (i in seq_len(nrow(S))) {
    s <- S[i]
    say(sprintf("    %-12s %+8d %-13s %5d %8.3f %7.3f %9.3f %9.3f %+8.3f",
                s$probe, s$dist_tss, s$set, s$n_lines,
                s$line_median, s$line_iqr, s$other_median, s$tcga_median, s$delta))
  }
  for (st in c("tss_proximal", "shore", "distal")) {
    ss <- S[set == st]
    if (nrow(ss))
      say(sprintf("    %-13s mean of probe medians: lines %.3f | pan-cancer %.3f | TCGA %.3f",
                  st, mean(ss$line_median, na.rm = TRUE),
                  mean(ss$other_median, na.rm = TRUE),
                  mean(ss$tcga_median, na.rm = TRUE)))
  }
  fwrite(S, file.path(RES_DIR, sprintf("44b_probe_summary_%s.csv", tolower(g))))
}

## ---- 6. Per-line core score and the verdict ----------------------------------
say("")
say("== 6. Per-line promoter-core score =============================")

line_scores <- function(g) {
  core  <- PROBES[gene == g & set == "tss_proximal" & probe %in% rownames(B), probe]
  shore <- PROBES[gene == g & set == "shore"        & probe %in% rownames(B), probe]
  dist  <- PROBES[gene == g & set == "distal"       & probe %in% rownames(B), probe]
  sc <- function(pr, cols) if (!length(pr)) rep(NA_real_, length(cols)) else
    colMeans(B[pr, cols, drop = FALSE], na.rm = TRUE)
  data.table(
    gene = g,
    geo_name = mel_names,
    subtype  = unname(sub_of[mel_names]),
    core  = sc(core,  mel_names),
    shore = sc(shore, mel_names),
    distal= sc(dist,  mel_names),
    n_core_probes = length(core)
  )
}

verdicts <- character(0)
for (g in c("RASSF9", "RASSF10")) {
  LS <- line_scores(g)
  v  <- LS$core[is.finite(LS$core)]

  med <- median(v)
  verdict <- if (med >= HYPER_CUT) "HYPERMETHYLATED" else
             if (med <= UNMETH_CUT) "UNMETHYLATED"   else "INTERMEDIATE"
  verdicts[g] <- verdict

  tcga_core <- mean(PROBES[gene == g & set == "tss_proximal", tcga_beta])

  say("")
  say("  ", g, " promoter core across ", length(v), " melanoma lines:")
  say(sprintf("    median %.3f | IQR %.3f | range %.3f - %.3f",
              med, IQR(v), min(v), max(v)))
  say(sprintf("    lines >= %.2f (hypermethylated): %d (%.1f%%)",
              HYPER_CUT, sum(v >= HYPER_CUT), 100 * mean(v >= HYPER_CUT)))
  say(sprintf("    lines <= %.2f (unmethylated)   : %d (%.1f%%)",
              UNMETH_CUT, sum(v <= UNMETH_CUT), 100 * mean(v <= UNMETH_CUT)))
  say(sprintf("    TCGA bulk-tumour core mean for reference: %.3f", tcga_core))
  say("    VERDICT: ", verdict)

  ## shore and distal as internal references
  say(sprintf("    shore  median %.3f | distal median %.3f",
              median(LS$shore, na.rm = TRUE), median(LS$distal, na.rm = TRUE)))

  ## by subtype, reported not pooled
  bys <- LS[is.finite(core), .(n = .N, median_core = median(core)), by = subtype][order(-n)]
  for (i in seq_len(nrow(bys)))
    say(sprintf("    %-40s n = %2d | median core %.3f",
                bys$subtype[i], bys$n[i], bys$median_core[i]))

  fwrite(LS[order(-core)], file.path(RES_DIR, sprintf("44b_line_scores_%s.csv", tolower(g))))
}

## ---- 7. The few expressing lines ---------------------------------------------
say("")
say("== 7. The handful of RASSF9-expressing lines ===================")

f_e9 <- file.path(RES_DIR, "42_rassf9_melanoma_lines.csv")
if (file.exists(f_e9)) {
  E <- fread(f_e9, showProgress = FALSE)
  EM <- merge(MAP[!is.na(field), .(model_id, geo_name)], E, by = "model_id")
  EM <- EM[geo_name %in% mel_names]
  core9 <- PROBES[gene == "RASSF9" & set == "tss_proximal" & probe %in% rownames(B), probe]
  EM[, core := colMeans(B[core9, geo_name, drop = FALSE], na.rm = TRUE)]

  pos <- EM[RASSF9 >= 1]
  say("  lines with BOTH expression and methylation: ", nrow(EM))
  say("  of those at TPM >= 1: ", nrow(pos))
  if (nrow(pos)) {
    say("")
    say(sprintf("    %-14s %10s %8s", "line", "log2TPM", "core"))
    for (i in seq_len(nrow(pos)))
      say(sprintf("    %-14s %10.3f %8.3f", pos$geo_name[i], pos$RASSF9[i], pos$core[i]))
    say(sprintf("    median core, expressing (n=%d) %.3f | non-expressing (n=%d) %.3f",
                nrow(pos), median(pos$core, na.rm = TRUE),
                nrow(EM) - nrow(pos), median(EM[RASSF9 < 1]$core, na.rm = TRUE)))
  }
  say("")
  say("  NO TEST IS PERFORMED HERE. With this few expressing lines the")
  say("  comparison is anecdote, reported for completeness only. It must not")
  say("  be described as evidence in either direction.")
} else {
  say("  42_rassf9_melanoma_lines.csv not found -- skipped")
}

## ---- 8. Figures ---------------------------------------------------------------
## RASSF9 and RASSF10 in separate files, per project convention.
say("")
say("== 8. Figures ==================================================")

for (g in c("RASSF9", "RASSF10")) {
  S <- fread(file.path(RES_DIR, sprintf("44b_probe_summary_%s.csv", tolower(g))),
             showProgress = FALSE)
  if (!nrow(S)) next

  pr <- PROBES[gene == g & probe %in% rownames(B)][order(dist_tss)]
  D  <- rbindlist(lapply(seq_len(nrow(pr)), function(i) {
    v <- B[pr$probe[i], mel_names]
    data.table(probe = pr$probe[i], dist_tss = pr$dist_tss[i], set = pr$set[i],
               beta = as.numeric(v))
  }))[is.finite(beta)]

  lv <- pr$probe
  D[,  probe := factor(probe, levels = lv)]
  S[,  probe := factor(probe, levels = lv)]
  D[,  set := factor(set, levels = c("tss_proximal", "shore", "distal"),
                     labels = c("promoter core", "promoter shore", "distal"))]

  p <- ggplot(D, aes(x = probe, y = beta)) +
    geom_jitter(width = 0.17, height = 0, size = 0.9, alpha = 0.45, colour = "grey35") +
    ## ymin/ymax given explicitly: geom_crossbar requires them, and relying on
    ## stat_summary to supply them implicitly would only fail AFTER the long
    ## streaming pass. Also avoids the deprecated `fatten` argument.
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", width = 0.5,
                 linewidth = 0.3, colour = "steelblue4") +
    geom_point(data = S, aes(x = probe, y = tcga_median),
               inherit.aes = FALSE, shape = 23, size = 2.2,
               fill = "firebrick", colour = "white", stroke = 0.4) +
    facet_grid(~ set, scales = "free_x", space = "free_x") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title    = paste0(g, " promoter methylation in melanoma cell lines"),
      subtitle = sprintf(paste0("GSE68379 (Illumina HM450) | %d non-uveal melanoma lines | ",
                                "grey = lines, blue bar = line median, red diamond = TCGA-SKCM median (n = 469)"),
                         length(mel_names)),
      x = NULL, y = "beta value"
    ) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          plot.title    = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 7.4, colour = "grey25"),
          plot.margin   = margin(10, 14, 8, 10))

  out <- file.path(FIG_DIR, sprintf("44b_%s_promoter.png", tolower(g)))
  ggsave(out, p, width = 8.0, height = 4.4, dpi = 400)
  say("  wrote: ", basename(out))
}

## ---- 9. Interpretation guardrails --------------------------------------------
say("")
say("== 9. How to read this =========================================")
say("  RASSF9 verdict: ", verdicts["RASSF9"])
say("")
say("  If INTERMEDIATE or UNMETHYLATED: the RASSF9 promoter is not uniformly")
say("  silenced by methylation in melanoma cells, yet the gene is off. Promoter")
say("  methylation is therefore not sufficient to explain its absence from the")
say("  malignant compartment, and the TCGA correlation is more plausibly driven")
say("  by the compartments that express RASSF9 at all. That is a result, not a")
say("  failure, and it is the honest answer to whether the relationship holds")
say("  'within melanoma cells'.")
say("")
say("  If HYPERMETHYLATED: consistent with the TCGA direction. Note it still")
say("  does not demonstrate a methylation-expression relationship in lines,")
say("  because there is no expression variance to relate it to.")
say("")
say("  LIMITATIONS to carry into the write-up:")
say("   1. Cell lines are not tumours. Cultured lines drift epigenetically,")
say("      and global hypermethylation is common in long-passage lines.")
say("   2. 56 of 154 DepMap melanoma models matched a GDSC sample by name.")
say("      Most of the rest are absent from the GDSC panel, but some may be")
say("      name-spelling losses.")
say("   3. Methylation (GDSC/Sanger) and expression (DepMap) come from")
say("      different laboratories and different cell stocks.")
say("   4. No correlation is reported, by design, because expression has no")
say("      usable variance in these lines.")
