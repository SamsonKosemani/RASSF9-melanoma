# ==============================================================================
# 26_visium_deconvolve.R
#
# Reproduce the published STdeconvolve topic deconvolution of GSE298774, using
# the exact parameters from Kreuger et al.'s Deconvolute_Visium.md.
#
# WHY, AFTER SCRIPT 25 FAILED. Script 25's per-spot marker score separated
# spots that differed up to eight-fold in sequencing depth, and depth
# residualisation left 13 of 18 sections still unbalanced. The reason is
# structural: counts-per-10,000 corrects the MAGNITUDE of a detected count but
# not the PROBABILITY of detection, and every marker that distinguishes nevus
# from melanoma here (PRAME, CDK2, MKI67) is sparse.
#
# Latent Dirichlet allocation does not have that weakness in the same way. It
# models the counts themselves and returns theta, a COMPOSITION per spot that
# sums to one, so a shallow and a deep spot from the same tissue get similar
# compositions. And topics are identified from beta -- each topic's
# distribution over genes -- which is not a per-spot measurement at all and so
# cannot be driven by per-spot depth.
#
# This is also not a new method of my own invention: it is the authors' own
# published pipeline at their own published K, which is why adopting it is not
# the retuning that script 25's stop rule forbade.
#
# TWO DELIBERATE DEVIATIONS, recorded because they are real.
#   1. They ran LDA on SCTransform-corrected counts (se@assays$SCT$counts).
#      This uses RAW INTEGER counts, which is what STdeconvolve documents as
#      its input and what cleanCounts expects. SCT counts are a corrected,
#      non-integer-derived quantity; feeding the raw matrix is the more
#      faithful use of the tool, at the cost of not matching them exactly.
#   2. They fitted Ks = 10:20 and then used models[["15"]]. This fits K = 15
#      only. Sweeping K and choosing our own would be a new analytic choice
#      requiring its own justification; taking their K is inheriting theirs.
#
# ALL SECTIONS ARE FITTED JOINTLY, as they did. That is what makes topics
# comparable across sections -- a per-section LDA would give 18 unrelated sets
# of topics and no way to pool patients.
#
# RUNTIME. LDA on tens of thousands of spots takes a long time: expect
# 30 minutes to a few hours. The fitted model is saved, so this script is run
# ONCE and script 27 reads its output.
#
# NEITHER RASSF9 NOR RASSF10 APPEARS IN THIS SCRIPT. The corpus is restricted
# by STdeconvolve's own overdispersion criteria, which are blind to our genes,
# and nothing here selects or inspects them.
#
# Output: data/processed/gse298774_lda_K15.rds   (the fitted model + theta/beta)
#         data/processed/gse298774_corpus.rds    (the corpus actually fitted)
#         results/geo/gse298774_lda_corpus_stats.csv
#         results/logs/26_visium_deconvolve.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(Matrix))

LOG <- file.path(PATHS$logs, "26_visium_deconvolve.log")
cat("GSE298774 deconvolution —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

EXT_DIR <- file.path(PATHS$raw, "geo", "GSE298774", "extracted")
PROC    <- PATHS$processed
RES_GEO <- file.path(PATHS$results, "geo")
invisible(lapply(c(PROC, RES_GEO), dir.create, recursive = TRUE,
                 showWarnings = FALSE))

set.seed(123)

# ---- Their parameters, verbatim ---------------------------------------------
MIN_LIB_SIZE  <- 100    # cleanCounts
MIN_READS     <- 1
MIN_DETECTED  <- 1
REMOVE_ABOVE  <- 1.0    # restrictCorpus
REMOVE_BELOW  <- 0.05
ALPHA         <- 0.05
PERC_RARE     <- 0.05   # fitLDA
K             <- 15     # the model they used
PERC_FILT     <- 0.05   # getBetaTheta
BETA_SCALE    <- 1
DROP_SECTIONS <- "A2"   # no barcodes or features in the deposit

# ==============================================================================
# 1. Packages
# ==============================================================================
hdr("1. Packages")

# STdeconvolve is NOT in Bioconductor 3.23 -- the first attempt failed with
# "package 'STdeconvolve' is not available for Bioconductor version '3.23'".
# The live source is the authors' GitHub repository, currently version 1.3.2,
# and it contains no src/ directory, so it is pure R and needs no compiler.
# Its dependencies split across CRAN and Bioconductor and must go in first,
# because install_github will not fetch Bioconductor packages by itself.
CRAN_DEPS <- c("topicmodels", "scatterpie", "viridis", "slam", "clue",
               "reshape2", "mgcv", "Matrix", "remotes")
BIOC_DEPS <- c("BiocParallel", "fgsea")

need_cran <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    log_line("  installing (CRAN): ", p, logfile = LOG)
    utils::install.packages(p, repos = "https://cloud.r-project.org")
  }
}
need_bioc <- function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    log_line("  installing (Bioconductor): ", p, logfile = LOG)
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

invisible(lapply(CRAN_DEPS, need_cran))
invisible(lapply(BIOC_DEPS, need_bioc))

miss <- c(CRAN_DEPS, BIOC_DEPS)[
  !vapply(c(CRAN_DEPS, BIOC_DEPS), requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  log_line("  [FLAG] dependencies still missing: ",
           paste(miss, collapse = ", "), logfile = LOG)
}

if (!requireNamespace("STdeconvolve", quietly = TRUE)) {
  log_line("  installing STdeconvolve from GitHub (JEFworks-Lab)", logfile = LOG)
  ok <- try(remotes::install_github("JEFworks-Lab/STdeconvolve",
                                    upgrade = "never",
                                    build_vignettes = FALSE,
                                    dependencies = FALSE), silent = TRUE)
  if (inherits(ok, "try-error")) {
    log_line("  [FLAG] GitHub install failed: ", as.character(ok), logfile = LOG)
  }
}

if (!requireNamespace("STdeconvolve", quietly = TRUE)) {
  log_line("", logfile = LOG)
  log_line("  STdeconvolve could not be installed. It is not in Bioconductor",
           logfile = LOG)
  log_line("  3.23 and the GitHub install did not succeed. Options, in order:",
           logfile = LOG)
  log_line("    1. install it by hand in R:", logfile = LOG)
  log_line('       remotes::install_github("JEFworks-Lab/STdeconvolve")',
           logfile = LOG)
  log_line("       and report whatever error it prints -- the message matters",
           logfile = LOG)
  log_line("       more than the failure", logfile = LOG)
  log_line("    2. if the blocker is a dependency, install that first and",
           logfile = LOG)
  log_line("       retry; the package itself is pure R and needs no compiler",
           logfile = LOG)
  log_line("    3. as a last resort the LDA can be fitted with topicmodels",
           logfile = LOG)
  log_line("       directly, but cleanCounts and restrictCorpus would then",
           logfile = LOG)
  log_line("       have to be reimplemented, which is a deviation from their",
           logfile = LOG)
  log_line("       published pipeline and would need stating as one.",
           logfile = LOG)
  stop("STdeconvolve is required for this script; see the log.")
}
log_line("  STdeconvolve ", as.character(utils::packageVersion("STdeconvolve")),
         logfile = LOG)
log_line("  topicmodels  ", as.character(utils::packageVersion("topicmodels")),
         logfile = LOG)

# ==============================================================================
# 2. Build one combined counts matrix across sections
# ==============================================================================
hdr("2. Combined matrix")

files <- list.files(EXT_DIR, full.names = TRUE)
if (!length(files)) stop("No extracted files in ", EXT_DIR, " -- run script 23.")

sec_of <- function(b) sub("^GSM\\d+_(.+?)_(aligned_fiducials|barcodes|features|matrix|scalefactors_json|tissue_hires_image|tissue_lowres_image|tissue_positions).*$", "\\1", b)
inv <- tibble::tibble(path = files, base = basename(files)) |>
  dplyr::mutate(section = sec_of(base),
                role = dplyr::case_when(
                  grepl("_matrix\\.mtx", base)     ~ "matrix",
                  grepl("_barcodes\\.tsv", base)   ~ "barcodes",
                  grepl("_features\\.tsv", base)   ~ "features",
                  grepl("_tissue_positions", base) ~ "positions",
                  TRUE ~ "other")) |>
  dplyr::filter(role != "other")
grab <- function(sec, role) {
  p <- inv$path[inv$section == sec & inv$role == role]
  if (!length(p)) NA_character_ else p[1]
}
SECTIONS <- sort(setdiff(unique(inv$section), DROP_SECTIONS))
SECTIONS <- SECTIONS[vapply(SECTIONS, function(s)
  !any(is.na(c(grab(s, "matrix"), grab(s, "barcodes"), grab(s, "features")))),
  logical(1))]
log_line("  sections: ", length(SECTIONS), logfile = LOG)

read_pos <- function(p) {
  if (is.na(p)) return(NULL)
  first <- readr::read_lines(p, n_max = 1)
  has_hdr <- grepl("barcode", first, ignore.case = TRUE)
  tb <- readr::read_csv(p, col_names = has_hdr, show_col_types = FALSE,
                        progress = FALSE)
  if (!has_hdr) {
    nm <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row", "pxl_col")
    names(tb)[seq_len(min(6, ncol(tb)))] <- nm[seq_len(min(6, ncol(tb)))]
  }
  bcol <- grep("barcode", names(tb), ignore.case = TRUE, value = TRUE)[1]
  icol <- grep("in_?tissue", names(tb), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(bcol) || is.na(icol)) return(NULL)
  tibble::tibble(barcode = as.character(tb[[bcol]]),
                 in_tissue = as.integer(tb[[icol]]))
}

mats <- list(); meta <- list()
for (sec in SECTIONS) {
  m <- Matrix::readMM(grab(sec, "matrix"))
  bc <- readr::read_tsv(grab(sec, "barcodes"), col_names = FALSE,
                        show_col_types = FALSE, progress = FALSE)[[1]]
  ft <- readr::read_tsv(grab(sec, "features"), col_names = FALSE,
                        show_col_types = FALSE, progress = FALSE)
  if (ncol(m) == nrow(ft) && nrow(m) == length(bc)) m <- Matrix::t(m)
  if (nrow(m) != nrow(ft) || ncol(m) != length(bc)) {
    log_line("  [FLAG] ", sec, ": dimension mismatch; skipped", logfile = LOG)
    next
  }
  sym <- toupper(as.character(ft[[min(2, ncol(ft))]]))

  pos <- read_pos(grab(sec, "positions"))
  keep <- rep(TRUE, length(bc))
  if (!is.null(pos)) {
    it <- pos$in_tissue[match(bc, pos$barcode)]
    keep <- !is.na(it) & it == 1L
  }
  m <- m[, keep, drop = FALSE]
  bc <- bc[keep]

  # Collapse duplicate symbols by summing, so one gene is one row.
  if (any(duplicated(sym))) {
    grp <- factor(sym, levels = unique(sym))
    agg <- Matrix::sparse.model.matrix(~ grp - 1)
    m <- Matrix::t(agg) %*% m
    rownames(m) <- levels(grp)
  } else {
    rownames(m) <- sym
  }
  colnames(m) <- paste0(sec, "|", bc)

  mats[[sec]] <- m
  meta[[sec]] <- tibble::tibble(spot = colnames(m), section = sec, barcode = bc)
  log_line(sprintf("  %-7s %6d in-tissue spots, %6d genes%s", sec, ncol(m),
                   nrow(m),
                   paste0("  [RASSF9 ",
                          ifelse("RASSF9" %in% rownames(m), "yes", "NO"),
                          ", RASSF10 ",
                          ifelse("RASSF10" %in% rownames(m), "yes", "NO"), "]")),
           logfile = LOG)
  gc(verbose = FALSE)
}
if (!length(mats)) stop("No section could be loaded.")

# ---- Gene sets differ between sections --------------------------------------
# The first run exposed this: A16, A21 and A23 carry 18,082 genes while the
# other 15 sections carry 37,072. Taking the UNION and padding with zeros, as
# that version did, records ~19,000 genes as "zero in every spot" for those
# three sections when the truth is "never measured" -- which is why sparsity
# came out at 6.52%. For the LDA that biases those sections; for a gene-level
# comparison later it would be worse, since a gene absent from the smaller
# annotation would contribute structural zeros and the result would be an
# artefact of annotation rather than biology.
#
# So the INTERSECTION is used. Every gene then has real data in every section.
# It is also probably closer to what was actually measured: script 24 found the
# 37,072-gene list contains genes with exactly zero counts everywhere (GAPDH,
# EEF1A1), because that list is the full annotation while only probe-panel
# genes carry data. The smaller list is plausibly the panel itself.
gene_sets <- lapply(mats, rownames)
sizes <- vapply(gene_sets, length, integer(1))
log_line("", logfile = LOG)
log_line("  distinct annotation sizes across sections: ",
         paste(sort(unique(sizes)), collapse = ", "), logfile = LOG)
if (length(unique(sizes)) > 1) {
  for (s in sort(unique(sizes))) {
    log_line("    ", s, " genes: ",
             paste(names(sizes)[sizes == s], collapse = ", "), logfile = LOG)
  }
}

genes_union <- Reduce(union, gene_sets)
genes <- sort(Reduce(intersect, gene_sets))
log_line("", logfile = LOG)
log_line("  union: ", length(genes_union), " genes | INTERSECTION USED: ",
         length(genes), logfile = LOG)
log_line("  dropped as not measured in every section: ",
         length(genes_union) - length(genes), logfile = LOG)

# Whether the genes of interest survive the intersection is a fact script 28
# depends on, so it is settled here and stated plainly.
for (g in c("RASSF9", "RASSF10")) {
  inter <- g %in% genes
  uni <- g %in% genes_union
  miss_in <- names(gene_sets)[!vapply(gene_sets, function(x) g %in% x,
                                      logical(1))]
  log_line("  ", g, ": in intersection = ", inter,
           " | in union = ", uni,
           if (length(miss_in)) paste0(" | absent from: ",
                                       paste(miss_in, collapse = ", ")) else
             " | present in all sections", logfile = LOG)
}
if (!all(c("RASSF9", "RASSF10") %in% genes)) {
  log_line("", logfile = LOG)
  log_line("  [IMPORTANT] a gene of interest is NOT measured in every section.",
           logfile = LOG)
  log_line("  Script 28 must restrict that gene's comparison to the sections",
           logfile = LOG)
  log_line("  where it is actually on the annotation, and report how many",
           logfile = LOG)
  log_line("  patients that leaves. It must NOT treat absence as zero.",
           logfile = LOG)
}

align <- function(m) m[genes, , drop = FALSE]
cnt <- do.call(cbind, lapply(mats, align))
rm(mats, gene_sets); gc(verbose = FALSE)
META <- dplyr::bind_rows(meta)

log_line("", logfile = LOG)
log_line("  combined: ", nrow(cnt), " genes x ", ncol(cnt), " spots",
         logfile = LOG)
log_line("  sparsity: ", sprintf("%.2f%%", 100 * Matrix::nnzero(cnt) /
                                   (as.numeric(nrow(cnt)) * ncol(cnt))),
         " non-zero (compare 6.52% with the zero-padded union)",
         logfile = LOG)

# ==============================================================================
# 3. cleanCounts and restrictCorpus, their parameters
# ==============================================================================
hdr("3. Corpus")

cleaned <- STdeconvolve::cleanCounts(counts = cnt,
                                     min.lib.size = MIN_LIB_SIZE,
                                     min.reads = MIN_READS,
                                     min.detected = MIN_DETECTED,
                                     verbose = TRUE)
log_line("  after cleanCounts: ", nrow(cleaned), " genes x ", ncol(cleaned),
         " spots", logfile = LOG)

corpus <- STdeconvolve::restrictCorpus(cleaned,
                                       removeAbove = REMOVE_ABOVE,
                                       removeBelow = REMOVE_BELOW,
                                       alpha = ALPHA,
                                       plot = FALSE,
                                       verbose = TRUE)
log_line("  after restrictCorpus: ", nrow(corpus), " genes x ", ncol(corpus),
         " spots", logfile = LOG)

# ---- Spots left empty by the gene restriction -------------------------------
# restrictCorpus keeps only the top overdispersed genes, and some spots carry
# no counts at all among those. fitLDA rejects them outright: "Each row (pixel)
# of `counts` needs to contain at least one non-zero entry". They are dropped
# here, and the number is reported rather than passed over -- a large loss would
# mean the corpus is too narrow for this tissue, which is worth knowing.
empty <- Matrix::colSums(corpus) == 0
log_line("  spots with zero counts across the retained genes: ", sum(empty),
         " (", sprintf("%.2f%%", 100 * mean(empty)), ")", logfile = LOG)
if (any(empty)) {
  lost <- META$section[match(colnames(corpus)[empty], META$spot)]
  tb <- sort(table(lost), decreasing = TRUE)
  log_line("  by section: ",
           paste(sprintf("%s=%d", names(tb), as.integer(tb)),
                 collapse = ", "), logfile = LOG)
  corpus <- corpus[, !empty, drop = FALSE]
  log_line("  corpus after dropping them: ", nrow(corpus), " genes x ",
           ncol(corpus), " spots", logfile = LOG)
}
if (mean(empty) > 0.10) {
  log_line("  [FLAG] more than 10% of spots were emptied by the gene",
           logfile = LOG)
  log_line("  restriction. The retained gene set may be too narrow to",
           logfile = LOG)
  log_line("  represent this tissue; treat the topics with caution.",
           logfile = LOG)
}

# Recorded, not acted on: whether our genes of interest happen to be in the
# corpus is a consequence of STdeconvolve's overdispersion criteria, which know
# nothing about them. It is logged so the fact is on the record either way, and
# it changes nothing downstream -- topic identification in script 27 uses
# lineage markers, and expression is read from the FULL matrix, not the corpus.
in_corpus <- c("RASSF9", "RASSF10") %in% rownames(corpus)
log_line("  (for the record, genes of interest in the fitted corpus: ",
         paste(c("RASSF9", "RASSF10")[in_corpus], collapse = ", "),
         if (!any(in_corpus)) "none" else "", ")", logfile = LOG)

stats <- tibble::tibble(
  stage = c("combined", "cleanCounts", "restrictCorpus"),
  genes = c(nrow(cnt), nrow(cleaned), nrow(corpus)),
  spots = c(ncol(cnt), ncol(cleaned), ncol(corpus)))
readr::write_csv(stats, file.path(RES_GEO, "gse298774_lda_corpus_stats.csv"))

saveRDS(list(corpus = corpus, meta = META),
        file.path(PROC, "gse298774_corpus.rds"))
log_line("  saved: data/processed/gse298774_corpus.rds", logfile = LOG)

# ==============================================================================
# 4. Fit LDA at K = 15
# ==============================================================================
hdr("4. LDA, K = 15")

log_line("  fitting -- this is the slow step, 30 minutes to a few hours",
         logfile = LOG)
log_line("  started: ", format(Sys.time()), logfile = LOG)

t0 <- Sys.time()
ldas <- STdeconvolve::fitLDA(t(as.matrix(corpus)),
                             Ks = K,
                             perc.rare.thresh = PERC_RARE,
                             plot = FALSE,
                             verbose = TRUE)
el <- difftime(Sys.time(), t0, units = "mins")
log_line("  finished in ", round(as.numeric(el), 1), " minutes", logfile = LOG)

# Raw fitLDA output goes to disk FIRST, before anything is done to it. The
# first version of this script lost a 196-minute fit to a syntax error that sat
# between fitLDA returning and the first saveRDS. Nothing expensive is held in
# memory unsaved again.
saveRDS(ldas, file.path(PROC, "gse298774_lda_raw.rds"))
log_line("  saved raw fit: data/processed/gse298774_lda_raw.rds", logfile = LOG)

mod <- if (!is.null(ldas$models[[as.character(K)]])) {
  ldas$models[[as.character(K)]]
} else {
  ldas$models[[1]]
}

res <- STdeconvolve::getBetaTheta(mod, perc.filt = PERC_FILT,
                                  betaScale = BETA_SCALE)
theta <- res$theta   # spots x topics, each row a composition
beta  <- res$beta    # topics x genes

log_line("", logfile = LOG)
log_line("  theta: ", nrow(theta), " spots x ", ncol(theta), " topics",
         logfile = LOG)
log_line("  beta:  ", nrow(beta), " topics x ", ncol(beta), " genes",
         logfile = LOG)

# How much of the tissue each topic accounts for, as a first look.
share <- colMeans(theta)
log_line("", logfile = LOG)
log_line("  mean proportion per topic:", logfile = LOG)
for (i in order(-share)) {
  log_line(sprintf("    topic %-4s %6.3f", colnames(theta)[i], share[i]),
           logfile = LOG)
}

saveRDS(list(model = mod, theta = theta, beta = beta, meta = META,
             params = list(K = K, min_lib_size = MIN_LIB_SIZE,
                           remove_above = REMOVE_ABOVE,
                           remove_below = REMOVE_BELOW, alpha = ALPHA,
                           perc_rare = PERC_RARE, perc_filt = PERC_FILT,
                           input = "raw integer counts, in-tissue spots",
                           sections = SECTIONS)),
        file.path(PROC, "gse298774_lda_K15.rds"))
log_line("", logfile = LOG)
log_line("  saved: data/processed/gse298774_lda_K15.rds", logfile = LOG)

# ==============================================================================
# 5. Next
# ==============================================================================
hdr("5. Next")

log_line("  Script 27 identifies which topics are which from beta -- the", logfile = LOG)
log_line("  topic-by-gene loadings -- using lineage markers, and validates the",
         logfile = LOG)
log_line("  melanocytic topics against the authors' reported CDK2 and PRAME",
         logfile = LOG)
log_line("  differences. Topic NUMBERING here will not match theirs: LDA is",
         logfile = LOG)
log_line("  stochastic, so their X15 (nevus), X3 and X4 (tumour) are not our",
         logfile = LOG)
log_line("  15, 3 and 4. Identification must come from the loadings, never",
         logfile = LOG)
log_line("  from the numbers.", logfile = LOG)
log_line("", logfile = LOG)
log_line("  The depth checks that failed in script 25 are applied again to the",
         logfile = LOG)
log_line("  topic-based labels. LDA should be far more robust to depth, but",
         logfile = LOG)
log_line("  'should be' is not a result, so it is tested.", logfile = LOG)

report_warnings(logfile = LOG)
write_session_info("26_visium_deconvolve")
log_line("\nFull log: ", LOG, logfile = NULL)
