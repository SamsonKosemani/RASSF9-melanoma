# ==============================================================================
# 26b_visium_lda_fit.R
#
# RECOVERY SCRIPT. Fits the LDA from the corpus script 26 already saved, and
# saves the fitted object immediately.
#
# WHY THIS EXISTS. Script 26's first run built the corpus, saved it, fitted the
# LDA in 195.9 minutes -- and then died on a syntax error in the two lines that
# followed fitLDA, before the first saveRDS. The fit was lost. The corpus was
# not: data/processed/gse298774_corpus.rds is on disk and is exactly the matrix
# that was fitted.
#
# So this script does only the expensive part, and does it defensively:
#   - loads the saved corpus rather than re-reading 18 matrices
#   - writes the RAW fitLDA object to disk the instant it returns, before any
#     downstream processing can fail
#   - then extracts theta and beta, in a second, cheap step that can be rerun
#     from the raw object at no cost if anything goes wrong
#
# NOTHING ABOUT THE MODEL CHANGES. Same corpus, same K, same parameters, same
# seed. This is a resumption, not a reanalysis.
#
# If gse298774_lda_raw.rds already exists, the fit is SKIPPED and the script
# goes straight to extraction. So it is safe to rerun.
#
# NEITHER RASSF9 NOR RASSF10 APPEARS IN THIS SCRIPT.
#
# Output: data/processed/gse298774_lda_raw.rds   (raw fitLDA output)
#         data/processed/gse298774_lda_K15.rds   (model + theta + beta + meta)
#         results/logs/26b_visium_lda_fit.log
#
# RUNTIME: about 200 minutes for the fit, seconds for extraction.
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages(library(Matrix))

LOG <- file.path(PATHS$logs, "26b_visium_lda_fit.log")
cat("GSE298774 LDA fit (recovery) —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

PROC    <- PATHS$processed
RES_GEO <- file.path(PATHS$results, "geo")

set.seed(123)

K          <- 15
PERC_RARE  <- 0.05
PERC_FILT  <- 0.05
BETA_SCALE <- 1

CORPUS_F <- file.path(PROC, "gse298774_corpus.rds")
RAW_F    <- file.path(PROC, "gse298774_lda_raw.rds")
OUT_F    <- file.path(PROC, "gse298774_lda_K15.rds")

# ==============================================================================
# 1. Packages and corpus
# ==============================================================================
hdr("1. Inputs")

if (!requireNamespace("STdeconvolve", quietly = TRUE)) {
  stop("STdeconvolve not available -- it installed during script 26's run, so ",
       "something has changed. Rerun script 26's section 1.")
}
log_line("  STdeconvolve ", as.character(utils::packageVersion("STdeconvolve")),
         logfile = LOG)

if (!file.exists(CORPUS_F)) {
  stop("Saved corpus not found at ", CORPUS_F,
       " -- rerun R/26_visium_deconvolve.R, which rebuilds it.")
}
CC <- readRDS(CORPUS_F)
corpus <- CC$corpus
META   <- CC$meta
log_line("  corpus: ", nrow(corpus), " genes x ", ncol(corpus), " spots",
         logfile = LOG)
log_line("  loaded from the file script 26 wrote, so this is the same matrix",
         logfile = LOG)
log_line("  that was fitted before -- not a rebuild", logfile = LOG)

# The sanity checks fitLDA itself enforces, verified here so a 200-minute run
# cannot fail on them at minute 196.
if (any(Matrix::colSums(corpus) == 0)) {
  stop("Corpus contains spots with no counts; fitLDA will reject it.")
}
if (any(Matrix::rowSums(corpus) == 0)) {
  stop("Corpus contains genes with no counts; fitLDA will reject it.")
}
log_line("  checked: no empty spots, no empty genes", logfile = LOG)

# ==============================================================================
# 2. Fit, and save before touching the result
# ==============================================================================
hdr("2. LDA, K = 15")

if (file.exists(RAW_F)) {
  log_line("  raw fit already on disk -- skipping the fit", logfile = LOG)
  log_line("  (delete ", basename(RAW_F), " to force a refit)", logfile = LOG)
  ldas <- readRDS(RAW_F)
} else {
  log_line("  fitting -- roughly 200 minutes, based on the first run",
           logfile = LOG)
  log_line("  started: ", format(Sys.time()), logfile = LOG)
  t0 <- Sys.time()

  ldas <- STdeconvolve::fitLDA(t(as.matrix(corpus)),
                               Ks = K,
                               perc.rare.thresh = PERC_RARE,
                               plot = FALSE,
                               verbose = TRUE)

  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  # Immediately. Nothing between the fit and the save.
  saveRDS(ldas, RAW_F)

  log_line("  finished in ", round(el, 1), " minutes", logfile = LOG)
  log_line("  SAVED: data/processed/gse298774_lda_raw.rds", logfile = LOG)
  log_line("  the fit is now on disk; everything below is cheap and rerunnable",
           logfile = LOG)
}

# ==============================================================================
# 3. Extract theta and beta
# ==============================================================================
hdr("3. Theta and beta")

mod <- if (!is.null(ldas$models[[as.character(K)]])) {
  ldas$models[[as.character(K)]]
} else {
  ldas$models[[1]]
}
log_line("  model class: ", paste(class(mod), collapse = "/"), logfile = LOG)

res   <- STdeconvolve::getBetaTheta(mod, perc.filt = PERC_FILT,
                                    betaScale = BETA_SCALE)
theta <- res$theta   # spots x topics, each row a composition summing to 1
beta  <- res$beta    # topics x genes

log_line("  theta: ", nrow(theta), " spots x ", ncol(theta), " topics",
         logfile = LOG)
log_line("  beta:  ", nrow(beta), " topics x ", ncol(beta), " genes",
         logfile = LOG)

# Row sums of theta should be 1. If they are not, something is wrong with the
# extraction and script 27's purity threshold would be meaningless.
rs <- range(rowSums(theta))
log_line(sprintf("  theta row sums range [%.4f, %.4f]", rs[1], rs[2]),
         logfile = LOG)
if (abs(rs[1] - 1) > 0.01 || abs(rs[2] - 1) > 0.01) {
  log_line("  [FLAG] theta rows do not sum to 1. Script 27 treats theta as a",
           logfile = LOG)
  log_line("  composition and thresholds on it, so this must be understood",
           logfile = LOG)
  log_line("  before continuing.", logfile = LOG)
}

share <- colMeans(theta)
log_line("", logfile = LOG)
log_line("  mean proportion of tissue per topic:", logfile = LOG)
for (i in order(-share)) {
  log_line(sprintf("    topic %-4s %6.3f", colnames(theta)[i], share[i]),
           logfile = LOG)
}

# ==============================================================================
# 4. Save in the shape script 27 expects
# ==============================================================================
hdr("4. Save")

saveRDS(list(model = mod, theta = theta, beta = beta, meta = META,
             params = list(K = K, perc_rare = PERC_RARE,
                           perc_filt = PERC_FILT,
                           corpus_genes = nrow(corpus),
                           corpus_spots = ncol(corpus),
                           input = "raw integer counts, in-tissue spots",
                           fitted_by = "26b_visium_lda_fit.R")),
        OUT_F)
log_line("  saved: data/processed/gse298774_lda_K15.rds", logfile = LOG)
log_line("  this is the file R/27_visium_topics.R reads", logfile = LOG)

report_warnings(logfile = LOG)
write_session_info("26b_visium_lda_fit")
log_line("\nFull log: ", LOG, logfile = NULL)
