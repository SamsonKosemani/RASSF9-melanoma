# ==============================================================================
# 14_nevus_vs_melanoma.R
#
# Does RASSF9 or RASSF10 expression change between melanocytic nevi and primary
# melanoma? Two INDEPENDENT patient cohorts, same question:
#
#   GSE112509  23 laser-microdissected nevi vs 57 primary melanomas
#              (DESeq2-normalised counts, Ensembl IDs with version suffix)
#   GSE98394   27 common acquired nevi vs 51 primary melanomas
#              (pre-normalised log expression, Ensembl IDs without version)
#
# PRE-SPECIFIED, fixed before any value was inspected:
#   - Primary test: Wilcoxon rank-sum, two-sided, alpha = 0.05
#   - Effect size: Cohen's d with 95% CI
#   - The two cohorts are treated as REPLICATION, not as multiple testing: each
#     is reported in full and no correction is applied across them. A result is
#     persuasive if it reproduces, not if one of the two reaches p < 0.05.
#   - Both genes reported regardless of outcome, per Gustavo's instruction to
#     report negative and discordant results.
#
# DIRECTIONAL CONTEXT (not a prediction for this axis). The TCGA analysis found
# RASSF9 HIGHER in metastatic than primary tumours and RASSF10 LOWER. That is a
# primary-vs-metastatic contrast; nevus-vs-melanoma is a different axis, so
# neither result licenses an expectation here. Stated so that whatever comes
# out is not retrofitted to a prediction afterwards.
#
# Input:  data/raw/geo/GSE112509/GSE112509_DESeq2_normalized_counts.txt.gz
#         data/raw/geo/GSE98394/GSE98394_expression.txt.gz
# Output: results/geo/nevus_vs_melanoma_values.csv
#         results/geo/nevus_vs_melanoma_tests.csv
#         figures/geo/nevus_vs_melanoma.{pdf,png}
#         results/logs/14_nevus_melanoma.log
# ==============================================================================

source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(patchwork)
})
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))

LOG <- file.path(PATHS$logs, "14_nevus_melanoma.log")
cat("Nevus vs melanoma —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

GEO_DIR <- file.path(PATHS$raw, "geo")
RES_GEO <- file.path(PATHS$results, "geo")
FIG_GEO <- file.path(PATHS$figures, "geo")
dir.create(RES_GEO, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_GEO, recursive = TRUE, showWarnings = FALSE)

ENS <- c(RASSF9 = "ENSG00000198774", RASSF10 = "ENSG00000189431")

options(timeout = 1800)

# ==============================================================================
# Helpers
# ==============================================================================

# Column names and GEO sample titles differ in punctuation between the two
# sources (LE.063.SP.064_N vs LE-063-SP-064_N), so matching is done on a
# punctuation-stripped key rather than on the literal string.
norm_key <- function(x) toupper(gsub("[^A-Za-z0-9]", "", as.character(x)))

cohens_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  s <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  (mean(x) - mean(y)) / s
}
cohens_d_ci <- function(d, n1, n2, level = 0.95) {
  se <- sqrt((n1 + n2) / (n1 * n2) + d^2 / (2 * (n1 + n2)))
  z <- qnorm(1 - (1 - level) / 2)
  c(lower = d - z * se, upper = d + z * se)
}

# Pull the two genes out of an expression matrix whose rownames are Ensembl IDs
# (version suffix tolerated).
extract_genes <- function(mat, label) {
  ids <- sub("\\..*$", "", rownames(mat))
  out <- list()
  for (g in names(ENS)) {
    idx <- which(ids == ENS[[g]])
    if (length(idx) == 0) {
      log_line("    [FLAG] ", g, " (", ENS[[g]], ") not found in ", label,
               logfile = LOG)
      next
    }
    if (length(idx) > 1) {
      log_line("    [note] ", g, " matched ", length(idx),
               " rows; summing them", logfile = LOG)
      v <- colSums(mat[idx, , drop = FALSE])
    } else {
      v <- as.numeric(mat[idx, ])
    }
    out[[g]] <- setNames(v, colnames(mat))
    log_line("    [ ok ] ", g, " extracted from ", label, logfile = LOG)
  }
  out
}

# ==============================================================================
# 1. GSE112509
# ==============================================================================
hdr("1. GSE112509 (23 nevi vs 57 primary melanomas)")

f1 <- list.files(file.path(GEO_DIR, "GSE112509"), pattern = "normalized_counts",
                 full.names = TRUE)
if (!length(f1)) stop("GSE112509 counts file not found in ", file.path(GEO_DIR, "GSE112509"))

m1 <- as.data.frame(readr::read_tsv(f1[1], show_col_types = FALSE))
rownames(m1) <- m1[[1]]; m1 <- as.matrix(m1[, -1, drop = FALSE])
log_line("  matrix: ", nrow(m1), " genes x ", ncol(m1), " samples", logfile = LOG)
log_line("  value range: ", round(min(m1, na.rm = TRUE), 2), " to ",
         round(max(m1, na.rm = TRUE), 1),
         " -> DESeq2 normalised counts, log2(x+1) applied", logfile = LOG)

e1 <- extract_genes(m1, "GSE112509")

# Groups: the series matrix carries "tissue: benign melanocytic nevi" /
# "primary melanoma". Match it to the count columns; the _N / _M suffix in the
# column names is used as an independent cross-check, not as the primary source.
gse1 <- getGEO("GSE112509", destdir = GEO_DIR, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
pd1 <- pData(gse1)
tis1 <- grep("tissue", names(pd1), value = TRUE, ignore.case = TRUE)
tis1 <- tis1[grepl("nevi|melanoma", tolower(sapply(pd1[tis1], function(z) paste(z, collapse = "")))) ][1]

map1 <- tibble::tibble(
  key   = norm_key(pd1$title),
  group = ifelse(grepl("nevi|nevus|naev", tolower(as.character(pd1[[tis1]]))),
                 "Nevus", "Melanoma")
)

d1 <- tibble::tibble(sample = names(e1[[1]]), key = norm_key(names(e1[[1]]))) |>
  dplyr::left_join(map1, by = "key") |>
  dplyr::mutate(suffix_group = dplyr::case_when(
    grepl("_N$", sample) ~ "Nevus",
    grepl("_M$", sample) ~ "Melanoma",
    TRUE ~ NA_character_))

log_line("\n  group assignment:", logfile = LOG)
log_line("    matched to series matrix: ", sum(!is.na(d1$group)), " of ",
         nrow(d1), logfile = LOG)
log_line("    from _N/_M suffix:        ", sum(!is.na(d1$suffix_group)),
         logfile = LOG)
agree <- sum(d1$group == d1$suffix_group, na.rm = TRUE)
both  <- sum(!is.na(d1$group) & !is.na(d1$suffix_group))
log_line("    the two agree on ", agree, " of ", both, " samples where both exist",
         logfile = LOG)
if (both > 0 && agree < both) {
  log_line("    [FLAG] disagreement -- series-matrix labels are used",
           logfile = LOG)
}
d1$group <- dplyr::coalesce(d1$group, d1$suffix_group)

for (g in names(e1)) d1[[g]] <- log2(e1[[g]][d1$sample] + 1)
d1$dataset <- "GSE112509"

tab1 <- as.data.frame(table(d1$group, useNA = "ifany"))
names(tab1) <- c("group", "n")
cat(capture.output(print(tab1, row.names = FALSE)), sep = "\n",
    file = LOG, append = TRUE)
print(tab1)

# ==============================================================================
# 2. GSE98394
# ==============================================================================
hdr("2. GSE98394 (27 nevi vs 51 primary melanomas)")

# GSE98394 deposits TWO matrices. The "expression" file is a pre-normalised
# log-scale matrix but covers only ~12,800 genes -- a filtered subset that does
# NOT include RASSF10. The "counts" file carries the full gene set. We therefore
# read both: a gene is taken from the expression file when present there, and
# otherwise from the counts file as log2 counts-per-million. Which source was
# used for each gene is reported, because the two are on different scales and
# that must not be hidden.

read_matrix <- function(path) {
  x <- as.data.frame(readr::read_tsv(path, show_col_types = FALSE))
  rownames(x) <- x[[1]]
  as.matrix(x[, -1, drop = FALSE])
}

f2e <- list.files(file.path(GEO_DIR, "GSE98394"), pattern = "expression",
                  full.names = TRUE)
f2c <- list.files(file.path(GEO_DIR, "GSE98394"), pattern = "counts",
                  full.names = TRUE)

m2e <- if (length(f2e)) read_matrix(f2e[1]) else NULL
m2c <- if (length(f2c)) read_matrix(f2c[1]) else NULL

if (!is.null(m2e)) {
  log_line("  expression file: ", nrow(m2e), " genes x ", ncol(m2e),
           " samples, range ", round(min(m2e, na.rm = TRUE), 2), " to ",
           round(max(m2e, na.rm = TRUE), 2), " (log scale)", logfile = LOG)
}
if (!is.null(m2c)) {
  log_line("  counts file:     ", nrow(m2c), " genes x ", ncol(m2c),
           " samples, max ", round(max(m2c, na.rm = TRUE), 0), logfile = LOG)
}

# --- per gene, pick a source and put it on a log scale ------------------------
e2 <- list(); src2 <- character(0)
for (g in names(ENS)) {
  got <- NULL
  if (!is.null(m2e)) {
    idx <- which(sub("\\..*$", "", rownames(m2e)) == ENS[[g]])
    if (length(idx)) {
      got <- colSums(m2e[idx, , drop = FALSE])
      src2[g] <- "expression file (log scale, as provided)"
    }
  }
  if (is.null(got) && !is.null(m2c)) {
    idx <- which(sub("\\..*$", "", rownames(m2c)) == ENS[[g]])
    if (length(idx)) {
      lib <- colSums(m2c, na.rm = TRUE)
      cpm <- colSums(m2c[idx, , drop = FALSE]) / lib * 1e6
      got <- log2(cpm + 1)
      src2[g] <- "counts file (log2 CPM, computed here)"
    }
  }
  if (is.null(got)) {
    log_line("    [FLAG] ", g, " not found in either GSE98394 matrix",
             logfile = LOG)
    next
  }
  e2[[g]] <- got
  log_line("    [ ok ] ", g, " from ", src2[g], logfile = LOG)
}
if (!length(e2)) stop("Neither gene could be extracted from GSE98394.")

# --- sample matching ---------------------------------------------------------
gse2 <- getGEO("GSE98394", destdir = GEO_DIR, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
pd2 <- pData(gse2)
tis2 <- grep("tissue", names(pd2), value = TRUE, ignore.case = TRUE)
tis2 <- tis2[grepl("nevus|nevi|melanoma",
                   tolower(sapply(pd2[tis2], function(z) paste(z, collapse = ""))))][1]
if (is.na(tis2)) stop("No tissue column found in GSE98394 metadata.")

# The JC### ids are embedded inside longer metadata strings, so an exact match
# fails. Score each column by how many expression headers appear as substrings.
cols2 <- names(pd2)
hdrs  <- names(e2[[1]])
score <- sapply(cols2, function(cc) {
  v <- as.character(pd2[[cc]])
  sum(vapply(hdrs, function(h) any(grepl(h, v, fixed = TRUE)), logical(1)))
})
best <- cols2[which.max(score)]
log_line("\n  best metadata column for sample ids: ", best,
         " (matches ", max(score), " of ", length(hdrs), " headers)",
         logfile = LOG)
log_line("  example value: ", substr(as.character(pd2[[best]])[1], 1, 90),
         logfile = LOG)

if (max(score) < length(hdrs) * 0.5) {
  log_line("  [FLAG] weak matching; metadata columns are: ",
           paste(head(cols2, 20), collapse = ", "), logfile = LOG)
  stop("Could not reliably match GSE98394 samples. Headers: ",
       paste(head(hdrs, 3), collapse = ", "))
}

grp2 <- ifelse(grepl("nevus|nevi|naev", tolower(as.character(pd2[[tis2]]))),
               "Nevus", "Melanoma")
vals2 <- as.character(pd2[[best]])

d2 <- tibble::tibble(sample = hdrs) |>
  dplyr::mutate(group = vapply(sample, function(h) {
    i <- which(grepl(h, vals2, fixed = TRUE))
    if (length(i) == 1) grp2[i] else NA_character_
  }, character(1)))

log_line("  matched to metadata: ", sum(!is.na(d2$group)), " of ", nrow(d2),
         logfile = LOG)

for (g in names(e2)) d2[[g]] <- as.numeric(e2[[g]][d2$sample])
d2$dataset <- "GSE98394"
d2 <- d2 |> dplyr::filter(!is.na(group))

tab2 <- as.data.frame(table(d2$group))
names(tab2) <- c("group", "n")
cat(capture.output(print(tab2, row.names = FALSE)), sep = "\n",
    file = LOG, append = TRUE)
print(tab2)

# ==============================================================================
# 3. Tests
# ==============================================================================
hdr("3. Nevus vs melanoma")

all_d <- dplyr::bind_rows(
  d1 |> dplyr::select(dataset, sample, group, dplyr::any_of(names(ENS))),
  d2 |> dplyr::select(dataset, sample, group, dplyr::any_of(names(ENS)))
) |> dplyr::filter(!is.na(group))

rows <- list()
for (ds in unique(all_d$dataset)) {
  log_line("\n  ", ds, ":", logfile = LOG)
  dd <- all_d |> dplyr::filter(dataset == ds)
  for (g in names(ENS)) {
    if (!g %in% names(dd)) next
    nv <- dd[[g]][dd$group == "Nevus"];    nv <- nv[!is.na(nv)]
    ml <- dd[[g]][dd$group == "Melanoma"]; ml <- ml[!is.na(ml)]
    if (length(nv) < 3 || length(ml) < 3) next

    w  <- suppressWarnings(wilcox.test(ml, nv, exact = FALSE))
    dd_eff <- cohens_d(ml, nv)
    ci <- cohens_d_ci(dd_eff, length(ml), length(nv))

    rows[[length(rows) + 1]] <- tibble::tibble(
      dataset = ds, gene = g,
      n_nevus = length(nv), n_melanoma = length(ml),
      median_nevus = median(nv), median_melanoma = median(ml),
      W = unname(w$statistic), p_value = w$p.value,
      cohens_d = dd_eff, d_lower = unname(ci[1]), d_upper = unname(ci[2]),
      direction = ifelse(median(ml) < median(nv),
                         "lower in melanoma", "higher in melanoma"))

    log_line(sprintf(
      "    %-8s n=%d nevi / %d melanoma  median %.2f vs %.2f  W=%.0f  p=%s  d=%+.3f [%+.2f, %+.2f]  (%s)",
      g, length(nv), length(ml), median(nv), median(ml), w$statistic,
      fmt_p(w$p.value), dd_eff, ci[1], ci[2],
      rows[[length(rows)]]$direction), logfile = LOG)
  }
}
res <- dplyr::bind_rows(rows)

# ==============================================================================
# 4. Replication check
# ==============================================================================
hdr("3b. Expression sources used")
log_line("  GSE112509: DESeq2-normalised counts, log2(x+1)", logfile = LOG)
for (g in names(src2)) log_line("  GSE98394 ", g, ": ", src2[g], logfile = LOG)
log_line("  The two cohorts are on different scales, so effect sizes are",
         logfile = LOG)
log_line("  comparable in direction but not in absolute magnitude.",
         logfile = LOG)

hdr("4. Does it replicate?")

for (g in names(ENS)) {
  r <- res |> dplyr::filter(gene == g)
  if (nrow(r) < 2) next
  same_dir <- length(unique(r$direction)) == 1
  both_sig <- all(r$p_value < 0.05)
  log_line(sprintf("  %-8s %s | %s",
                   g,
                   ifelse(same_dir, "same direction in both cohorts",
                          "DIRECTIONS DISAGREE"),
                   ifelse(both_sig, "significant in both",
                          paste0("significant in ", sum(r$p_value < 0.05),
                                 " of ", nrow(r)))), logfile = LOG)
}

log_line("\n  A finding that reproduces across two independent cohorts is",
         logfile = LOG)
log_line("  worth more than one that is significant in a single one. Where the",
         logfile = LOG)
log_line("  two disagree, report both and say so rather than choosing.",
         logfile = LOG)

# ==============================================================================
# 5. Figure
# ==============================================================================
hdr("5. Figure")

panel <- function(ds, g) {
  # A gene can be absent from one cohort's matrix (RASSF10 is not in the
  # GSE98394 files). Draw an explicit "not measured" placeholder rather than an
  # empty panel, so the gap is visible in the figure instead of being silent.
  if (!g %in% names(all_d) ||
      !any(all_d$dataset == ds & !is.na(all_d[[g]]))) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, size = 2.9, colour = PAL$grey_pt,
                 label = paste0(g, " not present\nin ", ds)) +
        theme_void()
    )
  }
  r <- res |> dplyr::filter(dataset == ds, gene == g)
  pd <- all_d |> dplyr::filter(dataset == ds, !is.na(.data[[g]]))
  pd$group <- factor(pd$group, levels = c("Nevus", "Melanoma"))
  ymax <- max(pd[[g]], na.rm = TRUE); ymin <- min(pd[[g]], na.rm = TRUE)
  pad <- (ymax - ymin) * 0.12

  ggplot(pd, aes(x = group, y = .data[[g]], fill = group)) +
    geom_violin(trim = TRUE, linewidth = 0.3, colour = NA, alpha = 0.30,
                scale = "width") +
    geom_jitter(width = 0.11, size = 0.6, alpha = 0.45, colour = PAL$grey_pt,
                stroke = 0) +
    geom_boxplot(width = 0.20, outlier.shape = NA, linewidth = 0.45,
                 colour = PAL$ink, fill = "white") +
    annotate("segment", x = 1, xend = 2, y = ymax + pad * 0.5,
             yend = ymax + pad * 0.5, linewidth = 0.45, colour = PAL$ink) +
    annotate("text", x = 1.5, y = ymax + pad * 0.7,
             label = if (nrow(r)) fmt_p_plot(r$p_value) else "", parse = TRUE,
             size = 2.9, vjust = 0, colour = PAL$ink) +
    scale_fill_manual(values = c("Nevus" = PAL$uv_low,
                                 "Melanoma" = PAL$uv_high), guide = "none") +
    coord_cartesian(ylim = c(ymin - pad * 0.3, ymax + pad * 1.6)) +
    labs(x = NULL, y = bquote(bold(italic(.(g))~"expression")),
         subtitle = paste0(ds, " · ", g)) +
    theme_pub()
}

fig <- (panel("GSE112509", "RASSF9")  | panel("GSE98394", "RASSF9")) /
       (panel("GSE112509", "RASSF10") | panel("GSE98394", "RASSF10")) +
  plot_annotation(
    title = "RASSF9 and RASSF10 in melanocytic nevi versus primary melanoma",
    caption = "Two independent cohorts. Wilcoxon rank-sum, two-sided.",
    tag_levels = "A") &
  theme(plot.tag = element_text(size = FIG$base_size + 3, face = "bold"),
        plot.caption = element_text(size = FIG$base_size - 1,
                                    colour = PAL$grey_pt, hjust = 0))

for (fmt in c("pdf", "png")) {
  f <- file.path(FIG_GEO, paste0("nevus_vs_melanoma.", fmt))
  if (file.exists(f)) file.remove(f)
  args <- list(filename = f, plot = fig, width = 160, height = 150,
               units = "mm")
  if (fmt == "pdf") args$device <- grDevices::cairo_pdf
  else { args$dpi <- 600; args$bg <- "white" }
  do.call(ggplot2::ggsave, args)
}
log_line("  Wrote: figures/geo/nevus_vs_melanoma.{pdf,png}", logfile = LOG)

# ==============================================================================
# 6. Output
# ==============================================================================
hdr("6. Output")

readr::write_csv(all_d, file.path(RES_GEO, "nevus_vs_melanoma_values.csv"))
readr::write_csv(res,   file.path(RES_GEO, "nevus_vs_melanoma_tests.csv"))
cat(capture.output(print(as.data.frame(res), row.names = FALSE, digits = 3)),
    sep = "\n", file = LOG, append = TRUE)
print(as.data.frame(res), digits = 3)

log_line("\n  Wrote: results/geo/nevus_vs_melanoma_values.csv", logfile = LOG)
log_line("  Wrote: results/geo/nevus_vs_melanoma_tests.csv", logfile = LOG)

report_warnings(LOG)
write_session_info("14_nevus_melanoma")
log_line("\nFull log: ", LOG, logfile = LOG)
