# ==============================================================================
# 07_download_methylation.R
#
# Purpose: download TCGA-SKCM Illumina HM450 methylation beta values, identify
# the probes annotated to RASSF9 and RASSF10, and report what is available.
# Download and inventory only -- no hypothesis is tested here.
#
# This is slow on a first run (large download). The script detects an existing
# download and skips to reporting, so re-running is cheap.
#
# PRE-SPECIFIED PROBE SELECTION (decisions D11-D12, recorded before any
# methylation-expression association was examined):
#
#   D11 Promoter-proximal probes are those whose UCSC RefGene annotation names
#       the gene AND whose RefGene group is TSS1500, TSS200, 5'UTR or 1stExon.
#       Gene-body probes are retained separately for a secondary analysis but
#       are not part of the promoter definition: body methylation often
#       correlates positively with expression and mixing the two would blur
#       the silencing hypothesis being tested.
#
#   D12 The promoter summary for each sample is the MEAN beta across its
#       promoter-proximal probes. Per-probe results are also reported, so the
#       summary never hides a divergent probe.
#
# Input:  data/processed/analysis_cohort.csv
# Output: data/raw/skcm_hm450.rds
#         data/processed/methylation_probes.csv   (probe annotation)
#         data/processed/methylation_betas.csv    (per sample, per probe)
#         results/logs/07_methylation_download.log
# ==============================================================================

source(here::here("R", "00_setup.R"))

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})
# S4Vectors masks several dplyr verbs; put dplyr back in front.
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE))

LOG <- file.path(PATHS$logs, "07_methylation_download.log")
cat("HM450 methylation download —", format(Sys.time()), "\n\n", file = LOG)

hdr <- function(txt) log_line("\n== ", txt, " ",
                              strrep("=", max(0, 60 - nchar(txt))), logfile = LOG)

GENES <- c("RASSF9", "RASSF10")
PROMOTER_GROUPS <- c("TSS1500", "TSS200", "5'UTR", "1stExon")

ME_FILE <- file.path(PATHS$raw, "skcm_hm450.rds")
GDC_DIR <- file.path(PATHS$raw, "GDCdata")

# ==============================================================================
# 1. Download or load
# ==============================================================================
hdr("1. Methylation data")

if (file.exists(ME_FILE)) {
  log_line("Found existing file: ", ME_FILE, logfile = LOG)
  log_line("Size: ", round(file.info(ME_FILE)$size / 1024^2, 1), " MB",
           logfile = LOG)
  log_line("Loading it instead of re-downloading.", logfile = LOG)
  me <- readRDS(ME_FILE)
} else {
  log_line("Querying the GDC for TCGA-SKCM HM450 beta values...", logfile = LOG)

  query <- GDCquery(
    project       = "TCGA-SKCM",
    data.category = "DNA Methylation",
    data.type     = "Methylation Beta Value",
    platform      = "Illumina Human Methylation 450"
  )

  n_files <- nrow(getResults(query))
  log_line("Files matched: ", n_files, logfile = LOG)
  if (n_files == 0) {
    stop("GDC query returned no HM450 files for TCGA-SKCM. Check connectivity ",
         "and that the platform string is still current.")
  }

  log_line("Downloading (this is the slow step)...", logfile = LOG)
  GDCdownload(query, directory = GDC_DIR, method = "api", files.per.chunk = 10)

  log_line("Assembling SummarizedExperiment...", logfile = LOG)
  me <- GDCprepare(query, directory = GDC_DIR)

  saveRDS(me, ME_FILE)
  log_line("Saved: ", ME_FILE, logfile = LOG)
}

log_line("\nProbes:  ", nrow(me), logfile = LOG)
log_line("Samples: ", ncol(me), logfile = LOG)

# ==============================================================================
# 2. Probe annotation
# ==============================================================================
hdr("2. Probe annotation")

# AMENDED D11 (see note below): the GDC methylation object is annotated on
# hg38 by sesame and carries probe coordinates but no RefGene region labels.
# The Bioconductor HM450 annotation package would supply region labels but is
# built on hg19, and mixing genome builds is a silent-error hazard. We
# therefore define the promoter by DISTANCE TO THE TRANSCRIPTION START SITE in
# hg38, taking gene coordinates from the same expression object used for the
# expression analysis. One build throughout, and an explicit criterion rather
# than an inherited one.

anno <- as.data.frame(rowRanges(me))
anno$probe <- rownames(me)

gene_col <- intersect(c("gene_HGNC", "gene", "UCSC_RefGene_Name"),
                      names(anno))[1]
if (is.na(gene_col)) {
  stop("No gene-annotation column found in the methylation object; columns: ",
       paste(names(anno), collapse = ", "))
}

chr_col <- if ("seqnames" %in% names(anno)) "seqnames" else "chrm_A"
pos_col <- if ("start" %in% names(anno)) "start" else "probeBeg"

log_line("Probe coordinates: ", chr_col, " / ", pos_col, " (hg38)", logfile = LOG)
log_line("Gene labels:       ", gene_col, logfile = LOG)

# --- gene coordinates, from the expression object ----------------------------
se_expr <- readRDS(file.path(PATHS$raw, "skcm_star_counts.rds"))
gr <- as.data.frame(rowRanges(se_expr))
gr$gene_name <- as.character(rowData(se_expr)$gene_name)

tss_tbl <- list()
for (g in GENES) {
  hit <- gr[gr$gene_name == g, , drop = FALSE]
  if (nrow(hit) == 0) {
    log_line("  [FLAG] ", g, " not found in the expression object", logfile = LOG)
    next
  }
  h <- hit[1, ]
  strand_chr <- as.character(h$strand)
  tss <- if (strand_chr == "-") h$end else h$start
  tss_tbl[[g]] <- tibble::tibble(
    gene = g, chr = as.character(h$seqnames), tss = tss,
    strand = strand_chr, gene_start = h$start, gene_end = h$end
  )
  log_line(sprintf("  %-8s %s:%s-%s (%s strand), TSS at %s", g,
                   h$seqnames, format(h$start, big.mark = ","),
                   format(h$end, big.mark = ","), strand_chr,
                   format(tss, big.mark = ",")), logfile = LOG)
}
tss_tbl <- dplyr::bind_rows(tss_tbl)

# ==============================================================================
# 3. Probe selection (D11, amended)
# ==============================================================================
hdr("3. Probe selection (D11, TSS-distance definition)")

# Promoter window: 1500 bp upstream to 500 bp downstream of the TSS. This is
# the conventional definition and matches the span covered by the TSS1500 and
# TSS200 categories of the RefGene annotation we could not use.
PROMOTER_UP   <- 1500
PROMOTER_DOWN <- 500

log_line("Promoter window: TSS ", -PROMOTER_UP, " to +", PROMOTER_DOWN,
         " bp (strand-aware)", logfile = LOG)
log_line("Gene-body probes are retained but analysed separately.", logfile = LOG)

# Semicolon-delimited, repeated per transcript. Match whole tokens so RASSF9
# can never match RASSF9P1 or similar.
has_gene <- function(x, gene) {
  vapply(strsplit(as.character(x), ";", fixed = TRUE),
         function(v) gene %in% trimws(v), logical(1))
}

probe_rows <- list()
for (g in GENES) {
  if (!g %in% tss_tbl$gene) next
  ti <- tss_tbl[tss_tbl$gene == g, ]

  # Probes labelled with the gene, OR lying within 10 kb of the locus (which
  # catches upstream promoter probes that carry no gene label).
  by_label <- has_gene(anno[[gene_col]], g)
  by_locus <- as.character(anno[[chr_col]]) == ti$chr &
    !is.na(anno[[pos_col]]) &
    anno[[pos_col]] > (min(ti$gene_start, ti$tss) - 10000) &
    anno[[pos_col]] < (max(ti$gene_end,   ti$tss) + 10000)

  idx <- which(by_label | by_locus)
  if (length(idx) == 0) {
    log_line("  [FLAG] ", g, ": no probes found", logfile = LOG)
    next
  }
  sub <- anno[idx, , drop = FALSE]

  # Signed distance from TSS, positive = downstream (into the gene).
  raw_d <- sub[[pos_col]] - ti$tss
  dist  <- if (ti$strand == "-") -raw_d else raw_d

  probe_rows[[g]] <- tibble::tibble(
    gene       = g,
    probe      = sub$probe,
    chr        = as.character(sub[[chr_col]]),
    pos        = sub[[pos_col]],
    dist_tss   = dist,
    labelled   = by_label[idx],
    region     = dplyr::case_when(
      dist >= -PROMOTER_UP & dist <= PROMOTER_DOWN ~ "Promoter",
      dist > PROMOTER_DOWN & sub[[pos_col]] >= min(ti$gene_start, ti$gene_end) &
        sub[[pos_col]] <= max(ti$gene_start, ti$gene_end)  ~ "Gene body",
      dist < -PROMOTER_UP                          ~ "Upstream",
      TRUE                                         ~ "Downstream"
    ),
    promoter   = dist >= -PROMOTER_UP & dist <= PROMOTER_DOWN,
    on_array   = sub$probe %in% rownames(me)
  ) |>
    dplyr::arrange(dist_tss)

  pr <- probe_rows[[g]]
  log_line(sprintf("\n  %-8s %d probes at the locus (%d gene-labelled)",
                   g, nrow(pr), sum(pr$labelled)), logfile = LOG)
  reg <- as.data.frame(table(pr$region))
  names(reg) <- c("region", "n")
  cat(capture.output(print(reg, row.names = FALSE)), sep = "\n",
      file = LOG, append = TRUE)

  prom_here <- pr |> dplyr::filter(promoter, on_array)
  if (nrow(prom_here) > 0) {
    log_line("    promoter probes (", nrow(prom_here), "):", logfile = LOG)
    for (k in seq_len(nrow(prom_here))) {
      log_line(sprintf("      %-12s %s:%s  %+d bp from TSS",
                       prom_here$probe[k], prom_here$chr[k],
                       format(prom_here$pos[k], big.mark = ","),
                       prom_here$dist_tss[k]), logfile = LOG)
    }
  } else {
    log_line("    [FLAG] no promoter probes within the window for ", g,
             logfile = LOG)
    near <- pr |> dplyr::arrange(abs(dist_tss)) |> head(3)
    log_line("    nearest probes: ",
             paste(sprintf("%s (%+d bp)", near$probe, near$dist_tss),
                   collapse = ", "), logfile = LOG)
  }
}

probe_tbl <- dplyr::bind_rows(probe_rows)
readr::write_csv(probe_tbl, file.path(PATHS$processed, "methylation_probes.csv"))

# ==============================================================================
# 4. Extract beta values and match to the cohort
# ==============================================================================
hdr("4. Beta values")

keep_probes <- probe_tbl |> dplyr::filter(on_array) |> dplyr::pull(probe) |> unique()
log_line("Extracting ", length(keep_probes), " probes across ", ncol(me),
         " samples", logfile = LOG)

betas <- assay(me)[keep_probes, , drop = FALSE]

cd <- colData(me)
bc <- if ("barcode" %in% names(cd)) as.character(cd$barcode) else colnames(me)

beta_long <- as.data.frame(t(betas))
beta_long$barcode <- bc
beta_long$patient <- tcga_patient(bc)
beta_long$sample_id <- substr(bc, 1, 16)

# Missingness per probe -- HM450 probes can fail in a subset of samples.
log_line("\nMissingness per probe:", logfile = LOG)
for (p in keep_probes) {
  miss <- mean(is.na(betas[p, ]))
  log_line(sprintf("  %-12s %.1f%% missing", p, 100 * miss), logfile = LOG)
}

# One record per patient, mirroring decision D9.
coh <- readr::read_csv(file.path(PATHS$processed, "analysis_cohort.csv"))
coh_samples <- substr(coh$rna_barcode, 1, 16)

beta_dedup <- beta_long |>
  dplyr::mutate(.same_sample = as.integer(sample_id %in% coh_samples)) |>
  dplyr::arrange(patient, dplyr::desc(.same_sample), barcode) |>
  dplyr::group_by(patient) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(-.same_sample)

log_line("\nMethylation samples: ", nrow(beta_long), " -> ", nrow(beta_dedup),
         " one per patient", logfile = LOG)

matched <- dplyr::inner_join(
  coh |> dplyr::select(patient),
  beta_dedup, by = "patient"
)
log_line("Patients with both expression and methylation: ", nrow(matched),
         " of ", nrow(coh), logfile = LOG)

same_vial <- mean(matched$sample_id %in% coh_samples)
log_line("Methylation from the same sample vial as RNA: ",
         round(100 * same_vial, 1), "%", logfile = LOG)

readr::write_csv(beta_dedup,
                 file.path(PATHS$processed, "methylation_betas.csv"))
log_line("\nWrote: data/processed/methylation_probes.csv", logfile = LOG)
log_line("Wrote: data/processed/methylation_betas.csv", logfile = LOG)

hdr("5. Warnings")
report_warnings(LOG)

write_session_info("07_methylation_download")

hdr("DONE")
log_line("Next: script 08 tests methylation against expression using the ",
         "probes inventoried above. Review the probe counts and regions in ",
         "this log before running it -- if a gene has only one or two usable ",
         "promoter probes, that is a real limitation and belongs in the ",
         "manuscript rather than being discovered later.", logfile = LOG)
log_line("Full log: ", LOG, logfile = LOG)
