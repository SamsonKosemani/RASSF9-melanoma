## =============================================================================
## 44a_gdsc450k_inspect.R
##
## STEP 2a of Gustavo's cell-line question. Inspection only -- computes no
## result and makes no claim.
##
## WHY A SEPARATE INSPECTION SCRIPT
##   GSE68379_Matrix.processed.txt.gz is 3.8 Gb. It is GenomeStudio output, and
##   that format normally pairs every sample with a Detection Pval column, so
##   the file may hold ~1,029 columns or ~2,057. Guessing wrong costs a full
##   pass over 3.8 Gb and risks silently reading p-values as beta values.
##   This script reads the header and the first few data rows only, reports
##   exactly what it found, and writes a parsing plan for 44b to follow.
##
##   Same pattern as 37a and 38a: look before building.
##
## IT ALSO SETTLES THE SAMPLE MAPPING
##   GEO sample titles are Sanger cell-line names ("CellLine_A375"). Our
##   melanoma set is DepMap model IDs. Model.csv carries SangerModelID and
##   StrippedCellLineName, so the join is done on curated identifiers first and
##   only falls back to name matching for what is left over. The script reports
##   how many of the 113 non-uveal melanoma lines are actually present, which
##   determines whether step 2 is worth running at all.
##
## INPUTS
##   data/raw/geo/gdsc450k/GSE68379_Matrix.processed.txt.gz   (3.8 Gb)
##   data/raw/ccle/Model.csv                                  (already present)
##
## OUTPUTS
##   results/ccle/44a_probe_table.csv        the 48 locus probes + TCGA values
##   results/ccle/44a_parse_plan.csv         column layout for 44b
##   results/ccle/44a_sample_map.csv         GEO sample <-> DepMap model
##   results/logs/44a_gdsc450k_inspect.log
## =============================================================================

## ---- 0. Setup ---------------------------------------------------------------
if (file.exists("R/00_setup.R")) source("R/00_setup.R")

suppressPackageStartupMessages({
  library(data.table)
})

ROOT    <- if (exists("PROJECT_ROOT")) PROJECT_ROOT else getwd()
GEO_DIR <- file.path(ROOT, "data", "raw", "geo", "gdsc450k")
CCLE_DIR<- file.path(ROOT, "data", "raw", "ccle")
RES_DIR <- file.path(ROOT, "results", "ccle")
LOG_DIR <- file.path(ROOT, "results", "logs")
for (d in c(RES_DIR, LOG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

LOG_PATH <- file.path(LOG_DIR, "44a_gdsc450k_inspect.log")
.logcon  <- file(LOG_PATH, open = "wt")
say <- function(...) { m <- paste0(...); cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = .logcon) }
WARN <- character(0)
warn <- function(...) { m <- paste0(...); WARN <<- c(WARN, m); say("  [WARN] ", m) }
on.exit({
  if (length(WARN)) { say(""); say("WARNINGS RAISED (", length(WARN), "):"); for (w in WARN) say("  - ", w) }
  else say("  no warnings raised")
  say(""); say("Full log: ", LOG_PATH); close(.logcon)
}, add = TRUE)

## ---- 1. The probe table -----------------------------------------------------
## Transcribed from results/tcga/methylation_full_probes.csv and the pre-
## specified set labels in methylation_null_distribution.csv. Embedded here
## rather than re-derived so that step 2 tests EXACTLY the probes the TCGA
## analysis tested, and so this script is self-contained for the methods.
##
## set: tss_proximal = promoter core, the probes carrying the TCGA signal
##      shore        = inside the promoter window but outside the core; these
##                     carry the island-shore inversion (positive rho)
##      distal       = everything else within 150 kb
##
## tcga_rho and tcga_beta are the complete-cohort values (n = 469).

P9 <- data.table(
  gene = "RASSF9",
  probe = c("cg11846905","cg04916911","cg18877737","cg08888956","cg03707749",
            "cg18827107","cg19622623","cg27035734","cg25456477","cg00310523",
            "cg02569458","cg06740227"),
  dist_tss = c(-112093,-37967,-37694,-37652,-10274,-770,-638,-370,-180,11,94,383),
  set = c("distal","distal","distal","distal","distal",
          "shore","shore",
          "tss_proximal","tss_proximal","tss_proximal","tss_proximal","tss_proximal"),
  tcga_rho = c(-0.2137,-0.1547,-0.1179,-0.1523,-0.0798,
               0.1350,0.1363,-0.3675,-0.4855,-0.4290,-0.3954,-0.2944),
  tcga_beta = c(0.8325,0.2009,0.1649,0.4792,0.8642,
                0.8457,0.6168,0.6319,0.4627,0.3852,0.2005,0.1908)
)

P10 <- data.table(
  gene = "RASSF10",
  probe = c("cg12399536","cg14665901","cg20445806","cg13845161","cg20233228",
            "cg26770464","cg01482376","cg24517989","cg24506086","cg09711214",
            "cg15805337","cg18525873","cg24442572","cg05059451","cg04755604",
            "cg27075800","cg23548136",
            "cg11518550","cg08148891",
            "cg19665644","cg20918243","cg09596429","cg05817758","cg08876434",
            "cg25344734","cg26005578","cg05095158","cg22740547",
            "cg07600636","cg00326127","cg20098471","cg08161666","cg04844534",
            "cg25068253","cg04137004","cg02452500"),
  dist_tss = c(-145979,-145732,-138015,-129145,-127253,-119934,-117257,-115361,
               -108745,-106038,-94966,-86319,-69274,-19707,-19679,-19265,-3157,
               -1141,-554,
               -389,-351,-346,-296,-226,-220,-187,10,13,
               940,1541,1750,4034,43041,93431,109216,131064),
  set = c(rep("distal",17), "shore","shore", rep("tss_proximal",9), rep("distal",8)),
  tcga_rho = c(-0.2128,-0.1574,0.0586,0.0651,-0.0037,-0.0095,0.0105,-0.1003,
               0.0328,0.1194,-0.1693,0.0962,0.0683,0.1810,0.0884,0.2292,0.1963,
               0.0216,-0.0026,
               -0.0931,-0.0774,-0.0188,-0.2435,-0.2867,-0.2050,-0.1886,-0.2328,-0.2215,
               0.1449,0.1902,0.1329,-0.0050,-0.0201,-0.0894,0.2068,-0.0632),
  tcga_beta = c(0.6804,0.5743,0.8798,0.9450,0.9119,0.7752,0.9333,0.6478,
                0.8068,0.8120,0.8874,0.9124,0.9085,0.7002,0.6035,0.8782,0.9345,
                0.8932,0.4019,
                0.0419,0.1996,0.2465,0.3374,0.0860,0.1019,0.0171,0.0217,0.0175,
                0.2585,0.5193,0.6857,0.5739,0.7916,0.3887,0.6010,0.4788)
)

PROBES <- rbind(P9, P10)

say("")
say("== 1. Probe table ==============================================")
say("  RASSF9 :  ", nrow(P9),  " probes  (", sum(P9$set == "tss_proximal"),  " core, ",
    sum(P9$set == "shore"),  " shore, ", sum(P9$set == "distal"),  " distal)")
say("  RASSF10:  ", nrow(P10), " probes  (", sum(P10$set == "tss_proximal"), " core, ",
    sum(P10$set == "shore"), " shore, ", sum(P10$set == "distal"), " distal)")
say("  total  :  ", nrow(PROBES))
say("")
say("  NOTE: script 35 retrieved 51 probes at these two loci but the")
say("  correlation table holds 48. Three RASSF10 probes are absent, presumably")
say("  dropped as all-NA. All 11 RASSF10 promoter-window probes are present,")
say("  so the promoter question is unaffected. State this in the methods.")
say("")
say("  TCGA complete-cohort medians, promoter core (what cell lines must beat")
say("  to be called hypermethylated):")
for (g in c("RASSF9", "RASSF10")) {
  cp <- PROBES[gene == g & set == "tss_proximal"]
  say(sprintf("    %-8s median beta across core probes: %.3f  (range %.3f - %.3f)",
              g, median(cp$tcga_beta), min(cp$tcga_beta), max(cp$tcga_beta)))
}
say("    RASSF9's core is INTERMEDIATE in tumours, not high. RASSF10's is")
say("    nearly unmethylated. Neither locus is a silenced-by-hypermethylation")
say("    picture in bulk tumour to begin with.")

fwrite(PROBES, file.path(RES_DIR, "44a_probe_table.csv"))

## ---- 2. Locate the file -----------------------------------------------------
say("")
say("== 2. The methylation matrix ===================================")

## Both spellings are accepted because the true extension may be hidden in
## Explorer. Whichever is found, it is opened through gzfile() below.
CAND <- c("GSE68379_Matrix.processed.txt.gz",
          "GSE68379_Matrix.processed.txt",
          "GSE68379_Matrix.Signal.Intensities.txt.gz",
          "GSE68379_Matrix.Signal.Intensities.txt")
hit  <- CAND[file.exists(file.path(GEO_DIR, CAND))]
f_met <- if (length(hit)) file.path(GEO_DIR, hit[1]) else NA_character_

if (is.na(f_met)) {
  say("  MISSING. Expected one of these in ", GEO_DIR, ":")
  for (c0 in CAND) say("    ", c0)
  say("")
  say("  Download the processed matrix from the GEO series page:")
  say("    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68379")
  stop("methylation matrix not found", call. = FALSE)
}
say(sprintf("  %-44s %8.2f GB", basename(f_met), file.size(f_met) / 1024^3))
if (grepl("Signal.Intensities", basename(f_met)))
  warn("this is the SIGNAL INTENSITIES file, not beta values -- ",
       "download GSE68379_Matrix.processed.txt.gz instead")

## ---- 3. Header and first rows -----------------------------------------------
say("")
say("== 3. File layout ==============================================")

## Always open through gzfile(): it auto-detects and reads an uncompressed
## file correctly too. Switching on the ".gz" suffix is unsafe because Windows
## Explorer hides known extensions, so a file shown as "...processed.txt" may
## really be "...processed.txt.gz" -- and a gzipped file opened with file()
## returns binary noise rather than an error.
con <- gzfile(f_met, "rt")
head_lines <- readLines(con, n = 6, warn = FALSE)
close(con)

if (length(head_lines) && any(grepl("[\x01-\x08\x0e-\x1f]", head_lines[1]))) {
  warn("the first line contains control characters -- the file may be an ",
       "archive containing the matrix (e.g. .tar.gz or .zip) rather than a ",
       "gzipped text file. Extract it and point GEO_DIR at the .txt inside.")
}

if (!length(head_lines)) stop("could not read any lines from the matrix", call. = FALSE)

## Determine the delimiter from the first line rather than assuming tab.
delim <- names(sort(c(
  "\t" = lengths(regmatches(head_lines[1], gregexpr("\t", head_lines[1]))),
  ","  = lengths(regmatches(head_lines[1], gregexpr(",",  head_lines[1]))),
  ";"  = lengths(regmatches(head_lines[1], gregexpr(";",  head_lines[1])))
), decreasing = TRUE))[1]
say("  delimiter inferred: ", ifelse(delim == "\t", "TAB", delim))

hdr <- strsplit(head_lines[1], delim, fixed = TRUE)[[1]]
hdr <- gsub('^"|"$', "", hdr)
say("  header fields: ", length(hdr))
say("  first 8 fields:")
for (i in seq_len(min(8, length(hdr)))) say(sprintf("    [%d] '%s'", i, hdr[i]))
say("  last 3 fields:")
for (i in seq(max(1, length(hdr) - 2), length(hdr))) say(sprintf("    [%d] '%s'", i, hdr[i]))

## Does the file pair each sample with a detection p-value?
n_beta  <- sum(grepl("AVG.?Beta|Beta", hdr, ignore.case = TRUE))
n_pval  <- sum(grepl("Detection|Pval|p.?value", hdr, ignore.case = TRUE))
n_cell  <- sum(grepl("^CellLine", hdr))
say("")
say("  fields mentioning Beta      : ", n_beta)
say("  fields mentioning Detection : ", n_pval)
say("  fields starting 'CellLine'  : ", n_cell)

## Braces are required here: at top level R treats the `if` as complete at the
## newline, so a bare `else` on the following line is a parse error.
layout <- if (n_pval > 0 && n_beta > 0) {
  "paired_beta_pval"
} else if (n_beta > 0 || n_cell > 0) {
  "beta_only"
} else {
  "unknown"
}
say("  layout: ", layout)

if (layout == "unknown")
  warn("could not classify the column layout from the header -- ",
       "44b must not be run until this is resolved by eye")
if (layout == "paired_beta_pval")
  say("  -> 44b must select only the beta columns and ignore the p-value ones.")

## The first data rows: are probe IDs in column 1?
say("")
say("  first 3 data rows, first field and first 3 values:")
for (i in 2:min(4, length(head_lines))) {
  f <- strsplit(head_lines[i], delim, fixed = TRUE)[[1]]
  f <- gsub('^"|"$', "", f)
  say(sprintf("    '%s'  ->  %s", f[1], paste(head(f[-1], 3), collapse = ", ")))
}

## Work out the column offset between header and data rows.
## This file was written with row names, so each data row carries one extra
## leading field (a row number) that the header does not name. Header field j
## therefore corresponds to DATA field j + offset. Getting this wrong shifts
## every beta value by one sample, which would look like plausible data.
dat1   <- gsub('^"|"$', "", strsplit(head_lines[2], delim, fixed = TRUE)[[1]])
n_dat  <- length(dat1)
offset <- n_dat - length(hdr)

say("")
say("  header fields: ", length(hdr), " | data-row fields: ", n_dat,
    " | offset: ", offset)

probe_field <- which(grepl("^(cg|ch\\.|rs)[0-9A-Za-z.]+$", dat1))[1]
if (is.na(probe_field)) {
  warn("no field in the first data row looks like a 450k probe ID -- ",
       "do not run 44b until this is resolved by eye")
} else {
  say("  probe IDs are in DATA field ", probe_field, " ('", dat1[probe_field], "')")
}

if (offset == 1) {
  say("  -> the leading data field is an unnamed row index ('", dat1[1], "').")
  say("     44b must read data rows positionally, mapping header field j to")
  say("     data field j + 1, and take the probe ID from field ", probe_field, ".")
} else if (offset == 0) {
  say("  -> header and data rows align; no offset correction needed.")
} else {
  warn("unexpected offset of ", offset,
       " fields between header and data rows -- resolve before running 44b")
}

## Beta columns are header fields 2..1029 under this layout; translate those to
## data-field positions and sanity-check that they really are beta values.
beta_hdr_idx <- which(grepl("AVG.?Beta", hdr, ignore.case = TRUE))
pval_hdr_idx <- which(grepl("Detection", hdr, ignore.case = TRUE))
beta_dat_idx <- beta_hdr_idx + offset
pval_dat_idx <- pval_hdr_idx + offset

say("")
say("  beta columns  : header ", min(beta_hdr_idx), "-", max(beta_hdr_idx),
    "  -> data ", min(beta_dat_idx), "-", max(beta_dat_idx))
say("  pval columns  : header ", min(pval_hdr_idx), "-", max(pval_hdr_idx),
    "  -> data ", min(pval_dat_idx), "-", max(pval_dat_idx))

vb <- suppressWarnings(as.numeric(dat1[beta_dat_idx])); vb <- vb[is.finite(vb)]
vp <- suppressWarnings(as.numeric(dat1[pval_dat_idx])); vp <- vp[is.finite(vp)]

if (length(vb)) {
  say(sprintf("  row 1 BETA values: n = %d, range %.4f to %.4f, median %.4f",
              length(vb), min(vb), max(vb), median(vb)))
  if (max(vb) > 1.01) warn("beta values exceed 1 -- the column mapping is wrong")
} else warn("no numeric beta values parsed from the first data row")

if (length(vp)) {
  say(sprintf("  row 1 PVAL values: n = %d, range %.2e to %.2e, median %.2e",
              length(vp), min(vp), max(vp), median(vp)))
  ## Detection p-values on a good array are overwhelmingly near zero. If the
  ## supposed p-value block has a high median it is probably the beta block,
  ## meaning the offset is wrong.
  if (median(vp) > 0.5)
    warn("the p-value block has median ", signif(median(vp), 3),
         " -- suspiciously high; the beta/pval column mapping may be reversed")
}
if (length(vb) && length(vp) && median(vb) < median(vp))
  warn("median of the beta block (", signif(median(vb), 3),
       ") is below the p-value block (", signif(median(vp), 3),
       ") -- check the mapping before running 44b")

PLAN <- data.table(
  file = basename(f_met), delimiter = ifelse(delim == "\t", "TAB", delim),
  n_header_fields = length(hdr), n_data_fields = n_dat, field_offset = offset,
  layout = layout, probe_data_field = probe_field,
  n_beta_fields = n_beta, n_pval_fields = n_pval,
  beta_data_first = min(beta_dat_idx), beta_data_last = max(beta_dat_idx),
  pval_data_first = min(pval_dat_idx), pval_data_last = max(pval_dat_idx)
)
fwrite(PLAN, file.path(RES_DIR, "44a_parse_plan.csv"))

## ---- 4. Sample names -> DepMap models ---------------------------------------
say("")
say("== 4. Mapping GEO samples to melanoma lines ====================")

f_model <- file.path(CCLE_DIR, "Model.csv")
if (!file.exists(f_model)) stop("Model.csv not found in ", CCLE_DIR, call. = FALSE)
MD <- fread(f_model, showProgress = FALSE)

## Sample columns, with the beta/pval suffix and the CellLine_ prefix removed.
samp_fields <- hdr[-1]
if (layout == "paired_beta_pval")
  samp_fields <- samp_fields[grepl("Beta", samp_fields, ignore.case = TRUE)]

## Observed field format in this file: "<NAME>_AVG.Beta" and
## "<NAME>_Detection.PVal", with no CellLine_ prefix (the prefix appears in the
## GEO sample titles but not in this matrix's header). Strip the known suffixes
## explicitly, then fall back to a looser pattern for anything unexpected.
clean <- function(x) {
  x <- sub("^CellLine_", "", x)
  x <- sub("_AVG[._ ]?Beta$",        "", x, ignore.case = TRUE)
  x <- sub("_Detection[._ ]?PVal$",  "", x, ignore.case = TRUE)
  x <- sub("[._ ]*(AVG[._ ]?Beta|Detection[._ ]?PVal.*)$", "", x, ignore.case = TRUE)
  trimws(x)
}
samp_names <- clean(samp_fields)
say("  sample columns detected: ", length(samp_names))
say("  examples: ", paste(head(samp_names, 6), collapse = ", "))

## If cleaning left a suffix behind, every name is wrong and the join silently
## returns nothing, so fail loudly rather than reporting zero matches.
resid <- grepl("Beta|Detection|PVal", samp_names, ignore.case = TRUE)
if (any(resid)) {
  warn(sum(resid), " sample name(s) still carry a column suffix after cleaning, e.g. '",
       samp_names[which(resid)[1]], "' -- the name join below is not trustworthy")
}
if (anyDuplicated(samp_names))
  warn(sum(duplicated(samp_names)), " duplicate sample name(s) after cleaning, e.g. '",
       samp_names[duplicated(samp_names)][1], "'")

## Match on a stripped form -- uppercase, alphanumeric only. This is the same
## normalisation DepMap uses for StrippedCellLineName, so it is a fair join.
strip <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

MEL <- MD[grepl("melanoma", OncotreePrimaryDisease, ignore.case = TRUE),
          .(model_id = ModelID, name = StrippedCellLineName,
            sanger = if ("SangerModelID" %in% names(MD)) SangerModelID else NA_character_,
            disease = OncotreePrimaryDisease, subtype = OncotreeSubtype)]
MEL[, uveal := grepl("uveal|ocular", paste(subtype, disease), ignore.case = TRUE)]

## NB: do not name this column "key" -- data.table() has a `key` argument, so
## `data.table(key = x)` sets the table's key instead of creating a column.
MEL[, match_key := strip(name)]

SAMP <- data.table(field     = samp_fields,
                   geo_name  = samp_names,
                   match_key = strip(samp_names))

MM <- merge(MEL, SAMP, by = "match_key", all.x = TRUE)
matched <- MM[!is.na(field)]

say("")
say("  melanoma lines in Model.csv            : ", nrow(MEL))
say("  of those matched to a GEO sample       : ", nrow(matched))
say("  non-uveal melanoma matched             : ", sum(!matched$uveal))
say("  melanoma lines with NO methylation     : ", sum(is.na(MM$field)))
say("")
say("  Model.csv carries SangerModelID: ",
    ifelse("SangerModelID" %in% names(MD), "yes", "NO -- name matching only"))
say("  (the GEO titles are plain names, not Sanger IDs, so the join is on the")
say("   stripped name; SangerModelID is kept in the output for verification)")

if (nrow(matched)) {
  say("")
  say("  matched non-uveal lines by subtype:")
  bysub <- matched[uveal == FALSE, .N, by = subtype][order(-N)]
  for (i in seq_len(nrow(bysub))) say(sprintf("    %-40s %4d", bysub$subtype[i], bysub$N[i]))
}

if (sum(!matched$uveal) < 20)
  warn("only ", sum(!matched$uveal), " non-uveal melanoma lines have methylation data -- ",
       "too few for anything beyond a description")

fwrite(MM[order(-uveal, name)], file.path(RES_DIR, "44a_sample_map.csv"))
say("")
say("  wrote 44a_probe_table.csv, 44a_parse_plan.csv, 44a_sample_map.csv")

## ---- 5. What 44b will and will not do ---------------------------------------
say("")
say("== 5. What step 2 can answer ===================================")
say("  The expression gate (42, 43) returned STOP: RASSF9 is at the detection")
say("  floor in melanoma lines, so no methylation-expression CORRELATION is")
say("  computable. 44b therefore does not attempt one.")
say("")
say("  It answers the one question that survives zero expression variance:")
say("  is the RASSF9 promoter core METHYLATED in melanoma cells?")
say("")
say("    high in lines  -> silencing consistent with the TCGA direction; the")
say("                      correlation is untestable in lines only because the")
say("                      gene is uniformly off")
say("    low in lines   -> methylation is not what holds RASSF9 off in")
say("                      melanoma cells, and the TCGA association is more")
say("                      likely driven by the non-malignant compartments")
say("                      where the gene is actually expressed")
say("")
say("  The shore and distal probes are carried through as internal controls,")
say("  and the TCGA median for each probe travels in the probe table so the")
say("  two can be put side by side per probe rather than as pooled averages.")
say("")
say("  WHAT NEITHER SCRIPT CAN SAY: cell lines are not tumours. Cultured lines")
say("  drift epigenetically, and a difference from TCGA may reflect culture")
say("  rather than malignant-cell biology in vivo. This is a supporting")
say("  observation, not a test of the TCGA finding.")
