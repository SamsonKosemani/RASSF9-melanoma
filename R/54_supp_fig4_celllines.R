## =============================================================================
## 54_supp_fig4_celllines.R
##
## Supplementary Figure 4 -- RASSF9 expression across cancer cell lines.
##
## THE CLAIM THIS FIGURE SUPPORTS (Results, cell-line paragraph):
##   "Basal RASSF9 expression was below 1 TPM in 109 of 113 non-uveal melanoma
##    lines, including all 15 explicitly annotated cutaneous melanoma lines."
##
##   That sentence currently cites a figure that was never drawn. The files in
##   figures/ccle/ are working diagnostics from 23 October at screen size, not
##   artwork. This draws the figure natively at 180 mm, as scripts 50, 52 and
##   53 do.
##
## PANELS
##   A  Every non-uveal melanoma line, ranked by RASSF9 expression, against the
##      1 TPM threshold. Shows the floor directly rather than asserting it.
##   B  The same lines split by annotated subtype, so the "all 15 cutaneous"
##      part of the claim is visible rather than taken on trust.
##
## ROBUSTNESS
##   The column names in the CCLE tables are not known to me, so section 1
##   prints an inventory of every candidate file and its columns BEFORE
##   anything is drawn, and every column used is resolved by pattern with a
##   loud failure if it is absent. A run either produces the figure or says
##   exactly what it could not find. Guessing column names is how scripts 44
##   and 53 broke.
##
## OUTPUT
##   figures/ccle/SupplementaryFigure4.{pdf,png}
##   results/logs/54_supp_fig4_celllines.log   (send this back)
## =============================================================================

if (file.exists("R/00_setup.R")) source("R/00_setup.R")
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

ROOT    <- if (exists("PROJECT_ROOT")) PROJECT_ROOT else getwd()
LOG_DIR <- file.path(ROOT, "results", "logs")
FIG_DIR <- file.path(ROOT, "figures", "ccle")
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_PATH <- file.path(LOG_DIR, "54_supp_fig4_celllines.log")
.logcon  <- file(LOG_PATH, open = "wt")
WARN <- character(0)
say  <- function(...) { m <- paste0(...); cat(m, "\n", sep = ""); cat(m, "\n", sep = "", file = .logcon) }
warn <- function(...) { m <- paste0(...); WARN <<- c(WARN, m); say("  [WARN] ", m) }
on.exit({
  say("")
  if (length(WARN)) { say("WARNINGS (", length(WARN), "):"); for (w in WARN) say("  - ", w) }
  say("Full log: ", LOG_PATH); close(.logcon)
}, add = TRUE)

## ---- journal constants, identical to scripts 50/52/53 -----------------------
WIDTH_MM <- 180; HEIGHT_MM <- 78
BASE_PT <- 7; AXIS_PT <- 6; LAB_PT <- 5.5; TAG_PT <- 9
COL_MEL <- "#C0392B"; COL_NEV <- "#1F6FB2"; COL_LO <- "#9A9A9A"; COL_INK <- "#2A2A2A"

theme_fig <- function() {
  theme_classic(base_size = BASE_PT) +
    theme(axis.text   = element_text(size = AXIS_PT, colour = COL_INK),
          axis.title  = element_text(size = BASE_PT, colour = COL_INK),
          axis.line   = element_line(linewidth = 0.3, colour = COL_INK),
          axis.ticks  = element_line(linewidth = 0.3, colour = COL_INK),
          legend.text = element_text(size = LAB_PT),
          legend.title = element_blank(),
          legend.key.size = unit(3, "mm"),
          plot.margin = margin(4, 6, 4, 4, unit = "pt"))
}

say("")
say("================================================================")
say(" Supplementary Figure 4 -- RASSF9 across cancer cell lines")
say("================================================================")

## ---- 1. Inventory -----------------------------------------------------------
say(""); say("== 1. Inventory of candidate tables ============================")
CC <- file.path(ROOT, "results", "ccle")
if (!dir.exists(CC)) stop("results/ccle not found under ", ROOT)
cands <- list.files(CC, pattern = "\\.csv$", full.names = TRUE)
say("  ", length(cands), " csv file(s) in results/ccle")
for (f in cands) {
  h <- tryCatch(names(fread(f, nrows = 1L, showProgress = FALSE)),
                error = function(e) NULL)
  n <- tryCatch(nrow(fread(f, select = 1L, showProgress = FALSE)),
                error = function(e) NA_integer_)
  say(sprintf("    %-42s %5s rows | %s", basename(f),
              ifelse(is.na(n), "?", n),
              ifelse(is.null(h), "unreadable", paste(h, collapse = ", "))))
}

find1 <- function(pat) {
  hit <- cands[grepl(pat, basename(cands), ignore.case = TRUE)]
  if (!length(hit)) NA_character_ else hit[1]
}
pick <- function(DT, pats, what, where) {
  for (p in pats) {
    m <- grep(p, names(DT), ignore.case = TRUE, value = TRUE)
    if (length(m)) return(m[1])
  }
  stop("could not find a column for ", what, " in ", basename(where),
       "\n  columns present: ", paste(names(DT), collapse = ", "),
       "\n  patterns tried:  ", paste(pats, collapse = ", "))
}

## Per-line melanoma table is the backbone. 46_line_level_data is the fallback
## because it carries lineage and expression for every line, melanoma included.
f_mel  <- find1("^42_rassf9_melanoma_lines")
f_line <- find1("^46_line_level_data")
f_sub  <- find1("^43_subtype_breakdown_rassf9")
say("")
say("  melanoma per-line table: ", ifelse(is.na(f_mel), "NOT FOUND", basename(f_mel)))
say("  line-level table:        ", ifelse(is.na(f_line), "NOT FOUND", basename(f_line)))
say("  subtype breakdown:       ", ifelse(is.na(f_sub), "NOT FOUND", basename(f_sub)))
if (is.na(f_mel) && is.na(f_line))
  stop("neither 42_rassf9_melanoma_lines.csv nor 46_line_level_data.csv is present; ",
       "nothing to draw from")

## ---- 2. Build the melanoma line table ---------------------------------------
say(""); say("== 2. Melanoma lines ===========================================")
src <- if (!is.na(f_mel)) f_mel else f_line
D   <- fread(src, showProgress = FALSE)
say("  using ", basename(src), " (", nrow(D), " rows)")

c_expr <- pick(D, c("^rassf9$", "rassf9.*(tpm|expr|log2)", "(tpm|expr|log2).*rassf9",
                    "^expr", "^tpm"), "RASSF9 expression", src)
c_name <- pick(D, c("^(cell_line_)?name$", "stripped", "^line$", "^model", "^cell"),
               "cell-line name", src)
say("  expression column: '", c_expr, "' | name column: '", c_name, "'")

## If the fallback table was used it holds every lineage, so melanoma must be
## selected explicitly rather than assumed.
if (is.na(f_mel)) {
  c_lin <- pick(D, c("^lineage$", "^disease$", "oncotree", "primary_disease"),
                "lineage", src)
  before <- nrow(D)
  D <- D[grepl("melanoma|skin", get(c_lin), ignore.case = TRUE) &
           !grepl("uveal|ocular|conjunctiv", get(c_lin), ignore.case = TRUE)]
  say("  filtered to non-uveal melanoma on '", c_lin, "': ", before, " -> ", nrow(D))
}

A <- D[, .(line = as.character(get(c_name)), expr = as.numeric(get(c_expr)))]
A <- A[is.finite(expr)]
say("  lines with finite expression: ", nrow(A))
if (nrow(A) < 20) stop("only ", nrow(A), " lines -- this is not the melanoma cohort")

## ---- 3. Which scale are we on? ----------------------------------------------
## The manuscript's threshold is 1 TPM. On a log2(TPM+1) axis that is log2(2)=1.
## Getting this wrong would move the threshold line and misstate the count, so
## it is decided from the data rather than assumed.
is_log2 <- grepl("log2", c_expr, ignore.case = TRUE) || max(A$expr, na.rm = TRUE) < 15
THRESH  <- if (is_log2) 1 else 1
say("")
say("== 3. Scale ====================================================")
say("  expression range: ", sprintf("%.3f to %.3f", min(A$expr), max(A$expr)))
say("  treated as ", ifelse(is_log2, "log2(TPM + 1)", "TPM"),
    " | threshold drawn at ", THRESH,
    ifelse(is_log2, "  (= 1 TPM on the log2(TPM+1) scale)", " TPM"))
if (!is_log2)
  warn("values look like raw TPM, not log2(TPM+1); check the axis label is right")

n_below <- sum(A$expr < THRESH); n_tot <- nrow(A)
say(sprintf("  below threshold: %d of %d (%.1f%%)", n_below, n_tot,
            100 * n_below / n_tot))
say("  the manuscript states 109 of 113")
if (n_tot != 113)
  warn("cohort is ", n_tot, " lines, the manuscript says 113 -- reconcile before use")
if (n_below != 109)
  warn("count below 1 TPM is ", n_below, ", the manuscript says 109 -- reconcile before use")

## ---- 4. Panel A : every line, ranked ----------------------------------------
A <- A[order(expr)]
A[, rank := seq_len(.N)]
A[, lo := expr < THRESH]

pA <- ggplot(A, aes(rank, expr, colour = lo)) +
  geom_hline(yintercept = THRESH, linetype = "dashed", linewidth = 0.3,
             colour = "grey45") +
  geom_point(size = 0.7) +
  scale_colour_manual(values = c(`TRUE` = COL_MEL, `FALSE` = COL_INK), guide = "none") +
  annotate("text", x = 1, y = THRESH, vjust = -0.7, hjust = 0, size = LAB_PT / .pt,
           colour = "grey35", label = "1 TPM") +
  annotate("text", x = n_tot, y = max(A$expr), hjust = 1, vjust = 1,
           size = LAB_PT / .pt, colour = COL_INK,
           label = sprintf("%d of %d lines below 1 TPM", n_below, n_tot)) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.12))) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.03))) +
  labs(x = "Melanoma cell lines, ranked",
       y = if (is_log2) expression(italic("RASSF9")*", log"[2]*"(TPM + 1)")
           else expression(italic("RASSF9")*" (TPM)")) +
  theme_fig()

## ---- 5. Panel B : by annotated subtype --------------------------------------
say(""); say("== 4. Panel B : subtype ========================================")
c_sub <- tryCatch(pick(D, c("^subtype$", "oncotree.*subtype", "subtype", "^disease_sub",
                            "^lineage_sub"), "subtype", src),
                  error = function(e) NA_character_)
pB <- NULL
if (!is.na(c_sub)) {
  B <- D[, .(line = as.character(get(c_name)), expr = as.numeric(get(c_expr)),
             sub = as.character(get(c_sub)))]
  B <- B[is.finite(expr)]
  B[is.na(sub) | !nzchar(trimws(sub)), sub := "Not annotated"]
  B[, lo := expr < THRESH]
  ## "Mucosal Melanoma of the Vulva/Vagina" is 36 characters and would eat a
  ## third of the panel. The site clause adds nothing the subtype name does
  ## not already carry, so it is dropped.
  B[, sub := sub("\\s+of the\\s+.*$", "", sub)]
  ## Keep every subtype as its own group. Folding "Melanoma NOS" into cutaneous
  ## was explicitly rejected in script 43 and the same rule holds here.
  cnt <- B[, .(n = .N, n_lo = sum(lo), med = median(expr)), by = sub][order(-n)]
  say("  subtypes found on '", c_sub, "':")
  for (i in seq_len(nrow(cnt)))
    say(sprintf("    %-34s %3d lines, %3d below 1 TPM, median %.3f",
                cnt$sub[i], cnt$n[i], cnt$n_lo[i], cnt$med[i]))
  cut_row <- cnt[grepl("cutane", sub, ignore.case = TRUE)]
  if (nrow(cut_row))
    say(sprintf("  cutaneous: %d of %d below 1 TPM (manuscript says all 15)",
                cut_row$n_lo[1], cut_row$n[1]))
  else
    warn("no subtype matching 'cutaneous' -- the manuscript's 15-line claim ",
         "cannot be shown in panel B from this column")

  ## Carry n in the label: with four rows there is room, and it puts the
  ## "all 15 cutaneous" claim on the axis instead of only in the text.
  cnt[, lab := sprintf("%s (n = %d)", sub, n)]
  lab_of <- setNames(cnt$lab, cnt$sub)
  B[, lab := lab_of[sub]]
  B[, lab := factor(lab, levels = unname(lab_of[cnt[order(med)]$sub]))]
  pB <- ggplot(B, aes(expr, lab)) +
    geom_vline(xintercept = THRESH, linetype = "dashed", linewidth = 0.3,
               colour = "grey45") +
    geom_jitter(aes(colour = lo), height = 0.17, width = 0, size = 0.7, alpha = 0.9) +
    stat_summary(fun = median, geom = "point", shape = 124, size = 2.4,
                 colour = COL_INK) +
    scale_colour_manual(values = c(`TRUE` = COL_MEL, `FALSE` = COL_INK), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.10))) +
    labs(x = if (is_log2) expression(italic("RASSF9")*", log"[2]*"(TPM + 1)")
             else expression(italic("RASSF9")*" (TPM)"),
         y = NULL) +
    theme_fig()
} else {
  warn("no subtype column in ", basename(src),
       " -- panel B cannot be drawn; the figure will be panel A alone")
}

## ---- 6. Assemble -------------------------------------------------------------
say(""); say("== 5. Assembling ===============================================")
FIG <- if (is.null(pB)) pA else
  (pA | pB) + plot_layout(widths = c(1, 1.1))
FIG <- FIG + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = TAG_PT, face = "bold", hjust = 0, vjust = 1))

save_fig <- function(path, ...) {
  tmp <- file.path(dirname(path), paste0(".tmp_", basename(path)))
  ok <- tryCatch({ ggsave(tmp, FIG, ...); TRUE }, error = function(e) { warn(conditionMessage(e)); FALSE })
  if (!ok) return(invisible(FALSE))
  if (file.exists(path)) try(file.remove(path), silent = TRUE)
  if (!file.rename(tmp, path)) { warn("could not replace ", basename(path),
                                      " -- is it open in a viewer?"); return(invisible(FALSE)) }
  say("  wrote ", basename(path)); invisible(TRUE)
}
save_fig(file.path(FIG_DIR, "SupplementaryFigure4.pdf"),
         width = WIDTH_MM, height = HEIGHT_MM, units = "mm", device = cairo_pdf)
save_fig(file.path(FIG_DIR, "SupplementaryFigure4.png"),
         width = WIDTH_MM, height = HEIGHT_MM, units = "mm", dpi = 600)

## ---- 7. Layout report --------------------------------------------------------
## Text in place of an image, as in script 53, so the figure can be checked
## from a log rather than by uploading a PNG.
say(""); say("== 6. Layout report ============================================")
report <- function() {
  say("")
  say("  6.1 Type sizes (pt), floor is 5.5")
  tp <- c(base = BASE_PT, axis = AXIS_PT, annotation = LAB_PT, tag = TAG_PT)
  for (n in names(tp))
    say(sprintf("    %-12s %4.1f  %s", n, tp[[n]],
                ifelse(tp[[n]] >= 5.5, "ok", "TOO SMALL")))
  if (any(tp < 5.5)) warn("a type size is below the 5.5 pt floor")

  say("")
  say("  6.2 Panel A data against the axis")
  say(sprintf("    expression %.3f to %.3f, threshold %.2f", min(A$expr), max(A$expr), THRESH))
  if (THRESH < min(A$expr) || THRESH > max(A$expr))
    warn("the 1 TPM line falls outside the data range, so it will not be visible")

  if (!is.null(pB)) {
    say("")
    say("  6.3 Panel B rows")
    h <- HEIGHT_MM * 0.78 / nrow(cnt)
    say(sprintf("    %d subtypes over ~%.0f mm -> %.1f mm per row  %s",
                nrow(cnt), HEIGHT_MM * 0.78, h,
                ifelse(h > 18, "AIRY", ifelse(h < 3.5, "CRAMPED", "ok"))))
    if (h > 18 || h < 3.5) warn("panel B row spacing is ", round(h, 1), " mm")
    w <- max(nchar(as.character(cnt$lab)))
    say("    longest subtype label: ", w, " characters")
    if (w > 28) warn("a subtype label is ", w, " characters and may crowd the panel")
  }

  say("")
  say("  6.4 Colour roles")
  say(sprintf("    %s red   below 1 TPM", COL_MEL))
  say(sprintf("    %s ink   at or above 1 TPM", COL_INK))
  say("    Red means 'below the detection threshold' in both panels, the same")
  say("    meaning in each, and it is not used for anything else here.")

  ## ASCII view of the rendered PNG, same method as script 53.
  say("")
  say("  6.5 Low-resolution view of the rendered figure")
  if (!requireNamespace("png", quietly = TRUE))
    try(utils::install.packages("png", repos = "https://cloud.r-project.org"), silent = TRUE)
  if (!requireNamespace("png", quietly = TRUE)) {
    say("    the 'png' package is unavailable, so the view is skipped"); return(invisible())
  }
  tf <- tempfile(fileext = ".png")
  ok <- tryCatch({ ggsave(tf, FIG, width = WIDTH_MM, height = HEIGHT_MM,
                          units = "mm", dpi = 100); TRUE }, error = function(e) FALSE)
  if (!ok || !file.exists(tf)) { say("    preview could not be rendered"); return(invisible()) }
  im <- try(png::readPNG(tf), silent = TRUE); unlink(tf)
  if (inherits(im, "try-error")) { say("    preview could not be read"); return(invisible()) }
  d <- dim(im)
  g <- if (length(d) == 2L) im else if (length(d) == 3L && d[3] >= 3L) {
         lum <- 0.299 * im[, , 1] + 0.587 * im[, , 2] + 0.114 * im[, , 3]
         if (d[3] >= 4L) lum * im[, , 4] + (1 - im[, , 4]) else lum
       } else if (length(d) == 3L) im[, , 1] else NULL
  if (is.null(g) || !is.matrix(g <- as.matrix(g))) {
    say("    the render had an unexpected shape, so the view is skipped"); return(invisible())
  }
  nr <- nrow(g); nc <- ncol(g); cols <- 104L
  rows <- max(10L, round(cols * (HEIGHT_MM / WIDTH_MM) / 2.05))
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
  ## Clamp by assignment: pmin/pmax take attributes from their FIRST argument
  ## and would silently drop the matrix dim, breaking the row loop below.
  idx <- 1L + floor(blk * length(ramp))
  idx[idx < 1L] <- 1L; idx[idx > length(ramp)] <- length(ramp)
  if (!is.matrix(idx)) dim(idx) <- dim(blk)
  say(sprintf("    %d x %d characters from the %d x %d px render", cols, rows, nc, nr))
  say("    |", strrep("-", cols), "|")
  for (i in seq_len(rows)) say("    |", paste0(ramp[idx[i, ]], collapse = ""), "|")
  say("    |", strrep("-", cols), "|")
  say(sprintf("    ink coverage %.1f%%", 100 * mean(blk < 0.92)))
  invisible()
}
tryCatch(report(), error = function(e)
  say("  [layout report failed: ", conditionMessage(e), " -- the figure itself is fine]"))

say("")
say("  Drawn at ", WIDTH_MM, " mm, the final width. Place at 100%.")
