#!/usr/bin/env Rscript

suppressMessages(library(data.table))
setwd("/mnt/d/ngs_rebuild")
dir.create("results/gene_lists", recursive = TRUE, showWarnings = FALSE)

ann_f <- Sys.glob("nfcore/results/bwa/merged_library/macs2/narrow_peak/consensus/*annotatePeaks.txt")
if (!length(ann_f)) stop("annotatePeaks.txt not found")
ann  <- fread(ann_f[1])
gcol <- grep("^Gene.?Name$", names(ann), ignore.case = TRUE, value = TRUE)[1]
if (is.na(gcol)) stop("no 'Gene Name' column found in ", ann_f[1])
map <- setNames(ann[[gcol]], ann[[1]])

dump <- function(bed, out) {
  ids <- fread(bed, header = FALSE)$V4
  g   <- sort(unique(na.omit(map[ids])))
  g   <- g[g != "" & !is.na(g)]
  writeLines(g, out)
  cat(sprintf("  %-22s %6d genes\n", basename(out), length(g)))
}
dump("results/diff_beds/E21_up.bed",     "results/gene_lists/E21_up_genes.txt")
dump("results/diff_beds/P11_up.bed",     "results/gene_lists/P11_up_genes.txt")
dump("results/diff_beds/all_tested.bed", "results/gene_lists/background_genes.txt")
cat("\nPaste E21_up / P11_up into g:Profiler, organism Rattus norvegicus,\n")
cat("with background_genes.txt as the custom background. Export as CSV.\n")
