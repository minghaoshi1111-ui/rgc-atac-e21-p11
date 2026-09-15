#!/usr/bin/env Rscript

suppressMessages({library(DESeq2); library(data.table)})
ROOT <- "/mnt/d/ngs_rebuild"; setwd(ROOT)
dir.create("results/diff_beds", recursive = TRUE, showWarnings = FALSE)

FC <- "nfcore/results/bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.featureCounts.txt"
fc  <- fread(FC, skip = 1)
cnt <- as.matrix(fc[, 7:ncol(fc)]); rownames(cnt) <- fc[[1]]
colnames(cnt) <- sub("\\.mLb.*$", "", basename(colnames(cnt)))
grp <- factor(sub("_.*$", "", colnames(cnt)), levels = c("P11","E21"))
dds <- DESeq(DESeqDataSetFromMatrix(cnt, DataFrame(grp), ~ grp))
res <- results(dds, contrast = c("grp","E21","P11"))

coord <- data.table(
  chr   = sub(";.*", "", fc$Chr),
  start = as.integer(sub(";.*", "", fc$Start)) - 1L,   # BED is 0-based
  end   = as.integer(sub(".*;", "", fc$End)),
  id    = fc[[1]])
setkey(coord, id)
coord <- coord[rownames(res)]
coord[, `:=`(lfc = res$log2FoldChange, padj = res$padj)]
coord <- coord[!is.na(padj)]

wr <- function(dt, f) {
  fwrite(dt[order(chr, start), .(chr, start, end, id)], f,
         sep = "\t", col.names = FALSE, quote = FALSE)
  cat(sprintf("  %-28s %7d regions\n", basename(f), nrow(dt)))
}
cat("padj < 0.05 and |log2FC| >= 1\n")
wr(coord[padj < 0.05 & lfc >=  1], "results/diff_beds/E21_up.bed")
wr(coord[padj < 0.05 & lfc <= -1], "results/diff_beds/P11_up.bed")
wr(coord[padj > 0.5 & abs(lfc) < 0.25], "results/diff_beds/static.bed")
wr(coord,                              "results/diff_beds/all_tested.bed")
