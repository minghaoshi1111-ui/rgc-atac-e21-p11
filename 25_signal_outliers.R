#!/usr/bin/env Rscript

# 25_signal_outliers.R
#This is a SENSITIVITY CHECK ONLY

suppressMessages({library(data.table); library(GenomicRanges)})

QT <- 0.999   

bc <- fread("counts/bin_counts.tsv")
cat("columns found:", paste(names(bc), collapse=" | "), "\n")

coordn <- names(bc)[1:3]                                  # chrom / start / end
samp   <- setdiff(names(bc), coordn)
samp   <- samp[vapply(bc[, ..samp], is.numeric, logical(1))]
cat("coordinates:", paste(coordn, collapse=","), "\n")
cat("samples    :", paste(samp,   collapse=","), "\n")
stopifnot(length(samp) >= 2)

M   <- as.matrix(bc[, ..samp])
cpm <- t(t(M) / colSums(M)) * 1e6       
thr <- apply(cpm, 2, function(v) quantile(v[v > 0], QT))
cat("\nper-sample 99.9th percentile CPM (non-zero bins):\n"); print(round(thr, 3))

flag <- rowSums(cpm > matrix(thr, nrow(cpm), ncol(cpm), byrow = TRUE)) == ncol(cpm)
cat(sprintf("\n[B1] bins total %d ; extreme in ALL %d samples: %d (%.4f%%) = %.2f Mb\n",
            nrow(bc), ncol(cpm), sum(flag), 100*mean(flag), sum(flag)*10000/1e6))

top <- head(order(-rowMeans(cpm)*flag), 12)
cat("\ntop flagged bins by mean CPM:\n")
print(data.table(chrom=bc[[1]][top], start=bc[[2]][top], end=bc[[3]][top],
                 meanCPM=round(rowMeans(cpm)[top],1), flagged=flag[top]))

bl <- reduce(GRanges(bc[[1]][flag], IRanges(bc[[2]][flag] + 1, bc[[3]][flag])))
fwrite(data.table(as.character(seqnames(bl)), start(bl)-1, end(bl)),
       "qc/signal_outlier_bins.bed", sep="\t", col.names=FALSE)
cat(sprintf("\nwrote qc/signal_outlier_bins.bed : %d merged intervals\n", length(bl)))

readbed <- function(f){ d <- fread(f, header=FALSE); GRanges(d$V1, IRanges(d$V2+1, d$V3)) }
chk <- function(f, tag){
  g <- readbed(f); n <- sum(overlapsAny(g, bl))
  cat(sprintf("%-5s %-32s %7d regions | %5d flagged (%.3f%%)\n",
              tag, basename(f), length(g), n, 100*n/length(g)))
  g[overlapsAny(g, bl)]
}
chk("peaks/consensus.named.bed",  "[B2]")
chk("results/robust_FDR05.bed",   "[B3]")
up <- chk("results/composite_up.bed", "[B4]")
dn <- chk("results/composite_dn.bed", "[B4]")

if (length(up) || length(dn)) {
  cat("\nflagged composite regions (these need individual inspection):\n")
  print(c(as.character(up), as.character(dn)))
} else cat("\nno composite region touches a signal-outlier bin\n")
