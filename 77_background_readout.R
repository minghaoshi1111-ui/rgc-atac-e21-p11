#!/usr/bin/env Rscript

#

suppressMessages(library(SummarizedExperiment))

Z <- readRDS("/mnt/d/ngs_rebuild/results/monalisa_rat.rds")            
G <- readRDS("/mnt/d/ngs_rebuild/results/monalisa_rat_genomebg.rds")   
stopifnot(identical(rowData(Z)$motif.name, rowData(G)$motif.name))
nmv <- rowData(Z)$motif.name

FAM <- list(ONECUT = c("ONECUT1","ONECUT2","ONECUT3","CUX1","CUX2"),
            ISL    = c("ISL1","ISL2","Isl1"),
            MEF2   = c("MEF2A","MEF2B","MEF2C","MEF2D"),
            NFI    = c("NFIA","NFIB","NFIC","NFIX","NFIC::TLX1"),
            ROR    = c("RORA","RORB","RORC"),
            CTCF   = c("CTCF","CTCFL"),
            CREB1  = c("CREB1"))

bin_bounds <- function(nm) {
  v  <- strsplit(gsub("[][()]", "", nm), ",", fixed = TRUE)
  lo <- suppressWarnings(as.numeric(vapply(v, `[`, "", 1)))
  hi <- suppressWarnings(as.numeric(vapply(v, `[`, "", 2)))
  list(lo = lo, hi = hi)
}
changing <- function(x) {
  b <- bin_bounds(colnames(x))
  drop <- which(!is.na(b$lo) & !is.na(b$hi) & b$lo < 0 & b$hi > 0)
  keep <- setdiff(seq_len(ncol(x)), drop)
  attr(keep, "dropped") <- colnames(x)[drop]
  keep
}

kZ <- changing(Z); kG <- changing(G)
cat(sprintf("bins: %d total | zeroBin run keeps %d (dropped: %s)\n",
            ncol(Z), length(kZ), paste(attr(kZ,"dropped"), collapse=", ")))
cat(sprintf("                | genome  run keeps %d (dropped: %s)\n\n",
            length(kG), paste(attr(kG,"dropped"), collapse=", ")))

best <- function(x, keep, idx) {
  P <- assay(x, "negLog10Padj")[idx, keep, drop = FALSE]
  E <- assay(x, "log2enr")[idx, keep, drop = FALSE]
  k <- which(P == max(P, na.rm = TRUE), arr.ind = TRUE)[1, ]
  list(motif = nmv[idx][k[1]], bin = colnames(P)[k[2]],
       nlp = P[k[1], k[2]], enr = E[k[1], k[2]])
}

cat(sprintf("%-8s %-11s %-16s %10s %9s | %-11s %-16s %10s %9s\n",
            "family","motif","zeroBin bin","-log10padj","log2enr",
            "motif","genome bin","-log10padj","log2enr"))
rows <- list()
for (f in names(FAM)) {
  idx <- which(toupper(nmv) %in% toupper(FAM[[f]]))
  if (!length(idx)) { cat(sprintf("%-8s  no member in the PWM set\n", f)); next }
  z <- best(Z, kZ, idx); g <- best(G, kG, idx)
  cat(sprintf("%-8s %-11s %-16s %10.2f %+9.2f | %-11s %-16s %10.2f %+9.2f\n",
              f, z$motif, z$bin, z$nlp, z$enr, g$motif, g$bin, g$nlp, g$enr))
  rows[[f]] <- data.frame(family=f,
    motif_zero=z$motif, bin_zero=z$bin, nlp_zero=z$nlp, enr_zero=z$enr,
    motif_gen =g$motif, bin_gen =g$bin, nlp_gen =g$nlp, enr_gen =g$enr)
}
tab <- do.call(rbind, rows)
write.table(tab, "/mnt/d/ngs_rebuild/results/background_comparison.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

