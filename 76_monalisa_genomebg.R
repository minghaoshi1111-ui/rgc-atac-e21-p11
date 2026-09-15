#!/usr/bin/env Rscript
suppressMessages({library(monaLisa); library(JASPAR2024); library(TFBSTools)
                  library(DESeq2); library(GenomicRanges); library(Biostrings)
                  library(Rsamtools); library(SummarizedExperiment); library(data.table); library(BiocParallel)})
dir.create("/mnt/d/ngs_rebuild/figures", showWarnings = FALSE)

FC <- "/mnt/d/ngs_rebuild/nfcore/results/bwa/merged_library/macs2/narrow_peak/consensus/consensus_peaks.mLb.clN.featureCounts.txt"
fc  <- fread(FC, skip = 1)
cnt <- as.matrix(fc[, 7:ncol(fc)]); rownames(cnt) <- fc[[1]]
colnames(cnt) <- sub("\\.mLb.*$", "", basename(colnames(cnt)))
stopifnot(ncol(cnt) == 4)
grp <- factor(sub("_.*$", "", colnames(cnt)), levels = c("P11","E21"))
cat("sample -> group:\n"); print(setNames(as.character(grp), colnames(cnt)))
stopifnot(all(table(grp) == 2))

dds <- DESeqDataSetFromMatrix(cnt, DataFrame(grp), ~ grp)
dds <- DESeq(dds)
res <- results(dds, contrast = c("grp","E21","P11"))
cat(sprintf("regions %d | padj<0.05 %d | |log2FC|>1 %d\n",
            nrow(res), sum(res$padj < .05, na.rm=TRUE),
            sum(abs(res$log2FoldChange) > 1, na.rm=TRUE)))

gr <- GRanges(sub(";.*","",fc$Chr),
              IRanges(as.integer(sub(";.*","",fc$Start)),
                      as.integer(sub(".*;","",fc$End))))
names(gr) <- fc[[1]]
gr <- gr[rownames(res)]
gr$lfc <- res$log2FoldChange
gr <- gr[!is.na(gr$lfc)]

fa <- FaFile("/mnt/d/ngs_rebuild/reference/genome.fa")
si <- seqinfo(fa)
seqlevels(gr, pruning.mode = "coarse") <- seqlevels(si)
seqinfo(gr) <- si[seqlevels(gr)]

n0 <- length(gr)
gr <- trim(resize(gr, width = 300, fix = "center"))
gr <- gr[width(gr) == 300]
cat(sprintf("regions after resize/trim: %d (dropped %d at sequence ends)\n", length(gr), n0 - length(gr)))
seqs <- getSeq(fa, gr)
names(seqs) <- names(gr)
stopifnot(length(seqs) == length(gr), all(width(seqs) == 300),
          !anyDuplicated(names(seqs)))
cat(sprintf("sequences: %d x 300 bp\n", length(seqs)))

bins <- bin(x = gr$lfc, binmode = "equalN", nElement = 5000, minAbsX = 1.0)
cat("\nbins:\n"); print(table(bins))
cat("zero bin index:", attr(bins, "bin0"), "\n\n")

# ---- motifs -----------------------------------------------------------------
pwms <- getMatrixSet(db(JASPAR2024()),
                     list(collection = "CORE", tax_group = "vertebrates",
                          matrixtype = "PWM"))
stopifnot(length(pwms) == 879)

gsq <- readDNAStringSet("/mnt/d/ngs_rebuild/reference/genome.fa")
names(gsq) <- sub(" .*", "", names(gsq))
cat(sprintf("genome loaded: %d sequences, %.2f Gb\n", length(gsq), sum(width(gsq))/1e9))

se <- calcBinnedMotifEnrR(seqs = seqs, bins = bins, pwmL = pwms,
                          background = "genome",
                          genome = gsq,
                          genome.oversample = 2,
                          BPPARAM = SerialParam(), verbose = TRUE)
saveRDS(se, "/mnt/d/ngs_rebuild/results/monalisa_rat_genomebg.rds")
cat("wrote results/monalisa_rat_genomebg.rds\n")

z <- readRDS("/mnt/d/ngs_rebuild/results/monalisa_rat.rds")   # zeroBin background
stopifnot(identical(rowData(se)$motif.name, rowData(z)$motif.name))

FAM <- list(ONECUT = c("ONECUT1","ONECUT2","ONECUT3"),
            ISL    = c("ISL1","ISL2","Isl1"),
            MEF2   = c("MEF2A","MEF2B","MEF2C","MEF2D"),
            NFI    = c("NFIA","NFIB","NFIC","NFIX","NFIC::TLX1"),
            CTCF   = c("CTCF","CTCFL"),
            CREB1  = c("CREB1"))

nmv <- rowData(se)$motif.name
best <- function(x, idx) {           # strongest bin for this family, that background
  P <- assay(x, "negLog10Padj")[idx, , drop = FALSE]
  E <- assay(x, "log2enr")[idx, , drop = FALSE]
  k <- which(abs(E) == max(abs(E), na.rm = TRUE), arr.ind = TRUE)[1, ]
  list(motif = nmv[idx][k[1]], bin = colnames(P)[k[2]],
       nlp = P[k[1], k[2]], enr = E[k[1], k[2]])
}

cat(sprintf("\n%-8s %-10s %-16s %10s %10s | %-16s %10s %10s\n",
            "family","motif","zeroBin bin","-log10padj","log2enr",
            "genome bin","-log10padj","log2enr"))
out <- list()
for (f in names(FAM)) {
  idx <- which(toupper(nmv) %in% toupper(FAM[[f]]))
  if (!length(idx)) { cat(sprintf("%-8s  no member in the PWM set\n", f)); next }
  bz <- best(z, idx); bg <- best(se, idx)
  cat(sprintf("%-8s %-10s %-16s %10.2f %+10.2f | %-16s %10.2f %+10.2f\n",
              f, bz$motif, bz$bin, bz$nlp, bz$enr, bg$bin, bg$nlp, bg$enr))
  out[[f]] <- data.frame(family = f, motif_zero = bz$motif, bin_zero = bz$bin,
                         nlp_zero = bz$nlp, enr_zero = bz$enr,
                         motif_gen = bg$motif, bin_gen = bg$bin,
                         nlp_gen = bg$nlp, enr_gen = bg$enr)
}
res_tab <- do.call(rbind, out)
write.table(res_tab, "/mnt/d/ngs_rebuild/results/background_comparison.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

