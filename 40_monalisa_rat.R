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

pwms <- getMatrixSet(db(JASPAR2024()),
                     list(collection = "CORE", tax_group = "vertebrates",
                          matrixtype = "PWM"))
stopifnot(length(pwms) == 879)


se <- calcBinnedMotifEnrR(seqs = seqs, bins = bins, pwmL = pwms,
                          background = "zeroBin",
                          BPPARAM = SerialParam(), verbose = TRUE)
saveRDS(se, "/mnt/d/ngs_rebuild/results/monalisa_rat.rds")

sel <- apply(assay(se, "negLog10Padj"), 1,
             function(x) max(abs(x), 0, na.rm = TRUE)) > 4.0
cat(sprintf("motifs passing -log10(padj) > 4 in any bin: %d of %d\n", sum(sel), nrow(se)))

pdf("/mnt/d/ngs_rebuild/figures/monalisa_rat.pdf", width = 8, height = 14)
plotMotifHeatmaps(x = se[sel, ], which.plots = c("log2enr","negLog10Padj"),
                  width = 2.0, cluster = TRUE, maxEnr = 2, maxSig = 10,
                  show_motif_GC = TRUE)
dev.off()

FAM <- list(ONECUT=c("ONECUT1","ONECUT2","ONECUT3"), ISL=c("ISL1","ISL2","Isl1"),
            MEF2=c("MEF2A","MEF2B","MEF2C","MEF2D"),
            NFI=c("NFIA","NFIB","NFIC","NFIX","NFIC::TLX1"),
            DBOX=c("DBP","TEF","HLF","NFIL3"), CTCF_CONTROL=c("CTCF","CTCFL"))
nmv <- rowData(se)$motif.name
E <- assay(se,"log2enr"); P <- assay(se,"negLog10Padj")
bl <- colnames(E)
cat("\nbins (left = P11-high, right = E21-high):", paste(bl, collapse=" | "), "\n")
for (f in names(FAM)) {
  idx <- which(toupper(nmv) %in% toupper(FAM[[f]]))
  if (!length(idx)) next
  cat(sprintf("\n%s\n", f))
  for (i in idx) {
    j <- which.max(abs(E[i,]))
    cat(sprintf("   %-12s strongest bin %-14s log2enr %+6.2f  -log10padj %6.2f\n",
                nmv[i], bl[j], E[i,j], P[i,j]))
  }
}
