#!/usr/bin/env Rscript
# 58_footprint_motifs.R

suppressMessages({library(JASPAR2024); library(TFBSTools)})
setwd("/mnt/d/ngs_rebuild")
dir.create("results", showWarnings = FALSE)

pfms <- getMatrixSet(db(JASPAR2024()),
                     list(collection = "CORE", tax_group = "vertebrates",
                          matrixtype = "PFM"))
nms <- vapply(pfms, TFBSTools::name, character(1))

consensus <- function(p) {
  m <- as.matrix(Matrix(p))
  paste(rownames(m)[apply(m, 2, which.max)], collapse = "")
}

want <- c(MEF2 = "MEF2C", NFI = "NFIC", Dbox = "DBP",
          ONECUT = "ONECUT2", ISL = "ISL1")

rows <- list(c("CTCF", "TTGGCCACTAGGGGGCGCTAT"))   # published, Zenodo 7443683
for (fam in names(want)) {
  hit <- which(nms == want[[fam]])
  if (!length(hit)) { cat("NOT FOUND in JASPAR:", want[[fam]], "\n"); next }
  s <- consensus(pfms[[hit[1]]])
  rows[[length(rows) + 1]] <- c(fam, s)
  cat(sprintf("  %-8s %-10s %s\n", fam, want[[fam]], s))
}

out <- do.call(rbind, rows)
write.table(out, "results/motif_to_pwm.rat.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("\nwrote results/motif_to_pwm.rat.tsv\n")
writeLines(readLines("results/motif_to_pwm.rat.tsv"))
