**Analysis scripts for "A Developmental Shift from ONECUT to MEF2 Motif Enrichment in Rat Retinal Ganglion Cells"**
Archived at Zenodo: https://doi.org/10.5281/zenodo.22761406

These are the nine scripts that produced the numbers in the paper, plus the seven files they read or write. Nothing has been tidied, re-pathed or re-run for the deposit: each script is the script that produced the result.

**Data**
No new data were generated. Everything starts from public records: GEO GSE163564 (with GSE163562 and GSE163563), BioProject PRJNA686776, SRA study SRP298616, genome mRatBN7.2 (GCA_015227675.2, RefSeq GCF_015227675.2), annotation Ensembl release 112. Alignment and peak calling were done with nf-core/atacseq as described in Methods 4.3; software versions are in Methods 4.11.

**Run order**

1. 51_diff_beds_1.R, then 54_gene_lists.R, then 55_gprofiler_api.py
2. 40_monalisa_rat.R, then 76_monalisa_genomebg.R, then 77_background_readout.R
3. 25_signal_outliers.R, then 50_chrombpnet_prep.sh, then chromBPNet training and footprinting (commands in Methods 4.6 and 4.8), then 58_footprint_motifs_1.R

**What each script does**
51_diff_beds_1.R reads the consensus count matrix and writes E21_up.bed, P11_up.bed, all_tested.bed and static.bed into results/diff_beds/.
54_gene_lists.R reads the peak annotation and those region files, and writes E21_up_genes.txt, P11_up_genes.txt and background_genes.txt into results/gene_lists/.
55_gprofiler_api.py sends those gene lists to g:Profiler and writes GO_E21_up.tsv, GO_P11_up.tsv and GO_query_provenance.json into results/.
40_monalisa_rat.R reads the count matrix and the genome FASTA and writes monalisa_rat.rds.
76_monalisa_genomebg.R repeats that run against a GC-matched genome background and writes monalisa_rat_genomebg.rds and background_comparison.tsv.
77_background_readout.R reads those three files and prints the values reported in Table 4 and Figure 2E. Writes no file.
25_signal_outliers.R reads a binned count table and writes qc/signal_outlier_bins.bed, the empirical exclusion set.
50_chrombpnet_prep.sh runs inside the chromBPNet container. It filters the per-sample peaks against that exclusion set and writes the splits and the GC-matched negatives used for training.
58_footprint_motifs_1.R builds motif_to_pwm.rat.tsv, the -pwm_f input for chrombpnet footprints.

**Where the data files came from**

consensus_peaks.mLb.clN.featureCounts.txt.gz and consensus_peaks.mLb.clN.annotatePeaks.txt.gz came from nfcore/results/bwa/merged_library/macs2/narrow_peak/consensus/. Both are gzipped to stay under the upload limit.
E21_up.bed, P11_up.bed, all_tested.bed and static.bed came from results/diff_beds/. 
signal_outlier_bins.bed came from qc/.

**Four things to know before running anything**
Paths are absolute and unchanged. Two scripts read the count matrix from /mnt/d/ngs_rebuild/nfcore/..., and the other R scripts set their working directory to /mnt/d/ngs_rebuild. Edit those lines for your own layout.
50_chrombpnet_prep.sh expects /workspace/data and /workspace/prep inside the container, so signal_outlier_bins.bed has to be copied to /workspace/data.
25_signal_outliers.R needs five files that are not here: counts/bin_counts.tsv for its main step, and peaks/consensus.named.bed, results/robust_FDR05.bed, results/composite_up.bed and results/composite_dn.bed for its closing overlap check. All five come back from the count matrix in this repository by the calls in Methods 4.3 and 4.4, or on request.
Other intermediate and derived files are not deposited. They are reproduced by running the scripts in the order above, and are available from the author on request.

**Citing this repository**
Cite the paper once it is published, and this repository by its archived release (see CITATION.cff). Until the paper is out, the repository release is the only citable record.

Licence
Code is MIT, see LICENSE. The derived tables and region files are CC BY 4.0.
