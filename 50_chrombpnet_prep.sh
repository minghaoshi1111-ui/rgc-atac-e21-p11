#!/usr/bin/env bash
#ran inside of runpod
set -eu

D=/workspace/data
P=/workspace/prep
mkdir -p "$P"

GENOME=$D/genome.fa
SIZES=$D/genome.canonical.sizes      # canonical chromosomes only
BL=$D/signal_outlier_bins.bed

echo "sanity checks"
echo "canonical chrom sizes"
cat "$SIZES"
echo
echo "chromosome names present in the peak files"
cut -f1 "$D"/E21.mRp.clN_peaks.narrowPeak | sort -u | tr '\n' ' '; echo
echo
echo "blacklist substitute: $(wc -l < "$BL") intervals, $(awk 'NR==1{print NF}' "$BL") columns "
head -3 "$BL"
echo
echo "peak file column count and it must be 10"
awk 'NR==1{print NF" columns; col10="$10}' "$D"/E21.mRp.clN_peaks.narrowPeak
echo

FIRST=$(cut -f1 "$SIZES" | head -1)
if [[ "$FIRST" == chr* ]]; then PFX="chr"; else PFX=""; fi
TEST_CHR="${PFX}1 ${PFX}3 ${PFX}6"
VAL_CHR="${PFX}8 ${PFX}20"
echo "chromosome prefix detected: '${PFX:-none}'"
echo "test chromosomes:       $TEST_CHR"
echo "validation chromosomes: $VAL_CHR"
echo

echo "=1. blacklist-filter peaks"
bedtools slop -i "$BL" -g "$SIZES" -b 1057 | sort -k1,1 -k2,2n > "$P/blacklist_slop.bed"
echo "slopped blacklist: $(wc -l < "$P/blacklist_slop.bed") intervals"

for S in E21 P11; do
  bedtools intersect -v \
    -a "$D/${S}.mRp.clN_peaks.narrowPeak" \
    -b "$P/blacklist_slop.bed" > "$P/${S}_peaks.narrowPeak"
  # keep only canonical chromosomes
  awk -v OFS='\t' 'NR==FNR{keep[$1];next} ($1 in keep)' \
      "$SIZES" "$P/${S}_peaks.narrowPeak" > "$P/${S}_peaks.filt.narrowPeak"
  mv "$P/${S}_peaks.filt.narrowPeak" "$P/${S}_peaks.narrowPeak"
  echo "$S: $(wc -l < "$D/${S}.mRp.clN_peaks.narrowPeak") -> $(wc -l < "$P/${S}_peaks.narrowPeak") peaks after blacklist + canonical filter"
done
echo

echo "2.shared chromosome fold"
chrombpnet prep splits \
  -c "$SIZES" \
  -tcr $TEST_CHR \
  -vcr $VAL_CHR \
  -op "$P/fold_0"
echo "--- fold_0.json ---"
cat "$P/fold_0.json"
echo

echo "3.GC-matched negatives"
for S in E21 P11; do
  echo "--- $S ---"
  chrombpnet prep nonpeaks \
    -g "$GENOME" \
    -p "$P/${S}_peaks.narrowPeak" \
    -c "$SIZES" \
    -fl "$P/fold_0.json" \
    -br "$BL" \
    -o "$P/${S}"
  echo "$S negatives: $(wc -l < "$P/${S}_negatives.bed")"
done

echo
echo "prep complete! "
ls -lh "$P"
