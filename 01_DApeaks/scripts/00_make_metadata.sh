#!/bin/bash
set -euo pipefail

BAM_DIR1=../../00_Peaks/00_inputs
BAM_DIR2=../../00_Peaks/00_inputs/merged_bams
PEAK_DIR=../../00_Peaks/01_outputs
OUT_DIR=../00_inputs
OUT_TSV="${OUT_DIR}/metadata.tsv"

TIMEPOINTS=(0h 6h 24h)
EXPECTED_REPS=2

HEADER="gene\tsampleID\ttimepoint\tbam_path\tmerged_bam_path\tpeak_path"
echo -e "$HEADER" > "$OUT_TSV"

for sample_dir in "$BAM_DIR1"/*/; do
    gene=$(basename "$sample_dir")
    gene_tsv="${OUT_DIR}/${gene}.metadata.tsv"
    echo -e "$HEADER" > "$gene_tsv"

    for tp in "${TIMEPOINTS[@]}"; do

        # Merged BAM — one per gene/timepoint
        merged_bam="${BAM_DIR2}/ATAC-${gene}-${tp}-merged.possort.bam"
        [[ -f "$merged_bam" ]] && merged_bam_val="$merged_bam" || merged_bam_val="NULL"

        # Peak file — prefer filtered over unfiltered, .gz over plain
        peak_val="NULL"
        for peak_candidate in \
            "${PEAK_DIR}/${gene}/${tp}/ATAC-${gene}-${tp}_peaks.filtered.narrowPeak.gz" \
            "${PEAK_DIR}/${gene}/${tp}/ATAC-${gene}-${tp}_peaks.filtered.narrowPeak" \
            "${PEAK_DIR}/${gene}/${tp}/ATAC-${gene}-${tp}_peaks.narrowPeak.gz" \
            "${PEAK_DIR}/${gene}/${tp}/ATAC-${gene}-${tp}_peaks.narrowPeak"; do
            if [[ -f "$peak_candidate" ]]; then
                peak_val="$peak_candidate"
                break
            fi
        done

        # One row per replicate — sampleID is the replicate directory name
        rep_count=0
        for rep_dir in "${BAM_DIR1}/${gene}/${tp}"/*/; do
            [[ -d "$rep_dir" ]] || continue
            sampleID=$(basename "$rep_dir")

            bams=("${rep_dir}"*.bam)
            [[ -f "${bams[0]}" ]] && bam_val="${bams[0]}" || bam_val="NULL"

            echo -e "${gene}\t${sampleID}\t${tp}\t${bam_val}\t${merged_bam_val}\t${peak_val}" | tee -a "$OUT_TSV" >> "$gene_tsv"
            rep_count=$((rep_count + 1))
        done

        # Pad with NULL rows for any missing replicates
        while [[ $rep_count -lt $EXPECTED_REPS ]]; do
            echo -e "${gene}\tNULL\t${tp}\tNULL\t${merged_bam_val}\t${peak_val}" | tee -a "$OUT_TSV" >> "$gene_tsv"
            rep_count=$((rep_count + 1))
        done
    done
done

echo "Written: $OUT_TSV and per-gene TSVs in ${OUT_DIR}/"
