#!/bin/bash
set -euo pipefail

SAMPLE=$1
TP=$2
BAM_OUTDIR=$3
PEAK_DIR=$4
EXCL_LIST=${5:-}
shift 5
BAMS=("$@")

NTH="${SLURM_CPUS_PER_TASK:-1}"
PREFIX="ATAC-${SAMPLE}-${TP}"
MERGED_BAM_possrt="${BAM_OUTDIR}/${PREFIX}-merged.possort.bam"
OUT_DIR="${PEAK_DIR}/${SAMPLE}/${TP}"
TAG="${OUT_DIR}/${PREFIX}.tagAlign.gz"

mkdir -p "$BAM_OUTDIR" "$OUT_DIR"

if [[ -n "${CONDA_INIT:-}" && -n "${CONDA_ENV:-}" ]]; then
    source "$CONDA_INIT"
    conda activate "$CONDA_ENV"
fi

#####------#####
echo; echo "##### Step 1: Merging and position-sorting BAMs for ${PREFIX}"

if [ ! -f "$MERGED_BAM_possrt" ]; then
    samtools merge -@ "$NTH" -u - "${BAMS[@]}" | samtools sort -@ "$NTH" -o "$MERGED_BAM_possrt"
    samtools index "$MERGED_BAM_possrt"
    echo "Written: $MERGED_BAM_possrt"
else
    echo "Skipping (already exists): $MERGED_BAM_possrt"
fi

#####------#####
echo; echo "##### Step 2: Converting to TagAlign for ${PREFIX}"
# bedpe requires name-sorted input; pipe through samtools sort -n on the fly to avoid writing nmsort bam

if [ ! -f "$TAG" ]; then
    BEDPE="${OUT_DIR}/${PREFIX}.bedpe.gz"

    LC_COLLATE=C samtools sort -@ "$NTH" -n "$MERGED_BAM_possrt" \
        | bedtools bamtobed -bedpe -mate1 -i /dev/stdin | gzip -nc > "$BEDPE"

    zcat -f "$BEDPE" | awk 'BEGIN{OFS="\t"}{
        printf "%s\t%s\t%s\tN\t1000\t%s\n%s\t%s\t%s\tN\t1000\t%s\n",
        $1,$2,$3,$9,$4,$5,$6,$10
    }' | gzip -nc > "$TAG"

    rm -f "$BEDPE"
    echo "Written: $TAG"
else
    echo "Skipping (already exists): $TAG"
fi

#####------#####
echo; echo "##### Step 3: Calling peaks for ${PREFIX}"

macs3 callpeak \
    -t "$TAG" \
    -f BED \
    -g hs \
    -n "$PREFIX" \
    --outdir "$OUT_DIR" \
    -p 0.01 \
    --keep-dup all \
    --call-summits \
    --nomodel \
    --shift -75 \
    --extsize 150

#####------#####
echo; echo "##### Step 4: Removing Exclusion List from Peaks (NOT summits) ${PREFIX}"

if [[ -n "$EXCL_LIST" ]]; then
    NARROWPEAK="${OUT_DIR}/${PREFIX}_peaks.narrowPeak"
    FILTERED="${OUT_DIR}/${PREFIX}_peaks.filtered.narrowPeak"
    bedtools intersect \
        -v \
        -a "$NARROWPEAK" \
        -b "$EXCL_LIST" > "$FILTERED"
    echo "Written: $FILTERED"
else
    echo "Skipping exclusion list filtering (no EXCL_LIST provided)"
fi

#####------#####
echo; echo "##### Step 5: Cleaning up ${PREFIX}"
rm "$TAG"
rm "${OUT_DIR}/${PREFIX}_peaks.xls"
gzip "${OUT_DIR}/*.narrowPeak"
gzip "${OUT_DIR}/*.bed"

echo "Done: ${SAMPLE} ${TP}"

