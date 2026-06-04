#!/bin/bash
set -euo pipefail

#file params
BAM_INDIR=../00_inputs/
BAM_OUTDIR=../00_inputs/merged_bams
PEAK_DIR=../01_outputs
EXCLUSION_LIST=../00_inputs/hg38_exclusion_ENCFF356LFX.bed.gz
SCRIPT_DIR=$(dirname "$(realpath "$0")")
LOG_DIR="${PEAK_DIR}/logs"

#env params
CONDA_INIT=/home/groups/engreitz/Software/anaconda3/etc/profile.d/conda.sh
CONDA_ENV=peak_callers

#specific run params
TIMEPOINTS=(0h 6h 24h)
GENES=$(ls "$BAM_INDIR") 
#SKIPGENES=(CHD4 YY1)

#GENES=(CHD4 YY1) 
SKIPGENES=()

mapfile -t FILTERED_GENES < <(
    printf '%s\n' "${GENES[@]}" |
    grep -vFf <(printf '%s\n' "${SKIPGENES[@]}")
)

#user checks
mkdir -p "$BAM_OUTDIR" "$LOG_DIR"

#main loop
for sample in "${FILTERED_GENES[@]}"; do

    for tp in "${TIMEPOINTS[@]}"; do
    
        bams=("${BAM_INDIR}/${sample}/${tp}"/*/ATAC-${sample}-${tp}-*.bam)

        if [[ ! -f "${bams[0]}" ]]; then
            echo "No BAMs found for ${sample} ${tp}, skipping"
            continue
        fi

        echo "Submitting ${sample} ${tp} BAM(s))"
        sbatch \
            --job-name="atac_${sample}_${tp}" \
            --output="${LOG_DIR}/${sample}_${tp}.%j.out" \
            --partition=engreitz \
            --cpus-per-task=6 \
            --mem=40GB \
            --time=12:00:00 \
            --export=ALL,CONDA_INIT="${CONDA_INIT}",CONDA_ENV="${CONDA_ENV}" \
            "${SCRIPT_DIR}/00_run_pipeline.sh" \
            "$sample" "$tp" "$BAM_OUTDIR" "$PEAK_DIR" "$EXCLUSION_LIST" "${bams[@]}"
    done
done
