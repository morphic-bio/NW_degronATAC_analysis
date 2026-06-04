#!/bin/bash
set -euo pipefail

#file params
OUTDIR=../01_outputs
META_DIR=../00_inputs
SCRIPT_DIR=$(dirname "$(realpath "$0")")
LOG_DIR="${OUTDIR}/logs"

mapfile -t GENES < <(ls "${META_DIR}"/*.metadata.tsv | xargs -n1 basename)
GENES=("${GENES[@]%.metadata.tsv}")
SKIPGENES=()

mapfile -t FILTERED_GENES < <(
    if [[ ${#SKIPGENES[@]} -gt 0 ]]; then
        printf '%s\n' "${GENES[@]}" |
        grep -vFf <(printf '%s\n' "${SKIPGENES[@]}")
    else
        printf '%s\n' "${GENES[@]}"
    fi
)

mkdir -p "$OUTDIR" "$LOG_DIR"

for sample in "${FILTERED_GENES[@]}"; do
  meta="${META_DIR}/${sample}.metadata.tsv"
  
  cat <<EOF > ${LOG_DIR}/${sample}_qsub.sh
#!/bin/bash
#SBATCH --job-name=${sample}
#SBATCH --output=${LOG_DIR}/${sample}.%j.out
#SBATCH --partition=engreitz
#SBATCH --cpus-per-task=6
#SBATCH --mem=80GB
#SBATCH --time=10:00:00
set -e

source /home/groups/engreitz/Software/anaconda3/etc/profile.d/conda.sh 
conda activate degron_env2
    
Rscript ${SCRIPT_DIR}/01_atac_pipeline.R \
  --meta ${meta} \
  -o ${OUTDIR} \
  --cores \${SLURM_CPUS_PER_TASK:-1} \
  --saveOnFail
    
[[ -n "\$SLURM_JOB_ID" ]] && scontrol show job "\$SLURM_JOB_ID"

EOF

echo sbatch ${LOG_DIR}/${sample}_qsub.sh

done

#   --skipMotifEnrich \
#   --skipGSEA \
#   --skipFixedPeaks 
