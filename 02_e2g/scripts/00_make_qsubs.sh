
QSUB_DIR=../00_inputs/qsubs
CONFIG_DIR=../00_inputs/configs
RESULTS_DIR=../01_outputs
BAM_INDIR=../../00_Peaks/00_inputs/merged_bams
PEAK_INDIR=../../00_Peaks/01_outputs

mkdir -p ${QSUB_DIR} ${CONFIG_DIR}

# Iterate over sample names (subdirectories of PEAK_INDIR)
for sample in $(ls ${PEAK_INDIR})
do

config_file=${CONFIG_DIR}/config_${sample}.tsv

# Validate all required files before writing anything
all_missing=()
valid_timepoints=()

for tp in 0h 6h 24h; do
    # Only check timepoints that have a peak directory
    [ ! -d "${PEAK_INDIR}/${sample}/${tp}" ] && continue

    bam="${BAM_INDIR}/ATAC-${sample}-${tp}-merged.possort.bam"
    bai="${bam}.bai"
    peak="${PEAK_INDIR}/${sample}/${tp}/ATAC-${sample}-${tp}_peaks.narrowPeak"
    summit="${PEAK_INDIR}/${sample}/${tp}/ATAC-${sample}-${tp}_summits.bed"

    [ ! -f "$bam" ]    && all_missing+=("$bam")
    [ ! -f "$bai" ]    && all_missing+=("$bai")
    [ ! -f "$peak" ]   && all_missing+=("$peak")
    [ ! -f "$summit" ] && all_missing+=("$summit")

    valid_timepoints+=("$tp")
done

if [ ${#all_missing[@]} -gt 0 ]; then
    echo "ERROR: Skipping ${sample} — missing files:"
    printf "  %s\n" "${all_missing[@]}"
    continue
fi

if [ ${#valid_timepoints[@]} -eq 0 ]; then
    echo "ERROR: No timepoint directories found for ${sample}, skipping"
    continue
fi

# All files present — write config and create peak symlinks
echo -e "biosample\tDHS\tATAC\tH3K27ac\tdefault_accessibility_feature\tHiC_file\tHiC_type\tHiC_resolution\talt_TSS\talt_genes" > ${config_file}

for tp in "${valid_timepoints[@]}"; do
    bam="${BAM_INDIR}/ATAC-${sample}-${tp}-merged.possort.bam"
    peak="${PEAK_INDIR}/${sample}/${tp}/ATAC-${sample}-${tp}_peaks.narrowPeak"

    echo -e "${sample}_${tp}\t\t${bam}\t\tATAC\thttps://s3.us-central-1.wasabisys.com/aiden-encode-hic-mirror/bifocals_iter2/tissues.hic\thic\t5000" >> ${config_file}

    # Pre-place peak symlinks so Snakemake skips MACS2 peak calling
    e2g_peak_dir="${RESULTS_DIR}${sample}/${sample}_${tp}/Peaks"
    mkdir -p "${e2g_peak_dir}"
    ln -sf "${peak}" "${e2g_peak_dir}/macs2_peaks.narrowPeak"
    touch "${e2g_peak_dir}/macs2_peaks.narrowPeak"
done

#note: there are explicit min_gb=20 for some steps, so mem needs to be at least more than that.
cat <<EOF > ${QSUB_DIR}/${sample}_e2g_qsub.sh
#!/bin/bash
#SBATCH --job-name=${sample}
#SBATCH --output=${QSUB_DIR}/${sample}.%j.out
#SBATCH --partition=engreitz
#SBATCH --cpus-per-task=6
#SBATCH --mem=80GB
#SBATCH --time=14:00:00
set -e

source /home/groups/engreitz/Software/anaconda3/etc/profile.d/conda.sh
conda activate encode_re2g

#####------#####
echo "Running E2G for ${sample}"

cd ../ENCODE_rE2G

snakemake \
    --snakefile workflow/Snakefile \
    --configfile config/config.yaml \
    --config ABC_BIOSAMPLES=${config_file} \
             results_dir=${RESULTS_DIR}${sample} \
    --use-conda \
    --conda-prefix /oak/stanford/groups/engreitz/Users/aeveritt/.conda \
    --cores \${SLURM_CPUS_PER_TASK:-1} \
    --rerun-incomplete \
    --keep-going \
    --rerun-triggers mtime

echo "Complete"
#####------#####

[[ -n "\$SLURM_JOB_ID" ]] && scontrol show job "\$SLURM_JOB_ID"

EOF

done
