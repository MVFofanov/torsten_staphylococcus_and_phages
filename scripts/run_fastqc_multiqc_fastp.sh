#!/usr/bin/env bash
#SBATCH --job-name=fastqc_fastp_multiqc
#SBATCH --partition=short
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=80
#SBATCH --mem=80G
#SBATCH --output=/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages/slurm_logs/%x/result_%x.%j.log

set -euo pipefail

# -------------------- paths & settings --------------------
wd="/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages"
RAW_DIR="${wd}/data/250710_Staph_Laboklin/X208SC25056733-Z01-F001/01.RawData/StaphIII12"
SAMPLE="StaphIII12"
THREADS="${SLURM_CPUS_PER_TASK:-80}"

FASTQC_BIN="/home/groups/VEO/tools/fastqc/v0.12.1/fastqc"

# fastp env
source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh
conda activate fastp_v0.23.4

# output dirs
RAW_QC="${wd}/fastqc_raw"
FASTP_DIR="${wd}/fastp"
TRIM_DIR="${FASTP_DIR}/trimmed"
TRIM_QC="${wd}/fastqc_trimmed"
MERGED_DIR="${wd}/trimmed_merged"
MERGED_QC="${wd}/fastqc_merged"
MQC_DIR="${wd}/multiqc"

mkdir -p "$RAW_QC" "$FASTP_DIR" "$TRIM_DIR" "$TRIM_QC" "$MERGED_DIR" "$MERGED_QC" "$MQC_DIR"

# -------------------- (optional) MD5 check --------------------
if [[ -f "${RAW_DIR}/MD5.txt" ]]; then
  ( cd "${RAW_DIR}" && md5sum -c MD5.txt )
fi

# -------------------- FastQC on raw --------------------
"$FASTQC_BIN" -t "$THREADS" -o "$RAW_QC" "${RAW_DIR}"/*_1.fq.gz "${RAW_DIR}"/*_2.fq.gz

# -------------------- fastp per pair --------------------
for R1 in "${RAW_DIR}"/*_1.fq.gz; do
  [[ -e "$R1" ]] || { echo "No R1 files found in ${RAW_DIR}"; exit 1; }
  R2="${R1/_1.fq.gz/_2.fq.gz}"
  [[ -f "$R2" ]] || { echo "Missing R2 for $(basename "$R1")"; exit 1; }

  # lane tag like L3 / L6 if present; fallback to basename stem
  base="$(basename "$R1")"
  lane="$(sed -n 's/.*_\(L[0-9]\+\)_.*/\1/p' <<<"$base")"
  [[ -n "$lane" ]] || lane="${base%%_1.fq.gz}"

  fastp \
    -i "$R1" -I "$R2" \
    -o "${TRIM_DIR}/${lane}_R1.trim.fq.gz" \
    -O "${TRIM_DIR}/${lane}_R2.trim.fq.gz" \
    --detect_adapter_for_pe \
    --trim_poly_g --trim_poly_x \
    -q 20 -u 30 -n 3 -l 50 \
    --thread 16 \
    --report_title "${SAMPLE} ${lane}" \
    --html "${FASTP_DIR}/${lane}_fastp.html" \
    --json "${FASTP_DIR}/${lane}_fastp.json"
done

# -------------------- FastQC on trimmed --------------------
"$FASTQC_BIN" -t "$THREADS" -o "$TRIM_QC" "${TRIM_DIR}"/*.trim.fq.gz

# -------------------- Merge lanes --------------------
cat $(ls -1 "${TRIM_DIR}"/*_R1.trim.fq.gz | sort) > "${MERGED_DIR}/${SAMPLE}_R1.trim.fq.gz"
cat $(ls -1 "${TRIM_DIR}"/*_R2.trim.fq.gz | sort) > "${MERGED_DIR}/${SAMPLE}_R2.trim.fq.gz"

# -------------------- FastQC on merged --------------------
"$FASTQC_BIN" -t "$THREADS" -o "$MERGED_QC" \
  "${MERGED_DIR}/${SAMPLE}_R1.trim.fq.gz" \
  "${MERGED_DIR}/${SAMPLE}_R2.trim.fq.gz"

# -------------------- MultiQC (your env) --------------------
# Activates MultiQC v1.15 env and runs it
source /home/groups/VEO/tools/multiQC/v1.15/bin/activate
multiqc -o "$MQC_DIR" "$RAW_QC" "$TRIM_QC" "$MERGED_QC" "$FASTP_DIR"

echo "Done."
