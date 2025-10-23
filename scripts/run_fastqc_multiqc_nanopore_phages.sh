#!/usr/bin/env bash
#SBATCH --job-name=fastqc_multiqc_nanopore_prophages_no_trim
#SBATCH --partition=short
#SBATCH --cpus-per-task=80
#SBATCH --mem=80G
#SBATCH --time=3:00:00
#SBATCH --output=/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages/slurm_logs/%x/result_%x.%j.log


set -euo pipefail

# ---- USER VARS (override with sbatch --export) ----
WD="/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages"
READS_DIR="${WD}/phages_analysis/ont_basecalls_20251007_160801/demux_no_trim"
OUT_DIR="${WD}/phages_analysis/ont_basecalls_20251007_160801/QC_fastqc_no_trim"

THREADS="${SLURM_CPUS_PER_TASK:-80}"

FastQC="/home/groups/VEO/tools/fastqc/v0.12.1/fastqc"

mkdir -p "$OUT_DIR"

# Find inputs
mapfile -t FILES < <(find "$READS_DIR" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | sort)
(( ${#FILES[@]} > 0 )) || { echo "[ERR] No FASTQ files in $READS_DIR"; exit 2; }

echo "[INFO] Running FastQC on ${#FILES[@]} files -> $OUT_DIR"
# FastQC can take multiple files at once; -t sets threads
# It will skip already-finished files (html/zip present).
$FastQC -t "$THREADS" -o "$OUT_DIR" "${FILES[@]}"

echo "[INFO] Running MultiQC"
source /home/groups/VEO/tools/multiQC/v1.15/bin/activate

multiqc -o "$OUT_DIR" "$OUT_DIR"

echo "[INFO] Done. Open: $OUT_DIR/multiqc_report_no_trim.html"
