#!/usr/bin/env bash
#SBATCH --job-name=nanoplot_qc_phage_reads_no_trim
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --output=/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages/slurm_logs/%x/result_%x.%j.log

set -euo pipefail

# ===== USER VARS (override with sbatch --export) =====
WD="/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages"
READS_DIR="${WD}/phages_analysis/ont_basecalls_20251007_160801/demux_no_trim"
OUT_ROOT="${WD}/phages_analysis/ont_basecalls_20251007_160801/NanoPlot_out_no_trim"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ===== Env: conda NanoPlot =====
source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh
conda activate nanoplot_v1.41.3

mkdir -p "$OUT_ROOT"

shopt -s nullglob
mapfile -t FILES < <(find "$READS_DIR" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | sort)
if (( ${#FILES[@]} == 0 )); then
  echo "[ERR] No .fastq or .fastq.gz files in: $READS_DIR" >&2
  exit 2
fi

echo "[INFO] Found ${#FILES[@]} files. Output root: $OUT_ROOT"

for f in "${FILES[@]}"; do
  base="$(basename "$f")"
  name="${base%.gz}"
  name="${name%.fastq}"
  outdir="$OUT_ROOT/$name"

  if [[ -s "$outdir/NanoStats.txt" ]]; then
    echo "[SKIP] $base -> already has NanoPlot stats"
    continue
  fi

  echo "[RUN ] $base -> $outdir"
  mkdir -p "$outdir"
  # Common, informative flags: per-read TSV stats + loglength plots
  NanoPlot --fastq "$f" -o "$outdir" -t "$THREADS" --tsv_stats --loglength
  echo "[DONE] $base"
done

echo "[INFO] All reports written under: $OUT_ROOT"
