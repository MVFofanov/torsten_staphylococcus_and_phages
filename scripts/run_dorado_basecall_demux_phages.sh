#!/usr/bin/env bash
#SBATCH --job-name=dorado_basecall_demux_phages
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --output=/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages/slurm_logs/%x/result_%x.%j.log

set -euo pipefail

# ------------------- USER PATHS -------------------
WD="${WD:-/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages}"
POD5_DIR="${POD5_DIR:-/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages/data/250715_StaphPhages_TS/250715_StaphPhages_TS/20250715_0915_MD-102423_FBA89536_6b18c2cc/pod5_skip}"
OUT="${OUT:-${WD}/phages_analysis/ont_basecalls_$(date +%Y%m%d_%H%M%S)}"

# Chemistry/model & kit
MODEL="${MODEL:-dna_r10.4.1_e8.2_400bps_sup@v4.3.0}"     # or MODEL=sup to auto-pick latest compatible
KIT="${KIT:-SQK-NBD114-96}"                               # exact barcode kit (Kit14)

# Dorado binary
DORADO="${DORADO:-/home/groups/VEO/tools/dorado/v1.1.0/bin/dorado}"

echo "[INFO] Job: $SLURM_JOB_NAME id=$SLURM_JOB_ID"
echo "[INFO] POD5_DIR: $POD5_DIR"
echo "[INFO] OUT:      $OUT"
echo "[INFO] MODEL:    $MODEL"
echo "[INFO] KIT:      $KIT"
echo "[INFO] DORADO:   $DORADO"

mkdir -p "$OUT"
ALL_FASTQ="$OUT/all.fastq.gz"
DEMUX_DIR="$OUT/demux_no_trim"            # keep barcode sequence
ASM_DIR="$OUT/demux_trimmed_for_assembly" # barcode/adapters/primers trimmed

# compressor
if command -v pigz >/dev/null 2>&1; then
  COMPRESS="pigz -c -p ${SLURM_CPUS_PER_TASK:-8}"
else
  COMPRESS="gzip -c"
fi

# Ensure model is present (safe to re-run)
"$DORADO" download --model "$MODEL" || true

# Pick device
if command -v nvidia-smi >/dev/null 2>&1; then DEV="cuda:all"; else DEV="cpu"; fi
echo "[INFO] Device: $DEV"

# ------------------- Step 1: Basecall all POD5 recursively (no trimming) -------------------
if [[ -s "$ALL_FASTQ" ]]; then
  echo "[SKIP] Basecalling exists: $ALL_FASTQ"
else
  echo "[RUN ] Basecalling -> $ALL_FASTQ"
  tmp="$ALL_FASTQ.partial"
  "$DORADO" basecaller "$MODEL" "$POD5_DIR" --recursive --device "$DEV" \
    --emit-fastq --no-trim \
    | ${COMPRESS} > "$tmp"
  mv "$tmp" "$ALL_FASTQ"
  echo "[DONE] Basecalling."
fi

# ------------------- Step 2: Demultiplex WITHOUT trimming (preserve barcodes) -------------
if [[ -d "$DEMUX_DIR" && ( -s "$DEMUX_DIR/barcode01.fastq" || -s "$DEMUX_DIR/barcode01.fastq.gz" ) ]]; then
  echo "[SKIP] Demux (no-trim) appears done: $DEMUX_DIR"
else
  echo "[RUN ] Demultiplex (no-trim) -> $DEMUX_DIR"
  mkdir -p "$DEMUX_DIR"
  "$DORADO" demux --kit-name "$KIT" --no-trim --emit-fastq "$ALL_FASTQ" --output-dir "$DEMUX_DIR" --emit-summary
  echo "[DONE] Demultiplex (no-trim)."
fi

# ------------------- Step 3 (new): Trim a COPY for assembly -------------------------------
# We keep original demuxed FASTQs untouched in $DEMUX_DIR, and produce assembly-ready FASTQs in $ASM_DIR.
if [[ -d "$ASM_DIR" && -s "$ASM_DIR/trim_done.flag" ]]; then
  echo "[SKIP] Trimming for assembly already done: $ASM_DIR"
else
  echo "[RUN ] Trimming demuxed reads for assembly -> $ASM_DIR"
  mkdir -p "$ASM_DIR"
  shopt -s nullglob
  for fq in "$DEMUX_DIR"/barcode*.fastq; do
    b="$(basename "$fq" .fastq)"
    # dorado trim expects kit name for Kit14 sets
    "$DORADO" trim "$fq" --sequencing-kit "$KIT" > "$ASM_DIR/${b}.trimmed.fastq"
  done
  # (optional) gzip the trimmed files
  if compgen -G "$ASM_DIR/*.trimmed.fastq" > /dev/null; then
    find "$ASM_DIR" -type f -name "*.trimmed.fastq" -print0 | xargs -0 -n1 -P "${SLURM_CPUS_PER_TASK:-8}" gzip -f
  fi
  touch "$ASM_DIR/trim_done.flag"
  echo "[DONE] Trimming for assembly."
fi

echo "[INFO] Results:"
echo "  Basecalled FASTQ           : $ALL_FASTQ"
echo "  Demux (barcodes preserved) : $DEMUX_DIR/*.fastq"
echo "  Assembly-ready FASTQs      : $ASM_DIR/*.trimmed.fastq.gz"
