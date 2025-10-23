#!/usr/bin/env bash
#SBATCH --job-name=barbell_test
#SBATCH --partition=interactive
#SBATCH --cpus-per-task=80
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages/slurm_logs/%x/result_%x.%j.log

set -euo pipefail
IFS=$'\n\t'

echo "Started: $(date)"
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-NA}"

# -------------------------
# USER CONFIG (yours)
# -------------------------
WD="${WD:-/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages}"

# Can be relative to WD (as you provided) or absolute:
RAW_FASTQ="${RAW_FASTQ:-phages_analysis/ont_pipeline_20251009_134551/demux_trimmed/6b18c2cc-16b3-463f-9d76-d8db7ce6d962_SQK-NBD114-96_barcode01.fastq}"

OUTROOT="${OUTROOT:-$WD/phages_analysis/ont_pipeline_20251009_134551}"
OUT_BARBELL="${OUT_BARBELL:-$OUTROOT/barbell_test}"

KIT="${KIT:-SQK-NBD114-96}"

# Threads (use all cores unless overridden)
THREADS="${SLURM_CPUS_PER_TASK:-80}"

# -------------------------
# barbell binary on cluster
# -------------------------
export CARGO_HOME=/home/groups/VEO/tools/cargo
export RUSTUP_HOME=/home/groups/VEO/tools/rust
#BARBELL=/home/groups/VEO/tools/barbell/v0.1.9/barbell/target/release/barbell
BARBELL="/home/groups/VEO/tools/barbell/v0.1.8/barbell/target/release/barbell"
# -------------------------
# Resolve paths & prep I/O
# -------------------------
# Make RAW_FASTQ absolute if user passed it relative to WD
if [[ "$RAW_FASTQ" = /* ]]; then
  IN_FASTQ="$RAW_FASTQ"
else
  IN_FASTQ="$WD/$RAW_FASTQ"
fi

mkdir -p "$OUT_BARBELL"

echo "[INFO] WD=$WD"
echo "[INFO] IN_FASTQ=$IN_FASTQ"
echo "[INFO] OUT_BARBELL=$OUT_BARBELL"
echo "[INFO] KIT=$KIT"
echo "[INFO] THREADS=$THREADS"
echo "[INFO] BARBELL=$BARBELL"

# -------------------------
# Run barbell (kit preset)
# -------------------------
# --maximize: favor yield; add/remove flags below as you prefer
# Other handy flags (optional):
#   --no-orientation         # keep original read orientation
#   --min-length 1000        # drop very short fragments
#   --only-side left         # trim only left side
#   --verbose                # more logging
"$BARBELL" kit \
  -k "$KIT" \
  -i "$IN_FASTQ" \
  -o "$OUT_BARBELL" \
  --maximize \
  -t "$THREADS"

echo "Finished: $(date)"
