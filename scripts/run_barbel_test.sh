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
# USER CONFIG
# -------------------------
WD="${WD:-/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages}"

RAW_FASTQ="${RAW_FASTQ:-phages_analysis/ont_pipeline_20251009_134551/demux_trimmed/6b18c2cc-16b3-463f-9d76-d8db7ce6d962_SQK-NBD114-96_barcode01.fastq}"

OUTROOT="${OUTROOT:-$WD/phages_analysis/ont_pipeline_20251009_134551}"
OUTDIR="${OUTDIR:-$OUTROOT/barbell_v0.2.0_test}"

THREADS="${SLURM_CPUS_PER_TASK:-80}"

# v0.1.8 binary (no 'kit' subcommand)
# export CARGO_HOME=/home/groups/VEO/tools/cargo
# export RUSTUP_HOME=/home/groups/VEO/tools/rust
# BARBELL="/home/groups/VEO/tools/barbell/v0.1.8/barbell/target/release/barbell"

export CARGO_HOME=/home/zo49sog/miniconda3/envs/barbell_rust_env/.cargo
export RUSTUP_HOME=/home/zo49sog/miniconda3/envs/barbell_rust_env/.rustup
BARBELL="/work/zo49sog/barbell/target/release/barbell"

# -------------------------
# Resolve paths & prep I/O
# -------------------------
if [[ "$RAW_FASTQ" = /* ]]; then
  IN_FASTQ="$RAW_FASTQ"
else
  IN_FASTQ="$WD/$RAW_FASTQ"
fi

mkdir -p "$OUTDIR"

echo "[INFO] IN_FASTQ=$IN_FASTQ"
echo "[INFO] OUTDIR=$OUTDIR"
echo "[INFO] THREADS=$THREADS"
echo "[INFO] BARBELL=$BARBELL"

# -------------------------
# 1) ANNOTATE (with kit)
# -------------------------
# v0.1.8 supports --kit in 'annotate'.
# Tip: if too few reads match, you can relax flanks a bit with: --flank-max-errors 5
# (but always inspect afterwards). :contentReference[oaicite:1]{index=1}
#  --kit SQK-NBD114-96 \
"$BARBELL" annotate \
  -i "$IN_FASTQ" \
  -o "$OUTDIR/anno.tsv" \
  -t "$THREADS"

# -------------------------
# 2) INSPECT (summaries + per-read pattern)
# -------------------------
"$BARBELL" inspect -i "$OUTDIR/anno.tsv" > "$OUTDIR/inspect_summary.txt"
"$BARBELL" inspect -i "$OUTDIR/anno.tsv" -o "$OUTDIR/pattern_per_read.tsv"

# -------------------------
# 3) FILTER (choose patterns and define cut positions)
# -------------------------
# Keep classic native-kit cases:
#  • Left barcode only (cut after left tag) 
#  • Left + right barcode (keep the insert: cut after left, before right)
# Pattern & cut syntax documented in PDF. 
cat > "$OUTDIR/filters.txt" <<'EOF'
Ftag[fw, *, @left(0..250), >>]
Ftag[fw, *, @left(0..250), >>]__Ftag[<<, rc, *, @right(0..250)]
EOF

# Produce filtered.tsv with cuts column (required for trim). :contentReference[oaicite:3]{index=3}
"$BARBELL" filter -i "$OUTDIR/anno.tsv" -f "$OUTDIR/filters.txt" -o "$OUTDIR/filtered.tsv"

# -------------------------
# 4) TRIM (write FASTQs)
# -------------------------
# Simplify filenames: keep only left label, ignore orientation → NB01.trimmed.fastq, etc. :contentReference[oaicite:4]{index=4}
mkdir -p "$OUTDIR/trimmed"
"$BARBELL" trim \
  -i "$OUTDIR/filtered.tsv" \
  -r "$IN_FASTQ" \
  -o "$OUTDIR/trimmed" \
  --no-orientation \
  --only-side left

echo "Finished: $(date)"
