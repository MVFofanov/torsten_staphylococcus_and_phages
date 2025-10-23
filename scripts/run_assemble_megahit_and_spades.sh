#!/usr/bin/env bash
#SBATCH --job-name=assembly_megahit_and_spades
#SBATCH --partition=short
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=80
#SBATCH --mem=80G
#SBATCH --output=/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages/slurm_logs/%x/%x_%j.log
#SBATCH --open-mode=append

set -euo pipefail

# ---------------- paths & inputs ----------------
WD="/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages"
R1="${WD}/trimmed_merged/StaphIII12_R1.trim.fq.gz"
R2="${WD}/trimmed_merged/StaphIII12_R2.trim.fq.gz"

MEGAHIT_ENV="megahit_v1.2.9"
CONDA_SETUP="/vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh"
SPADES_BIN="/home/groups/VEO/tools/SPAdes/v3.15.5/spades.py"

THREADS="${SLURM_CPUS_PER_TASK:-80}"
MEM_GB="80"   # for SPAdes -m

# ---------------- preflight checks ----------------
echo "[INFO] Job: ${SLURM_JOB_NAME:-local}  ID: ${SLURM_JOB_ID:-NA}"
echo "[INFO] Using ${THREADS} threads; ${MEM_GB} GB mem"
echo "[INFO] WD: $WD"
[[ -r "$R1" && -r "$R2" ]] || { echo "[ERROR] R1/R2 not found"; exit 1; }

mkdir -p "${WD}/spades"

# ---------------- MEGAHIT ----------------
echo "[STEP] MEGAHIT"
# activate conda env with megahit
# shellcheck disable=SC1090
source "$CONDA_SETUP"
conda activate "$MEGAHIT_ENV"

# echo "[INFO] MEGAHIT version:"
# megahit -v || true
# echo "[INFO] MEGAHIT help (truncated):"
# megahit -h | head -n 20 || true

# time megahit \
#   -1 "$R1" -2 "$R2" \
#   -o "${WD}/megahit" \
#   --min-contig-len 5000 \
#   -t "$THREADS"

# # symlink a friendly name
# ln -sf "${WD}/megahit/final.contigs.fa" "${WD}/megahit/contigs.fa" || true

# conda deactivate

# ---------------- SPAdes ----------------
echo "[STEP] SPAdes"
echo "[INFO] SPAdes path: $SPADES_BIN"
"$SPADES_BIN" --version || true
echo "[INFO] SPAdes help (truncated):"
"$SPADES_BIN" -h | head -n 30 || true

# Use --isolate + --careful for bacterial isolates; cov-cutoff auto is typical
time "$SPADES_BIN" \
  -1 "$R1" -2 "$R2" \
  -o "${WD}/spades" \
  -t "$THREADS" -m "$MEM_GB" \
  --isolate --cov-cutoff auto

# convenience links
ln -sf "${WD}/spades/contigs.fasta"   "${WD}/spades/contigs.fa"   || true
ln -sf "${WD}/spades/scaffolds.fasta" "${WD}/spades/scaffolds.fa" || true

# ---------------- post-run summary ----------------
echo "[DONE] Assemblies finished."
echo "[OUT] MEGAHIT contigs: ${WD}/megahit/contigs.fa"
echo "[OUT] SPAdes contigs : ${WD}/spades/contigs.fa"
echo "[OUT] SPAdes scaffolds: ${WD}/spades/scaffolds.fa"

# quick size stats
for f in "${WD}/megahit/contigs.fa" "${WD}/spades/contigs.fa" "${WD}/spades/scaffolds.fa"; do
  if [[ -f "$f" ]]; then
    echo "==> $(basename "$f")"
    awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {seqlen+=length($0)} END {print seqlen}' "$f" \
      | awk '{sum+=$1; if($1>max)max=$1; n++;} END {if(n){printf "  contigs=%d  total_bp=%d  max_contig=%d\n", n, sum, max} else {print "  empty"}}'
  fi
done
