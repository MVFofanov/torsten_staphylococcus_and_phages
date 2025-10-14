#!/usr/bin/env bash
#SBATCH --job-name=assemble_phages_flye_canu_raven
#SBATCH --partition=interactive
#SBATCH --cpus-per-task=95
#SBATCH --mem=350G
#SBATCH --time=12:00:00
#SBATCH --output=/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages/slurm_logs/%x/result_%x.%j.log

set -euo pipefail

# =========================
# USER CONFIG
# =========================
WD="${WD:-/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages}"
OUTROOT="${OUTROOT:-$WD/phages_analysis/ont_pipeline_20251009_134551}"   # same root used earlier
PER_PHAGE_RAW="${PER_PHAGE_RAW:-$OUTROOT/demux_trimmed_per_phage}"       # where M*.trimmed.fastq.gz live

# Include only these phages (space-separated). Default targets M1..M5.
PHAGES="${PHAGES:-M1 M2 M3 M4 M5}"

# Assembly parameters
THREADS="${SLURM_CPUS_PER_TASK:-95}"
GENOME_SIZE="${GENOME_SIZE:-200k}"      # for Canu (e.g., 50k..500k for typical phages)
FLYE_READ_MODE="${FLYE_READ_MODE:---nano-raw}"  # --nano-raw or --nano-hq (use hq if you’re confident in Q20+ reads)
FLYE_MIN_OVERLAP="${FLYE_MIN_OVERLAP:-2000}"
FLYE_ASM_COV="${FLYE_ASM_COV:-100}"
MIN_CONTIG="${MIN_CONTIG:-1000}"         # for QUAST reporting

# Tool envs/paths
# Flye
FLYE_ACTIVATE="${FLYE_ACTIVATE:-/vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh}"
FLYE_ENV="${FLYE_ENV:-flye_v2.9.2}"

# Canu (binary)
CANU_BIN="${CANU_BIN:-/home/groups/VEO/tools/canu/v2.2/canu-2.2/bin/canu}"

# Raven
RAVEN_ACTIVATE="${RAVEN_ACTIVATE:-/vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh}"
RAVEN_ENV="${RAVEN_ENV:-raven_v1.8.3}"

# QUAST
QUAST_CMD="${QUAST_CMD:-python3 /home/groups/VEO/tools/quast/v5.2.0/quast.py}"

# =========================
# STEP CONTROL
# =========================
# Steps: flye, canu, raven, quast
ONLY_STEPS="${ONLY_STEPS:-}"   # e.g. "flye,quast"
SKIP_STEPS="${SKIP_STEPS:-}"   # e.g. "canu"

_contains_step() { local needle="$1" hay="$2"; [[ ",$hay," == *",$needle,"* ]]; }
should_run() {
  local step="$1"
  if [[ -n "$ONLY_STEPS" ]] && ! _contains_step "$step" "$ONLY_STEPS"; then return 1; fi
  if _contains_step "$step" "$SKIP_STEPS"; then return 1; fi
  return 0
}

# =========================
# OUTPUT DIRS
# =========================
FLYE_ROOT="${OUTROOT}/phage_assembly_flye"
CANU_ROOT="${OUTROOT}/phage_assembly_canu"
RAVEN_ROOT="${OUTROOT}/phage_assembly_raven"
QUAST_ROOT="${OUTROOT}/phage_assembly_quast"

mkdir -p "$FLYE_ROOT" "$CANU_ROOT" "$RAVEN_ROOT" "$QUAST_ROOT"

echo "[INFO] OUTROOT: $OUTROOT"
echo "[INFO] PER_PHAGE_RAW: $PER_PHAGE_RAW"
echo "[INFO] PHAGES: $PHAGES"
echo "[INFO] GENOME_SIZE: $GENOME_SIZE"
echo "[INFO] THREADS: $THREADS"
echo "[INFO] Steps: ONLY_STEPS='${ONLY_STEPS}' SKIP_STEPS='${SKIP_STEPS}'"

# =========================
# HELPERS
# =========================
activate_flye() { source "$FLYE_ACTIVATE"; conda activate "$FLYE_ENV"; }
activate_raven() { source "$RAVEN_ACTIVATE"; conda activate "$RAVEN_ENV"; }

find_fastq_for_phage() {
  local ph="$1"
  # prefer exact Mx.trimmed.fastq.gz; fallback to glob if needed
  local fq="${PER_PHAGE_RAW}/${ph}.trimmed.fastq.gz"
  if [[ -s "$fq" ]]; then
    printf "%s" "$fq"
    return 0
  fi
  # fallback: any Mx*.trimmed.fastq.gz (shouldn’t usually be needed)
  local first
  first="$(ls -1 "${PER_PHAGE_RAW}/${ph}"*.trimmed.fastq.gz 2>/dev/null | head -n1 || true)"
  if [[ -n "$first" ]]; then
    printf "%s" "$first"
    return 0
  fi
  return 1
}

# =========================
# MAIN
# =========================
for ph in $PHAGES; do
  echo "=========================="
  echo "[PHAGE] $ph"
  fq="$(find_fastq_for_phage "$ph" || true)"
  if [[ -z "${fq:-}" ]]; then
    echo "[WARN] Cannot find trimmed FASTQ for $ph under ${PER_PHAGE_RAW}. Skipping."
    continue
  fi
  echo "[INFO] Reads: $fq"

  # ---------- Flye ----------
  if should_run flye; then
    outdir="${FLYE_ROOT}/${ph}"
    flye_ctg="${outdir}/assembly.fasta"
    if [[ -s "$flye_ctg" ]]; then
      echo "[SKIP] Flye exists: $flye_ctg"
    else
      echo "[RUN ] Flye → $outdir"
      mkdir -p "$outdir"
      activate_flye
      flye $FLYE_READ_MODE "$fq" --out-dir "$outdir" --threads "$THREADS" --min-overlap "$FLYE_MIN_OVERLAP" --asm-coverage "$FLYE_ASM_COV" || {
        echo "[WARN] Flye failed for $ph"; }
    fi
  else
    echo "[SKIP] Step flye (per switch)"
  fi

  # ---------- Canu ----------
  if should_run canu; then
    outdir="${CANU_ROOT}/${ph}"
    canu_prefix="${ph}"
    canu_ctg="${outdir}/${canu_prefix}.contigs.fasta"
    if [[ -s "$canu_ctg" ]]; then
      echo "[SKIP] Canu exists: $canu_ctg"
    else
      echo "[RUN ] Canu → $outdir"
      mkdir -p "$outdir"
      # Canu expects uncompressed or compressed fastq is fine; it handles .gz
      "$CANU_BIN" \
        -p "$canu_prefix" -d "$outdir" \
        genomeSize="$GENOME_SIZE" \
        useGrid=false \
        maxThreads="$THREADS" \
        -nanopore "$fq" || {
          echo "[WARN] Canu failed for $ph"; }
    fi
  else
    echo "[SKIP] Step canu (per switch)"
  fi

  # ---------- Raven ----------
  if should_run raven; then
    outdir="${RAVEN_ROOT}/${ph}"
    raven_ctg="${outdir}/raven.fasta"
    if [[ -s "$raven_ctg" ]]; then
      echo "[SKIP] Raven exists: $raven_ctg"
    else
      echo "[RUN ] Raven → $outdir"
      mkdir -p "$outdir"
      activate_raven
      # Raven can write to stdout, we capture to file
      raven -t "$THREADS" "$fq" > "${raven_ctg}.partial" 2> "${outdir}/raven.log" || {
        echo "[WARN] Raven failed for $ph"; rm -f "${raven_ctg}.partial"; }
      [[ -s "${raven_ctg}.partial" ]] && mv "${raven_ctg}.partial" "$raven_ctg"
    fi
  else
    echo "[SKIP] Step raven (per switch)"
  fi

  # ---------- QUAST ----------
  if should_run quast; then
    echo "[RUN ] QUAST for $ph"
    qout="${QUAST_ROOT}/${ph}"
    mkdir -p "$qout"

    # Collect assemblies that exist
    flye_ctg="${FLYE_ROOT}/${ph}/assembly.fasta"
    canu_ctg="${CANU_ROOT}/${ph}/${ph}.contigs.fasta"
    raven_ctg="${RAVEN_ROOT}/${ph}/raven.fasta"

    declare -a asm_paths=()
    declare -a asm_labels=()

    if [[ -s "$flye_ctg" ]]; then asm_paths+=("$flye_ctg"); asm_labels+=("flye"); fi
    if [[ -s "$canu_ctg" ]]; then asm_paths+=("$canu_ctg"); asm_labels+=("canu"); fi
    if [[ -s "$raven_ctg" ]]; then asm_paths+=("$raven_ctg"); asm_labels+=("raven"); fi

    if (( ${#asm_paths[@]} == 0 )); then
      echo "[WARN] No assemblies found for $ph; skipping QUAST."
    else
      # join labels with commas
      IFS=, read -r -a labels_csv <<< "${asm_labels[*]/%/,}"
      labels_joined="$(printf "%s," "${asm_labels[@]}")"; labels_joined="${labels_joined%,}"

      echo "[INFO] QUAST assemblies: ${asm_paths[*]}"
      $QUAST_CMD \
        --threads "$THREADS" \
        --min-contig "$MIN_CONTIG" \
        --labels "$labels_joined" \
        -o "$qout" \
        "${asm_paths[@]}" || {
          echo "[WARN] QUAST failed for $ph"; }
    fi
  else
    echo "[SKIP] Step quast (per switch)"
  fi

done

echo "[DONE] All requested steps completed."
echo "Outputs:"
echo "  Flye : ${FLYE_ROOT}/<phage>/assembly.fasta"
echo "  Canu : ${CANU_ROOT}/<phage>/<phage>.contigs.fasta"
echo "  Raven: ${RAVEN_ROOT}/<phage>/raven.fasta"
echo "  QUAST: ${QUAST_ROOT}/<phage>/report.html"
