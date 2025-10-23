#!/usr/bin/env bash
#SBATCH --job-name=trim_filter_qc_reads_per_phage
#SBATCH --partition=short
#SBATCH --cpus-per-task=80
#SBATCH --mem=80G
#SBATCH --time=3:00:00
#SBATCH --output=/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages/slurm_logs/%x/result_%x.%j.log

set -euo pipefail

# --- user inputs (edit or pass with sbatch --export=ALL,VAR=VAL) ---
WD="/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages"
IN_DIR="${IN_DIR:-${WD}/phages_analysis/ont_basecalls_20251007_160801/demux_no_trim_per_phage_raw}"
OUT_ROOT="${OUT_ROOT:-${IN_DIR%/}/..}"  # parent of IN_DIR
KIT="${KIT:-SQK-NBD114-96}"
DORADO="${DORADO:-/home/groups/VEO/tools/dorado/v1.1.0/bin/dorado}"

# filtering options (set FILTER=none to skip)
FILTER="${FILTER:-filtlong}"        # filtlong|chopper|none
MINLEN="${MINLEN:-1000}"            # min length for filtering
KEEPPC="${KEEPPC:-95}"              # filtlong keep percent
QCHOP="${QCHOP:-10}"                # chopper q cutoff

# QC options
RUN_NANOPLOT="${RUN_NANOPLOT:-1}"         # 1 to enable NanoPlot on trimmed
RUN_FASTQC_TRIM="${RUN_FASTQC_TRIM:-1}"   # 1 to FastQC+MultiQC trimmed
RUN_FASTQC_FILTER="${RUN_FASTQC_FILTER:-1}" # 1 to FastQC+MultiQC filtered
NANOPLOT_THREADS="${SLURM_CPUS_PER_TASK:-80}"
THREADS="${SLURM_CPUS_PER_TASK:-80}"

# Tool paths (your site)
FastQC="${FastQC:-/home/groups/VEO/tools/fastqc/v0.12.1/fastqc}"
MultiQC_ACTIVATE="${MultiQC_ACTIVATE:-/home/groups/VEO/tools/multiQC/v1.15/bin/activate}"

# --- env helpers ---
if command -v pigz >/dev/null 2>&1; then
  COMPRESS="pigz -c -p ${SLURM_CPUS_PER_TASK:-8}"
  DECOMP="pigz -dc"
else
  COMPRESS="gzip -c"
  DECOMP="gzip -dc"
fi

TRIM_DIR="${OUT_ROOT}/per_phage_trimmed"
FILT_DIR="${OUT_ROOT}/per_phage_filtered"
NANO_DIR="${OUT_ROOT}/per_phage_nanoplot_trimmed"

# FastQC/MultiQC outputs
FQCTRIM_DIR="${OUT_ROOT}/QC_fastqc_trimmed"
FQCFILT_DIR="${OUT_ROOT}/QC_fastqc_filtered"

mkdir -p "$TRIM_DIR" "$FILT_DIR"

echo "[INFO] IN_DIR : $IN_DIR"
echo "[INFO] TRIM_DIR -> $TRIM_DIR"
echo "[INFO] FILT_DIR -> $FILT_DIR"
echo "[INFO] KIT    : $KIT"
echo "[INFO] FILTER : $FILTER (minlen=$MINLEN keep%=$KEEPPC q=$QCHOP)"

shopt -s nullglob
mapfile -t PHAGE_FASTQS < <(find "$IN_DIR" -maxdepth 1 -type f -name "M*.raw.fastq.gz" | sort)
(( ${#PHAGE_FASTQS[@]} > 0 )) || { echo "[ERR] No M*.raw.fastq.gz in $IN_DIR"; exit 2; }

# --- loop phages: trim -> (optional) filter -> (optional) NanoPlot ---

for IN in "${PHAGE_FASTQS[@]}"; do
  base="$(basename "$IN" .raw.fastq.gz)"   # e.g. M1
  OUT_TRIM="$TRIM_DIR/${base}.trimmed.fastq.gz"
  OUT_FILT="$FILT_DIR/${base}.filtered.fastq.gz"

  # 1) TRIM (adapters/primers/barcodes)
  source /home/zo49sog/miniconda3/etc/profile.d/conda.sh && conda activate nanopore_analysis

  if [[ -s "$OUT_TRIM" ]]; then
    echo "[SKIP] Trim exists: $OUT_TRIM"
  else
    echo "[RUN ] Trim $base -> $OUT_TRIM"
    $DECOMP "$IN" \
      | "$DORADO" trim - --sequencing-kit "$KIT" --emit-fastq \
      | $COMPRESS > "${OUT_TRIM}.partial"
    mv "${OUT_TRIM}.partial" "$OUT_TRIM"
    echo "[DONE] Trim $base"
  fi

  # 2) FILTER (optional)
  if [[ "$FILTER" == "none" ]]; then
    echo "[INFO] Skipping filter for $base"
  else
    if [[ -s "$OUT_FILT" ]]; then
      echo "[SKIP] Filter exists: $OUT_FILT"
    else
      echo "[RUN ] Filter $base ($FILTER) -> $OUT_FILT"
      case "$FILTER" in
        filtlong)
          filtlong --min_length "$MINLEN" --keep_percent "$KEEPPC" "$OUT_TRIM" \
            | $COMPRESS > "${OUT_FILT}.partial"
          ;;
        chopper)
          $DECOMP "$OUT_TRIM" | chopper -q "$QCHOP" \
            | $COMPRESS > "${OUT_FILT}.partial"
          ;;
        *)
          echo "[ERR] Unknown FILTER=$FILTER"; exit 1;;
      esac
      mv "${OUT_FILT}.partial" "$OUT_FILT"
      echo "[DONE] Filter $base"
    fi
  fi

  # 3) NanoPlot on trimmed (optional quick QC)
  source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh && conda activate nanoplot_v1.41.3

  if [[ "$RUN_NANOPLOT" == "1" ]]; then
    mkdir -p "$NANO_DIR"
    OUT_NANO="$NANO_DIR/${base}"
    if [[ -s "$OUT_NANO/NanoStats.txt" ]]; then
      echo "[SKIP] NanoPlot exists for $base"
    else
      echo "[RUN ] NanoPlot trimmed $base"

      NanoPlot --fastq "$OUT_TRIM" -o "$OUT_NANO" -t "$NANOPLOT_THREADS" --tsv_stats --loglength || true
    fi
  fi
done

# -------------------- FastQC + MultiQC on TRIMMED --------------------
if [[ "$RUN_FASTQC_TRIM" == "1" ]]; then
  echo "[INFO] FastQC on trimmed -> $FQCTRIM_DIR"
  mkdir -p "$FQCTRIM_DIR"
  mapfile -t TRIM_FILES < <(find "$TRIM_DIR" -maxdepth 1 -type f -name "*.fastq.gz" | sort)
  if (( ${#TRIM_FILES[@]} > 0 )); then
    if [[ -x "$FastQC" ]]; then
      "$FastQC" -t "$THREADS" -o "$FQCTRIM_DIR" "${TRIM_FILES[@]}"
    else
      echo "[ERR] FastQC not found at: $FastQC"; exit 1
    fi
    # MultiQC
    if [[ -f "$MultiQC_ACTIVATE" ]]; then
      # shellcheck disable=SC1090
      source "$MultiQC_ACTIVATE"
    fi
    multiqc -o "$FQCTRIM_DIR" "$FQCTRIM_DIR"
    echo "[INFO] Trimmed MultiQC: $FQCTRIM_DIR/multiqc_report.html"
  else
    echo "[WARN] No trimmed files found for FastQC."
  fi
fi

# -------------------- FastQC + MultiQC on FILTERED -------------------
if [[ "$RUN_FASTQC_FILTER" == "1" ]]; then
  echo "[INFO] FastQC on filtered -> $FQCFILT_DIR"
  mkdir -p "$FQCFILT_DIR"
  mapfile -t FILT_FILES < <(find "$FILT_DIR" -maxdepth 1 -type f -name "*.fastq.gz" | sort)
  if (( ${#FILT_FILES[@]} > 0 )); then
    if [[ -x "$FastQC" ]]; then
      "$FastQC" -t "$THREADS" -o "$FQCFILT_DIR" "${FILT_FILES[@]}"
    else
      echo "[ERR] FastQC not found at: $FastQC"; exit 1
    fi
    # MultiQC
    if [[ -f "$MultiQC_ACTIVATE" ]]; then
      # shellcheck disable=SC1090
      source "$MultiQC_ACTIVATE"
    fi
    multiqc -o "$FQCFILT_DIR" "$FQCFILT_DIR"
    echo "[INFO] Filtered MultiQC: $FQCFILT_DIR/multiqc_report.html"
  else
    echo "[WARN] No filtered files found for FastQC."
  fi
fi

echo "[INFO] Done. Trimmed: $TRIM_DIR ; Filtered: $FILT_DIR"
