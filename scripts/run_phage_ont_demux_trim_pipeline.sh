#!/usr/bin/env bash
#SBATCH --job-name=phage_ont_demux_trim_filter_qc_pipeline
#SBATCH --partition=interactive
#SBATCH --cpus-per-task=80
#SBATCH --mem=90G
#SBATCH --time=12:00:00
#SBATCH --output=/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages/slurm_logs/%x/result_%x.%j.log

set -euo pipefail

# =========================
# USER CONFIG
# =========================
WD="${WD:-/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages}"

# You already have basecalled reads (no GPU needed here):
IN=${IN:-"${WD}/phages_analysis/ont_basecalls_20251007_160801/all.fastq.gz"}

# OUTROOT can be fixed or timestamped
#OUTROOT="${OUTROOT:-$WD/phages_analysis/ont_pipeline_$(date +%Y%m%d_%H%M%S)}"
OUTROOT="${OUTROOT:-$WD/phages_analysis/ont_pipeline_20251009_134551}"

# Dorado demux settings (demux doesn't need a basecalling model, but we keep a placeholder)
MODEL="${MODEL:-sup}"
KIT="${KIT:-SQK-NBD114-96}"
BARCODE_BOTH_ENDS="${BARCODE_BOTH_ENDS:-0}"

# Phage map TSV (phage<TAB>barcode01,barcode02)
PHAGE_MAP="${PHAGE_MAP:-$WD/data/phage_barcode_map.tsv}"

# Binaries / envs
DORADO="${DORADO:-/home/groups/VEO/tools/dorado/v1.1.0/bin/dorado}"
FASTQC_BIN="${FASTQC_BIN:-/home/groups/VEO/tools/fastqc/v0.12.1/fastqc}"
MULTIQC_ACTIVATE="${MULTIQC_ACTIVATE:-/home/groups/VEO/tools/multiQC/v1.15/bin/activate}"
NANO_ENV_ACTIVATE="${NANO_ENV_ACTIVATE:-/vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh}"
NANOPLOT_ENV="${NANOPLOT_ENV:-nanoplot_v1.41.3}"
ONT_ENV_ACTIVATE="${ONT_ENV_ACTIVATE:-/home/zo49sog/miniconda3/etc/profile.d/conda.sh}"
ONT_ENV_NAME="${ONT_ENV_NAME:-nanopore_analysis}"   # should contain porechop_abi, filtlong, chopper, etc.

# Filtering options
FILTER="${FILTER:-filtlong}"   # filtlong|chopper|none
MINLEN="${MINLEN:-1000}"
KEEPPC="${KEEPPC:-95}"
QCHOP="${QCHOP:-10}"

# QC toggles
RUN_NANOPLOT="${RUN_NANOPLOT:-1}"
RUN_FASTQC_TRIM="${RUN_FASTQC_TRIM:-1}"
RUN_FASTQC_FILT="${RUN_FASTQC_FILT:-1}"

# Downstream set selection: porechop | trimmed
DOWNSTREAM_SET="${DOWNSTREAM_SET:-porechop}"

THREADS="${SLURM_CPUS_PER_TASK:-80}"

# =========================
# OUTPUT PATHS
# =========================
mkdir -p "$OUTROOT"
ALL_FASTQ="${ALL_FASTQ:-$IN}"

DEMUX_DIR="$OUTROOT/demux_trimmed"                 # per-barcode, trimmed by dorado demux
PER_PHAGE_RAW="$OUTROOT/demux_trimmed_per_phage"   # concatenated per-phage (already dorado-trimmed)
PORECHOP_DIR="$OUTROOT/per_phage_porechop"         # per-phage porechop_abi outputs + logs
PER_PHAGE_FILT="$OUTROOT/per_phage_filtered"
NANOPLOT_DIR="$OUTROOT/per_phage_nanoplot"
FQCTRIM_DIR="$OUTROOT/QC_fastqc_trimmed"
FQCFILT_DIR="$OUTROOT/QC_fastqc_filtered"

mkdir -p "$DEMUX_DIR" "$PER_PHAGE_RAW" "$PORECHOP_DIR" "$PER_PHAGE_FILT"

# compressors
if command -v pigz >/dev/null 2>&1; then
  COMPRESS="pigz -c -p ${THREADS}"
  DECOMP="pigz -dc"
else
  COMPRESS="gzip -c"
  DECOMP="gzip -dc"
fi

echo "[INFO] Job: $SLURM_JOB_NAME id=$SLURM_JOB_ID"
echo "[INFO] Using basecalled reads: $ALL_FASTQ"
echo "[INFO] OUTROOT:   $OUTROOT"
echo "[INFO] KIT:       $KIT   (BARCODE_BOTH_ENDS=$BARCODE_BOTH_ENDS)"
echo "[INFO] PHAGE_MAP: $PHAGE_MAP"

# small FASTQ check (gzip OK + FASTQ record structure looks sane)
is_fastq_gz () {
  local f="$1"
  # 1) gzip integrity
  gzip -t "$f" 2>/dev/null || return 1
  # 2) 4-line block shape and markers (@ on 1st, + on 3rd lines)
  zcat -f "$f" 2>/dev/null | awk '
    NR%4==1 && $0 !~ /^@/ {exit 1}
    NR%4==3 && $0 !~ /^\+/ {exit 1}
    END{ if (NR%4!=0) exit 1 }'
}

# # =========================
# # STEP 2: Demultiplex WITH trimming → per-barcode FASTQ
# # =========================
# if compgen -G "$DEMUX_DIR/*barcode*.fastq" > /dev/null || compgen -G "$DEMUX_DIR/*barcode*.fastq.gz" > /dev/null; then
#   echo "[SKIP] Demux-trim outputs present in $DEMUX_DIR"
# else
#   echo "[RUN ] Demux (TRIM ON) → $DEMUX_DIR"
#   DEMUX_CMD=( "$DORADO" demux --kit-name "$KIT" --emit-fastq "$ALL_FASTQ" --output-dir "$DEMUX_DIR" --emit-summary )
#   if [[ "$BARCODE_BOTH_ENDS" == "1" ]]; then
#     DEMUX_CMD+=( --barcode-both-ends )
#   fi
#   "${DEMUX_CMD[@]}"
#   echo "[DONE] Demux (trimmed)."
# fi

# # =========================
# # STEP 3: CONCAT per-phage (trimmed per-barcode → M*.trimmed.fastq.gz)
# # =========================
# echo "[RUN ] Concatenate per-phage → $PER_PHAGE_RAW"
# while IFS=$'\t' read -r phage barcs; do
#   # skip blanks / header / comments
#   [[ -z "${phage:-}" || -z "${barcs:-}" ]] && continue
#   [[ "${phage}" == "phage" || "${phage:0:1}" == "#" ]] && continue

#   mapfile -t files < <(
#     IFS=',' read -ra arr <<< "$barcs"
#     for bc in "${arr[@]}"; do
#       bc="${bc//[[:space:]]/}"
#       find "$DEMUX_DIR" -maxdepth 1 -type f \
#         \( -name "*_${bc}.fastq.gz" -o -name "*_${bc}.fastq" -o -name "${bc}.fastq.gz" -o -name "${bc}.fastq" \)
#     done | sort
#   )

#   echo "[DEBUG] $phage barcodes: $barcs"
#   printf '[DEBUG]   %s\n' "${files[@]:-<none>}"

#   if (( ${#files[@]} == 0 )); then
#     echo "[WARN] No per-barcode FASTQs for $phage ($barcs). Skipping."
#     continue
#   fi

#   out="$PER_PHAGE_RAW/${phage}.trimmed.fastq.gz"
#   if [[ -s "$out" ]]; then
#     echo "[SKIP] $phage exists: $out"
#     continue
#   fi

#   echo "[RUN ] $phage: concat ${#files[@]} → $out"
#   {
#     for f in "${files[@]}"; do
#       case "$f" in
#         *.fastq.gz) $DECOMP "$f" ;;
#         *.fastq)    cat "$f" ;;
#         *)          echo "[WARN] Skip unknown extension: $f" >&2 ;;
#       esac
#     done
#   } | $COMPRESS > "$out"

#   if ! is_fastq_gz "$out"; then
#     echo "[WARN] FASTQ structure check failed for $phage: $out (continuing)"
#   fi
# done < "$PHAGE_MAP"
# echo "[DONE] Per-phage concatenation (trimmed)."

# =========================
# STEP 3.5: Porechop-ABI (adapter/artefact discovery + trimming) per-phage
# =========================
echo "[RUN ] porechop_abi per-phage → $PORECHOP_DIR (logs per phage)"

# use the env that has porechop_abi installed
source "$ONT_ENV_ACTIVATE"; conda activate "$ONT_ENV_NAME"

# make loop robust if no files match
shopt -s nullglob

found_any=0
for INPHAGE in "$PER_PHAGE_RAW"/M*.trimmed.fastq.gz; do
  found_any=1
  [[ -f "$INPHAGE" ]] || continue

  # derive phage name, e.g. M1 from M1.trimmed.fastq.gz
  bn="$(basename "$INPHAGE")"
  ph="${bn%.trimmed.fastq.gz}"

  if [[ ! -s "$INPHAGE" ]]; then
    echo "[WARN] Input is missing or empty: $(basename "$INPHAGE")"
    continue
  fi

  ph_outdir="$PORECHOP_DIR/$ph"
  mkdir -p "$ph_outdir"

  OUT_FASTQ="$ph_outdir/${ph}.porechop.fastq.gz"
  LOG="$ph_outdir/porechop_abi.log"

  if [[ -s "$OUT_FASTQ" ]]; then
    echo "[SKIP] porechop_abi exists for $ph: $OUT_FASTQ"
    continue
  fi

  # quieter logs
  PORECHOP_VERBOSITY="${PORECHOP_VERBOSITY:-1}"

  echo "[RUN ] porechop_abi $ph"
  porechop_abi \
    -i "$INPHAGE" \
    -o "$OUT_FASTQ" \
    --format fastq.gz \
    -t "$THREADS" \
    -v "$PORECHOP_VERBOSITY" \
    --min_split_read_size 2000 \
    >"$LOG" 2>&1 || { echo "[WARN] porechop_abi failed for $ph (see $LOG)"; continue; }

  # verify gz integrity & FASTQ structure (uses your is_fastq_gz helper)
  if ! is_fastq_gz "$OUT_FASTQ"; then
    echo "[WARN] porechop_abi output failed FASTQ check for $ph: $OUT_FASTQ (continuing)"
  fi

  # --- concise stats -> $ph_outdir/porechop_abi.summary.tsv and a global _summary.tsv ---
  SUM="$ph_outdir/porechop_abi.summary.tsv"
  {
    # extract totals like: "16,494 / 344,434 reads had adapters trimmed from their start (153,974 bp removed)"
    tot_reads="$(grep -Eo ' / [0-9,]+ reads had adapters trimmed from their start' "$LOG" | head -n1 | awk '{gsub(",","",$3); print $2}' | tr -d ' /')"
    start_line="$(grep -E 'reads had adapters trimmed from their start ' "$LOG" | head -n1)"
    end_line="$(grep -E 'reads had adapters trimmed from their end ' "$LOG"   | head -n1)"

    start_n="$(printf '%s\n' "$start_line" | awk '{gsub(",","",$1); print $1}')"
    start_bp="$(printf  '%s\n' "$start_line" | sed -E 's/.*\(([0-9,]+) bp removed\).*/\1/' | tr -d ,)"
    end_n="$(printf   '%s\n' "$end_line"   | awk '{gsub(",","",$1); print $1}')"
    end_bp="$(printf  '%s\n' "$end_line"   | sed -E 's/.*\(([0-9,]+) bp removed\).*/\1/' | tr -d ,)"

    # crude estimate of split events = lines showing a detected middle adapter (those “read coords:” examples)
    splits="$(grep -c 'read coords:' "$LOG" || true)"

    # barcode tallies seen in the log (e.g. BC01, BC01_rev, etc.)
    # format as counts like "BC01:123;BC01_rev:45"
    bc_counts="$(grep -Eo 'BC[0-9]{2}(_rev)?' "$LOG" \
        | sort | uniq -c \
        | awk '{printf("%s:%s;", $2,$1)}' \
        | sed 's/;$//')"

    # header + row (one row file per phage)
    echo -e "phage\ttotal_reads\tstart_trimmed\tstart_bp_removed\tend_trimmed\tend_bp_removed\tsplit_events\tbarcode_counts"
    echo -e "${ph}\t${tot_reads:-NA}\t${start_n:-0}\t${start_bp:-0}\t${end_n:-0}\t${end_bp:-0}\t${splits:-0}\t${bc_counts:-NA}"
  } > "$SUM"

  # also append to a run-level summary
  GLOBAL_SUM="$PORECHOP_DIR/_summary.tsv"
  if [[ ! -s "$GLOBAL_SUM" ]]; then
    head -n1 "$SUM" > "$GLOBAL_SUM"
  fi
  tail -n1 "$SUM" >> "$GLOBAL_SUM"

done
shopt -u nullglob

echo "[DONE] porechop_abi."


# =========================
# SELECT DOWNSTREAM INPUTS (for filtering + "trimmed" QC)
# =========================
shopt -s nullglob
declare -a SRC_FOR_FILTER_LIST TRIMMED_FOR_QC_LIST
if [[ "$DOWNSTREAM_SET" == "porechop" ]]; then
  SRC_FOR_FILTER_LIST=( "$PORECHOP_DIR"/M*/M*.porechop.fastq.gz )
  TRIMMED_FOR_QC_LIST=( "${SRC_FOR_FILTER_LIST[@]}" )
  echo "[INFO] DOWNSTREAM_SET=porechop → use porechop outputs."
else
  SRC_FOR_FILTER_LIST=( "$PER_PHAGE_RAW"/M*.trimmed.fastq.gz )
  TRIMMED_FOR_QC_LIST=( "${SRC_FOR_FILTER_LIST[@]}" )
  echo "[INFO] DOWNSTREAM_SET=trimmed → use dorado-trimmed concatenations."
fi
if (( ${#SRC_FOR_FILTER_LIST[@]} == 0 )); then
  echo "[WARN] No inputs found for downstream in set: $DOWNSTREAM_SET"
fi

# =========================
# STEP 4: (OPTIONAL) FILTER per-phage → per_phage_filtered
# =========================
if [[ "$FILTER" == "none" ]]; then
  echo "[INFO] Skipping filtering (FILTER=none)"
else
  echo "[RUN ] Filtering ($FILTER) → $PER_PHAGE_FILT (from set: $DOWNSTREAM_SET)"
  source "$ONT_ENV_ACTIVATE"; conda activate "$ONT_ENV_NAME"
  mkdir -p "$PER_PHAGE_FILT"

  for IN in "${SRC_FOR_FILTER_LIST[@]}"; do
    [[ -f "$IN" ]] || continue
    fname="$(basename "$IN")"
    case "$fname" in
      *.porechop.fastq.gz) base="${fname%.porechop.fastq.gz}";;
      *.trimmed.fastq.gz)  base="${fname%.trimmed.fastq.gz}";;
      *)                   base="${fname%.fastq.gz}";;
    esac
    OUT_FILT="$PER_PHAGE_FILT/${base}.filtered.fastq.gz"
    [[ -s "$OUT_FILT" ]] && { echo "[SKIP] Filter exists: $OUT_FILT"; continue; }

    case "$FILTER" in
      filtlong)
        filtlong --min_length "$MINLEN" --keep_percent "$KEEPPC" "$IN" \
          | $COMPRESS > "${OUT_FILT}.partial"
        ;;
      chopper)
        $DECOMP "$IN" | chopper -q "$QCHOP" | $COMPRESS > "${OUT_FILT}.partial"
        ;;
      *) echo "[ERR] Unknown FILTER=$FILTER"; exit 1 ;;
    esac
    mv "${OUT_FILT}.partial" "$OUT_FILT"
    echo "[DONE] Filtered: $OUT_FILT"
  done
  echo "[DONE] Filtering."
fi

# =========================
# STEP 5: QC — NanoPlot (on selected 'trimmed' set); FastQC+MultiQC (trimmed & filtered)
# =========================
# NanoPlot on selected 'trimmed' set
if [[ "$RUN_NANOPLOT" == "1" ]]; then
  mkdir -p "$NANOPLOT_DIR"
  if [[ -f "$NANO_ENV_ACTIVATE" ]]; then source "$NANO_ENV_ACTIVATE"; conda activate "$NANOPLOT_ENV" || true; fi
  for IN in "${TRIMMED_FOR_QC_LIST[@]}"; do
    [[ -f "$IN" ]] || continue
    fname="$(basename "$IN")"
    case "$fname" in
      *.porechop.fastq.gz) base="${fname%.porechop.fastq.gz}";;
      *.trimmed.fastq.gz)  base="${fname%.trimmed.fastq.gz}";;
      *)                   base="${fname%.fastq.gz}";;
    esac
    OUT_NANO="$NANOPLOT_DIR/${base}"
    [[ -s "$OUT_NANO/NanoStats.txt" ]] && { echo "[SKIP] NanoPlot exists: $base"; continue; }
    echo "[RUN ] NanoPlot ($DOWNSTREAM_SET) $base"
    NanoPlot --fastq "$IN" -o "$OUT_NANO" -t "$THREADS" --tsv_stats --loglength || true
  done
fi

# FastQC/MultiQC on selected 'trimmed' set
if [[ "$RUN_FASTQC_TRIM" == "1" ]]; then
  mkdir -p "$FQCTRIM_DIR"
  if (( ${#TRIMMED_FOR_QC_LIST[@]} > 0 )); then
    "$FASTQC_BIN" -t "$THREADS" -o "$FQCTRIM_DIR" "${TRIMMED_FOR_QC_LIST[@]}"
    if [[ -f "$MULTIQC_ACTIVATE" ]]; then source "$MULTIQC_ACTIVATE"; fi
    multiqc -o "$FQCTRIM_DIR" "$FQCTRIM_DIR"
  else
    echo "[WARN] No inputs found for FastQC (trimmed set: $DOWNSTREAM_SET)."
  fi
fi

# FastQC/MultiQC on filtered
if [[ "$RUN_FASTQC_FILT" == "1" ]]; then
  mkdir -p "$FQCFILT_DIR"
  mapfile -t FILT_FILES < <(find "$PER_PHAGE_FILT" -maxdepth 1 -type f -name "*.fastq.gz" | sort)
  if (( ${#FILT_FILES[@]} > 0 )); then
    "$FASTQC_BIN" -t "$THREADS" -o "$FQCFILT_DIR" "${FILT_FILES[@]}"
    if [[ -f "$MULTIQC_ACTIVATE" ]]; then source "$MULTIQC_ACTIVATE"; fi
    multiqc -o "$FQCFILT_DIR" "$FQCFILT_DIR"
  else
    echo "[INFO] No filtered FASTQs found for FastQC (maybe FILTER=none)."
  fi
fi

# =========================
# DONE
# =========================
echo "[INFO] Outputs:"
echo "  Demux (trimmed)         : $DEMUX_DIR/*_barcodeXX.fastq"
echo "  Per-phage (trimmed)     : $PER_PHAGE_RAW/M*.trimmed.fastq.gz"
echo "  Porechop per-phage      : $PORECHOP_DIR/M*/M*.porechop.fastq.gz (+ logs)"
echo "  Per-phage filtered      : $PER_PHAGE_FILT/M*.filtered.fastq.gz (if enabled)"
echo "  NanoPlot (selected set) : $NANOPLOT_DIR/*/NanoPlot-report.html"
echo "  FastQC+MultiQC          : $FQCTRIM_DIR , $FQCFILT_DIR"
echo "  DOWNSTREAM_SET          : $DOWNSTREAM_SET (porechop|trimmed)"
