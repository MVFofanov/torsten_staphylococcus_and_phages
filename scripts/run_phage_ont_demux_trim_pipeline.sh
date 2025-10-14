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

# Porechop verbosity
PORECHOP_VERBOSITY="${PORECHOP_VERBOSITY:-1}"

# =========================
# STEP CONTROL
# =========================
# Run only these steps (comma-separated) or leave empty to allow all:
ONLY_STEPS="${ONLY_STEPS:-}"          # e.g. "filter,qc_trim,qc_filt"
# Skip these steps regardless (comma-separated):
SKIP_STEPS="${SKIP_STEPS:-}"          # e.g. "demux,concat,porechop"

# Canonical step names:
# demux, concat, porechop, filter, nanoplot, qc_trim, qc_filt, summarize
_contains_step() {  # _contains_step "name" "list,a,b"
  local needle="$1"; local hay="$2"
  [[ ",$hay," == *",$needle,"* ]]
}
should_run() {
  local step="$1"
  if [[ -n "$ONLY_STEPS" ]] && ! _contains_step "$step" "$ONLY_STEPS"; then
    return 1
  fi
  if _contains_step "$step" "$SKIP_STEPS"; then
    return 1
  fi
  return 0
}

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

echo "[INFO] Job: ${SLURM_JOB_NAME:-NA} id=${SLURM_JOB_ID:-NA}"
echo "[INFO] Using basecalled reads: $ALL_FASTQ"
echo "[INFO] OUTROOT:   $OUTROOT"
echo "[INFO] KIT:       $KIT   (BARCODE_BOTH_ENDS=$BARCODE_BOTH_ENDS)"
echo "[INFO] PHAGE_MAP: $PHAGE_MAP"

# =========================
# HELPERS
# =========================
is_fastq_gz () {
  local f="$1"
  gzip -t "$f" 2>/dev/null || return 1
  zcat -f "$f" 2>/dev/null | awk '
    NR%4==1 && $0 !~ /^@/ {exit 1}
    NR%4==3 && $0 !~ /^\+/ {exit 1}
    END{ if (NR%4!=0) exit 1 }'
}

# SAFE DEFAULTS for maps (avoid nounset issues)
declare -A B2P=()     # barcode -> phage
declare -A PH_SEEN=() # phage -> 1

# =========================
# MAPPING BUILDERS
# =========================
build_barcode_to_phage() {
  B2P=()
  while IFS=$'\t' read -r phage barcs; do
    [[ -z "${phage:-}" || -z "${barcs:-}" ]] && continue
    [[ "$phage" == "phage" || "${phage:0:1}" == "#" ]] && continue
    IFS=',' read -ra arr <<< "$barcs"
    for bc in "${arr[@]}"; do
      bc="$(echo "$bc" | tr '[:upper:]' '[:lower:]' | tr -d '[:space:]')"
      if   [[ "$bc" =~ ^barcode[0-9]{2}$ ]];       then norm="$bc"
      elif [[ "$bc" =~ ^bc([0-9]{2})$ ]];          then norm="barcode${BASH_REMATCH[1]}"
      elif [[ "$bc" =~ ^([0-9]{1,2})$ ]];          then norm=$(printf "barcode%02d" "${BASH_REMATCH[1]}")
      else
        echo "[WARN] Unrecognized barcode token in PHAGE_MAP: '$bc' for $phage" >&2
        continue
      fi
      B2P["$norm"]="$phage"
    done
  done < "$PHAGE_MAP"
}

seed_phages_from_map() {
  PH_SEEN=()
  while IFS=$'\t' read -r phage barcs; do
    [[ -z "${phage:-}" || -z "${barcs:-}" ]] && continue
    [[ "$phage" == "phage" || "${phage:0:1}" == "#" ]] && continue
    PH_SEEN["$phage"]=1
  done < "$PHAGE_MAP"
}

seed_phages_from_existing_outputs() {
  while IFS= read -r f; do
    bn="$(basename "$f")"
    if [[ "$bn" =~ ^([^.]*)\.(porechop|trimmed)\.fastq\.gz$ ]]; then
      PH_SEEN["${BASH_REMATCH[1]}"]=1
    fi
  done < <(find "$PORECHOP_DIR" "$PER_PHAGE_RAW" -maxdepth 2 -type f -name "*.fastq.gz" 2>/dev/null || true)
}

# Parse porechop log → one TSV line
_parse_porechop_log() {
  local log="$1" phage="$2" barcode="$3"
  local start_line end_line split_line start_n end_n start_bp end_bp splits tot_reads bc_counts

  start_line="$(grep -E 'reads had adapters trimmed from their start ' "$log" | head -n1 || true)"
  end_line="$(  grep -E 'reads had adapters trimmed from their end '   "$log" | head -n1 || true)"

  start_n="$(printf '%s\n' "$start_line" | awk '{gsub(",","",$1); print $1+0}')"
  end_n="$(  printf '%s\n' "$end_line"   | awk '{gsub(",","",$1); print $1+0}')"

  start_bp="$(printf '%s\n' "$start_line" | sed -E 's/.*\(([0-9,]+) bp removed\).*/\1/;t;d' | tr -d ,)"
  end_bp="$(  printf '%s\n' "$end_line"   | sed -E 's/.*\(([0-9,]+) bp removed\).*/\1/;t;d' | tr -d ,)"

  split_line="$(grep -E '^[[:space:]]*[0-9,]+ / [0-9,]+ reads were split based on middle adapters' "$log" | head -n1 || true)"
  if [[ -n "$split_line" ]]; then
    splits="$(printf '%s\n' "$split_line" | sed -E 's/^[[:space:]]*([0-9,]+) \/ .*/\1/;t;d' | tr -d ,)"
    tot_reads="$(printf '%s\n' "$split_line" | sed -E 's/^[[:space:]]*[0-9,]+ \/ ([0-9,]+) .*/\1/;t;d' | tr -d ,)"
  else
    splits=0
    if [[ -n "$start_line" ]]; then
      tot_reads="$(printf '%s\n' "$start_line" | sed -E 's/.* \/ ([0-9,]+) reads .*/\1/;t;d' | tr -d ,)"
    else
      tot_reads="NA"
    fi
  fi

  bc_counts="$(grep -Eo 'BC[0-9]{2}(_rev)?' "$log" 2>/dev/null \
      | sort | uniq -c \
      | awk '{printf("%s:%s;", $2,$1)}' \
      | sed 's/;$//' )"
  [[ -z "$bc_counts" ]] && bc_counts="NA"

  echo -e "${phage}\t${barcode}\t${tot_reads}\t${start_n:-0}\t${start_bp:-0}\t${end_n:-0}\t${end_bp:-0}\t${splits:-0}\t${bc_counts}"
}

# =========================
# Build maps up-front (so later steps never crash)
# =========================
build_barcode_to_phage
seed_phages_from_map
seed_phages_from_existing_outputs

# =========================
# STEP 2: Demultiplex WITH trimming → per-barcode FASTQ
# =========================
if should_run demux; then
  if compgen -G "$DEMUX_DIR/*barcode*.fastq*" > /dev/null; then
    echo "[SKIP] Demux-trim outputs present in $DEMUX_DIR"
  else
    echo "[RUN ] Demux (TRIM ON) → $DEMUX_DIR"
    DEMUX_CMD=( "$DORADO" demux --kit-name "$KIT" --emit-fastq "$ALL_FASTQ" --output-dir "$DEMUX_DIR" --emit-summary )
    if [[ "$BARCODE_BOTH_ENDS" == "1" ]]; then
      DEMUX_CMD+=( --barcode-both-ends )
    fi
    "${DEMUX_CMD[@]}"
    echo "[DONE] Demux (trimmed)."
  fi
else
  echo "[SKIP] Step demux (per switch)"
fi

# =========================
# STEP 3: CONCAT per-phage (trimmed per-barcode → M*.trimmed.fastq.gz)
# =========================
if should_run concat; then
  echo "[RUN ] Concatenate per-phage → $PER_PHAGE_RAW"
  while IFS=$'\t' read -r phage barcs; do
    [[ -z "${phage:-}" || -z "${barcs:-}" ]] && continue
    [[ "${phage}" == "phage" || "${phage:0:1}" == "#" ]] && continue

    mapfile -t files < <(
      IFS=',' read -ra arr <<< "$barcs"
      for bc in "${arr[@]}"; do
        bc="${bc//[[:space:]]/}"
        find "$DEMUX_DIR" -maxdepth 1 -type f \
          \( -name "*_${bc}.fastq.gz" -o -name "*_${bc}.fastq" -o -name "${bc}.fastq.gz" -o -name "${bc}.fastq" \)
      done | sort
    )
    if (( ${#files[@]} == 0 )); then
      echo "[WARN] No per-barcode FASTQs for $phage ($barcs). Skipping."
      continue
    fi

    out="$PER_PHAGE_RAW/${phage}.trimmed.fastq.gz"
    if [[ -s "$out" ]]; then
      echo "[SKIP] $phage exists: $out"
      continue
    fi
    echo "[RUN ] $phage: concat ${#files[@]} → $out"
    { for f in "${files[@]}"; do
        case "$f" in
          *.fastq.gz) $DECOMP "$f" ;;
          *.fastq)    cat "$f" ;;
          *)          echo "[WARN] Skip unknown extension: $f" >&2 ;;
        esac
      done
    } | $COMPRESS > "${out}.partial"
    mv "${out}.partial" "$out"

    if ! is_fastq_gz "$out"; then
      echo "[WARN] FASTQ structure check failed for $phage: $out (continuing)"
    fi
  done < "$PHAGE_MAP"
  echo "[DONE] Per-phage concatenation (trimmed)."
else
  echo "[SKIP] Step concat (per switch)"
fi

# =========================
# STEP 3.5: Porechop-ABI per-barcode → concat per-phage + summaries
# =========================
if should_run porechop; then
  echo "[RUN ] porechop_abi per-barcode → $PORECHOP_DIR (logs per barcode)"
  source "$ONT_ENV_ACTIVATE"; conda activate "$ONT_ENV_NAME"

  mkdir -p "$PORECHOP_DIR"

  # Find demuxed per-barcode FASTQs (already dorado-trimmed)
  shopt -s nullglob
  mapfile -t BC_FASTQS < <(
    find "$DEMUX_DIR" -maxdepth 1 -type f \( \
        -name '*_barcode[0-9][0-9].fastq'      -o \
        -name '*_barcode[0-9][0-9].fastq.gz'   -o \
        -name '*barcode[0-9][0-9].fastq'       -o \
        -name '*barcode[0-9][0-9].fastq.gz' \
      \) | sort
  )
  if (( ${#BC_FASTQS[@]} == 0 )); then
    echo "[WARN] No per-barcode FASTQs found in $DEMUX_DIR. Did you run demux? File pattern expects '*barcodeNN.fastq[.gz]'"
  fi

  # per-run barcode summary (one row per processed barcode FASTQ)
  RUN_BC_SUM="$PORECHOP_DIR/_per_barcode_summary.tsv"
  echo -e "phage\tbarcode\ttotal_reads\tstart_trimmed\tstart_bp_removed\tend_trimmed\tend_bp_removed\tsplit_reads\tbarcode_counts" > "$RUN_BC_SUM"

  declare -A local_seen_ph=()

  for INBC in "${BC_FASTQS[@]}"; do
    # derive 'barcodeNN' token from filename
    if [[ "$(basename "$INBC")" =~ (barcode[0-9]{2}) ]]; then
      bc="${BASH_REMATCH[1],,}"  # lower-case
    else
      echo "[WARN] Skip (cannot extract barcodeNN): $(basename "$INBC")"
      continue
    fi

    ph="${B2P[$bc]-}"   # empty if missing
    if [[ -z "$ph" ]]; then
      echo "[WARN] $bc has no mapping in $PHAGE_MAP; skipping."
      continue
    fi
    PH_SEEN["$ph"]=1
    local_seen_ph["$ph"]=1

    ph_bc_dir="$PORECHOP_DIR/$ph/barcodes"
    mkdir -p "$ph_bc_dir"

    OUT_FASTQ="$ph_bc_dir/${ph}.${bc}.porechop.fastq.gz"
    LOG="$ph_bc_dir/${ph}.${bc}.porechop.log"

    if [[ -s "$OUT_FASTQ" && -s "$LOG" ]]; then
      echo "[SKIP] porechop_abi exists: $(basename "$OUT_FASTQ")"
    else
      echo "[RUN ] porechop_abi ${ph}/${bc}"
      porechop_abi \
        -i "$INBC" \
        -o "$OUT_FASTQ" \
        --format fastq.gz \
        -t "$THREADS" \
        -v "$PORECHOP_VERBOSITY" \
        --min_split_read_size 2000 \
        >"$LOG" 2>&1 \
        || { echo "[WARN] porechop_abi failed for ${ph}/${bc} (see $LOG)"; continue; }

      if ! is_fastq_gz "$OUT_FASTQ"; then
        echo "[WARN] porechop_abi output failed FASTQ check: $OUT_FASTQ"
      fi
    fi

    # per-barcode summary (with phage column)
    SUM_BC="${OUT_FASTQ%.fastq.gz}.summary.tsv"
    if [[ ! -s "$SUM_BC" ]]; then
      echo -e "phage\tbarcode\ttotal_reads\tstart_trimmed\tstart_bp_removed\tend_trimmed\tend_bp_removed\tsplit_reads\tbarcode_counts" > "$SUM_BC"
      _parse_porechop_log "$LOG" "$ph" "$bc" >> "$SUM_BC"
    fi

    # append/update to run-level barcode summary (avoid duplicates)
    tail -n1 "$SUM_BC" >> "$RUN_BC_SUM"
  done

  # Concat the *porechop* per-barcode outputs per phage; build per-phage roll-up
  for ph in "${!local_seen_ph[@]}"; do
    ph_dir="$PORECHOP_DIR/$ph"
    ph_cat="$ph_dir/${ph}.porechop.fastq.gz"
    mkdir -p "$ph_dir"

    mapfile -t PH_BARCODE_FASTQS < <(find "$ph_dir/barcodes" -maxdepth 1 -type f -name "${ph}.barcode*.porechop.fastq.gz" | sort)
    if (( ${#PH_BARCODE_FASTQS[@]} == 0 )); then
      echo "[WARN] No porechopped barcodes to concat for $ph"
      continue
    fi

    if [[ -s "$ph_cat" ]]; then
      echo "[SKIP] phage concat exists: $ph_cat"
    else
      echo "[RUN ] Concat porechop outputs → ${ph_cat}"
      { for f in "${PH_BARCODE_FASTQS[@]}"; do pigz -dc "$f"; done; } | $COMPRESS > "${ph_cat}.partial"
      mv "${ph_cat}.partial" "$ph_cat"
      if ! is_fastq_gz "$ph_cat"; then
        echo "[WARN] Concat FASTQ check failed for $ph: $ph_cat"
      fi
    fi

    # phage roll-up summary by summing per-barcode summaries
    PH_BC_SUMS=( "$ph_dir"/barcodes/*.summary.tsv )
    PH_SUM="$ph_dir/porechop_abi.summary.tsv"
    awk -F'\t' -v ph="$ph" '
      BEGIN{ OFS="\t";
             print "phage\ttotal_reads\tstart_trimmed\tstart_bp_removed\tend_trimmed\tend_bp_removed\tsplit_reads\tbarcode_counts" }
      FNR==1 && NR>1 { next }
      NR>1 {
        # cols: 1=phage 2=barcode 3=total_reads 4=start_n 5=start_bp 6=end_n 7=end_bp 8=split 9=bc_counts
        if ($1==ph) {
          tr+=($3=="NA"||$3==""?0:$3)
          st+=($4==""?0:$4)
          sb+=($5==""?0:$5)
          et+=($6==""?0:$6)
          eb+=($7==""?0:$7)
          sp+=($8==""?0:$8)
          if ($9!="NA" && $9!="") { bc=(bc!=""?bc";":"")$9 }
        }
      }
      END{
        tot = (tr>0?tr:"NA")
        print ph, tot, st+0, sb+0, et+0, eb+0, sp+0, (bc!=""?bc:"NA")
      }' "${PH_BC_SUMS[@]}" > "$PH_SUM"
  done
  echo "[DONE] porechop_abi stage."
else
  echo "[SKIP] Step porechop (per switch)"
fi

# =========================
# SELECT DOWNSTREAM INPUTS (for filtering + 'trimmed' QC)
# =========================
shopt -s nullglob
declare -a SRC_FOR_FILTER_LIST TRIMMED_FOR_QC_LIST
if [[ "$DOWNSTREAM_SET" == "porechop" ]]; then
  SRC_FOR_FILTER_LIST=( "$PORECHOP_DIR"/M*/M*.porechop.fastq.gz )
  TRIMMED_FOR_QC_LIST=( "${SRC_FOR_FILTER_LIST[@]}" )
  echo "[INFO] DOWNSTREAM_SET=porechop → use per-phage porechop outputs."
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
if should_run filter; then
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
else
  echo "[SKIP] Step filter (per switch)"
fi

# =========================
# STEP 5: QC — NanoPlot (on selected 'trimmed' set)
# =========================
if should_run nanoplot; then
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
  else
    echo "[SKIP] NanoPlot (RUN_NANOPLOT=0)"
  fi
else
  echo "[SKIP] Step nanoplot (per switch)"
fi

# =========================
# STEP 6: FastQC/MultiQC on selected 'trimmed' set
# =========================
if should_run qc_trim; then
  if [[ "$RUN_FASTQC_TRIM" == "1" ]]; then
    mkdir -p "$FQCTRIM_DIR"
    if (( ${#TRIMMED_FOR_QC_LIST[@]} > 0 )); then
      "$FASTQC_BIN" -t "$THREADS" -o "$FQCTRIM_DIR" "${TRIMMED_FOR_QC_LIST[@]}"
      if [[ -f "$MULTIQC_ACTIVATE" ]]; then source "$MULTIQC_ACTIVATE"; fi
      multiqc -o "$FQCTRIM_DIR" "$FQCTRIM_DIR"
    else
      echo "[WARN] No inputs found for FastQC (trimmed set: $DOWNSTREAM_SET)."
    fi
  else
    echo "[SKIP] FastQC/MultiQC trimmed (RUN_FASTQC_TRIM=0)"
  fi
else
  echo "[SKIP] Step qc_trim (per switch)"
fi

# =========================
# STEP 7: FastQC/MultiQC on filtered
# =========================
if should_run qc_filt; then
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
  else
    echo "[SKIP] FastQC/MultiQC filtered (RUN_FASTQC_FILT=0)"
  fi
else
  echo "[SKIP] Step qc_filt (per switch)"
fi

# =========================
# STEP 8: Summaries (run-level rollups)
# =========================
if should_run summarize; then
  echo "[RUN ] Building run-level summaries"
  # build a per-phage run-level table (one line per phage) from existing per-phage summaries
  RUN_PH_SUM="$PORECHOP_DIR/_per_phage_summary.tsv"
  echo -e "phage\ttotal_reads\tstart_trimmed\tstart_bp_removed\tend_trimmed\tend_bp_removed\tsplit_reads\tbarcode_counts" > "$RUN_PH_SUM"
  while IFS= read -r ph; do
    if [[ -s "$PORECHOP_DIR/$ph/porechop_abi.summary.tsv" ]]; then
      tail -n1 "$PORECHOP_DIR/$ph/porechop_abi.summary.tsv" >> "$RUN_PH_SUM"
    fi
  done < <(printf "%s\n" "${!PH_SEEN[@]}" | sort)
  echo "[DONE] Summaries."
else
  echo "[SKIP] Step summarize (per switch)"
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
