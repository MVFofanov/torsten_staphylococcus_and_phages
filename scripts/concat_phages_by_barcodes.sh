#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./concat_by_phage.sh <READS_DIR> <MAP_TSV> <OUT_DIR>
#
# Example:
#   ./concat_by_phage.sh \
#     /work/.../ont_basecalls_20251007_160801/demux_no_trim \
#     phage_barcode_map.tsv \
#     /work/.../per_phage_raw

READS_DIR="${1:?Give READS_DIR (folder with per-barcode FASTQs)}"
MAP_TSV="${2:?Give mapping TSV (columns: phage<TAB>barcode_list)}"
OUT_DIR="${3:?Give OUT_DIR for per-phage outputs}"

mkdir -p "$OUT_DIR"

# pick compressor/decompressor
if command -v pigz >/dev/null 2>&1; then
  COMPRESS="pigz -c -p ${SLURM_CPUS_PER_TASK:-8}"
  DECOMP_GZ="pigz -dc"
else
  COMPRESS="gzip -c"
  DECOMP_GZ="gzip -dc"
fi

echo "[INFO] Reads: $READS_DIR"
echo "[INFO] Map  : $MAP_TSV"
echo "[INFO] Out  : $OUT_DIR"

# read mapping; skip header if present
# expected format:
# phage<TAB>barcode01,barcode02
# phage names like M1, M2, ...
# barcodes like barcode01,barcode02 (zero-padded)
tail -n +1 "$MAP_TSV" | while IFS=$'\t' read -r phage barcs; do
  # skip empty/comment lines
  [[ -z "${phage:-}" || -z "${barcs:-}" ]] && continue
  [[ "${phage:0:1}" == "#" || "${phage}" == "phage" ]] && continue

  # build list of files for this phage
  mapfile -t files < <(
    IFS=',' read -ra arr <<< "$barcs"
    for bc in "${arr[@]}"; do
      bc="${bc//[[:space:]]/}"                      # strip spaces
      # match both plain and gz fastq
      # Dorado demux typically names: ..._barcodeXX.fastq
      # adjust the glob if your names differ
      find "$READS_DIR" -maxdepth 1 -type f -name "*_${bc}.fastq.gz" -o -name "*_${bc}.fastq" 2>/dev/null
    done | sort
  )

  if (( ${#files[@]} == 0 )); then
    echo "[WARN] No FASTQs found for $phage ($barcs). Skipping."
    continue
  fi

  out="$OUT_DIR/${phage}.raw.fastq.gz"
  if [[ -s "$out" ]]; then
    echo "[SKIP] $phage already exists: $out"
    continue
  fi

  echo "[RUN ] Concatenating ${#files[@]} files for $phage -> $out"
  # stream all inputs (mix of .fastq and .fastq.gz) into one gzipped output
  {
    for f in "${files[@]}"; do
      case "$f" in
        *.fastq.gz)  $DECOMP_GZ "$f" ;;
        *.fastq)     cat "$f" ;;
        *)           echo "[WARN] Skipping unknown extension: $f" >&2 ;;
      esac
    done
  } | $COMPRESS > "$out"

  echo "[DONE] $phage: $(du -h "$out" | cut -f1) -> $out"
done

echo "[INFO] Finished. Outputs in: $OUT_DIR"
