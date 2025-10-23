# paths
WD="/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages"
POD5_DIR="/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages/data/250715_StaphPhages_TS/250715_StaphPhages_TS/20250715_0915_MD-102423_FBA89536_6b18c2cc/pod5_skip"
OUT="${WD}/phages_analysis/ont_basecalls_$(date +%Y%m%d_%H%M%S)"
MODEL="dna_r10.4.1_e8.2_400bps_sup@v4.3.0"      # R10.4.1 SUP
KIT="SQK-NBD114-96"                              # <- put YOUR barcode kit here

DORADO=" /home/groups/VEO/tools/dorado/v1.1.0/bin/dorado"
mkdir -p "$OUT"

# pick device
if command -v nvidia-smi >/dev/null; then DEV="cuda:all"; else DEV="cpu"; fi

# 1) Basecall ALL pod5s in the folder (recursive)
"$DORADO" download --model "$MODEL" || true
"$DORADO" basecaller "$MODEL" "$POD5_DIR" --recursive --device "$DEV" --emit-fastq \
  | gzip > "$OUT/all.fastq.gz"

# 2) Demultiplex with the exact kit name
mkdir -p "$OUT/demux"
"$DORADO" demux --kit-name "$KIT" --emit-fastq "$OUT/all.fastq.gz" --output-dir "$OUT/demux"

# results:
#   $OUT/demux/barcode01.fastq
#   $OUT/demux/barcode02.fastq
#   $OUT/demux/unclassified.fastq
