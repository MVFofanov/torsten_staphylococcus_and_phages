#!/usr/bin/env bash
#SBATCH --job-name=compare_staph_assemblies
#SBATCH --partition=short
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=80
#SBATCH --mem=80G
#SBATCH --output=/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages/slurm_logs/%x/%x_%j.log
#SBATCH --open-mode=append

set -euo pipefail

# ======== paths ========
WD="/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages"

# Your inputs
REF="${WD}/data/ncbi_dataset/data/GCF_016026635.1/GCF_016026635.1_ASM1602663v1_genomic.fna"
MEGA="${WD}/megahit/contigs.fa"
SPAD="${WD}/spades/scaffolds.fa"
CONS="${WD}/data/StaphIII12_consensus.fna"

# Reads for coverage checks
R1="${WD}/trimmed_merged/StaphIII12_R1.trim.fq.gz"
R2="${WD}/trimmed_merged/StaphIII12_R2.trim.fq.gz"

# Tools
CONDA_SETUP="/vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh"
# mummer needs gnuplot for PNG plots; we keep paths explicit below
SAMTOOLS="/home/groups/VEO/tools/samtools/v1.17/bin/samtools"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
RUN_TAG="$(date +%Y%m%d_%H%M%S)"
OUT="${WD}/assembly_comparison/${RUN_TAG}"
mkdir -p "${OUT}"/{quast,mummer,fastani,bbmap}
echo "[INFO] Output: ${OUT}"

# ======== env activation (for python if needed) ========
# shellcheck disable=SC1090
source "$CONDA_SETUP"

# ======== collect assemblies that exist ========
declare -a ASMS=()
declare -a LABELS=()

add_if_exist () {
  local f="$1"; local label="$2"
  if [[ -s "$f" ]]; then
    ASMS+=("$f")
    LABELS+=("$label")
    echo "[INFO] Added $label: $f"
  else
    echo "[WARN] Missing $label: $f (skipping)"
  fi
}

[[ -s "$REF" ]] || { echo "[ERROR] Reference not found: $REF"; exit 1; }
add_if_exist "$MEGA" "megahit"
add_if_exist "$SPAD" "spades"
add_if_exist "$CONS" "consensus"

if ((${#ASMS[@]}==0)); then
  echo "[ERROR] No assemblies to compare"; exit 1
fi

# ======== 1) QUAST multi-assembly vs reference ========
echo "[STEP] QUAST"
# Headless backend for matplotlib
export MPLBACKEND=Agg
labels_csv="$(IFS=,; echo "${LABELS[*]}")"

python3 /home/groups/VEO/tools/quast/v5.2.0/quast.py \
  -r "$REF" \
  --min-contig 500 \
  --threads "$THREADS" \
  --labels "$labels_csv" \
  -o "${OUT}/quast" \
  "${ASMS[@]}"

echo "[INFO] QUAST report:"
echo "  ${OUT}/quast/report.html"
echo "  ${OUT}/quast/transposed_report.tsv"

# ======== 2) Whole-genome alignments & dotplots (MUMmer4) ========
echo "[STEP] MUMmer (nucmer + mummerplot)"
if ! command -v gnuplot >/dev/null 2>&1; then
  echo "[WARN] gnuplot not found; mummerplot PNGs may not be generated (coords.txt will still be produced)."
fi
# mummerplot has a perl dep; activate if needed
source /vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh && conda activate perl_v5.32.1

for i in "${!ASMS[@]}"; do
  asm="${ASMS[$i]}"; label="${LABELS[$i]}"
  prefix="${OUT}/mummer/${label}_vs_ref"
  echo "[INFO] Aligning $label -> ref"
  /home/groups/VEO/tools/mummer/v4.1/bin/nucmer --maxmatch -l 100 -c 500 -p "$prefix" "$REF" "$asm"
  /home/groups/VEO/tools/mummer/v4.1/bin/delta-filter -1 "$prefix.delta" > "$prefix.1delta"
  /home/groups/VEO/tools/mummer/v4.1/bin/show-coords -rcl "$prefix.1delta" > "$prefix.coords.txt"
  /home/groups/VEO/tools/mummer/v4.1/bin/mummerplot --png --large --fat -p "$prefix" "$prefix.1delta" 2>/dev/null || true

  # Also a many-to-many view to catch paralogs/repeats:
  /home/groups/VEO/tools/mummer/v4.1/bin/delta-filter -r -q "$prefix.delta" > "$prefix.rq.delta"
  /home/groups/VEO/tools/mummer/v4.1/bin/mummerplot --png --large --fat -p "${prefix}.rq" "$prefix.rq.delta" 2>/dev/null || true
done

# ======== 3) ANI vs reference (fastANI) ========
echo "[STEP] fastANI"
for i in "${!ASMS[@]}"; do
  asm="${ASMS[$i]}"; label="${LABELS[$i]}"
  /home/groups/VEO/tools/fastANI/v1.33/fastANI \
    --query "$asm" \
    --ref "$REF" \
    --threads "$THREADS" \
    --fragLen 3000 \
    -o "${OUT}/fastani/${label}_vs_ref.fastani.tsv" || true
done

# ======== 4) Read mapping to each assembly (coverage, potential misassemblies) ========
echo "[STEP] Read mapping (BBMap -> BAM via samtools sort)"
for i in "${!ASMS[@]}"; do
  asm="${ASMS[$i]}"; label="${LABELS[$i]}"
  echo "[INFO] Mapping reads to $label"
  workdir="${OUT}/bbmap/${label}"
  mkdir -p "$workdir"
  pushd "$workdir" >/dev/null

  # Stream SAM to samtools sort -> BAM; keep covstats & mapping stats
  /home/groups/VEO/tools/bbmap/v39.06/bbmap.sh \
      ref="$asm" in="$R1" in2="$R2" t="$THREADS" nodisk \
      covstats="${label}.covstats.txt" statsfile="${label}.stats.txt" \
      out=stdout.sam 2> "${label}.bbmap.log" \
    | "$SAMTOOLS" sort -@ "$THREADS" -o "${label}.sorted.bam" -

  "$SAMTOOLS" index "${label}.sorted.bam"

  # quick mean depth per contig
  awk 'NR>1{sum+=$2*$3; len+=$2} END{if(len>0) printf("Approx_mean_depth: %.2f\n", sum/len)}' \
      "${label}.covstats.txt" > "${label}.depth_summary.txt" || true

  popd >/dev/null
done

# ======== 5) Summary table from QUAST + ANI + depth ========
echo -e "label\tN50\tNGA50\tGenome_fraction(%)\tMisassemblies\tANI(%)\tApprox_mean_depth" > "${OUT}/summary.tsv"
QT="${OUT}/quast/transposed_report.tsv"

get_quast_metric () {
  # Args: metric label
  local metric="$1"; local label="$2"
  awk -v l="$label" -v m="$metric" -F'\t' '
    NR==1 { for(i=2;i<=NF;i++) if($i==l) {idx=i} next }
    $1==m { if(idx>0) print $idx; else print "NA" }
  ' "$QT"
}

for label in "${LABELS[@]}"; do
  n50=$(get_quast_metric "N50" "$label")
  nga50=$(get_quast_metric "NGA50" "$label")
  gf=$(get_quast_metric "Genome fraction (%)" "$label")
  mis=$(get_quast_metric "# misassemblies" "$label")
  ani="NA"
  if [[ -s "${OUT}/fastani/${label}_vs_ref.fastani.tsv" ]]; then
    ani=$(awk '{print $3; exit}' "${OUT}/fastani/${label}_vs_ref.fastani.tsv" 2>/dev/null || echo "NA")
  fi
  depth="NA"
  if [[ -s "${OUT}/bbmap/${label}/${label}.depth_summary.txt" ]]; then
    depth=$(awk '{print $2}' "${OUT}/bbmap/${label}/${label}.depth_summary.txt" 2>/dev/null || echo "NA")
  fi
  echo -e "${label}\t${n50}\t${nga50}\t${gf}\t${mis}\t${ani}\t${depth}" >> "${OUT}/summary.tsv"
done

echo "[DONE] Comparison finished."
echo "Open QUAST HTML: ${OUT}/quast/report.html"
echo "See dotplots:    ${OUT}/mummer/*png"
echo "ANI tables:      ${OUT}/fastani/*fastani.tsv"
echo "Coverage stats:  ${OUT}/bbmap/*/*.covstats.txt"
echo "Summary TSV:     ${OUT}/summary.tsv"
