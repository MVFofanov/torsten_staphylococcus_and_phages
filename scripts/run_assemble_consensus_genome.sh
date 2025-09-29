#!/usr/bin/env bash
#SBATCH --job-name=assemble_consensus_genome
#SBATCH --partition=short
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=80
#SBATCH --mem=80G
#SBATCH --output=/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages/slurm_logs/%x/%x_%j.log

set -euo pipefail

source /home/zo49sog/miniconda3/etc/profile.d/conda.sh && conda activate samtools_env

# ==================== paths ====================
WD=/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages
REF=${WD}/data/ncbi_dataset/data/GCF_016026635.1/GCF_016026635.1_ASM1602663v1_genomic.fna
R1=${WD}/trimmed_merged/StaphIII12_R1.trim.fq.gz
R2=${WD}/trimmed_merged/StaphIII12_R2.trim.fq.gz
OUT=${WD}/ref_consensus_StaphIII12

THREADS=${SLURM_CPUS_PER_TASK:-80}

mkdir -p "$OUT"

# ==================== pipeline ====================
echo "[INFO] Indexing reference"
bwa-mem2 index "$REF" || bwa index "$REF"
samtools faidx "$REF"

# 1) align -> name sort (for fixmate)
echo "[INFO] Mapping reads"
bwa-mem2 mem -t "$THREADS" "$REF" "$R1" "$R2" \
  | samtools sort -@ "$THREADS" -n -o "$OUT/aln.namesorted.bam" -

# 2) fix mate info (required for accurate markdup on PE)
samtools fixmate -@ "$THREADS" -m "$OUT/aln.namesorted.bam" "$OUT/aln.fixmate.bam"

# 3) coordinate sort
samtools sort -@ "$THREADS" -o "$OUT/aln.coordsorted.bam" "$OUT/aln.fixmate.bam"
samtools index "$OUT/aln.coordsorted.bam"

# 4) mark duplicates
samtools markdup -@ "$THREADS" -s "$OUT/aln.coordsorted.bam" "$OUT/aln.mkdup.bam"
samtools index "$OUT/aln.mkdup.bam"

# quick QC
samtools flagstat "$OUT/aln.mkdup.bam" > "$OUT/flagstat.txt"
samtools coverage "$OUT/aln.mkdup.bam" > "$OUT/coverage.txt"
samtools idxstats "$OUT/aln.mkdup.bam" > "$OUT/idxstats.txt"

# ---------- low-coverage mask (depth < 5) ----------
echo "[INFO] Generating low coverage mask (depth<5)"
samtools depth -aa "$OUT/aln.mkdup.bam" \
  | awk '$3<5{print $1"\t"$2-1"\t"$2}' \
  | bedtools merge -i - > "$OUT/lowcov_lt5.bed"

# ---------- variant calling (haploid, with MQ/BQ filters) ----------
echo "[INFO] Calling variants"
bcftools mpileup -Ou -f "$REF" -q 20 -Q 20 "$OUT/aln.mkdup.bam" \
| bcftools call -mv --ploidy 1 -Ou \
| bcftools filter -s LowQual -e 'QUAL<20 || DP<5' \
  -Oz -o "$OUT/variants.vcf.gz"

bcftools index "$OUT/variants.vcf.gz"

# ---------- consensus ----------
echo "[INFO] Building consensus"
bcftools consensus -f "$REF" -m "$OUT/lowcov_lt5.bed" "$OUT/variants.vcf.gz" \
  > "$OUT/StaphIII12_consensus.fna"

echo "[INFO] Done. Consensus genome at:"
echo "  $OUT/StaphIII12_consensus.fna"
