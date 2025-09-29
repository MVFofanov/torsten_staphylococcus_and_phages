#!/usr/bin/env bash
#SBATCH --job-name=genomad_bakta_padloc
#SBATCH --partition=short
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=80
#SBATCH --mem=80G
#SBATCH --output=/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages/slurm_logs/%x/%x_%j.log

set -euo pipefail

# ===== paths =====
wd="/work/zo49sog/crassvirales/torsten_staphylococcus_and_phages"
threads=80
consensus_genome="${wd}/data/StaphIII12_consensus.fna"
file_name=$(basename "${consensus_genome%.fna}")

genomad_db="/veodata/03/databases/geNomad/v1.3"
genomad_outdir="${wd}/genomad_${file_name}"

bakta_db="/veodata/03/databases/bakta/v20240516/db"
bakta_outdir="${wd}/bakta_${file_name}"

infernal_outdir="${wd}/infernal_${file_name}"
infernal_db="/home/zo49sog/miniconda3/envs/padloc_env/data/cm/padlocdb.cm"

crisprdetect_outdir="${wd}/crisprdetect_${file_name}"
padloc_outdir="${wd}/padloc_${file_name}"

run_infernal="/home/zo49sog/miniconda3/envs/padloc_env/bin/bin/run-infernal"
run_crisprdetect="/home/zo49sog/miniconda3/envs/padloc_env/bin/bin/run-crisprdetect"

# # ===== genomad =====
# source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh
# conda activate genNomad_v20230721

# genomad end-to-end --cleanup "${consensus_genome}" "${genomad_outdir}" "${genomad_db}"
# conda deactivate

# # ===== bakta =====
# source /vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh
# conda activate bakta_v1.9.3

# bakta --db "${bakta_db}" --threads "${threads}" --output "${bakta_outdir}" \
#       --prefix "${file_name}" ${consensus_genome}
# conda deactivate

# ===== padloc + helpers =====
source /home/zo49sog/miniconda3/etc/profile.d/conda.sh
# conda activate padloc_env
conda activate padloc

# Step 1: run-infernal
mkdir -p "${infernal_outdir}" "${crisprdetect_outdir}" "${padloc_outdir}"

# "${run_infernal}" --input "${consensus_genome}" \
#              --output "${infernal_outdir}/${file_name}_ncrna.tblout"

# cmsearch \
#   -Z 10 \
#   --FZ 500 \
#   --acc \
#   --noali \
#   --tblout "${infernal_outdir}/${file_name}_ncrna.tblout" \
#   "${infernal_db}" \
#   "${consensus_genome}"

sed 's/\t/ /18g' "${infernal_outdir}/${file_name}_ncrna.tblout" > "${infernal_outdir}/${file_name}_ncrna.tblout.formatted"

# # Step 2: run-crisprdetect
# "${run_crisprdetect}" --input "${consensus_genome}" \
#                  --output "${crisprdetect_outdir}/${file_name}_crispr"

# Step 3: padloc
# --fna "${consensus_genome}" 
# --faa "${bakta_outdir}/${file_name}.faa" --gff "${bakta_outdir}/${file_name}.gff3"
# --ncrna "${infernal_outdir}/${file_name}_ncrna.tblout.formatted"
# --ncrna "${infernal_outdir}/${file_name}_ncrna.tblout"
padloc --fna "${consensus_genome}" \
       --crispr "${crisprdetect_outdir}/${file_name}_crispr.gff" \
       --cpu "${threads}" \
       --force \
       --outdir "${padloc_outdir}"

conda deactivate

echo "[INFO] All analyses finished successfully."
