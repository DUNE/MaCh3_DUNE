#!/bin/bash
#SBATCH --job-name=ProcessMCMC
#SBATCH --account=def-nilic
#SBATCH --output=ProcessMCMC_%j.out
#SBATCH --error=ProcessMCMC_%j.err
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G

# ---------------- Modules ---------------- #
module purge
module load StdEnv/2020
module load cmake/3.18.4
module load gcc/9.3.0
module load python/3.9.6
module load root/6.26.06

# ---------------- Setup environment ---------------- #
source /scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/build/bin/setup.MaCh3.sh
source /scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/build/bin/setup.MaCh3DUNE.sh

# ---------------- Paths ---------------- #
BIN_DIR=/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/build/bin

CONFIG_YAML=/scratch/abipeake/Off_axis_chains/newadaptiveonaxis/cfg/newadaptiveonaxis_chain_0_job_0.yaml
OUTPUT_ROOT=/scratch/abipeake/Off_axis_chains/newadaptiveonaxis/newadaptiveonaxis_chain_0_job_0.root

# ---------------- Run ProcessMCMC ---------------- #
echo "Running ProcessMCMC..."
srun "$BIN_DIR/ProcessMCMC" "$CONFIG_YAML" "$OUTPUT_ROOT"

echo "ProcessMCMC finished."
