#!/bin/bash
#SBATCH --job-name=MCMCWorkflow
#SBATCH --account=def-nilic
#SBATCH --output=MCMCWorkflow_%j.out
#SBATCH --error=MCMCWorkflow_%j.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G

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
INPUT_CHAIN=/scratch/abipeake/August_tesing/monday_noadaption_q0q3_ptpzenu2/monday_noadaption_q0q3_ptpzenu2_chain_0_job_0.root
INPUT_ROOT=/scratch/abipeake/August_tesing/monday_noadaption_q0q3_ptpzenu2/monday_noadaption_q0q3_ptpzenu2_chain_0_job_0_MCMC_Diag.root
CONFIG_YAML=/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/TutorialConfigs/TutoriaDiagConfigq0q3.yaml

WORKDIR=$(dirname "$INPUT_ROOT")
cd "$WORKDIR" || { echo "Failed to cd to $WORKDIR"; exit 1; }

# ---------------- Step: PlotMCMCDiag ---------------- #
#echo "Running PlotMCMCDiag on $INPUT_ROOT..."
#srun "$BIN_DIR/PlotMCMCDiag" "$INPUT_ROOT"

# ---------------- Optional: ProcessMCMC ---------------- #
# If you want to process the chain into posteriors, uncomment:
srun "$BIN_DIR/ProcessMCMC" "$CONFIG_YAML" "$INPUT_CHAIN"

echo "Workflow completed. All outputs are in $WORKDIR"
