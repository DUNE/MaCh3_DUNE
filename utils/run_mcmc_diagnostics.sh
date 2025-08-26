#!/bin/bash
#SBATCH --job-name=MCMCWorkflow
#SBATCH --account=def-nilic
#SBATCH --output=MCMCWorkflow_%j.out
#SBATCH --error=MCMCWorkflow_%j.err
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1        # adjust if your code is multithreaded
#SBATCH --mem=80G

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
# Path to your binaries
BIN_DIR=/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/build/bin

# Input ROOT file and YAML config
INPUT_ROOT=/scratch/abipeake/newChains_August/EventRates_Elep_erec_afteradaptive_total/EventRates_Elep_erec_afteradaptive_total_chain_0_job_0.root         #scratch/abipeake/August_tesing/tuesday_newconfig_adaptive_enubias_ptpzenurec/tuesday_newconfig_adaptive_enubias_ptpzenurec_chain_0_job_0.root
CONFIG_YAML=/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/TutorialConfigs/TutorialDiagConfig.yaml

# Directory where the input ROOT file lives (outputs go here)
WORKDIR=$(dirname "$INPUT_ROOT")
cd "$WORKDIR" || { echo "Failed to cd to $WORKDIR"; exit 1; }

# DiagMCMC output ROOT file
OUTPUT_DIAG=${INPUT_ROOT%.root}_MCMC_Diag.root

# ---------------- Step 1: DiagMCMC ---------------- #
echo "Running DiagMCMC..."
srun "$BIN_DIR/DiagMCMC" "$INPUT_ROOT" "$CONFIG_YAML"

# ---------------- Step 2: PlotMCMCDiag ---------------- #
echo "Running PlotMCMCDiag on $OUTPUT_DIAG..."
srun "$BIN_DIR/PlotMCMCDiag" "$OUTPUT_DIAG"

# ---------------- Step 3: ProcessMCMC ---------------- #
echo "Running ProcessMCMC..."
srun "$BIN_DIR/ProcessMCMC" "$CONFIG_YAML" "$INPUT_ROOT"

echo "Workflow completed. All outputs are in $WORKDIR"
