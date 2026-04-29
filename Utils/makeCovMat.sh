#!/bin/bash
#SBATCH --job-name=MakeCovMatrix
#SBATCH --account=def-deborahh
#SBATCH --output=MakeCovMatrix_%j.out
#SBATCH --error=MakeCovMatrix_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=400G

# ---------------- Modules ---------------- #
module purge
module load StdEnv/2020
module load cmake/3.18.4
module load gcc/9.3.0
module load python/3.9.6
module load root/6.26.06

# ---------------- Setup environment ---------------- #
source /scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/build/bin/setup.MaCh3.sh
source /scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/build/bin/setup.MaCh3DUNE.sh

# ---------------- Paths ---------------- #
BASE_DIR=/scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE
BIN=$BASE_DIR/build/Apps/MakeCovarianceMatrix

EVENT_CFG=$BASE_DIR/Configs/BeamOffAxis/EventRate.yml
SYST_CFG=$BASE_DIR/Configs/BeamOffAxis/NDDetSysts.yml

OUTPUT=normalbinning_covmat.root

# ---------------- Run ---------------- #
echo "[INFO] Starting covariance matrix generation"
echo "[INFO] Using configs:"
echo "  Event: $EVENT_CFG"
echo "  Systs: $SYST_CFG"

$BIN $EVENT_CFG $SYST_CFG 1000000 $OUTPUT

echo "[INFO] Finished"