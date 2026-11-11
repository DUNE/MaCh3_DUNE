#!/bin/bash
#SBATCH --job-name=MakeCovMatrix
#SBATCH --account=rpp-dune
#SBATCH --output=/scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/Fractional_covmat_offaxis/MakeCovMatrix_%A_%a.out
#SBATCH --error=/scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/Fractional_covmat_offaxis//MakeCovMatrix_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --array=0-50

module purge
module load StdEnv/2020
module load cmake/3.18.4
module load gcc/9.3.0
module load python/3.9.6
module load root/6.26.06

export MACH3_SEED=$((12345 + SLURM_ARRAY_TASK_ID))

source /scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/build/bin/setup.MaCh3.sh
source /scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/build/bin/setup.MaCh3DUNE.sh

BASE_DIR=/scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE
BIN=$BASE_DIR/build/Apps/MakeCovarianceMatrix

EVENT_CFG=$BASE_DIR/Configs/BeamOffAxis/Fit_new.yaml
SYST_CFG=$BASE_DIR/Configs/BeamOffAxis/NDDetSysts.yml

NTHROWS=50000

OUTPUT=/scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/Fractional_covmat_offaxis/covmat_010825_fractionalcovmat_offaxis_${SLURM_ARRAY_TASK_ID}.root

echo "[INFO] Running array job ${SLURM_ARRAY_TASK_ID}"

$BIN $EVENT_CFG $SYST_CFG $NTHROWS $OUTPUT