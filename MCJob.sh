#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --job-name=FHCsSelecTwen
#SBATCH --account=def-deborahh
#SBATCH --cpus-per-task=16
#SBATCH --mem=25G
#SBATCH --output=DUNEChain2.out
#SBATCH --error=DUNEError2.out

module restore

cd /home/acarney/scratch/MaCh3_DUNE/

source ../MaCh3/build/bin/setup.MaCh3.sh
source build/bin/setup.MaCh3DUNE.sh

./build/bin/Fit Configs/EventRates_Beam.yaml #General:FittingAlgorithm:DelayedMCMC