#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name=MCChain
#SBATCH --account=def-deborahh
#SBATCH --cpus-per-task=16
#SBATCH --mem=25G
#SBATCH --output=DUNEChain.out
#SBATCH --error=DUNEError.out

module restore

cd /home/acarney/scratch/MaCh3_DUNE/

source build/bin/setup.MaCh3DUNE.sh

./build/bin/Fit Configs/EventRates_Beam.yaml