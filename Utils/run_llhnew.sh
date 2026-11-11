#!/bin/bash
#SBATCH --job-name=LikelihoodScan
#SBATCH --account=rpp-dune
#SBATCH --output=LikelihoodScan_%j.out
#SBATCH --error=LikelihoodScan_%j.err
#SBATCH --time=5:00:00
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

# ---------------- Move to project directory ---------------- #
cd /scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE || exit 1

# ---------------- Run ---------------- #
./build/Apps/LikelihoodScan /scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/Configs/BeamOffAxis/Fit_new.yaml