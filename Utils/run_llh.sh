#!/bin/bash
#SBATCH --job-name=SigmaVar
#SBATCH --account=rpp-dune
#SBATCH --output=/scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/SigmaVarTemplate/LikelihoodScan_%j.out
#SBATCH --error=/scratch/abipeake/MaCh3DUNE_LukesVersion/MaCh3_DUNE/SigmaVarTemplate/LikelihoodScan_%j.err
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G

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
#./build/Apps/LikelihoodScan /home/abipeake/scratch/July2026/Testwithflatpriors_July29th/cfg/Testwithflatpriors_July29th_chain_0_job_0notflat.yaml
./build/Apps/SigmaVariation /home/abipeake/scratch/MaCh3DUNE_LukesVersion/MaCh3_DUNE/Configs/BeamOffAxis/Fit.yaml #/home/abipeake/scratch/MaCh3DUNE_LukesVersion/MaCh3_DUNE/Configs/BeamOffAxis/FitFD.yaml #/home/abipeake/scratch/MaCh3DUNE_LukesVersion/MaCh3_DUNE/Configs/BeamOffAxis/Fit.yaml