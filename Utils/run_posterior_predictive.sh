#!/bin/bash
#SBATCH --account=def-nilic
#SBATCH --job-name=posterior_q0q3
#SBATCH --output=logs/posterior_q0q3_%j.out
#SBATCH --error=logs/posterior_q0q3_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G

# Load modules or environment if needed
# module load gcc cmake root
# source setup_maCh3.sh   # if you have a custom environment script

# Move to your build directory (if necessary)
cd $SLURM_SUBMIT_DIR

# Run your code
./build/Apps/Abi_Posteriorpredicitive \
  /scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/CERN_PosteriorPredictives/PosteriorPredictive_q0q3_Erecyrec.yaml
