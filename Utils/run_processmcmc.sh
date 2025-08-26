#!/bin/bash
#SBATCH --job-name=ProcessMCMC
#SBATCH --account=def-nilic
#SBATCH --output=ProcessMCMC_%j.out
#SBATCH --error=ProcessMCMC_%j.err
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1          # change to 8 if your code is parallel
#SBATCH --mem=80G                  # increase memory to avoid OOM

module purge
module load StdEnv/2020
module load cmake/3.18.4
module load gcc/9.3.0              # or intel/2020.1.217 if required
module load python/3.9.6
module load root/6.26.06

source build/bin/setup.MaCh3.sh
source build/bin/setup.MaCh3DUNE.sh

srun ./build/bin/DiagMCMC \
    /scratch/abipeake/August_tesing/monday_noadaption_q0q3_ptpzenu2/monday_noadaption_q0q3_ptpzenu2_chain_0_job_0.root \
    /scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/TutorialConfigs/TutoriaDiagConfigq0q3.yaml

#/build/bin/PlotMCMCDiag /scratch/abipeake/August_tesing/monday_noadaption_q0q3_ptpzenu2/monday_noadaption_q0q3_ptpzenu2_chain_0_job_0_MCMC_Diag.root
#./build/bin/ProcessMCMC configs/TutorialConfigs/TutorialDiagConfig.yaml /scratch/abipeake/August_tesing/monday_noadaption_q0q3_ptpzenu2/monday_noadaption_q0q3_ptpzenu2_chain_0_job_0.root
#./build/bin/ProcessMCMC configs/TutorialConfigs/TutorialDiagConfig.yaml new_nonadaptive_enuenubias_eleperc_chain_0_job_0.root
#./build/bin/GetPostfitParamPlots /scratch/abipeake/August_tesing/Tuesday_Enubiasrel_enu_erec_yrec_withthrows/Tuesday_Enubiasrel_enu_erec_yrec_withthrows_chain_0_job_0_Process.root
#./build/bin/PlotMCMCDiag /scratch/abipeake/August_tesing/new_nonadaptive_enuenubias_eleperc/new_nonadaptive_enuenubias_eleperc_chain_0_job_0_MCMC_Diag.root
#/build/bin/ProcessMCMC configs/TutorialConfigs/TutorialDiagConfig.yaml /scratch/abipeake/August_tesing/new_nonadaptive_enuenubias_eleperc/new_nonadaptive_enuenubias_eleperc_chain_0_job_0.root