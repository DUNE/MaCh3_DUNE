#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --job-name=ProcSelec
#SBATCH --account=def-deborahh
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --output=DUNEDiag.out
#SBATCH --error=DUNEDiagE.out

module restore

cd /home/acarney/scratch/
source MaCh3_DUNE/build/bin/setup.MaCh3DUNE.sh
./MaCh3/build/bin/ProcessMCMC MaCh3_DUNE/Configs/EventRates_Beam.yaml MaCh3_DUNE/FHCBothReg/FHCsRegVal1SelecPars.root

#./MaCh3/build/bin/DiagMCMC MaCh3_DUNE/FHCnueReg/RegVal0.1FHCnue.root MaCh3_DUNE/Configs/EventRates_Beam.yaml
#./MaCh3/build/bin/PlotMCMCDiag MaCh3_DUNE/RegVal1_MCMC_Diag.root MaCh3_DUNE/"RegVal1Diag"

# cd /home/acarney/scratch/MaCh3_DUNE
# source build/bin/setup.MaCh3DUNE.sh
# ./build/bin/PostPredDUNE Configs/EventRates_Beam.yaml General:OutputFile:PostPredFHCsRegVal1SelecPars.root