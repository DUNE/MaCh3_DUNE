#!/bin/bash
#SBATCH --job-name=MCMCWorkflow
#SBATCH --account=rpp-nilic
#SBATCH --output=MCMCWorkflow_%j.out
#SBATCH --error=MCMCWorkflow_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8    # adjust if your code is multithreaded
#SBATCH --mem=400G

# ---------------- Modules ---------------- #
module purge
module load StdEnv/2020
module load cmake/3.18.4
module load gcc/9.3.0
module load python/3.9.6
module load root/6.26.06

# ---------------- Setup environment ---------------- #
source /scratch/abipeake/MaCh3DUNE_OffAxisnew/MaCh3_DUNE/build/bin/setup.MaCh3.sh
source /scratch/abipeake/MaCh3DUNE_OffAxisnew/MaCh3_DUNE/build/bin/setup.MaCh3DUNE.sh

# ---------------- Paths ---------------- #
# Path to your binaries
Apps_DIR=/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/build/Apps/AbiP/Plotting/


CONFIG_YAML=/home/abipeake/scratch/March2026/Monday9thMarch_OnAxis_templateparams_postadaptive/cfg/posterior.yaml #/scratch/abipeake/October_chains/Tuesday_041125_enubias_newattempt_withinitialproposal/cfg/FD_Posteriorpredictive.yaml #/scratch/abipeake/October_chains/Tuesday_041125_enubias_newattempt_withinitialproposal/cfg/FD_Posteriorpredictive.yaml
 #/scratch/abipeake/October_chains/Tuesday_041125_enubias_newattempt_withinitialproposal/cfg/FD_Posteriorpredictive.yaml #/scratch/abipeake/October_chains/Tuesday_041125_enubias_newattempt_adaptingfor4millionseteps_robbinsmonroe/cfg/FD_Postpred.yaml #/scratch/abipeake/October_chains/Saturday_enubiasenu_newchain/cfg/FD_Postpred.yaml #/scratch/abipeake/October_chains/Saturday_enubiasenu_newchain/cfg/FD_Postpred.yaml #/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/TestingnoRM/FDpostpred.yaml #/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/TestingnoRM/Tuesday_041125_enubias_newattempt_withinitialproposal_chain_0_job_0copy_recovered.root #/scratch/abipeake/October_chains/Tuesday_041125_enubias_newattempt_adaptingfor4millionseteps_robbinsmonroe/cfg/FD_Postpred.yaml #/scratch/abipeake/October_chains/Saturday_EnubiasEnu_with12moffaxis/cfg/FD_Postpred.yaml #/scratch/abipeake/October_chains/Saturday_enubiasenu_newchain/cfg/FD_Postpred.yaml #/scratch/abipeake/October_chains/Saturday_enubiasenu_newchain/cfg/FD_Postpred.yaml #/scratch/abipeake/Off_axis_chains_new/ACTUALLYIS_ENUBIAS_ELEPEREC_postadaptive_offaxis_v2/cfg/FDPostpred.yaml #/scratch/abipeake/October_chains/THURS_newadaptive_enubias/cfg/FD_Posterior.yaml #/scratch/abipeake/Off_axis_chains_new/ACTUALLYIS_ENUBIAS_ELEPEREC_postadaptive_offaxis_v2/cfg/FDPostpred.yaml #/scratch/abipeake/October_chains/Wed_EnuWnubias_adaptivefit/cfg/FD_Postpredictive.yaml #/scratch/abipeake/Off_axis_chains_new/afteradaptive_q0q3_ptpzenu/cfg/FD_postpred.yaml #/scratch/abipeake/October_chains/Wed_EnuWnubias_adaptivefit/cfg/FD_Postpredictive.yaml #/scratch/abipeake/Off_axis_chains_new/ACTUALLYIS_ENUBIAS_ELEPEREC_postadaptive_offaxis_v2/cfg/FDPostpred.yaml #/scratch/abipeake/Off_axis_chains_new/Onaxis_EventRates_Elep_erec_afteradaptive/cfg/FD_posteriorpredictive.yaml #/scratch/abipeake/Off_axis_chains_new/RunPlan_10_90_aadapted_part2/cfg/FD_Postpred.yaml #/scratch/abipeake/Off_axis_chains_new/Onaxis_EventRates_Elep_erec_afteradaptive/cfg/FD_posteriorpredictive.yaml #/scratch/abipeake/Off_axis_chains_new/ACTUALLYIS_ENUBIAS_ELEPEREC_postadaptive_offaxis_v2/cfg/FDPostpred.yaml #/scratch/abipeake/Off_axis_chains_new/RunPlan_10_90_aadapted_part2/cfg/FD_Postpred.yaml #/scratch/abipeake/October_chains/EventRates_1D_30bins_adapted_new2/cfg/FDPostpred.yaml #/scratch/abipeake/October_chains/EventRates_1D_30bins_adapted_new2/cfg/FDPostpred.yaml #/scratch/abipeake/October_chains/EventRates_1D_70bins_afteradapted_wed/cfg/FD_postpred.yaml #/scratch/abipeake/October_chains/EventRates_1D_50bins_afteradapted_postadap2/cfg/FDPostpred.yaml #/scratch/abipeake/Off_axis_chains_new/Onaxis_EventRates_Elep_erec_afteradaptive/cfg/FD_posteriorpredictive.yaml #/scratch/abipeake/Off_axis_chains_new/RunPlan_10_90_aadapted_part2/cfg/FD_Postpred.yaml #/scratch/abipeake/Off_axis_chains_new/ACTUALLYIS_ENUBIAS_ELEPEREC_postadaptive_offaxis_v2/cfg/FDPostpred.yaml #/scratch/abipeake/Off_axis_chains_new/RunPlan_10_90_aadapted_part2/cfg/FD_Postpred.yaml #/scratch/abipeake/Off_axis_chains_new/Onaxis_EventRates_Elep_erec_afteradaptive/cfg/FD_posteriorpredictive.yaml #/scratch/abipeake/Off_axis_chains_new/ACTUALLYIS_ENUBIAS_ELEPEREC_postadaptive_offaxis_v2/cfg/FDPostpred.yaml #/scratch/abipeake/Off_axis_chains_new/RunPlan_10_90_afteradapted/cfg/FD_Postpred.yaml #/scratch/abipeake/Off_axis_chains_new/RunPlan_90_10_afteradapted/cfg/FD_Postpred.yaml #/scratch/abipeake/Off_axis_chains_new/Onaxis_EventRates_Elep_erec_afteradaptive/cfg/FD_posteriorpredictive.yaml #/scratch/abipeake/October_chains/EventRates_1D_30bins_afteradapted_afteradapwed/cfg/FD_Postpred.yaml #/scratch/abipeake/Off_axis_chains_new/ACTUALLYIS_ENUBIAS_ELEPEREC_postadaptive_offaxis_v2/cfg/FDPostpred.yaml  #/scratch/abipeake/October_chains/EventRates_1D_30bins_afteradapted_afteradapwed/cfg/FD_Postpred.yaml #/scratch/abipeake/Off_axis_chains_new/Onaxis_EventRates_Elep_erec_afteradaptive/cfg/FD_posteriorpredictive.yaml #/scratch/abipeake/Off_axis_chains_new/RunPlan_90_10_afteradapted/cfg/FD_Postpred.yaml #/scratch/abipeake/Off_axis_chains_new/ACTUALLYIS_ENUBIAS_ELEPEREC_postadaptive_offaxis_v2/cfg/FDPostpred.yaml #/scratch/abipeake/Off_axis_chains_new/ACTUALLYIS_ENUBIAS_ELEPEREC_postadaptive_offaxis_v2/cfg/FDPostpred.yaml #/scratch/abipeake/October_chains/EventRates_1D_70bins_afteradapted_wed/cfg/FD_postpred.yaml #/scratch/abipeake/October_chains/EventRates_1D_70bins_afteradapted_wed/cfg/FD_postpred.yaml #/scratch/abipeake/October_chains/EventRates_1D_30bins_adapted_new2/cfg/FDPostpred.yaml #/scratch/abipeake/Off_axis_chains_new/ACTUALLYIS_ENUBIAS_ELEPEREC_postadaptive_offaxis_v2/cfg/FDPostpred.yaml #/scratch/abipeake/Off_axis_chains_new/TUesday_afteradaptive_Elep_Erec_EnuEnubias_offaxis/cfg/FD_Postpred.yaml #/scratch/abipeake/Off_axis_chains_new/Onaxis_EventRates_Elep_erec_afteradaptive/cfg/FD_posteriorpredictive.yaml
 #/scratch/abipeake/Off_axis_chains_new/TUesday_afteradaptive_Elep_Erec_EnuEnubias_offaxis/cfg/FD_Postpred.yaml #/scratch/abipeake/Off_axis_chains_new/Onaxis_EventRates_Elep_erec_afteradaptive/cfg/FD_posteriorpredictive.yaml

echo "Running Posterior Predictive code"
srun "$Apps_DIR/PosteriorPredictive_forND" "$CONFIG_YAML"


