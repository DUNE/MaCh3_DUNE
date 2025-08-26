#!/bin/bash
#set -x

#
# Submit job script, bash style
# Because I prefer sed and shizz
#
# This script will run multiple MCMC chains for a SINGLE fit
# on the Emerald cluster. Typically when the jobs finish you can combine
# or reduce the .root files into a single file for analysis.
#
# Current good clusters:  STFC Emerald, 
#                         Compute Canada (Cedar, Graham, Beluga, Niagara, Guillimin, Hades, Helios)
#                         Imperial HPC
#                         Imperial HEP (heppc105, heppc205, lt2gpu00)

# A function to calculate the time from seconds to HH:MM:SS
function SecToH
{
  if [[ "$#" -ne 1 ]]; then
    echo "Second to hour convertor needs one argument"
    exit -1
  fi

  num=$1
  ss=00
  mm=00
  hh=00
  if ((num>59)); then
    ((ss=num%60))
    ((num=num/60))
    if ((num>59)); then
      ((mm=num%60))
      ((num=num/60))
      if ((num>23)); then
        ((hh=num))
      else
        ((hh=num))
      fi
    else
      ((mm=num))
    fi
  else
    ((ss=num))
  fi

  # Now set the variables to 00 form
  ss=$(printf "%.2i" $ss)
  mm=$(printf "%.2i" $mm)
  hh=$(printf "%.2i" $hh)

  echo "$hh:$mm:$ss"
}

# Add this at the top of your script
timestamp=$(date +"%Y%m%d_%H%M%S")

if [[ "$#" -ne 6 ]]; then
  echo "I need 6 parameters, you gave $#"
  echo "./submit_DUNE_niagara.sh  EXECUTABLE  CONFIG_FILE  FIT_NAME  N_CHAINS  N_STEPS JOB"
  exit -1
fi

# First argument is executable
EXE=$1

# Second argument is base config file
CFG=$2

# Third argument is output name
FitName=$3

# Fourth argument is number of jobs
NCHAINS=$4
case $NCHAINS in
  ''|[!0-9]*) echo "N_CHAINS = $NCHAINS NOT A NUMBER"; exit -1 ;;
esac

# Fifth argument is number of steps
NSTEPS=$5
case $NSTEPS in
  ''|[!0-9]*) echo "N_STEPS = $NSTEPS NOT A NUMBER"; exit -1 ;;
esac

# Fifth argument is number of steps
JOB=$6
case $JOB in
  ''|[!0-9]*) echo "JOB = $JOB NOT A NUMBER"; exit -1 ;;
esac


# Number of threads to use for multi-threading jobs
# Change these if you want!
# Number of threads, memory, and GPUs
NTHREADS=20
RAM=10
NGPU=0
SEC_PER_STEP=0.05


###########################################################
# DO NOT EDIT BELOW HERE
###########################################################

# Different settings for the different clusters
ScratchDir="/scratch/$USER/"
##################


##################
# Set up stuff depending on which cluster we're on
# Emerald
IsEmerald=false
# Compute Canada Guilliumin
IsCC=false
# Compute Canada Helios
IsHelios=false
# Imperial HEP
IsICHEPPC=false
# Imperial HEP
IsICHEPQ=false
# Imperial HPC
IsICHPC=false
# RHUL Linappserv2
IsLinapp=false
# Cedar
IsCedar=false
# Niagara
IsNiagara=false
IsFir=false

echo "---------------"
# Emerald
if [[ $HOSTNAME == *gpu.rl.ac.uk ]]; then
  IsEmerald=true
  echo "Found Emerald cluster, setting scratch"
  ScratchDir="/work/scratch"

  # ComputeCanada Guillimin
elif [[ $HOSTNAME == lg-1r1[47]-n[0][1-4] ]]; then
  IsCC=true
  echo "Found Guilliumin cluster, setting scratch"
  ScratchDir="/gs/project/sab-064-aa/${USER}/scratch"
  # Guillimin is very fast!
  SEC_PER_STEP=0.05

  # ComputeCanada Helios
elif [[ $HOSTNAME == helios[0-9] ]]; then
  IsHelios=true
  echo "Found Helios cluster, setting scratch"
  ScratchDir="/rap/sab-064-aa/${USER}/scratch"
  # Use less threads on Helios because how the nodes are set up
  NTHREADS=5
  # Also Helios is just slow...
  SEC_PER_STEP=0.06

elif [[ $HOSTNAME == cedar* || $HOSTNAME == login1 || $HOSTNAME == gra-login* || $HOSTNAME == beluga[0-9]* ]]; then
  IsCedar=true
  echo "Found Cedar, Graham or Beluga cluster, setting scratch"
  ScratchDir="${SCRATCH}"
  SEC_PER_STEP=1.0
  # Currently a bug in SLURM scheduler which screws up RAM usage
  # Cedar people have been notified, recommended upping RAM usage
  RAM=7

  # Niagara
elif [[ $HOSTNAME == nia-login* ]]; then
  IsNiagara=true
  echo "Found Niagara cluster, setting scratch"
  ScratchDir="${SCRATCH}"
  SEC_PER_STEP=0.16
  # Niagara automatically gives 40 threads per job, so use them all!
  NTHREADS=40

  # Imperial heppc GPU
elif [[ $HOSTNAME == @(heppc105|heppc205).hep.ph.ic.ac.uk ]]; then
  IsICHEPPC=true
  echo "Found Imperial HEP PC, setting scratch"
  ScratchDir="/vols/t2k/users/${USER}/scratch/MaCh3"

  # Imperial lt2gpu00 GPU
elif [[ $HOSTNAME == "lt2gpu00.hep.ph.ic.ac.uk" ]]; then
  IsICHEPQ=true
  echo "Found Imperial HEP Queue, setting scratch"
  ScratchDir="/vols/t2k/users/${USER}/scratch/MaCh3"
  SEC_PER_STEP=0.05

  # Imperial HPC
elif [[ $HOSTNAME == *[0-9]-internal ]]; then
  IsICHPC=true
  echo "Found Imperial HPC cluster, setting scratch"
  ScratchDir="$WORK/scratch/${USER}"
  SEC_PER_STEP=0.05

elif [[ $HOSTNAME == "linappserv2.pp.rhul.ac.uk" ]]; then
  IsLinapp=true
  echo "Found RHUL Linapp, setting scratch"
  ScratchDir="/home/${USER}/s3/OA2019/EBDial/EBDial"
  SEC_PER_STEP=0.6

else 
  echo "Couldn't find correct cluster, so can't set scratch"
  echo "Edit me please"
  echo "Running on: $HOSTNAME"
  exit -1
fi

##################
# End cluster-specifics


# Number of threads, memory, and GPUs
NTHREADS=20
RAM=7
NGPU=0
SEC_PER_STEP=0.3

# Number of steps must be defined before this block
# Example: NSTEPS=2000000

##################
echo "You've given me $NTHREADS CPU threads and ${NGPU} GPUs with $NSTEPS steps"
echo "---------------"
echo "Setting walltime automatically, assuming $SEC_PER_STEP s/step and 30 min start-up..."

# Calculate walltime in seconds
WALLTIME=$(echo "${SEC_PER_STEP} * ${NSTEPS} + 1800" | bc)  # 1800 = 30 min setup buffer

# Optional: Cap max walltime to 8 hours (28800 sec) for better scheduling
MAX_WALLTIME=28800  #86400 #28800 #86400 #28800
if (( $(echo "$WALLTIME > $MAX_WALLTIME" | bc -l) )); then
  echo "⚠️  Walltime capped at 8 hours to improve queue time"
  WALLTIME=$MAX_WALLTIME
fi

# Convert to integer seconds
WALLTIME_INT=$(printf "%.0f" "$WALLTIME")

# Convert seconds to HH:MM:SS format
SecToH() {
  local total_seconds=$1
  local hours=$((total_seconds / 3600))
  local minutes=$(((total_seconds % 3600) / 60))
  local seconds=$((total_seconds % 60))
  printf "%02d:%02d:%02d\n" $hours $minutes $seconds
}
WALLTIME_CC=$(SecToH "$WALLTIME_INT")
WALLTIME_EM="${WALLTIME_CC%:*}"  # HH:MM for Emerald (if needed)

# Output for logs
echo "Final walltime estimate:"
echo " - In seconds     : $WALLTIME"
echo " - In HH:MM:SS    : $WALLTIME_CC"

# Check if suggested WALLTIME is OK with user
#while true; do
#  read -p "Walltime = ${WALLTIME_CC}, OK (y/n)?" yn
#  case $yn in
#    [Yy]* ) break;;
#    [Nn]* ) exit;;
#    * ) echo "Please answer...";;
#  esac
#done
#echo "---------------"
##################
# End walltime check



# Find MaCh3 folder and setup script
# Assuming we're running in MaCh3 folder
##################
MaCh3Dir=$(pwd -P)
echo $MaCh3Dir | grep -q 'MaCh3'
greprc=$?
if [ $greprc -ne 0 ]; then
  echo "Coulnd't find MaCh3 in current pwd"
  echo "pwd = $(pwd -P)"
  echo "Have to quit, sorry"
  exit -1
fi
MaCh3Dir=${MaCh3Dir%%/clusters}

# Check that the executable exists
#if [ ! -e ${EXE} ]; then
#  echo "Did not find executable ${EXE} in ${MaCh3Dir}..."
#  echo "Have you built MaCh3 yet?"
#  exit -1
#fi
##################


##################
# Check that the config file exists
if [ ! -e ${CFG} ]; then
  echo "Did not find config file ${CFG}"
  exit -1
fi
##################


##################
# Output directory
#OutputDir=${ScratchDir}/"FinalChains/"${FitName}
OutputDir=${ScratchDir}/"Off_axis_chains/"${FitName}

counter=0
# Check to see if this path exists, display a warning if exists
#while [ -d $OutputDir ]; do
#  echo "Output directory $OutputDir already exists"
#  while true; do
#    read -p 'Add number to directory name? (y/n) _OR_ Delete current directory (d)' yn
#    case $yn in
#      [Yy]* ) break;;
#      [Nn]* ) exit;;
#      [Dd]* ) rm -rfv $OutputDir; break;;
#      * ) echo "Please answer...";;
#    esac
#  done
#  OutputDir+=_$counter
#  ((counter++))
#done

# Make a temp folder that will store all the intermediate scripts, configs and logs
scriptDir=${OutputDir}/scripts
cfgDir=${OutputDir}/cfg
logDir=${OutputDir}/log

# if this is the first job then make the Output Directories
#if [[ $JOB -eq 0 ]]
#then
 # mkdir -p $OutputDir
 # mkdir -p $scriptDir
  #mkdir -p $cfgDir
 ## mkdir -p $logDir
#fi
##################
# Only create the directories if they don't already exist
if [[ ! -d $scriptDir || ! -d $cfgDir || ! -d $logDir ]]; then
  echo "Creating necessary output directories in $OutputDir"
  mkdir -p "$scriptDir" "$cfgDir" "$logDir"
else
  echo "Output directories already exist. Reusing: $OutputDir"
fi



##################
# Special treatment of picky task spooler at Imperial
# Essentially set the path, logging directory, socket and simultaneous slots
if [[ $IsICHEP == true ]]; then
  # Export the task spooler
  export PATH=/vols/t2k/users/$USER/software/ts-1.0/bin:${PATH}
  # Change the TMPDIR env variable which we write to
  export TMPDIR=${logDir}
  export TS_SOCKET=/tmp/ts.socket
  export TS_SLOTS=1
fi
##################


##################
# Loop over each job and make the config file and submission files
for ((i = 0 ; i < $NCHAINS ; i++)); do
  JobName=${FitName}_chain_${i}_job_${JOB}
  cfgFileName=${cfgDir}/${JobName}.yaml

  echo "cfgDir='$cfgDir'"
  echo "CFG='$CFG'"
  echo "OutputDir='$OutputDir'"
  echo "cfgFileName='$cfgFileName'"


  # Copy the "template" config file to new location
  cp "$CFG" "$cfgFileName"

  # Replace the output name
  sed -i "s|^\( *OutputFile:\).*|\1 \"${OutputDir}/${JobName}.root\"|g" "$cfgFileName"
  grep "OutputFile" "$cfgFileName"

  # Replace NSteps in YAML file (preserve comments and spacing)
  sed -i "s/^\([[:space:]]*NSteps:[[:space:]]*\)[0-9]\+\(.*\)\$/\1$NSTEPS\2/" "$cfgFileName"
  grep "NSteps" "$cfgFileName"
done


  #############################################
  #if you're continuing from the previous chain
  #############################################
  if [[ $JOB -gt 0 ]]
  then
    OLD_JOB=$(( $JOB-1 ))
    sed -i "s|^  StartFromPos.*|  StartFromPos: true|g" $cfgFileName
    OldJobName=${FitName}_chain_${i}_job_${OLD_JOB}.root
    sed -i "s|^  PosFileName.*|  PosFileName: \"${OutputDir}/${OldJobName}\"|g" $cfgFileName
    if [[ -e "${OutputDir}/${OldJobName}" ]]
    then
      echo "Previous chain exists!"
      echo "${OutputDir}/${OldJobName}" 
    else
      echo "Cannot find previous chain ${OutputDir}/${OldJobName}"
      exit -1
    fi
  fi



  # Now make the .sh script that we submit to the cluster
  ScriptFileName=${scriptDir}/${JobName}.sh

  # Write the temporary submission file
  executable="${EXE} ${cfgFileName}"
  echo "#!/bin/bash" > $ScriptFileName

  # Write the OMP_NUM_THREADS variable
  echo "export OMP_NUM_THREADS=${NTHREADS}" >> $ScriptFileName
  # Save the MaCh3 directory into the file
  echo "export MACH3=${MaCh3Dir}" >> $ScriptFileName

  # Have the script cat the config so we have a log of the exact config file that was run
  echo "cat ${cfgFileName}" >> $ScriptFileName
  echo "cd \${MACH3}" >> $ScriptFileName
  # Load up cluster defaults
  echo "source setup.sh" >> $ScriptFileName
  echo "source setup_dune_cedar_env.sh" >> $ScriptFileName
  echo "source build/bin/setup.MaCh3.sh" >> $ScriptFileName
  echo "source build/bin/setup.MaCh3DUNE.sh" >> $ScriptFileName
  echo "PATH=../build/src:$PATH" >> $ScriptFileName
  # Run exec
  echo "${executable}" >> $ScriptFileName
  # Make executable
  chmod 744 $ScriptFileName

  # Now temp bash file is written and we can submit the contents
  stdout="-o ${logDir}/${JobName}.log"
  stderr="-e ${logDir}/${JobName}.err"

  # For Emerald we use bsub
  if [[ ${IsEmerald} == true ]]; then
    bsubmit="bsub"
    jobid="-J ${JobName}"
    # General options
    # e.g. WALLTIME etc
    bsubopt="-n ${NTHREADS} -N -W ${WALLTIME_EM}"

    # The final command to submit
    bsubmit="${bsubmit} ${bsubopt} ${stdout} ${stderr} ${jobid} ${ScriptFileName}"
    eval $bsubmit

    # For ComputeCanada we use qsub
  elif [[ ${IsCC} == true ]]; then
    qsubmit="qsub"
    RAMMB=$(echo "${RAM} * 1000" | bc)
    RAMMB=${RAMMB%%.*}
    qsubopt="-l nodes=1:ppn=${NTHREADS}:gpus=${NGPU} -l pmem=${RAMMB}mb -l walltime=${WALLTIME_CC}"
    qsubmit="${qsubmit} ${qsubopt} ${stdout} ${stderr} ${ScriptFileName}"
    eval $qsubmit

  elif [[ ${IsCedar} == true ]]; then
    qsubmit="sbatch"
    # Need to specify the user, walltime, gpu, cpu
    RAMMB=$(echo "${RAM} * 1000" | bc)
    RAMMB=${RAMMB%%.*}
    echo "Asking for ${RAMMB}M of RAM"
    if [[ $JOB -gt 15 ]]
    then
      OLD_JOB=$(( $JOB-1 ))
      OldJobName=${FitName}_chain_${i}_job_${OLD_JOB}
      #PrevJobPID=$(squeue -u liban --format="%.18i %.9P %.60j %.8u %.8T %.10M %.9l %.6D %R" | grep "${OldJobName}" | awk '{print $1}')
      #qsubopt="--account=rpp-nilic --time=${WALLTIME_CC} --cpus-per-task=${NTHREADS} --mem=${RAMMB}M --dependency=afterany:${PrevJobPID}"
      #echo "This job will wait for job with name: ${OldJobName} with PID: "${PrevJobPID}""
      PrevJobPID=$(squeue -u liban --format="%.18i %.9P %.60j %.8u %.8T %.10M %.9l %.6D %R" | grep "${OldJobName}" | awk '{print $1}')
      if [[ -n "$PrevJobPID" ]]; then
        qsubopt="--account=def-nilic --time=${WALLTIME_CC} --cpus-per-task=${NTHREADS} --mem=${RAMMB}M --dependency=afterany:${PrevJobPID}"
        echo "This job will wait for job with name: ${OldJobName} with PID: ${PrevJobPID}"
      else
  echo "⚠️  No previous job found for ${OldJobName}; submitting without dependency"
  qsubopt="--account=def-nilic --time=${WALLTIME_CC} --cpus-per-task=${NTHREADS} --mem=${RAMMB}M"
fi

    else
      qsubopt="--account=def-nilic --time=${WALLTIME_CC} --cpus-per-task=${NTHREADS} --mem=${RAMMB}M"
    fi
    # Also different way of specifiying output
    stdout="--output ${logDir}/${JobName}_%j.log"
    stderr="--error ${logDir}/${JobName}_%j.err"
    qsubmit="${qsubmit} ${qsubopt} ${stdout} ${stderr} ${ScriptFileName}"
    eval $qsubmit

  elif [[ ${IsNiagara} == true ]]; then
    qsubmit="sbatch"
    # Need to specify the user, walltime, gpu, cpu
    RAMMB=$(echo "${RAM} * 1000" | bc)
    RAMMB=${RAMMB%%.*}
    # Need some special requests on Niagara
    qsubopt="--account=def-deborahh --time=${WALLTIME_CC} --cpus-per-task=${NTHREADS} --nodes=1"
    # Also different way of specifiying output
    stdout="--output ${logDir}/${JobName}_%j.log"
    stderr="--error ${logDir}/${JobName}_%j.err"
    qsubmit="${qsubmit} ${qsubopt} ${stdout} ${stderr} ${ScriptFileName}"
    eval $qsubmit

    # On helios you need -A to specify submitter and you also need to specify the queue
    # Also use msub
    # Seems like we can't submit more than 5 threads? prefer k20 queue because we get more CPU
  elif [[ ${IsHelios} == true ]]; then
    qsubmit="msub"
    qsubopt="-l nodes=1:gpus=1 -l walltime=${WALLTIME_CC} -A sab-064-aa -N ${JobName}"
    qsubmit="${qsubmit} ${qsubopt} ${stdout} ${stderr} ${ScriptFileName}"
    eval $qsubmit

    # For Imperial heppc we use TaskSpooler
    # Needs to be submitted from GPU machine
  elif [[ ${IsICHEPPC} == true ]]; then
    # Give a good name (-L), split output into stdout and stderr (-E)
    taskspool="ts"
    taskopt="-L ${JobName} -E"

    # The final command to submit
    taskspool="${taskspool} ${taskopt} ${ScriptFileName}"
    eval $taskspool

    # For Imperial lt2gpu00 we use qsub now
    # Needs to be submitted from GPU machine
  elif [[ ${IsICHEPQ} == true ]]; then
    qsubmit="qsub"
    # Send to GPU queue, choose one GPU, set walltime to HH:MM:SS, request 4GB, request multi-thread
    qsubopt="-q gpu.q -l gpu=1 -l h_rt=${WALLTIME_CC} -l h_vmem=${RAM}G -pe hep.pe ${NTHREADS}"
    # Put together the submit command
    qsubmit="${qsubmit} ${qsubopt} ${stdout} ${stderr} ${ScriptFileName}"
    eval $qsubmit

  elif [[ ${IsICHPC} == true ]]; then
    qsubmit="qsub"
    RAMMB=$(echo "${RAM} * 1000" | bc)
    RAMMB=${RAMMB%%.*}
    qsubopt="-l select=1:ncpus=${NTHREADS}:ngpus=${NGPU}:mem=${RAMMB}mb -q gpgpu -l walltime=${WALLTIME_CC}"
    # Can also give gpu_type=P100 for spanking new P100 cards
    qsubmit="${qsubmit} ${qsubopt} ${stdout} ${stderr} ${ScriptFileName}"
    eval $qsubmit

 elif [[ ${IsLinapp} == true ]]; then
    qsubmit="qsub"
    RAMMB=$(echo "${RAM} * 1000" | bc)
    RAMMB=${RAMMB%%.*}
    #qsubopt="-l select=1:ncpus=${NTHREADS}:ngpus=${NGPU}:mem=${RAMMB}mb -qlong -l walltime=${WALLTIME_CC}"
    qsubopt="-q long -l  pmem=8gb,pvmem=8gb -l walltime=${WALLTIME_CC}"
    qsubmit="${qsubmit} ${qsubopt} ${stdout} ${stderr} ${ScriptFileName}"
    eval $qsubmit

  fi

done
##################

echo "All ${NCHAINS} chains submitted (This is Job number $JOB for these chains!)"
echo "Used ${CFG} template config, doing ${FitName} with ${NSTEPS} steps each"
echo "Submitted with GPU and ${NTHREADS} CPUs"