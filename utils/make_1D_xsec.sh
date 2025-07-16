#!/bin/bash

# User-defined variables
RECO_MIN=0.0
RECO_MAX=6.0
RECO_BINS=60
OUTPUT_FILE="RecoNeutrinoEnergy_1Dbinning.yaml"
PARAM_LIST_FILE="RecoNeutrinoEnergy_1Dbinning_param_list.txt"

# Calculate bin width
RECO_STEP=$(echo "($RECO_MAX - $RECO_MIN) / $RECO_BINS" | bc -l)

# Start writing to the output file
echo "Systematics:" > $OUTPUT_FILE
echo "ParameterName, reco_energy_min, reco_energy_max" > $PARAM_LIST_FILE

# Generate systematic blocks
for ((i=0; i<$RECO_BINS; i++)); do
    RECO_LOW=$(echo "$RECO_MIN + $i * $RECO_STEP" | bc -l)
    RECO_HIGH=$(echo "$RECO_LOW + $RECO_STEP" | bc -l)

    PARAM_NAME="reco_energy_$((i+1))"

    cat <<EOL >> $OUTPUT_FILE
- Systematic:
    Sample_Name: ["ND"]
    Error: 0.5
    FlatPrior: false
    KinematicCuts:
    - RecoNeutrinoEnergy:
      - $RECO_LOW
      - $RECO_HIGH
    Names:
      FancyName: $PARAM_NAME
      ParameterName: $PARAM_NAME
    ParameterBounds:
    - 0.0
    - 4.0
    ParameterValues:
      Generated: 1.0
      PreFitValue: 1.0
    StepScale:
      MCMC: 1.0
    Type: Norm
    ParameterGroup: Xsec
EOL

    echo "$PARAM_NAME, $RECO_LOW, $RECO_HIGH" >> $PARAM_LIST_FILE
done

echo "Config file generated: $OUTPUT_FILE"
echo "Parameter list generated: $PARAM_LIST_FILE"
