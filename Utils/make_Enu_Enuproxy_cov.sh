#!/bin/bash

# User-defined variables
ENE_MIN=0.0
ENE_MAX=10.0
EREC_MINUS_ENE_MIN=-60.0
EREC_MINUS_ENE_MAX=0.0
ENE_BINS=50
EREC_MINUS_ENE_BINS=60
OUTPUT_FILE="TrueNeutrinoEnergy_ERecProxy_minus_Enu_0.0_10.0GeV_fullgrid_smallerbins.yaml"
PARAM_LIST_FILE="TrueNeutrinoEnergy_ERecProxy_minus_Enu_0.0_10.0GeV_parameter_list.txt"

# Calculate bin widths
ENE_STEP=$(echo "($ENE_MAX - $ENE_MIN) / $ENE_BINS" | bc -l)
EREC_MINUS_ENE_STEP=$(echo "($EREC_MINUS_ENE_MAX - $EREC_MINUS_ENE_MIN) / $EREC_MINUS_ENE_BINS" | bc -l)

# Start writing to the output file
echo "Systematics:" > $OUTPUT_FILE
echo "ParameterName, ene_min, ene_max, erec_minus_ene_min, erec_minus_ene_max" > $PARAM_LIST_FILE

# Generate systematic blocks
COUNT=1
for ((i=0; i<$ENE_BINS; i++)); do
    for ((j=0; j<$EREC_MINUS_ENE_BINS; j++)); do
        ENE_LOW=$(echo "$ENE_MIN + $i * $ENE_STEP" | bc -l)
        ENE_HIGH=$(echo "$ENE_LOW + $ENE_STEP" | bc -l)
        EREC_MINUS_ENE_LOW=$(echo "$EREC_MINUS_ENE_MIN + $j * $EREC_MINUS_ENE_STEP" | bc -l)
        EREC_MINUS_ENE_HIGH=$(echo "$EREC_MINUS_ENE_LOW + $EREC_MINUS_ENE_STEP" | bc -l)

        PARAM_NAME="enu_erec_minus_enu_$COUNT"

        cat <<EOL >> $OUTPUT_FILE
- Systematic:
    Sample_Name: ["ND"]
    Error: 0.5
    FlatPrior: false
    KinematicCuts:
    - TrueNeutrinoEnergy:
      - $ENE_LOW
      - $ENE_HIGH
    - ERecProxy_minus_Enu:
      - $EREC_MINUS_ENE_LOW
      - $EREC_MINUS_ENE_HIGH
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

        echo "$PARAM_NAME, $ENE_LOW, $ENE_HIGH, $EREC_MINUS_ENE_LOW, $EREC_MINUS_ENE_HIGH" >> $PARAM_LIST_FILE
        COUNT=$((COUNT+1))
    done
done

echo "Config file generated: $OUTPUT_FILE"
echo "Parameter list generated: $PARAM_LIST_FILE"
