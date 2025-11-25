#!/bin/bash

# Non-uniform bin edges for TrueNeutrinoEnergy (ENE) and isRelativeEnubias (EREC_MINUS_ENE)
#ENE_BINS=(0.0 1.0 1.25 1.5 1.75 2.0 2.25 2.5 2.75 3.0 3.25 3.5 3.75 4.0 5.0 6.0 10.0)
#Enubiasbins=(-2.0 -0.6 -0.581 -0.5595 -0.538 -0.5165 -0.495 -0.4735 -0.452 -0.4305 -0.409 -0.3875 -0.366 -0.3445 -0.323 -0.3015 -0.28 -0.2585 -0.237 -0.2155 -0.194 -0.1725 -0.151 -0.1295 -0.108 -0.0865 -0.065 -0.0435 -0.022 0.0 0.1)

Enubiasbins=(-2.0 -1.3 -0.6 -0.5905 -0.581 -0.57025 -0.5595 -0.54875 -0.538 -0.52725 -0.5165 -0.50575 -0.495 -0.48425 -0.4735 -0.46275 -0.452 -0.44125 -0.4305 -0.41975 -0.409 -0.39825 -0.3875 -0.37675 -0.366 -0.35525 -0.3445 -0.33375 -0.323 -0.31225 -0.3015 -0.29075 -0.28 -0.26925 -0.2585 -0.24775 -0.237 -0.22625 -0.2155 -0.20475 -0.194 -0.18325 -0.1725 -0.16175 -0.151 -0.14025 -0.1295 -0.11875 -0.108 -0.09725 -0.0865 -0.07575 -0.065 -0.05425 -0.0435 -0.03275 -0.022 -0.011 0.0 0.05 0.1)

ENE_BINS=(0.0 0.5 1.0 1.125 1.25 1.375 1.5 1.625 1.75 1.875 2.0 2.125 2.25 2.375 2.5 2.625 2.75 2.875 3.0 3.125 3.25 3.375 3.5 3.625 3.75 3.875 4.0 4.5 5.0 5.5 6.0 8.0 10.0)


#ENE_BINS=(0.0 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.25 2.5 2.75 3.0 3.25 3.5 3.75 4.0 5 10.0)
#EREC_MINUS_ENE_BINS=(-1.0 -0.75 -0.5 -0.49902296 -0.47777235 -0.4579873 -0.43820225 -0.41988276 -0.40156326 -0.38470933 -0.3678554 -0.35246702 -0.33634587 -0.32169028 -0.30630191 -0.29237909 -0.2777235 -0.26380068 -0.24987787 -0.23668784 -0.22276502 -0.20957499 -0.19638495 -0.18319492 -0.17000489 -0.15681485 -0.14362482 -0.13043478 -0.11651197 -0.10332193 -0.08939912 -0.07547631 -0.06082071 -0.04616512 -0.03150953 -0.01612115 0.0 0.005)
OUTPUT_FILE="Enu_Enubiasxsec_enubiasandnuetwiceassmall.yaml"
PARAM_LIST_FILE="Enu_Enubiasxsec_params.txt"

# Start writing to the output files
echo "Systematics:" > $OUTPUT_FILE
echo "ParameterName, ene_min, ene_max, erec_minus_ene_min, erec_minus_ene_max" > $PARAM_LIST_FILE

# Generate systematic blocks using the non-uniform bin edges
COUNT=1
for ((i=0; i<${#ENE_BINS[@]}-1; i++)); do
    for ((j=0; j<${#Enubiasbins[@]}-1; j++)); do
        ENE_LOW=${ENE_BINS[$i]}
        ENE_HIGH=${ENE_BINS[$((i+1))]}
        EREC_MINUS_ENE_LOW=${Enubiasbins[$j]}
        EREC_MINUS_ENE_HIGH=${Enubiasbins[$((j+1))]}

        PARAM_NAME="enu_enubias_$COUNT"

        cat <<EOL >> $OUTPUT_FILE
- Systematic:
    Sample_Name: ["ND*", "FD*"]
    Error: 0.5
    FlatPrior: false
    KinematicCuts:
    - TrueNeutrinoEnergy:
      - $ENE_LOW
      - $ENE_HIGH
    - Enubias:
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
