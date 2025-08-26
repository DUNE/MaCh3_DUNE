#!/bin/bash

# Non-uniform bin edges for TrueNeutrinoEnergy (ENE) and isRelativeEnubias (EREC_MINUS_ENE)
ENE_BINS_ORIG=(0.0 0.9184172 1.52906693 2.02735711 2.48656571 2.94088911 3.4147533 3.95701026 4.64582316 5.77918906 10.0)
EREC_MINUS_ENE_BINS_ORIG=(-1.0 -0.60527601 -0.51794333 -0.46268197 -0.42222277 -0.39015144 -0.36350757 -0.34130435 -0.32156815 -0.30380557 \
-0.28801661 -0.27370787 -0.26087934 -0.24854421 -0.2376893 -0.22683439 -0.21696629 -0.208085 -0.19920371 -0.19032242 \
-0.18242794 -0.17453346 -0.16713239 -0.15973131 -0.15282364 -0.14591597 -0.13950171 -0.13308745 -0.12667318 -0.12025892 \
-0.11433806 -0.10792379 -0.10200293 -0.09608207 -0.09016121 -0.08424035 -0.07831949 -0.07239863 -0.06697118 -0.06105032 \
-0.05512946 -0.0492086 -0.04279433 -0.03687347 -0.03095261 -0.02453835 -0.01812408 -0.01121641 -0.00430874 0.00259893 0.01)


ENE_BINS=(0.0 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.25 2.5 2.75 3.0 3.25 3.5 3.75 4.0 5 10.0)
EREC_MINUS_ENE_BINS=(-1.0 -0.75 -0.5 -0.49902296 -0.47777235 -0.4579873 -0.43820225 -0.41988276 -0.40156326 -0.38470933 -0.3678554 -0.35246702 -0.33634587 -0.32169028 -0.30630191 -0.29237909 -0.2777235 -0.26380068 -0.24987787 -0.23668784 -0.22276502 -0.20957499 -0.19638495 -0.18319492 -0.17000489 -0.15681485 -0.14362482 -0.13043478 -0.11651197 -0.10332193 -0.08939912 -0.07547631 -0.06082071 -0.04616512 -0.03150953 -0.01612115 0.0 0.005)
OUTPUT_FILE="TrueNeutrinoEnergy_Enubias_relative_20thaugbinning.yaml"
PARAM_LIST_FILE="TrueNeutrinoEnergy_Enubias_relative_parameter__20thaugbinning.txt"

# Start writing to the output files
echo "Systematics:" > $OUTPUT_FILE
echo "ParameterName, ene_min, ene_max, erec_minus_ene_min, erec_minus_ene_max" > $PARAM_LIST_FILE

# Generate systematic blocks using the non-uniform bin edges
COUNT=1
for ((i=0; i<${#ENE_BINS[@]}-1; i++)); do
    for ((j=0; j<${#EREC_MINUS_ENE_BINS[@]}-1; j++)); do
        ENE_LOW=${ENE_BINS[$i]}
        ENE_HIGH=${ENE_BINS[$((i+1))]}
        EREC_MINUS_ENE_LOW=${EREC_MINUS_ENE_BINS[$j]}
        EREC_MINUS_ENE_HIGH=${EREC_MINUS_ENE_BINS[$((j+1))]}

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
    - isRelativeEnubias:
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
