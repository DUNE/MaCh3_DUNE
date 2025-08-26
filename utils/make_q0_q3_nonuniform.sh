#!/bin/bash

# Define bin edges
Q0_BINS=(0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 \
         0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 \
         1.00 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.45 \
         1.50 1.55 1.60 1.65 1.70 1.75 1.80 1.85 1.90 1.95 \
         2.00 2.05 2.10 2.15 2.20 2.25 2.30 2.35 2.40 2.45 \
         2.50 2.55 2.60 2.65 2.70 2.75 2.80 2.85 2.90 2.95 \
         3.00 3.05 3.10 3.15 3.20 3.25 3.30 3.35 3.40 3.45 \
         3.50 3.55 3.60 3.65 3.70 3.75 3.80 3.85 3.90 4.40 \
         4.90 5.00)

Q3_BINS=(0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 \
         0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 \
         1.00 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.45 \
         1.50 1.55 1.60 1.65 1.70 1.75 1.80 1.85 1.90 1.95 \
         2.00 2.05 2.10 2.15 2.20 2.25 2.30 2.35 2.40 2.45 \
         2.50 2.55 2.60 2.65 2.70 2.75 2.80 2.85 2.90 2.95 \
         3.00 3.05 3.10 3.15 3.20 3.25 3.30 3.35 3.40 3.45 \
         3.50 3.55 3.60 3.65 3.70 3.75 3.80 3.85 3.90 4.40 \
         4.90 5.00)

OUTPUT_FILE="q0_q3_systematics.yaml"
PARAM_LIST_FILE="q0_q3_param_list.txt"

# Start writing to the output files
echo "Systematics:" > $OUTPUT_FILE
echo "ParameterName, q0_min, q0_max, q3_min, q3_max" > $PARAM_LIST_FILE

COUNT=1
for ((i=0; i<${#Q0_BINS[@]}-1; i++)); do
    for ((j=0; j<${#Q3_BINS[@]}-1; j++)); do
        
        Q0_LOW=${Q0_BINS[$i]}
        Q0_HIGH=${Q0_BINS[$((i+1))]}
        Q3_LOW=${Q3_BINS[$j]}
        Q3_HIGH=${Q3_BINS[$((j+1))]}
        
        # Apply the triangular mask
        if (( $(echo "$Q0_LOW <= $Q3_LOW" | bc -l) )); then
            
            PARAM_NAME="q0_q3_bin_$COUNT"
            
            cat <<EOL >> $OUTPUT_FILE
- Systematic:
    Sample_Name: ["ND"]
    Error: 0.5
    FlatPrior: false
    KinematicCuts:
    - q0:
      - $Q0_LOW
      - $Q0_HIGH
    - q3:
      - $Q3_LOW
      - $Q3_HIGH
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
            
            echo "$PARAM_NAME, $Q0_LOW, $Q0_HIGH, $Q3_LOW, $Q3_HIGH" >> $PARAM_LIST_FILE
            COUNT=$((COUNT+1))
        fi
    done
done

echo "Config file generated: $OUTPUT_FILE"
echo "Parameter list generated: $PARAM_LIST_FILE"
