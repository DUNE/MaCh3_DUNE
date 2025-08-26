#!/bin/bash

# User-defined variables
Q0_MIN=0.0
Q0_MAX=5.0
Q3_MIN=0.0
Q3_MAX=5.0
Q0_BINS=100
Q3_BINS=100
OUTPUT_FILE="q0q3_0.0_5.0GeV_triangle_10e4bins.yaml"
PARAM_LIST_FILE="q0q3_0.0_5.0GeV_triangle_parameter_list_10e4bins.txt"

# Calculate bin widths
Q0_STEP=$(echo "($Q0_MAX - $Q0_MIN) / $Q0_BINS" | bc -l)
Q3_STEP=$(echo "($Q3_MAX - $Q3_MIN) / $Q3_BINS" | bc -l)

# Start writing to the output file
echo "Systematics:" > $OUTPUT_FILE
echo "ParameterName, q0_min, q0_max, q3_min, q3_max" > $PARAM_LIST_FILE

# Generate systematic blocks
COUNT=1
for ((i=0; i<$Q0_BINS; i++)); do
    for ((j=0; j<$Q3_BINS; j++)); do
        Q0_LOW=$(echo "$Q0_MIN + $i * $Q0_STEP" | bc -l)
        Q0_HIGH=$(echo "$Q0_LOW + $Q0_STEP" | bc -l)
        Q3_LOW=$(echo "$Q3_MIN + $j * $Q3_STEP" | bc -l)
        Q3_HIGH=$(echo "$Q3_LOW + $Q3_STEP" | bc -l)

        # Only include bins where q3 >= q0
        COMPARE=$(echo "$Q3_LOW >= $Q0_LOW" | bc)
        if [ "$COMPARE" -eq 1 ]; then
            PARAM_NAME="q0q3_$COUNT"
            
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
