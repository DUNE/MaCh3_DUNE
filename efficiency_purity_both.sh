#!/bin/bash

# Input YAML files

EVENTRATES_CONFIG="/scratch/abipeake/MaCh3_DUNE_FD/MaCh3_DUNE/configs/EventRates_BeamVD.yaml"

INPUT_FILE="/scratch/abipeake/MaCh3_DUNE_FD/MaCh3_DUNE/configs/Samples/FD_VD/FD_VD_FHC.yaml"
OUTPUT_PREFIX="modifiedcuts_"
EVENTRATES_CONFIG="/scratch/abipeake/MaCh3_DUNE_FD/MaCh3_DUNE/configs/EventRates_BeamVD.yaml"

ROOT_FILE="efficiency_purity_results_nuenumucut_2flavour.root"
ROOT_SCRIPT="store_results.C"
CSV_FILE="efficiency_purity_results_nuenumucut_2flavour.csv"

# Initialize the ROOT file
echo "void store_results() { TFile f(\"$ROOT_FILE\", \"RECREATE\"); TTree t(\"EfficiencyPurity\", \"Efficiency and Purity Data\"); }" > "$ROOT_SCRIPT"
root -l -b -q "$ROOT_SCRIPT"

# Loop through values from 0.4 to 0.7 in increments of 0.1
for nu_cut in $(seq 0.1 0.1 1.0); do
    OUTPUT_FILE="configs/Samples/test/${OUTPUT_PREFIX}${nu_cut}_nue.yaml"

    cp $INPUT_FILE $OUTPUT_FILE

    # Modify YAML file with new nue_cut value
    sed -i "s|^nue_cut.*|nue_cut: $nu_cut|g" $OUTPUT_FILE
    sed -i "s|^numu_cut.*|numu_cut: $nu_cut|g" $OUTPUT_FILE
    #sed "s/^ nue_cut: .*/nue_cut: $nue_cut/" "$INPUT_FILE" > "$OUTPUT_FILE"
    echo "Created: $OUTPUT_FILE with nu_cut = $nu_cut"

    # Update EventRates config to use the modified YAML file
    #sed "s|DUNESamples: \[.*\]|DUNESamples: [\"$OUTPUT_FILE\"]|" "$EVENTRATES_CONFIG" > "updated_eventrates.yaml"
    #sed -i "s|DUNESamples.*|DUNESamples: ["$OUTPUT_FILE"]|g"  "$EVENTRATES_CONFIG"
    sed -i "s|DUNESamples.*|DUNESamples: [\"$OUTPUT_FILE\"]|g" "$EVENTRATES_CONFIG"
    echo "Updated EventRates config with $OUTPUT_FILE"

    grep "DUNESamples" "$EVENTRATES_CONFIG"


    # Debugging: Check nue_cut value in modified YAML
    grep "nue_cut" "$OUTPUT_FILE"
    echo "DEBUG: EVENTRATES_CONFIG='$EVENTRATES_CONFIG'"

    # Run EventRates and check for errors
    #OUTPUT=$(./build/src/EventRates updated_eventrates.yaml 2>&1)
    OUTPUT=$(./build/src/EventRates "$EVENTRATES_CONFIG" 2>&1)

    if [[ $? -ne 0 ]]; then
        echo "ERROR: EventRates execution failed for nue_cut=$nue_cut"
        continue
    fi

    echo "DEBUG: Full EventRates Output"
    echo "$OUTPUT"

    # Extract relevant values
    eventsthatpasscut=$(echo "$OUTPUT" | grep "eventsthatpassedcut =" | awk -F '= ' '{print $2}' | tail -n1)
    total_true_ccnue=$(echo "$OUTPUT" | grep "total_true_ccnue  =" | awk -F '= ' '{print $2}' | tail -n1)
    total_events_incut=$(echo "$OUTPUT" | grep "total_events_incut =" | awk -F '= ' '{print $2}' | tail -n1)

    echo "Extracted values: eventsthatpasscut=$eventsthatpasscut, total_true_ccnue=$total_true_ccnue, total_events_incut=$total_events_incut"

    
done