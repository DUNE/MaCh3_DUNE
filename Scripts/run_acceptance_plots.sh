#!/bin/bash

# Loop through all files in ./outputs/ that match the pattern
for file in ./outputs/Projections*.root; do
    # Extract the base name (e.g., ProjectionsXYZ.root → XYZ)
    filename=$(basename "$file")             # ProjectionsXYZ.root
    suffix="${filename#Projections}"         # XYZ.root
    stem="${suffix%.root}"                   # XYZ

    # Construct input and output paths
    input="./outputs/Projections${stem}.root"
    output="outputs/AcceptancePlots${stem}.pdf"

    # Run the ROOT macro
    root -l -q -b "scripts/makeAcceptanceCorrectionPlots.C(\"$input\", \"$output\")"
done

