#!/bin/bash

# Loop through all files in ./outputs/ that match the pattern
for file in ./Outputs/projections_outputs/Projections*.root; do
    # Extract the base name (e.g., ProjectionsXYZ.root → XYZ)
    filename=$(basename "$file")             # ProjectionsXYZ.root
    suffix="${filename#Projections}"         # XYZ.root
    stem="${suffix%.root}"                   # XYZ

    # Construct input and output paths
    input="./Outputs/projections_outputs/Projections${stem}.root"
    output="Outputs/acceptance_plots/AcceptancePlots${stem}.pdf"

    # Run the ROOT macro
    root -l -q -b "Scripts/makeAcceptanceCorrectionPlots.C(\"$input\", \"$output\")"
done

