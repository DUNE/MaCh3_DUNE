#!/bin/bash

CONFIG_TEMPLATE="Configs/EventRates_Beam_NDGAr_q0q3.yaml"
CONFIG_DIR="Configs/GeneratedConfigs"
OUTPUT_DIR="Outputs/projections_outputs/loop"
mkdir -p "$CONFIG_DIR"
mkdir -p "$OUTPUT_DIR"

# List of Enu ranges
declare -a RANGES=(
  "0.00 1.50"
  "1.50 1.75"
  "1.75 2.00"
  "2.00 2.25"
  "2.25 2.50"
  "2.50 2.75"
  "2.75 3.00"
  "3.00 3.25"
  "3.25 3.50"
  "3.50 4.00"
  "4.00 4.50"
  "4.50 5.00"
)

for range in "${RANGES[@]}"; do
  ENUMIN=$(echo "$range" | awk '{printf "%.2f", $1}')
  ENUMAX=$(echo "$range" | awk '{printf "%.2f", $2}')

  TAG="Enu_${ENUMIN}to${ENUMAX}"
  CONFIG="$CONFIG_DIR/EventRates_${TAG}.yaml"
  OUTPUT_FILE="${OUTPUT_DIR}/Projections_CC_${TAG}.root"
  LOG_FILE="output_${TAG}.txt"

  # Copy the template
  cp "$CONFIG_TEMPLATE" "$CONFIG"

  # Update OutputFile in 'General' section
  yq eval ".General.OutputFile = \"${OUTPUT_FILE}\"" -i "$CONFIG"

  # Find index of the Enu cut in the GeneralKinematicCuts list
  # ENU_INDEX=$(yq eval '.GeneralKinematicCuts | map(.Name == "Enu") | index(true)' "$CONFIG")
  ENU_INDEX=$(yq eval '.GeneralKinematicCuts | to_entries | map(select(.value.Name == "Enu")) | .[0].key' "$CONFIG")

  if [[ "$ENU_INDEX" == "null" ]]; then
    echo "Error: Could not find 'Enu' cut in $CONFIG"
    continue
  fi

  # Update Enu Range
  yq eval ".GeneralKinematicCuts[$ENU_INDEX].Range = [${ENUMIN}, ${ENUMAX}]" -i "$CONFIG"

  echo "Launching: $CONFIG -> $OUTPUT_FILE"
  nohup Projections "$CONFIG" > "$LOG_FILE" 2>&1 &
done

echo "All jobs launched."

