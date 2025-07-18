#!/bin/bash

CONFIG_TEMPLATE="Configs/EventRates_Beam_NDGAr_q0q3.yaml"
CONFIG_DIR="Configs/GeneratedConfigs"
OUTPUT_DIR="Outputs/projections_outputs/loop/"
mkdir -p "$CONFIG_DIR"

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
  ENUMIN=$(echo $range | awk '{print $1}')
  ENUMAX=$(echo $range | awk '{print $2}')

  TAG="Enu_${ENUMIN}to${ENUMAX}"
  CONFIG="$CONFIG_DIR/EventRates_${TAG}.yaml"
  OUTPUT_FILE="${OUTPUT_DIR}/Projections_CC_${TAG}.root"
  LOG_FILE="output_${TAG}.txt"

  # Copy the base config and modify it
  cp "$CONFIG_TEMPLATE" "$CONFIG"

  # Replace the Enu range
  sed -i -E "s/(Name: \"Enu\".*\n.*Range: )\[[^]]+\]/\1[${ENUMIN},${ENUMAX}]/" "$CONFIG"

  # Replace the OutputFile line
  sed -i -E "s|(OutputFile: ).*|\1\"${OUTPUT_FILE}\"|" "$CONFIG"

  echo "Launching: $CONFIG -> $OUTPUT_FILE"
  nohup Projections "$CONFIG" > "$LOG_FILE" 2>&1 &
done

