#!/usr/bin/env bash

set -x

./build/bin/$EventRates CIValidations/CIInputs/CIFitConfig.yaml | grep '^\['"$EventRates"'\.cpp\]' | tee CIValidations/CIOutputs/BeamEventRates_New.txt

diff CIValidations/CIOutputs/BeamEventRates.txt CIValidations/CIOutputs/BeamEventRates_New.txt > diff_output.txt
DIFF_EXIT=$?
cat diff_output.txt

if [ $DIFF_EXIT -gt 0 ]; then
  echo "Differences found:"
  exit 1
else
  echo "Workflow run successfully, no differences found!"
fi
