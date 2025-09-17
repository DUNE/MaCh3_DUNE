#!/bin/bash

bash Scripts/q0q3loop/loop_projections.sh

echo "Waiting for all background jobs to finish..."

while pgrep -f "Projections" > /dev/null; do
  echo "Jobs still running... sleeping for 10s"
  sleep 10
done

echo "All projections completed."

echo "Running ROOT macro for roughness analysis..."
root -l -b -q 'Scripts/q0q3loop/test_acceptance_metrics.C'

echo "All done!"

