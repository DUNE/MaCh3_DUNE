#!/usr/bin/env bash
#
# run_all_splines.sh
#
# Batch-run plot_splines.py over every ROOT file in
# Analysis_splines/OnplusOffAxis/, using one shared YAML config,
# and save each file's plots into its own subfolder under an
# output directory (named after the off-axis position extracted
# from the filename).
#
# Usage
# -----
#   ./run_all_splines.sh /path/to/config.yaml [input_dir] [output_dir] [extra plot_splines.py args...]
#
# Examples
# --------
#   ./run_all_splines.sh config.yaml
#   ./run_all_splines.sh config.yaml Analysis_splines/OnplusOffAxis plots
#   ./run_all_splines.sh config.yaml Analysis_splines/OnplusOffAxis plots --format png --group Xsec
#
# Requires plot_splines.py to be in the same directory as this script,
# or set PLOT_SCRIPT below to point at it.

set -euo pipefail

YAML="${1:?Usage: $0 <config.yaml> [input_dir] [output_dir] [extra args...]}"
shift

INPUT_DIR="${1:-Analysis_splines/OnplusOffAxis}"
[[ $# -ge 1 ]] && shift

OUTPUT_DIR="${1:-plots}"
[[ $# -ge 1 ]] && shift

# Anything left on the command line is passed straight through to
# plot_splines.py (e.g. --format png --group Xsec --cmap viridis)
EXTRA_ARGS=("$@")

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLOT_SCRIPT="${PLOT_SCRIPT:-$SCRIPT_DIR/plot_splines.py}"

if [[ ! -f "$YAML" ]]; then
    echo "ERROR: YAML config not found: $YAML" >&2
    exit 1
fi
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: input directory not found: $INPUT_DIR" >&2
    exit 1
fi
if [[ ! -f "$PLOT_SCRIPT" ]]; then
    echo "ERROR: plot_splines.py not found at $PLOT_SCRIPT (set PLOT_SCRIPT env var to override)" >&2
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

shopt -s nullglob
root_files=("$INPUT_DIR"/*.root)
shopt -u nullglob

if [[ ${#root_files[@]} -eq 0 ]]; then
    echo "ERROR: no .root files found in $INPUT_DIR" >&2
    exit 1
fi

echo "YAML config : $YAML"
echo "Input dir   : $INPUT_DIR"
echo "Output dir  : $OUTPUT_DIR"
echo "Found ${#root_files[@]} ROOT file(s)."
echo "Extra args  : ${EXTRA_ARGS[*]:-<none>}"
echo

n_ok=0
n_fail=0
fail_list=()

for root_file in "${root_files[@]}"; do
    base="$(basename "$root_file" .root)"
    # e.g. OffAxis13.75m_numuCC_numode_splines -> OffAxis13.75m
    label="${base%%_*}"

    subdir="$OUTPUT_DIR/$label"
    mkdir -p "$subdir"

    echo "=================================================================="
    echo "Processing: $root_file"
    echo "  -> label : $label"
    echo "  -> outdir: $subdir"
    echo "=================================================================="

    if python3 "$PLOT_SCRIPT" \
        --yaml "$YAML" \
        --root "$root_file" \
        --outdir "$subdir" \
        "${EXTRA_ARGS[@]}"; then
        n_ok=$((n_ok + 1))
    else
        n_fail=$((n_fail + 1))
        fail_list+=("$root_file")
    fi
    echo
done

echo "=================================================================="
echo "Batch complete."
echo "  Succeeded: $n_ok"
echo "  Failed:    $n_fail"
if [[ $n_fail -gt 0 ]]; then
    echo "  Failed files:"
    for f in "${fail_list[@]}"; do
        echo "    - $f"
    done
fi
echo "  All outputs under: $OUTPUT_DIR/<label>/"