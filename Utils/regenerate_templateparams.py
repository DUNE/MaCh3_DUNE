#!/usr/bin/env python3
"""
regenerate_templateparams.py

Regenerates a MaCh3 templateparams.yaml-style systematics file, keeping the
existing TrueNeutrinoEnergy slicing (e0..eN) exactly as-is, but replacing the
Enubias bin edges (eb0..ebM) in every energy slice with a new set -- e.g. the
quantile-based edges computed from the MC distribution.

Usage:
    python regenerate_templateparams.py \
        --old templateparams.yaml \
        --new-enubias-edges enubias_quantile_edges.txt \
        --out templateparams_quantile.yaml

The --new-enubias-edges file can be either:
  (a) a plain comma/space-separated list of N+1 edge values on one line, or
  (b) the "# Plain edge list: # [...]" line as printed directly by the
      mcmc_diagnostics tool's WriteQuantileEdgesToFile output -- this script
      will find and parse that line automatically if present.

Everything else in each systematic block (Error, FlatPrior, ParameterBounds,
ParameterValues, StepScale, Type, ParameterGroup, SampleNames) is copied
verbatim from the corresponding OLD block for that energy slice, so only the
Enubias bounds and the eb-index (and therefore FancyName) change.
"""

import argparse
import re
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("This script requires PyYAML. Install with: pip install pyyaml")


def parse_edges_from_file(path: str):
    """Parse either a plain edge list line, or a bare comma/space separated
    list of numbers, from the given file."""
    text = Path(path).read_text()

    # Look for a "# Plain edge list:\n# [...]" style block first.
    match = re.search(r"#\s*\[([^\]]+)\]", text)
    if match:
        raw = match.group(1)
    else:
        # Fall back: treat the whole file as a list of numbers.
        raw = text

    tokens = re.split(r"[,\s]+", raw.strip())
    edges = [float(t) for t in tokens if t not in ("", "#")]
    if len(edges) < 2:
        sys.exit(f"Could not parse at least 2 edge values from {path}")
    return edges


def load_old_yaml(path: str):
    with open(path) as f:
        return yaml.safe_load(f)


def group_by_energy_slice(systematics):
    """
    Group the flat list of systematic blocks by their TrueNeutrinoEnergy
    bounds, preserving the order energy slices first appear in the file.
    Returns: list of (e_lo, e_hi, representative_block) tuples, in order.
    The representative_block is the FIRST eb-block found for that energy
    slice, used as a template for all fields except Enubias bounds/FancyName.
    """
    seen = {}
    order = []
    for entry in systematics:
        sys_block = entry["Systematic"]
        cuts = sys_block["KinematicCuts"]

        e_bounds = None
        for cut in cuts:
            if "TrueNeutrinoEnergy" in cut:
                e_bounds = tuple(cut["TrueNeutrinoEnergy"])
                break
        if e_bounds is None:
            continue  # skip anything without a TrueNeutrinoEnergy cut

        if e_bounds not in seen:
            seen[e_bounds] = sys_block
            order.append(e_bounds)

    return [(e_lo, e_hi, seen[(e_lo, e_hi)]) for (e_lo, e_hi) in order]


def build_new_systematics(energy_slices, new_enubias_edges):
    """
    Build the new flat list of systematic blocks: for every energy slice
    (in the order given) and every new Enubias bin, emit one block whose
    non-Enubias fields are copied from that slice's representative block.
    """
    new_systs = []
    n_e = len(energy_slices)
    n_eb = len(new_enubias_edges) - 1

    for e_idx, (e_lo, e_hi, template_block) in enumerate(energy_slices):
        for eb_idx in range(n_eb):
            eb_lo = new_enubias_edges[eb_idx]
            eb_hi = new_enubias_edges[eb_idx + 1]

            block = {
                "SampleNames": template_block.get("SampleNames", ["*"]),
                "Error": template_block.get("Error", 0.5),
                "FlatPrior": template_block.get("FlatPrior", False),
                "KinematicCuts": [
                    {"TrueNeutrinoEnergy": [e_lo, e_hi]},
                    {"Enubias": [eb_lo, eb_hi]},
                ],
                "Names": {"FancyName": f"tmplw_e{e_idx}_eb{eb_idx}"},
                "ParameterBounds": template_block.get("ParameterBounds", [0.0, 4.0]),
                "ParameterValues": template_block.get(
                    "ParameterValues", {"Generated": 1.0, "PreFitValue": 1.0}
                ),
                "StepScale": template_block.get("StepScale", {"MCMC": 1.0}),
                "Type": template_block.get("Type", "Norm"),
                "ParameterGroup": template_block.get("ParameterGroup", "Xsec"),
            }
            new_systs.append({"Systematic": block})

    return new_systs


# ---------------------------------------------------------------------------
# Custom YAML dumper: match the original file's block style, list-item order,
# and indentation as closely as possible.
# ---------------------------------------------------------------------------
class OrderedDumper(yaml.SafeDumper):
    def ignore_aliases(self, data):
        # Force fully explicit output (no &id.../*id... anchors), since the
        # original templateparams.yaml has no anchors and we don't want to
        # rely on yaml-cpp correctly resolving PyYAML-style aliases.
        return True


def represent_dict_in_order(dumper, data):
    return dumper.represent_mapping("tag:yaml.org,2002:map", data.items())


OrderedDumper.add_representer(dict, represent_dict_in_order)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old", required=True, help="Existing templateparams.yaml")
    ap.add_argument("--new-enubias-edges", required=True,
                     help="File containing the new Enubias bin edges "
                          "(plain list, or the 'Plain edge list' line from "
                          "WriteQuantileEdgesToFile output)")
    ap.add_argument("--out", required=True, help="Output path for the new YAML")
    args = ap.parse_args()

    print(f"[Info] Loading existing YAML: {args.old}")
    old_data = load_old_yaml(args.old)
    old_systs = old_data.get("Systematics", [])
    if not old_systs:
        sys.exit("No 'Systematics' key found in the old YAML -- wrong file?")

    print(f"[Info] Parsing new Enubias edges from: {args.new_enubias_edges}")
    new_edges = parse_edges_from_file(args.new_enubias_edges)
    print(f"[Info] Parsed {len(new_edges) - 1} new Enubias bins: {new_edges}")

    energy_slices = group_by_energy_slice(old_systs)
    print(f"[Info] Found {len(energy_slices)} TrueNeutrinoEnergy slices in the old YAML.")

    new_systs = build_new_systematics(energy_slices, new_edges)
    print(f"[Info] Built {len(new_systs)} new systematic blocks "
          f"({len(energy_slices)} energy slices x {len(new_edges) - 1} Enubias bins).")

    out_data = {"Systematics": new_systs}

    with open(args.out, "w") as f:
        f.write("---\n")
        yaml.dump(out_data, f, Dumper=OrderedDumper,
                   default_flow_style=False, sort_keys=False)

    print(f"[Info] Wrote new YAML to: {args.out}")
    print("[Info] Done. Sanity-check a few blocks by eye before using this in a fit run.")


if __name__ == "__main__":
    main()