#!/usr/bin/env python
"""
Write the (TrueE x RecoE) template-weight systematics YAML, skipping bins
with no events in any sample.

Usage:
  python make_xsec_bin_systs.py binmask_bins.csv > xsec_bins.yaml
  python make_xsec_bin_systs.py binmask_bins.csv --threshold 1 > xsec_bins.yaml
  python make_xsec_bin_systs.py --all > xsec_bins_full.yaml    # old behaviour

binmask_bins.csv comes from find_empty_bins.py. A bin is kept if its summed
event rate is > --threshold (default 0, i.e. any events at all).

Parameter names keep their ORIGINAL indices (tmplw_e{ie}_eb{ieb}), so a given
name always refers to the same bin whatever is removed.
"""
import argparse
import csv
import sys

enu_bins = [0, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.575, 1.65, 1.725, 1.8, 1.875, 1.95, 2.025, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.225, 3.35, 3.475, 3.6, 3.75, 3.9, 4.075, 4.25, 4.425, 4.6, 4.8, 5, 5.5, 6, 7, 8, 10]       # TrueNeutrinoEnergy
enubias_bins = [0, 0.2, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.575, 1.65, 1.725, 1.8, 1.875, 1.95, 2.025, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3, 3.1, 3.225, 3.35, 3.475, 3.6, 3.75, 3.9, 4.075, 4.25, 4.425, 4.6, 4.8, 5, 5.5, 6, 7, 8, 10]  # RecoNeutrinoEnergy

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("bincsv", nargs="?", help="binmask_bins.csv from find_empty_bins.py")
ap.add_argument("--threshold", type=float, default=0.0,
                help="keep a bin only if its summed events are > this (default 0)")
ap.add_argument("--all", action="store_true", help="write every bin (no filtering)")
a = ap.parse_args()
if not a.all and not a.bincsv:
    ap.error("give binmask_bins.csv, or --all")

# Event rate per (true bin, reco bin), checked against the edges above so a
# CSV made with a different binning can't silently select the wrong bins.
rate = {}
if not a.all:
    close = lambda x, y: abs(float(x) - y) < 1e-6
    with open(a.bincsv) as f:
        for r in csv.DictReader(f):
            ie, ieb = int(r["true_bin"]), int(r["reco_bin"])
            if not (close(r["true_lo"], enu_bins[ie]) and close(r["true_hi"], enu_bins[ie + 1])
                    and close(r["reco_lo"], enubias_bins[ieb]) and close(r["reco_hi"], enubias_bins[ieb + 1])):
                sys.exit(f"Bin edges in {a.bincsv} don't match this script "
                         f"(true bin {ie}, reco bin {ieb}). Was it made with the same binning?")
            rate[(ie, ieb)] = float(r["total_events"])
    expected = (len(enu_bins) - 1) * (len(enubias_bins) - 1)
    if len(rate) != expected:
        sys.exit(f"{a.bincsv} has {len(rate)} bins, expected {expected}.")

print("""---
Systematics:""")

nkept = 0
for ie in range(len(enu_bins) - 1):
    for ieb in range(len(enubias_bins) - 1):
        if not a.all and rate[(ie, ieb)] <= a.threshold:
            continue
        nkept += 1
        pnm = f"tmplw_e{ie}_eb{ieb}"

        print(f"""- Systematic:
    SampleNames: ["*"]
    Error:  0.5
    FlatPrior: false
    KinematicCuts:
    - TrueNeutrinoEnergy:
      - {enu_bins[ie]}
      - {enu_bins[ie+1]}
    - RecoNeutrinoEnergy:
      - {enubias_bins[ieb]}
      - {enubias_bins[ieb+1]}
    Names:
      FancyName: {pnm}
    ParameterBounds:
    - 0.0
    - 4.0
    ParameterValues:
      Generated: 1.0
      PreFitValue: 1.0
    StepScale:
      MCMC: 1.0
    Type: Norm
    ParameterGroup: Xsec""")

total = (len(enu_bins) - 1) * (len(enubias_bins) - 1)
print(f"Wrote {nkept} of {total} parameters", file=sys.stderr)