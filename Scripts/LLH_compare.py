#!/usr/bin/env python3

import sys
import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

if len(sys.argv) != 4:
    print("Usage: compare_llh_scans.py file1.root file2.root output.pdf")
    sys.exit(1)

f1 = ROOT.TFile.Open(sys.argv[1])
f2 = ROOT.TFile.Open(sys.argv[2])
out = sys.argv[3]

# Oscillation parameters to keep
keep_osc_param = [
    "sin2th_12",
    "sin2th_23",
    "sin2th_13",
    "delm2_12",
    "delm2_23",
    "delta_cp"
]

# Sample types to keep
keep_samples = [
    "fc",
    "pc"
]

def is_osc_param(name):
    """
    Select only oscillation parameter histograms
    containing FC or PC.
    """
    name_low = name.lower()

    has_param = any(p in name_low for p in keep_osc_param)
    has_sample = any(s in name_low for s in keep_samples)

    return has_param and has_sample


def get_hists(directory):
    hists = {}

    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        name = key.GetName()

        if obj.InheritsFrom("TH1"):
            if is_osc_param(name):
                hists[name] = obj.Clone(name + "_clone")
                hists[name].SetDirectory(0)

        elif obj.InheritsFrom("TDirectory"):
            hists.update(get_hists(obj))

    return hists


h1s = get_hists(f1)
h2s = get_hists(f2)

print(f"Found {len(h1s)} matching histograms in file1")
print(f"Found {len(h2s)} matching histograms in file2")

c = ROOT.TCanvas("c", "LLH comparison", 1200, 800)

c.Print(out + "[")

for name in sorted(h1s):

    if name not in h2s:
        print(f"Skipping {name} (not found in file2)")
        continue

    h1 = h1s[name]
    h2 = h2s[name]

    ymin = min(h1.GetMinimum(), h2.GetMinimum())
    ymax = max(h1.GetMaximum(), h2.GetMaximum())

    margin = 0.05 * (ymax - ymin + 1e-9)

    h1.SetLineColor(ROOT.kRed)
    h1.SetLineWidth(2)

    h2.SetLineColor(ROOT.kBlue)
    h2.SetLineWidth(2)

    h1.SetMinimum(0)
    h1.SetMaximum(ymax + margin)

    h2.SetMinimum(0)
    h2.SetMaximum(ymax + margin)

    c.Clear()

    h1.SetTitle(f"Sample LLH;{name};-2(ln L_{{sample}})")

    h1.Draw("HIST")
    h2.Draw("HIST SAME")

    leg = ROOT.TLegend(0.15, 0.75, 0.40, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)

    leg.AddEntry(h1, "2D Binning", "l")
    leg.AddEntry(h2, "3D Binning", "l")

    leg.Draw()

    c.Print(out)

c.Print(out + "]")

print("Saved:", out)

f1.Close()
f2.Close()
