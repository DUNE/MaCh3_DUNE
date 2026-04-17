#!/usr/bin/env python3

import sys
import ROOT

ROOT.gStyle.SetOptStat(0)

if len(sys.argv) != 4:
    print("Usage: compare_llh_scans.py file1.root file2.root output.pdf")
    sys.exit(1)

f1 = ROOT.TFile.Open(sys.argv[1])
f2 = ROOT.TFile.Open(sys.argv[2])
out = sys.argv[3]

ROOT.gROOT.SetBatch(True)

keep_osc_param = [
    "sin2th_12_sam",
    "sin2th_23_sam",
    "sin2th_13_sam",
    "delm2_12_sam",
    "delm2_23_sam",
    "delta_cp_sam"
]

def is_osc_param(name):
    name_low = name.lower()
    return any(k in name_low for k in keep_osc_param)

def get_hists(directory):
    hists = {}

    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        name =  key.GetName()

        if obj.InheritsFrom("TH1"):
            if is_osc_param(name):
                hists[name] = obj.Clone(name + "_clone")

        elif obj.InheritsFrom("TDirectory"):
            hists.update(get_hists(obj))

    return hists


h1s = get_hists(f1)
h2s = get_hists(f2)

c = ROOT.TCanvas("c", "LLH comparison", 1200, 800)
c.Print(out + "[")

for name in h1s:
    if name not in h2s:
        continue

    h1 = h1s[name]
    h2 = h2s[name]

    ymin = min(h1.GetMinimum(), h2.GetMinimum())
    ymax = max(h1.GetMaximum(), h2.GetMaximum())

    margin = 0.05 * (ymax - ymin + 1e-9)

    # plot style
    h1.SetLineColor(ROOT.kRed)
    h1.SetLineWidth(2)

    h2.SetLineColor(ROOT.kBlue)
    h2.SetLineWidth(2)

   # h1.SetMinimum(ymin - margin)
    h1.SetMinimum(0)
    h1.SetMaximum(ymax + margin)
    h2.SetMinimum(0)
   # h2.SetMinimum(ymin - margin)
    h2.SetMaximum(ymax + margin)

    c.Clear()

    h1.Draw("HIST")
    h2.Draw("HIST SAME")
    
    h1.SetTitle("Sample LLH"+ ";" + name+ ";-2(ln L_{sample})")
    

    leg = ROOT.TLegend(0.15, 0.75, 0.35, 0.89)
    leg.AddEntry(h1, "2D Binning", "l")
    leg.AddEntry(h2, "3D Binning", "l")
    leg.Draw()

    c.Print(out)

c.Print(out + "]")

print("Saved:", out)
