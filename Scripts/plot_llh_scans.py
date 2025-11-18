#!/usr/bin/python
#
# plot the results of "LLH_scan"
#

import os
import sys
import fnmatch
import ROOT
from ROOT import TFile, TIter, TKey, TCanvas, TGraph

if (len(sys.argv) != 3):
    print "Useage: " + str(sys.argv[0]) + " /path/to/llh_scan.root /path/to/output"
    quit(1)

### main 
data_loc = sys.argv[1]
output_loc = sys.argv[2] #+ "/"
out_name = "llh_scans"

print "-------------------------------------------"
ROOT.gROOT.SetBatch(True)

file1 = TFile(data_loc)
#file1.ls()
nextkey = TIter(file1.GetListOfKeys())

canv = TCanvas("llh_scans","") #,1366,768)
canv.Print(output_loc + out_name + ".ps[")

ppp = 6 # Plots Per Page
pc = 1  # Plot Counter
canv.Divide(2,3)

for key in nextkey:
    if (str(key.GetClassName()) == "TGraph"): 
        
        print key.GetName()
        llh_plot = key.ReadObj()
        canv.cd(pc)
        ROOT.gPad.SetGrid()

        if (fnmatch.fnmatch(key.GetName(), "*flux*")):
            llh_plot.SetMaximum(1000)
        else:
            llh_plot.SetMaximum(100)
    
        llh_plot.SetMinimum(93)
        llh_plot.SetTitle(key.GetName())
        llh_plot.SetLineColor(4)
        llh_plot.SetLineStyle(3)
        llh_plot.SetMarkerColor(4)
        llh_plot.SetMarkerStyle(7)
        llh_plot.SetMarkerSize(1)
        llh_plot.Draw("APC")
        
        pc = pc + 1
        if (pc > ppp):
            canv.Print(output_loc + out_name + ".ps")
            pc = 1

canv.Print(output_loc + out_name + ".ps]")

os.system("ps2pdf " + output_loc + out_name + ".ps " + output_loc + out_name + ".pdf")
print "finished!"        
