#!/usr/bin/env python

import sys
import math
import os
import logging
import getopt
import fnmatch

#sys.argv.append( '-b-' )
import ROOT
from ROOT import gROOT, TCanvas, TPad, TH1D, TH2D, TStyle, TFile, TApplication, TString, TChain, gSystem, gDirectory, TGraph, TMultiGraph, TH1F, TF1, AddressOf, TPaveText

from random import choice
from array import array

if len(sys.argv) != 2:
        print "Usage: " + sys.argv[0] + " /path/of/files/"
        sys.exit(1)

def histGausFit(gainDist):
	func = TF1("gaus","gaus",0,1)
	gainDist.Fit("gaus","RQ")
	const = func.GetParameter(0)
	mean  = func.GetParameter(1)
	sigma = func.GetParameter(2)
	return const, mean, sigma

print "summing files from " + sys.argv[1]
#hpost = TH1D("hpost","hpost",100,0.5,1)

# config ################
t_off = 1.7 # y axis title offset
webloc = "." # CHANGE THIS TO WHEREVER YOU WANT THE PLOTS TO BE SAVED
# maybe add it as an argument
nburn = 1#0000 # adjust burn in steps
binning = 100
nflux = 49
ndet = 19
nxsec = 23
nNDdet = 39
nDiv_delm23 = 5 # number of divisions on xaxis for delm23
dud_check = True #False # check for files whose osc pars dont change

# histograms ############

h_flux = []
h_xsec = []
h_skdet = []
h_nddet = []

#########################

totalfile = sys.argv[1]
#outname = sys.argv[2]

#print totalfile

#tfile = TFile(totalfile);

#hpost = TH1D("hpost","hpost",100,0.5,1)

code = "struct mcmc_obj {"
code += "double sin2th_23;"
code += "double sin2th_13;"
code += "double delta_cp;"
code += "double delm2_12;"
code += "double delm2_23;"
code += "double accProb;"
code += "double LogL;"
code += "};"

#sin2th_12       = 0
# sin2th_23       = 0
# sin2th_13       = 0
# delm2_12        = 0
# delm2_23        = 0
# delta_cp

ROOT.gROOT.SetBatch(True)
gROOT.ProcessLine(code)

from ROOT import mcmc_obj
mcmc_tree = mcmc_obj()

# osc params #################

sin12_u = 0.5
sin12_l = 0.2
h_sin12 = TH1D("sin12","",binning,sin12_l, sin12_u)

sin23_u = 0.7
sin23_l = 0.3
h_sin23 = TH1D("sin23","",binning, sin23_l, sin23_u)

sin13_u = 0.25
sin13_l = 0.0
h_sin13 = TH1D("sin13","",binning, sin13_l, sin13_u)
h_sin13.SetFillStyle(3003)
h_sin13.SetFillColor(4)
h_sin13.SetLineColor(4)
h_sin13.SetLineStyle(1)

delm12_u = 8.5E-5 #7.50005E-5
delm12_l = 6.25E-5 #7.49995E-5
h_delm12 = TH1D("delm12","",binning, delm12_l, delm12_u)

delm23_u = 2.85E-3 #2.392E-3
delm23_l = 2E-3 #2.388E-3
h_delm23 = TH1D("delm23","",binning, delm23_l, delm23_u)

deltacp_u = 3.14
deltacp_l = -3.14
h_deltacp = TH1D("deltacp","",binning, deltacp_l, deltacp_u)

h_13_cp = TH2D("13_cp","",binning, sin13_l, sin13_u, binning, deltacp_l, deltacp_u)
h_23_cp = TH2D("23_cp","",binning, sin23_l, sin23_u, binning, deltacp_l, deltacp_u)
h_d23_23 = TH2D("d23_23","",binning, delm23_l, delm23_u, binning, sin23_l, sin23_u)
h_d23_13 = TH2D("d23_23","",binning, delm23_l, delm23_u, binning, sin13_l, sin13_u)
h_23_13 = TH2D("23_13","",binning, sin23_l, sin23_u, binning, sin13_l, sin13_u)
h_d23_cp = TH2D("d23_cp","",binning, delm23_l, delm23_u, binning, deltacp_l, deltacp_u)

#misc
h_accProb = TH1D("accProb","",binning, 0, 1.1)
h_LogL = TH1D("LogL","",binning, 100, 200)
h_13_LogL = TH2D("13logl","",binning, 0, 0.25, binning, 150, 160)
h_23_LogL = TH2D("23logl","",binning, 0.9, 1, binning, 150, 160)

h_projection = TH1D("proj_13","",binning/2,0, 0.25)
llh_13 = TH1D("llh13","",10000,0,10000)

# bad file list
bad_files = []

# file loop
dir_list = os.listdir(totalfile)

check = 0

for ttfile in dir_list:
	if (fnmatch.fnmatch(ttfile,"*.root") is not True):
		continue

	statinfo = os.stat(totalfile + "/" + ttfile)
	if (statinfo.st_size < 2000):
		continue

	tfile = TFile(totalfile + "/" + ttfile)
	ttree = tfile.Get("posteriors")

	ttree.SetBranchAddress("sin2th_13", AddressOf(mcmc_tree, "sin2th_13"))
	ttree.SetBranchAddress("sin2th_23", AddressOf(mcmc_tree, "sin2th_23"))
	ttree.SetBranchAddress("delm2_12", AddressOf(mcmc_tree, "delm2_12"))
	ttree.SetBranchAddress("delm2_23", AddressOf(mcmc_tree, "delm2_23"))
	ttree.SetBranchAddress("delta_cp", AddressOf(mcmc_tree, "delta_cp"))
	ttree.SetBranchAddress("accProb", AddressOf(mcmc_tree, "accProb"))
	ttree.SetBranchAddress("LogL", AddressOf(mcmc_tree, "LogL"))
	
	if(dud_check):
		ttree.GetEntry(100)
		check = mcmc_tree.sin2th_23
		print "sin23 " + str(check)
		if(float(check) >= float(1)):
			print "GREATER THAN 1"
			bad_files.append(str(ttfile))
			continue
	
		#if (int(check) is int(1)):
	
		"""
	for lol in range(10000):
		ttree.GetEntry(lol)
		print "llh " + str(lol) + " " + str(mcmc_tree.LogL)
		llh_13.Fill(float(lol),float(mcmc_tree.LogL))

	canvg = TCanvas()
	llh_13.GetXaxis().SetTitle("MCMC Step Number")
	llh_13.GetYaxis().SetTitle("Neg. Log Likelihood")
	llh_13.GetYaxis().SetTitleOffset(1.3)
	llh_13.Draw()
	canvg.SaveAs(webloc + "/llhsteps.png")
	break
#check = 0
	
#else:
	#	break
"""
	for e in range(nburn, ttree.GetEntries()):
		ttree.GetEntry(e)
		h_13_cp.Fill(mcmc_tree.sin2th_13, mcmc_tree.delta_cp)
		h_23_cp.Fill(mcmc_tree.sin2th_23, mcmc_tree.delta_cp)
		h_d23_23.Fill(mcmc_tree.delm2_23, mcmc_tree.sin2th_23)
		h_d23_13.Fill(mcmc_tree.delm2_23, mcmc_tree.sin2th_13)
		h_23_13.Fill(mcmc_tree.sin2th_23, mcmc_tree.sin2th_13)
		h_d23_cp.Fill(mcmc_tree.delm2_23, mcmc_tree.delta_cp)

		h_sin13.Fill(mcmc_tree.sin2th_13)
		h_accProb.Fill(mcmc_tree.accProb)
		h_LogL.Fill(mcmc_tree.LogL)
		h_13_LogL.Fill(mcmc_tree.sin2th_13, mcmc_tree.LogL)
		h_23_LogL.Fill(mcmc_tree.sin2th_23, mcmc_tree.LogL)
		
	############################################

	hpost = TH1D("hpost","hpost",100,0.5,1)

	# sin12                                                 
	hpost.SetBins(binning, sin12_l, sin12_u)
	ttree.Project("hpost", "sin2th_12", "", "", ttree.GetEntries(), nburn)
	hpost.SetTitle("sin2th_12")
	h_sin12.Add(hpost)

	# sin23                                         
	hpost.SetBins(binning, sin23_l, sin23_u)
	ttree.Project("hpost", "sin2th_23", "", "", ttree.GetEntries(), nburn)
	hpost.SetTitle("sin2th_23")
	h_sin23.Add(hpost)

	# sin13                                                                    
	#hpost.SetBins(binning, sin13_l, sin13_u)
	#ttree.Project("hpost", "sin2th_13", "", "", ttree.GetEntries(), nburn)
	#hpost.SetTitle("sin2th_13")
	#h_sin13.Add(hpost)

	# delm12                 
	hpost.SetBins(binning, delm12_l, delm12_u)
	ttree.Project("hpost", "delm2_12", "", "", ttree.GetEntries(), nburn)
	hpost.SetTitle("delm2_12")
	h_delm12.Add(hpost)
	
	# delm23                                              
	hpost.SetBins(binning, delm23_l, delm23_u)
	ttree.Project("hpost", "delm2_23", "", "", ttree.GetEntries(), nburn)
	hpost.SetTitle("delm2_23")
	h_delm23.Add(hpost)

	# deltacp                                    
	hpost.SetBins(binning, deltacp_l, deltacp_u)
	ttree.Project("hpost", "delta_cp", "", "", ttree.GetEntries(), nburn)
	hpost.SetTitle("delta_cp")
	h_deltacp.Add(hpost)

# owaru ################

ROOT.gStyle.SetOptStat(0)

nudge = 2.5
fcanv = TCanvas("","",int(600*nudge),int(600*nudge))
#fcanv.SetGrid()
#fcanv.Print(outname + ".ps[")
fcanv.Divide(4,4)

fcanv.cd(8)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
#h_sin12.SetTitle("sin2th_12")
h_sin12.GetXaxis().SetTitle("sin^{2}(#theta_{12})")
h_sin12.Draw()

fcanv.cd(11)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
#h_sin13.SetTitle("sin2th_13")
h_sin13.GetXaxis().SetTitle("sin^{2}(#theta_{13})")
h_sin13.Draw()
#fcanv.Print(outname + ".ps")

fcanv.cd(1)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
ROOT.gPad.SetBottomMargin(0.15)
#h_delm23.SetTitle("delm23")
h_delm23.GetXaxis().SetNdivisions(nDiv_delm23)
h_delm23.GetXaxis().SetTitle("#Deltam^{2}_{23}")
a,b,c = histGausFit(h_delm23)

tbox_delm23 = TPaveText(0.125, 0.7, 0.425, 0.87, "NDC")
tbox_delm23.SetFillColor(10)
m = TString("mean - " + str("%.3e"%b))
s = TString("sigma - " + str("%.3e"%c))
tbox_delm23.AddText(str(m))
tbox_delm23.AddText(str(s))

h_delm23.Draw()
tbox_delm23.Draw()

ROOT.gPad.Update()
#fcanv.Print(outname + ".ps")

fcanv.cd(4)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
#h_delm12.SetTitle("delm12")
h_delm12.GetXaxis().SetTitle("#Deltam^{2}_{12}")
a,b,c = histGausFit(h_delm12)

tbox_delm12 = TPaveText(0.125, 0.7, 0.425, 0.87, "NDC")
tbox_delm12.SetFillColor(10)
m = TString("mean - " + str("%.3e"%b))
s = TString("sigma - " + str("%.3e"%c))
tbox_delm12.AddText(str(m))
tbox_delm12.AddText(str(s))

h_delm12.Draw()
tbox_delm12.Draw()
ROOT.gPad.Update()

#fcanv.Print(outname + ".ps")

ROOT.gStyle.SetOptStat("e")
ROOT.gStyle.SetPalette(1)

fcanv.cd(15)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
#h_13_cp.GetYaxis().SetTitleOffset(t_off)
h_13_cp.GetXaxis().SetTitle("sin^{2}(#theta_{13})")
h_13_cp.GetYaxis().SetTitle("#delta_{cp}")
h_13_cp.Draw("COLZ")

fcanv.cd(14)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
#h_23_cp.GetYaxis().SetTitleOffset(t_off)
h_23_cp.GetXaxis().SetTitle("sin^{2}(#theta_{23})")
h_23_cp.GetYaxis().SetTitle("#delta_{cp}")
h_23_cp.Draw("COLZ")

fcanv.cd(5)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
#ROOT.gPad.SetLogz()
h_d23_23.GetYaxis().SetTitleOffset(t_off)
h_d23_23.GetXaxis().SetNdivisions(nDiv_delm23)
h_d23_23.GetXaxis().SetTitle("#Deltam^{2}_{23}")
h_d23_23.GetYaxis().SetTitle("sin^{2}(#theta_{23})")
h_d23_23.Draw("COLZ")

fcanv.cd(6)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
#h_sin23.SetTitle("sin2th_23")
h_sin23.GetXaxis().SetTitle("sin^{2}(#theta_{23})")
h_sin23.Draw()


fcanv.cd(9)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
#ROOT.gPad.SetLogz()
h_d23_13.GetYaxis().SetTitleOffset(t_off)
h_d23_13.GetXaxis().SetNdivisions(nDiv_delm23)
h_d23_13.GetXaxis().SetTitle("#Deltam^{2}_{23}")
h_d23_13.GetYaxis().SetTitle("sin^{2}(#theta_{13})")
h_d23_13.Draw("COLZ")

fcanv.cd(10)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
#ROOT.gPad.SetLogz()
h_23_13.GetYaxis().SetTitleOffset(t_off)
h_23_13.GetXaxis().SetTitle("sin^{2}(#theta_{23})")
h_23_13.GetYaxis().SetTitle("sin^{2}(#theta_{13})")
h_23_13.Draw("COLZ")

fcanv.cd(13)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()
h_d23_cp.GetXaxis().SetNdivisions(nDiv_delm23)
h_d23_cp.GetXaxis().SetTitle("#Deltam^{2}_{23}")
h_d23_cp.GetYaxis().SetTitle("#delta_{cp}")
h_d23_cp.Draw("COLZ")

fcanv.cd(16)
ROOT.gPad.SetGridx()
ROOT.gPad.SetGridy()

h_deltacp.GetXaxis().SetTitle("#delta_{cp}")
#h_deltacp.SetTitle("delta_cp")
h_deltacp.Draw()


fcanv.SaveAs(webloc + "/in_progress.png")

ROOT.gStyle.SetOptStat(1101)

canv2 = TCanvas("l","",375*5, 375)
canv2.Divide(5)
canv2.cd(1)
h_accProb.Draw()
canv2.cd(2)
h_LogL.Draw()
canv2.cd(3)
ROOT.gPad.SetGrid()
h_13_LogL.Draw("COLZ")
canv2.cd(4)
ROOT.gPad.SetGrid()
h_23_LogL.Draw("COLZ")
canv2.cd(5)
h_projection = h_23_13.ProjectionY("proj",95)
h_projection.Rebin()
ROOT.gPad.SetGrid();
h_projection.SetTitle("Projected #theta_{13}")
h_projection.Draw()
#h_13.Draw("SAME")
canv2.SaveAs(webloc + "/fit_misc.png")

#canvg = TCanvas()
#llh_13.Draw()
#canvg.SaveAs(webloc + "/llhsteps.png")
#os.system("ps2pdf " + outname + ".ps")
#os.system("cp " + outname + ".pdf " + str(webloc))

# bad files
print "\nthe following files look bad:"
for f in bad_files:
	print f
