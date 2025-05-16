#!/usr/bin/env python

import sys
import math
import os
import logging
import getopt
import fnmatch

import ROOT
from ROOT import gROOT, TCanvas, TH1D, TH2D, TStyle, TFile, TApplication, TString, TChain, gSystem, gDirectory, TGraph, TMultiGraph, TH1F, TF1

from random import choice
from array import array

if len(sys.argv) != 3:
        print "Usage: " + sys.argv[0] + " /path/of/files/ outname"
        sys.exit(1)

print "summing files from " + sys.argv[1]
#hpost = TH1D("hpost","hpost",100,0.5,1)

# config ################
webloc = "/user2/rcalland/www/MaCh3"
nburn = 0 # adjust burn in steps
binning = 100
nflux = 32
ndet = 19
nxsec = 23
nNDdet = 39

# histograms ############

h_flux = []
h_xsec = []
h_skdet = []
h_nddet = []

#########################

totalfile = sys.argv[1]
outname = sys.argv[2]

#print totalfile

#tfile = TFile(totalfile);

#hpost = TH1D("hpost","hpost",100,0.5,1)

# print flux params #####
flux_pre = "b_"
fupper = 2
flower = 0

for i in range(nflux + 1):
        histo = TH1D(flux_pre + str(i),"",binning, flower, fupper)
	h_flux.append(histo)
	
ROOT.gStyle.SetOptStat(0)

# print det #############
ndet_pre = "skd_comb_"
skupper = 5
sklower = -5

for i in range(ndet +1):
        histo = TH1D(ndet_pre + str(i),"",binning, sklower, skupper)
        h_skdet.append(histo)

# print xsec ##############
xsec_pre = "xsec_"
xsecupper = 5
xseclower = -5

for i in range(nxsec +1):
	histo = TH1D(xsec_pre + str(i),"",binning, xseclower, xsecupper)
        h_xsec.append(histo)

# print xsec ##############                                                                                                                                               
ndd_pre = "ndd_"
nddupper = 2
nddlower = 0

for i in range(nNDdet +1):
	histo = TH1D(ndd_pre + str(i),"",binning, nddlower, nddupper)
        h_nddet.append(histo)


# osc params #################

sin12_u = 1
sin12_l = 0.5
h_sin12 = TH1D("sin12","",binning,sin12_l, sin12_u)

sin23_u = 1
sin23_l = 0.9
h_sin23 = TH1D("sin23","",binning, sin23_l, sin23_u)

sin13_u = 0.15
sin13_l = 0.05
h_sin13 = TH1D("sin13","",binning, sin13_l, sin13_u)

delm12_u = 7.5001E-5
delm12_l = 7.4999E-5
h_delm12 = TH1D("delm12","",binning, delm12_l, delm12_u)

delm23_u = 2.4E-3
delm23_l = 2.385E-3
h_delm23 = TH1D("delm23","",binning, delm23_l, delm23_u)

deltacp_u = 1
deltacp_l = -1
h_deltacp = TH1D("deltacp","",binning, deltacp_l, deltacp_u)

# file loop
dir_list = os.listdir(totalfile)

for ttfile in dir_list:
	if (fnmatch.fnmatch(ttfile,"*.root") is not True):
		continue
	tfile = TFile(totalfile + "/" + ttfile)
	ttree = tfile.Get("posteriors")
	hpost = TH1D("hpost","hpost",100,0.5,1)

	
	# flux
	for i in range(nflux + 1):
		name = flux_pre + str(i)
		hpost.SetBins(binning, flower, fupper)
		ttree.Project("hpost", name, "", "", ttree.GetEntries(), nburn)
		hpost.SetTitle(name)
		h_flux[i].Add(hpost)
	
	for i in range(ndet + 1):
		name = ndet_pre + str(i)
		hpost.SetBins(binning, sklower, skupper)
		ttree.Project("hpost", name, "", "", ttree.GetEntries(), nburn)
		hpost.SetTitle(name)
		h_skdet[i].Add(hpost)

	for i in range(nxsec + 1):
	       	name = xsec_pre + str(i)
		hpost.SetBins(binning, xseclower, xsecupper)
		ttree.Project("hpost", name, "", "", ttree.GetEntries(), nburn)
		hpost.SetTitle(name)
		h_xsec[i].Add(hpost)

	for i in range(nNDdet + 1):
		name = ndd_pre + str(i)
		hpost.SetBins(binning, nddlower, nddupper)
		ttree.Project("hpost", name, "", "", ttree.GetEntries(), nburn)
		hpost.SetTitle(name)
		h_nddet[i].Add(hpost)

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
	hpost.SetBins(binning, sin13_l, sin13_u)
	ttree.Project("hpost", "sin2th_13", "", "", ttree.GetEntries(), nburn)
	hpost.SetTitle("sin2th_13")
	h_sin13.Add(hpost)

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

fcanv = TCanvas()
fcanv.Print(outname + ".ps[")

for i in range(nflux + 1):
	h_flux[i].SetTitle(flux_pre + str(i))
	h_flux[i].Draw()
	fcanv.Print(outname + ".ps")

for i in range(ndet + 1):
	h_skdet[i].SetTitle(ndet_pre + str(i))
	h_skdet[i].Draw()
	fcanv.Print(outname + ".ps")

for i in range(nxsec + 1):
	h_xsec[i].SetTitle(xsec_pre + str(i))
	h_xsec[i].Draw()
	fcanv.Print(outname + ".ps")

for i in range(nNDdet + 1):
	h_nddet[i].SetTitle(ndd_pre + str(i))
	h_nddet[i].Draw()
	fcanv.Print(outname + ".ps")

h_sin12.SetTitle("sin2th_12")
h_sin12.Draw()
fcanv.Print(outname + ".ps")

h_sin23.SetTitle("sin2th_23")
h_sin23.Draw()
fcanv.Print(outname + ".ps")

h_sin13.SetTitle("sin2th_13")
h_sin13.Draw()
fcanv.Print(outname + ".ps")

h_delm23.SetTitle("delm23")
h_delm23.Draw()
fcanv.Print(outname + ".ps")

h_delm12.SetTitle("delm12")
h_delm12.Draw()
fcanv.Print(outname + ".ps")

h_deltacp.SetTitle("delta_cp")
h_deltacp.Draw()
fcanv.Print(outname + ".ps")

fcanv.Print(outname + ".ps]")


os.system("ps2pdf " + outname + ".ps")
os.system("cp " + outname + ".pdf " + str(webloc))
