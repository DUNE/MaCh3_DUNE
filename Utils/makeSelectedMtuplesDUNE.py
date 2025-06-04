import os, sys
from ROOT import *
import argparse


if __name__ == "__main__":
    
    # Get arguments (ie. SK file name)
    parser = argparse.ArgumentParser(description='Splits an SK file into five files, containing the numu-,nue-, nue1pi- numucc1pi- selected events')
    parser.add_argument('filename', nargs=1)
    args = parser.parse_args()
    skfile = args.filename[0]

    print "Working on file ", skfile

    # Get SK tree
    f_sk = TFile(skfile, "OPEN")
    if not f_sk:
        print "Error: could not open file ", skfile
        sys.exit()
    sk_tree = f_sk.Get("caf")
    if not sk_tree:
        sk_tree = f_sk.Get("h1")

    # Get entries in SK tree
    n_total = sk_tree.GetEntries()
    print n_total, " events in tree."

    # Get the SK filename (without '.root')
    # skfile = skfile[:-5]

    # Make new filenames for files containing only the numu/nue trees
    numu_numu_file  = skfile[:-13]+"_numu_x_numu.root"
    numu_nue_file   = skfile[:-13]+"_numu_x_nue.root"
    numu_nutau_file = skfile[:-13]+"_numu_x_nutau.root"
    nue_nue_file    = skfile[:-13]+"_nue_x_nue.root"
    nue_numu_file   = skfile[:-13]+"_nue_x_numu.root"
    nue_nutau_file  = skfile[:-13]+"_nue_x_nutau.root"
    numubar_numubar_file  = skfile[:-13]+"_numubar_x_numubar.root"
    numubar_nuebar_file   = skfile[:-13]+"_numubar_x_nuebar.root"
    numubar_nutaubar_file = skfile[:-13]+"_numubar_x_nutaubar.root"
    nuebar_nuebar_file    = skfile[:-13]+"_nuebar_x_nuebar.root"
    nuebar_numubar_file   = skfile[:-13]+"_nuebar_x_numubar.root"
    nuebar_nutaubar_file  = skfile[:-13]+"_nuebar_x_nutaubar.root"
    print "Making new files: "
    print numu_numu_file
    print numu_nue_file
    print numu_nutau_file
    print nue_nue_file
    print nue_numu_file
    print nue_nutau_file
    print numubar_numubar_file
    print numubar_nuebar_file
    print numubar_nutaubar_file
    print nuebar_nuebar_file
    print nuebar_numubar_file
    print nuebar_nutaubar_file
    print "Check: number of events in original SK tree = ", n_total

    # Make numu_numu tree and save to file
    f_numu_numu = TFile(numu_numu_file,"RECREATE")
    numu_numu_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==14 and sk_tree.nuPDG==14: #numu_numu selection
            numu_numu_tree.Fill()

    n_numu_numu = numu_numu_tree.GetEntries()
    print "Check: number of events in numu_numu tree = ", n_numu_numu	
    numu_numu_tree.Write("caf")
    f_numu_numu.Close()
    if n_numu_numu == 0:
       os.remove(numu_numu_file)


    # Make numu_nue tree and save to file
    f_numu_nue = TFile(numu_nue_file,"RECREATE")
    numu_nue_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==14 and sk_tree.nuPDG==12: #numu_nue selection
            numu_nue_tree.Fill()

    n_numu_nue = numu_nue_tree.GetEntries()
    print "Check: number of events in numu_nue tree = ", n_numu_nue

    numu_nue_tree.Write("caf")
    f_numu_nue.Close()
    if n_numu_nue == 0:
       os.remove(numu_nue_file)

    # Make numu_nutau tree and save to file
    f_numu_nutau = TFile(numu_nutau_file,"RECREATE")
    numu_nutau_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==14 and sk_tree.nuPDG==16: #numu_nutau selection
            numu_nutau_tree.Fill()

    n_numu_nutau = numu_nutau_tree.GetEntries()
    print "Check: number of events in numu_nutau tree = ", n_numu_nutau

    numu_nutau_tree.Write("caf")
    f_numu_nutau.Close()
    if n_numu_nutau == 0:
       os.remove(numu_nutau_file)
    
    
    # Make nue_nue tree and save to file
    f_nue_nue = TFile(nue_nue_file,"RECREATE")
    nue_nue_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==12 and sk_tree.nuPDG==12: #nue_nue selection
            nue_nue_tree.Fill()

    n_nue_nue = nue_nue_tree.GetEntries()
    print "Check: number of events in nue_nue tree = ", n_nue_nue

    nue_nue_tree.Write("caf")
    f_nue_nue.Close()
    if n_nue_nue == 0:
       os.remove(nue_nue_file)

    # Make nue_numu tree and save to file
    f_nue_numu = TFile(nue_numu_file,"RECREATE")
    nue_numu_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==12 and sk_tree.nuPDG==14: #nue_numu selection
            nue_numu_tree.Fill()

    n_nue_numu = nue_numu_tree.GetEntries()
    print "Check: number of events in nue_numu tree = ", n_nue_numu

    nue_numu_tree.Write("caf")
    f_nue_numu.Close()
    if n_nue_numu == 0:
       os.remove(nue_numu_file)

    # Make nue_nutau tree and save to file
    f_nue_nutau = TFile(nue_nutau_file,"RECREATE")
    nue_nutau_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==12 and sk_tree.nuPDG==16: #nue_nutau selection
            nue_nutau_tree.Fill()

    n_nue_nutau = nue_nutau_tree.GetEntries()
    print "Check: number of events in nue_nutau tree = ", n_nue_nutau

    nue_nutau_tree.Write("caf")
    f_nue_nutau.Close()
    if n_nue_nutau == 0:
       os.remove(nue_nutau_file)
    
    # Make numubar_numubar tree and save to file
    f_numubar_numubar = TFile(numubar_numubar_file,"RECREATE")
    numubar_numubar_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-14 and sk_tree.nuPDG==-14: #numubar_numubar selection
            numubar_numubar_tree.Fill()

    n_numubar_numubar = numubar_numubar_tree.GetEntries()
    print "Check: number of events in numubar_numubar tree = ", n_numubar_numubar

    numubar_numubar_tree.Write("caf")
    f_numubar_numubar.Close()
    if n_numubar_numubar == 0:
       os.remove(numubar_numubar_file)

    # Make numubar_nuebar tree and save to file
    f_numubar_nuebar = TFile(numubar_nuebar_file,"RECREATE")
    numubar_nuebar_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-14 and sk_tree.nuPDG==-12: #numubar_nuebar selection
            numubar_nuebar_tree.Fill()

    n_numubar_nuebar = numubar_nuebar_tree.GetEntries()
    print "Check: number of events in numubar_nuebar tree = ", n_numubar_nuebar

    numubar_nuebar_tree.Write("caf")
    f_numubar_nuebar.Close()
    if n_numubar_nuebar == 0:
       os.remove(numubar_nuebar_file)

    # Make numubar_nutaubar tree and save to file
    f_numubar_nutaubar = TFile(numubar_nutaubar_file,"RECREATE")
    numubar_nutaubar_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-14 and sk_tree.nuPDG==-16: #numubar_nutaubar selection
            numubar_nutaubar_tree.Fill()

    n_numubar_nutaubar = numubar_nutaubar_tree.GetEntries()
    print "Check: number of events in numubar_nutaubar tree = ", n_numubar_nutaubar

    numubar_nutaubar_tree.Write("caf")
    f_numubar_nutaubar.Close()
    if n_numubar_nutaubar == 0:
       os.remove(numubar_nutaubar_file)
    
    
    # Make nuebar_nuebar tree and save to file
    f_nuebar_nuebar = TFile(nuebar_nuebar_file,"RECREATE")
    nuebar_nuebar_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-12 and sk_tree.nuPDG==-12: #nuebar_nuebar selection
            nuebar_nuebar_tree.Fill()

    n_nuebar_nuebar = nuebar_nuebar_tree.GetEntries()
    print "Check: number of events in nuebar_nuebar tree = ", n_nuebar_nuebar

    nuebar_nuebar_tree.Write("caf")
    f_nuebar_nuebar.Close()
    if n_nuebar_nuebar == 0:
       os.remove(nuebar_nuebar_file)

    # Make nuebar_numubar tree and save to file
    f_nuebar_numubar = TFile(nuebar_numubar_file,"RECREATE")
    nuebar_numubar_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-12 and sk_tree.nuPDG==-14: #nuebar_numubar selection
            nuebar_numubar_tree.Fill()

    n_nuebar_numubar = nuebar_numubar_tree.GetEntries()
    print "Check: number of events in nuebar_numubar tree = ", n_nuebar_numubar

    nuebar_numubar_tree.Write("caf")
    f_nuebar_numubar.Close()
    if n_nuebar_numubar == 0:
       os.remove(nuebar_numubar_file)

    # Make nuebar_nutaubar tree and save to file
    f_nuebar_nutaubar = TFile(nuebar_nutaubar_file,"RECREATE")
    nuebar_nutaubar_tree = sk_tree.CloneTree(0)

    for i in xrange(0,n_total):
        sk_tree.GetEntry(i)
        if sk_tree.nuPDGunosc==-12 and sk_tree.nuPDG==-16: #nuebar_nutaubar selection
            nuebar_nutaubar_tree.Fill()

    n_nuebar_nutaubar = nuebar_nutaubar_tree.GetEntries()
    print "Check: number of events in nuebar_nutaubar tree = ", n_nuebar_nutaubar

    nuebar_nutaubar_tree.Write("caf")
    f_nuebar_nutaubar.Close()
    if n_nuebar_nutaubar == 0:
       os.remove(nuebar_nutaubar_file)
    



    print "Sum of event selections = ", n_numu_numu + n_numu_nue + n_numu_nutau + n_nue_nue + n_nue_numu + n_nue_nutau + n_numubar_numubar + n_numubar_nuebar + n_numubar_nutaubar +  n_nuebar_nuebar + n_nuebar_numubar + n_nuebar_nutaubar

    print "This should be equal to the total number of events: ", n_total 
