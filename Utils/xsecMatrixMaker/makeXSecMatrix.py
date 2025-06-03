#!/usr/bin/python

import xml.etree.ElementTree as ET
import ROOT
import math
import sys

# LOL DOCUMENTATION
# ACK, Oct 2016
#
# Produces a MaCh3 compatible matrix for cross section systematics
# First argument is the xml file
# Second argument is the output root file
# Code is hugely hacked and ugly, but totes functional
#
# ETA - Oct 2021
# Adding in a few elements to read in generic xsec kinematic strings that
# then gets converted within MaCh3 to get the kinematic variable for an event
# which then gets cut on

if len(sys.argv) != 3:
  print("Sorry, I need two arguments")
  print("./makeXSecMatrix.py input.xml output.root")
  sys.exit()

tree = ET.parse(sys.argv[1])
fromxml = tree.getroot()

maxelements=320
nummodes=12


xsec_param_names        = ROOT.TObjArray()
xsec_param_prior_unnorm = ROOT.TVectorD(maxelements)
xsec_param_nom_unnorm   = ROOT.TVectorD(maxelements)
xsec_param_prior        = ROOT.TVectorD(maxelements)
xsec_param_nom          = ROOT.TVectorD(maxelements)
xsec_param_lb           = ROOT.TVectorD(maxelements)
xsec_param_ub           = ROOT.TVectorD(maxelements)
xsec_error              = ROOT.TVectorD(maxelements)
xsec_stepscale          = ROOT.TVectorD(maxelements)
xsec_param_id           = ROOT.TMatrixD(maxelements, 2)

xsec_fd_spline_names    = ROOT.TObjArray()
xsec_fd_spline_modes    = ROOT.TObjArray()

xsec_nd_spline_names    = ROOT.TObjArray()

xsec_norm_modes         = ROOT.TObjArray()
xsec_norm_elements      = ROOT.TObjArray()
xsec_norm_nupdg         = ROOT.TObjArray()
xsec_norm_prod_nupdg         = ROOT.TObjArray()
xsec_norm_horncurrents         = ROOT.TObjArray()
#ETA - adding in generic kinematic type
xsec_norm_kinematic_type = ROOT.TObjArray()
xsec_norm_kinematic_lb = ROOT.TVectorD(maxelements) 
xsec_norm_kinematic_ub = ROOT.TVectorD(maxelements)


corr_dicts = []
corr_names = []

kinematic_dicts = []


nelem = 0
for child in fromxml:

    # Get the name attribute
    name = ROOT.TObjString(child.attrib['name'])


    # Add the name to the TObjArray
    xsec_param_names.AddLast(name)
    corr_names.append(name)

    # Get the prior from attribute
    prior = float(child.attrib['prior'])
    xsec_param_prior_unnorm[nelem]=float(prior)

    # Get the nominal from attribute
    nom = float(child.attrib['nom'])
    xsec_param_nom_unnorm[nelem]=float(nom)

    # Get the lower bound from attribute
    lb = float(child.attrib['lb'])
    ub = float(child.attrib['ub'])
    error = float(child.attrib['error'])
    renorm = child.attrib['renorm']

    stepscale = float(child.attrib['stepscale'])
    xsec_stepscale[nelem] = float(stepscale)

    # Read in what to of detector the parameter applies to
    # This is meant to be used wtih binary operator in MaCh3, e.g. 25 means apply both to ND280 and SK, 1 means apply only to ND280, etc
    # Is used because for some analyses we used parameters only for numu samples at SK, which had detid 8
    # Could be switched to simple int comparisons...
    xsec_param_id[nelem][1]=int(child.attrib['detid'])

    # If we want to renormalize so that the NEUT default = 1
    if (int(renorm) == 1):
        xsec_param_nom[nelem] = nom/nom
        xsec_param_prior[nelem] = prior/nom
        xsec_error[nelem] = error/nom
        if (lb != -9999):
            xsec_param_lb[nelem] = lb/nom
        else:
            xsec_param_lb[nelem] = lb
        if (ub!= 9999):
            xsec_param_ub[nelem] = ub/nom
        else:
            xsec_param_ub[nelem] = ub
    else:
        xsec_param_nom[nelem] = nom
        xsec_param_prior[nelem] = prior
        xsec_param_lb[nelem] = lb
        xsec_param_ub[nelem] = ub
        xsec_error[nelem] = error

    # Find the correlations
    correlations = child.findall('correlation')
    dict = {}
    for cor in correlations:
        dict[cor.attrib['par']]=float(cor.text)
    corr_dicts.append(dict)

    # Find what type of parameter we're dealing with
    type = str(child.attrib['type'])


    # Spline parameters (e.g. MAQE)

    #All variables need the spline par entry


    if(type=="spline"):
        xsec_param_id[nelem][0]=int(child.attrib['splineind'])
        # Get the SK spline file name attribute
        print(xsec_param_id[nelem][1])

        if 'fd_spline_name' in child.attrib:
          fd_spline_name = ROOT.TObjString(child.attrib['fd_spline_name'])
        else: 
          fd_spline_name=ROOT.TObjString("")
        xsec_fd_spline_names.AddAtAndExpand(fd_spline_name, nelem)
        if "_" in str(fd_spline_name):
          print("ERROR : ", fd_spline_name, " contains _ (this is a problem for MaCh3!)")
          sys.exit()


        # Find the mode which this parameter might apply to
        spline_modes = child.findall('fd_mode')
        if(len(spline_modes)==0):
          #print "Found no mode list for ",name," filling with blank which means apply to all"
          dummyvec = ROOT.TVectorD(0)
          xsec_fd_spline_modes.AddAtAndExpand(dummyvec,nelem)
        else:
          for mode in spline_modes:
            if mode.text == None or mode.tail == None:
              #print "Found empty mode list for",name," filling with blank which means apply to all"
              mode.text= ""
            # Split multiple entries up
            modelist = mode.text.split()
            # Make a TVectorD containing each mode
            modevec = ROOT.TVectorD(len(modelist))
            for i in range(len(modelist)):
              modevec[i]=int(modelist[i])
            xsec_fd_spline_modes.AddAtAndExpand(modevec, nelem)

        if 'nd_spline_name' in child.attrib:
          nd_spline_name = ROOT.TObjString(child.attrib['nd_spline_name'])
        else: 
          nd_spline_name=ROOT.TObjString("")
        xsec_nd_spline_names.AddAtAndExpand(nd_spline_name, nelem)


    #End Spline pars

    # Normalisation parameter (e.g. 2p2h neutrino normalisation)
    elif(type=="norm"):
        xsec_param_id[nelem][0]=-1
    # Functional parameter (e.g. BeRPA)
    elif(type=="function"):
        xsec_param_id[nelem][0]=-2
    else:
        print("Wrong parameter type!")
        quit()

    # Find the mode which this parameter might apply to
    normmodes = child.findall('mode')
    if(len(normmodes)==0 and type=="norm"):
        print("Found no mode list for ",name," filling with blank which means apply to all")
        dummyvec = ROOT.TVectorD(0)
        xsec_norm_modes.AddAtAndExpand(dummyvec,nelem)
    else:
      for mode in normmodes:
        if mode.text == None or mode.tail == None:
          print("Found empty mode list for",name," filling with blank which means apply to all")
          mode.text= ""
        # Split multiple entries up
        modelist = mode.text.split()
        # Make a TVectorD containing each mode
        modevec = ROOT.TVectorD(len(modelist))
        for i in range(len(modelist)):
            modevec[i]=int(modelist[i])
        xsec_norm_modes.AddAtAndExpand(modevec, nelem)

    # Find the target element
    normelements = child.findall('element')
    if(len(normelements)==0 and type=="norm"):
        print("Found no target list for ",name," filling with blank which means apply to all")
        dummyvec = ROOT.TVectorD(0)
        xsec_norm_elements.AddAtAndExpand(dummyvec,nelem)
    else:
      for elem in normelements:
        if  elem.text == None or elem.tail == None:
          print("Found empty target list for ",name," filling with blank which means apply to all")
          elem.text = ""
        elemlist = elem.text.split()
        elemvec = ROOT.TVectorD(len(elemlist))
        for i in range(len(elemlist)):
            elemvec[i]=int(elemlist[i])
        xsec_norm_elements.AddAtAndExpand(elemvec, nelem)

    # Find the mode which this parameter might apply to
    horncurrents = child.findall('horncurrent')
    if(len(horncurrents)==0 and type=="norm"):
        print("Found no horn current list for ",name," filling with blank which means apply to all")
        dummyvec = ROOT.TVectorD(0)
        xsec_norm_horncurrents.AddAtAndExpand(dummyvec,nelem)
    else:
      for horncurrent in horncurrents:
        if horncurrent.text == None or horncurrent.tail == None:
          print("Found empty horncurrent list for",name," filling with blank which means apply to all")
          horncurrent.text= ""
        # Split multiple entries up
        horncurrentlist = horncurrent.text.split()
        # Make a TVectorD containing each mode
        horncurrentvec = ROOT.TVectorD(len(horncurrentlist))
        for i in range(len(horncurrentlist)):
            horncurrentvec[i]=int(horncurrentlist[i])
        xsec_norm_horncurrents.AddAtAndExpand(horncurrentvec, nelem)

    # Find the neutrino type
    normpdg = child.findall('nupdg')
    if(len(normpdg)==0 and type=="norm"):
        print("Found no nupdg list for ",name," filling with blank which means apply to all")
        dummyvec = ROOT.TVectorD(0)
        xsec_norm_nupdg.AddAtAndExpand(dummyvec,nelem)
    else:
      for pdg in normpdg:
        if pdg.text == None or pdg.tail == None:
          print("Found empty nupdg list for ",name," filling with blank which means apply to all")
          pdg.text = ""
        pdglist = pdg.text.split()
        pdgvec = ROOT.TVectorD(len(pdglist))
        for i in range(len(pdglist)):
            pdgvec[i]=int(pdglist[i])
        xsec_norm_nupdg.AddAtAndExpand(pdgvec, nelem)

    normprodpdg = child.findall('prod_nupdg')
    if(len(normprodpdg)==0 and type=="norm"):
        print("Found no prod_nupdg list for ",name," filling with blank which means apply to all")
        dummyvec = ROOT.TVectorD(0)
        xsec_norm_prod_nupdg.AddAtAndExpand(dummyvec,nelem)
    else:
      for prodpdg in normprodpdg:
        if prodpdg.text == None or prodpdg.tail == None:
          print("Found empty prod_nupdg list for ",name," filling with blank which means apply to all")
          prodpdg.text = ""
        prodpdglist = prodpdg.text.split()
        prodpdgvec = ROOT.TVectorD(len(prodpdglist))
        for i in range(len(prodpdglist)):
            prodpdgvec[i]=int(prodpdglist[i])
        xsec_norm_prod_nupdg.AddAtAndExpand(prodpdgvec, nelem)

	#Find if there are any kinematic cuts
    kin_cuts = child.findall('kinematic_cut')
    if(len(kin_cuts)!=0 and type=="norm"):
	  #print "Didn't find any kinematic_cuts"  
      print("FOUND kinematic cut")
      for var in kin_cuts:
        xsec_norm_kinematic_type.AddLast(ROOT.TObjString(var.attrib['var']))
        bnds_list = var.text.split()
        if(len(bnds_list) !=2):
          print("ERROR: there should only be two bounds for a kinamatic cut; the upper and lower") 
          sys.exit()
        else:
          xsec_norm_kinematic_lb[nelem] = float(bnds_list[0])
          print("Found lower bound of ", xsec_norm_kinematic_lb[nelem], " and put it in element ", nelem)
          xsec_norm_kinematic_ub[nelem] = float(bnds_list[1])
          print("Found upper bound of ", xsec_norm_kinematic_ub[nelem])

          print("Name of kinematic var is ",var.attrib['var'])
    elif(type=="norm"):
      print("NO KINEMATIC CUT FOUND!!")
      xsec_norm_kinematic_type.AddLast(ROOT.TObjString(""))


    nelem+=1 #Only have this line once par par

# Resize the vectors to a reasonable number of elements
xsec_param_nom_unnorm.ResizeTo(nelem)
xsec_param_nom.ResizeTo(nelem)
xsec_param_prior_unnorm.ResizeTo(nelem)
xsec_param_prior.ResizeTo(nelem)
xsec_param_lb.ResizeTo(nelem)
xsec_param_ub.ResizeTo(nelem)
xsec_stepscale.ResizeTo(nelem)
xsec_param_id.ResizeTo(nelem, 2)
xsec_norm_kinematic_ub.ResizeTo(nelem)
xsec_norm_kinematic_lb.ResizeTo(nelem)

xsec_cov = ROOT.TMatrixD(nelem, nelem)

# Loop over the elements
for i in range(nelem):

    # The diagonal entries in the covariance is simply the error squared
    xsec_cov[i][i] = xsec_error[i]*xsec_error[i]

    for item in corr_dicts[i]:
        # Script will fail if not in list of parameters; this is desired behavior. Could do more neatly, but meh
        index = corr_names.index(item)
        corr1 = corr_dicts[i][item]
        corr2 = 0;
        if (str(corr_names[i]) in list(corr_dicts[index].keys())):
            corr2 = corr_dicts[index][str(corr_names[i])]
        if (round(corr1, 12) == round(corr2, 12)):
            xsec_cov[i][index] = xsec_cov[index][i] = round(corr1, 12) * xsec_error[i] * xsec_error[index]
        elif (corr2==0):
            print("Warning: no reciprocal correlation between "+str(item)+" and "+str(corr_names[i])+" proceeding with non-reciprocal")
            xsec_cov[i][index]=xsec_cov[index][i] = corr1 * xsec_error[i] * xsec_error[index]
        else:
            print("Correlations between "+str(item)+" and "+str(corr_names[i])+" don't match. Quitting!")
            print("corr 1 = ", corr1, "corr 2 = ", corr2)
            quit()
        

# Make the TH2D covariance matrix
hcov = ROOT.TH2D("hcov", "", nelem, 0, nelem, nelem, 0, nelem)
# Set the content
for i in range(nelem):
    for j in range(nelem):
        if (xsec_cov[i][j] >= 0):
            hcov.SetBinContent(i+1, j+1, math.sqrt(xsec_cov[i][j]))
        else:
            hcov.SetBinContent(i+1, j+1, -1*math.sqrt(abs(xsec_cov[i][j])))

for i in range(nelem):
    hcov.GetXaxis().SetBinLabel(i+1, str(corr_names[i]))
    hcov.GetYaxis().SetBinLabel(i+1, str(corr_names[i]))

# Set the covariance matrix pretty
hcov.SetMaximum(1)
hcov.SetMinimum(-1)
hcov.GetXaxis().LabelsOption("v")
hcov.GetXaxis().SetLabelSize(0.03)
hcov.GetYaxis().SetLabelSize(0.03)
hcov.GetZaxis().SetTitle("#sqrt{|V_{ij}|}#times sign(V_{ij})")

# Write the output file
file = ROOT.TFile(sys.argv[2], "RECREATE")
xsec_param_names.Write("xsec_param_names", 1)
xsec_fd_spline_names.Write("fd_spline_names", 1);
xsec_fd_spline_modes.Write("fd_spline_modes",1);
xsec_nd_spline_names.Write("nd_spline_names", 1);
xsec_param_nom.Write("xsec_param_nom")
xsec_param_nom_unnorm.Write("xsec_param_nom_unnorm")
xsec_param_prior.Write("xsec_param_prior")
xsec_param_prior_unnorm.Write("xsec_param_prior_unnorm")
xsec_param_lb.Write("xsec_param_lb")
xsec_param_ub.Write("xsec_param_ub")
xsec_param_id.Write("xsec_param_id")
xsec_stepscale.Write("xsec_stepscale")
xsec_norm_kinematic_type.Write("xsec_norm_kinematic_type", 1)
xsec_norm_kinematic_lb.Write("xsec_norm_kinematic_lb")
xsec_norm_kinematic_ub.Write("xsec_norm_kinematic_ub")
xsec_norm_modes.Write("xsec_norm_modes", 1)
xsec_norm_horncurrents.Write("xsec_norm_horncurrents", 1)
xsec_norm_elements.Write("xsec_norm_elements", 1)
xsec_norm_nupdg.Write("xsec_norm_nupdg", 1)
xsec_norm_prod_nupdg.Write("xsec_norm_prod_nupdg", 1)
xsec_cov.Write("xsec_cov")
hcov.Write("hcov")
file.Close()

#filename=sys.argv[1];
#filename = filename.replace('.xml', '.tex')
#print filename
#texfile = open(filename, 'w')

#for i in range(nelem):
#    snomu =  "%.2f" % xsec_param_nom_unnorm[i]
#    snom =  "%.2f" % xsec_param_nom[i]
#    sprioru =  "%.2f" % xsec_param_prior_unnorm[i]
#    sprior =  "%.2f" % xsec_param_prior[i]
#    serroru = ""
#    if(snomu==snom):
#        serroru =  "%.2f" % xsec_error[i];
#    else:    
#        serroru =  "%.2f" % (xsec_error[i]*xsec_param_nom[i])
#    serror =  "%.2f" % xsec_error[i]
#    slb =  "%.2f" % xsec_param_lb[i]
#    sub =  "%.2f" % xsec_param_ub[i]
#    texfile.write(sprioru+"&"+serroru+"&"+snomu+"&"+sprior+"&"+serror+"&"+slb+"&"+sub+"\\\\\n")

