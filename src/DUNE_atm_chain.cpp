#include <iostream>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>

#include "samplePDFDUNE/samplePDFDUNEBase.h"

#include "mcmc/mcmc.h"


int main(int argc, char * argv[]) {

  // ----------------------- OPTIONS ---------------------------------------- //

  if(argc == 1){
    std::cout << "Usage: bin/mini_MCMC config.cfg" << std::endl;
    return 1;
  }

  manager *fitMan = new manager(argv[1]);


  std::string  XsecMatrixFile = fitMan->raw()["General"]["Systematics"]["XsecCovFile"].as<std::string>(); 
  std::string  XsecMatrixName = fitMan->raw()["General"]["Systematics"]["XsecCovName"].as<std::string>();
  std::string  OscMatrixFile = fitMan->raw()["General"]["Systematics"]["OscCovFile"].as<std::string>(); 
  std::string  OscMatrixName = fitMan->raw()["General"]["Systematics"]["OscCovName"].as<std::string>();
 
  // Asimov fit
  bool asimovfit = true;


  // ----------------------- COVARIANCE AND SAMPLEPDF OBJECTS ---------------------------------------- //

  gStyle -> SetPalette(1);

  // make file to save plots

  covarianceXsec *xsec = new covarianceXsec(XsecMatrixName.c_str(), XsecMatrixFile.c_str()) ;

  xsec->setParameters();

  xsec->setStepScale(0.01);

  // Oscillation covariance
  // covarianceOsc *osc = new covarianceOsc("osc_cov","inputs/oscillation_covariance_6par_nondouble_PDG2019.root");
  covarianceOsc *osc = new covarianceOsc("osc_cov","inputs/osc_covariance_DUNE_PDG2021_v1.root");
  osc->setStepScale(0.045);

  // oscpars from manager in order:
  // sin2th_12, sin2th_23, sin2th_13, delm2_12, delm2_23, delta_cp
  std::vector<double> oscpars{0.307,0.528,0.0218,7.53e-5, 2.509e-3,-1.601}; // Asimov A

  //  oscpars[1] = 0.3988888;

  std::cout<<"Using these oscillation parameters: ";
  for(unsigned ipar=0;ipar<oscpars.size();ipar++) std::cout<<" "<<oscpars.at(ipar);
  std::cout << std::endl;

  // Ask config file whether to use reactor constraint
  //bool useRC =  fitMan -> getRC() ;
  //std::cout << "use reactor prior is : " << useRC << std::endl ;
  //  osc->useReactorPrior(useRC); // this is hard coded inside, and is bad

  // Use prior for 12 parameters only
  //osc->setEvalLikelihood(0,false);
  osc->setFlipDeltaM23(true);
  osc->setEvalLikelihood(1,false);
  osc->setEvalLikelihood(2,false);
  //osc->setEvalLikelihood(3,false);
  osc->setEvalLikelihood(4,false);
  osc->setEvalLikelihood(5,false);
  osc->setParameters(oscpars);
  osc->acceptStep();
  //osc->setStepScale(fitMan->getOscStepScale());

  //  osc->setFlipDeltaM23(false);



  std::cout << "-----------------------------------------------------------" << std::endl;
  std::cout << "Creating Oscillator object" << std::endl;
  std::string OscillatorCfgName = fitMan->raw()["General"]["OscillatorConfigName"].as<std::string>();
  Oscillator* Oscill = new Oscillator(OscillatorCfgName);

  Oscill->FillOscillogram(osc->getPropPars(),25.0,0.5);
  
  std::vector<samplePDFDUNEBase*> pdfs;
  samplePDFDUNEBase *numu_pdf = new samplePDFDUNEBase(1.3628319e+23, "configs/AtmSample_numuselec.yaml", xsec);
  samplePDFDUNEBase *nue_pdf = new samplePDFDUNEBase(1.3628319e+23, "configs/AtmSample_nueselec.yaml", xsec);

  numu_pdf->SetOscillator(Oscill);
  nue_pdf->SetOscillator(Oscill);

  pdfs.push_back(numu_pdf);
  pdfs.push_back(nue_pdf);

  osc -> setParameters(oscpars);
  std::cout << "oscpars[0] = " << (osc -> getPropPars())[0] << std::endl
	    << "oscpars[1] = " << (osc -> getPropPars())[1] << std::endl
	    << "oscpars[2] = " << (osc -> getPropPars())[2] << std::endl
	    << "oscpars[3] = " << (osc -> getPropPars())[3] << std::endl
	    << "oscpars[4] = " << (osc -> getPropPars())[4] << std::endl
	    << "oscpars[5] = " << (osc -> getPropPars())[5] << std::endl;

  //###########################################################################################################
  // Set covariance objects equal to output of previous chain
  
  std::vector<double> oscparstarts;
  std::map<TString,std::vector<double> > parstarts;
  int lastStep = 0;

 

  if(fitMan->raw()["MCMC"]["StartFromPos"].as<bool>()) {//Start from values at the end of an already run chain
    //Read in paramter names and values from file
    std::cout << "MCMC getting starting position from " << fitMan->raw()["MCMC"]["PosFileName"].as<std::string>() << std::endl;
    TFile *infile = new TFile(fitMan->raw()["MCMC"]["PosFileName"].as<std::string>().c_str(), "READ");
    TTree *posts = (TTree*)infile->Get("posteriors");
    TObjArray* brlis = (TObjArray*)posts->GetListOfBranches();
    int nbr = brlis->GetEntries();
    TString branch_names[nbr];
    double branch_vals[nbr];
    int step_val;
    for (int i = 0; i < nbr; ++i) {
      TBranch *br = (TBranch*)brlis->At(i);
      TString bname = br->GetName();
      branch_names[i] = bname;
      if(branch_names[i] == "step") {
        posts->SetBranchAddress("step",&step_val);
        continue;
      }
      std::cout << " * Loading " << bname << std::endl;
      posts->SetBranchAddress(branch_names[i], &branch_vals[i]);
    }
    posts->GetEntry(posts->GetEntries()-1);
    std::map<TString,double> init_pars;
    for (int i = 0; i < nbr; ++i) {
      init_pars.insert( std::pair<TString, double>(branch_names[i], branch_vals[i]));
    }
    infile->Close();
    delete infile;
    
    //Make vectors of parameter value for each covariance
    std::vector<TString> covtypes;
    covtypes.push_back("xsec");
 
    for(unsigned icov=0;icov<covtypes.size();icov++){
      std::vector<double> covparstarts;
      std::map<TString, double>::const_iterator it;
      int iPar=0;
      while(it != init_pars.end()){
  	it = init_pars.find(covtypes[icov]+"_"+TString::Format("%d",iPar));
  	if (it != init_pars.end()) {
  	  covparstarts.push_back(it->second);
  	}
  	iPar++;
      }
      if(covparstarts.size()!=0) parstarts.insert(std::pair<TString,std::vector<double> >(covtypes[icov],covparstarts));
      else std::cout<<"Did not find any parameters in posterior tree for: "<<covtypes[icov]<<std::endl<<"assuming previous chain didn't use them"<<std::endl;
    }

    std::map<TString, double>::const_iterator itt;

    // set the oscillation parameters
    itt = init_pars.find("sin2th_12");
    oscparstarts.push_back(itt->second);
    itt = init_pars.find("sin2th_23");
    oscparstarts.push_back(itt->second);
    itt = init_pars.find("sin2th_13");
    oscparstarts.push_back(itt->second);
    itt = init_pars.find("delm2_12");
    oscparstarts.push_back(itt->second);
    itt = init_pars.find("delm2_23");
    oscparstarts.push_back(itt->second);
    itt = init_pars.find("delta_cp");
    oscparstarts.push_back(itt->second);

    lastStep = step_val;

  }

  // Set to nominal pars
  vector<double> xsecpar = xsec->getNominalArray();  
  xsec->setParameters(xsecpar);
  //  xsec->Print();

  numu_pdf->reweight(osc->getPropPars());
  TH1D *numu_asimov = (TH1D*)numu_pdf->get1DHist()->Clone("numu_asimov");
  nue_pdf->reweight(osc->getPropPars());
  TH1D *nue_asimov = (TH1D*)nue_pdf->get1DHist()->Clone("nue_asimov");

  TH2D *numu_asimov_2d = (TH2D*)numu_pdf->get2DHist()->Clone("numu_asimov_2d");
  TH2D *nue_asimov_2d = (TH2D*)nue_pdf->get2DHist()->Clone("nue_asimov_2d");

  // Print event rates to check
  std::cout << "-------- DUNE event rates for Asimov fit (Asimov fake data) ------------" << std::endl;
  std::cout << "mu like:   " << numu_asimov->Integral() << std::endl;
  std::cout << "e like:    " << nue_asimov->Integral() << std::endl;

  numu_pdf->addData(numu_asimov_2d); 
  nue_pdf->addData(nue_asimov_2d); 

	
  bool statsonly = fitMan->raw()["MCMC"]["StatOnly"].as<bool>(); 

    //###########################################################################################################

  // Back to actual nominal for fit
  // If starting from end values of a previous chain set everything to the values from there
  // and do acceptStep() to update fParCurr with these values
  //
  if(fitMan->raw()["MCMC"]["StartFromPos"].as<bool>()) {
    osc->setParameters(oscparstarts);
    osc->acceptStep();
    if(parstarts.find("xsec")!=parstarts.end()) {
      xsec->setParameters(parstarts["xsec"]);
      xsec->acceptStep();
    }
    //else {xsec->setParameters();}
  }
  else {
    //xsec->setParameters();
    //start chain from random position
    if(!statsonly)
    	xsec->throwParameters();
    osc->throwParameters();
  }


    //###########################################################################################################
  // MCMC

  mcmc *markovChain = new mcmc(fitMan);

  //numu_fds->Write();

  // set up
  int NSteps = fitMan->raw()["MCMC"]["NSTEPS"].as<int>(); 
  markovChain->setChainLength(NSteps);
  markovChain->addOscHandler(osc);
  if(lastStep > 0) markovChain->setInitialStepNumber(lastStep+1);

  // add samples
  markovChain->addSamplePDF(numu_pdf);
  markovChain->addSamplePDF(nue_pdf);


  // add systematic objects
  if (!statsonly) {
    markovChain->addSystObj(xsec);
  }
  else {std::cout << "STATS ONLY!" << std::endl;}
  
  // run!
  markovChain->runMCMC();

  return 0;
 }

