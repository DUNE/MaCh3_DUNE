#include <ctime>
#include <iostream>
#include <vector>
#include <sys/time.h>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TStopwatch.h>

#include "samplePDFDUNE/samplePDFDUNEBase.h"
#include "samplePDFDUNE/samplePDFDUNEBaseND.h"
#include "mcmc/mcmc.h"

int main(int argc, char **argv)
{

  #ifdef MULTITHREAD
    std::cout << "MUTLI THREADING IS ON" << std::endl;
  #else
    std::cout << "MUTLI THREADING IS OFF" << std::endl;
  #endif

  #ifdef CPU_ONLY
    std::cout << "ONLY CPU OSC CALC" << std::endl;
  #else
    std::cout << "GPU OSC CALC IS ON" << std::endl;
  #endif

  // Read config
  manager *fitMan = new manager(argv[1]);

  std::string  XsecMatrixFile = fitMan->raw()["General"]["Systematics"]["XsecCovFile"].as<std::string>(); 
  std::string  XsecMatrixName = fitMan->raw()["General"]["Systematics"]["XsecCovName"].as<std::string>();
  std::string  OscMatrixFile = fitMan->raw()["General"]["Systematics"]["OscCovFile"].as<std::string>(); 
  std::string  OscMatrixName = fitMan->raw()["General"]["Systematics"]["OscCovName"].as<std::string>();
  double FDPOT = fitMan->raw()["General"]["FDPOT"].as<double>(); 
  double NDPOT = fitMan->raw()["General"]["NDPOT"].as<double>(); 
  
  //###########################################################################################################
  // Covariance Objects

  covarianceXsec *xsec;
  covarianceOsc *osc;

  xsec = new covarianceXsec(XsecMatrixName.c_str(), XsecMatrixFile.c_str());


  // Setting flat priors based on XSECPARAMFLAT list in configuration file 

  bool statonly = fitMan->raw()["MCMC"]["StatOnly"].as<bool>();;

  std::vector<double> oscpars = fitMan->raw()["General"]["OscParams"].as<std::vector<double>>();

  osc = new covarianceOsc(OscMatrixName.c_str(), OscMatrixFile.c_str());

  if (!(oscpars.size()==6 || oscpars.size()==7))
    {
      std::cout<<"Input osc pars not of right size, there should be six entries (or seven if setting beta)"<<std::endl;
      std::cout<<"oscpars.size() = " << oscpars.size() << std::endl;
      exit(1);
    }

  std::cout<<"Using these oscillation parameters: ";
  for(unsigned ipar=0;ipar<oscpars.size();ipar++)
    std::cout<<" "<<oscpars.at(ipar);
  std::cout << std::endl;

  osc->setFlipDeltaM23(true);

  // Use prior for 12 parameters only
  //osc->setEvalLikelihood(0,false);
  osc->setEvalLikelihood(1,false);
  osc->setEvalLikelihood(2,false);
  //osc->setEvalLikelihood(3,false);
  osc->setEvalLikelihood(4,false);
  osc->setEvalLikelihood(5,false);

  osc->setParameters(oscpars);
  osc->acceptStep();


  xsec->setStepScale(fitMan->raw()["General"]["Systematics"]["XsecStepScale"].as<double>());
  osc->setStepScale(fitMan->raw()["General"]["Systematics"]["OscStepScale"].as<double>());

  bool addFD = fitMan->raw()["General"]["IncludeFD"].as<bool>();
  bool addND = fitMan->raw()["General"]["IncludeND"].as<bool>();

  if (!addFD && !addND) {std::cerr << "[ERROR:] You've chosen NOT to include FD or ND samples in the config file... you need to add something!" << std::endl; throw;}


  std::vector<samplePDFFDBase*> SamplePDFs;
  
  if(addFD) { 
    samplePDFDUNEBase *numu_pdf = new samplePDFDUNEBase(FDPOT, "configs/SamplePDFDune_FHC_numuselec.yaml", xsec);
    SamplePDFs.push_back(numu_pdf);
    samplePDFDUNEBase *nue_pdf = new samplePDFDUNEBase(FDPOT, "configs/SamplePDFDune_FHC_nueselec.yaml", xsec);
    SamplePDFs.push_back(nue_pdf);
    samplePDFDUNEBase *numubar_pdf = new samplePDFDUNEBase(FDPOT, "configs/SamplePDFDune_RHC_numuselec.yaml", xsec);
    SamplePDFs.push_back(numubar_pdf);
    samplePDFDUNEBase *nuebar_pdf = new samplePDFDUNEBase(FDPOT, "configs/SamplePDFDune_RHC_nueselec.yaml", xsec);
    SamplePDFs.push_back(nuebar_pdf);
  }
 
  if(addND) {
    samplePDFDUNEBaseND * numu_cc_nd_pdf = new samplePDFDUNEBaseND(NDPOT, "configs/SamplePDFDuneND_FHC_CCnumuselec.yaml", xsec);
    SamplePDFs.push_back(numu_cc_nd_pdf);
    samplePDFDUNEBaseND * numubar_cc_nd_pdf = new samplePDFDUNEBaseND(NDPOT, "configs/SamplePDFDuneND_RHC_CCnumuselec.yaml", xsec);
    SamplePDFs.push_back(numubar_cc_nd_pdf);
  }


  for(unsigned sample_i=0; sample_i < SamplePDFs.size() ; ++sample_i) {
    SamplePDFs[sample_i] -> useNonDoubledAngles(true);
    SamplePDFs[sample_i] -> SetupOscCalc(osc->GetPathLength(), osc->GetDensity());
  }


  //###########################################################################################################
  // Set covariance objects equal to output of previous chain
  
  std::vector<double> oscparstarts;
  std::map<TString,std::vector<double> > parstarts;
  int lastStep = 0;

 
  //Start from values at the end of an already run chain
  if(fitMan->raw()["MCMC"]["StartFromPos"].as<bool>()) {
    //Read in paramter names and values from file
    std::cout << "MCMC getting starting position from " << fitMan->raw()["MCMC"]["PosFileName"].as<std::string>() << std::endl;
    TFile *infile = TFile::Open(fitMan->raw()["MCMC"]["PosFileName"].as<std::string>().c_str(), "READ");
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
  
  // set to nominal
  xsec->setParameters();

  //###########################################################################################################
  // Apply reweight and find event rates

  std::cout << "Nominal Oscillated Event rates: " << std::endl;

  for(unsigned sample_i=0; sample_i < SamplePDFs.size() ; ++sample_i) {
    SamplePDFs[sample_i] -> reweight(osc->getPropPars()); // Get Nominal Predictions 

	std::string name = SamplePDFs[sample_i] -> GetSampleName();
    TString NameTString = TString(name.c_str());

	if (SamplePDFs[sample_i] -> GetBinningOpt() == 1) {
      TH1D *Nominal_1D = (TH1D*)SamplePDFs[sample_i]->get1DHist()->Clone(NameTString+"_nominal");
      std::cout << name.c_str() << ": " << Nominal_1D->Integral() << std::endl;
	}
      
	else if (SamplePDFs[sample_i] -> GetBinningOpt() == 2) {
      TH2D *Nominal_2D = (TH2D*)SamplePDFs[sample_i]->get2DHist()->Clone(NameTString+"_nominal");
      std::cout << name.c_str() << ": " << Nominal_2D->Integral() << std::endl;
	}

	else {
	  std::cerr << "[ERROR:] Unknown binning option for sample, don't know whether to make 1D or 2D hist! : " << name.c_str() << " : " << SamplePDFs[sample_i] -> GetBinningOpt() << std::endl; throw;
	}

  }


  std::cout << "-----------------------" << std::endl; 

  // Set to Asimov values
  // If wanting to do an Asimov fit to non-nominal values of xsec params set them here!
  // e.g. xsecpar[i] = ....

  xsec->setParameters();


  // Print Asimov event rates to check
  std::cout << "-------- Event rates for Asimov Data  ------------" << std::endl;

  for(unsigned sample_i=0; sample_i < SamplePDFs.size() ; ++sample_i) {
    SamplePDFs[sample_i] -> reweight(osc->getPropPars()); // Get Nominal Predictions 

	std::string name = SamplePDFs[sample_i] -> GetSampleName();
    TString NameTString = TString(name.c_str());

	if (SamplePDFs[sample_i] -> GetBinningOpt() == 1) {
      TH1D *Asimov_1D = (TH1D*)SamplePDFs[sample_i]->get1DHist()->Clone(NameTString+"_asimov");
      std::cout << name.c_str() << ": " << Asimov_1D->Integral() << std::endl;
      SamplePDFs[sample_i] -> addData(Asimov_1D); 
	}
      
	if (SamplePDFs[sample_i] -> GetBinningOpt() == 2) {
      TH2D *Asimov_2D = (TH2D*)SamplePDFs[sample_i]->get2DHist()->Clone(NameTString+"_asimov");
      std::cout << name.c_str() << ": " << Asimov_2D->Integral() << std::endl;
      SamplePDFs[sample_i] -> addData(Asimov_2D); 
	}

  }

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
    else {xsec->setParameters();}
  }
  else {
    xsec->setParameters();
  }

  //###########################################################################################################
  // MCMC

  mcmc *markovChain = new mcmc(fitMan);

  // set up
  int NSteps = fitMan->raw()["MCMC"]["NSTEPS"].as<int>(); 
  markovChain->setChainLength(NSteps);
  markovChain->addOscHandler(osc);
  if(lastStep > 0) markovChain->setInitialStepNumber(lastStep+1);

  // add samples

  for(unsigned sample_i=0; sample_i < SamplePDFs.size() ; ++sample_i) {
    markovChain->addSamplePDF(SamplePDFs[sample_i]);
  }

  if(!fitMan->raw()["MCMC"]["StartFromPos"].as<bool>()) {
    //start chain from random position unless starting from a previous chain
    xsec->proposeStep();
    osc->proposeStep();
  }

  // add systematic objects
  if (!statonly) {
    markovChain->addSystObj(xsec);
  }
  else {std::cout << "STAT ONLY!" << std::endl;}
  
  // run!
  markovChain->runMCMC();

  return 0;
}
