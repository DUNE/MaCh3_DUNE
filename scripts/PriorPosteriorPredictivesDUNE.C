#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>

#include "samplePDFDUNE/samplePDFDUNEBase.h"
#include "samplePDFDUNE/samplePDFDUNEBaseND.h"
#include "manager/manager.h"


std::string getNameNoExt(std::string name, std::string ext)  
{                                                            
  std::size_t pos ;                                          
  pos = name.find(ext);                                      
  name = name.substr(0,pos);                                 
  return name ;                                              
}                                                            

void saveCanvas(TCanvas* canvas, std::string name, std::string legend)                                                                  
{                                                            
  name = getNameNoExt(name, ".root") ;                       
  name = name + legend + ".root" ;                     
  canvas -> SaveAs(name.c_str()) ;                           

  name = getNameNoExt(name, ".root") ;                       
  name = name + ".png" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           

  name = getNameNoExt(name, ".png") ;                        
  name = name + ".pdf" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           

  name = getNameNoExt(name, ".pdf") ;                        
  name = name + ".eps" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           

} 


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

  double FDPOT = fitMan->raw()["General"]["FDPOT"].as<double>(); 
  double NDPOT = fitMan->raw()["General"]["NDPOT"].as<double>(); 

  // ----------------------- COVARIANCE AND SAMPLEPDF OBJECTS ---------------------------------------- //

  gStyle -> SetPalette(1);

  // make file to save the variation of xsec params to
  std::string OutputFileName = fitMan->raw()["General"]["Output"]["FileName"].as<std::string>();
  TFile *OutputFile = new TFile(OutputFileName.c_str() , "RECREATE");

  //initialise xsec
  covarianceXsec *xsec = new covarianceXsec(XsecMatrixName.c_str(), XsecMatrixFile.c_str());

  std::cout << "---------- Printing off nominal parameter values ----------" << std::endl;
  std::cout << "Cross section parameters:" << std::endl;
  xsec->printNominal();

  bool XsecParsAtGen = fitMan->raw()["General"]["Systematics"]["XsecAtGen"].as<bool>();

  std::cout << "---------- Finished printing nominal parameter values ----------" << std::endl;

  covarianceOsc *osc = new covarianceOsc(OscMatrixName.c_str(), OscMatrixFile.c_str());

  std::vector<double> oscpars = fitMan->raw()["General"]["OscParams"].as<std::vector<double>>();

  std::cout<<"Using these oscillation parameters: ";
  for(unsigned ipar=0;ipar<oscpars.size();ipar++) std::cout<<" "<<oscpars.at(ipar);
  std::cout << std::endl;
  //osc->setFlipDeltaM23(true);

  // Ask config file whether to use reactor constraint
  //std::cout << "use reactor prior is : " << useRC << std::endl ;

  // Use prior for 12 parameters only
  //osc->setEvalLikelihood(0,false);
  ///osc->setEvalLikelihood(1,false);
  //osc->setEvalLikelihood(2,false);
  //osc->setEvalLikelihood(3,false);
  //osc->setEvalLikelihood(4,false);
  //osc->setEvalLikelihood(5,false);
  // This line gives a crash and stack trace...
  //osc->setParameters(oscpars);
  //osc->acceptStep();

  //bool addFD = fitMan->raw()["General"]["IncludeFD"].as<bool>();
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
    samplePDFDUNEBaseND * FHC_numuCCND_pdf = new samplePDFDUNEBaseND(NDPOT, "configs/SamplePDFDuneND_FHC_CCnumuselec.yaml", xsec);
    SamplePDFs.push_back(FHC_numuCCND_pdf);
    samplePDFDUNEBaseND * RHC_numuCCND_pdf = new samplePDFDUNEBaseND(NDPOT, "configs/SamplePDFDuneND_RHC_CCnumuselec.yaml", xsec);
    SamplePDFs.push_back(RHC_numuCCND_pdf);
  }

  // Oscillated
  osc -> setParameters(oscpars);
  std::cout << "oscpars[0] = " << (osc -> getPropPars())[0] << std::endl
	<< "oscpars[1] = " << (osc -> getPropPars())[1] << std::endl
	<< "oscpars[2] = " << (osc -> getPropPars())[2] << std::endl
	<< "oscpars[3] = " << (osc -> getPropPars())[3] << std::endl
	<< "oscpars[4] = " << (osc -> getPropPars())[4] << std::endl
	<< "oscpars[5] = " << (osc -> getPropPars())[5] << std::endl;

  // unosc
  std::vector<double> oscpars_un(oscpars);
  oscpars_un[0] = 0;
  oscpars_un[1] = 0;
  oscpars_un[2] = 0;
  //oscpars_un[3] = 0;
  //oscpars_un[4] = 0;

  //Setup the cross-section parameters
  //This should get the prior values.
  std::vector<double> XsecParVals = xsec->getNominalArray();

  if(XsecParsAtGen){
	TFile* XsecFile = new TFile(XsecMatrixFile.c_str(), "READ");
	TVectorD* XsecGeneratedParamArray = (TVectorD*)XsecFile->Get("xsec_param_nom");
	std::cout << "Setting xsec systs to their generated values " << std::endl;
	for (unsigned param_i = 0 ; param_i < XsecParVals.size() ; ++param_i) {
	  std::cout << "Generated value for param " << param_i << " is " << (*XsecGeneratedParamArray)(param_i) << std::endl;
	  XsecParVals[param_i] = (*XsecGeneratedParamArray)(param_i);
	  std::cout << "Set parameter " << param_i << " to value " << XsecParVals[param_i] << std::endl;
	}
  }
  else{
	std::cout << "Keeping xsec parameters at their prior values" << std::endl;
  }

  xsec->setParameters(XsecParVals);
  xsec->setStepScale(fitMan->raw()["General"]["Systematics"]["XsecStepScale"].as<double>());
  //osc->setStepScale(fitMan->raw()["General"]["Systematics"]["OscStepScale"].as<double>());


  //Some place to store the histograms
  std::vector<TH1D*> oscillated_hists;
  std::vector<TH1D*> unoscillated_hists;
  std::vector<std::string> sample_names;

  for (unsigned sample_i = 0 ; sample_i < SamplePDFs.size() ; ++sample_i) {

	std::string name = SamplePDFs[sample_i]->GetSampleName();

	sample_names.push_back(name);
	TString NameTString = TString(name.c_str());
	// Unoscillated
	osc -> setParameters(oscpars_un);
	osc -> acceptStep();

	//SamplePDFs[sample_i] -> SetupOscCalc(osc->GetPathLength(), osc->GetDensity());
	//SamplePDFs[sample_i] -> reweight(osc->getPropPars());
	TH1D *sample_unosc = (TH1D*)SamplePDFs[sample_i] -> get1DHist() -> Clone(NameTString+"_unosc");
	unoscillated_hists.push_back(sample_unosc);

	OutputFile -> cd();
	sample_unosc			-> Write(NameTString+"_unosc");

	osc -> setParameters(oscpars);
	osc -> acceptStep();

	//SamplePDFs[sample_i] -> reweight(osc->getPropPars());
	TH1D *sample_osc = (TH1D*)SamplePDFs[sample_i] -> get1DHist()->Clone(NameTString+"_osc");
	oscillated_hists.push_back(sample_osc);

	OutputFile -> cd();
	sample_osc			-> Write(NameTString+"_osc");

  }

  //Now print out some event rates, we'll make a nice latex table at some point 
  std::cout << "Integrals of nominal hists: " << std::endl;
  for (unsigned sample_i = 0; sample_i < SamplePDFs.size() ; ++sample_i) {
	std::cout << sample_names[sample_i].c_str() << " unosc:      " << unoscillated_hists[sample_i]-> Integral() << std::endl;
	std::cout << sample_names[sample_i].c_str() << "   osc:      " << oscillated_hists[sample_i]-> Integral() << std::endl; 
	std::cout << " " << std::endl;
  }
  std::cout << "~~~~~~~~~~~~~~~~" << std::endl;

  /////////////////
  //
  // This is where the Posterior/Prior throwing happens
  //
  //////////////////

  TRandom3* rnd = new TRandom3(100);
  unsigned int NToys = fitMan->raw()["PostPriorPredict"]["NToys"].as<int>();

  bool OnlyPrior = fitMan->raw()["PostPriorPredict"]["OnlyPrior"].as<int>();

  if(OnlyPrior) 
    {std::cout << "You've chosen to produce a Prior Predicitive Distribution" << std::endl}

  else 
    {std::cout << "You've chosen to produce a Posterior Predicitive Distribution" << std::endl}

  //if Posterior predictive read in MCMC chain specified in the exec config 
  if(!OnlyPrior) {




  TH1D* Hist1D;
  for(unsigned int toy=0; toy < NToys; toy++) {

    if( toy % (NToys/10) == 0){cout << "Throw " << toy << " / " << NToys << endl;}


    // For Prior Predictives throw parameters according to their prior uncertainty
    if (OnlyPrior) {
	  xsec->throwParameters();
	}

	  osc -> setParameters(oscpars);
	  osc -> acceptStep();
	
      for (unsigned sample_i = 0; sample_i < SamplePDFs.size() ; ++sample_i) {
	    //SamplePDFs[sample_i] -> reweight(osc->getPropPars());
		Hist1D =  (TH1D*)SamplePDFs[sample_i] -> get1DHist();
		TString WriteName = SamplePDFs[sample_i]->GetSampleName()+"_";
		WriteName +=toy;
        Hist1D->Write(WriteName);
	  }
  }

  OutputFile->Close();
  return 0;
}