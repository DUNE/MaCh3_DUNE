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
#include "samplePDFDUNE/samplePDFDUNEBaseND.h"

#include "mcmc/mcmc.h"


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


  // Read config
  manager *fitMan = new manager(argv[1]);

  std::string  XsecMatrixFile = fitMan->raw()["General"]["Systematics"]["XsecCovFile"].as<std::string>(); 
  std::string  XsecMatrixName = fitMan->raw()["General"]["Systematics"]["XsecCovName"].as<std::string>();
  std::string  OscMatrixFile = fitMan->raw()["General"]["Systematics"]["OscCovFile"].as<std::string>(); 
  std::string  OscMatrixName = fitMan->raw()["General"]["Systematics"]["OscCovName"].as<std::string>();
  double FDPOT = fitMan->raw()["General"]["FDPOT"].as<double>(); 
  double NDPOT = fitMan->raw()["General"]["NDPOT"].as<double>(); 
 
  std::string OutFileName = fitMan->raw()["General"]["Output"]["FileName"].as<std::string>();
  TFile *OutFile = TFile::Open(OutFileName.c_str() , "RECREATE");
  
  //###########################################################################################################
  // Covariance Objects

  covarianceXsec *xsec;
  covarianceOsc *osc;

  xsec = new covarianceXsec(XsecMatrixName.c_str(), XsecMatrixFile.c_str());


  // Setting flat priors based on XSECPARAMFLAT list in configuration file 

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

  vector<double> xsecpar = xsec->getNominalArray();  
  // Set to Asimov values
  // If wanting to do an Asimov fit to non-nominal values of xsec params set them here!
  // e.g. xsecpar[i] = ....

  xsec->setParameters();


  // Print Asimov event rates to check
  std::cout << "-------- Event rates for Asimov Data  ------------" << std::endl;
  
  std::vector<std::string> sample_names;
  for(unsigned sample_i=0; sample_i < SamplePDFs.size() ; ++sample_i) {
    SamplePDFs[sample_i] -> useNonDoubledAngles(true);
    SamplePDFs[sample_i] -> SetupOscCalc(osc->GetPathLength(), osc->GetDensity());

    SamplePDFs[sample_i] -> reweight(osc->getPropPars()); // Get Nominal Predictions 

	std::string name = SamplePDFs[sample_i] -> GetSampleName();
    TString NameTString = TString(name.c_str());
	sample_names.push_back(name);

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

  TCanvas *sys_canv = new TCanvas("sys_canv","",1200,600);
  
  // Save pars before modifying
  std::vector<double> tpar = xsecpar;

  unsigned int n_points = fitMan->raw()["LLHScans"]["ParPoints"].as<int>();
  int n_sigma = 3;

  // Loop over parameters
  for (unsigned  par=0; par < xsecpar.size(); par++) {

    std::vector<std::string> histnames;
    std::vector<std::string> histtitles;
    std::vector<TH1D*> llh_hists;

    xsecpar = xsec->getNominalArray();
    double lowerb = xsec->getNominal(par) - (n_sigma+0.01) * sqrt((*xsec->getCovMatrix())(par,par));
    double upperb = xsec->getNominal(par) + (n_sigma+0.01) * sqrt((*xsec->getCovMatrix())(par,par));
 
	 // Create Histogram for each sample + systematic penalty + total sample	
    for(unsigned sample_i=0; sample_i < SamplePDFs.size() ; ++sample_i) {
	  std::string histname = "xsec_" + std::to_string(par) + "_llh_" + sample_names[sample_i];
	  std::string histtitle = "xsec_" + std::to_string(par) + "_" + sample_names[sample_i];
      TH1D *hScan = new TH1D(histname.c_str(), histtitle.c_str(), n_points, lowerb, upperb);
	  llh_hists.push_back(hScan);
	}
	
    std::string histname_systpen = "xsec_" + std::to_string(par) + "_llh_syst";
    std::string histtitle_systpen = "xsec_" + std::to_string(par) + "_syst";
    TH1D *hScan_pen = new TH1D(histname_systpen.c_str(), histtitle_systpen.c_str(), n_points, lowerb, upperb);
	llh_hists.push_back(hScan_pen);

    std::string histname_total = "xsec_" + std::to_string(par) + "_llh_total_sample";
    std::string histtitle_total = "xsec_" + std::to_string(par) + "_total_sample";
    TH1D *hScan_total = new TH1D(histname_total.c_str(), histtitle_total.c_str(), n_points, lowerb, upperb);
	llh_hists.push_back(hScan_total);
	
    xsecpar[par] = xsec->getNominal(par)-n_sigma*sqrt((*xsec->getCovMatrix())(par,par));
 
	// Increment in sigma units
    double dsigma = (2*n_sigma)/((double)n_points-1);
	std::cout << "dsigma = " << dsigma << std::endl;

	// Loop over parameter values
    for (unsigned val =0; val < n_points; val++) {
      double totalllh = 0;
      
	  xsec->setParameters(xsecpar);
	  std::cout << "Val = " << xsecpar[par] << std::endl;

	  // Calc LLH for each sample
      for(unsigned sample_i=0; sample_i < SamplePDFs.size() ; ++sample_i) {
        SamplePDFs[sample_i]  -> reweight(osc -> getPropPars());
	    llh_hists[sample_i]->Fill(xsecpar[par], 2 * SamplePDFs[sample_i] -> GetLikelihood());	
        totalllh += 2 * SamplePDFs[sample_i] -> GetLikelihood();
		std::cout << "LLH for sample " << sample_i << " = " << 2 * SamplePDFs[sample_i] -> GetLikelihood() << std::endl;
      } //end of sample loop

	  // Calc Penalty LLH
	  llh_hists[SamplePDFs.size()]->Fill(xsecpar[par], 2 * (xsec->GetLikelihood()));
	  // Total sample LLH
	  llh_hists[SamplePDFs.size()+1]->Fill(xsecpar[par], 2 * totalllh);

      xsecpar[par] += dsigma*sqrt((*xsec->getCovMatrix())(par,par)); // Increment parameter
    } //end of value loop

    OutFile->cd();

	for(unsigned hist=0; hist < llh_hists.size(); hist++)
	{
      llh_hists[hist]->Write();
	}

    std::cout << "Finished xsec param " << par << std::endl;
  } //end of parameter loop

  std::cout << "Finished LLH Scans!" << std::endl;
  return 0;
}


