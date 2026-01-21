#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>

#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>

#include "Fitters/mcmc.h"
#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"

//CS YAML "Variations" node expects to be given each parameter to vary with :
// - "Name": name of the parameter as written in the CovObjs YAML file 
// - "OscParDefault" array: values that oscillation parameters should be given as a default for this specific parameter variation 
//    ex : put sterile mixing angles to non zero values to test variations of sterile cp phases
//    If "OscParDefault" not given, will use the "OscillationParameters" values specified in "General" node 
// - "VarValues" array: values we want the parameter to take

//TODO: Consider merging with SigmaVariations app at some point

int main(int argc, char * argv[]) {
  if(argc == 1){
    MACH3LOG_ERROR("Usage: bin/EventRatesDUNEBeam config.cfg");
    return 1;
  }
  manager* FitManager = new manager(argv[1]);

  
  //###############################################################################################################################
  //Create samplePDFFD objects
  
  ParameterHandlerGeneric* xsec = nullptr;

  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  std::vector<double> oscpars = FitManager->raw()["General"]["OscillationParameters"].as<std::vector<double>>();
  
  //###############################################################################################################################
  //Perform reweight and print total integral

  MACH3LOG_INFO("=======================================================");
  for(SampleHandlerFD* Sample: DUNEPdfs){
    Sample->Reweight();
    MACH3LOG_INFO("Event rate for {} : {:<5.2f}", Sample->GetTitle(), Sample->GetMCHist(Sample->GetNDim())->Integral());
  }
  
  //###############################################################################################################################
  //DB Can't use the core sigma variations as it's entirely set up around the concept of multiple selections per samplePDF object
  //   Thats not the case in the FD code, which has one selection per samplePDF object
  //   Consequently have to write out own code
  
  std::vector<ParameterHandlerBase*> CovObjs;
  CovObjs.emplace_back(xsec);

  MACH3LOG_INFO("=======================================================");

  std::string OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();
  TFile* File = TFile::Open(OutputFileName.c_str(),"RECREATE");

  for (ParameterHandlerBase* CovObj: CovObjs) {
    MACH3LOG_INFO("Starting Variations for covarianceBase Object: {}",CovObj->GetName());
    
    int nPars = CovObj->GetNumParams();

    for (int iPar=0;iPar<nPars;iPar++) {
      std::string ParName = CovObj->GetParName(iPar);

      for (auto const &param : FitManager->raw()["Variations"]) {

        std::string VarName = (param["Name"].as<std::string>());

        if(ParName == VarName) {

          MACH3LOG_INFO("\tParameter : {:<30}",ParName);

          if(!param["OscParDefault"]){ //if specific default values not specified for the parameter then use global default ones
            CovObj->SetParameters(oscpars);
          }
          else {
            CovObj->SetParameters((param["OscParDefault"].as<std::vector<double>>()));
          }

          std::vector<double> valVariations = (param["VarValues"].as<std::vector<double>>());

          File->cd();
          File->mkdir(ParName.c_str());
          File->cd(ParName.c_str());
      
          for (size_t iSigVar=0;iSigVar<valVariations.size();iSigVar++) {
	    double VarVal = valVariations[iSigVar];
	    
	    MACH3LOG_INFO("\t\tParameter Value : {:<10.7f}",VarVal);
            CovObj->SetParProp(iPar,VarVal);
	    
	    for (size_t iSample=0;iSample<DUNEPdfs.size();iSample++) {
	      std::string SampleName = DUNEPdfs[iSample]->GetTitle();
	      
	      File->cd(ParName.c_str());
	      if (iSigVar == 0) {
		File->mkdir((ParName+"/"+SampleName).c_str());
	      }
	      File->cd((ParName+"/"+SampleName).c_str());
	      
	      DUNEPdfs[iSample]->Reweight();
              TH1* Hist = DUNEPdfs[iSample]->GetMCHist(DUNEPdfs[iSample]->GetNDim());
	      MACH3LOG_INFO("\t\t\tSample : {:<30} - Integral : {:<10}",SampleName,Hist->Integral());
	      Hist->Write(Form("Variation_%.2e",VarVal));
	    }
          }

        }
      }
    }

    MACH3LOG_INFO("=======================================================");
  }
  

  File->Close();
  
}
