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

#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Fitters/MaCh3Factory.h"

int main(int argc, char * argv[]) {
  auto FitManager = MaCh3ManagerFactory(argc, argv);

  //###############################################################################################################################

  //DB Sigma variations in units of each parameters Sigma
  std::vector<double> sigmaVariations = {-3, -1, 0, 1, 3};

  //###############################################################################################################################
  //Create samplePDFFD objects
  
  ParameterHandlerGeneric* xsec = nullptr;

  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  //###############################################################################################################################
  //Perform reweight and print total integral

    TFile* Osc = TFile::Open("EventRates/FHCsTwentyIndChanOsc.root"); // Loading in histograms
  TFile* Unosc = TFile::Open("EventRates/FHCsTwentyIndChanUnosc.root");

  TIter next(Osc->GetListOfKeys()); // Get list of different items within oscillated data histograms (48 in total, 12 channels in 4 samples)
  TKey* key;  // Initialise

  std::vector<std::string> KeyNames;
  std::vector<int> SelectedKeys = {0, 13, 2, 15, 16, 17, 6, 7, 20, 9, 22, 11};
  int ParIndex = 0.0;
  //int KeyIndex = 12.0;

  for(int i = 0; i < xsec->GetNumParams(); i++){ // For every param in the xsec group
    if(xsec->IsParFromGroup(i, "EParam")){ // If param is from our energy normalisation parameter group
      ParIndex = i; // Set the index of first parameter
      break;
    }
  }

  while((key = (TKey*)next())) {
    KeyNames.push_back(key->GetName()); 
  }

  //for (int index = 0; index < 12; index++) {
  for (int index : SelectedKeys) {
    const std::string& name = KeyNames[index];
    auto HistoOsc = Osc->Get<TH1D>(name.c_str()); // Getting names of histograms
    auto HistoUnosc = Unosc->Get<TH1D>(name.c_str()); 
    int NumBins = HistoOsc->GetNbinsX();
    for(int j = 1; j <= NumBins; j++) { // For each energy bin, starting from 1 to avoid the overflow bin
      double BinSizeOsc = HistoOsc->GetBinContent(j); // Get the number of events in specific energy bin
      double BinSizeUnosc = HistoUnosc->GetBinContent(j);
      double Param = (BinSizeUnosc == 0.0 ? 0.0 : BinSizeOsc/BinSizeUnosc);
      xsec->SetPar(ParIndex, Param);
      if (BinSizeUnosc == 0.0) xsec->ToggleFixParameter(ParIndex);
      ParIndex++; // Increment parameter index to keep amending in sequence
    }
  }

  MACH3LOG_INFO("=======================================================");
  for(SampleHandlerFD* handler: DUNEPdfs){
    handler->Reweight();
    for (int iSample=0;iSample<handler->GetNsamples();iSample++) {
      MACH3LOG_INFO("Event rate for {} : {:<5.2f}", handler->GetSampleTitle(iSample), handler->GetMCHist(iSample)->Integral());
    }
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
      std::string ParName = CovObj->GetParFancyName(iPar);
      double VarInit = CovObj->GetParInit(iPar);
      double VarSigma = CovObj->GetDiagonalError(iPar);
      
      MACH3LOG_INFO("\tParameter : {:<30} - Variations around value : {:<10.7f} , in units of 1 Sigma : {:<10.7f}",ParName,VarInit,VarSigma);

      File->cd();
      File->mkdir(ParName.c_str());
      File->cd(ParName.c_str());
      
      for (size_t iSigVar=0;iSigVar<sigmaVariations.size();iSigVar++) {
        double VarVal = VarInit + sigmaVariations[iSigVar] * VarSigma;
        if (VarVal < CovObj->GetLowerBound(iPar)) VarVal = CovObj->GetLowerBound(iPar);
        if (VarVal > CovObj->GetUpperBound(iPar)) VarVal = CovObj->GetUpperBound(iPar);

        MACH3LOG_INFO("\t\tVariation {:<5.3f} - Parameter Value : {:<10.7f}",
                      sigmaVariations[iSigVar], VarVal);
        CovObj->SetParProp(iPar, VarVal);

        for (auto handler : DUNEPdfs) {
          for (int iSample = 0; iSample < handler->GetNsamples(); iSample++) {
            std::string SampleName = handler->GetSampleTitle(iSample);

            File->cd(ParName.c_str());
            if (iSigVar == 0) {
              File->mkdir((ParName + "/" + SampleName).c_str());
            }
            File->cd((ParName + "/" + SampleName).c_str());

            handler->Reweight();
            TH1 *Hist = handler->GetMCHist(iSample);
            MACH3LOG_INFO("\t\t\tSample : {:<30} - Integral : {:<10}", SampleName,
                          Hist->Integral());

            Hist->Write(Form("Variation_%i", (int)iSigVar));
          }
        }
      }

      CovObj->SetParProp(iPar, VarInit);
    }

    MACH3LOG_INFO("=======================================================");
  }

  File->Close();
}
