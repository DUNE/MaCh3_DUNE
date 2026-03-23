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
#include <TFile.h>
#include <TKey.h>

#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Fitters/MaCh3Factory.h"
#include "Fitters/PredictiveThrower.h"

int main(int argc, char * argv[]) {

  // if(argc == 1){
  //   MACH3LOG_ERROR("Usage: bin/EventRatesDUNEBeam config.cfg");
  //   return 1;
  // }

  auto FitManager = MaCh3ManagerFactory(argc, argv);

  // #########################################################################
  // Produce sample PDF objects

  ParameterHandlerGeneric* xsec = nullptr;

  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  // #########################################################################
  // Load in and create all relevant energy parameters

  TFile* Osc = TFile::Open("TrueIndChanOsc.root");
  TFile* Unosc = TFile::Open("TrueIndChanUnosc.root");

  TIter next(Osc->GetListOfKeys()); 
  TKey* key; 
  int ParIndex = 0.0;

  for(int i = 0; i < xsec->GetNumParams(); i++){
    if( xsec->IsParFromGroup(i, "EParam")){
      ParIndex = i;
      break;
    }
  }

  while ((key = (TKey*)next())) { 

    auto HistoOsc = Osc->Get<TH1D>(key->GetName());
    auto HistoUnosc = Unosc->Get<TH1D>(key->GetName());

    int NumBins = HistoOsc->GetNbinsX();

    for(int j = 1; j <= NumBins; j++) {
      double BinSizeOsc = HistoOsc->GetBinContent(j);
      double BinSizeUnosc = HistoUnosc->GetBinContent(j);

      double Param;
      if(BinSizeUnosc == 0) {
        Param = 0;
        xsec->SetPar(ParIndex, Param);
      }
      else {
        Param = BinSizeOsc / BinSizeUnosc; 
        xsec->SetPar(ParIndex, Param);
      }
        
      ParIndex++;
    }
  }

  for(int k = 0; k < xsec->GetNumParams(); k++){
    if((xsec->GetParProp(k) == 0) && xsec->IsParFromGroup(k, "EParam")){
      xsec->ToggleFixParameter(k);
    }
  }

  // #########################################################################
  // Perform reweight, print total integral and set the data

  std::vector<TH1*> DUNEHists;
  for(auto Sample : DUNEPdfs){
    Sample->Reweight();
    DUNEHists.push_back(Sample->GetMCHist(Sample->GetNDim()));
    MACH3LOG_INFO("Event rate for {} : {:<5.2f}", Sample->GetTitle(), Sample->GetMCHist(Sample->GetNDim())->Integral());

    if (Sample->GetNDim() == 1) {
      Sample->AddData((TH1D*)DUNEHists.back());
    } else if (Sample->GetNDim() == 2) {
      Sample->AddData((TH2D*)DUNEHists.back());
    }
  }

  std::unique_ptr<PredictiveThrower> MaCh3Fitter = std::make_unique<PredictiveThrower>(FitManager.get());

  // #########################################################################
  // Add samples and run the Posterior Predictive

  MaCh3Fitter->AddSystObj(xsec);

  for(auto Sample : DUNEPdfs){
    MaCh3Fitter->AddSampleHandler(Sample);
  }

  MaCh3Fitter->ProduceToys();
  MaCh3Fitter->RunPredictiveAnalysis();

}