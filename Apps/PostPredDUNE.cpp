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

  auto FitManager = MaCh3ManagerFactory(argc, argv);

  // #########################################################################
  // Produce sample PDF objects

  ParameterHandlerGeneric* xsec = nullptr;

  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  // #########################################################################
  // Load in and create all relevant energy parameters

  TFile* Osc = TFile::Open("EventRates/CCIndChanOsc.root"); // Loading in histograms
  TFile* Unosc = TFile::Open("EventRates/CCIndChanUnosc.root");

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

  // for (int i = 0; i < 12; i++) {
  //   key = (TKey*) next();
  // }

  // while ((key = (TKey*)next())) { // Go through all keys in sequence
  //   if(KeyIndex == 24.0) break;
  //   KeyIndex++;
  //   auto HistoOsc = Osc->Get<TH1D>(key->GetName()); // Getting names of histograms
  //   auto HistoUnosc = Unosc->Get<TH1D>(key->GetName()); 
  //   int NumBins = HistoOsc->GetNbinsX(); // Find number of bins (number of energy normalisation parameters for this histogram)
  //   for(int j = 1; j <= NumBins; j++) { // For each energy bin, starting from 1 to avoid the overflow bin
  //     double BinSizeOsc = HistoOsc->GetBinContent(j); // Get the number of events in specific energy bin
  //     double BinSizeUnosc = HistoUnosc->GetBinContent(j);
  //     double Param;
  //     if(BinSizeUnosc == 0) { // If no unoscillated data, set param to 0
  //       Param = 0;
  //       xsec->SetPar(ParIndex, Param);
  //       xsec->ToggleFixParameter(ParIndex); // Fix these params at 0
  //     }
  //     else { // If unoscillated data, calculate ratio between these as needed to induce oscillation
  //       Param = BinSizeOsc / BinSizeUnosc; 
  //       xsec->SetPar(ParIndex, Param);
  //     } 
  //     ParIndex++; // Increment parameter index to keep amending in sequence
  //   }
  // }

  // for(int k = 0; k < xsec->GetNumParams(); k++){ // For every param in the xsec group
  //   if((xsec->GetParProp(k) == 0) && xsec->IsParFromGroup(k, "EParam")){ // If it has a value of 0 and is an energy normalisation parameter
  //     xsec->ToggleFixParameter(k); // Fix these params at 0
  //   }
  // }

  // #########################################################################
  // Perform reweight, print total integral and set the data

  std::vector<std::string> sample_names; // Create vector to store names of each sample
  std::vector<TH1*> DUNEHists;
  std::vector<std::string> MySamples{"FHC_numu", "FHC_nue"}; 
  
  TFile* PMNSData = TFile::Open("EventRates/OscPMNSNoNC.root"); // Load in the data we want to fit our Posterior Predictive to
  for (auto handler : DUNEPdfs) {
    for (unsigned iSample = 0; iSample < handler->GetNsamples(); ++iSample) {
    
    std::string name = handler->GetSampleTitle(iSample);
    sample_names.push_back(name);
    TString NameTString = TString(name.c_str());
    TString MyNameTString = TString(MySamples[iSample].c_str());

    handler->Reweight();

    TString HistName = "hRecoNeutrinoEnergy" + MyNameTString; // Name the histograms as they appear in PMNSData
    //TString HistName = "hFD_" + MyNameTString;
    TH1D* blarbHist = PMNSData->Get<TH1D>(HistName); // Get the histogram from our data
    TH1D* CloneHist = (TH1D*) blarbHist->Clone(); // Create clone of data
    CloneHist->SetDirectory(nullptr);

    DUNEHists.push_back(static_cast<TH1*>(CloneHist->Clone(NameTString+"_unosc"))); // Add our data histograms to the DUNE histograms
    //DUNEHists.push_back(Sample->GetMCHist(Sample->GetNDim()));
    //MACH3LOG_INFO("Event rate for {} : {:<5.2f}", Sample->GetTitle(), Sample->GetMCHist(Sample->GetNDim())->Integral());

      if (handler->GetNDim(iSample) == 1){
        handler->AddData(iSample, static_cast<TH1D*>(DUNEHists.back()));
      } else if (handler->GetNDim(iSample) == 2){
        handler->AddData(iSample, static_cast<TH2D*>(DUNEHists.back()));
      }
      
      else {
        MACH3LOG_ERROR("Unsupported number of dimensions > 2 - Quitting"); 
        throw MaCh3Exception(__FILE__ , __LINE__ );
      }

      MACH3LOG_INFO("Integrals of nominal hists: ");
      MACH3LOG_INFO("{} : {}",name.c_str(),DUNEHists.back()->Integral());
      MACH3LOG_INFO("--------------");
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
