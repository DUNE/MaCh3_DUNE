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

#include "Fitters/MaCh3Factory.h"
#include "Samples/MaCh3DUNEFactory.h"

int main(int argc, char * argv[]) {

  auto FitManager = MaCh3ManagerFactory(argc, argv);
  auto OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();

  ParameterHandlerGeneric* xsec = nullptr;

  // #########################################################################
  // Produce sample PDF objects

  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  // #########################################################################
  // Load in and create all relevant energy parameters

  // std::vector<double> Params;
  // std::vector<double> AvgParams(480, 0.0);
  // std::vector<int> Count(480, 0);

  TFile* Osc = TFile::Open("CCIndChanOsc.root"); // Loading in histograms
  TFile* Unosc = TFile::Open("CCIndChanUnosc.root");

  TIter next(Osc->GetListOfKeys()); // Get list of different items within oscillated data histograms (48 in total, 12 channels in 4 samples)
  TKey* key;  // Initialise
  int ParIndex = 0.0;
  int KeyIndex = 0.0;

  for(int i = 0; i < xsec->GetNumParams(); i++){ // For every param in the xsec group
    if(xsec->IsParFromGroup(i, "EParam")){ // If param is from our energy normalisation parameter group
      ParIndex = i; // Set the index of first parameter
      break;
    }
  }

  while ((key = (TKey*)next())) { // Go through all keys in sequence
    if(KeyIndex == 12.0) break;
    KeyIndex++;
    auto HistoOsc = Osc->Get<TH1D>(key->GetName()); // Getting names of histograms
    auto HistoUnosc = Unosc->Get<TH1D>(key->GetName()); 
    int NumBins = HistoOsc->GetNbinsX(); // Find number of bins (number of energy normalisation parameters for this histogram)
    for(int j = 1; j <= NumBins; j++) { // For each energy bin, starting from 1 to avoid the overflow bin
      double BinSizeOsc = HistoOsc->GetBinContent(j); // Get the number of events in specific energy bin
      double BinSizeUnosc = HistoUnosc->GetBinContent(j);
      double Param;
      if(BinSizeUnosc == 0) { // If no unoscillated data, set param to 0
        Param = 0;
        xsec->SetPar(ParIndex, Param);
        //Params.push_back(Param);
      }
      else { // If unoscillated data, calculate ratio between these as needed to induce oscillation
        Param = BinSizeOsc / BinSizeUnosc; 
        xsec->SetPar(ParIndex, Param);
        //Params.push_back(Param);
      } 
      ParIndex++; // Increment parameter index to keep amending in sequence
    }
  }
  // for(int p = 0; p < 1920; p++){
  //   int sample = p / 480;
  //   int channel = (p % 480) / 40;
  //   int bin = p % 40;
  //   int Index = channel * 40 + bin;
  //   AvgParams[Index] += Params[p];
  //   Count[Index] += 1;
  // }
  // for(int i = 0; i < 480; i++){
  //   AvgParams[i] /= Count[i];
  //   xsec->SetPar(ParIndex, AvgParams[i]);
  //   ParIndex++;
  // }

  for(int k = 0; k < xsec->GetNumParams(); k++){ // For every param in the xsec group
    if((xsec->GetParProp(k) == 0) && xsec->IsParFromGroup(k, "EParam")){ // If it has a value of 0 and is an energy normalisation parameter
      xsec->ToggleFixParameter(k); // Fix these params at 0
    }
  }
  
  // #########################################################################
  // Perform reweight, print total integral and set the data

  //Some place to store the histograms
  std::vector<TH1*> PredictionHistograms;
  std::vector<std::string> sample_names;

  auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();

  TFile* PMNSData = TFile::Open("OscPMNSNoNC.root"); // Loading in the oscillated data we want to fit to
  for (unsigned sample_i = 0 ; sample_i < DUNEPdfs.size() ; ++sample_i) {
    
    std::string name = DUNEPdfs[sample_i]->GetTitle();
    sample_names.push_back(name);
    TString NameTString = TString(name.c_str());
    
    DUNEPdfs[sample_i] -> Reweight();

    TString HistName = "hRecoNeutrinoEnergy" + NameTString; // Name the histograms as they appear in PMNSData
    TH1D* blarbHist = PMNSData->Get<TH1D>(HistName); // Get the histogram from our data
    TH1D* CloneHist = (TH1D*) blarbHist->Clone(); // Create clone of data
    CloneHist->SetDirectory(nullptr);

    PredictionHistograms.push_back(static_cast<TH1*>(CloneHist->Clone(NameTString+"_unosc"))); // Add our data histograms to the DUNE histograms
    //PredictionHistograms.push_back(static_cast<TH1*>(DUNEPdfs[sample_i]->GetMCHist(DUNEPdfs[sample_i]->GetNDim())->Clone(NameTString+"_unosc"))); // Use this for PMNS data?

    if (DUNEPdfs[sample_i]->GetNDim() == 1){
      DUNEPdfs[sample_i]->AddData(static_cast<TH1D*>(PredictionHistograms[sample_i]));
    } else if (DUNEPdfs[sample_i]->GetNDim() == 2){
      DUNEPdfs[sample_i]->AddData(static_cast<TH2D*>(PredictionHistograms[sample_i]));
    }
    
    else {
      MACH3LOG_ERROR("Unsupported number of dimensions > 2 - Quitting"); 
      throw MaCh3Exception(__FILE__ , __LINE__ );
    } 
  }
  
  //Now print out some event rates, we'll make a nice latex table at some point 
  for (unsigned iPDF = 0; iPDF < DUNEPdfs.size() ; ++iPDF) {
    MACH3LOG_INFO("Integrals of nominal hists: ");
    MACH3LOG_INFO("{} : {}",sample_names[iPDF].c_str(),PredictionHistograms[iPDF]->Integral());
    MACH3LOG_INFO("--------------");
  }
  
  //###########################################################################################################
  //MCMC

  auto MaCh3Fitter = MaCh3FitterFactory(FitManager.get());

  bool StartFromPreviousChain = GetFromManager(FitManager->raw()["General"]["StartFromPos"], false);
  //Start chain from random position unless continuing a chain
  if(!StartFromPreviousChain){
    if (!GetFromManager(FitManager->raw()["General"]["StatOnly"], false)) {
      xsec->ThrowParameters();
    }
  }

  //Add systematic objects
  MaCh3Fitter->AddSystObj(xsec);

  if (StartFromPreviousChain) {
    std::string PreviousChainPath = FitManager->raw()["General"]["PosFileName"].as<std::string>();
    MACH3LOG_INFO("MCMC getting starting position from: {}",PreviousChainPath);
    MaCh3Fitter->StartFromPreviousFit(PreviousChainPath);
  }
  
  //Add samples
  for(auto Sample : DUNEPdfs){
    MaCh3Fitter->AddSampleHandler(Sample);
  }

  
  //Run fit
  MaCh3Fitter->RunMCMC();

  //Writing the memory usage at the end to eventually spot some nasty leak
  MACH3LOG_WARN("\033[0;31mCurrent Total RAM usage is {:.2f} GB\033[0m", MaCh3Utils::getValue("VmRSS") / 1048576.0);
  MACH3LOG_WARN("\033[0;31mOut of Total available RAM {:.2f} GB\033[0m", MaCh3Utils::getValue("MemTotal") / 1048576.0);

  return 0;
}
