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
#include <yaml-cpp/yaml.h>
#include <fstream>
#include <string>
#include <map>

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

  std::vector<double> EParams;

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
        EParams.push_back(Param);
      }
      else { // If unoscillated data, calculate ratio between these as needed to induce oscillation
        Param = BinSizeOsc / BinSizeUnosc; 
        xsec->SetPar(ParIndex, Param);
        EParams.push_back(Param);
      } 
      ParIndex++; // Increment parameter index to keep amending in sequence
    }
  }

  const int NumBins = 40;
  const double EMin = 0.0;
  const double EMax = 10.0;
  const double BinWidth = (EMax - EMin) / NumBins;

  std::map<int, std::string> NuType = {{12, "nue"}, {14, "numu"}, {16, "nutau"}, {-12, "nuebar"}, {-14, "numubar"}, {-16, "nutaubar"}}; // Creating map between neutrino "value" and "name"
  std::vector<std::pair<int,int>> OscChan = {{14, 14}, {12, 12}, {-14, -14}, {-12, -12}, {14, 12}, {-14, -12}, {12, 14}, {-12, -14}, {14, 16}, {12, 16}, {-14, -16}, {-12, -16}}; // Creating all oscillation channels as pairs of neutrino "value"

  YAML::Node systematics(YAML::NodeType::Sequence);

  int NewIndex = 0;

  for(const auto& chan : OscChan){
    int nu_unosc = chan.first;
    int nu_osc   = chan.second;
    std::string unosc_name = NuType[nu_unosc]; 
    std::string osc_name   = NuType[nu_osc];
    for(int b = 0; b < NumBins; b++){
        double ELow = EMin + b * BinWidth;
        double EHigh = ELow + BinWidth;
        char lowbuf[16], highbuf[16];
        snprintf(lowbuf,  sizeof(lowbuf),  "%.2f", ELow);
        snprintf(highbuf, sizeof(highbuf), "%.2f", EHigh);
        std::string pname = unosc_name + "_" + osc_name + "_" + lowbuf + "_" + highbuf;
        double Value = EParams[NewIndex];        
        YAML::Node block;
        block["Systematic"]["SampleNames"].push_back("BeamFD");
        block["Systematic"]["Error"] = 1.0;
        block["Systematic"]["FlatPrior"] = true;
        if(Value == 0){
            block["Systematic"]["FixParam"] = true;
        }
        else{
            block["Systematic"]["FixParam"] = false;
        } 

        block["Systematic"]["KinematicCuts"].push_back(YAML::Node());
        block["Systematic"]["KinematicCuts"][0]["TrueNeutrinoEnergy"].push_back(ELow);
        block["Systematic"]["KinematicCuts"][0]["TrueNeutrinoEnergy"].push_back(EHigh);

        block["Systematic"]["Names"]["FancyName"] = pname;
        block["Systematic"]["Names"]["ParameterName"] = pname;

        block["Systematic"]["NeutrinoFlavourUnosc"].push_back(nu_unosc);
        block["Systematic"]["NeutrinoFlavour"].push_back(nu_osc);

        block["Systematic"]["ParameterBounds"].push_back(0.0);
        block["Systematic"]["ParameterBounds"].push_back(1.0);

        block["Systematic"]["ParameterGroup"] = "EParam";

        block["Systematic"]["ParameterValues"]["Generated"] = Value;
        block["Systematic"]["ParameterValues"]["PreFitValue"] = Value;

        block["Systematic"]["StepScale"]["MCMC"] = 1.0;
        block["Systematic"]["Type"] = "Norm";

        systematics.push_back(block);

        NewIndex++;
    }
  }
  YAML::Node root;
  root["Systematics"] = systematics;

  std::ofstream fout("EParams.yaml");
  fout << root;
  fout.close();

  std::cout << "Wrote " << systematics.size() << " parameters to EParams.yaml" << std::endl;
}
  //###########################################################################################################

