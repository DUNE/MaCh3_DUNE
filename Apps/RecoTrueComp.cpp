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
#include <filesystem>
#include <TFile.h>
#include <TTree.h>

#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Fitters/MaCh3Factory.h"

int main(int argc, char * argv[]) {
  MaCh3Utils::MaCh3Usage(argc, argv);
  auto fitMan = MaCh3ManagerFactory(argc, argv);

  //###############################################################################################################################
  //Create SampleHandlerFD objects
  
  ParameterHandlerGeneric* xsec = nullptr;
  
  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

  //###############################################################################################################################
  //Perform reweight and print total integral

  namespace fs = std::filesystem;
  std::string folder = "Inputs/DUNE_CAF_files/";
  std::string prefix = "FD_RHC";
  std::string suffix = "numuselec.root";
  std::string suffix1 = "numuselec.root";
  std::string sample = "RHC_numu";
  std::vector<std::string> files;
  double RecoNeutrinoEnergy;
  double TrueNeutrinoEnergy;
  TH2D* hRecoTrue = new TH2D("hRecoTrue", ("Reco vs True Neutrino Energy for " + sample).c_str(), 120, 0, 30, 120, 0, 30);

  for (const auto& entry : fs::directory_iterator(folder)) {
    std::string name = entry.path().filename().string();

    if (name.rfind(prefix, 0) != 0) continue;
    bool match_suffix = (name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0) || (name.compare(name.size() - suffix1.size(), suffix1.size(), suffix1) == 0);
    if (!match_suffix) continue;
    files.push_back(entry.path().string());
  }

  for (const auto& filename : files) {
    std::cout << "Processing file: " << filename << std::endl;

    auto ends_with = [](const std::string& s, const std::string& suffix) {
      return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
    };

    bool is_nue = ends_with(filename, suffix);
    bool is_numu = ends_with(filename, suffix1);
    
    TFile* CAFFile = TFile::Open(filename.c_str());
    TTree* CAFTree = CAFFile->Get<TTree>("caf");

    CAFTree->SetBranchAddress("Ev", &TrueNeutrinoEnergy);

    if (is_nue) {
      CAFTree->SetBranchAddress("Ev_reco_nue", &RecoNeutrinoEnergy);
    } else if (is_numu) {
      CAFTree->SetBranchAddress("Ev_reco_numu", &RecoNeutrinoEnergy);
    }

    Long64_t nEntries = CAFTree->GetEntries();

    for (Long64_t i = 0; i < nEntries; i++) {
      CAFTree->GetEntry(i);
      hRecoTrue->Fill(TrueNeutrinoEnergy, RecoNeutrinoEnergy);
    }
  }

  auto Canvas = new TCanvas();
  Canvas->Draw();
  Canvas->cd();
  Canvas->SetLogz();
  hRecoTrue->Draw("COLZ");
  hRecoTrue->SetStats(0);
  hRecoTrue->GetXaxis()->SetTitle("True Neutrino Energy (GeV)"); 
  hRecoTrue->GetYaxis()->SetTitle("Reco Neutrino Energy (GeV)");
  hRecoTrue->GetZaxis()->SetTitle("Number of Events");
  hRecoTrue->GetZaxis()->SetTitleOffset(1.5);
  Canvas->SetRightMargin(0.15);
  Canvas->Update();
  Canvas->Print((sample + "RecoTrueComp.pdf").c_str());

  TFile outFile((sample + "RecoTrueComp.root").c_str(), "RECREATE");
  hRecoTrue->Write();
  outFile.Close();

}