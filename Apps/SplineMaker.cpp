#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>

#include <TH1D.h>
#include <TH3F.h>
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
#include "Samples/SampleHandlerBeamOffAxis.h"

void ParseXsecYaml(std::unique_ptr<manager>& FitManager, std::vector<std::string>& ParamNames, std::vector<std::string>& SplineNames, std::vector<std::vector<int>>& SplineModes) {
  std::string SplineYAML = FitManager->raw()["General"]["Systematics"]["SplineMakingCov"].as<std::string>();
  YAML::Node XsecYaml = M3OpenConfig(SplineYAML);
  for (auto const &param : XsecYaml["Systematics"]) {
    ParamNames.push_back(param["Systematic"]["Names"]["FancyName"].as<std::string>());
    SplineNames.push_back(param["Systematic"]["SplineInformation"]["SplineName"].as<std::string>());
    std::vector<int> Modes = param["Systematic"]["SplineInformation"]["Mode"].as<std::vector<int>>();
    SplineModes.push_back(Modes);
  }
}

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
  
  int nshifts = 7;
  std::vector<double> SigmaKnots{-3., -2., -1., 0., 1., 2., 3.};
  std::vector<double> TrueEBinning{0., 0.5, 1.0, 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 1000.};
  size_t NTrueEbins = TrueEBinning.size() - 1;

  
  for (auto Sample : DUNEPdfs){

    auto SampleOffAxis = dynamic_cast<SampleHandlerBeamOffAxis*>(Sample);
    if(!SampleOffAxis){
      std::cout << "Casting Failed!" << std::endl; 
    }
    std::vector<std::vector<std::vector<std::vector<TH2D*>>>> BinnedWeights;
    std::string SampleName = Sample->GetTitle();
    MaCh3Modes* Modes = Sample->GetMaCh3Modes();

    //Sample Binning
    std::vector<double> BinEdgesX = Sample->ReturnKinematicParameterBinning(Sample->GetXBinVarName());
    int NBinsX =  static_cast<int>(BinEdgesX.size()) - 1;


    std::vector<double> BinEdgesY = Sample->ReturnKinematicParameterBinning(Sample->GetYBinVarName()); 
    int NBinsY =  static_cast<int>(BinEdgesY.size()) - 1;

    //Make Template Binning Histogram
    TH3F* BinningTemplate = new TH3F("dev_tmp_0_0", "dev_tmp_0_0", NTrueEbins, TrueEBinning.data(), NBinsX, BinEdgesX.data(), NBinsY, BinEdgesY.data());

    std::vector<std::string> ParamNames;
    std::vector<std::string> SplineNames;
    std::vector<std::vector<int>> SplineModes;
    ParseXsecYaml(fitMan, ParamNames, SplineNames, SplineModes);

    size_t NSplineParams = ParamNames.size();

    BinnedWeights = SampleOffAxis->GetBinnedWeights(ParamNames, SplineModes, TrueEBinning);

    //Make Splines
    
    std::string SplineFileName = SampleName + "_splines.root";
    auto SplineFile = std::unique_ptr<TFile>(TFile::Open(SplineFileName.c_str(), "RECREATE"));
    SplineFile->cd();

    for (size_t iParam = 0; iParam < NSplineParams; iParam++) {
      for (size_t mode = 0; mode < SplineModes[iParam].size(); mode++) {
        for (size_t TrueEbin = 0; TrueEbin < NTrueEbins; TrueEbin++) {
          for (int xbin = 0; xbin < NBinsX; xbin++) {
            for (int ybin = 0; ybin < NBinsY; ybin++) {
              TGraph* WeightGraph = new TGraph();
	      for (int knot = 0; knot < nshifts; knot++) {
	        WeightGraph->SetPoint(knot, SigmaKnots[knot], BinnedWeights[iParam][knot][mode][TrueEbin]->GetBinContent(xbin+1, ybin+1));
	      }
	      char SplineName[1000];
	      sprintf(SplineName, "dev_%s_%s_sp_%d_%d_%d", SplineNames[iParam].c_str(), Modes->GetSplineSuffixFromMaCh3Mode(mode).c_str(), TrueEbin, xbin, ybin);
	      TSpline3* Spline = new TSpline3(SplineName, WeightGraph);
	      delete WeightGraph;
	      Spline->Write(SplineName);
	    }
	  }
	}
      }
    }

    BinningTemplate->Write();
    SplineFile->Close();
  }

}
