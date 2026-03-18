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
  std::vector<double> TrueEBinning{0.0, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 5.0, 6.0, 10.0, 1200};   //{0., 0.5, 1.0, 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 6., 1000.};
  
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
    std::vector<double> BinEdgesX = SampleOffAxis->ReturnKinematicParameterBinning(Sample->GetXBinVarName());
    int NBinsX =  static_cast<int>(BinEdgesX.size()) - 1;


    std::vector<double> BinEdgesY = SampleOffAxis->ReturnKinematicParameterBinning(Sample->GetYBinVarName());
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

    std::string SplineFileName = SampleName + "_splines_offaxis.root";
    auto SplineFile = std::unique_ptr<TFile>(TFile::Open(SplineFileName.c_str(), "RECREATE"));
    SplineFile->cd();

    for (size_t iParam = 0; iParam < NSplineParams; iParam++) {
      for (size_t mode = 0; mode < SplineModes[iParam].size(); mode++) {
        for (size_t TrueEbin = 0; TrueEbin < NTrueEbins; TrueEbin++) {
          for (int xbin = 0; xbin < NBinsX; xbin++) {
            for (int ybin = 0; ybin < NBinsY; ybin++) {
              TGraph* WeightGraph = new TGraph();
	      for (int knot = 0; knot < nshifts; knot++) {
          auto* hist = BinnedWeights[iParam][knot][mode][TrueEbin];

          // if (!hist) {
          //     std::cout << "NULL histogram at "
          //               << "iParam=" << iParam
          //               << ", knot=" << knot
          //               << ", mode=" << mode
          //               << ", TrueEbin=" << TrueEbin
          //               << std::endl;
          // }
          // std::cout << "Integral = " << hist->Integral() << std::endl;

	        WeightGraph->SetPoint(knot, SigmaKnots[knot], BinnedWeights[iParam][knot][mode][TrueEbin]->GetBinContent(xbin+1, ybin+1) > 0. ? BinnedWeights[iParam][knot][mode][TrueEbin]->GetBinContent(xbin+1, ybin+1) : 1.);
          double weight = BinnedWeights[iParam][knot][mode][TrueEbin]->GetBinContent(xbin+1, ybin+1);

          // std::cout << "[DEBUG] "
          //           << "iParam=" << iParam
          //           << ", mode=" << mode
          //           << ", TrueEbin=" << TrueEbin
          //           << ", xbin=" << xbin
          //           << ", ybin=" << ybin
          //           << ", knot=" << knot
          //           << ", sigma=" << SigmaKnots[knot]
          //           << ", weight=" << weight
          //           << std::endl;

	      }
        // Check if spline is flat (all weights are 1.0)
        bool isFlat = true;
        for (int knot = 0; knot < nshifts; knot++) {
          double w = BinnedWeights[iParam][knot][mode][TrueEbin]->GetBinContent(xbin+1, ybin+1);
          if (std::abs(w - 1.0) > 1e-6 && w > 0.) { isFlat = false; break; }
        }
        // if (isFlat) {
        //   std::cout << "[FLAT SPLINE] "
        //             << "Param=" << ParamNames[iParam]
        //             << ", Mode=" << Modes->GetSplineSuffixFromMaCh3Mode(SplineModes[iParam][mode])
        //             << ", TrueEbin=" << TrueEbin
        //             << ", erec=" << xbin
        //             << ", leprec=" << ybin
        //             << std::endl;

        // }
        // std::cout << "Enurec bin edges: ";
        // for (auto e : BinEdgesX) std::cout << e << " ";
        // std::cout << std::endl;
        // std::cout << "RecoLep bin edges: ";
        // for (auto e : BinEdgesY) std::cout << e << " ";
        // std::cout << std::endl;
        
	      char SplineName[1000];
	      sprintf(SplineName, "dev_%s_%s_sp_%d_%d_%d", SplineNames[iParam].c_str(), Modes->GetSplineSuffixFromMaCh3Mode(SplineModes[iParam][mode]).c_str(), TrueEbin, xbin, ybin);
	      TSpline3* Spline = new TSpline3(SplineName, WeightGraph);
	      delete WeightGraph;
	      Spline->Write(SplineName);
	    }
	  }
	}
      }
    }

    // After all splines are written, print summary of non-flat splines
std::vector<std::string> nonFlatSummary;

for (size_t iParam = 0; iParam < NSplineParams; iParam++) {
  for (size_t mode = 0; mode < SplineModes[iParam].size(); mode++) {
    for (size_t TrueEbin = 0; TrueEbin < NTrueEbins; TrueEbin++) {
      for (int xbin = 0; xbin < NBinsX; xbin++) {
        for (int ybin = 0; ybin < NBinsY; ybin++) {

          bool isFlat = true;
          for (int knot = 0; knot < nshifts; knot++) {
            double w = BinnedWeights[iParam][knot][mode][TrueEbin]->GetBinContent(xbin+1, ybin+1);
            if (w > 0. && std::abs(w - 1.0) > 1e-6) { isFlat = false; break; }
          }

          if (!isFlat) {
            double xLow = BinEdgesX[xbin];
            double xHigh = BinEdgesX[xbin+1];
            double yLow = BinEdgesY[ybin];
            double yHigh = BinEdgesY[ybin+1];
            double eLow = TrueEBinning[TrueEbin];
            double eHigh = TrueEBinning[TrueEbin+1];

            std::ostringstream oss;
            oss << "Param=" << ParamNames[iParam]
                << " Mode=" << Modes->GetSplineSuffixFromMaCh3Mode(SplineModes[iParam][mode])
                << " TrueE=[" << eLow << "," << eHigh << "]"
                << " Erec=[" << xLow << "," << xHigh << "]"
                << " LepRec=[" << yLow << "," << yHigh << "]"
                << " weights: [";
            for (int k = 0; k < nshifts; k++) {
              oss << BinnedWeights[iParam][k][mode][TrueEbin]->GetBinContent(xbin+1, ybin+1);
              if (k < nshifts-1) oss << ", ";
            }
            oss << "]";
            nonFlatSummary.push_back(oss.str());
          }
        }
      }
    }
  }
}

std::cout << "\n========== NON-FLAT SPLINE SUMMARY (" << nonFlatSummary.size() << " total) ==========" << std::endl;
for (auto const& s : nonFlatSummary) std::cout << s << std::endl;
std::cout << "================================================================\n" << std::endl;

    BinningTemplate->Write();
    SplineFile->Close();
  }

}
