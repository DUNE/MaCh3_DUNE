#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include <TCanvas.h>
#include <TColor.h>
#include <TH1D.h>
#include <TH3F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRint.h>
#include <TStyle.h>

#include "Fitters/MaCh3Factory.h"
#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/SampleHandlerBeamOffAxis.h"
#include "Samples/BeamOffAxis/Utility.h"
#include "Samples/StructsDUNE.h"

void ParseXsecYaml(std::unique_ptr<Manager> &FitManager,
                   std::vector<std::string> &ParamNames,
                   std::vector<std::string> &SplineNames,
                   std::vector<std::vector<int>> &SplineModes) {
  std::string SplineYAML =
      FitManager->raw()["General"]["Systematics"]["SplineMakingCov"]
          .as<std::string>();
  YAML::Node XsecYaml = M3OpenConfig(SplineYAML);
  for (auto const &param : XsecYaml["Systematics"]) {
    ParamNames.push_back(
        param["Systematic"]["Names"]["FancyName"].as<std::string>());
    SplineNames.push_back(param["Systematic"]["SplineInformation"]["SplineName"]
                              .as<std::string>());
    std::vector<int> Modes =
        param["Systematic"]["SplineInformation"]["Mode"].as<std::vector<int>>();
    SplineModes.push_back(Modes);
  }
}

int main(int argc, char *argv[]) {
  MaCh3Utils::MaCh3Usage(argc, argv);
  auto fitMan = MaCh3ManagerFactory(argc, argv);

  // ###############################################################################################################################
  // Create SampleHandlerFD objects

  ParameterHandlerGeneric *xsec = nullptr;

  std::vector<SampleHandlerFD *> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

  // ###############################################################################################################################
  // Perform reweight and print total integral

  int nshifts = 7;
  std::vector<double> SigmaKnots{-3., -2., -1., 0., 1., 2., 3.};
  std::vector<double> TrueEBinning{
      0.0, 1.0,  1.25, 1.5,  1.75, 2.0, 2.25, 2.5,  2.75,
      3.0, 3.25, 3.5,  3.75, 4.0,  5.0, 6.0,  10.0, 1200}; //{0.,
                                                           // 0.5, 1.0, 1.5, 2.,
                                                           // 2.5,
                                                           // 3., 3.5, 4., 4.5,
                                                           // 5., 6.,
                                                           // 1000.};

  size_t NTrueEbins = TrueEBinning.size() - 1;

  for (auto Sample : DUNEPdfs) {
    auto SampleOffAxis =
        dynamic_cast<dune::beamoffaxis::SampleHandlerBeamOffAxis *>(Sample);
    if (!SampleOffAxis) {
      throw "Casting Failed!";
    }
    MaCh3Modes *Modes = Sample->GetMaCh3Modes();

    for (unsigned iSample = 0; iSample < Sample->GetNsamples(); ++iSample) {

      std::vector<std::vector<std::vector<std::vector<std::unique_ptr<TH2D>>>>> BinnedWeights;
      std::string SubSampleName = Sample->GetSampleTitle(iSample);

      // Sample Binning
      std::vector<double> BinEdgesX =
          SampleOffAxis->ReturnKinematicParameterBinning(
              iSample, Sample->GetXBinVarName(iSample));
      int NBinsX = static_cast<int>(BinEdgesX.size()) - 1;

      std::vector<double> BinEdgesY =
          SampleOffAxis->ReturnKinematicParameterBinning(
              iSample, Sample->GetYBinVarName(iSample));
      int NBinsY = static_cast<int>(BinEdgesY.size()) - 1;

      // Make Template Binning Histogram
      TH3F *BinningTemplate = new TH3F(
          "dev_tmp_0_0", "dev_tmp_0_0", NTrueEbins, TrueEBinning.data(), NBinsX,
          BinEdgesX.data(), NBinsY, BinEdgesY.data());

      std::vector<std::string> ParamNames;
      std::vector<std::string> SplineNames;
      std::vector<std::vector<int>> SplineModes;
      ParseXsecYaml(fitMan, ParamNames, SplineNames, SplineModes);

      size_t NSplineParams = ParamNames.size();

      BinnedWeights = dune::beamoffaxis::GetBinnedWeights(
          *SampleOffAxis, iSample, ParamNames, SplineModes, TrueEBinning);

      // Make Splines

      std::string SplineFileName = SubSampleName + "_splines.root";
      auto SplineFile = std::unique_ptr<TFile>(
          TFile::Open(SplineFileName.c_str(), "RECREATE"));
      SplineFile->cd();

      for (size_t iParam = 0; iParam < NSplineParams; iParam++) {
        for (size_t mode = 0; mode < SplineModes[iParam].size(); mode++) {
          for (size_t TrueEbin = 0; TrueEbin < NTrueEbins; TrueEbin++) {
            for (int xbin = 0; xbin < NBinsX; xbin++) {
              for (int ybin = 0; ybin < NBinsY; ybin++) {
                auto WeightGraph = std::make_unique<TGraph>();
                for (int knot = 0; knot < nshifts; knot++) {
                  auto &hist = BinnedWeights[iParam][knot][mode][TrueEbin];

                  if (!hist) {
                    std::cout << "NULL histogram at "
                              << "iParam=" << iParam << ", knot=" << knot
                              << ", mode=" << mode << ", TrueEbin=" << TrueEbin
                              << std::endl;
                  }
                  std::cout << "Integral = " << hist->Integral() << std::endl;

                  WeightGraph->SetPoint(
                      knot, SigmaKnots[knot],
                      BinnedWeights[iParam][knot][mode][TrueEbin]
                                  ->GetBinContent(xbin + 1, ybin + 1) > 0.
                          ? BinnedWeights[iParam][knot][mode][TrueEbin]
                                ->GetBinContent(xbin + 1, ybin + 1)
                          : 1.);
                  double weight = BinnedWeights[iParam][knot][mode][TrueEbin]
                                      ->GetBinContent(xbin + 1, ybin + 1);

                  std::cout << "[DEBUG] "
                            << "iParam=" << iParam << ", mode=" << mode
                            << ", TrueEbin=" << TrueEbin << ", xbin=" << xbin
                            << ", ybin=" << ybin << ", knot=" << knot
                            << ", sigma=" << SigmaKnots[knot]
                            << ", weight=" << weight << std::endl;
                }

                auto SplineName =
                    fmt::format("dev_{}_{}_sp_{}_{}_{}", SplineNames[iParam],
                                Modes->GetSplineSuffixFromMaCh3Mode(
                                    SplineModes[iParam][mode]),
                                TrueEbin, xbin, ybin);

                auto Spline = std::make_unique<TSpline3>(SplineName.c_str(),
                                                         WeightGraph.get());
                SplineFile->WriteTObject(Spline.get(), SplineName.c_str());
              }
            }
          }
        }
      }

      BinningTemplate->Write();
      SplineFile->Close();
    }
  }
}
