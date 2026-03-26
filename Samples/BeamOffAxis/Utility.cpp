#include "Samples/BeamOffAxis/Utility.h"

#include "Samples/SampleHandlerBeamOffAxis.h"

#include <cmath>

double ERecQE(int nupdg, bool isCC, double el_MeV, double theta_lep_rad) {
  constexpr double V = 0;        // 0 binding energy for now
  constexpr double mn = 939.565; // neutron mass
  constexpr double mp = 938.272; // proton mass
  double mN_eff = mn - V;
  double mN_oth = mp;

  if (nupdg < 0) { // if anti-neutrino, swap target/out masses
    mN_eff = mp - V;
    mN_oth = mn;
  }

  double ml = 0;

  switch (std::abs(nupdg)) {
  case 12:
    ml = isCC ? 0.511 : 0;
    break;
  case 14:
    ml = isCC ? 105.66 : 0;
    break;
  case 16:
    ml = isCC ? 1777.0 : 0;
    break;
  default:
    std::cerr << "Warning: Unexpected PDG code " << nupdg
              << " passed to ml lambda.\n";
    assert(false && "Unexpected neutrino PDG code in ml lambda");
  }

  double pl = std::sqrt(el_MeV * el_MeV - ml * ml); // momentum of lepton

  return (2 * mN_eff * el_MeV - ml * ml + mN_oth * mN_oth - mN_eff * mN_eff) /
         (2 * (mN_eff - el_MeV + pl * std::cos(theta_lep_rad)));
}

namespace dune::beamoffaxis {

std::vector<std::vector<std::vector<std::vector<std::unique_ptr<TH2D>>>>>
GetBinnedWeights(SampleHandlerBeamOffAxis &sample, int iSubSample,
                 std::vector<std::string> ParamNames,
                 std::vector<std::vector<int>> ParamModes,
                 std::vector<double> TrueEBins) {

  std::cout << fmt::format("Making Binned Weights for Sample: {}/{}",
                           sample.GetName(), sample.GetSampleTitle(iSubSample))
            << std::endl;

  // use sample.GetNDim(iSubSample) to work outif 1D or 2D

  // Vector Structure:
  // Parameter<Knot<Mode<TrueE<TH2D>>>
  std::vector<std::vector<std::vector<std::vector<std::unique_ptr<TH2D>>>>>
      histVec, NomVec;
  int nshifts = 7;

  // True Energy Binning
  int NTrueEBins = static_cast<int>(TrueEBins.size()) - 1;
  TH1D *TrueEbinning =
      new TH1D("Template True E binning", "", NTrueEBins, TrueEBins.data());

  // Sample Binning
  std::vector<double> BinEdgesX = sample.ReturnKinematicParameterBinning(
      iSubSample, sample.GetXBinVarName(iSubSample));
  int NBinsX = static_cast<int>(BinEdgesX.size()) - 1;

  std::vector<double> BinEdgesY = sample.ReturnKinematicParameterBinning(
      iSubSample, sample.GetYBinVarName(iSubSample));
  int NBinsY = static_cast<int>(BinEdgesY.size()) - 1;

  // Setup Histograms
  histVec.resize(ParamNames.size());
  NomVec.resize(ParamNames.size());

  for (size_t iParam = 0; iParam < ParamNames.size(); iParam++) {
    histVec[iParam].resize(nshifts);
    NomVec[iParam].resize(nshifts);
    for (int shift = 0; shift < nshifts; shift++) {
      histVec[iParam][shift].resize(ParamModes[iParam].size());
      NomVec[iParam][shift].resize(ParamModes[iParam].size());
      for (size_t mode = 0; mode < ParamModes[iParam].size(); mode++) {
        histVec[iParam][shift][mode].resize(NTrueEBins);
        NomVec[iParam][shift][mode].resize(NTrueEBins);
        for (int b_etrue = 0; b_etrue < NTrueEBins; b_etrue++) {

          histVec[iParam][shift][mode][b_etrue] = std::make_unique<TH2D>(
              Form("bn_p:%zu_s:%d_m:%zu_e:%d", iParam, shift, mode, b_etrue),
              "", NBinsX, BinEdgesX.data(), NBinsY, BinEdgesY.data());

          NomVec[iParam][shift][mode][b_etrue] = std::make_unique<TH2D>(
              Form("nom_p:%zu_s:%d_m:%zu_e:%d", iParam, shift, mode, b_etrue),
              "", NBinsX, BinEdgesX.data(), NBinsY, BinEdgesY.data());
        }
      }
    }
  }

  // TChain to store MC which contains weight information
  TChain CAFChain("cafTree");

  for (const std::string &filename :
       sample.SampleDetails[iSubSample].mc_files) {
    MACH3LOG_INFO("Adding file to TChains: {}", filename);

    // HH: Check whether the file exists, see
    // https://root.cern/doc/master/classTChain.html#a78a896924ac6c7d3691b7e013bcbfb1c
    if (!CAFChain.Add(filename.c_str(), -1)) {
      MACH3LOG_ERROR("Could not add file {} to TChain, please check the file "
                     "exists and is readable",
                     filename);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  MACH3LOG_INFO("Number of entries in CAF TChain: {}", CAFChain.GetEntries());

  CAFChain.SetBranchStatus("*", 0);

  std::vector<std::vector<double>> weightArr(ParamNames.size(),
                                             std::vector<double>(nshifts));
  std::vector<double> cvweight(ParamNames.size());

  for (size_t iParam = 0; iParam < ParamNames.size(); iParam++) {
    std::string WeightBranchName = "wgt_" + ParamNames[iParam];
    std::string CVBranchName = ParamNames[iParam] + "_cvwgt";
    CAFChain.SetBranchStatus(WeightBranchName.c_str(), 1);
    CAFChain.SetBranchStatus(CVBranchName.c_str(), 1);

    CAFChain.SetBranchAddress(WeightBranchName.c_str(),
                              weightArr[iParam].data());
    CAFChain.SetBranchAddress(CVBranchName.c_str(), &cvweight[iParam]);
  }

  int NEvents = static_cast<int>(CAFChain.GetEntries());

  for (int iEvent = 0; iEvent < NEvents; iEvent++) {
    CAFChain.GetEntry(iEvent);
    double x_var = sample.ReturnKinematicParameter(
        sample.GetXBinVarName(iSubSample), iEvent);
    double y_var = sample.ReturnKinematicParameter(
        sample.GetYBinVarName(iSubSample), iEvent);

    // skip if event does not pass selections
    if (!sample.IsEventSelected(iSubSample, iEvent)) {
      continue;
    }

    int TrueEbin =
        TrueEbinning->FindBin(sample.DUNEMCEvents[iEvent].truth.nu.e) - 1;
    if (iEvent < 10) {
      std::cout << "Event " << iEvent
                << " truth.nu.e=" << sample.DUNEMCEvents[iEvent].truth.nu.e
                << " truth.lep.e=" << sample.DUNEMCEvents[iEvent].truth.lep.e
                << " TrueEbin=" << TrueEbin << std::endl;
    }
    std::cout << "TrueE bin edges: ";
    for (auto e : TrueEBins)
      std::cout << e << " ";
    std::cout << std::endl;

    for (size_t iParam = 0; iParam < ParamNames.size(); iParam++) {

      std::string FancyName = ParamNames[iParam];
      std::string WeightBranchName = "wgt_" + FancyName;
      std::string CVBranchName = FancyName + "_cvwgt";

      for (size_t mode = 0; mode < ParamModes[iParam].size(); mode++) {

        if (ParamModes[iParam][mode] ==
            sample.DUNEMCEvents[iEvent].truth.mach3_mode) {

          // double weightArr[1000];
          // double cvweight;
          //  Vector Structure:
          //  Parameter<Knot<ETrue<InteractionMode<TH2D>>>

          // CAFChain.SetBranchAddress(WeightBranchName.c_str(), &weightArr);
          // CAFChain.SetBranchAddress(CVBranchName.c_str(), &cvweight);
          // CAFChain.GetEntry(iEvent);

          // Debug print for weights being read
          std::cout << "Event " << iEvent << " | Param: " << ParamNames[iParam]
                    << " | Mode: " << ParamModes[iParam][mode]
                    << " | CVWeight: " << cvweight[iParam]
                    << " | TrueBin: " << TrueEbin << " | Weights: [";
          for (int s = 0; s < nshifts; s++) {
            std::cout << weightArr[iParam][s];
            if (s < nshifts - 1)
              std::cout << ", ";
          }
          std::cout << "]" << std::endl;

          for (int shift = 0; shift < nshifts; shift++) {
            histVec[iParam][shift][mode][TrueEbin]->Fill(
                x_var, y_var, weightArr[iParam][shift]);
            NomVec[iParam][shift][mode][TrueEbin]->Fill(x_var, y_var, 1.0);
          }
        }
      }
    }
  }

  // Summary of filled histograms before taking ratio
  for (size_t iParam = 0; iParam < ParamNames.size(); iParam++) {

    for (size_t mode = 0; mode < ParamModes[iParam].size(); mode++) {
      for (int b_etrue = 0; b_etrue < NTrueEBins; b_etrue++) {
        for (int shift = 0; shift < nshifts; shift++) {
          int nBinsX = histVec[iParam][shift][mode][b_etrue]->GetNbinsX();
          int nBinsY = histVec[iParam][shift][mode][b_etrue]->GetNbinsY();
          for (int ix = 1; ix <= nBinsX; ix++) {
            for (int iy = 1; iy <= nBinsY; iy++) {
              std::cout << "Param: " << ParamNames[iParam]
                        << " | Shift: " << shift
                        << " | Mode: " << ParamModes[iParam][mode]
                        << " | TrueEBin: " << b_etrue << " | xbin: " << ix
                        << " | ybin: " << iy << " | Content (weighted): "
                        << histVec[iParam][shift][mode][b_etrue]->GetBinContent(
                               ix, iy)
                        << " | Content (nominal): "
                        << NomVec[iParam][shift][mode][b_etrue]->GetBinContent(
                               ix, iy)
                        << std::endl;
            }
          }
        }
      }
    }
  }

  // Take the ratio
  for (size_t iParam = 0; iParam < ParamNames.size(); iParam++) {
    for (int shift = 0; shift < nshifts; shift++) {
      for (size_t mode = 0; mode < ParamModes[iParam].size(); mode++) {
        for (int b_etrue = 0; b_etrue < NTrueEBins; b_etrue++) {
          histVec[iParam][shift][mode][b_etrue]->Divide(
              NomVec[iParam][shift][mode][b_etrue].get());
        }
      }
    }
  }

  std::cout << "finished making binned weights" << std::endl;
  // return final vector of response histograms
  return histVec;
}

} // namespace dune::beamoffaxis
