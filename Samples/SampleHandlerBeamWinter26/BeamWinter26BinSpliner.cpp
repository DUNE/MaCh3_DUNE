#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/SampleHandlerBeamWinter26/SampleHandlerBeamWinter26.h"

#include "Fitters/Algorithms/MaCh3Factory.h"

#include "TH3.h"
#include "TSpline.h"

static const std::string SplineName = "Winter26Xsec";

void EmitXSecConfig(
    std::vector<std::pair<std::string, std::vector<float>>> const
        &EventSplineMetadata) {
  std::ofstream fout("cfg.yml");

  fout << "---\nSystematics:";

  std::string group = "Xsec";
  for (int pi = 0; pi < EventSplineMetadata.size(); ++pi) {
    auto [pname, pknots] = EventSplineMetadata[pi];

    for (auto &c : pname) {
      if (c == '_') {
        c = '-';
      }
    }

    if (pname.substr(0, 6) == "Target") {
      group = "Flux";
    }

    fout << fmt::format(
        R"""(
- Systematic:
    SampleNames: ["Beam*", ]
    Error: 1.0
    FlatPrior: false
    Names:
      FancyName: {0}
      ParameterName: {0}
    ParameterBounds:
    - -9999.0
    - 9999.0
    ParameterValues:
      Generated: 0.0
      PreFitValue: 0.0
    SplineInformation:
      Mode: [ 0 ]
      SplineName: {1}
    StepScale: 0.4
    Type: Spline
    ParameterGroup: {2}
)""",
        pname, SplineName, group);
    break;
  }
}

int main(int argc, char *argv[]) {

  auto sh = std::make_unique<dune::beam::SampleHandlerBeamWinter26<true>>(
      argv[1], nullptr, nullptr);

  auto &modes = *sh->GetMaCh3Modes();

  auto &Binning = sh->GetBinning();
  if (!Binning->IsUniform(0)) {
    throw std::runtime_error("Cannot work with nonuniform binning.");
  }
  int ndims = Binning->GetNDim(0);

  EmitXSecConfig(sh->EventSplineMetadata);

  // would be nice to check that these are the right projections, but
  // BinningHandler doesn't provide that info in the public interface
  auto mode_bins = Binning->GetBinEdges(0, 0);
  auto enu_bins = Binning->GetBinEdges(0, 1);
  auto xbins = Binning->GetBinEdges(0, 2);
  std::vector<double> ybins = {0, 1};

  TFile fout("fout.root", "RECREATE");
  switch (ndims - 2) {
  case 1: {
    TH2F *BinningHist =
        new TH2F("dev_tmp_0_0", "dev_tmp_0_0", enu_bins.size() - 1,
                 enu_bins.data(), xbins.size() - 1, xbins.data());
    fout.WriteTObject(BinningHist, "dev_tmp_0_0");
    break;
  }
  case 2: {
    ybins = Binning->GetBinEdges(0, 3);
    TH3F *BinningHist = new TH3F(
        "dev_tmp_0_0", "dev_tmp_0_0", enu_bins.size() - 1, enu_bins.data(),
        xbins.size() - 1, xbins.data(), ybins.size() - 1, ybins.data());
    fout.WriteTObject(BinningHist, "dev_tmp_0_0");
    break;
  }
  default:
    throw std::runtime_error(
        fmt::format("Cannot work with {} dim binning.", ndims));
  }

  bool havey = ybins.size() > 2;

  sh->Reweight();
  std::vector<double> nominal = sh->GetMCArray();

  for (int pi = 0; pi < sh->EventSplineMetadata.size(); ++pi) {
    auto const &[pname, pknots] = sh->EventSplineMetadata[pi];

    std::vector<TGraph> splines(Binning->GetNBins());
    std::vector<std::string> splnm(Binning->GetNBins());

    for (int ki = 0; ki < pknots.size(); ++ki) {
      sh->SetEventResponseWeights(pi, ki);
      sh->Reweight();
      auto const &response = sh->GetMCArray();

      std::vector<int> binidxs(ndims);
      for (int yi = 0; yi < ybins.size() - 1; ++yi) {

        if (havey) {
          binidxs[3] = yi;
        }

        for (int mi = 0; mi < mode_bins.size() - 1; ++mi) {
          binidxs[0] = mi;
          for (int ei = 0; ei < enu_bins.size() - 1; ++ei) {
            binidxs[1] = ei;
            for (int xi = 0; xi < xbins.size() - 1; ++xi) {
              binidxs[2] = xi;

              auto gbin = Binning->GetBinSafe(0, binidxs);
              double wght = response[gbin] / nominal[gbin];

              // if (response[gbin] <= 0 ||
              //     !std::isnormal(response[gbin])) { // skip empty bins
              //   continue;
              // }

              splines[gbin].SetPoint(ki, pknots[ki],
                                     std::isnormal(wght) ? wght : 1);

              if (!ki) { // name the spline
                std::string y = "";
                if (havey) {
                  y = fmt::format("_{}", yi);
                }
                splnm[gbin] =
                    fmt::format("dev_{}_{}_sp_{}_{}{}", SplineName,
                                modes.GetSplineSuffixFromMaCh3Mode(mi), ei, xi,
                                y);
              }
            }
          }
        }
      }
    }

    size_t written_sp = 0;
    for (int si = 0; si < splines.size(); ++si) {
      // if (!splines[si].GetN()) {
      //   continue;
      // }
      TSpline3 spline(splnm[si].c_str(), splines[si].GetX(),
                      splines[si].GetY(), splines[si].GetN());
      fout.WriteTObject(&spline, splnm[si].c_str());
      written_sp++;
    }

    fmt::print(
        "Done {}, wrote {}/{} splines which covered bins containing events.\n",
        pname, written_sp, splines.size());
    break;
  }
  fout.Close();
}
