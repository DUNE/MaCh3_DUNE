#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/SampleHandlerBeamWinter26/SampleHandlerBeamWinter26.h"

#include "Fitters/Algorithms/MaCh3Factory.h"

#include "TH3.h"

int main(int argc, char *argv[]) {

  auto sh = std::make_unique<dune::beam::SampleHandlerBeamWinter26<true>>(
      argv[1], nullptr, nullptr);

  auto const &modes = *sh->GetMaCh3Modes();

  bool dummy;
  auto enu_bins = BuildBinEdgesFromNode(
      YAML::Load("{ linspace: { nb: 5, low: 0, up: 10 } }"), dummy);

  auto &Binning = sh->GetBinning();
  if (!Binning->IsUniform(0)) {
    throw std::runtime_error("Cannot work with nonuniform binning.");
  }
  int ndims = Binning->GetNDim(0);

  TFile fout("fout.root", "RECREATE");
  switch (ndims) {
  case 1: {
    auto xbins = Binning->GetBinEdges(0, 0);
    TH2F *BinningHist =
        new TH2F("dev_tmp.0.0", "dev_tmp.0.0", enu_bins.size() - 1,
                 enu_bins.data(), xbins.size() - 1, xbins.data());
    break;
  }
  case 2: {
    auto xbins = Binning->GetBinEdges(0, 0);
    auto ybins = Binning->GetBinEdges(0, 1);
    TH3F *BinningHist = new TH3F(
        "dev_tmp.0.0", "dev_tmp.0.0", enu_bins.size() - 1, enu_bins.data(),
        xbins.size() - 1, xbins.data(), ybins.size() - 1, ybins.data());
    break;
  }
  default:
    throw std::runtime_error(
        fmt::format("Cannot work with {} dim binning.", ndims));
  }

  sh->UnsetWeights();
  sh->Reweight();
  auto h = sh->GetMCHist(0);
  auto nomdir = fout.mkdir("nominal");
  auto spdir = fout.mkdir("splines");
  nomdir->WriteTObject(h, "hist");

  std::map<int, std::map<int, std::unique_ptr<TH1>>> mode_enu_nominals;

  for (int pi = 0; pi < sh->EventSplineMetadata.size(); ++pi) {
    auto const &[pname, pknots] = sh->EventSplineMetadata[pi];
    auto pdir = spdir->mkdir(pname.c_str());
    for (int im = 0; im < modes.GetNModes(); ++im) {

      sh->UnsetWeights();
      sh->SelectMode(im);
      sh->Reweight();
      if (sh->GetMCHist(0)->Integral() == 0) {
        continue;
      }
      auto mpdir = pdir->mkdir(modes.GetMaCh3ModeName(im).c_str());
      if (!pi) {
        auto h = sh->GetMCHist(0);
        nomdir->WriteTObject(h, fmt::format("nominal_mode{}", im).c_str());
      }

      for (int ib = 0; ib < (enu_bins.size() - 1); ++ib) {

        sh->UnsetWeights();
        sh->SelectTrueEnu(enu_bins[ib], enu_bins[ib + 1]);
        sh->SelectMode(im);
        sh->Reweight();
        if (sh->GetMCHist(0)->Integral() == 0) {
          continue;
        }
        if (!pi) {
          auto h = sh->GetMCHist(0);
          auto hname = fmt::format("nominal_mode{}_enu{}", im, ib);
          nomdir->WriteTObject(h, hname.c_str());
          mode_enu_nominals[im][ib] =
              std::unique_ptr<TH1>(static_cast<TH1 *>(h->Clone(hname.c_str())));
          mode_enu_nominals[im][ib]->SetDirectory(nullptr);
        }

        std::vector<TGraph> m3_splines;

        for (int ki = 0; ki < pknots.size(); ++ki) {

          sh->SetEventResponseWeights(pi, ki);
          sh->Reweight();
          auto h = sh->GetMCHist(0);
          if (h->Integral() != 0) {
            fmt::print("\rWriting spline for param {} {}/{}, for mode: {}, enu "
                       "bin: {:.3f} -- {:.3f}, knot: {}/{}",
                       pname, pi + 1, sh->EventSplineMetadata.size(),
                       modes.GetMaCh3ModeName(im), enu_bins[ib],
                       enu_bins[ib + 1], ki + 1, pknots.size());
            auto hname = fmt::format("mode{}_enu{}_{}_ki{}", im, ib, pname, ki);
            mpdir->WriteTObject(h, hname.c_str());

            switch (ndims) {
            case 1: {
              for (int xi = 0; xi < h->GetXaxis()->GetNbins(); ++xi) {
                if (m3_splines.size() <= xi) {
                  m3_splines.emplace_back();
                  auto gname =
                      fmt::format("{}_mode_{}_enu_{}_x_{}", pname, im, ib, xi);
                  m3_splines.back().SetName(gname.c_str());
                }
                m3_splines[xi].SetPoint(
                    ki, pknots[ki],
                    h->GetBinContent(xi + 1) /
                        mode_enu_nominals[im][ib]->GetBinContent(xi + 1));
              }
            }
            case 2: {
              auto h2 = static_cast<TH2 const *>(h);
            }
            }
          }
        }

        for (auto const &g : m3_splines) {
          mpdir->WriteTObject(&g, g.GetName());
        }
      }
    }
    fmt::print("\n");
    break;
  }
  fout.Close();
}
