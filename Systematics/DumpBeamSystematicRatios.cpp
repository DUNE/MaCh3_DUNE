#include "TFile.h"
#include "TH1.h"

#include "Systematics/Beam.h"

int main() {

  TFile out("BeamSystematicRatios.root", "RECREATE");

  for (auto const &spec : std::vector<std::pair<std::string, int>>{
           {"numu", 14}, {"nue", 12}, {"numubar", -14}, {"nuebar", -12}}) {
    std::vector<TH1D *> foc, hp;

    for (auto const pn : dune::syst::GetFluxFocussingParamNames()) {
      foc.push_back(
          new TH1D(Form("%s_numode_fd_%s", spec.first.c_str(), pn.c_str()), "",
                   1000, 0, 10));
    }

    for (auto const pn : dune::syst::GetFluxHadProdParamNames()) {
      hp.push_back(
          new TH1D(Form("%s_numode_fd_%s", spec.first.c_str(), pn.c_str()), "",
                   1000, 0, 10));
    }

    for (int i = 0; i < 1000; ++i) {
      double enu = (i + 0.5) * 10.0 / 1000.0;
      auto ws =
          dune::syst::GetFluxVariationWeights(spec.second, enu, true, true);
      for (int j = 0; j < ws.first.size(); ++j) {
        foc[j]->SetBinContent(i + 1, ws.first[j]);
      }
      for (int j = 0; j < ws.second.size(); ++j) {
        hp[j]->SetBinContent(i + 1, ws.second[j]);
      }
    }
  }
  out.Write();
  out.Close();
}
