#include "Systematics/Beam.h"

#include "duneanafluxtools/FluxWeighter.h"

namespace dune::syst {
std::pair<std::vector<float>, std::vector<float>>
GetFluxVariationWeights(int nu_pdg, double enu_true_GeV, bool is_FD,
                        bool is_numode, double off_axis_pos_m) {

  auto const &fw = FluxWeighter::Get();
  auto nucfg = fw.GetNuConfig(nu_pdg, !is_FD, is_numode);

  auto flux_focussing_systbin =
      fw.GetFocussingBin(nu_pdg, enu_true_GeV, off_axis_pos_m, nucfg);

  auto flux_hadprod_systbin =
      fw.GetFocussingBin(nu_pdg, enu_true_GeV, off_axis_pos_m, nucfg);

  std::pair<std::vector<float>, std::vector<float>> weights;
  for (size_t i = 0; i < fw.GetNFocussingParams(); i++) {
    weights.first.push_back(
        float(fw.GetFluxFocussingWeight(i, 1, flux_focussing_systbin, nucfg)));
  }
  for (size_t i = 0; i < fw.GetNHadProdPCAComponents(); i++) {
    weights.second.push_back(
        float(fw.GetFluxHadProdWeight(i, 1, flux_hadprod_systbin, nucfg)));
  }
  return weights;
}

std::vector<std::string> GetFluxFocussingParamNames() {
  auto const &fw = FluxWeighter::Get();
  size_t nfocus_par = fw.GetNFocussingParams();
  std::vector<std::string> focussing_par_names;
  for (size_t focussing_par = 0; focussing_par < nfocus_par; focussing_par++) {
    focussing_par_names.push_back(fw.GetFocussingParamName(focussing_par));
  }
  return focussing_par_names;
}

std::vector<std::string> GetFluxHadProdParamNames() {
  auto const &fw = FluxWeighter::Get();
  size_t Nhadprod_par = fw.GetNHadProdPCAComponents();
  std::vector<std::string> hadprod_par_names;
  for (size_t hadprod_par = 0; hadprod_par < Nhadprod_par; hadprod_par++) {
    hadprod_par_names.push_back("Flux_HadProd_Param_" +
                                std::to_string(hadprod_par));
  }
  return hadprod_par_names;
}

} // namespace dune::syst
