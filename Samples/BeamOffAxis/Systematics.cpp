#include "Samples/BeamOffAxis/Systematics.h"

#include "Systematics/Flux/OffAxisFluxUncertaintyHelper.h"

namespace dune::beamoffaxis {

void EMEnergyResolution(const double *par_val, EventInfo &ev) {
  ev.varied_reco.e_pi0 = std::max(
      0.0, ev.reco.e_pi0 + (*par_val) * (ev.truth.had.e_pi0 - ev.reco.e_pi0));

  if (std::abs(ev.truth.lep.pdg) == 11) {
    ev.varied_reco.e_lep = std::max(
        0.0, ev.reco.e_lep + (*par_val) * (ev.truth.lep.e - ev.reco.e_lep));
  }
}

void MuonEnergyScale(const double *par_val, EventInfo &ev) {
  if (std::abs(ev.truth.lep.pdg) == 13) {
    ev.varied_reco.e_lep = std::max(0.0, ev.reco.e_lep * (1 + *par_val));
  }
}

void MuonEnergyResolution(const double *par_val, EventInfo &ev) {
  if (std::abs(ev.truth.lep.pdg) == 13) {
    ev.varied_reco.e_lep = std::max(
        0.0, ev.reco.e_lep + (*par_val) * (ev.truth.lep.e - ev.reco.e_lep));
  }
}

void NeutronEnergyResolution(const double *par_val, EventInfo &ev) {
  ev.varied_reco.e_neutron =
      std::max(0.0, ev.reco.e_neutron + (*par_val) * (ev.truth.had.e_neutron -
                                                    ev.reco.e_neutron));
}

void ChargedHadronEnergyResolution(const double *par_val, EventInfo &ev) {
  ev.varied_reco.e_proton =
      std::max(0.0, ev.reco.e_proton +
                      (*par_val) * (ev.truth.had.e_proton - ev.reco.e_proton));
  ev.varied_reco.e_piplus =
      std::max(0.0, ev.reco.e_piplus +
                      (*par_val) * (ev.truth.had.e_piplus - ev.reco.e_piplus));
  ev.varied_reco.e_piminus =
      std::max(0.0, ev.reco.e_piminus + (*par_val) * (ev.truth.had.e_piminus -
                                                    ev.reco.e_piminus));
}

std::pair<std::vector<float>, std::vector<float>>
GetFluxVariationRatios(int nu_pdg, double enu_true_GeV, double off_axis_pos_m,
                       bool is_numode) {
  auto nucfg =
      OffAxisFluxUncertaintyHelper::Get().GetNuConfig(nu_pdg, true, is_numode);

  auto flux_focussing_systbin =
      OffAxisFluxUncertaintyHelper::Get().GetFocussingBin(
          nu_pdg, enu_true_GeV, off_axis_pos_m, nucfg);

  auto flux_hadprod_systbin =
      OffAxisFluxUncertaintyHelper::Get().GetFocussingBin(
          nu_pdg, enu_true_GeV, off_axis_pos_m, nucfg);

  std::pair<std::vector<float>, std::vector<float>> ratios;
  for (size_t i = 0;
       i < OffAxisFluxUncertaintyHelper::Get().GetNFocussingParams(); i++) {
    ratios.first.push_back(
        float(OffAxisFluxUncertaintyHelper::Get().GetFluxFocussingRatio(
            i, flux_focussing_systbin, nucfg)));
  }
  for (size_t i = 0;
       i < OffAxisFluxUncertaintyHelper::Get().GetNHadProdPCAComponents();
       i++) {
    ratios.second.push_back(
        float(OffAxisFluxUncertaintyHelper::Get().GetFluxHadProdRatio(
            i, flux_hadprod_systbin, nucfg)));
  }
  return ratios;
}

void PrintFluxParameterNames() {

  std::cout << "=== Focusing Flux Params ===\n";
  for (size_t i = 0;
       i < OffAxisFluxUncertaintyHelper::Get().GetNFocussingParams(); i++) {
    std::cout << OffAxisFluxUncertaintyHelper::Get().GetFocussingParamName(i)
              << "\n";
  }

  std::cout << "=== Hadron-Production PCA Params ===\n";
  for (size_t i = 0;
       i < OffAxisFluxUncertaintyHelper::Get().GetNHadProdPCAComponents();
       i++) {
    std::cout << "Flux_HadProd_Param_" << i << "\n";
  }
}

} // namespace dune::beamoffaxis
