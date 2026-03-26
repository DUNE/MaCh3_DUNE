#include "Samples/BeamOffAxis/Systematics.h"

#include "Systematics/Flux/OffAxisFluxUncertaintyHelper.h"

#include <cmath>

namespace dune::beamoffaxis {

inline double EnergyScaleVariation(double const *par_vals, double e,
                                   double sqrte) {
  double e_prime =
      e * ((1 + par_vals[0]) + par_vals[1] * sqrte) + (par_vals[2] * sqrte);
  // std::cout << "e = " << e << ", e_prime = " << e_prime << std::endl;
  return (e_prime < 0) ? 0 : e_prime - e;
}

void EnergyScales(std::vector<double> const &par_vals, EventInfo &ev) {
  auto MuonContained_vals = par_vals.data();
  auto MuonTracked_vals = MuonContained_vals + 3;
  auto EM_vals = MuonTracked_vals + 3;
  auto ChgHad_vals = EM_vals + 3;
  auto Neutron_vals = ChgHad_vals + 3;
  auto Total_vals = Neutron_vals + 3;

  if (std::abs(ev.truth.lep.pdg) == 13) {
    if (ev.reco.muonlike_contained) {
      ev.varied_reco.e_lep += EnergyScaleVariation(
          MuonContained_vals, ev.reco.e_lep, ev.syst.sqrt_e.lep);
    } else if (ev.reco.muonlike_tracker) {
      ev.varied_reco.e_lep += EnergyScaleVariation(
          MuonTracked_vals, ev.reco.e_lep, ev.syst.sqrt_e.lep);
    }
  }

  ev.varied_reco.e_pi0 +=
      EnergyScaleVariation(EM_vals, ev.reco.e_pi0, ev.syst.sqrt_e.pi0);
  if (std::abs(ev.truth.lep.pdg) == 11) {
    ev.varied_reco.e_lep +=
        EnergyScaleVariation(EM_vals, ev.reco.e_lep, ev.syst.sqrt_e.lep);
  }

  ev.varied_reco.e_proton += EnergyScaleVariation(ChgHad_vals, ev.reco.e_proton,
                                                  ev.syst.sqrt_e.proton);
  ev.varied_reco.e_piplus += EnergyScaleVariation(ChgHad_vals, ev.reco.e_piplus,
                                                  ev.syst.sqrt_e.piplus);
  ev.varied_reco.e_piminus += EnergyScaleVariation(
      ChgHad_vals, ev.reco.e_piminus, ev.syst.sqrt_e.piminus);

  ev.varied_reco.e_neutron += EnergyScaleVariation(
      Neutron_vals, ev.reco.e_neutron, ev.syst.sqrt_e.neutron);

  if (std::abs(ev.truth.lep.pdg) == 11) {
    ev.varied_reco.e_lep +=
        EnergyScaleVariation(Total_vals, ev.reco.e_lep, ev.syst.sqrt_e.lep);
  }
  ev.varied_reco.e_had +=
      EnergyScaleVariation(Total_vals, ev.reco.e_had, ev.syst.sqrt_e.had);
}

void ParticleEnergyResolutions(std::vector<double> const &par_vals,
                               EventInfo &ev) {
  auto const &Muon_val = par_vals[0];
  auto const &EM_val = par_vals[1];
  auto const &ChgHad_val = par_vals[2];
  auto const &Neutron_val = par_vals[3];

  if (std::abs(ev.truth.lep.pdg) == 13) {
    auto e_lep_bias = (ev.truth.lep.e - ev.reco.e_lep);
    ev.varied_reco.e_lep += Muon_val * e_lep_bias;
  }

  auto e_pi0_bias = (ev.truth.had.e_pi0 - ev.reco.e_pi0);
  ev.varied_reco.e_pi0 += EM_val * e_pi0_bias;

  if (std::abs(ev.truth.lep.pdg) == 11) {
    auto e_lep_bias = (ev.truth.lep.e - ev.reco.e_lep);
    ev.varied_reco.e_lep += EM_val * e_lep_bias;
  }

  auto e_proton_bias = (ev.truth.had.e_proton - ev.reco.e_proton);
  ev.varied_reco.e_proton += ChgHad_val * e_proton_bias;
  auto e_piplus_bias = (ev.truth.had.e_piplus - ev.reco.e_piplus);
  ev.varied_reco.e_piplus += ChgHad_val * e_piplus_bias;
  auto e_piminus_bias = (ev.truth.had.e_piminus - ev.reco.e_piminus);
  ev.varied_reco.e_piminus += ChgHad_val * e_piminus_bias;

  auto e_neutron_bias = (ev.truth.had.e_neutron - ev.reco.e_neutron);
  ev.varied_reco.e_neutron += Neutron_val * e_neutron_bias;
}

void CalculateVariedCompositeQuantities(EventInfo &ev) {
  auto const &reco = ev.reco;
  auto &varreco = ev.varied_reco;

  auto e_lep_shift = (varreco.e_lep < 0) ? 0 : (varreco.e_lep - reco.e_lep);
  auto e_had_shift = (varreco.e_had < 0) ? 0 : (varreco.e_had - reco.e_had);
  auto e_proton_shift =
      (varreco.e_proton < 0) ? 0 : (varreco.e_proton - reco.e_proton);
  auto e_piplus_shift =
      (varreco.e_piplus < 0) ? 0 : (varreco.e_piplus - reco.e_piplus);
  auto e_piminus_shift =
      (varreco.e_piminus < 0) ? 0 : (varreco.e_piminus - reco.e_piminus);
  auto e_neutron_shift =
      (varreco.e_neutron < 0) ? 0 : (varreco.e_neutron - reco.e_neutron);

  varreco.enu = reco.enu + e_lep_shift + e_had_shift + e_proton_shift +
                e_piplus_shift + e_piminus_shift + e_neutron_shift;
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

void UpdateFluxFocussingWeight(std::vector<double> const &par_vals,
                               EventInfo &ev) {
  for (int i = 0; i < par_vals.size(); ++i) {
    ev.syst.flux.total_weight *=
        1 + (par_vals[i] * ev.syst.flux.focussing_ratio[i]);
  }
}
void UpdateFluxHadProdWeight(std::vector<double> const &par_vals,
                             EventInfo &ev) {
  for (int i = 0; i < par_vals.size(); ++i) {
    ev.syst.flux.total_weight *=
        1 + (par_vals[i] * ev.syst.flux.hadprod_ratio[i]);
  }
}

std::vector<std::string> GetFluxFocussingParamNames() {
  size_t nfocus_par = OffAxisFluxUncertaintyHelper::Get().GetNFocussingParams();
  std::vector<std::string> focussing_par_names;
  for (size_t focussing_par = 0; focussing_par < nfocus_par; focussing_par++) {
    focussing_par_names.push_back(
        OffAxisFluxUncertaintyHelper::Get().GetFocussingParamName(
            focussing_par));
  }
  return focussing_par_names;
}
std::vector<std::string> GetFluxHadProdParamNames() {
  size_t Nhadprod_par =
      OffAxisFluxUncertaintyHelper::Get().GetNHadProdPCAComponents();
  std::vector<std::string> hadprod_par_names;
  for (size_t hadprod_par = 0; hadprod_par < Nhadprod_par; hadprod_par++) {
    hadprod_par_names.push_back("Flux_HadProd_Param_" +
                                std::to_string(hadprod_par));
  }
  return hadprod_par_names;
}

void PrintFluxParameterNames() {

  std::cout << "=== Focussing Flux Params ===\n";
  for (auto const &name : GetFluxFocussingParamNames()) {
    std::cout << name << "\n";
  }

  std::cout << "=== Hadron-Production PCA Params ===\n";
  for (auto const &name : GetFluxHadProdParamNames()) {
    std::cout << name << "\n";
  }
}

} // namespace dune::beamoffaxis
