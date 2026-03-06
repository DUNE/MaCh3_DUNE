#pragma once

#include "Samples/BeamOffAxis/EventInfo.h"

#include <utility>
#include <vector>

namespace dune::beamoffaxis {

enum EDetSyst {
  kEMEnergyResolution = 0,
  kMuonEnergyScale,
  kMuonEnergyResolution,
  kNeutronEnergyResolution,
  kChargedHadronEnergyResolution,
  kNDetSysts
};

void EMEnergyResolution(const double *par_val, EventInfo &ev);
void MuonEnergyScale(const double *par_val, EventInfo &ev);
void MuonEnergyResolution(const double *par_val, EventInfo &ev) ;
void NeutronEnergyResolution(const double *par_val, EventInfo &ev);
void ChargedHadronEnergyResolution(const double *par_val, EventInfo &ev);

inline auto GetFluxFocussingWeighter(int focussing_par) {
  return [=](double const *par_val, EventInfo const &ev) {
    return 1 + (*par_val) * ev.syst.flux.focussing_ratio[focussing_par];
  };
}
inline auto GetFluxHadronProductionWeighter(int hadprod_par) {
  return [=](double const *par_val, EventInfo const &ev) {
    return 1 + (*par_val) * ev.syst.flux.hadprod_ratio[hadprod_par];
  };
}

std::pair<std::vector<float>, std::vector<float>>
GetFluxVariationRatios(int nu_pdg, double enu_true_GeV, double off_axis_pos_m,
                       bool is_numode);

void PrintFluxParameterNames();

} // namespace dune::beamoffaxis
