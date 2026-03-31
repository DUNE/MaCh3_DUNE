#pragma once

#include "Samples/BeamOffAxis/EventInfo.h"

#include <string>
#include <utility>
#include <vector>

namespace dune::beamoffaxis {

void EnergyScales(std::vector<double> const &par_vals, EventInfo &ev);
void ParticleEnergyResolutions(std::vector<double> const &par_vals,
                               EventInfo &ev);
void CalculateVariedCompositeQuantities(EventInfo &ev);

std::pair<std::vector<float>, std::vector<float>>
GetFluxVariationRatios(int nu_pdg, double enu_true_GeV, double off_axis_pos_m,
                       bool is_numode);

std::vector<std::string> GetFluxFocussingParamNames();
std::vector<std::string> GetFluxHadProdParamNames();

void UpdateFluxFocussingWeight(std::vector<double> const &par_vals,
                               EventInfo &ev);
void UpdateFluxHadProdWeight(std::vector<double> const &par_vals,
                             EventInfo &ev);
void PrintFluxParameterNames();

void MissingProtonFD(double const &par_val, EventInfo &ev) {
  ev.varied_truth.enurec_hadavailable_missed -= par_val * ev.truth.had.e_proton;

  ev.varied_reco.enu -= par_val * ev.truth.had.e_proton;
  ev.varied_reco.e_had -= par_val * ev.truth.had.e_proton;
  ev.varied_reco.e_proton -= par_val * ev.truth.had.e_proton;
}

} // namespace dune::beamoffaxis
