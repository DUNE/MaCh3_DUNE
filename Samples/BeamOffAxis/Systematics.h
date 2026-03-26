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

} // namespace dune::beamoffaxis
