#pragma once

#include <utility>
#include <vector>

namespace dune::beamoffaxis {

enum FuncParEnum {
  kTotalEScaleND = 0,
  kTotalEScaleND_mu,
  kEMResND,
  kEScaleMuSpectND,
  kMuonRes_ND,
  kNRes_ND,
  kHadRes_ND,
  kNFuncPars
};

std::pair<std::vector<float>, std::vector<float>>
GetFluxVariationRatios(int nu_pdg, double enu_true_GeV, double off_axis_pos_m,
                       bool is_numode);

void PrintFluxParameterNames();

} // namespace dune::beamoffaxis
