#pragma once

#include <string>
#include <utility>
#include <vector>

namespace dune::syst {

std::vector<std::string> GetFluxFocussingParamNames();
std::vector<std::string> GetFluxHadProdParamNames();

std::pair<std::vector<float>, std::vector<float>>
GetFluxVariationWeights(int nu_pdg, double enu_true_GeV, bool is_FD,
                        bool is_numode, double off_axis_pos_m = 0);

} // namespace dune::syst
