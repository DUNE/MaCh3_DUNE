#pragma once

#include "TH2D.h"

#include <vector>
#include <memory>
#include <string>

double ERecQE(int nupdg, bool isCC, double el_MeV, double theta_lep_rad);

namespace dune::beamoffaxis {

class SampleHandlerBeamOffAxis;

std::vector<std::vector<std::vector<std::vector<std::unique_ptr<TH2D>>>>>
GetBinnedWeights(SampleHandlerBeamOffAxis &sample, int iSubSample,
                 std::vector<std::string> ParamNames,
                 std::vector<std::vector<int>> ParamModes,
                 std::vector<double> TrueEBins);

} // namespace dune::beamoffaxis
