#pragma once

#include "Samples/BeamFDStandardRecord/EventInfo.h"

#include "TTree.h"

namespace dune::beamfd {

float GetPOT(TTree &);
std::vector<EventInfo> ReadEvents(TTree &, float downsamplefraction = 0);

} // namespace dune::beamfd
