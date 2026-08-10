#pragma once

#include "Samples/BeamFDStandardRecord/EventInfo.h"

#include "TTree.h"

namespace dune::beamfd {

double GetPOT(TTree &);
std::vector<EventInfo> ReadEvents(TTree &);

} // namespace dune::beamfd
