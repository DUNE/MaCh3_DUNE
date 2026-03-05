#pragma once

#include "Samples/BeamOffAxis/EventInfo.h"

#include "TTree.h"

namespace dune::beamoffaxis {

double GetPOT(TTree &);
std::vector<EventInfo> ReadEvents(TTree &);

} // namespace dune::beamoffaxis
