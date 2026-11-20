#pragma once

#include "Samples/SampleHandlerBeamWinter26/EventInfo.h"

#include "TDirectory.h"

namespace dune::beam {

float GetPOT(TDirectory & dir);
std::pair<std::vector<EventInfo>, size_t>
ReadEvents(TDirectory & dir, float downsamplefraction = 0);

} // namespace dune::beam
