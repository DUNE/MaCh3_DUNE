#pragma once

#include "Samples/SampleHandlerBeamWinter26/EventInfo.h"

#include "TDirectory.h"

namespace dune::beam {

float GetPOT(TDirectory &dir);

std::pair<std::vector<EventInfo>, size_t>
ReadEvents(TDirectory &dir, float downsamplefraction = 0);

std::pair<std::vector<EventInfoSplines>, size_t>
ReadEventsWithSplines(TDirectory &dir, float downsamplefraction = 0);

template <bool READ_EVENT_SPLINES>
inline auto ReadEvents(TDirectory &dir, float downsamplefraction = 0) {
  if constexpr (!READ_EVENT_SPLINES) {
    return ReadEvents(dir, downsamplefraction);
  } else {
    return ReadEventsWithSplines(dir, downsamplefraction);
  }
}

} // namespace dune::beam
