#pragma once

#include "Samples/BeamFDStandardRecord/EventInfo.h"

#include "Manager/MaCh3Logger.h"
#include "Manager/MaCh3Exception.h"

#include <string>
#include <unordered_map>

namespace dune::beamfd {

// below is ugly, but lets us define it only once and get the enum and both
// maps https://en.wikipedia.org/wiki/X_macro
#define LIST_OF_VARIABLES                                                      \
  X(TrueNeutrinoEnergy)                                                        \
  X(OscillationChannel)                                                        \
  X(IsCC)                                                                      \
  X(Mode)                                                                      \
  X(ELepRec)                                                                   \
  X(EHadRec)                                                                   \
  X(RecoNeutrinoEnergy)

#define X(a) k##a,

/// @brief Enum to identify kinematics
enum KinematicTypes { LIST_OF_VARIABLES };

#undef X
#define X(a) {#a, k##a},
const std::unordered_map<std::string, int> KinematicParametersDUNE = {
    LIST_OF_VARIABLES};

#undef X
#define X(a) {k##a, #a},
const std::unordered_map<int, std::string> ReversedKinematicParametersDUNE = {
    LIST_OF_VARIABLES};

#undef X
#undef LIST_OF_VARIABLES

inline const double *ResolveKinematicEventMember(KinematicTypes KinPar,
                                                 EventInfo const &ev) {
  switch (KinPar) {
  case kTrueNeutrinoEnergy:
    return &ev.truth.nu.e;
  case kOscillationChannel:
    return &ev.truth.mach3_mode;
  case kIsCC:
    return &ev.truth.is_cc;
  case kMode:
    return &ev.truth.mach3_mode;

  case kELepRec:
    return &ev.reco.e_lep;
  case kEHadRec:
    return &ev.reco.e_had;
  case kRecoNeutrinoEnergy:
    return &ev.reco.enu;

  default:
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type...");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

} // namespace dune::beamfd
