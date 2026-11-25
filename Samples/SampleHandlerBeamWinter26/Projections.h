#pragma once

#include "Samples/SampleHandlerBeamWinter26/EventInfo.h"

#include "Manager/MaCh3Exception.h"
#include "Manager/MaCh3Logger.h"

#include <string>
#include <unordered_map>

namespace dune::beam {

// below is ugly, but lets us define it only once and get the enum and both
// maps https://en.wikipedia.org/wiki/X_macro
#define LIST_OF_VARIABLES                                                      \
  X(TrueNeutrinoEnergy)                                                        \
  X(TargetNucleus)                                                             \
  X(IsCC)                                                                      \
  X(Mode)                                                                      \
  X(RecoSample)                                                                \
  X(RecoNeutrinoEnergy)                                                        \
  X(RecoLeptonEnergy)                                                          \
  X(OscillationChannel)

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
  static double dummy = 0;
  switch (KinPar) {
  case kTrueNeutrinoEnergy:
    return &ev.truth.nu.e_GeV;
  case kTargetNucleus:
    return &ev.truth.tgt_a;
  case kIsCC:
    return &ev.truth.is_cc;
  case kMode:
    return &ev.truth.mach3_mode;

  case kRecoSample:
    return &ev.reco.sample;
  case kRecoNeutrinoEnergy:
    return &ev.reco.e_nu_GeV;
  case kRecoLeptonEnergy:
    return &ev.reco.e_mu_GeV;

  case kOscillationChannel:
    return &dummy;

  default:
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type...");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

} // namespace dune::beam
