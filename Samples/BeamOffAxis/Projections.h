#pragma once

#include "Samples/BeamOffAxis/EventInfo.h"

#include <string>
#include <unordered_map>

namespace dune::beamoffaxis {

// below is ugly, but lets us define it only once and get the enum and both
// maps https://en.wikipedia.org/wiki/X_macro
#define LIST_OF_VARIABLES                                                      \
  X(TrueNeutrinoEnergy)                                                        \
  X(RecoNeutrinoEnergy)                                                        \
  X(TrueXPos)                                                                  \
  X(TrueYPos)                                                                  \
  X(TrueZPos)                                                                  \
  X(TrueW)                                                                     \
  X(Mode)                                                                      \
  X(OscillationChannel)                                                        \
  X(ELepRec)                                                                   \
  X(Enubias)                                                                   \
  X(IsCC)                                                                      \
  X(OffAxisPosition)                                                           \
  X(ENuRes)                                                                    \
  X(ELepRes)                                                                   \
  X(EHadRes)                                                                   \
  X(EEMRes)                                                                    \
  X(EChgHadRes)                                                                \
  X(ENeutronRes)
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
  case kMode:
    return &ev.truth.mach3_mode;
  case kOscillationChannel:
    return &ev.truth.mach3_mode;
  case kIsCC:
    return &ev.truth.is_cc;

  case kTrueNeutrinoEnergy:
    return &ev.truth.nu.e;

  case kTrueXPos:
    return &ev.truth.vtx.detcoords_cm[0];
  case kTrueYPos:
    return &ev.truth.vtx.detcoords_cm[1];
  case kTrueZPos:
    return &ev.truth.vtx.detcoords_cm[2];
  case kOffAxisPosition:
    return &ev.truth.vtx.off_axis_pos_m;

  case kTrueW:
    return &ev.truth.kine.invariant_mass;
  case kEnubias:
    return &ev.truth.kine.enurec_hadavailable_missed;

  case kRecoNeutrinoEnergy:
    return &ev.varied_reco.enu;
  case kELepRec:
    return &ev.varied_reco.e_lep;

  case kENuRes:
    return &ev.varied_res.enu;
  case kELepRes:
    return &ev.varied_res.e_lep;
  case kEHadRes:
    return &ev.varied_res.e_had;
  case kEEMRes:
    return &ev.varied_res.e_EM;
  case kEChgHadRes:
    return &ev.varied_res.e_ChgHad;
  case kENeutronRes:
    return &ev.varied_res.e_neutron;
  default:
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type...");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

} // namespace dune::beamoffaxis
