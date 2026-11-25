#pragma once

#include "Math/Point3D.h"

#include <array>
#include <cstddef>
#include <vector>

namespace dune::beam {

// for storing event information read out of inputs
struct CAFEventInfo {
  struct Truth {

    int generator_mode;
    double tgt_a;
    double mach3_mode;
    double is_cc;

    struct Neutrino {
      int pdg, pdg_unosc;
      double e;
    } nu;

  } truth;

  struct Reconstructed {

    double sample;

    double e_nu;

    ROOT::Math::XYZPoint vtx_pos_cm;

  } reco;
};

// for storing any additional derived information needed by M3
struct EventInfo : public CAFEventInfo {

  size_t tag_id;
  bool is_numode;
  int sample;

  struct SystInfo {
  } syst;

  struct Weights {
    double pot;
  } weights;
};

struct EventInfoSplines : public EventInfo {

  std::vector<std::vector<float>> knot_responses;

};

} // namespace dune::beam
