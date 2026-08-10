#pragma once

#include <array>
#include <vector>
#include <cstddef>

namespace dune::beamfd {

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

    enum ESample { kRejected = 0, kNuMuCCLike, kNuECCLike, kNCLike };
    double sample;

    double e_nu;

    std::array<double, 3> vtx_pos_cm;

  } reco;
};

struct EventInfo : public CAFEventInfo {

  size_t tag_id;
  bool is_numode;
  int sample;

  struct SystInfo {
    struct Flux {
      std::vector<float> focussing_weights;
      std::vector<float> hadprod_weights;

      double total_weight;
    } flux;

  } syst;

  struct Weights {
    double pot;
  } weights;
};

} // namespace dune::beamfd
