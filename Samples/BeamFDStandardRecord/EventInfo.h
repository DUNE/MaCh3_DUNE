#pragma once

#include <array>
#include <vector>

namespace dune::beamfd {

struct CAFEventInfo {
  struct Truth {

    int generator_mode;
    double mach3_mode;
    double is_cc;

    struct Neutrino {
      int pdg, pdg_unosc;
      double e;
    } nu;

  } truth;

  struct Reconstructed {
    double enu;

    double e_lep;
    double e_had;

  } reco;
};

struct EventInfo : public CAFEventInfo {

  int subsample;
  int is_numode;

  struct SystInfo {
    struct Flux {
      std::vector<float> focussing_ratio;
      std::vector<float> hadprod_ratio;

      double total_weight;
    } flux;

  } syst;

  struct Weights {
    double pot;
  } weights;
};

} // namespace dune::beamfd
