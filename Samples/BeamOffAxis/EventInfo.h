#pragma once

#include <array>
#include <vector>

namespace dune::beamoffaxis {

struct CAFEventInfo {
  struct Truth {

    int target_a;

    int generator_mode;
    double mach3_mode;
    double is_cc;

    struct Neutrino {
      int pdg, pdg_unosc;
      double e;
    } nu;

    struct Lepton {
      int pdg;
      double e;
      double theta_nu;
    } lep;

    struct Hadron {

      double e_proton;
      double e_neutron;
      double e_piplus;
      double e_piminus;
      double e_pi0;

      double e_available;

      int n_proton;
      int n_neutron;
      int n_piplus;
      int n_piminus;
      int n_pi0;

    } had;

    struct Vertex {
      std::array<double, 3> detcoords_cm;
      double off_axis_pos_m;
    } vtx;

    struct Kinematics {
      double Q2;
      double invariant_mass;
      double bjorken_X;
      double bjorken_Y;

      double enurec_QE, enurec_hadav, enurec_hadavailable_missed;

    } kine;

  } truth;

  struct Reconstructed {
    int is_muonlike;
    int muonlike_contained;
    int muonlike_tracker;

    double enu;

    double e_lep;
    double e_had;

    double e_proton;
    double e_neutron;
    double e_pi0;
    double e_piplus;
    double e_piminus;
  } reco;
};

struct EventInfo : public CAFEventInfo {

  // A copy of the reconstructed info for systematically varying
  CAFEventInfo::Reconstructed varied_reco;

  int subsample;
  int is_numode;

  struct SystInfo {
    struct Flux {
      std::vector<float> focussing_ratio;
      std::vector<float> hadprod_ratio;

      double total_weight;
    } flux;

    struct {
      double lep;
      double had;

      double proton;
      double neutron;
      double pi0;
      double piplus;
      double piminus;
    } sqrt_e;
  } syst;

  struct Weights {
    double pot;
  } weights;
};

} // namespace dune::beamoffaxis
