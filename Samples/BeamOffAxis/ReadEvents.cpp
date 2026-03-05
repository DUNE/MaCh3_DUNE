#include "Samples/BeamOffAxis/ReadEvents.h"

#include "Samples/BeamOffAxis/Utility.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <cmath>

namespace dune::beamoffaxis {

double GetPOT(TTree &tree) {
  double pot = 0;
  TTreeReader metardr(&tree);
  TTreeReaderValue<double> entry_pot(metardr, "pot");
  while (metardr.Next()) {
    pot += *entry_pot;
  }
  return pot;
}

std::vector<EventInfo> ReadEvents(TTree &tree) {

  // Reco Variables
  TTreeReader caf_reader(&tree);

  TTreeReaderValue<double> Ev_reco(caf_reader, "Ev_reco");
  TTreeReaderValue<double> Elep_reco(caf_reader, "Elep_reco");

  TTreeReaderValue<int> reco_numu(caf_reader, "reco_numu");
  TTreeReaderValue<int> muon_contained(caf_reader, "muon_contained");
  TTreeReaderValue<int> muon_tracker(caf_reader, "muon_tracker");

  TTreeReaderValue<double> eRecoP(caf_reader, "eRecoP");
  TTreeReaderValue<double> eRecoPip(caf_reader, "eRecoPip");
  TTreeReaderValue<double> eRecoPim(caf_reader, "eRecoPim");
  TTreeReaderValue<double> eRecoPi0(caf_reader, "eRecoPi0");
  TTreeReaderValue<double> eRecoN(caf_reader, "eRecoN");

  TTreeReaderValue<double> vtx_x(caf_reader, "vtx_x");
  TTreeReaderValue<double> vtx_y(caf_reader, "vtx_y");
  TTreeReaderValue<double> vtx_z(caf_reader, "vtx_z");
  TTreeReaderValue<double> det_x(caf_reader, "det_x");

  // Truth Variables
  TTreeReaderValue<int> mode(caf_reader, "mode");
  TTreeReaderValue<double> Ev(caf_reader, "Ev");

  TTreeReaderValue<double> LepE(caf_reader, "LepE");
  TTreeReaderValue<double> LepNuAngle(caf_reader, "LepNuAngle");
  TTreeReaderValue<double> LepMomX(caf_reader, "LepMomX");
  TTreeReaderValue<double> LepMomY(caf_reader, "LepMomY");
  TTreeReaderValue<double> LepMomZ(caf_reader, "LepMomZ");
  TTreeReaderValue<double> NuMomX(caf_reader, "NuMomX");
  TTreeReaderValue<double> NuMomY(caf_reader, "NuMomY");
  TTreeReaderValue<double> NuMomZ(caf_reader, "NuMomZ");
  TTreeReaderValue<double> eP(caf_reader, "eP");
  TTreeReaderValue<double> ePip(caf_reader, "ePip");
  TTreeReaderValue<double> ePim(caf_reader, "ePim");
  TTreeReaderValue<double> ePi0(caf_reader, "ePi0");
  TTreeReaderValue<double> eN(caf_reader, "eN");
  TTreeReaderValue<double> eOther(caf_reader, "eOther");

  TTreeReaderValue<double> Q2(caf_reader, "Q2");
  TTreeReaderValue<double> W(caf_reader, "W");
  TTreeReaderValue<double> X(caf_reader, "X");
  TTreeReaderValue<double> Y(caf_reader, "Y");

  TTreeReaderValue<int> nipi0(caf_reader, "nipi0");

  TTreeReaderValue<int> isCC(caf_reader, "isCC");

  TTreeReaderValue<int> nP(caf_reader, "nP");
  TTreeReaderValue<int> nN(caf_reader, "nN");
  TTreeReaderValue<int> nPi0(caf_reader, "nipi0");
  TTreeReaderValue<int> nPim(caf_reader, "nipim");
  TTreeReaderValue<int> nPip(caf_reader, "nipip");
  TTreeReaderValue<int> niem(caf_reader, "niem");

  TTreeReaderValue<int> nuPDG(caf_reader, "nuPDG");
  TTreeReaderValue<int> nuPDGunosc(caf_reader, "nuPDGunosc");
  TTreeReaderValue<int> LepPDG(caf_reader, "LepPDG");

  std::vector<EventInfo> events(tree.GetEntries());

  size_t ev_it = 0;
  while (caf_reader.Next()) {

    auto &ev = events[ev_it++];

    ev.truth.target_a = 40;

    ev.truth.generator_mode = *mode;
    ev.truth.is_cc = *isCC;

    ev.truth.nu.pdg = *nuPDG;
    ev.truth.nu.pdg_unosc = *nuPDGunosc;
    ev.truth.nu.e = *Ev;

    ev.truth.lep.pdg = *LepPDG;
    ev.truth.lep.e = *LepE;
    ev.truth.lep.theta_nu = *LepNuAngle;

    ev.truth.had.e_proton = *eP;
    ev.truth.had.e_neutron = *eN;
    ev.truth.had.e_piplus = *ePip;
    ev.truth.had.e_piminus = *ePim;
    ev.truth.had.e_pi0 = *ePi0;

    if (std::isfinite(*eP) && std::isfinite(*ePip) && std::isfinite(*ePim) &&
        std::isfinite(*ePi0) && std::isfinite(*eOther)) {
      ev.truth.had.e_available =
          *eP + *ePip + *ePim + *ePi0 + *eOther + (*nipi0) * 0.1349;
    } else {
      ev.truth.had.e_available = -999;
    }

    ev.truth.had.n_proton = *nP;
    ev.truth.had.n_neutron = *nN;
    ev.truth.had.n_piplus = *nPip;
    ev.truth.had.n_piminus = *nPim;
    ev.truth.had.n_pi0 = *nPi0;

    ev.truth.vtx.detcoords_cm = {*vtx_x, *vtx_y, *vtx_z};
    ev.truth.vtx.off_axis_pos_m =
        (ev.truth.vtx.detcoords_cm[0] + *vtx_x) / 100.0;

    ev.truth.kine.Q2 = *Q2;
    ev.truth.kine.invariant_mass = *W;

    if (!std::isnormal(*W)) {
      ev.truth.kine.invariant_mass = -999;
    }

    ev.truth.kine.bjorken_X = *X;
    ev.truth.kine.bjorken_Y = *Y;

    ev.truth.kine.enurec_QE = ERecQE(ev.truth.nu.pdg, ev.truth.is_cc,
                                     ev.truth.lep.e, ev.truth.lep.theta_nu);

    if (ev.truth.had.e_available >= 0) {
      ev.truth.kine.enurec_hadav = (ev.truth.lep.e + ev.truth.had.e_available);
      ev.truth.kine.enurec_hadavailable_missed =
          ev.truth.kine.enurec_hadav - ev.truth.nu.e;
    } else {
      ev.truth.kine.enurec_hadav = -999;
      ev.truth.kine.enurec_hadavailable_missed = -999;
    }

    ev.reco.is_muonlike = *reco_numu;
    ev.reco.muonlike_contained = *muon_contained;
    ev.reco.muonlike_tracker = *muon_tracker;

    ev.reco.enu = *Ev_reco;
    ev.reco.e_lep = *Elep_reco;

    ev.reco.e_proton = *eRecoP;
    ev.reco.e_neutron = *eRecoN;
    ev.reco.e_pi0 = *eRecoPip;
    ev.reco.e_piplus = *eRecoPim;
    ev.reco.e_piminus = *eRecoPi0;
  }

  return events;
}

} // namespace dune::beamoffaxis
