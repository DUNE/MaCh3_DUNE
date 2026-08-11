#include "Samples/BeamFDStandardRecord/ReadEvents.h"

#include "Manager/Manager.h"

_MaCh3_Safe_Include_Start_ //{
#include "duneanaobj/StandardRecord/StandardRecord.h"
_MaCh3_Safe_Include_End_ //}

#include "duneanasel/common/TruthTools.h"
#include "duneanasel/fd/beam/Observables.h"
#include "duneanasel/fd/beam/Selections.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <cmath>
#include <random>

namespace dune::beamfd {

  float GetPOT(TTree & tree) {
    double pot = 0;
    TTreeReader metardr(&tree);
    TTreeReaderValue<double> entry_pot(metardr, "pot");
    while (metardr.Next()) {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnull-dereference"
      pot += *entry_pot;
#pragma GCC diagnostic pop
    }
    return float(pot);
  }

  std::vector<EventInfo> ReadEvents(TTree & tree, float downsamplefraction) {

    std::random_device r;
    std::default_random_engine e1(r());
    std::uniform_real_distribution<> dis(0, 1);

    // Reco Variables
    TTreeReader caf_reader(&tree);

    TTreeReaderValue<caf::StandardRecord> sr(caf_reader, "rec");

    std::vector<EventInfo> events;

    size_t sr_it = 0;
    while (caf_reader.Next()) {

      if (dis(e1) < downsamplefraction) {
        continue;
      }

      if (!sr->common.ixn.pandora.size()) { // empty event
        continue;
      }

      auto const &ixn = sr->common.ixn.pandora[0];
      EventInfo ev;

      auto sample = sel::beam::FD1::numode::ApplySelection(ixn);
      ev.reco.sample = sample;
      ev.reco.e_nu = proj::beam::FD1::ENuReco(ixn, sample);
      auto const &det_pos_cm = proj::beam::det_pos_cm(ixn);
      ev.reco.vtx_pos_cm = {det_pos_cm.x, det_pos_cm.y, det_pos_cm.z};

      auto const &tixn = sr->mc.nu.front();

      ev.truth.generator_mode = tixn.mode;
      ev.truth.is_cc = tixn.iscc;
      ev.truth.tgt_a = 40; // should do properly

      ev.truth.nu.pdg = tixn.pdg;
      ev.truth.nu.pdg_unosc = tixn.pdgorig;
      ev.truth.nu.e = tixn.E;

      if (ev.reco.sample > 0) {
        events.push_back(ev);
      }

      sr_it++;
    }
    return events;
  }
} // namespace dune::beamfd
