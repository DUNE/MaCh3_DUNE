#include "Samples/BeamFDStandardRecord/ReadEvents.h"

#include "Manager/Manager.h"

_MaCh3_Safe_Include_Start_ //{
#include "duneanaobj/StandardRecord/StandardRecord.h"
_MaCh3_Safe_Include_End_ //}

#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include <cmath>

namespace dune::beamfd {

  double GetPOT(TTree &) {
    double pot = 0;
    // TTreeReader metardr(&tree);
    // TTreeReaderValue<double> entry_pot(metardr, "pot");
    // while (metardr.Next()) {
    //   pot += *entry_pot;
    // }
    return pot;
  }

  std::vector<EventInfo> ReadEvents(TTree & tree) {

    // Reco Variables
    TTreeReader caf_reader(&tree);

    TTreeReaderValue<caf::StandardRecord> sr(caf_reader, "rec");

    std::vector<EventInfo> events(tree.GetEntries());

    size_t ev_it = 0;
    while (caf_reader.Next()) {

      auto &ev = events[ev_it++];
      (void)ev;

      // ev.truth.is_cc = *isCC;

      // ev.truth.nu.pdg = *nuPDG;
      // ev.truth.nu.pdg_unosc = *nuPDGunosc;
      // ev.truth.nu.e = *Ev;

      // ev.reco.enu = std::max(0.0, *Ev_reco);
      // ev.reco.e_lep = std::max(0.0, *Elep_reco);
      // ev.reco.e_had = std::max(0.0, ev.reco.enu - ev.reco.e_lep);
    }
    return events;
  }
} // namespace dune::beamfd
