#include "Samples/SampleHandlerBeamWinter26/ReadEvents.h"

#include "Manager/Manager.h"

_MaCh3_Safe_Include_Start_ //{
#include "duneanaobj/StandardRecord/StandardRecord.h"
    _MaCh3_Safe_Include_End_ //}

#include "duneanasel/skimmers/MaCh3Winter26SkimModels.h"

#include "duneanasel/common/TruthTools.h"
#include "duneanasel/fd/beam/Observables.h"
#include "duneanasel/fd/beam/Selections.h"

#include "ROOT/RNTuple.hxx"
#include "ROOT/RNTupleReader.hxx"

#include <cmath>
#include <random>

    namespace dune::beam {

  float GetPOT(TDirectory & dir) {

    mach3winter26::Meta meta;

    auto rdr = ROOT::RNTupleReader::Open(std::move(meta.model),
                                         *dir.Get<ROOT::RNTuple>("meta"));

    rdr->LoadEntry(0);

    return float(*meta.pot_processed);
  }

  std::pair<std::vector<EventInfo>, size_t> ReadEvents(
      TDirectory & dir, float downsamplefraction) {

    std::random_device r;
    std::default_random_engine e1(r());
    std::uniform_real_distribution<> dis(0, 1);

    mach3winter26::Event event;

    auto rdr = ROOT::RNTupleReader::Open(std::move(event.model),
                                         *dir.Get<ROOT::RNTuple>("event"));

    std::vector<EventInfo> events;

    for (auto entryId : *rdr) {
      rdr->LoadEntry(entryId);

      if (dis(e1) < downsamplefraction) {
        continue;
      }

      EventInfo ev;

      ev.reco.sample = *event.reco.sample;
      ev.reco.e_nu = *event.reco.enu_GeV;
      ev.reco.vtx_pos_cm = *event.reco.vtx_pos_cm;

      ev.truth.generator_mode = *event.truth.mode;
      ev.truth.is_cc = *event.truth.is_cc;
      ev.truth.tgt_a = 40; // should do properly

      ev.truth.nu.pdg = *event.truth.pdgnu;
      ev.truth.nu.pdg_unosc = *event.truth.pdgnu_original;
      ev.truth.nu.e = *event.truth.enu_GeV;

      events.push_back(ev);
    }
    return {events, rdr->GetNEntries()};
  }
} // namespace dune::beam
