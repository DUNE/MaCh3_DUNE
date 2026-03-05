#include "Samples/BeamOffAxis/Systematics.h"

#include "Systematics/Flux/OffAxisFluxUncertaintyHelper.h"

namespace dune::beamoffaxis {

// void TotalEScaleND(const double *par, EventInfo &ev) {

//   double enu_rec = ev.reco.shift.enu;
//   double e_lep = reco.ELep;

//   if (std::abs(ev.LepPDG) == 11) {
//     ev.reco.shift.enu *= (1.0 + *par);
//     ev.reco.shift.e_lep *= (1.0 + *par);
//   }

//   erec += (*par) * reco.EHad;

//   shift.erec = std::max(0, erec);
//   shift.ELep = std::max(0, ELep);
// };

// RegisterIndividualFunctionalParameter(
//     "TotalEScaleND_mu", kTotalEScaleND_mu,
//     [this](const double *par, std::size_t iEvent) {
//       const auto &reco = DUNEMCEvents[iEvent].reco;
//       auto &shift = DUNEMCEvents[iEvent].shift;

//       if (std::abs(DUNEMCEvents[iEvent].LepPDG) == 13 &&
//           (reco.muon_contained || reco.muon_tracker)) {

//         double ELep = reco.ELep * (1.0 + *par);
//         double erec = DUNEMCEvents[iEvent].rw_erec + (*par) * reco.ELep;

//         if (ELep < 0.0 || erec < 0.0)
//           throw std::runtime_error("Negative muon energy");

//         shift.ELep = ELep;
//         shift.erec = erec;
//       }
//     });
// RegisterIndividualFunctionalParameter(
//     "EMResND", kEMResND, [this](const double *par, std::size_t iEvent) {
//       const auto &reco = DUNEMCEvents[iEvent].reco;
//       const auto &truth = DUNEMCEvents[iEvent].truth;
//       auto &shift = DUNEMCEvents[iEvent].shift;

//       double recoPi0 = reco.ePi0;
//       if (recoPi0 < 0.0)
//         recoPi0 = 0.0;

//       double ePi0 = recoPi0 + (*par) * (truth.ePi0 - recoPi0);
//       if (ePi0 < 0.0)
//         ePi0 = 0.0;

//       shift.ePi0 = ePi0;
//       shift.erec += ePi0 - recoPi0;

//       if (std::abs(DUNEMCEvents[iEvent].LepPDG) == 11) {
//         double ELep = reco.ELep + (*par) * (truth.LepE - reco.ELep);
//         if (ELep < 0.0)
//           ELep = 0.0;

//         shift.ELep = ELep;
//         shift.erec += ELep - reco.ELep;
//       }
//     });

// RegisterIndividualFunctionalParameter(
//     "EScaleMuSpectND", kEScaleMuSpectND,
//     [this](const double *par, std::size_t iEvent) {
//       const auto &reco = DUNEMCEvents[iEvent].reco;
//       auto &shift = DUNEMCEvents[iEvent].shift;

//       if (std::abs(DUNEMCEvents[iEvent].LepPDG) == 13 && reco.muon_tracker) {

//         double ELep = reco.ELep * (1.0 + *par);
//         double erec = DUNEMCEvents[iEvent].rw_erec + (*par) * reco.ELep;

//         if (ELep < 0.0 || erec < 0.0)
//           throw std::runtime_error("Negative muon energy");

//         shift.ELep = ELep;
//         shift.erec = erec;
//       }
//     });

// RegisterIndividualFunctionalParameter(
//     "MuonRes_ND", kMuonRes_ND, [this](const double *par, std::size_t iEvent) {
//       const auto &reco = DUNEMCEvents[iEvent].reco;
//       const auto &truth = DUNEMCEvents[iEvent].truth;
//       auto &shift = DUNEMCEvents[iEvent].shift;

//       if (std::abs(DUNEMCEvents[iEvent].LepPDG) == 13) {

//         double dE = (*par) * (truth.LepE - reco.ELep);
//         double ELep = reco.ELep + dE;
//         double erec = DUNEMCEvents[iEvent].rw_erec + dE;

//         if (ELep < 0.0)
//           ELep = 0.0;
//         if (erec < 0.0)
//           erec = 0.0;

//         shift.ELep = ELep;
//         shift.erec = erec;
//       }
//     });

// RegisterIndividualFunctionalParameter(
//     "NRes_ND", kNRes_ND, [this](const double *par, std::size_t iEvent) {
//       const auto &reco = DUNEMCEvents[iEvent].reco;
//       const auto &truth = DUNEMCEvents[iEvent].truth;
//       auto &shift = DUNEMCEvents[iEvent].shift;

//       double eN = reco.eN + (*par) * (truth.eN - reco.eN);

//       if (eN < 0.0) {
//         // Set shift to safe values
//         shift.eN = 0.0;
//         shift.erec = 0.0;
//       } else {
//         // Normal computation
//         shift.eN = eN;
//         shift.erec += eN - reco.eN;
//       }
//     });

// RegisterIndividualFunctionalParameter(
//     "HadRes_ND", kHadRes_ND, [this](const double *par, std::size_t iEvent) {
//       const auto &reco = DUNEMCEvents[iEvent].reco;
//       const auto &truth = DUNEMCEvents[iEvent].truth;
//       auto &shift = DUNEMCEvents[iEvent].shift;

//       double dE = (*par) * ((truth.ePim - reco.ePim) +
//                             (truth.ePip - reco.ePip) + (truth.eP - reco.eP));

//       shift.erec = DUNEMCEvents[iEvent].rw_erec + dE;
//       shift.eP = reco.eP + (*par) * (truth.eP - reco.eP);
//       shift.ePip = reco.ePip + (*par) * (truth.ePip - reco.ePip);
//       shift.ePim = reco.ePim + (*par) * (truth.ePim - reco.ePim);
//     });

// /// Fake Data Syst
// RegisterIndividualFunctionalParameter(
//     "NuWroFakeDataWeight", kNuWroFakeDataWeight,
//     [this](const double *par, std::size_t iEvent) {
//       DUNEMCEvents[iEvent].flux_w *=
//           (((*par) * (this->NuWroFakeDataWeight(iEvent) - 1.0)) + 1);
//     });

// // don't register flux parameters if we're not using them.
// if (ParHandler->GetNumParFromGroup("Flux")) {

//   for (size_t par_it = 0;
//        par_it < OffAxisFluxUncertaintyHelper::Get().GetNFocussingParams();
//        par_it++) {
//     RegisterIndividualFunctionalParameter(
//         OffAxisFluxUncertaintyHelper::Get().GetFocussingParamName(par_it),
//         int(kNFuncPars + par_it),
//         [this, par_it](const double *par, std::size_t iEvent) {
//           if (DUNEMCEvents[iEvent].syst.flux_focussing_ratio.size() <= par_it) {
//             return;
//           }
//           DUNEMCEvents[iEvent].flux_w *=
//               (1 +
//                (*par) * DUNEMCEvents[iEvent].syst.flux_focussing_ratio[par_it]);
//         });
//   }

//   for (size_t par_it = 0;
//        par_it < OffAxisFluxUncertaintyHelper::Get().GetNHadProdPCAComponents();
//        par_it++) {
//     RegisterIndividualFunctionalParameter(
//         "Flux_HadProd_Param_" + std::to_string(par_it),
//         int(kNFuncPars +
//             OffAxisFluxUncertaintyHelper::Get().GetNFocussingParams() + par_it),
//         [this, par_it](const double *par, std::size_t iEvent) {
//           if (DUNEMCEvents[iEvent].syst.flux_hadprod_ratio.size() <= par_it) {
//             return;
//           }
//           DUNEMCEvents[iEvent].flux_w *=
//               (1 +
//                (*par) * DUNEMCEvents[iEvent].syst.flux_hadprod_ratio[par_it]);
//         });
//   }
// }

std::pair<std::vector<float>, std::vector<float>>
GetFluxVariationRatios(int nu_pdg, double enu_true_GeV, double off_axis_pos_m,
                       bool is_numode) {
  auto nucfg =
      OffAxisFluxUncertaintyHelper::Get().GetNuConfig(nu_pdg, true, is_numode);

  auto flux_focussing_systbin =
      OffAxisFluxUncertaintyHelper::Get().GetFocussingBin(
          nu_pdg, enu_true_GeV, off_axis_pos_m, nucfg);

  auto flux_hadprod_systbin =
      OffAxisFluxUncertaintyHelper::Get().GetFocussingBin(
          nu_pdg, enu_true_GeV, off_axis_pos_m, nucfg);

  std::pair<std::vector<float>, std::vector<float>> ratios;
  for (size_t i = 0;
       i < OffAxisFluxUncertaintyHelper::Get().GetNFocussingParams(); i++) {
    ratios.first.push_back(
        float(OffAxisFluxUncertaintyHelper::Get().GetFluxFocussingRatio(
            i, flux_focussing_systbin, nucfg)));
  }
  for (size_t i = 0;
       i < OffAxisFluxUncertaintyHelper::Get().GetNHadProdPCAComponents();
       i++) {
    ratios.second.push_back(
        float(OffAxisFluxUncertaintyHelper::Get().GetFluxHadProdRatio(
            i, flux_hadprod_systbin, nucfg)));
  }
  return ratios;
}

void PrintFluxParameterNames() {

  std::cout << "=== Focusing Flux Params ===\n";
  for (size_t i = 0;
       i < OffAxisFluxUncertaintyHelper::Get().GetNFocussingParams(); i++) {
    std::cout << OffAxisFluxUncertaintyHelper::Get().GetFocussingParamName(i)
              << "\n";
  }

  std::cout << "=== Hadron-Production PCA Params ===\n";
  for (size_t i = 0;
       i < OffAxisFluxUncertaintyHelper::Get().GetNHadProdPCAComponents();
       i++) {
    std::cout << "Flux_HadProd_Param_" << i << "\n";
  }
}

} // namespace dune::beamoffaxis
