#include "Samples/SampleHandlerBeamOffAxis.h"

#include "Samples/BeamOffAxis/ReadEvents.h"
#include "Samples/BeamOffAxis/Systematics.h"

#include <iostream>

#include <fstream>

namespace dune::beamoffaxis {

  SampleHandlerBeamOffAxis::SampleHandlerBeamOffAxis(
    std::string mc_version_, ParameterHandlerGeneric *ParHandler_,
    const std::shared_ptr<OscillationHandler> &Oscillator_)
    : SampleHandlerFD(mc_version_, ParHandler_, Oscillator_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  Initialise();

  unweighted_mc = GetUnweightedMCRate();
}

void SampleHandlerBeamOffAxis::Init() {
  subsample_analysispot.resize(GetNsamples());
  subsample_is_numode.resize(GetNsamples());

  for (int iSubSample = 0; iSubSample < GetNsamples(); iSubSample++) {
    auto const &sample_conf = SampleManager->raw()[GetSampleTitle(iSubSample)];
    subsample_analysispot[iSubSample] =
        Get<double>(sample_conf["POT"], __FILE__, __LINE__);
    subsample_is_numode[iSubSample] =
        Get<bool>(sample_conf["is_numode"], __FILE__, __LINE__);
  }

  if (SampleManager->raw()["InputFiles"]["CovarianceMatrix"]) {
    auto cvmx_details = Get<std::vector<std::string>>(
        SampleManager->raw()["InputFiles"]["CovarianceMatrix"], __FILE__,
        __LINE__);

    TFile cvmx_file(cvmx_details[0].c_str(), "READ");
    auto *rcvmx = cvmx_file.Get<TMatrixD>(cvmx_details[1].c_str());

    cvmx = Eigen::Map<Eigen::MatrixXd>(rcvmx->GetMatrixArray(),
                                       rcvmx->GetNrows(), rcvmx->GetNcols());

    MACH3LOG_INFO("Using ND Covariance Matrix({},{}):", cvmx.rows(),
                  cvmx.cols());
    std::stringstream ss;
    ss << cvmx;
    MACH3LOG_INFO("\n{}", ss.str());
  }
}

void SampleHandlerBeamOffAxis::SetupSplines() {
  if (!ParHandler) {
    return;
  }

  ///@todo move all of the spline setup into core

  int num_splines = 0;
  for (int iSubSample = 0; iSubSample < int(SampleDetails.size());
       iSubSample++) {
    num_splines += ParHandler->GetNumParamsFromSampleName(
        GetSampleTitle(iSubSample), kSpline);
  }

  if (num_splines > 0) {
    MACH3LOG_INFO(
        "Found {} splines for this sample so I will create a spline object",
        num_splines);
    SplineHandler = std::unique_ptr<BinnedSplineHandler>(
        new BinnedSplineHandlerDUNE(ParHandler, Modes.get()));
    InitialiseSplineObject();
  } else {
    MACH3LOG_INFO(
        "Found {} splines for this sample so I will not load or "
        "evaluate splines",
        ParHandler->GetNumParamsFromSampleName(SampleHandlerName, kSpline));
    SplineHandler = nullptr;
  }
}

void SampleHandlerBeamOffAxis::RegisterFunctionalParameters() {
  if (!ParHandler) {
    return;
  }

  RegisterIndividualFunctionalParameter(DUNEMCEvents,
                                        {
                                            "ContainedMuonEnergyScale",
                                            "ContainedMuonSqrtEnergyScale",
                                            "ContainedMuonInvSqrtEnergyScale",
                                            "TrackedMuonEnergyScale",
                                            "TrackedMuonSqrtEnergyScale",
                                            "TrackedMuonInvSqrtEnergyScale",
                                            "EMEnergyScale",
                                            "EMSqrtEnergyScale",
                                            "EMInvSqrtEnergyScale",
                                            "ChgHadEnergyScale",
                                            "ChgHadSqrtEnergyScale",
                                            "ChgHadInvSqrtEnergyScale",
                                            "NeutronEnergyScale",
                                            "NeutronSqrtEnergyScale",
                                            "NeutronInvSqrtEnergyScale",
                                            "TotalEnergyScale",
                                            "TotalSqrtEnergyScale",
                                            "TotalInvSqrtEnergyScale",
                                        },
                                        EnergyScales);

  RegisterIndividualFunctionalParameter(
      DUNEMCEvents,
      {"MuonEnergyResolution", "EMEnergyResolution", "ChgHadEnergyResolution",
       "NeutronEnergyResolution"},
      ParticleEnergyResolutions);

  if (ParHandler->GetNumParFromGroup("Flux")) {
    RegisterIndividualFunctionalParameter(
        DUNEMCEvents, GetFluxFocussingParamNames(), UpdateFluxFocussingWeight);

    RegisterIndividualFunctionalParameter(
        DUNEMCEvents, GetFluxHadProdParamNames(), UpdateFluxHadProdWeight);
  }

  RegisterIndividualFunctionalParameter(DUNEMCEvents, "MissingProtonFD",
                                        MissingProtonFD);
}

Eigen::VectorXd SampleHandlerBeamOffAxis::GetUnweightedMCRate() {
  Eigen::VectorXd mc = Eigen::VectorXd::Zero(Binning->GetNBins());

  for (int sample_i = 0; sample_i < GetNsamples(); ++sample_i) {
    const int ndim = GetNDim(sample_i);

    if (ndim == 1) {
      auto mc_hist = Get1DVarHist(sample_i, GetXBinVarName(sample_i), StoredSelection[sample_i], 1);
      if (!mc_hist) throw MaCh3Exception(__FILE__, __LINE__);

      for (int i = 0; i < mc_hist->GetNbinsX(); ++i) {
        const int global_bin = Binning->GetGlobalBinSafe(sample_i, {i});
        mc(global_bin) = mc_hist->GetBinContent(i + 1);
      }
    } else if (ndim == 2) {
      auto mc_hist = Get2DVarHist(sample_i, GetXBinVarName(sample_i), GetYBinVarName(sample_i),
                                        StoredSelection[sample_i], 1);
      if (!mc_hist) throw MaCh3Exception(__FILE__, __LINE__);

      for (int j = 0; j < mc_hist->GetNbinsY(); ++j) {
        for (int i = 0; i < mc_hist->GetNbinsX(); ++i) {
          const int global_bin = Binning->GetGlobalBinSafe(sample_i, {i, j});
          mc(global_bin) = mc_hist->GetBinContent(i + 1, j + 1);
        }
      }
    } else {
      MACH3LOG_ERROR("GetUnweightedMCRate: sample {} has {} dimensions, only 1D or 2D supported",
                     sample_i, ndim);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  return mc;
}

void SampleHandlerBeamOffAxis::ResetShifts(int iEvent) {
  auto &ev = DUNEMCEvents[iEvent];
  ev.varied_reco = ev.reco;

  // resolution variables
  ev.varied_res.enu = ev.reco.enu - ev.truth.nu.e;
  ev.varied_res.e_lep = ev.reco.e_lep - ev.truth.lep.e;
  ev.varied_res.e_had =
      (ev.reco.enu - ev.reco.e_lep) - (ev.truth.nu.e - ev.truth.lep.e);

  ev.varied_res.e_EM = (ev.reco.e_pi0 - ev.truth.had.e_pi0);
  if (std::abs(ev.truth.lep.pdg) == 11) {
    ev.varied_res.e_EM += ev.varied_res.e_lep;
  }
  ev.varied_res.e_ChgHad = (ev.reco.e_proton - ev.truth.had.e_proton) +
                           (ev.reco.e_piplus - ev.truth.had.e_piplus) +
                           (ev.reco.e_piminus - ev.truth.had.e_piminus);
  ev.varied_res.e_neutron = ev.reco.e_neutron - ev.truth.had.e_neutron;

  ev.varied_truth.enurec_hadavailable_missed =
      ev.truth.kine.enurec_hadavailable_missed;

  // flux weights
  ev.syst.flux.total_weight = 1.0;
}

void SampleHandlerBeamOffAxis::FinaliseShifts(int iEvent) {
  CalculateVariedCompositeQuantities(DUNEMCEvents[iEvent]);
}

void SampleHandlerBeamOffAxis::AddAdditionalWeightPointers() {
  for (size_t i = 0; i < DUNEMCEvents.size(); ++i) {
    MCEvents[i].total_weight_pointers.push_back(&(DUNEMCEvents[i].weights.pot));
    MCEvents[i].total_weight_pointers.push_back(
        &(DUNEMCEvents[i].syst.flux.total_weight));
  }
}

int SampleHandlerBeamOffAxis::SetupExperimentMC() {

  MACH3LOG_INFO(
      "-------------------------------------------------------------------");

  bool do_flux_systematics =
      ParHandler && ParHandler->GetNumParFromGroup("Flux");

  for (int iSubSample = 0; iSubSample < int(SampleDetails.size());
       iSubSample++) {
    MACH3LOG_INFO("-- subsample[{}]: {} ", iSubSample,
                  SampleDetails[iSubSample].SampleTitle);

    TChain MetaChain("meta");
    TChain CAFChain("cafTree");
    for (const std::string &filename : SampleDetails[iSubSample].mc_files) {
      if (filename.empty()) {
        MACH3LOG_INFO("-- -- Skipping empty filename entry");
        continue;
    }
      MACH3LOG_INFO("-- -- Adding file to TChain: {}", filename);
      if (!CAFChain.Add(filename.c_str(), -1)) {
        MACH3LOG_ERROR("Could not add file {} to TChain, please check the file "
                       "exists and is readable",
                       filename);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      MetaChain.Add(filename.c_str(), -1);
    }

    double subsample_cafpot = GetPOT(MetaChain);
    auto sample_evs = ReadEvents(CAFChain);

    // fix up any analysis specific information
    int n_contained = 0, n_tracked = 0, n_muon_neither = 0, n_not_muon = 0;
    for (auto &ev : sample_evs) {

      ev.subsample = iSubSample;
      ev.is_numode = subsample_is_numode[iSubSample];

      ev.weights.pot = subsample_analysispot[iSubSample] / subsample_cafpot;
     // std::cout<< "pot scaling = " << ev.weights.pot << std::endl;

      ev.truth.mach3_mode =
          Modes->GetModeFromGenerator(std::abs(ev.truth.generator_mode));
      if (!ev.truth.is_cc) {
        // Account for no ability to distinguish CC/NC
        ev.truth.mach3_mode += 14;
      }
      if (ev.truth.mach3_mode > 15) {
        // Account for no NCSingleKaon
        ev.truth.mach3_mode -= 1;
      }

      ev.syst.flux.total_weight = 1;
      if (do_flux_systematics) {
        std::tie(ev.syst.flux.focussing_ratio, ev.syst.flux.hadprod_ratio) =
            GetFluxVariationRatios(ev.truth.nu.pdg, ev.truth.nu.e,
                                   ev.truth.vtx.off_axis_pos_m, true);
      }
      // Add this block at the end of the loop body:
    if (std::abs(ev.truth.lep.pdg) == 13) {
        if (ev.reco.muonlike_contained)       n_contained++;
        else if (ev.reco.muonlike_tracker)    n_tracked++;
        else                                   n_muon_neither++;
    } else {
        n_not_muon++;
    }
    }
    // After the loop, print the results:
MACH3LOG_INFO("Subsample [{}] muon classification:", iSubSample);
MACH3LOG_INFO("  true muon + contained:  {}", n_contained);
MACH3LOG_INFO("  true muon + tracked:    {}", n_tracked);
MACH3LOG_INFO("  true muon + neither:    {}", n_muon_neither);
MACH3LOG_INFO("  not true muon (pdg!=13): {}", n_not_muon);

    DUNEMCEvents.reserve(DUNEMCEvents.size() + CAFChain.GetEntries());
    std::copy(sample_evs.begin(), sample_evs.end(),
              std::back_inserter(DUNEMCEvents));
  }

  return int(DUNEMCEvents.size());
}

const double *
SampleHandlerBeamOffAxis::GetPointerToKinematicParameter(KinematicTypes KinPar,
                                                         int iEvent) {
  return ResolveKinematicEventMember(KinPar, DUNEMCEvents[iEvent]);
}

const double *SampleHandlerBeamOffAxis::GetPointerToKinematicParameter(
    std::string KinematicParameter, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(
      ReturnKinematicParameterFromString(KinematicParameter));
  return GetPointerToKinematicParameter(KinPar, iEvent);
}

const double *SampleHandlerBeamOffAxis::GetPointerToKinematicParameter(
    double KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return GetPointerToKinematicParameter(KinPar, iEvent);
}

double SampleHandlerBeamOffAxis::ReturnKinematicParameter(int KinematicVariable,
                                                          int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iEvent);
}

double SampleHandlerBeamOffAxis::ReturnKinematicParameter(
    std::string KinematicParameter, int iEvent) {
  return *GetPointerToKinematicParameter(KinematicParameter, iEvent);
}

void SampleHandlerBeamOffAxis::SetupFDMC() {
  size_t iEvent = 0;
  for (auto const &ev : DUNEMCEvents) {
    MCEvents[iEvent].Target = ev.truth.target_a;
    MCEvents[iEvent].mode = ev.truth.mach3_mode;
    MCEvents[iEvent].isNC = !ev.truth.is_cc;

    MCEvents[iEvent].enu_true = ev.truth.nu.e;
    MCEvents[iEvent].nupdg = ev.truth.nu.pdg;
    MCEvents[iEvent].nupdgUnosc = ev.truth.nu.pdg_unosc;

    MCEvents[iEvent].NominalSample = ev.subsample;

    iEvent++;
  }
}

void SampleHandlerBeamOffAxis::BuildRegularisationMatrix(
    ParameterHandlerRegularised *RegParHandler, double lambda) {

  MACH3LOG_INFO("Building Regularisation matrix "
                "with penalty term = {}", lambda);

  struct TemplateParameterBins { //struct for each template parameter (Etrue,enubias)
    int index;
    double enu_lo, enu_hi;
    double enubias_lo, enubias_hi;
  };

  std::vector<TemplateParameterBins> template_parameters; //vector of all the template parameters

  auto normPars = RegParHandler->GetNormParsFromSampleName(SampleHandlerName); //find them all
  for (auto const &np : normPars) {
    if (!RegParHandler->IsParFromGroup(np.index, "Xsec")) continue;

    double enu_lo = -1e10, enu_hi = 1e10;
    double enubias_lo = -1e10, enubias_hi = 1e10;

    for (size_t iVar = 0; iVar < np.KinematicVarStr.size(); ++iVar) {
      if (np.KinematicVarStr[iVar] == "TrueNeutrinoEnergy") {
        enu_lo = np.Selection[iVar][0][0];
        enu_hi = np.Selection[iVar][0][1];
      } else if (np.KinematicVarStr[iVar] == "Enubias") {
        enubias_lo = np.Selection[iVar][0][0];
        enubias_hi = np.Selection[iVar][0][1];
      }
    }

    template_parameters.push_back({np.index, enu_lo, enu_hi, enubias_lo, enubias_hi});
  }

  if (template_parameters.empty()) {
    MACH3LOG_WARN("BuildRegularisationMatrix: no Xsec Norm parameters found, "
                  "regularisation will be zero");
    return;
  }

  MACH3LOG_INFO("  Found {} Xsec Normalisation parameters", template_parameters.size());

  auto round6 = [](double x) { return std::round(x * 1e6) / 1e6; };

  std::set<double> enu_edges_set, enubias_edges_set;
  for (auto const &p : template_parameters) {
    enu_edges_set.insert(round6(p.enu_lo));
    enubias_edges_set.insert(round6(p.enubias_lo));
  }

  std::vector<double> enu_edges(enu_edges_set.begin(), enu_edges_set.end());
  std::vector<double> enubias_edges(enubias_edges_set.begin(), enubias_edges_set.end());

  int nEnu     = int(enu_edges.size());
  int nEnubias = int(enubias_edges.size());

  MACH3LOG_INFO("  Found {} True Neutrino Energy bins and {} Enubias bins", nEnu, nEnubias);

  auto edge_index = [&](std::vector<double> const &edges, double val) {
    auto it = std::lower_bound(edges.begin(), edges.end(), round6(val));
    if (it == edges.end()) {
      MACH3LOG_ERROR("Could not find bin edge {} in sorted edges", val);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    return int(std::distance(edges.begin(), it));
  };

  struct IndexOfTemplateParam { int index, iEnuTrue, iEnuTruebias; };
  std::vector<IndexOfTemplateParam> indexed;
  indexed.reserve(template_parameters.size());

  for (auto const &p : template_parameters) {
    indexed.push_back({p.index,
                       edge_index(enu_edges,     p.enu_lo),
                       edge_index(enubias_edges, p.enubias_lo)});
  }

  // Build first-difference penalty matrix R = D^T D
  int nTotal = RegParHandler->GetNumParams();
  Eigen::MatrixXd R = Eigen::MatrixXd::Zero(nTotal, nTotal);

  std::map<std::pair<int,int>, int> bin_to_global;
  for (auto const &ip : indexed) {
    bin_to_global[{ip.iEnuTrue, ip.iEnuTruebias}] = ip.index;
  }

  int n_pairs = 0;
  for (int iEnuTruebias = 0; iEnuTruebias < nEnubias; ++iEnuTruebias) {
    for (int iEnuTrue = 0; iEnuTrue < nEnu - 1; ++iEnuTrue) {
      auto it_lo = bin_to_global.find({iEnuTrue,     iEnuTruebias});
      auto it_hi = bin_to_global.find({iEnuTrue + 1, iEnuTruebias});

      if (it_lo == bin_to_global.end() || it_hi == bin_to_global.end()) continue;

      int gi = it_lo->second;
      int gj = it_hi->second;

      R(gi, gi) += 1.0;
      R(gj, gj) += 1.0;
      R(gi, gj) -= 1.0;
      R(gj, gi) -= 1.0;
      ++n_pairs;
    }
  }


RegParHandler->SetPenalty(
    [R = std::move(R), lambda, nTotal](std::vector<double> const &propVal) {
      std::size_t counter = 0;
        Eigen::VectorXd shift(nTotal);
        for (int i = 0; i < nTotal; ++i) {
            shift(i) = propVal[i] - 1.0;
        }
        double penalty = 0.5 * lambda * double(shift.transpose() * R * shift);
        ++counter;
        if (counter % 100000 == 0) {
            MACH3LOG_INFO(
                "Regularisation penalty at call {} = {:.4f}",
                counter, penalty);
        }
        //MACH3LOG_INFO("Regularisation penalty = {:.4f}", pen);
        return penalty;
    });
}

} // namespace dune::beamoffaxis
