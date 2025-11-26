#include "Samples/SampleHandlerBeamOffAxis.h"

#include "Systematics/Flux/OffAxisFluxUncertaintyHelper.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"

SampleHandlerBeamOffAxis::SampleHandlerBeamOffAxis(
    std::string mc_version_, ParameterHandlerGeneric *ParHandler_,
    const std::shared_ptr<OscillationHandler> &Oscillator_)
    : SampleHandlerFD(mc_version_, ParHandler_, Oscillator_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;

  Initialise();
}

SampleHandlerBeamOffAxis::~SampleHandlerBeamOffAxis() {}

void SampleHandlerBeamOffAxis::Init() {
  if (CheckNodeExists(SampleManager->raw(), "DUNESampleBools", "isFHC")) {
    isFHC = SampleManager->raw()["DUNESampleBools"]["isFHC"].as<double>();
  } else {
    MACH3LOG_ERROR("Did not find DUNESampleBools:isFHC in {}, please add this",
                   SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

void SampleHandlerBeamOffAxis::SetupSplines() {

  ///@todo move all of the spline setup into core
  if (ParHandler->GetNumParamsFromSampleName(SampleName, kSpline) > 0) {
    MACH3LOG_INFO(
        "Found {} splines for this sample so I will create a spline object",
        ParHandler->GetNumParamsFromSampleName(SampleName, kSpline));
    SplineHandler = std::unique_ptr<BinnedSplineHandler>(
        new BinnedSplineHandlerDUNE(ParHandler, Modes.get()));
    InitialiseSplineObject();
  } else {
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or "
                  "evaluate splines",
                  ParHandler->GetNumParamsFromSampleName(SampleName, kSpline));
    SplineHandler = nullptr;
  }
}

void SampleHandlerBeamOffAxis::RegisterFunctionalParameters() {
  MACH3LOG_INFO("Registering functional parameters");
  // This function manually populates the map of functional parameters
  // Maps the name of the functional parameter to the pointer of the function

  RegisterIndividualFunctionalParameter(
      "TotalEScaleND", kTotalEScaleND,
      [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].shift.erec +=
            (*par) * dunemcSamples[iEvent].reco.EHad;
      });

  RegisterIndividualFunctionalParameter(
      "TotalEScaleND_sqrt", kTotalEScaleND_sqrt,
      [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].shift.erec +=
            (*par) * dunemcSamples[iEvent].reco.EHad *
            dunemcSamples[iEvent].syst.EHad_sqrt;
      });

  RegisterIndividualFunctionalParameter(
      "TotalEScaleND_invsqrt", kTotalEScaleND_invsqrt,
      [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].shift.erec +=
            (*par) * dunemcSamples[iEvent].syst.EHad_sqrt;
      });

  RegisterIndividualFunctionalParameter(
      "TotalEScaleND_mu", kTotalEScaleND_mu,
      [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].shift.erec +=
            (*par) * dunemcSamples[iEvent].reco.ELep;
      });

  RegisterIndividualFunctionalParameter(
      "TotalEScaleND_musqrt", kTotalEScaleND_musqrt,
      [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].shift.erec +=
            (*par) * dunemcSamples[iEvent].reco.ELep *
            dunemcSamples[iEvent].syst.ELep_sqrt;
      });

  RegisterIndividualFunctionalParameter(
      "TotalEScaleND_muinvsqrt", kTotalEScaleND_muinvsqrt,
      [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].shift.erec +=
            (*par) * dunemcSamples[iEvent].syst.ELep_sqrt;
      });

  RegisterIndividualFunctionalParameter(
      "TotalEScaleND_had", kTotalEScaleND_had,
      [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].shift.erec +=
            (*par) * dunemcSamples[iEvent].reco.sum_ehad;
      });

  RegisterIndividualFunctionalParameter(
      "TotalEScaleND_hadsqrt", kTotalEScaleND_hadsqrt,
      [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].shift.erec +=
            (*par) * dunemcSamples[iEvent].reco.sum_ehad *
            dunemcSamples[iEvent].syst.sum_ehad_sqrt;
      });

  RegisterIndividualFunctionalParameter(
      "TotalEScaleND_hadinvsqrt", kTotalEScaleND_hadinvsqrt,
      [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].shift.erec +=
            (*par) * dunemcSamples[iEvent].syst.sum_ehad_sqrt;
      });

  RegisterIndividualFunctionalParameter(
      "TotalEScaleND_EM", kTotalEScaleND_EM,
      [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].shift.erec +=
            (*par) * (dunemcSamples[iEvent].reco.ePi0 -
                      dunemcSamples[iEvent].truth.ePi0);
      });

  RegisterIndividualFunctionalParameter(
      "TotalEScaleND_EM_invsqrt", kTotalEScaleND_EMinvsqrt,
      [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].shift.erec +=
            (*par) *
            (dunemcSamples[iEvent].reco.ePi0 -
             dunemcSamples[iEvent].truth.ePi0) *
            sqrt((dunemcSamples[iEvent].reco.ePi0 -
                  dunemcSamples[iEvent].truth.ePi0));
      });

  RegisterIndividualFunctionalParameter(
      "TotalEScaleND_EM_sqrt", kTotalEScaleND_EMsqrt,
      [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].shift.erec +=
            (*par) * sqrt((dunemcSamples[iEvent].reco.ePi0 -
                           dunemcSamples[iEvent].truth.ePi0));
      });

  RegisterIndividualFunctionalParameter(
      "MuonRes_ND", kMuonRes_ND, [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].shift.erec +=
            (*par) * (dunemcSamples[iEvent].truth.LepE -
                      dunemcSamples[iEvent].reco.ELep);
      });

  RegisterIndividualFunctionalParameter(
      "NRes_ND", kNRes_ND, [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].shift.erec +=
            (*par) *
            (dunemcSamples[iEvent].truth.eN - dunemcSamples[iEvent].reco.eN);
      });

  RegisterIndividualFunctionalParameter(
      "HadRes_ND", kHadRes_ND, [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].shift.erec +=
            (*par) * ((dunemcSamples[iEvent].truth.ePim -
                       dunemcSamples[iEvent].reco.ePim) +
                      ((dunemcSamples[iEvent].truth.ePip -
                        dunemcSamples[iEvent].reco.ePip)) +
                      ((dunemcSamples[iEvent].truth.eP -
                        dunemcSamples[iEvent].reco.eP)));
      });

  for (size_t par_it = 0;
       par_it < OffAxisFluxUncertaintyHelper::Get().GetNFocussingParams();
       par_it++) {
    RegisterIndividualFunctionalParameter(
        OffAxisFluxUncertaintyHelper::Get().GetFocussingParamName(par_it),
        int(kNFuncPars + par_it),
        [this, par_it](const double *par, std::size_t iEvent) {
          dunemcSamples[iEvent].flux_w *=
              (1 + (*par) *
                       dunemcSamples[iEvent].syst.flux_focussing_ratio[par_it]);
        });
  }

  for (size_t par_it = 0;
       par_it < OffAxisFluxUncertaintyHelper::Get().GetNHadProdPCAComponents();
       par_it++) {
    RegisterIndividualFunctionalParameter(
        "Flux_HadProd_Param_" + std::to_string(par_it),
        int(kNFuncPars +
            OffAxisFluxUncertaintyHelper::Get().GetNFocussingParams() + par_it),
        [this, par_it](const double *par, std::size_t iEvent) {
          dunemcSamples[iEvent].flux_w *=
              (1 +
               (*par) * dunemcSamples[iEvent].syst.flux_hadprod_ratio[par_it]);
        });
  }

  MACH3LOG_INFO("Finished registering functional parameters");
}

// HH: Reset the shifted values to the original values
void SampleHandlerBeamOffAxis::resetShifts(int iEvent) {
  dunemcSamples[iEvent].shift.erec = dunemcSamples[iEvent].rw_erec;
  dunemcSamples[iEvent].flux_w = 1.0;
}
// =================================

void SampleHandlerBeamOffAxis::SetupWeightPointers() {
  for (size_t i = 0; i < dunemcSamples.size(); ++i) {
    MCSamples[i].total_weight_pointers.push_back(&(dunemcSamples[i].pot_s));
    MCSamples[i].total_weight_pointers.push_back(&(dunemcSamples[i].norm_s));
    MCSamples[i].total_weight_pointers.push_back(MCSamples[i].osc_w_pointer);
    MCSamples[i].total_weight_pointers.push_back(
        &(dunemcSamples[i].rw_berpaacvwgt));
    MCSamples[i].total_weight_pointers.push_back(&(dunemcSamples[i].flux_w));
    MCSamples[i].total_weight_pointers.push_back(&(MCSamples[i].xsec_w));
  }
}

double ERecQE(int nupdg, bool isCC, double el, double leptheta) {
  constexpr double V = 0;        // 0 binding energy for now
  constexpr double mn = 939.565; // neutron mass
  constexpr double mp = 938.272; // proton mass
  double mN_eff = mn - V;
  double mN_oth = mp;

  if (nupdg < 0) { // if anti-neutrino, swap target/out masses
    mN_eff = mp - V;
    mN_oth = mn;
  }

  // this is funky, but don't be scared, it defines an annonymous function
  // in place that grabs the lepton mass in MeV when given the neutrino PDG
  // and whether the interaction was CC or NC and then immediately calls it.
  // It's basically a generalisation of the ternary operator.
  double ml = [=]() {
    switch (std::abs(nupdg)) {
    case 12:
      return isCC ? 0.511 : 0;
    case 14:
      return isCC ? 105.66 : 0;
    case 16:
      return isCC ? 1777.0 : 0;
    default:
      std::cerr << "Warning: Unexpected PDG code " << nupdg
                << " passed to ml lambda.\n";
      assert(false && "Unexpected neutrino PDG code in ml lambda");
      return 0.0;
    }
  }();

  double pl = std::sqrt(el * el - ml * ml); // momentum of lepton

  return (2 * mN_eff * el - ml * ml + mN_oth * mN_oth - mN_eff * mN_eff) /
         (2 * (mN_eff - el + pl * std::cos(leptheta)));
}

int SampleHandlerBeamOffAxis::SetupExperimentMC() {

  MACH3LOG_INFO(
      "-------------------------------------------------------------------");
  TChain CAFChain("caf");
  TChain CAFMetaChain("meta");
  for (size_t iSample = 0; iSample < mc_files.size(); iSample++) {
    MACH3LOG_INFO("Adding file to TChains: {}", mc_files[iSample]);

    // HH: Check whether the file exists, see
    // https://root.cern/doc/master/classTChain.html#a78a896924ac6c7d3691b7e013bcbfb1c
    if (!CAFChain.Add(mc_files[iSample].c_str(), -1)) {
      MACH3LOG_ERROR("Could not add file {} to TChain, please check the file "
                     "exists and is readable",
                     mc_files[iSample]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    CAFMetaChain.Add(mc_files[iSample].c_str(), -1);
  }

  MACH3LOG_INFO("Number of entries in CAF TChain: {}", CAFChain.GetEntries());

  TTreeReader chrdr(&CAFChain);

  // Reco Variables
  TTreeReaderValue<double> erec(chrdr, "erec");
  TTreeReaderValue<double> Elep_reco(chrdr, "Elep_reco");

  TTreeReaderValue<int> reco_numu(chrdr, "reco_numu");
  TTreeReaderValue<int> muon_contained(chrdr, "muon_contained");
  TTreeReaderValue<int> muon_tracker(chrdr, "muon_tracker");

  TTreeReaderValue<double> eRecoP(chrdr, "eRecoP");
  TTreeReaderValue<double> eRecoPip(chrdr, "eRecoPip");
  TTreeReaderValue<double> eRecoPim(chrdr, "eRecoPim");
  TTreeReaderValue<double> eRecoPi0(chrdr, "eRecoPi0");
  TTreeReaderValue<double> eRecoN(chrdr, "eRecoN");

  TTreeReaderValue<double> vtx_x(chrdr, "vtx_x");
  TTreeReaderValue<double> vtx_y(chrdr, "vtx_y");
  TTreeReaderValue<double> vtx_z(chrdr, "vtx_z");
  TTreeReaderValue<double> det_x(chrdr, "det_x");

  // Truth Variables
  TTreeReaderValue<int> mode(chrdr, "mode");
  TTreeReaderValue<double> ev(chrdr, "ev");

  TTreeReaderValue<double> LepE(chrdr, "LepE");
  TTreeReaderValue<double> LepNuAngle(chrdr, "LepNuAngle");
  TTreeReaderValue<double> LepMomX(chrdr, "LepMomX");
  TTreeReaderValue<double> LepMomY(chrdr, "LepMomY");
  TTreeReaderValue<double> LepMomZ(chrdr, "LepMomZ");
  TTreeReaderValue<double> NuMomX(chrdr, "NuMomX");
  TTreeReaderValue<double> NuMomY(chrdr, "NuMomY");
  TTreeReaderValue<double> NuMomZ(chrdr, "NuMomZ");
  TTreeReaderValue<double> eP(chrdr, "eP");
  TTreeReaderValue<double> ePip(chrdr, "ePip");
  TTreeReaderValue<double> ePim(chrdr, "ePim");
  TTreeReaderValue<double> ePi0(chrdr, "ePi0");
  TTreeReaderValue<double> eN(chrdr, "eN");
  TTreeReaderValue<double> eOther(chrdr, "eOther");

  TTreeReaderValue<int> nipi0(chrdr, "nipi0");

  TTreeReaderValue<int> isCC(chrdr, "isCC");
  TTreeReaderValue<int> nuPDG(chrdr, "nuPDG");
  TTreeReaderValue<int> nuPDGunosc(chrdr, "nuPDGunosc");

  dunemcSamples.resize(CAFChain.GetEntries());

  // HH: A map to keep track of negative energies
  std::unordered_map<std::string, int> negative_counts;
  // Initialize the negative counts for each energy variable
  negative_counts["rw_erec_had"] = 0;
  negative_counts["rw_erec_lep"] = 0;
  negative_counts["rw_eRecoN"] = 0;
  negative_counts["rw_eRecoPi0"] = 0;
  negative_counts["rw_sum_ehad"] = 0;

  size_t iEvent = 0;
  while (chrdr.Next()) {

    // fill dunemc_base
    dunemcSamples[iEvent].Target = kTarget_Ar;

    dunemcSamples[iEvent].nupdg = *nuPDG;
    dunemcSamples[iEvent].nupdgUnosc = *nuPDGunosc;

    dunemcSamples[iEvent].rw_isCC = *isCC;

    dunemcSamples[iEvent].OscChannelIndex = static_cast<double>(
        GetOscChannel(OscChannels, dunemcSamples[iEvent].nupdgUnosc,
                      dunemcSamples[iEvent].nupdg));

    dunemcSamples[iEvent].rw_erec = *erec;
    dunemcSamples[iEvent].rw_etru = *ev;
    dunemcSamples[iEvent].flux_w = 1.0;

    int M3Mode = Modes->GetModeFromGenerator(std::abs(*mode));
    if (!*isCC) {
      M3Mode += 14; // Account for no ability to distinguish CC/NC
    }
    if (M3Mode > 15) {
      M3Mode -= 1; // Account for no NCSingleKaon
    }
    dunemcSamples[iEvent].mode = M3Mode;

    // fill dunemc_beamoffaxis
    dunemcSamples[iEvent].truth.LepE = *LepE;
    dunemcSamples[iEvent].truth.LepNuAngle = *LepNuAngle;
    dunemcSamples[iEvent].truth.eP = *eP;
    dunemcSamples[iEvent].truth.ePip = *ePip;
    dunemcSamples[iEvent].truth.ePim = *ePim;
    dunemcSamples[iEvent].truth.ePi0 = *ePi0;
    dunemcSamples[iEvent].truth.eN = *eN;
    dunemcSamples[iEvent].truth.eHad_av =
        *eP + *ePip + *ePim + *ePi0 + *eOther + (*nipi0 * 0.1349);
    dunemcSamples[iEvent].truth.enu_bias = dunemcSamples[iEvent].truth.LepE +
                                           dunemcSamples[iEvent].truth.eHad_av -
                                           dunemcSamples[iEvent].rw_etru;

    dunemcSamples[iEvent].truth.det_x = *det_x;
    dunemcSamples[iEvent].truth.vtx_x = *vtx_x;
    dunemcSamples[iEvent].truth.vtx_y = *vtx_y;
    dunemcSamples[iEvent].truth.vtx_z = *vtx_z;

    dunemcSamples[iEvent].reco.numu = *reco_numu;
    dunemcSamples[iEvent].reco.muon_contained = *muon_contained;
    dunemcSamples[iEvent].reco.muon_tracker = *muon_tracker;

    dunemcSamples[iEvent].reco.ELep = *Elep_reco;
    dunemcSamples[iEvent].reco.EHad =
        (dunemcSamples[iEvent].rw_erec - dunemcSamples[iEvent].reco.ELep);
    dunemcSamples[iEvent].reco.sum_ehad = *eRecoP + *eRecoPip + *eRecoPim;
    dunemcSamples[iEvent].reco.y =
        (dunemcSamples[iEvent].reco.EHad / dunemcSamples[iEvent].rw_erec);

    dunemcSamples[iEvent].reco.eP = *eRecoP;
    dunemcSamples[iEvent].reco.ePip = *eRecoPip;
    dunemcSamples[iEvent].reco.ePim = *eRecoPim;
    dunemcSamples[iEvent].reco.ePi0 = *eRecoPi0;
    dunemcSamples[iEvent].reco.eN = *eRecoN;

    dunemcSamples[iEvent].shift.erec = dunemcSamples[iEvent].rw_erec;

    dunemcSamples[iEvent].syst.EHad_sqrt =
        std::sqrt(dunemcSamples[iEvent].reco.EHad);
    dunemcSamples[iEvent].syst.ELep_sqrt =
        std::sqrt(dunemcSamples[iEvent].reco.ELep);
    dunemcSamples[iEvent].syst.sum_ehad_sqrt =
        std::sqrt(dunemcSamples[iEvent].reco.sum_ehad);

    dunemcSamples[iEvent].truth.ERec_QE =
        ERecQE(dunemcSamples[iEvent].nupdg, dunemcSamples[iEvent].rw_isCC,
               dunemcSamples[iEvent].truth.LepE,
               dunemcSamples[iEvent].truth.LepNuAngle);

    dunemcSamples[iEvent].truth.off_axis_pos_m = *det_x * 100 + *vtx_x;

    auto nucfg = OffAxisFluxUncertaintyHelper::Get().GetNuConfig(
        dunemcSamples[iEvent].nupdg, true, isFHC, false);
    auto flux_focussing_systbin =
        OffAxisFluxUncertaintyHelper::Get().GetFocussingBin(
            dunemcSamples[iEvent].nupdg, dunemcSamples[iEvent].rw_etru,
            dunemcSamples[iEvent].truth.off_axis_pos_m, nucfg);

    auto flux_hadprod_systbin =
        OffAxisFluxUncertaintyHelper::Get().GetFocussingBin(
            dunemcSamples[iEvent].nupdg, dunemcSamples[iEvent].rw_etru,
            dunemcSamples[iEvent].truth.off_axis_pos_m, nucfg);

    for (size_t i = 0;
         i < OffAxisFluxUncertaintyHelper::Get().GetNFocussingParams(); i++) {
      dunemcSamples[iEvent].syst.flux_focussing_ratio.push_back(
          OffAxisFluxUncertaintyHelper::Get().GetFluxFocussingRatio(
              i, flux_focussing_systbin, nucfg));
    }
    for (size_t i = 0;
         i < OffAxisFluxUncertaintyHelper::Get().GetNHadProdPCAComponents();
         i++) {
      dunemcSamples[iEvent].syst.flux_hadprod_ratio.push_back(
          OffAxisFluxUncertaintyHelper::Get().GetFluxHadProdRatio(
              i, flux_hadprod_systbin, nucfg));
    }

    // HH: Give a warning if any negative energies were found
    for (const auto &[key, count] : negative_counts) {
      if (count > 0) {
        MACH3LOG_WARN("Found {} negative values for {} in sample {}", count,
                      key, GetSampleName());
      }
    }
  }

  return int(dunemcSamples.size());
}

const double *
SampleHandlerBeamOffAxis::GetPointerToKinematicParameter(KinematicTypes KinPar,
                                                         int iEvent) {
  switch (KinPar) {
  case kTrueNeutrinoEnergy:
    return &dunemcSamples[iEvent].rw_etru;
  case kRecoNeutrinoEnergy:
    return &dunemcSamples[iEvent].shift.erec;
  case kTrueXPos:
    return &dunemcSamples[iEvent].truth.vtx_x;
  case kTrueYPos:
    return &dunemcSamples[iEvent].truth.vtx_y;
  case kTrueZPos:
    return &dunemcSamples[iEvent].truth.vtx_z;
  case kM3Mode:
    return &dunemcSamples[iEvent].mode;
  case kIsFHC:
    return &isFHC;
  case kELepRec:
    return &dunemcSamples[iEvent].reco.ELep;
  case kEnubias:
    return &dunemcSamples[iEvent].truth.enu_bias;
  default:
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type...");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
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
  for (int iEvent = 0; iEvent < int(GetNEvents()); ++iEvent) {
    MCSamples[iEvent].rw_etru = &(dunemcSamples[iEvent].rw_etru);
    MCSamples[iEvent].mode = &(dunemcSamples[iEvent].mode);
    MCSamples[iEvent].Target = &(dunemcSamples[iEvent].Target);
    MCSamples[iEvent].isNC = !(dunemcSamples[iEvent].rw_isCC);
    MCSamples[iEvent].nupdg = &(dunemcSamples[iEvent].nupdg);
    MCSamples[iEvent].nupdgUnosc = &(dunemcSamples[iEvent].nupdgUnosc);
  }
}
