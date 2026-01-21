#include "Samples/SampleHandlerBeamOffAxis.h"

#include "Systematics/Flux/OffAxisFluxUncertaintyHelper.h"

#include "Utils/FakeDataGenerators/FDTDRFDS/WeightToNuWro.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <iostream>

#include <fstream>

SampleHandlerBeamOffAxis::SampleHandlerBeamOffAxis(
    std::string mc_version_, ParameterHandlerGeneric *ParHandler_,
    const std::shared_ptr<OscillationHandler> &Oscillator_)
    : SampleHandlerFD(mc_version_, ParHandler_, Oscillator_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;

  Initialise();
}

double SampleHandlerBeamOffAxis::CalculatePOT() {
  TChain calc_pot_chain("meta"); // Use correct tree name

  std::string pot_branch = "pot"; // Use the correct branch name

  for (size_t i = 0; i < mc_files.size(); ++i) {
    calc_pot_chain.AddFile(mc_files[i].c_str());
  }

  std::cout << "Added " << calc_pot_chain.GetListOfFiles()->GetEntries()
            << " files to POT chain" << std::endl;

  std::cout << "calc_pot_chain.GetEntries() = " << calc_pot_chain.GetEntries()
            << std::endl;

  // Check if the branch exists before proceeding
  if (!calc_pot_chain.GetBranch(pot_branch.c_str())) {
    std::cerr << "Error: Branch " << pot_branch << " not found in the tree!"
              << std::endl;
    return 0.0;
  }

  double pot_value = 0.0;
  calc_pot_chain.SetBranchAddress(pot_branch.c_str(), &pot_value);

  double sum_pot = 0.0;
  Long64_t nEntries = calc_pot_chain.GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    calc_pot_chain.GetEntry(i);
    sum_pot += pot_value;
  }

  std::cout << "Summed POT: " << sum_pot << std::endl;
  return sum_pot;
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
  if (CheckNodeExists(SampleManager->raw(), "POT")) {
    pot = SampleManager->raw()["POT"].as<double>();
  } else {
    MACH3LOG_ERROR("POT not defined in {}, please add this!",
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

  // RegisterIndividualFunctionalParameter(
  //     "TotalEScaleND", kTotalEScaleND,
  //     [this](const double *par, std::size_t iEvent) {
  //       if (dunemcSamples[iEvent].shift.erec < 0.)
  //         dunemcSamples[iEvent].shift.erec = 0.;
  //       if (dunemcSamples[iEvent].reco.ELep < 0.)
  //         dunemcSamples[iEvent].reco.ELep = 0.;

  //       if (abs(dunemcSamples[iEvent].LepPDG == 11)) {
  //         dunemcSamples[iEvent].shift.erec *= (1.0 + (*par));
  //         dunemcSamples[iEvent].reco.ELep *= (1.0 + (*par));
  //       }

  //       dunemcSamples[iEvent].shift.erec +=
  //           (*par) * dunemcSamples[iEvent].reco.EHad;
  //       // dunemcSamples[iEvent].reco.ELep += (*par) *
  //       // dunemcSamples[iEvent].reco.EHad;
  //     });

  // RegisterIndividualFunctionalParameter(
  //     "TotalEScaleND_mu", kTotalEScaleND_mu,
  //     [this](const double *par, std::size_t iEvent) {
  //       if (abs(dunemcSamples[iEvent].LepPDG == 13) &&
  //           (dunemcSamples[iEvent].reco.muon_contained == 1 ||
  //            dunemcSamples[iEvent].reco.muon_tracker == 1)) {
  //         dunemcSamples[iEvent].shift.erec +=
  //             (*par) * (dunemcSamples[iEvent].reco.ELep);
  //         dunemcSamples[iEvent].reco.ELep *= (1.0 + (*par));
  //       }
  //     });

  // RegisterIndividualFunctionalParameter(
  //     "EMResND", kEMResND, [this](const double *par, std::size_t iEvent) {
  //       if (dunemcSamples[iEvent].reco.ePi0 < 0.0)
  //         dunemcSamples[iEvent].reco.ePi0 = 0.0;
  //       dunemcSamples[iEvent].shift.erec +=
  //           (*par) * (dunemcSamples[iEvent].truth.ePi0 -
  //                     dunemcSamples[iEvent].reco.ePi0);
  //       dunemcSamples[iEvent].shift.etrue +=
  //           (*par) * (dunemcSamples[iEvent].truth.ePi0 -
  //                     dunemcSamples[iEvent].reco.ePi0);
  //       dunemcSamples[iEvent].reco.ePi0 +=
  //           (*par) * (dunemcSamples[iEvent].truth.ePi0 -
  //                     dunemcSamples[iEvent].reco.ePi0);
  //       if (abs(dunemcSamples[iEvent].LepPDG == 11)) {
  //         dunemcSamples[iEvent].shift.erec +=
  //             (*par) * (dunemcSamples[iEvent].truth.LepE -
  //                       dunemcSamples[iEvent].reco.ELep);
  //         dunemcSamples[iEvent].reco.ELep +=
  //             (*par) * (dunemcSamples[iEvent].truth.LepE -
  //                       dunemcSamples[iEvent].reco.ELep);
  //         dunemcSamples[iEvent].shift.etrue +=
  //             (*par) * (dunemcSamples[iEvent].truth.LepE -
  //                       dunemcSamples[iEvent].reco.ELep);
  //       }
  //     });

  // RegisterIndividualFunctionalParameter(
  //     "EScaleMuSpectND", kEScaleMuSpectND,
  //     [this](const double *par, std::size_t iEvent) {
  //       if (abs(dunemcSamples[iEvent].LepPDG == 13) &&
  //           dunemcSamples[iEvent].reco.muon_tracker == 1) {
  //         dunemcSamples[iEvent].shift.erec +=
  //             (*par) * (dunemcSamples[iEvent].reco.ELep);
  //         dunemcSamples[iEvent].reco.ELep *= ((*par) + 1.0);
  //       }
  //     });

  // RegisterIndividualFunctionalParameter(
  //     "TotalEScaleND_hadsqrt", kTotalEScaleND_hadsqrt,
  //     [this](const double *par, std::size_t iEvent) {
  //       dunemcSamples[iEvent].shift.erec +=
  //           (*par) * dunemcSamples[iEvent].syst.sum_ehad_sqrt;
  //     });

  // RegisterIndividualFunctionalParameter(
  //     "TotalEScaleND_hadinvsqrt", kTotalEScaleND_hadinvsqrt,
  //     [this](const double *par, std::size_t iEvent) {
  //       dunemcSamples[iEvent].shift.erec +=
  //           (*par) * dunemcSamples[iEvent].syst.sum_ehad_sqrt /
  //           dunemcSamples[iEvent].reco.sum_ehad;
  //     });

  // RegisterIndividualFunctionalParameter(
  //     "TotalEScaleND_EM", kTotalEScaleND_EM,
  //     [this](const double *par, std::size_t iEvent) {
  //       dunemcSamples[iEvent].shift.erec +=
  //           (*par) * (dunemcSamples[iEvent].reco.ePi0 -
  //                     dunemcSamples[iEvent].truth.ePi0);
  //     });

  // RegisterIndividualFunctionalParameter(
  //     "TotalEScaleND_EM_invsqrt", kTotalEScaleND_EMinvsqrt,
  //     [this](const double *par, std::size_t iEvent) {
  //       dunemcSamples[iEvent].shift.erec +=
  //           (*par) *
  //           sqrt((dunemcSamples[iEvent].reco.ePi0 -
  //                 dunemcSamples[iEvent].truth.ePi0)) /
  //           ((dunemcSamples[iEvent].reco.ePi0 -
  //             dunemcSamples[iEvent].truth.ePi0));
  //     });

  // RegisterIndividualFunctionalParameter(
  //     "TotalEScaleND_EM_sqrt", kTotalEScaleND_EMsqrt,
  //     [this](const double *par, std::size_t iEvent) {
  //       dunemcSamples[iEvent].shift.erec +=
  //           (*par) * sqrt((dunemcSamples[iEvent].reco.ePi0 -
  //                          dunemcSamples[iEvent].truth.ePi0));
  //     });

  // RegisterIndividualFunctionalParameter(
  //     "MuonRes_ND", kMuonRes_ND, [this](const double *par, std::size_t iEvent) {
  //       if (abs(dunemcSamples[iEvent].LepPDG == 13)) {
  //         dunemcSamples[iEvent].shift.erec +=
  //             (*par) * (dunemcSamples[iEvent].truth.LepE -
  //                       dunemcSamples[iEvent].reco.ELep);
  //         dunemcSamples[iEvent].reco.ELep +=
  //             (*par) * (dunemcSamples[iEvent].truth.LepE -
  //                       dunemcSamples[iEvent].reco.ELep);
  //       }
  //     });

  // RegisterIndividualFunctionalParameter(
  //     "NRes_ND", kNRes_ND, [this](const double *par, std::size_t iEvent) {
  //       if (dunemcSamples[iEvent].reco.eN < 0.0)
  //         dunemcSamples[iEvent].reco.eN = 0.0;
  //       dunemcSamples[iEvent].shift.erec +=
  //           (*par) *
  //           (dunemcSamples[iEvent].truth.eN - dunemcSamples[iEvent].reco.eN);
  //       dunemcSamples[iEvent].reco.eN +=
  //           (*par) *
  //           (dunemcSamples[iEvent].truth.eN - dunemcSamples[iEvent].reco.eN);
  //     });

  // RegisterIndividualFunctionalParameter(
  //     "HadRes_ND", kHadRes_ND, [this](const double *par, std::size_t iEvent) {
  //       dunemcSamples[iEvent].shift.erec +=
  //           (*par) * ((dunemcSamples[iEvent].truth.ePim -
  //                      dunemcSamples[iEvent].reco.ePim) +
  //                     ((dunemcSamples[iEvent].truth.ePip -
  //                       dunemcSamples[iEvent].reco.ePip)) +
  //                     ((dunemcSamples[iEvent].truth.eP -
  //                       dunemcSamples[iEvent].reco.eP)));
  //       dunemcSamples[iEvent].reco.eP +=
  //           (*par) *
  //           (dunemcSamples[iEvent].truth.eP - dunemcSamples[iEvent].reco.eP);
  //       dunemcSamples[iEvent].reco.ePip +=
  //           (*par) * (dunemcSamples[iEvent].truth.ePip -
  //                     dunemcSamples[iEvent].reco.ePip);
  //       dunemcSamples[iEvent].reco.ePim +=
  //           (*par) * (dunemcSamples[iEvent].truth.ePim -
  //                     dunemcSamples[iEvent].reco.ePim);
  //     });
  RegisterIndividualFunctionalParameter(
      "TotalEScaleND", kTotalEScaleND,
      [this](const double *par, std::size_t iEvent) {

      if ( dunemcSamples[iEvent].shift.erec < 0.)  dunemcSamples[iEvent].shift.erec = 0.;
      if (dunemcSamples[iEvent].reco.ELep < 0.) dunemcSamples[iEvent].reco.ELep = 0.;

      if(abs(dunemcSamples[iEvent].LepPDG == 11)){
        dunemcSamples[iEvent].shift.erec *= (1.0 + (*par));

        dunemcSamples[iEvent].reco.ELep *= (1.0 + (*par));
      }
        dunemcSamples[iEvent].shift.erec += (*par) * dunemcSamples[iEvent].reco.EHad;
        //dunemcSamples[iEvent].reco.ELep += (*par) * dunemcSamples[iEvent].reco.EHad;


      });




  RegisterIndividualFunctionalParameter(
      "TotalEScaleND_mu", kTotalEScaleND_mu,
      [this](const double *par, std::size_t iEvent) {

        if(abs(dunemcSamples[iEvent].LepPDG == 13) && (dunemcSamples[iEvent].reco.muon_contained == 1 ||dunemcSamples[iEvent].reco.muon_tracker==1) ){
          dunemcSamples[iEvent].shift.erec += (*par) * (dunemcSamples[iEvent].reco.ELep);
           //dunemcSamples[iEvent].shift.etrue += (*par) * (dunemcSamples[iEvent].reco.ELep);
          dunemcSamples[iEvent].reco.ELep *= (1.0 + (*par));
        }

      });

    RegisterIndividualFunctionalParameter(
      "EMResND", kEMResND,
      [this](const double *par, std::size_t iEvent) {
        if(dunemcSamples[iEvent].reco.ePi0 < 0.0) dunemcSamples[iEvent].reco.ePi0 = 0.0;
        dunemcSamples[iEvent].shift.erec += (*par) * (dunemcSamples[iEvent].truth.ePi0 - dunemcSamples[iEvent].reco.ePi0);

        dunemcSamples[iEvent].reco.ePi0 += (*par) * (dunemcSamples[iEvent].truth.ePi0 - dunemcSamples[iEvent].reco.ePi0);
        if(abs(dunemcSamples[iEvent].LepPDG == 11)){
          dunemcSamples[iEvent].shift.erec += (*par) * (dunemcSamples[iEvent].truth.LepE - dunemcSamples[iEvent].reco.ELep);
          dunemcSamples[iEvent].reco.ELep += (*par) * (dunemcSamples[iEvent].truth.LepE - dunemcSamples[iEvent].reco.ELep);


        }
      }
    );



    RegisterIndividualFunctionalParameter(
      "EScaleMuSpectND", kEScaleMuSpectND,
      [this](const double *par, std::size_t iEvent) {
        if(abs(dunemcSamples[iEvent].LepPDG == 13) && dunemcSamples[iEvent].reco.muon_tracker==1){
          dunemcSamples[iEvent].shift.erec += (*par) * (dunemcSamples[iEvent].reco.ELep);
           dunemcSamples[iEvent].reco.ELep *= ((*par) + 1.0);


        }
    });


  RegisterIndividualFunctionalParameter(
      "MuonRes_ND", kMuonRes_ND, [this](const double *par, std::size_t iEvent) {
        if(abs(dunemcSamples[iEvent].LepPDG == 13)){
        dunemcSamples[iEvent].shift.erec += (*par) * (dunemcSamples[iEvent].truth.LepE - dunemcSamples[iEvent].reco.ELep);
        dunemcSamples[iEvent].reco.ELep += (*par) * (dunemcSamples[iEvent].truth.LepE - dunemcSamples[iEvent].reco.ELep);
        }
      });

  RegisterIndividualFunctionalParameter(
      "NRes_ND", kNRes_ND, [this](const double *par, std::size_t iEvent) {
        if(dunemcSamples[iEvent].reco.eN < 0.0) dunemcSamples[iEvent].reco.eN = 0.0;
        dunemcSamples[iEvent].shift.erec += (*par) * (dunemcSamples[iEvent].truth.eN - dunemcSamples[iEvent].reco.eN);
        dunemcSamples[iEvent].reco.eN += (*par) * (dunemcSamples[iEvent].truth.eN - dunemcSamples[iEvent].reco.eN);
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
        dunemcSamples[iEvent].reco.eP += (*par) *(dunemcSamples[iEvent].truth.eP -
                        dunemcSamples[iEvent].reco.eP) ;
        dunemcSamples[iEvent].reco.ePip += (*par) *(dunemcSamples[iEvent].truth.ePip -
                        dunemcSamples[iEvent].reco.ePip) ;
        dunemcSamples[iEvent].reco.ePim += (*par) * (dunemcSamples[iEvent].truth.ePim -
                       dunemcSamples[iEvent].reco.ePim);


      });


  /// Fake Data Syst
  RegisterIndividualFunctionalParameter(
      "NuWroFakeDataWeight", kNuWroFakeDataWeight,
      [this](const double *par, std::size_t iEvent) {
        dunemcSamples[iEvent].flux_w *=
            (((*par) * (this->NuWroFakeDataWeight(iEvent) - 1.0)) + 1);
      });

  // don't register flux parameters if we're not using them.
  if (ParHandler->GetNumParFromGroup("Flux")) {
    for (size_t par_it = 0;
         par_it < OffAxisFluxUncertaintyHelper::Get().GetNFocussingParams();
         par_it++) {
      RegisterIndividualFunctionalParameter(
          OffAxisFluxUncertaintyHelper::Get().GetFocussingParamName(par_it),
          int(kNFuncPars + par_it),
          [this, par_it](const double *par, std::size_t iEvent) {
            dunemcSamples[iEvent].flux_w *=
                (1 +
                 (*par) *
                     dunemcSamples[iEvent].syst.flux_focussing_ratio[par_it]);
          });
    }

    for (size_t par_it = 0;
         par_it <
         OffAxisFluxUncertaintyHelper::Get().GetNHadProdPCAComponents();
         par_it++) {
      RegisterIndividualFunctionalParameter(
          "Flux_HadProd_Param_" + std::to_string(par_it),
          int(kNFuncPars +
              OffAxisFluxUncertaintyHelper::Get().GetNFocussingParams() +
              par_it),
          [this, par_it](const double *par, std::size_t iEvent) {
            dunemcSamples[iEvent].flux_w *=
                (1 + (*par) *
                         dunemcSamples[iEvent].syst.flux_hadprod_ratio[par_it]);
          });
    }
  }

  MACH3LOG_INFO("Finished registering functional parameters");
}

void PrintFluxParameterNames() {
  auto &helper = OffAxisFluxUncertaintyHelper::Get();

  std::cout << "=== Focusing Flux Params ===\n";
  for (size_t i = 0; i < helper.GetNFocussingParams(); i++) {
    std::cout << helper.GetFocussingParamName(i) << "\n";
  }

  std::cout << "=== Hadron-Production PCA Params ===\n";
  for (size_t i = 0; i < helper.GetNHadProdPCAComponents(); i++) {
    std::cout << "Flux_HadProd_Param_" << i << "\n";
  }
}

// HH: Reset the shifted values to the original values
void SampleHandlerBeamOffAxis::resetShifts(int iEvent) {
  dunemcSamples[iEvent].shift.erec = dunemcSamples[iEvent].rw_erec;
  dunemcSamples[iEvent].reco.ELep = dunemcSamples[iEvent].reco.ELep;
  dunemcSamples[iEvent].flux_w = 1.0;
}
// =================================

void SampleHandlerBeamOffAxis::SetupWeightPointers() {
  // std::cout<< "WEIGHT POINTERS = " << std::endl;
  for (size_t i = 0; i < dunemcSamples.size(); ++i) {
    MCSamples[i].total_weight_pointers.push_back(&(dunemcSamples[i].pot_s));
    MCSamples[i].total_weight_pointers.push_back(&(dunemcSamples[i].norm_s));
    MCSamples[i].total_weight_pointers.push_back(MCSamples[i].osc_w_pointer);
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
  double ml = [](int _nupdg, bool _isCC) {
    switch (std::abs(_nupdg)) {
    case 12:
      return _isCC ? 0.511 : 0;
    case 14:
      return _isCC ? 105.66 : 0;
    case 16:
      return _isCC ? 1777.0 : 0;
    default:
      std::cerr << "Warning: Unexpected PDG code " << _nupdg
                << " passed to ml lambda.\n";
      assert(false && "Unexpected neutrino PDG code in ml lambda");
      return 0.0;
    }
  }(nupdg, isCC);

  double pl = std::sqrt(el * el - ml * ml); // momentum of lepton

  return (2 * mN_eff * el - ml * ml + mN_oth * mN_oth - mN_eff * mN_eff) /
         (2 * (mN_eff - el + pl * std::cos(leptheta)));
}

int SampleHandlerBeamOffAxis::SetupExperimentMC() {

  MACH3LOG_INFO(
      "-------------------------------------------------------------------");
  TChain CAFChain("cafTree");
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
  }

  MACH3LOG_INFO("Number of entries in CAF TChain: {}", CAFChain.GetEntries());

  // Calculate POT from the CAF files
  double cafpot = CalculatePOT(); // sums the 'pot' branch in 'meta' tree

  if (cafpot <= 0) {
    MACH3LOG_WARN("CalculatePOT() returned <= 0. Forcing cafpot = 1.");
    cafpot = 1.0; // avoid division by zero
  }

  // Compute POT scaling factor to normalize event weights
  double pot_scaling = pot / cafpot;

  std::cout << "POT in yaml file =  " << pot << ", Summed CAF POT: " << cafpot
            << ", Scaling factor = pot/cafpot: " << pot_scaling << std::endl;

  // Reco Variables
  TTreeReader chrdr(&CAFChain);

  TTreeReaderValue<double> Ev_reco(chrdr, "Ev_reco"); // Ev_reco
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
  TTreeReaderValue<double> Ev(chrdr, "Ev");

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

  TTreeReaderValue<double> Q2(chrdr, "Q2");
  TTreeReaderValue<double> W(chrdr, "W");
  TTreeReaderValue<double> X(chrdr, "X");
  TTreeReaderValue<double> Y(chrdr, "Y");

  TTreeReaderValue<int> nipi0(chrdr, "nipi0");

  TTreeReaderValue<int> isCC(chrdr, "isCC");

  TTreeReaderValue<int> nP(chrdr, "nP");
  TTreeReaderValue<int> nN(chrdr, "nN");
  TTreeReaderValue<int> nPi0(chrdr, "nipi0");
  TTreeReaderValue<int> nPim(chrdr, "nipim");
  TTreeReaderValue<int> nPip(chrdr, "nipip");
  TTreeReaderValue<int> niem(chrdr, "niem");

  TTreeReaderValue<int> nuPDG(chrdr, "nuPDG");
  TTreeReaderValue<int> nuPDGunosc(chrdr, "nuPDGunosc");
  TTreeReaderValue<int> LepPDG(chrdr, "LepPDG");

  dunemcSamples.resize(CAFChain.GetEntries());

  // HH: A map to keep track of negative energies
  std::unordered_map<std::string, int> negative_counts;
  // Initialize the negative counts for each energy variable
  negative_counts["rw_erec_had"] = 0;
  negative_counts["rw_erec_lep"] = 0;
  negative_counts["rw_eRecoN"] = 0;
  negative_counts["rw_eRecoPi0"] = 0;
  negative_counts["rw_sum_ehad"] = 0;

  /////////////////////////////////////

  size_t iEvent = 0;
  // std::cout<< "just before while..." << std::endl;


  bool do_flux_systs = ParHandler->GetNumParFromGroup("Flux");

  while (chrdr.Next()) {

    dunemcSamples[iEvent].norm_s = 1.0;
    dunemcSamples[iEvent].xsec_w = 1.0;
    dunemcSamples[iEvent].pot_s = pot_scaling;

    dunemcSamples[iEvent].Target = kTarget_Ar;

    dunemcSamples[iEvent].nupdg = *nuPDG;
    dunemcSamples[iEvent].nupdgUnosc = *nuPDGunosc;
    dunemcSamples[iEvent].LepPDG = *LepPDG;

    dunemcSamples[iEvent].rw_isCC = *isCC;

    dunemcSamples[iEvent].OscChannelIndex = 0;
    // static_cast<double>(
    //     GetOscChannel(OscChannels, dunemcSamples[iEvent].nupdgUnosc,
    //                   dunemcSamples[iEvent].nupdg));

    dunemcSamples[iEvent].rw_erec = *Ev_reco;
    dunemcSamples[iEvent].rw_etru = *Ev;
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
    dunemcSamples[iEvent].truth.LepE_sqrt = sqrt(*LepE);
    dunemcSamples[iEvent].truth.LepE_invsqrt = 1.0 / sqrt(*LepE);
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

    dunemcSamples[iEvent].truth.Q2 = *Q2;
    dunemcSamples[iEvent].truth.W = *W;
    dunemcSamples[iEvent].truth.X = *X;
    dunemcSamples[iEvent].truth.Y = *Y;

    dunemcSamples[iEvent].truth.nP = *nP;
    dunemcSamples[iEvent].truth.nN = *nN;
    dunemcSamples[iEvent].truth.nPip = *nPip;
    dunemcSamples[iEvent].truth.nPim = *nPim;
    dunemcSamples[iEvent].truth.nPi0 = *nPi0;

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

    dunemcSamples[iEvent].truth.off_axis_pos_m = (*det_x + *vtx_x) / 100.0;

    if (do_flux_systs) {
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
    }

    // HH: Give a warning if any negative energies were found
    for (const auto &[key, count] : negative_counts) {
      if (count > 0) {
        MACH3LOG_WARN("Found {} negative values for {} in sample {}", count,
                      key, GetSampleName());
      }
    }

    ++iEvent; // for pot?
  }

  if (do_flux_systs) {
    PrintFluxParameterNames();
  }
  return int(dunemcSamples.size());
}

// Fake Data Studies
float SampleHandlerBeamOffAxis::NuWroFakeDataWeight(std::size_t iEvent) {
  auto const &ev = dunemcSamples[iEvent];

  // struct features {
  //   int isCC;
  //   float Ev;
  //   float LepE;
  //   float LepNuAngle;
  //   float Q2;
  //   float W;
  //   float X;
  //   float Y;
  //   int nP;
  //   int nN;
  //   int nipip;
  //   int nipim;
  //   int nipi0;
  //   int niem;
  //   float eP;
  //   float eN;
  //   float ePip;
  //   float ePim;
  //   float ePi0;
  //   int isFD;
  //   int isFHC;
  //   int nuPDG;
  // };

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnarrowing"
#pragma GCC diagnostic ignored "-Wfloat-conversion"

  return WeightToNuWro::Weight(WeightToNuWro::features{
      ev.rw_isCC,    ev.rw_etru,    ev.truth.LepE, ev.truth.LepNuAngle,
      ev.truth.Q2,   ev.truth.W,    ev.truth.X,    ev.truth.Y,
      ev.truth.nP,   ev.truth.nN,   ev.truth.nPip, ev.truth.nPim,
      ev.truth.nPi0, ev.truth.niem, ev.truth.eP,   ev.truth.eN,
      ev.truth.ePip, ev.truth.ePim, ev.truth.ePi0, 0,
      isFHC,         ev.nupdg});
#pragma GCC diagnostic pop
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
  case kMode:
    return &dunemcSamples[iEvent].mode;
  case kIsFHC:
    return &isFHC;
  case kELepRec:
    return &dunemcSamples[iEvent].reco.ELep;
  case kEnubias:
    return &dunemcSamples[iEvent].truth.enu_bias;
  case kisCC:
    return &dunemcSamples[iEvent].rw_isCC;
  case kOscillationChannel:
    return &dunemcSamples[iEvent].OscChannelIndex;
  case kOffAxisPosition:
    return &dunemcSamples[iEvent].truth.off_axis_pos_m;
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

std::vector<std::vector<std::vector<std::vector<TH2D *>>>>
SampleHandlerBeamOffAxis::GetBinnedWeights(
    std::vector<std::string> ParamNames,
    std::vector<std::vector<int>> ParamModes, std::vector<double> TrueEBins) {

  std::cout << "Making Binned Weights" << std::endl;
  // Vector Structure:
  // Parameter<Knot<Mode<TrueE<TH2D>>>

  std::vector<std::vector<std::vector<std::vector<TH2D *>>>> histVec;
  std::vector<std::vector<std::vector<std::vector<TH2D *>>>> NomVec;
  int nshifts = 7;

  // True Energy Binning
  int NTrueEBins = static_cast<int>(TrueEBins.size()) - 1;
  TH1D *TrueEbinning =
      new TH1D("Template True E binning", "", NTrueEBins, TrueEBins.data());

  // Sample Binning
  std::vector<double> BinEdgesX =
      ReturnKinematicParameterBinning(GetXBinVarName());
  int NBinsX = static_cast<int>(BinEdgesX.size()) - 1;

  std::vector<double> BinEdgesY =
      ReturnKinematicParameterBinning(GetYBinVarName());
  int NBinsY = static_cast<int>(BinEdgesY.size()) - 1;

  // Setup Histograms
  histVec.resize(ParamNames.size());
  NomVec.resize(ParamNames.size());
  for (size_t iParam = 0; iParam < ParamNames.size(); iParam++) {
    histVec[iParam].resize(nshifts);
    NomVec[iParam].resize(nshifts);
    for (int shift = 0; shift < nshifts; shift++) {
      histVec[iParam][shift].resize(ParamModes[iParam].size());
      NomVec[iParam][shift].resize(ParamModes[iParam].size());
      for (size_t mode = 0; mode < ParamModes[iParam].size(); mode++) {
        histVec[iParam][shift][mode].resize(NTrueEBins);
        NomVec[iParam][shift][mode].resize(NTrueEBins);
        for (int b_etrue = 0; b_etrue < NTrueEBins; b_etrue++) {

          histVec[iParam][shift][mode][b_etrue] = new TH2D(
              Form("bn_p:%zu_s:%d_m:%zu_e:%d", iParam, shift, mode, b_etrue),
              "", NBinsX, BinEdgesX.data(), NBinsY, BinEdgesY.data());

          NomVec[iParam][shift][mode][b_etrue] = new TH2D(
              Form("nom_p:%zu_s:%d_m:%zu_e:%d", iParam, shift, mode, b_etrue),
              "", NBinsX, BinEdgesX.data(), NBinsY, BinEdgesY.data());
        }
      }
    }
  }

  // TChain to store MC which contains weight information
  TChain *MCData = new TChain("cafTree");

  // Read all files into TChain
  for (size_t iSample = 0; iSample < mc_files.size(); iSample++) {
    MACH3LOG_INFO("Adding file to TChains: {}", mc_files[iSample]);

    // HH: Check whether the file exists, see
    // https://root.cern/doc/master/classTChain.html#a78a896924ac6c7d3691b7e013bcbfb1c
    if (!MCData->Add(mc_files[iSample].c_str(), -1)) {
      MACH3LOG_ERROR("Could not add file {} to TChain, please check the file "
                     "exists and is readable",
                     mc_files[iSample]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    // MCData->Add(mc_files[iSample].c_str(), -1);
  }

  MACH3LOG_INFO("Number of entries in CAF TChain: {}", MCData->GetEntries());

  MCData->SetBranchStatus("*", 0);
  int NEvents = static_cast<int>(MCData->GetEntries());

  for (int iEvent = 0; iEvent < NEvents; iEvent++) {
    double x_var = ReturnKinematicParameter(GetXBinVarName(), iEvent);
    double y_var = ReturnKinematicParameter(GetYBinVarName(), iEvent);

    // skip if event does not pass selections
    if (!IsEventSelected(iEvent)) {
      continue;
    }

    int TrueEbin = TrueEbinning->FindBin(dunemcSamples[iEvent].rw_etru) - 1;

    for (size_t iParam = 0; iParam < ParamNames.size(); iParam++) {

      std::string FancyName = ParamNames[iParam];
      std::string WeightBranchName = "wgt_" + FancyName;
      std::string CVBranchName = FancyName + "_cvwgt";

      for (size_t mode = 0; mode < ParamModes[iParam].size(); mode++) {

        if (ParamModes[iParam][mode] == dunemcSamples[iEvent].mode) {

          double weightArr[1000];
          double cvweight;
          // Vector Structure:
          // Parameter<Knot<ETrue<InteractionMode<TH2D>>>

          MCData->SetBranchAddress(WeightBranchName.c_str(), &weightArr);
          MCData->SetBranchAddress(CVBranchName.c_str(), &cvweight);
          MCData->GetEntry(iEvent);

          for (int shift = 0; shift < nshifts; shift++) {
            histVec[iParam][shift][mode][TrueEbin]->Fill(x_var, y_var,
                                                         weightArr[shift]);
            NomVec[iParam][shift][mode][TrueEbin]->Fill(x_var, y_var, 1.0);
          }
        }
      }
    }
  }

  // Take the ratio
  for (size_t iParam = 0; iParam < ParamNames.size(); iParam++) {
    for (int shift; shift < nshifts; shift++) {
      for (size_t mode = 0; mode < ParamModes[iParam].size(); mode++) {
        for (int b_etrue = 0; b_etrue < NTrueEBins; b_etrue++) {
          histVec[iParam][shift][mode][b_etrue]->Divide(
              NomVec[iParam][shift][mode][b_etrue]);
        }
      }
    }
  }

  std::cout << "finished making binned weights" << std::endl;
  // return final vector of response histograms
  return histVec;
}
