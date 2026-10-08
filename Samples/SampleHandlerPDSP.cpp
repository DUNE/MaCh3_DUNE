#include "SampleHandlerPDSP.h"

#include <TFile.h>
#include <TKey.h>
#include <TString.h>

#include <algorithm>
#include <cstdint>
#include <cmath>

namespace {

// SplitMix64 gives each event a stable pseudo-random value without keeping a
// mutable RNG in the likelihood calculation.
std::uint64_t SplitMix64(std::uint64_t value) {
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

double UniformFromHash(const std::uint64_t value) {
  // Use the upper 53 bits and offset by half a bin so log(u) is always safe.
  return (static_cast<double>(value >> 11U) + 0.5) *
         (1.0 / 9007199254740992.0);
}

double FixedGaussianDeviate(const std::uint64_t eventIndex) {
  const double u1 = UniformFromHash(SplitMix64(2U * eventIndex));
  const double u2 = UniformFromHash(SplitMix64(2U * eventIndex + 1U));
  constexpr double twoPi = 6.28318530717958647692;
  return std::sqrt(-2.0 * std::log(u1)) * std::cos(twoPi * u2);
}

double PionMeanDEDX(const double kineticEnergy) {
  if (!std::isfinite(kineticEnergy) || kineticEnergy <= 0.0) return 0.0;

  // Liquid-argon Bethe--Bloch constants used by the analysis production.
  constexpr double rho = 1.39;
  constexpr double K = 0.307075;
  constexpr double Z = 18.0;
  constexpr double atomicMass = 39.948;
  constexpr double excitationEnergy = 188.0e-6;
  constexpr double electronMass = 0.51099895;
  constexpr double pionMass = 139.57039;
  constexpr double densityC = 5.2146;
  constexpr double densityY0 = 0.2;
  constexpr double densityY1 = 3.0;
  constexpr double densityA = 0.19559;
  constexpr double densityK = 3.0;

  const double gamma = kineticEnergy / pionMass + 1.0;
  const double beta2 = 1.0 - 1.0 / (gamma * gamma);
  if (beta2 <= 0.0) return 0.0;

  const double betaGamma = std::sqrt(beta2) * gamma;
  const double massRatio = electronMass / pionMass;
  const double wMax = 2.0 * electronMass * betaGamma * betaGamma /
      (1.0 + 2.0 * massRatio * gamma + massRatio * massRatio);

  const double y = std::log10(betaGamma);
  double densityCorrection = 0.0;
  const double delta0 = 2.0 * std::log(10.0) * y - densityC;
  if (y >= densityY1) {
    densityCorrection = delta0;
  } else if (y >= densityY0) {
    densityCorrection =
        delta0 + densityA * std::pow(densityY1 - y, densityK);
  }

  const double normalisation = rho * K * Z / (atomicMass * beta2);
  const double logTerm = 0.5 * std::log(
      2.0 * electronMass * gamma * gamma * beta2 * wMax /
      (excitationEnergy * excitationEnergy));
  const double dedx = normalisation *
      (logTerm - beta2 - 0.5 * densityCorrection);
  return std::isfinite(dedx) && dedx > 0.0 ? dedx : 0.0;
}

double PropagatePionKE(const double initialKE, const double length,
                       const int numberOfSteps) {
  if (!std::isfinite(initialKE) || initialKE <= 0.0) return 0.0;
  if (!std::isfinite(length) || length <= 0.0 || numberOfSteps <= 0) {
    return initialKE;
  }

  const double stepLength = length / numberOfSteps;
  double kineticEnergy = initialKE;
  for (int step = 0; step < numberOfSteps && kineticEnergy > 0.0; ++step) {
    const double dedx = PionMeanDEDX(kineticEnergy);
    if (dedx <= 0.0) break;
    kineticEnergy = std::max(0.0, kineticEnergy - dedx * stepLength);
  }
  return kineticEnergy;
}

} // namespace

// ************************************************
SampleHandlerPDSP::SampleHandlerPDSP(const std::string& config_name, ParameterHandlerGeneric* parameter_handler)
    : SampleHandlerBase(config_name, parameter_handler), MCGlobalScale(1.0) {
// ************************************************
  KinematicParameters = &KinematicParametersPDSP;
  ReversedKinematicParameters = &ReversedKinematicParametersPDSP;

  Initialise();
}

// ************************************************
SampleHandlerPDSP::~SampleHandlerPDSP() {
// ************************************************

}

// ************************************************
TH1* SampleHandlerPDSP::GetDataHistogramFromInputs(const int Sample) const {
// ************************************************
  if (Sample < 0 || Sample >= static_cast<int>(SampleDetails.size())) {
    MACH3LOG_ERROR("Requested data histogram for invalid sample index {}", Sample);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  const std::string expectedName = SampleDetails[Sample].SampleTitle + "_DataHist";

  std::vector<std::string> candidateNames = {expectedName};
  if (SampleDetails[Sample].SampleTitle == "PDSP_Abs") {
    candidateNames.emplace_back("absorption_DataHist");
  } else if (SampleDetails[Sample].SampleTitle == "PDSP_CEx") {
    candidateNames.emplace_back("charge_exchange_DataHist");
  } else if (SampleDetails[Sample].SampleTitle == "PDSP_Pip") {
    candidateNames.emplace_back("pion_production_DataHist");
  } else if (SampleDetails[Sample].SampleTitle == "PDSP_Uncategorised") {
    candidateNames.emplace_back("uncategorised_DataHist");
  }

  for (const auto& fileGroup : SampleDetails[Sample].mc_files) {
    for (const auto& fileName : fileGroup) {
      TFile inputFile(fileName.c_str(), "READ");
      if (inputFile.IsZombie()) {
        MACH3LOG_ERROR("Could not open input file while looking for data histogram: {}", fileName);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      for (const auto& candidateName : candidateNames) {
        TH1* dataHist = nullptr;
        inputFile.GetObject(candidateName.c_str(), dataHist);
        if (dataHist != nullptr) {
          TH1* clone = static_cast<TH1*>(dataHist->Clone(expectedName.c_str()));
          clone->SetDirectory(nullptr);
          MACH3LOG_INFO("Loaded data histogram '{}' from {} for sample '{}'",
                        candidateName, fileName, SampleDetails[Sample].SampleTitle);
          return clone;
        }
      }

      TIter nextKey(inputFile.GetListOfKeys());
      while (TKey* key = static_cast<TKey*>(nextKey())) {
        const TString keyName = key->GetName();
        if (!keyName.EndsWith("_DataHist")) {
          continue;
        }

        TObject* object = key->ReadObj();
        if (!object->InheritsFrom(TH1::Class())) {
          delete object;
          continue;
        }

        TH1* clone = static_cast<TH1*>(static_cast<TH1*>(object)->Clone(expectedName.c_str()));
        clone->SetDirectory(nullptr);
        delete object;
        MACH3LOG_INFO("Loaded data histogram '{}' from {} for sample '{}'",
                      keyName.Data(), fileName, SampleDetails[Sample].SampleTitle);
        return clone;
      }
    }
  }

  MACH3LOG_ERROR("Could not find '{}' in the input files for sample '{}'",
                 expectedName, SampleDetails[Sample].SampleTitle);
  MACH3LOG_ERROR("Expected each PDSP input to contain a histogram ending in '_DataHist'.");
  throw MaCh3Exception(__FILE__, __LINE__);
}

// ************************************************
void SampleHandlerPDSP::Init() {
// ************************************************
  MCGlobalScale = GetFromManager<double>(SampleManager->raw()["MCGlobalScale"], 1.0, __FILE__, __LINE__);
  MACH3LOG_INFO("PDSP MC global scale: {}", MCGlobalScale);
}

// ************************************************
void SampleHandlerPDSP::SetupSplines() {
// ************************************************

}

// ************************************************
void SampleHandlerPDSP::AddAdditionalWeightPointers() {
// ************************************************
  for (std::size_t iEvent = 0; iEvent < MCEvents.size(); ++iEvent) {
    MCEvents[iEvent].total_weight_pointers.push_back(&MCGlobalScale);
    MCEvents[iEvent].total_weight_pointers.push_back(&PDSPSamples[iEvent].BeamMomentumWeight);
  }
}

// ************************************************
void SampleHandlerPDSP::InititialiseData() {
// ************************************************
  Reweight();
  for (int iSample = 0; iSample < GetNSamples(); ++iSample) {
    AddData(iSample, GetMCArray(iSample));
  }
}

void SampleHandlerPDSP::CleanMemoryBeforeFit() {
  CleanVector(PDSPPlottingSamples);
}

// ************************************************
int SampleHandlerPDSP::SetupExperimentMC() {
// ************************************************

  // *** Get number of events from each file
  TChain* _Chain = new TChain("FlatTree_VARS");
  for(size_t iSample = 0; iSample < SampleDetails.size(); iSample++)
  {
    for (const auto& fileGroup : SampleDetails[iSample].mc_files) {
      for (const auto& filename : fileGroup) {
        _Chain->Add(filename.c_str());
      }
    }
  }

  int nEntries = static_cast<int>(_Chain->GetEntries());
  delete _Chain;
  // ***

  // Set size of data vectors
  PDSPSampleMetaData.resize(nEntries);
  PDSPSamples.resize(nEntries);
  PDSPPlottingSamples.resize(nEntries);

  int TotalEventCounter = 0;
  // loop over all Samples
  for(size_t iSample = 0; iSample < SampleDetails.size(); iSample++)
  {
    // loop over all samples in a file
    for(size_t iFile = 0; iFile < SampleDetails[iSample].mc_files.size(); iFile++)
    {
      for (const auto& fileName : SampleDetails[iSample].mc_files[iFile]) {
      MACH3LOG_INFO("-------------------------------------------------------------------");
      MACH3LOG_INFO("input file: {}", fileName);

      TFile* _sampleFile = new TFile(fileName.c_str(), "READ");
      TTree* _data = static_cast<TTree*>(_sampleFile->Get("FlatTree_VARS"));      

      if(_data){
        MACH3LOG_INFO("Found \"FlatTree_VARS\" tree in {}", fileName);
        MACH3LOG_INFO("With number of entries: {}", _data->GetEntries());
      } else{
        MACH3LOG_ERROR("Could not find \"FlatTree_VARS\" tree in {}", fileName);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      _data->SetBranchStatus("*", false);
      
      // Truth variables
      double trueKEIni;
      double trueKEInt;
      double trueEndZ;
      bool true_abs;
      bool true_cex;
      bool true_pip;
      bool true_decay;

      _data->SetBranchStatus("KE_ff_true", true);
      _data->SetBranchAddress("KE_ff_true", &trueKEIni);

      _data->SetBranchStatus("KE_int_true", true);
      _data->SetBranchAddress("KE_int_true", &trueKEInt);

      _data->SetBranchStatus("track_length_true", true);
      _data->SetBranchAddress("track_length_true", &trueEndZ);

      _data->SetBranchStatus("exclusive_process_absorption", true);
      _data->SetBranchAddress("exclusive_process_absorption", &true_abs);

      _data->SetBranchStatus("exclusive_process_charge_exchange", true);
      _data->SetBranchAddress("exclusive_process_charge_exchange", &true_cex);

      _data->SetBranchStatus("exclusive_process_pion_production", true);
      _data->SetBranchAddress("exclusive_process_pion_production", &true_pip);

      _data->SetBranchStatus("exclusive_process_decay", true);
      _data->SetBranchAddress("exclusive_process_decay", &true_decay);

      // Reco variables
      double recoKEFF;
      double recoKEIni;
      double recoKEInt;
      double recoEndZ;
      double recoPinst;
      double recoTrackLength;

      _data->SetBranchStatus("KE_ff_reco", true);
      _data->SetBranchAddress("KE_ff_reco", &recoKEFF);

      _data->SetBranchStatus("KE_init_reco", true);
      _data->SetBranchAddress("KE_init_reco", &recoKEIni);
      
      _data->SetBranchStatus("KE_int_reco", true);
      _data->SetBranchAddress("KE_int_reco", &recoKEInt);

      _data->SetBranchStatus("end_z_reco", true);
      _data->SetBranchAddress("end_z_reco", &recoEndZ);

      _data->SetBranchStatus("P_inst_reco", true);
      _data->SetBranchAddress("P_inst_reco", &recoPinst);

      _data->SetBranchStatus("track_length_reco", true);
      _data->SetBranchAddress("track_length_reco", &recoTrackLength);

      for (int i = 0; i < _data->GetEntries(); ++i) { // Loop through tree (events)
        _data->GetEntry(i);

        PDSPSampleMetaData[TotalEventCounter].SampleIndex = static_cast<int>(iSample);

        PDSPSamples[TotalEventCounter].TrueKEIni = trueKEIni;
        PDSPSamples[TotalEventCounter].TrueKEInt = trueKEInt;
        PDSPSamples[TotalEventCounter].TrueEndZ = trueEndZ;
        PDSPSamples[TotalEventCounter].RecoKEFF = recoKEFF;
        PDSPSamples[TotalEventCounter].RecoKEFFShifted = recoKEFF;
        PDSPSamples[TotalEventCounter].RecoKEIni = recoKEIni;
        PDSPSamples[TotalEventCounter].RecoKEInt = recoKEInt;
        PDSPSamples[TotalEventCounter].RecoEndZ = recoEndZ;
        PDSPSamples[TotalEventCounter].RecoPinst = recoPinst;
        PDSPSamples[TotalEventCounter].RecoPinstShifted = recoPinst;
        PDSPSamples[TotalEventCounter].RecoTrackLength = recoTrackLength;
        PDSPSamples[TotalEventCounter].RecoTrackLengthShifted = recoTrackLength;
        PDSPSamples[TotalEventCounter].TrackLengthSmearZ =
            FixedGaussianDeviate(static_cast<std::uint64_t>(TotalEventCounter));
        PDSPSamples[TotalEventCounter].RecoKEIniShifted = recoKEIni;
        PDSPSamples[TotalEventCounter].RecoKEIntShifted = recoKEInt;


        bool isPion = true_abs == 1 || true_cex == 1 || true_pip == 1 || true_decay == 1;
        int mode;
        if(trueEndZ > 220 && isPion) {
          mode = 4; // escaping pions
        }else if(true_abs == 1) {
          mode = 0;
        }else if(true_cex == 1) {
          mode = 1;
        }else if(true_pip == 1) {
          mode = 2;
        }else if(true_decay == 1) {
          mode = 3;
        }else {
          mode = 999;
        }
        // MACH3LOG_INFO("abs | {}, cex | {}, spip | {}, pip | {}, mode | {}", true_abs, true_cex, true_pip, mode);
        PDSPSamples[TotalEventCounter].Mode = mode;


        //? redundant?
        PDSPPlottingSamples[TotalEventCounter].TrueKEIni = trueKEIni;
        PDSPPlottingSamples[TotalEventCounter].TrueKEInt = trueKEInt;
        PDSPPlottingSamples[TotalEventCounter].RecoKEIni = recoKEIni;
        PDSPPlottingSamples[TotalEventCounter].RecoKEInt = recoKEInt;

        TotalEventCounter++;
      }
      _sampleFile->Close();
      delete _sampleFile;
      MACH3LOG_INFO("Initialised file: {}/{}", iFile, iSample);
      }
    }
  }
  return nEntries;
}

double SampleHandlerPDSP::ReturnKinematicParameter(const int KinematicVariable, const int iEvent) const {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iEvent);
}

const double* SampleHandlerPDSP::GetPointerToKinematicParameter(const int KinPar, const int iEvent) const {
  switch(KinPar) {
    case kTrueKEIni:
      return &PDSPSamples[iEvent].TrueKEIni;
    case kTrueKEInt:
      return &PDSPSamples[iEvent].TrueKEInt;
    case kRecoKEIni:
      return &PDSPSamples[iEvent].RecoKEIniShifted;
    case kRecoKEInt:
      return &PDSPSamples[iEvent].RecoKEIntShifted;
    case kMode:
      return &PDSPSamples[iEvent].Mode;
    case kOscChannel:
      return &PDSPSamples[iEvent].OscillationChannel;
    case kTargetNucleus:
      return &PDSPSamples[iEvent].Target;
    case kTrueEndZ:
      return &PDSPSamples[iEvent].TrueEndZ;
    case kRecoEndZ:
      return &PDSPSamples[iEvent].RecoEndZ;
    case kRecoPinst:
      return &PDSPSamples[iEvent].RecoPinstShifted;
    case kRecoTrackLength:
      return &PDSPSamples[iEvent].RecoTrackLengthShifted;
    default:
      MACH3LOG_ERROR("Unrecognized Kinematic Parameter type: {}", KinPar);
      throw MaCh3Exception(__FILE__, __LINE__);
  }
}

void SampleHandlerPDSP::SetupMC() {
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); ++iEvent) {
    MCEvents[iEvent].NominalSample = PDSPSampleMetaData[iEvent].SampleIndex;
  }
}

void SampleHandlerPDSP::RegisterFunctionalParameters() {
  MACH3LOG_INFO("Registering PDSP systematic response functions");

  // 1.2% beam-instrument momentum scale.  Shift KE at the TPC front face;
  // FinaliseShifts derives KE_init and KE_int from the completed KE_ff shift.
  RegisterIndividualFunctionalParameter(
      PDSPSamples, "BeamMomentumMeasurement",
      [](const double& theta, PDSPMCInfo& event) {
        const double shiftedP = event.RecoPinst * (1.0 + 0.012 * theta);
        const auto kineticEnergy = [](double p) {
          constexpr double mass = 139.57039;
          return std::sqrt(p * p + mass * mass) - mass;
        };
        const double deltaKE = kineticEnergy(shiftedP) - kineticEnergy(event.RecoPinst);
        event.RecoPinstShifted = shiftedP;
        event.RecoKEFFShifted += deltaKE;
      });

  // Beam sideband Gaussian.  Apply this spectrum correction before the
  // upstream-energy response, matching the analysis correction sequence.
  // A single standard-normal envelope parameter moves all three fit
  // coefficients coherently by their quoted one-sigma errors.
  RegisterIndividualFunctionalParameter(
      PDSPSamples, "BeamMomentumReweight",
      [](const double& theta, PDSPMCInfo& event) {
        const double p0 = 1.38 + 0.07 * theta;
        const double p1 = 2000.0 + 9.0 * theta;
        const double p2 = 150.0 + 9.0 * theta;
        if (p0 <= 0.0 || p2 <= 0.0) {
          event.BeamMomentumWeight = 1.8;
          return;
        }
        const double pull = (event.RecoPinstShifted - p1) / p2;
        const double ratio = p0 * std::exp(-0.5 * pull * pull);
        event.BeamMomentumWeight = (!std::isfinite(ratio) || ratio <= 0.0)
            ? 1.8
            : static_cast<M3::float_t>(std::min(1.8, 1.0 / ratio));
      });

  // Vary the externally fitted upstream correction relative to the central
  // correction already present in the nominal reconstructed energies.
  RegisterIndividualFunctionalParameter(
      PDSPSamples, "UpstreamEnergyCorrection",
      [](const double& theta, PDSPMCInfo& event) {
        constexpr double p0 = -25.0;
        constexpr double p1 = 0.23;
        constexpr double p2 = 2.75e-3;
        constexpr double e0 = 1.0;
        constexpr double e1 = 0.03;
        constexpr double e2 = 0.06e-3;
        const auto kineticEnergy = [](double p) {
          constexpr double mass = 139.57039;
          return std::sqrt(p * p + mass * mass) - mass;
        };
        const double shiftedX = kineticEnergy(event.RecoPinstShifted);
        // KE_ff_reco already contains the central upstream correction.  Use
        // the same current beam energy for the central and varied functions,
        // so theta=0 is exactly an identity even if the beam scale has moved.
        const double nominal = p0 + p1 * std::exp(p2 * shiftedX);
        const double varied = (p0 + theta * e0) +
                              (p1 + theta * e1) * std::exp((p2 + theta * e2) * shiftedX);
        const double deltaCorrection = varied - nominal;
        event.RecoKEFFShifted -= deltaCorrection;
      });

  // track_length_reco starts at the TPC front face.  This functional changes
  // only that length; FinaliseShifts propagates KE_ff over the shifted length.
  RegisterIndividualFunctionalParameter(
      PDSPSamples, "TrackLengthResolution",
      [](const double& theta, PDSPMCInfo& event) {
        const double scale = std::max(
            0.01, 1.0 + 0.026 * theta * event.TrackLengthSmearZ);
        event.RecoTrackLengthShifted = event.RecoTrackLength * scale;
      });

  MACH3LOG_INFO("Finished registering PDSP systematic response functions");
}

void SampleHandlerPDSP::ResetShifts(const int iEvent) {
  auto& event = PDSPSamples[iEvent];
  event.RecoPinstShifted = event.RecoPinst;
  event.RecoTrackLengthShifted = event.RecoTrackLength;
  event.RecoKEFFShifted = event.RecoKEFF;
  event.RecoKEIniShifted = event.RecoKEIni;
  event.RecoKEIntShifted = event.RecoKEInt;
  event.BeamMomentumWeight = 1.0;
}

void SampleHandlerPDSP::FinaliseShifts(const int iEvent) {
  auto& event = PDSPSamples[iEvent];

  // The reconstructed FV begins 30 cm downstream of the TPC front face.
  // Match the analysis production: 25 steps to the FV boundary and 50 steps
  // over the full reconstructed trajectory to the interaction point.
  constexpr double fiducialVolumeStart = 30.0;
  event.RecoKEIniShifted = PropagatePionKE(
      event.RecoKEFFShifted, fiducialVolumeStart, 25);
  event.RecoKEIntShifted = PropagatePionKE(
      event.RecoKEFFShifted, event.RecoTrackLengthShifted, 50);
}
