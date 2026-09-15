#include "SampleHandlerPDSP.h"

#include <TFile.h>
#include <TKey.h>
#include <TString.h>

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
  for (auto& sample : MCEvents) {
    sample.total_weight_pointers.push_back(&MCGlobalScale);
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
      double recoKEIni;
      double recoKEInt;
      double recoEndZ;

      _data->SetBranchStatus("KE_ff_reco", true);
      _data->SetBranchAddress("KE_ff_reco", &recoKEIni);
      
      _data->SetBranchStatus("KE_int_reco", true);
      _data->SetBranchAddress("KE_int_reco", &recoKEInt);

      _data->SetBranchStatus("end_z_reco", true);
      _data->SetBranchAddress("end_z_reco", &recoEndZ);

      for (int i = 0; i < _data->GetEntries(); ++i) { // Loop through tree (events)
        _data->GetEntry(i);

        PDSPSampleMetaData[TotalEventCounter].SampleIndex = static_cast<int>(iSample);

        PDSPSamples[TotalEventCounter].TrueKEIni = trueKEIni;
        PDSPSamples[TotalEventCounter].TrueKEInt = trueKEInt;
        PDSPSamples[TotalEventCounter].TrueEndZ = trueEndZ;
        PDSPSamples[TotalEventCounter].RecoKEIni = recoKEIni;
        PDSPSamples[TotalEventCounter].RecoKEInt = recoKEInt;
        PDSPSamples[TotalEventCounter].RecoEndZ = recoEndZ;


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
      return &PDSPSamples[iEvent].RecoKEIni;
    case kRecoKEInt:
      return &PDSPSamples[iEvent].RecoKEInt;
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
  MACH3LOG_INFO("No functional parameters");

}
