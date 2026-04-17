#include "SampleHandlerPDSP.h"

// ************************************************
SampleHandlerPDSP::SampleHandlerPDSP(const std::string& config_name, ParameterHandlerGeneric* parameter_handler)
                                             : SampleHandlerFD(config_name, parameter_handler) {
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
void SampleHandlerPDSP::Init() {
// ************************************************

}

// ************************************************
void SampleHandlerPDSP::SetupSplines() {
// ************************************************

}

// ************************************************
void SampleHandlerPDSP::AddAdditionalWeightPointers() {
// ************************************************
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
    for (const std::string& filename : SampleDetails[iSample].mc_files) {
      _Chain->Add(filename.c_str());
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
      auto fileName = SampleDetails[iSample].mc_files[iFile];
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
      double trueKEInt;
      bool true_abs;
      bool true_cex;
      bool true_spip;
      bool true_pip;

      _data->SetBranchStatus("KE_int_true", true);
      _data->SetBranchAddress("KE_int_true", &trueKEInt);

      _data->SetBranchStatus("exclusive_process_absorption", true);
      _data->SetBranchAddress("exclusive_process_absorption", &true_abs);

      _data->SetBranchStatus("exclusive_process_charge_exchange", true);
      _data->SetBranchAddress("exclusive_process_charge_exchange", &true_cex);

      _data->SetBranchStatus("exclusive_process_single_pion_production", true);
      _data->SetBranchAddress("exclusive_process_single_pion_production", &true_spip);

      _data->SetBranchStatus("exclusive_process_pion_production", true);
      _data->SetBranchAddress("exclusive_process_pion_production", &true_pip);

      // Reco variables
      double recoKEInt;
      // int reco_abs;
      // int reco_cex;
      // int reco_spip;
      // int reco_pip;

      _data->SetBranchStatus("KE_int_reco", true);
      _data->SetBranchAddress("KE_int_reco", &recoKEInt);

      for (int i = 0; i < _data->GetEntries(); ++i) { // Loop through tree (events)
        _data->GetEntry(i);

        PDSPSampleMetaData[TotalEventCounter].SampleIndex = static_cast<int>(iSample);

        PDSPSamples[TotalEventCounter].TrueKEInt = trueKEInt;
        PDSPSamples[TotalEventCounter].RecoKEInt = recoKEInt;

        int mode;
        if(true_abs == 1) {
          mode = 0;
        }else if(true_cex == 1) {
          mode = 1;
        }else if(true_spip == 1 || true_pip == 1) {
          mode = 2;
        }else {
          mode = -1;
        }
        // MACH3LOG_INFO("abs | {}, cex | {}, spip | {}, pip | {}, mode | {}", true_abs, true_cex, true_spip, true_pip, mode);
        PDSPSamples[TotalEventCounter].Mode = mode;


        //? redundant?
        PDSPPlottingSamples[TotalEventCounter].TrueKEInt = trueKEInt;
        PDSPPlottingSamples[TotalEventCounter].RecoKEInt = recoKEInt;

        TotalEventCounter++;
      }
      _sampleFile->Close();
      delete _sampleFile;
      MACH3LOG_INFO("Initialised file: {}/{}", iFile, iSample);
    }
  }
  return nEntries;
}

double SampleHandlerPDSP::ReturnKinematicParameter(KinematicTypes KinPar, int iEvent) {
  const double* paramPointer = GetPointerToKinematicParameter(KinPar, iEvent);
  return *paramPointer;
}

double SampleHandlerPDSP::ReturnKinematicParameter(int KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(std::round(KinematicVariable));
  return ReturnKinematicParameter(KinPar, iEvent);
}

double SampleHandlerPDSP::ReturnKinematicParameter(std::string KinematicParameter, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return ReturnKinematicParameter(KinPar, iEvent);
}

const double* SampleHandlerPDSP::GetPointerToKinematicParameter(KinematicTypes KinPar, int iEvent) {
  switch (KinPar) {
    case kTrueKEInt:
      return &PDSPSamples[iEvent].TrueKEInt;
    case kRecoKEInt:
      return &PDSPSamples[iEvent].RecoKEInt;
    case kMode: // required to work with SampleHandlerFD
      return &PDSPSamples[iEvent].Mode;
    case kOscChannel: // required to work with SampleHandlerFD
      return &PDSPSamples[iEvent].OscillationChannel;
    default:
      MACH3LOG_ERROR("Unrecognized Kinematic Parameter type: {}", static_cast<int>(KinPar));
      throw MaCh3Exception(__FILE__, __LINE__);
  }
}

const double* SampleHandlerPDSP::GetPointerToKinematicParameter(double KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(std::round(KinematicVariable));
  return GetPointerToKinematicParameter(KinPar, iEvent);
}

const double* SampleHandlerPDSP::GetPointerToKinematicParameter(std::string KinematicParameter, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return GetPointerToKinematicParameter(KinPar, iEvent);
}

void SampleHandlerPDSP::SetupFDMC() { 
  for (unsigned int iEvent = 0; iEvent < GetNEvents(); ++iEvent) {
    MCSamples[iEvent].NominalSample = PDSPSampleMetaData[iEvent].SampleIndex;
  }
}

void SampleHandlerPDSP::RegisterFunctionalParameters() {
  MACH3LOG_INFO("No functional parameters");

}
