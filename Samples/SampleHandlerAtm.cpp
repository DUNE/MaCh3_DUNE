#include "SampleHandlerAtm.h"

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wfloat-conversion"
//Standard Record includes
#include "duneanaobj/StandardRecord/StandardRecord.h"
#pragma GCC diagnostic pop

SampleHandlerAtm::SampleHandlerAtm(std::string mc_version_, ParameterHandlerGeneric* xsec_cov_, const std::shared_ptr<OscillationHandler>&  Oscillator_) : SampleHandlerFD(mc_version_, xsec_cov_, Oscillator_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  
  Initialise();
}

SampleHandlerAtm::~SampleHandlerAtm() {
}

void SampleHandlerAtm::Init() {
  std::vector<std::string> EnabledSamples = Get<std::vector<std::string>>(SampleManager->raw()["Samples"], __FILE__ , __LINE__);
  IsELike.resize(GetNsamples());
  
  ExposureScaling = Get<double>(SampleManager->raw()["AnalysisOptions"]["ExposureScaling"],__FILE__,__LINE__);
  for(int iSample=0;iSample<GetNsamples();iSample++){
    const std::string TempTitle = EnabledSamples[iSample];
    IsELike[iSample] = Get<int>(SampleManager->raw()[TempTitle]["SampleOptions"]["IsELike"],__FILE__,__LINE__);
  }
  
}

void SampleHandlerAtm::SetupSplines() {
  SplineHandler = nullptr;
}

void SampleHandlerAtm::AddAdditionalWeightPointers() {
  for (size_t i = 0; i < dunemcSamples.size(); ++i) {
    MCSamples[i].total_weight_pointers.push_back(&(dunemcSamples[i].flux_w));
    MCSamples[i].total_weight_pointers.push_back(&(ExposureScaling));
  }  
}

int SampleHandlerAtm::SetupExperimentMC() {
  int CurrErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  
  caf::StandardRecord* sr = new caf::StandardRecord();

  TChain* Chain = new TChain("cafTree");
  std::vector<size_t> FileIndexToSample;
  for (size_t iSample=0;iSample<SampleDetails.size();iSample++) {
    for (const std::string& Filename : SampleDetails[iSample].mc_files) {
      FileIndexToSample.push_back(iSample);
      int ChainAddCheck = Chain->Add(Filename.c_str(), -1);
      if(ChainAddCheck == 0){
        MACH3LOG_ERROR("Could not add file {} to TChain, please check the file exists and is readable", Filename);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
  }
  int nChainEntries = static_cast<int>(Chain->GetEntries());
  
  Chain->SetBranchStatus("*", 1);
  Chain->SetBranchAddress("rec", &sr);

  //================================================================================================
  //Count how many "sensible" entries we have

  MACH3LOG_INFO("Counting sensible entries...");
  
  int nSensibleEntries = 0;
  for (int iChainEntry=0;iChainEntry<nChainEntries;iChainEntry++) {
    Chain->GetEntry(iChainEntry);

    if ((iChainEntry % (nChainEntries/10))==0) {
      MACH3LOG_INFO("\tProcessing entry: {}/{}",iChainEntry,nChainEntries);
    }
    
    if(sr->common.ixn.pandora.size() != 1) {
      continue;
    }

    const size_t SampleIndex = FileIndexToSample[static_cast<size_t>(Chain->GetTreeNumber())];
    
    TVector3 RecoNuMomentumVector;
    double RecoENu;
    if (IsELike[SampleIndex]) {
      RecoENu = sr->common.ixn.pandora[0].Enu.e_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.heshw.X(),sr->common.ixn.pandora[0].dir.heshw.Y(),sr->common.ixn.pandora[0].dir.heshw.Z())).Unit();
    } else {
      RecoENu = sr->common.ixn.pandora[0].Enu.lep_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.lngtrk.X(),sr->common.ixn.pandora[0].dir.lngtrk.Y(),sr->common.ixn.pandora[0].dir.lngtrk.Z())).Unit();      
    }
    double RecoCZ = -RecoNuMomentumVector.Y(); // +Y in CAF files translates to +Z in typical CosZ
    if (std::isnan(RecoCZ)) {
      continue;
    }
    if (std::isnan(RecoENu)) {
      continue;
    }
    
    nSensibleEntries += 1;
  }

  MACH3LOG_INFO("Found {} sensible events out of the {} MC entries",nSensibleEntries,nChainEntries);
  
  dunemcSamples.resize(nSensibleEntries);

  //================================================================================================
  //Now load "sensible" events
  
  int iEvent = 0;
  for (int iChainEntry=0;iChainEntry<nChainEntries;iChainEntry++) {
    Chain->GetEntry(iChainEntry);

    if ((iChainEntry % (nChainEntries/10))==0) {
      MACH3LOG_INFO("\tProcessing entry: {}/{}",iChainEntry,nChainEntries);
    }
    
    if(sr->common.ixn.pandora.size() != 1) {
      MACH3LOG_WARN("Skipping entry {}/{} -> Number of neutrino slices found in event: {}",iChainEntry,nChainEntries,sr->common.ixn.pandora.size());
      continue;
    }

    const size_t SampleIndex = FileIndexToSample[static_cast<size_t>(Chain->GetTreeNumber())];
    
    TVector3 RecoNuMomentumVector;
    double RecoENu;
    if (IsELike[SampleIndex]) {
      RecoENu = sr->common.ixn.pandora[0].Enu.e_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.heshw.X(),sr->common.ixn.pandora[0].dir.heshw.Y(),sr->common.ixn.pandora[0].dir.heshw.Z())).Unit();
    } else {
      RecoENu = sr->common.ixn.pandora[0].Enu.lep_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.lngtrk.X(),sr->common.ixn.pandora[0].dir.lngtrk.Y(),sr->common.ixn.pandora[0].dir.lngtrk.Z())).Unit();      
    }
    double RecoCZ = -RecoNuMomentumVector.Y(); // +Y in CAF files translates to +Z in typical CosZ
    if (std::isnan(RecoCZ)) {
      MACH3LOG_WARN("Skipping entry {}/{} -> Reconstructed Cosine Z is NAN",iChainEntry,nChainEntries);
      continue;
    }
    if (std::isnan(RecoENu)) {
      MACH3LOG_WARN("Skipping entry {}/{} -> Reconstructed Neutrino Energy is NAN",iChainEntry,nChainEntries);
      continue;
    }

    if (iEvent >= nSensibleEntries) {
      MACH3LOG_ERROR("EventIndex ({}) greater than number of SensibleEntries found ({})",iEvent,nSensibleEntries);
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    dunemcSamples[iEvent].rw_erec = RecoENu;
    dunemcSamples[iEvent].rw_theta = RecoCZ;
    
    auto& OscillationChannels = SampleDetails[SampleIndex].OscChannels;
    std::string CurrFileName = Chain->GetCurrentFile()->GetName();
    dunemcSamples[iEvent].SampleIndex = static_cast<unsigned int>(SampleIndex);
    dunemcSamples[iEvent].nupdgUnosc = GetInitPDGFromFileName(CurrFileName);
    dunemcSamples[iEvent].nupdg = GetFinalPDGFromFileName(CurrFileName);
    dunemcSamples[iEvent].OscChannelIndex = static_cast<double>(GetOscChannel(OscillationChannels, dunemcSamples[iEvent].nupdgUnosc, dunemcSamples[iEvent].nupdg));
    
    int M3Mode = Modes->GetModeFromGenerator(std::abs(sr->mc.nu[0].mode));
    if (!sr->mc.nu[0].iscc) M3Mode += 14; //Account for no ability to distinguish CC/NC
    if (M3Mode > 15) M3Mode -= 1; //Account for no NCSingleKaon
    dunemcSamples[iEvent].mode = M3Mode;
    
    dunemcSamples[iEvent].rw_isCC = sr->mc.nu[0].iscc;
    dunemcSamples[iEvent].Target = kTarget_Ar;
    
    dunemcSamples[iEvent].rw_etru = static_cast<double>(sr->mc.nu[0].E);

    TVector3 TrueNuMomentumVector = (TVector3(sr->mc.nu[0].momentum.X(),sr->mc.nu[0].momentum.Y(),sr->mc.nu[0].momentum.Z())).Unit();
    dunemcSamples[iEvent].rw_truecz = -TrueNuMomentumVector.Y(); // +Y in CAF files translates to +Z in typical CosZ

    dunemcSamples[iEvent].flux_w = sr->mc.nu[0].genweight;
    
    iEvent += 1;
  }

  if (iEvent != nSensibleEntries) {
    MACH3LOG_ERROR("Mismatch in number of counted 'sensible' entries vs number of events loaded");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  //================================================================================================
  delete Chain;
  gErrorIgnoreLevel = CurrErrorLevel;

  return nSensibleEntries;
}

void SampleHandlerAtm::SetupFDMC() {
  for(int iEvent = 0 ;iEvent < int(GetNEvents()) ; ++iEvent) {
    MCSamples[iEvent].rw_etru = &(dunemcSamples[iEvent].rw_etru);
    MCSamples[iEvent].mode = &(dunemcSamples[iEvent].mode);
    MCSamples[iEvent].Target = &(dunemcSamples[iEvent].Target);    
    MCSamples[iEvent].isNC = !dunemcSamples[iEvent].rw_isCC;
    MCSamples[iEvent].nupdg = &(dunemcSamples[iEvent].nupdg);
    MCSamples[iEvent].nupdgUnosc = &(dunemcSamples[iEvent].nupdgUnosc);
    MCSamples[iEvent].NominalSample = dunemcSamples[iEvent].SampleIndex;
    
    MCSamples[iEvent].rw_truecz = &(dunemcSamples[iEvent].rw_truecz);
  }
}

const double* SampleHandlerAtm::GetPointerToKinematicParameter(KinematicTypes KinPar, int iEvent) {
  double* KinematicValue;

  switch (KinPar) {
  case kTrueNeutrinoEnergy:
    KinematicValue = &(dunemcSamples[iEvent].rw_etru);
    break;
  case kRecoNeutrinoEnergy:
    KinematicValue = &(dunemcSamples[iEvent].rw_erec);
    break;
  case kTrueCosZ:
    KinematicValue = &(dunemcSamples[iEvent].rw_truecz);
    break;
  case kRecoCosZ:
    KinematicValue = &(dunemcSamples[iEvent].rw_theta);
    break;
  case kOscChannel:
    KinematicValue = &(dunemcSamples[iEvent].OscChannelIndex);
    break;
  case kMode:
    KinematicValue = &(dunemcSamples[iEvent].mode);
    break;
  default:
    MACH3LOG_ERROR("Unknown KinPar: {}",static_cast<int>(KinPar));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  return KinematicValue;
}

const double* SampleHandlerAtm::GetPointerToKinematicParameter(double KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return GetPointerToKinematicParameter(KinPar,iEvent);
}

const double* SampleHandlerAtm::GetPointerToKinematicParameter(std::string KinematicParameter, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return GetPointerToKinematicParameter(KinPar,iEvent);
}

double SampleHandlerAtm::ReturnKinematicParameter(int KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iEvent);
}

double SampleHandlerAtm::ReturnKinematicParameter(std::string KinematicParameter, int iEvent) {
  return *GetPointerToKinematicParameter(KinematicParameter, iEvent);
}

std::vector<double> SampleHandlerAtm::ReturnKinematicParameterBinning(std::string KinematicParameterStr) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameterStr));
  return ReturnKinematicParameterBinning(KinPar);
}

std::vector<double> SampleHandlerAtm::ReturnKinematicParameterBinning(KinematicTypes KinPar)  {
  (void)KinPar;
  std::vector<double> ReturnVec;
  return ReturnVec;
}
