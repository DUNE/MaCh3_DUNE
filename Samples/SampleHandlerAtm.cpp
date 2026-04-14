#include "SampleHandlerAtm.h"

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wfloat-conversion"
//Standard Record includes
#include "duneanaobj/StandardRecord/StandardRecord.h"
#pragma GCC diagnostic pop

SampleHandlerAtm::SampleHandlerAtm(std::string mc_version_, ParameterHandlerGeneric* xsec_cov_, const std::shared_ptr<OscillationHandler>&  Oscillator_) : SampleHandlerBase(mc_version_, xsec_cov_, Oscillator_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  
  Initialise();
}

SampleHandlerAtm::~SampleHandlerAtm() {
}

void SampleHandlerAtm::Init() {
  std::vector<std::string> EnabledSamples = Get<std::vector<std::string>>(SampleManager->raw()["Samples"], __FILE__ , __LINE__);
  IsELike.resize(GetNSamples());
  
  ExposureScaling = Get<double>(SampleManager->raw()["AnalysisOptions"]["ExposureScaling"],__FILE__,__LINE__);
  for(int iSample=0;iSample<GetNSamples();iSample++){
    const std::string TempTitle = EnabledSamples[iSample];
    IsELike[iSample] = Get<int>(SampleManager->raw()[TempTitle]["SampleOptions"]["IsELike"],__FILE__,__LINE__);
  }
}

// ************************************************
void SampleHandlerAtm::InititialiseData()
{
  // ************************************************
  // Reweight MC to match
  Reweight();
  // set asimov data
  for (int iSample = 0; iSample < GetNSamples(); iSample++)
  {
    AddData(iSample, GetMCArray(iSample));
  }
}

void SampleHandlerAtm::SetupSplines() {
  SplineHandler = nullptr;
}

void SampleHandlerAtm::AddAdditionalWeightPointers() {
  for (size_t i = 0; i < dunemcSamples.size(); ++i) {
    MCEvents[i].total_weight_pointers.push_back(&(dunemcSamples[i].flux_w));
    MCEvents[i].total_weight_pointers.push_back(&(ExposureScaling));
  }  
}

int SampleHandlerAtm::SetupExperimentMC() {
  int CurrErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  
  caf::StandardRecord* sr = new caf::StandardRecord();

  TChain* Chain = new TChain("cafTree");
  std::vector<size_t> FileIndexToSample;
  for (size_t iSample=0;iSample<SampleDetails.size();iSample++) {
    for (const std::vector<std::string>& files : SampleDetails[iSample].mc_files) {
      for (const std::string& filename : files){
        FileIndexToSample.push_back(iSample);
        int ChainAddCheck = Chain->Add(filename.c_str(), -1);
        if(ChainAddCheck == 0){
          MACH3LOG_ERROR("Could not add file {} to TChain, please check the file exists and is readable", filename);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
      }
    }
  }
  int nChainEntries = static_cast<int>(Chain->GetEntries());
  
  Chain->SetBranchStatus("*", 1);
  Chain->SetBranchAddress("rec", &sr);

  //================================================================================================
  
  //First set struct length to maximum number of events from the MC
  dunemcSamples.resize(nChainEntries);

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
    
    dunemcSamples[iEvent].enu_true = static_cast<double>(sr->mc.nu[0].E);

    TVector3 TrueNuMomentumVector = (TVector3(sr->mc.nu[0].momentum.X(),sr->mc.nu[0].momentum.Y(),sr->mc.nu[0].momentum.Z())).Unit();
    dunemcSamples[iEvent].coszenith_true = -TrueNuMomentumVector.Y(); // +Y in CAF files translates to +Z in typical CosZ

    dunemcSamples[iEvent].flux_w = sr->mc.nu[0].genweight;
    
    iEvent += 1;
  }

  //Now resize to number of events which passed "sensible" cuts
  dunemcSamples.resize(iEvent);

  //================================================================================================
  delete Chain;
  gErrorIgnoreLevel = CurrErrorLevel;

  return iEvent;
}

void SampleHandlerAtm::SetupMC() {
  for(int iEvent = 0 ;iEvent < int(GetNEvents()) ; ++iEvent) {
    MCEvents[iEvent].enu_true = dunemcSamples[iEvent].enu_true;
    MCEvents[iEvent].isNC = !dunemcSamples[iEvent].rw_isCC;
    MCEvents[iEvent].nupdg = dunemcSamples[iEvent].nupdg;
    MCEvents[iEvent].nupdgUnosc = dunemcSamples[iEvent].nupdgUnosc;
    MCEvents[iEvent].NominalSample = dunemcSamples[iEvent].SampleIndex;
    
    MCEvents[iEvent].coszenith_true = dunemcSamples[iEvent].coszenith_true;
  }
}

const double* SampleHandlerAtm::GetPointerToKinematicParameter(KinematicTypes KinPar, int iEvent) {
  double* KinematicValue;

  switch (KinPar) {
  case kTrueNeutrinoEnergy:
    KinematicValue = &(dunemcSamples[iEvent].enu_true);
    break;
  case kRecoNeutrinoEnergy:
    KinematicValue = &(dunemcSamples[iEvent].rw_erec);
    break;
  case kTrueCosZ:
    KinematicValue = &(dunemcSamples[iEvent].coszenith_true);
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

const double* SampleHandlerAtm::GetPointerToKinematicParameter(const int KinematicVariable, const int iEvent) const {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return GetPointerToKinematicParameter(KinPar,iEvent);
}

double SampleHandlerAtm::ReturnKinematicParameter(const int KinematicVariable, const int iEvent) const {
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iEvent);
}
