#include "SampleHandlerAtm.h"

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wfloat-conversion"
//Standard Record includes
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#pragma GCC diagnostic pop

SampleHandlerAtm::SampleHandlerAtm(std::string mc_version_, ParameterHandlerGeneric* xsec_cov_, const std::shared_ptr<OscillationHandler>&  Oscillator_) : SampleHandlerFD(mc_version_, xsec_cov_, Oscillator_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  
  Initialise();
}

SampleHandlerAtm::~SampleHandlerAtm() {
}

void SampleHandlerAtm::Init() {
  // dunemcSamples.resize(nSamples,dunemc_base());
  
  IsELike = Get<bool>(SampleManager->raw()["SampleOptions"]["IsELike"],__FILE__,__LINE__);
  ExposureScaling = Get<double>(SampleManager->raw()["SampleOptions"]["ExposureScaling"],__FILE__,__LINE__);
  fInputFile = Get<std::string>(SampleManager->raw()["SampleOptions"]["InputFile"],__FILE__,__LINE__);

  if (SampleManager->raw()["SampleOptions"]["InputSplines"]) {
    fInputSplines = Get<std::string>(SampleManager->raw()["SampleOptions"]["InputSplines"],__FILE__,__LINE__);
  } else {
    fInputSplines = "";
  }
  fSampleId = Get<uint>(SampleManager->raw()["SampleOptions"]["SampleId"],__FILE__,__LINE__);
}

void SampleHandlerAtm::SetupSplines() {
  ///@todo move all of the spline setup into core
  if(ParHandler->GetNumParamsFromSampleName(SampleName, kSpline) > 0){
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", ParHandler->GetNumParamsFromSampleName(SampleName, kSpline));
    auto SplineFactory = SplineHandlerFactoryDUNE(ParHandler,Modes.get(), dunemcSamples, fInputSplines, SampleName);

    SplineHandler = std::move(SplineFactory.GetSplineHandler());
    if (SplineFactory.GetSplineType() == kBinned){
      InitialiseSplineObject(); //Running the "normal" initialisation for binned splines
    }
  }
  else{
    MACH3LOG_INFO("Found no spline for this sample so I will not load or evaluate splines");
    SplineHandler = nullptr;
  }
  
  return;
}

void SampleHandlerAtm::SetupWeightPointers() {
  for (size_t i = 0; i < dunemcSamples.size(); ++i) {
    MCSamples[i].total_weight_pointers.push_back(&(dunemcSamples[i].flux_w));
    MCSamples[i].total_weight_pointers.push_back(MCSamples[i].osc_w_pointer);
    MCSamples[i].total_weight_pointers.push_back(&(MCSamples[i].xsec_w));
    MCSamples[i].total_weight_pointers.push_back(&(ExposureScaling));
  }
  
}

int SampleHandlerAtm::SetupExperimentMC() {
  int CurrErrorLevel = gErrorIgnoreLevel; 
  gErrorIgnoreLevel = kFatal;
  
  TFile *f = TFile::Open(fInputFile.c_str(),"READ");
  if (!f || f->IsZombie()) {
    MACH3LOG_ERROR("Could not open input CAF file: {}",fInputFile);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  TTree *cafTree, *weightsTree;
  f->GetObject("cafTree",cafTree);
  if (!cafTree) {
    MACH3LOG_ERROR("Could not find cafTree in input CAF file: {}",fInputFile);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  f->GetObject("weights",weightsTree);
  if (!weightsTree) {
    MACH3LOG_ERROR("Could not find weights tree in input CAF file: {}",fInputFile);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  uint sample_id;
  double xsec_w, flux_nue_w, flux_numu_w;
  weightsTree->SetBranchAddress("sample_id",&sample_id);
  weightsTree->SetBranchAddress("xsec",&xsec_w);
  weightsTree->SetBranchAddress("flux_nue",&flux_nue_w);
  weightsTree->SetBranchAddress("flux_numu",&flux_numu_w);


  caf::StandardRecordProxy* sr = new caf::StandardRecordProxy(cafTree, "rec");

  int nEntries = static_cast<int>(cafTree->GetEntries());

  //Define Oscillation Channels here for the first file to avoid multiple lookups
 
  for (int iEvent=0;iEvent<nEntries;iEvent++) {

    weightsTree->GetEntry(iEvent);
    if (sample_id != fSampleId) continue;
    
    cafTree->LoadTree(iEvent); //Only loads the tree without reading the entire entry

    if ((iEvent % (nEntries/10))==0) {
      MACH3LOG_INFO("\tProcessing event: {}/{}",iEvent,nEntries);
    }
    
    if(sr->common.ixn.pandora.size() != 1) {
      // MACH3LOG_WARN("Skipping event {}/{} -> Number of neutrino slices found in event: {}",iEvent,nEntries,sr->common.ixn.pandora.size());
      continue;
    }
    
    TVector3 RecoNuMomentumVector;
    double RecoENu;
    if (IsELike) {
      RecoENu = sr->common.ixn.pandora[0].Enu.e_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.heshw.x,sr->common.ixn.pandora[0].dir.heshw.y,sr->common.ixn.pandora[0].dir.heshw.z)).Unit();
    } else {
      RecoENu = sr->common.ixn.pandora[0].Enu.lep_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.lngtrk.x,sr->common.ixn.pandora[0].dir.lngtrk.y,sr->common.ixn.pandora[0].dir.lngtrk.z)).Unit();      
    }
    double RecoCZ = -RecoNuMomentumVector.Y(); // +Y in CAF files translates to +Z in typical CosZ
    if (std::isnan(RecoCZ)) {
      // MACH3LOG_WARN("Skipping event {}/{} -> Reconstructed Cosine Z is NAN",iEvent,nEntries);
      continue;
    }
    if (std::isnan(RecoENu)) {
      // MACH3LOG_WARN("Skipping event {}/{} -> Reconstructed Neutrino Energy is NAN",iEvent,nEntries);
      continue;
    }

    struct dunemc_base currentEventFromNuMu;

    currentEventFromNuMu.rw_erec = RecoENu;
    currentEventFromNuMu.rw_theta = RecoCZ;
    
    int M3Mode = Modes->GetModeFromGenerator(std::abs(sr->mc.nu[0].mode));
    if (!sr->mc.nu[0].iscc) M3Mode += 14; //Account for no ability to distinguish CC/NC
    if (M3Mode > 15) M3Mode -= 1; //Account for no NCSingleKaon
    currentEventFromNuMu.mode = M3Mode;
    
    currentEventFromNuMu.rw_isCC = sr->mc.nu[0].iscc;
    currentEventFromNuMu.Target = kTarget_Ar;
    
    currentEventFromNuMu.rw_etru = static_cast<double>(sr->mc.nu[0].E);

    TVector3 TrueNuMomentumVector = (TVector3(sr->mc.nu[0].momentum.x,sr->mc.nu[0].momentum.y,sr->mc.nu[0].momentum.z)).Unit();
    currentEventFromNuMu.rw_truecz = -TrueNuMomentumVector.Y(); // +Y in CAF files translates to +Z in typical CosZ

    currentEventFromNuMu.nupdg = sr->mc.nu[0].pdg;

    currentEventFromNuMu.flux_w = xsec_w*flux_numu_w;
    currentEventFromNuMu.nupdgUnosc = (currentEventFromNuMu.nupdg > 0) ? 14 : -14;
    currentEventFromNuMu.OscChannelIndex = static_cast<double>(GetOscChannel(OscChannels, currentEventFromNuMu.nupdgUnosc, currentEventFromNuMu.nupdg));
    currentEventFromNuMu.eid = iEvent;

    struct dunemc_base currentEventFromNuE = currentEventFromNuMu;;
    currentEventFromNuE.nupdgUnosc = (currentEventFromNuE.nupdg > 0) ? 12 : -12;
    currentEventFromNuE.flux_w = xsec_w*flux_nue_w;
    currentEventFromNuE.OscChannelIndex = static_cast<double>(GetOscChannel(OscChannels, currentEventFromNuE.nupdgUnosc, currentEventFromNuE.nupdg));
    //Debug flux printout
    // std::cout << "Event " << iEvent 
    //           << ": nupdg = " << currentEventFromNuE.nupdg 
    //           << ", nupdgUnosc = " << currentEventFromNuE.nupdgUnosc 
    //           << ", flux_nue_w = " << flux_nue_w 
    //           << ", flux_numu_w = " << flux_numu_w 
    //           << ", xsec_w = " << xsec_w 
    //           << ", total flux weight nue = " << currentEventFromNuE.flux_w
    //           << ", total flux weight numu = " << currentEventFromNuMu.flux_w
    //           << std::endl;
    dunemcSamples.emplace_back(std::move(currentEventFromNuMu));
    dunemcSamples.emplace_back(std::move(currentEventFromNuE));
  }

  delete cafTree;
  delete weightsTree;
  delete f;
  delete sr;
  gErrorIgnoreLevel = CurrErrorLevel;

  return dunemcSamples.size();
}

void SampleHandlerAtm::SetupFDMC() {
  for(int iEvent = 0 ;iEvent < int(GetNEvents()) ; ++iEvent) {
    MCSamples[iEvent].rw_etru = &(dunemcSamples[iEvent].rw_etru);
    MCSamples[iEvent].mode = &(dunemcSamples[iEvent].mode);
    MCSamples[iEvent].Target = &(dunemcSamples[iEvent].Target);    
    MCSamples[iEvent].isNC = !dunemcSamples[iEvent].rw_isCC;
    MCSamples[iEvent].nupdg = &(dunemcSamples[iEvent].nupdg);
    MCSamples[iEvent].nupdgUnosc = &(dunemcSamples[iEvent].nupdgUnosc);

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
