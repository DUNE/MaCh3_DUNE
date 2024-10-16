#include "samplePDFDUNEAtm.h"

//Standard Record includes
#include "StandardRecord.h"

//ROOT includes
#include "TError.h"

samplePDFDUNEAtm::samplePDFDUNEAtm(std::string mc_version_, covarianceXsec* xsec_cov_) : samplePDFFDBase(mc_version_, xsec_cov_) {
  Initialise();
}

samplePDFDUNEAtm::~samplePDFDUNEAtm() {
}

void samplePDFDUNEAtm::Init() {
  for (int i=0;i<nSamples;i++) {
    struct dunemc_base obj = dunemc_base();
    dunemcSamples.push_back(obj);
  }
  
  IsELike = SampleManager->raw()["SampleBools"]["IsELike"].as<bool>();
}

void samplePDFDUNEAtm::SetupSplines() {
  splineFile = nullptr;
}

void samplePDFDUNEAtm::SetupWeightPointers() {
  for (int i = 0; i < (int)dunemcSamples.size(); ++i) {
    for (int j = 0; j < dunemcSamples[i].nEvents; ++j) {
      MCSamples[i].ntotal_weight_pointers[j] = 3;
      MCSamples[i].total_weight_pointers[j] = new const double*[MCSamples[i].ntotal_weight_pointers[j]];
      MCSamples[i].total_weight_pointers[j][0] = MCSamples[i].osc_w_pointer[j];
      MCSamples[i].total_weight_pointers[j][1] = &(dunemcSamples[i].flux_w[j]);
      MCSamples[i].total_weight_pointers[j][2] = &(MCSamples[i].xsec_w[j]);
    }
  }
}

int samplePDFDUNEAtm::setupExperimentMC(int iSample) {
  int CurrErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  
  caf::StandardRecord* sr = new caf::StandardRecord();
  dunemc_base* duneobj = &dunemcSamples[iSample];

  std::string FileName = mtuple_files[iSample];
  std::cout << "Reading File:" << FileName << std::endl;
  TFile* File = TFile::Open(FileName.c_str());
  if (!File || File->IsZombie()) {
    std::cerr << "Did not find File:" << FileName << std::endl;
    throw;
  }
  TTree* Tree = File->Get<TTree>("cafTree");
  if (!Tree){
    std::cerr << "Did not find Tree::cafTree in File:" << FileName << std::endl;
    throw;
  }
  
  Tree->SetBranchStatus("*", 1);
  Tree->SetBranchAddress("rec", &sr);

  duneobj->nEvents = Tree->GetEntries();
  duneobj->norm_s = 1;
  duneobj->pot_s = 1;

  duneobj->nutype = sample_nutype[iSample];
  duneobj->oscnutype = sample_oscnutype[iSample];
  duneobj->signal = sample_signal[iSample];
  
  duneobj->mode.resize(duneobj->nEvents);
  duneobj->rw_isCC.resize(duneobj->nEvents);
  duneobj->Target.resize(duneobj->nEvents);
  
  duneobj->rw_etru.resize(duneobj->nEvents);
  duneobj->rw_truecz.resize(duneobj->nEvents);
  duneobj->flux_w.resize(duneobj->nEvents);
  duneobj->rw_erec.resize(duneobj->nEvents);
  duneobj->rw_theta.resize(duneobj->nEvents);
 
  for (int iEvent=0;iEvent<duneobj->nEvents;iEvent++) {
    Tree->GetEntry(iEvent);

    if ((iEvent % (duneobj->nEvents/10))==0) {
      std::cout << "\tProcessing event: " << iEvent << "/" << duneobj->nEvents << std::endl;
    }

    duneobj->mode[iEvent] = SIMBMode_ToMaCh3Mode(sr->mc.nu[0].mode,sr->mc.nu[0].iscc);
    duneobj->rw_isCC[iEvent] = sr->mc.nu[0].iscc;
    duneobj->Target[iEvent] = 40;
    
    duneobj->rw_etru[iEvent] = (double)(sr->mc.nu[0].E);

    TVector3 TrueNuMomentumVector = (TVector3(sr->mc.nu[0].momentum.X(),sr->mc.nu[0].momentum.Y(),sr->mc.nu[0].momentum.Z())).Unit();
    duneobj->rw_truecz[iEvent] = TrueNuMomentumVector.Y();

    duneobj->flux_w[iEvent] = sr->mc.nu[0].genweight;

    TVector3 RecoNuMomentumVector;
    if (IsELike) {
      duneobj->rw_erec[iEvent] = sr->common.ixn.pandora[0].Enu.e_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.heshw.X(),sr->common.ixn.pandora[0].dir.heshw.Y(),sr->common.ixn.pandora[0].dir.heshw.Z())).Unit();
    } else {
      duneobj->rw_erec[iEvent] = sr->common.ixn.pandora[0].Enu.lep_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.lngtrk.X(),sr->common.ixn.pandora[0].dir.lngtrk.Y(),sr->common.ixn.pandora[0].dir.lngtrk.Z())).Unit();      
    }
    duneobj->rw_theta[iEvent] = RecoNuMomentumVector.Y();
    
  }

  delete Tree;
  delete File;

  gErrorIgnoreLevel = CurrErrorLevel;
  
  return duneobj->nEvents;
}

void samplePDFDUNEAtm::setupFDMC(int iSample) {
  dunemc_base *duneobj = &(dunemcSamples[iSample]);
  fdmc_base *fdobj = &(MCSamples[iSample]);  

  //a bit of a hack, make sure that this is only set if you're an atmoshperic object
  //i.e. this doesn't live in samplePDFFDBase
  fdobj->rw_truecz = new const double*[fdobj->nEvents];
  
  fdobj->nutype = duneobj->nutype;
  fdobj->oscnutype = duneobj->oscnutype;
  fdobj->signal = duneobj->signal;
  fdobj->SampleDetID = SampleDetID;
  
  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent) {
    fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
    fdobj->isNC[iEvent] = !duneobj->rw_isCC[iEvent];
    fdobj->Target[iEvent] = &(duneobj->Target[iEvent]);
    
    fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
    fdobj->rw_truecz[iEvent] = &(duneobj->rw_truecz[iEvent]);
  }
}

double const &samplePDFDUNEAtm::ReturnKinematicParameterByReference(int KinematicParameter, int iSample, int iEvent) {

  switch (KinematicParameter) {
  case kTrueNeutrinoEnergy:
    return (dunemcSamples[iSample].rw_etru[iEvent]);
  case kRecoNeutrinoEnergy:
    return (dunemcSamples[iSample].rw_erec[iEvent]);
  case kTrueCosZ:
    return (dunemcSamples[iSample].rw_truecz[iEvent]);
  case kRecoCosZ:
    return (dunemcSamples[iSample].rw_theta[iEvent]);
  default:
    std::cerr << "Unknown KinematicParameter:" << KinematicParameter << std::endl;
    throw;
  }
}

double samplePDFDUNEAtm::ReturnKinematicParameter(int KinematicParameter, int iSample, int iEvent) {
  return ReturnKinematicParameterByReference(KinematicParameter, iSample, iEvent);
}

std::vector<double> samplePDFDUNEAtm::ReturnKinematicParameterBinning(int KinematicParameterStr)  {
}

int samplePDFDUNEAtm::ReturnKinematicParameterFromString(std::string KinematicParameterStr) {
  int ReturnVal;

  if (KinematicParameterStr == "TrueNeutrinoEnergy") {ReturnVal = kTrueNeutrinoEnergy;}
  else if (KinematicParameterStr == "RecoNeutrinoEnergy") {ReturnVal = kRecoNeutrinoEnergy;}
  else if (KinematicParameterStr == "TrueCosineZ") {ReturnVal = kTrueCosZ;}
  else if (KinematicParameterStr == "RecoCosineZ") {ReturnVal = kRecoCosZ;}
  else {
    std::cerr << "KinematicParameterStr: " << KinematicParameterStr << " not found" << std::endl;
    throw;
  }
  
  return ReturnVal;
}

std::string samplePDFDUNEAtm::ReturnStringFromKinematicParameter(int KinematicParameter) {
  return "";
}
