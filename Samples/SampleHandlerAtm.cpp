#include "SampleHandlerAtm.h"

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wfloat-conversion"
//Standard Record includes
#include "duneanaobj/StandardRecord/StandardRecord.h"
#pragma GCC diagnostic pop

#include <Eigen/Dense>
#include "DUNEUtils.h"

SampleHandlerAtm::SampleHandlerAtm(std::string mc_version_, ParameterHandlerGeneric* xsec_cov_, const std::shared_ptr<OscillationHandler>&  Oscillator_) : SampleHandlerFD(mc_version_, xsec_cov_, Oscillator_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  
  Initialise();
}

SampleHandlerAtm::~SampleHandlerAtm() {
}

void SampleHandlerAtm::Init() {
  IsELike = Get<bool>(SampleManager->raw()["SampleOptions"]["IsELike"],__FILE__,__LINE__);
  ExposureScaling = Get<double>(SampleManager->raw()["SampleOptions"]["ExposureScaling"],__FILE__,__LINE__);

  if (SampleManager->raw()["SampleOptions"]["EigenFile"]) {
    EigenInputFile = Get<std::string>(SampleManager->raw()["SampleOptions"]["EigenFile"],__FILE__,__LINE__);
    EigenInputFileMD5Sum = Get<std::string>(SampleManager->raw()["SampleOptions"]["EigenFileMD5Sum"],__FILE__,__LINE__);
  } else {
    EigenInputFile = "";
  }
}

void SampleHandlerAtm::SetupSplines() {
  SplineHandler = nullptr;
}

void SampleHandlerAtm::SetupWeightPointers() {
  for (size_t i = 0; i < dunemcSamples.size(); ++i) {
    MCSamples[i].total_weight_pointers.push_back(&(dunemcSamples[i].flux_w));
    MCSamples[i].total_weight_pointers.push_back(MCSamples[i].osc_w_pointer);
    MCSamples[i].total_weight_pointers.push_back(&(MCSamples[i].xsec_w));
    MCSamples[i].total_weight_pointers.push_back(&(ExposureScaling));
  }  
}

int SampleHandlerAtm::ReadFromEigen() {
  Eigen::MatrixXd Matrix;

  MACH3LOG_INFO("Reading Eigen matrix from:{}",EigenInputFile);
  CompareMD5Sum(EigenInputFileMD5Sum,EigenInputFile);
  ReadEigenMatrixFromFile(EigenInputFile.c_str(),Matrix);

  int nEntries = static_cast<int>(Matrix.rows());
  dunemcSamples.resize(nEntries);

  for (size_t iEvent=0;iEvent<dunemcSamples.size();++iEvent) {
    dunemcSamples[iEvent].Target = static_cast<int>(Matrix(iEvent,Variables::Target));
    dunemcSamples[iEvent].nupdg = static_cast<int>(Matrix(iEvent,Variables::nupdg));
    dunemcSamples[iEvent].nupdgUnosc = static_cast<int>(Matrix(iEvent,Variables::nupdgUnosc));
    dunemcSamples[iEvent].rw_isCC = static_cast<int>(Matrix(iEvent,Variables::rw_isCC));
    dunemcSamples[iEvent].OscChannelIndex = Matrix(iEvent,Variables::OscChannelIndex);
    dunemcSamples[iEvent].rw_erec = Matrix(iEvent,Variables::rw_erec);
    dunemcSamples[iEvent].rw_etru = Matrix(iEvent,Variables::rw_etru);
    dunemcSamples[iEvent].flux_w = Matrix(iEvent,Variables::flux_w);
    dunemcSamples[iEvent].mode = Matrix(iEvent,Variables::mode);
    dunemcSamples[iEvent].rw_theta = Matrix(iEvent,Variables::rw_theta);
    dunemcSamples[iEvent].rw_truecz = Matrix(iEvent,Variables::rw_truecz);
  }

  return nEntries;
}

void SampleHandlerAtm::TransferToEigen(std::string FileName) {
  Eigen::MatrixXd Matrix = Eigen::MatrixXd(dunemcSamples.size(),Variables::nVariables);

  for (size_t iEvent=0;iEvent<dunemcSamples.size();++iEvent) {
    Matrix(iEvent,Variables::Target) = dunemcSamples[iEvent].Target;
    Matrix(iEvent,Variables::nupdg) = dunemcSamples[iEvent].nupdg;
    Matrix(iEvent,Variables::nupdgUnosc) = dunemcSamples[iEvent].nupdgUnosc;
    Matrix(iEvent,Variables::rw_isCC) = dunemcSamples[iEvent].rw_isCC;
    Matrix(iEvent,Variables::OscChannelIndex) = dunemcSamples[iEvent].OscChannelIndex;
    Matrix(iEvent,Variables::rw_erec) = dunemcSamples[iEvent].rw_erec;
    Matrix(iEvent,Variables::rw_etru) = dunemcSamples[iEvent].rw_etru;
    Matrix(iEvent,Variables::flux_w) = dunemcSamples[iEvent].flux_w;
    Matrix(iEvent,Variables::mode) = dunemcSamples[iEvent].mode;
    Matrix(iEvent,Variables::rw_theta) = dunemcSamples[iEvent].rw_theta;
    Matrix(iEvent,Variables::rw_truecz) = dunemcSamples[iEvent].rw_truecz;
  }

  MACH3LOG_INFO("Writing Eigen matrix to:{}",FileName);
  WriteEigenMatrixToFile(FileName.c_str(),Matrix);
}

int SampleHandlerAtm::SetupExperimentMC() {
  // If the Eigen input file is define, read from that. Otherwise read from CAF file
  if (EigenInputFile != "") {
    int nEntries = ReadFromEigen();
    return nEntries;
  }
  
  int CurrErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  
  caf::StandardRecord* sr = new caf::StandardRecord();

  TChain* Chain = new TChain("cafTree");
  for (size_t iSample=0;iSample<mc_files.size();iSample++) {
    Chain->Add(mc_files[iSample].c_str());
  }
  
  Chain->SetBranchStatus("*", 1);
  Chain->SetBranchAddress("rec", &sr);

  int nEntries = static_cast<int>(Chain->GetEntries());
  dunemcSamples.resize(nEntries);
 
  for (int iEvent=0;iEvent<nEntries;iEvent++) {
    Chain->GetEntry(iEvent);

    if ((iEvent % (nEntries/10))==0) {
      MACH3LOG_INFO("\tProcessing event: {}/{}",iEvent,nEntries);
    }
    
    if(sr->common.ixn.pandora.size() != 1) {
      MACH3LOG_WARN("Skipping event {}/{} -> Number of neutrino slices found in event: {}",iEvent,nEntries,sr->common.ixn.pandora.size());
      continue;
    }
    
    TVector3 RecoNuMomentumVector;
    double RecoENu;
    if (IsELike) {
      RecoENu = sr->common.ixn.pandora[0].Enu.e_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.heshw.X(),sr->common.ixn.pandora[0].dir.heshw.Y(),sr->common.ixn.pandora[0].dir.heshw.Z())).Unit();
    } else {
      RecoENu = sr->common.ixn.pandora[0].Enu.lep_calo;
      RecoNuMomentumVector = (TVector3(sr->common.ixn.pandora[0].dir.lngtrk.X(),sr->common.ixn.pandora[0].dir.lngtrk.Y(),sr->common.ixn.pandora[0].dir.lngtrk.Z())).Unit();      
    }
    double RecoCZ = -RecoNuMomentumVector.Y(); // +Y in CAF files translates to +Z in typical CosZ
    if (std::isnan(RecoCZ)) {
      MACH3LOG_WARN("Skipping event {}/{} -> Reconstructed Cosine Z is NAN",iEvent,nEntries);
      continue;
    }
    if (std::isnan(RecoENu)) {
      MACH3LOG_WARN("Skipping event {}/{} -> Reconstructed Neutrino Energy is NAN",iEvent,nEntries);
      continue;
    }
    dunemcSamples[iEvent].rw_erec = RecoENu;
    dunemcSamples[iEvent].rw_theta = RecoCZ;
    
    std::string CurrFileName = Chain->GetCurrentFile()->GetName();
    dunemcSamples[iEvent].nupdgUnosc = GetInitPDGFromFileName(CurrFileName);
    dunemcSamples[iEvent].nupdg = GetFinalPDGFromFileName(CurrFileName);
    dunemcSamples[iEvent].OscChannelIndex = static_cast<double>(GetOscChannel(OscChannels, dunemcSamples[iEvent].nupdgUnosc, dunemcSamples[iEvent].nupdg));
    
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
  }

  delete Chain;
  gErrorIgnoreLevel = CurrErrorLevel;

  return nEntries;
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
