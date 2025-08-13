#include "SampleHandlerAtm.h"
#include "TSpline.h"
#include "TFile.h"
#include <map>

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wfloat-conversion"
//Standard Record includes
#include "duneanaobj/StandardRecord/StandardRecord.h"
#pragma GCC diagnostic pop

////////// Muyuan He: adding flux spline handler for atm systematics //////////

// Flux Spline Handler for Atmospheric Systematics, local use only, not declared in .h ---
// Store each flavor's energy and coszen splines
struct AtmFluxSplinePair {
  TSpline3* spline_E = nullptr;
  TSpline3* spline_Cos = nullptr;
};
// add a global std::map to associate flavor string to its spline pair
std::map<int, AtmFluxSplinePair> flux_splines;

// Function to load flux splines from file
TSpline3* LoadSpline(const std::string& filename, const std::string& splinename) {
  TFile* file = TFile::Open(filename.c_str(), "READ");
  if (!file || file->IsZombie()) { // Check if file opened successfully
    MACH3LOG_ERROR("Failed to open spline file: {}", filename);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  TSpline3* spline = dynamic_cast<TSpline3*>(file->Get(splinename.c_str())); // Retrieve the spline
  if (!spline) {
    MACH3LOG_ERROR("Failed to retrieve spline from file: {}", splinename, filename);
    file->Close();
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  // Clone the spline to avoid issues when closing the file
  TSpline3* cloned_spline = dynamic_cast<TSpline3*>(spline->Clone());
  file->Close();
  
  return cloned_spline;
}
// add a function to load both E and Cos splines for a given flavor
void LoadFluxSplinesForFlavor(const std::string& file, int pdg, const std::string& tag) {
  AtmFluxSplinePair sp; // create a new spline pair
  sp.spline_E   = LoadSpline(file, "atmflux_" + tag + "_E");
  sp.spline_Cos = LoadSpline(file, "atmflux_" + tag + "_Cos");
  if (sp.spline_E && sp.spline_Cos) {
    flux_splines[pdg] = sp; // store in the global map
    MACH3LOG_INFO("Loaded flux splines for flavor PDG {} with tag {}", pdg, tag);
  } else {
    MACH3LOG_ERROR("Failed to load flux splines for flavor PDG {} with tag {}", pdg, tag);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}
// add a master function to load all splines based on the config
void LoadAllFluxSplines(const YAML::Node& sample_node) {
  if (!sample_node["FluxSplineFile"] || !sample_node["FluxFlavors"]) { // Check if section exists
    MACH3LOG_INFO("[LoadAllFluxSplines] No FluxSplineFile or FluxFlavors section found in YAML. Skipping spline loading.");
    return;
  }

  std::string spline_file = Get<std::string>(sample_node["FluxSplineFile"], __FILE__, __LINE__);
  const YAML::Node& flavors = sample_node["FluxFlavors"];
  for (const auto& flavor_node : flavors) {
    int pdg = Get<int>(flavor_node["PDG"], __FILE__, __LINE__);
    std::string tag = Get<std::string>(flavor_node["Tag"], __FILE__, __LINE__);
    LoadFluxSplinesForFlavor(spline_file, pdg, tag);
  }
}








////////// End of flux spline handler //////////






SampleHandlerAtm::SampleHandlerAtm(std::string mc_version_, ParameterHandlerGeneric* xsec_cov_, const std::shared_ptr<OscillationHandler>&  Oscillator_) : SampleHandlerFD(mc_version_, xsec_cov_, Oscillator_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;

  LoadAllFluxSplines(SampleManager->raw());
  
  Initialise();
}

SampleHandlerAtm::~SampleHandlerAtm() {
}

void SampleHandlerAtm::Init() {
  dunemcSamples.resize(nSamples,dunemc_base());
  
  IsELike = Get<bool>(SampleManager->raw()["SampleOptions"]["IsELike"],__FILE__,__LINE__);
  ExposureScaling = Get<double>(SampleManager->raw()["SampleOptions"]["ExposureScaling"],__FILE__,__LINE__);
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

int SampleHandlerAtm::SetupExperimentMC() {
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
    
    // Choose direction source based on IsELike
    double x, y, z;
    if (IsELike) {
      x = sr->common.ixn.pandora[0].dir.heshw.X();
      y = sr->common.ixn.pandora[0].dir.heshw.Y();
      z = sr->common.ixn.pandora[0].dir.heshw.Z();
    } else {
      x = sr->common.ixn.pandora[0].dir.lngtrk.X();
      y = sr->common.ixn.pandora[0].dir.lngtrk.Y();
      z = sr->common.ixn.pandora[0].dir.lngtrk.Z();
    }

    TVector3 raw_dir(x, y, z);
    if (raw_dir.Mag2() < 1e-6) { // safeguard against zero vector
      MACH3LOG_WARN("Skipping event {} -> zero direction vector", iEvent);
      continue;
    }

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

std::vector<double> SampleHandlerAtm::ReturnKinematicParameterBinning(std::string KinematicParameterStr) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameterStr));
  return ReturnKinematicParameterBinning(KinPar);
}

std::vector<double> SampleHandlerAtm::ReturnKinematicParameterBinning(KinematicTypes KinPar)  {
  (void)KinPar;
  std::vector<double> ReturnVec;
  return ReturnVec;
}
