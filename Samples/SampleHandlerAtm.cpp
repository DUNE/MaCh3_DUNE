#include "SampleHandlerAtm.h"

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wfloat-conversion"
//Standard Record includes
#include "duneanaobj/StandardRecord/StandardRecord.h"
#pragma GCC diagnostic pop

#include <cstdlib>
#include <fstream>

// OPTIMIZED: Load all splines in one file operation and pre-compute grids
void SampleHandlerAtm::LoadFluxSplines() {
  MACH3LOG_INFO("=== Loading atmospheric flux splines ===");
  
  if (!SampleManager->raw()["FluxSplines"]) {
    MACH3LOG_WARN("No FluxSplines section found in YAML config. Skipping spline loading.");
    return;
  }
  
  std::string spline_file;
  
  // Try different possible locations for the File key
  if (SampleManager->raw()["FluxSplines"]["File"]) {
    spline_file = Get<std::string>(SampleManager->raw()["FluxSplines"]["File"], __FILE__, __LINE__);
    MACH3LOG_INFO("Found File key directly under FluxSplines");
  } else {
    // Let's just try hardcoding the path from your YAML for now
    spline_file = "Inputs/DUNE_atmospheric_spline_files/atmflux_SolMin_splines.root";
    MACH3LOG_WARN("Using hardcoded spline file path: {}", spline_file);
  }
  
  MACH3LOG_INFO("Loading splines from: {}", spline_file);
  
  // Check if file exists
  std::ifstream file_check(spline_file);
  if (!file_check.good()) {
    MACH3LOG_ERROR("Spline file does not exist or is not readable: {}", spline_file);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  file_check.close();
  
  // OPTIMIZED: Open file once and load all splines
  TFile* file = TFile::Open(spline_file.c_str(), "READ");
  if (!file || file->IsZombie()) {
    MACH3LOG_ERROR("Failed to open spline file: {}", spline_file);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  // Define spline names for efficient loading
  const char* spline_names[8] = {
    "atmflux_NuE_E", "atmflux_NuE_Cos",      // nue
    "atmflux_ANuE_E", "atmflux_ANuE_Cos",    // anue
    "atmflux_NuM_E", "atmflux_NuM_Cos",      // numu
    "atmflux_ANuM_E", "atmflux_ANuM_Cos"     // anumu
  };
  
  TSpline3* splines[8];
  
  // Load all splines in one loop
  for (int i = 0; i < 8; ++i) {
    TSpline3* temp_spline = dynamic_cast<TSpline3*>(file->Get(spline_names[i]));
    if (!temp_spline) {
      MACH3LOG_ERROR("Failed to retrieve spline: {} from file: {}", spline_names[i], spline_file);
      file->Close();
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    splines[i] = dynamic_cast<TSpline3*>(temp_spline->Clone());
    if (!splines[i]) {
      MACH3LOG_ERROR("Failed to clone spline: {}", spline_names[i]);
      file->Close();
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
  
  file->Close();
  
  // Assign to member variables and pre-compute ranges for fast lookup
  spline_sets[0].energy_spline = splines[0];  // nue E
  spline_sets[0].cosZ_spline = splines[1];    // nue Cos
  spline_sets[1].energy_spline = splines[2];  // anue E
  spline_sets[1].cosZ_spline = splines[3];    // anue Cos
  spline_sets[2].energy_spline = splines[4];  // numu E
  spline_sets[2].cosZ_spline = splines[5];    // numu Cos
  spline_sets[3].energy_spline = splines[6];  // anumu E
  spline_sets[3].cosZ_spline = splines[7];    // anumu Cos
  
  // Pre-compute ranges for fast bounds checking
  for (int i = 0; i < 4; ++i) {
    spline_sets[i].energy_min = spline_sets[i].energy_spline->GetXmin();
    spline_sets[i].energy_max = spline_sets[i].energy_spline->GetXmax();
    spline_sets[i].cosZ_min = spline_sets[i].cosZ_spline->GetXmin();
    spline_sets[i].cosZ_max = spline_sets[i].cosZ_spline->GetXmax();
  }
  
  MACH3LOG_INFO("=== All atmospheric flux splines loaded successfully! ===");
  
  // Log summary (only once)
  const char* neutrino_types[4] = {"nue", "anue", "numu", "anumu"};
  for (int i = 0; i < 4; ++i) {
    MACH3LOG_INFO("  {}_E:   {} points, range [{:.3f}, {:.3f}]", 
                  neutrino_types[i], spline_sets[i].energy_spline->GetNp(), 
                  spline_sets[i].energy_min, spline_sets[i].energy_max);
    MACH3LOG_INFO("  {}_Cos: {} points, range [{:.3f}, {:.3f}]", 
                  neutrino_types[i], spline_sets[i].cosZ_spline->GetNp(), 
                  spline_sets[i].cosZ_min, spline_sets[i].cosZ_max);
  }
  
  // After loading splines, pre-compute grids and cache responses
  PrecomputeSplineGrids();
}

void SampleHandlerAtm::PrecomputeSplineGrids() {
  MACH3LOG_INFO("Pre-computing spline grids for fast lookup...");
  
  // Define grid resolution (balance memory vs accuracy)
  const int ENERGY_BINS = 200;  // Adjust based on your needs
  const int COSZ_BINS = 100;    // Adjust based on your needs
  
  const char* neutrino_names[4] = {"nue", "anue", "numu", "anumu"};
  
  for (int i = 0; i < 4; ++i) {
    const SplineSet& splines = spline_sets[i];
    SplineGrid& grid = spline_grids[i];
    
    // Set up grid parameters
    grid.energy_min = splines.energy_min;
    grid.energy_max = splines.energy_max;
    grid.cosZ_min = splines.cosZ_min;
    grid.cosZ_max = splines.cosZ_max;
    grid.energy_bins = ENERGY_BINS;
    grid.cosZ_bins = COSZ_BINS;
    grid.energy_step = (grid.energy_max - grid.energy_min) / ENERGY_BINS;
    grid.cosZ_step = (grid.cosZ_max - grid.cosZ_min) / COSZ_BINS;
    
    // Create 2D histogram for fast lookup
    std::string hist_name = Form("spline_grid_%s", neutrino_names[i]);
    grid.energy_cosZ_grid = new TH2D(hist_name.c_str(), hist_name.c_str(),
                                     ENERGY_BINS, grid.energy_min, grid.energy_max,
                                     COSZ_BINS, grid.cosZ_min, grid.cosZ_max);
    
    // Pre-compute all grid points (this is the expensive part, done once)
    MACH3LOG_INFO("Computing {} grid ({} x {} = {} points)...", 
                  neutrino_names[i], ENERGY_BINS, COSZ_BINS, ENERGY_BINS * COSZ_BINS);
    
    for (int e_bin = 1; e_bin <= ENERGY_BINS; ++e_bin) {
      double energy = grid.energy_cosZ_grid->GetXaxis()->GetBinCenter(e_bin);
      
      for (int cz_bin = 1; cz_bin <= COSZ_BINS; ++cz_bin) {
        double cosZ = grid.energy_cosZ_grid->GetYaxis()->GetBinCenter(cz_bin);
        
        // Evaluate splines at this point
        double energy_weight = 1.0;
        double cosZ_weight = 1.0;
        
        if (energy >= splines.energy_min && energy <= splines.energy_max) {
          energy_weight = splines.energy_spline->Eval(energy);
        }
        if (cosZ >= splines.cosZ_min && cosZ <= splines.cosZ_max) {
          cosZ_weight = splines.cosZ_spline->Eval(cosZ);
        }
        
        double combined_weight = energy_weight * cosZ_weight;
        grid.energy_cosZ_grid->SetBinContent(e_bin, cz_bin, combined_weight - 1.0); // Store response
      }
    }
    
    MACH3LOG_INFO("Grid for {} completed", neutrino_names[i]);
  }
  
  MACH3LOG_INFO("All spline grids pre-computed successfully!");
}

void SampleHandlerAtm::CacheAllSplineResponses() {
  MACH3LOG_INFO("Caching spline responses for all events...");
  
  // FIX: Use nEvents member variable instead of declaring new one
  size_t numEvents = dunemcSamples.size();
  cached_spline_responses.resize(numEvents);
  
  for (size_t iEvent = 0; iEvent < numEvents; ++iEvent) {
    double energy = dunemcSamples[iEvent].rw_etru;
    double cosZ = dunemcSamples[iEvent].rw_truecz;
    
    // Cache response for all neutrino types
    for (int i = 0; i < 4; ++i) {
      const SplineGrid& grid = spline_grids[i];
      
      // Fast histogram lookup (much faster than spline evaluation)
      if (energy >= grid.energy_min && energy <= grid.energy_max &&
          cosZ >= grid.cosZ_min && cosZ <= grid.cosZ_max) {
        int bin = grid.energy_cosZ_grid->FindBin(energy, cosZ);
        cached_spline_responses[iEvent][i] = grid.energy_cosZ_grid->GetBinContent(bin);
      } else {
        cached_spline_responses[iEvent][i] = 0.0; // Outside range
      }
    }
    
    // Progress logging
    if (iEvent % (numEvents/10) == 0) {
      MACH3LOG_INFO("Cached responses for {}/{} events", iEvent, numEvents);
    }
  }
  
  responses_cached = true;
  MACH3LOG_INFO("All {} event responses cached!", numEvents);
}

// Ultra-fast flux weight calculation using cached values
double SampleHandlerAtm::GetCachedFluxWeight(int iEvent, int pdg_index) const {
  if (!responses_cached || pdg_index < 0 || pdg_index >= 4) {
    return 0.0;
  }
  return cached_spline_responses[iEvent][pdg_index];
}

// OPTIMIZED: Fast flux weight calculation with pre-computed lookup
double SampleHandlerAtm::CalculateFluxWeight(int iEvent) {
  int pdg = dunemcSamples[iEvent].nupdgUnosc;
  
  // OPTIMIZED: Fast PDG to index mapping using switch (no string operations)
  int idx;
  switch(pdg) {
    case 12:  idx = 0; break;  // nue
    case -12: idx = 1; break;  // anue
    case 14:  idx = 2; break;  // numu
    case -14: idx = 3; break;  // anumu
    default: return 0.0;       // Unknown PDG, return 0 response
  }
  
  // Use cached response (extremely fast lookup)
  return GetCachedFluxWeight(iEvent, idx);
}

// OPTIMIZED: Simplified AtmFluxShift with minimal operations
void SampleHandlerAtm::AtmFluxShift(const double* par, std::size_t iSample, std::size_t iEvent) {
  (void)iSample; // suppress unused parameter warning
  
  // Calculate the spline response for this event (now ultra-fast)
  double spline_response = CalculateFluxWeight(static_cast<int>(iEvent));
  
  // Apply the systematic: Weight = Original * (1.0 + Dial_Value * Spline_Response)
  dunemcSamples[iEvent].flux_w = original_flux_weights[iEvent] * (1.0 + (*par) * spline_response);
}

// Enhanced constructor initialization
SampleHandlerAtm::SampleHandlerAtm(std::string mc_version_, ParameterHandlerGeneric* xsec_cov_, const std::shared_ptr<OscillationHandler>&  Oscillator_) : SampleHandlerFD(mc_version_, xsec_cov_, Oscillator_) {
  
  // OPTIMIZED: Initialize spline sets array instead of individual pointers
  for (int i = 0; i < 4; ++i) {
    spline_sets[i].energy_spline = nullptr;
    spline_sets[i].cosZ_spline = nullptr;
    spline_sets[i].energy_min = 0.0;
    spline_sets[i].energy_max = 0.0;
    spline_sets[i].cosZ_min = 0.0;
    spline_sets[i].cosZ_max = 0.0;
    
    spline_grids[i].energy_cosZ_grid = nullptr;
    spline_grids[i].energy_min = 0.0;
    spline_grids[i].energy_max = 0.0;
    spline_grids[i].cosZ_min = 0.0;
    spline_grids[i].cosZ_max = 0.0;
    spline_grids[i].energy_bins = 0;
    spline_grids[i].cosZ_bins = 0;
    spline_grids[i].energy_step = 0.0;
    spline_grids[i].cosZ_step = 0.0;
  }
  
  responses_cached = false;
  
  // Initialize flux parameters
  param_atmflux_nue = 1.0;
  param_atmflux_anue = 1.0;
  param_atmflux_numu = 1.0;
  param_atmflux_anumu = 1.0;
  
  LoadFluxSplines();

  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;
  
  Initialise();
}

// OPTIMIZED: Enhanced destructor to clean up spline sets and grids
SampleHandlerAtm::~SampleHandlerAtm() {
  MACH3LOG_INFO("Cleaning up atmospheric flux splines and grids");
  
  // Clean up all splines and grids in the sets
  for (int i = 0; i < 4; ++i) {
    delete spline_sets[i].energy_spline;
    delete spline_sets[i].cosZ_spline;
    delete spline_grids[i].energy_cosZ_grid;
    
    spline_sets[i].energy_spline = nullptr;
    spline_sets[i].cosZ_spline = nullptr;
    spline_grids[i].energy_cosZ_grid = nullptr;
  }
}

void SampleHandlerAtm::Init() {
  std::cout << "MUYUAN DEBUG: SampleHandlerAtm::Init() called!" << std::endl;

  dunemcSamples.resize(nSamples,dunemc_base());
  
  IsELike = Get<bool>(SampleManager->raw()["SampleOptions"]["IsELike"],__FILE__,__LINE__);
  ExposureScaling = Get<double>(SampleManager->raw()["SampleOptions"]["ExposureScaling"],__FILE__,__LINE__);
}

void SampleHandlerAtm::SetupSplines() {
  SplineHandler = nullptr;
}

void SampleHandlerAtm::RegisterFunctionalParameters() {
  MACH3LOG_INFO("Registering atmospheric flux functional parameters");
  
  RegisterIndividualFunctionalParameter("atmflux_0", 0,
    [this](const double* par, std::size_t iEvent) { 
      this->AtmFluxShift(par, 0, iEvent);  // Pass 0 for iSample
    }
  );
  
  MACH3LOG_INFO("Registered atmflux_0 functional parameter");
}

void SampleHandlerAtm::SetupWeightPointers() {
  MACH3LOG_INFO("=== SampleHandlerAtm::SetupWeightPointers() called ===");
  
  // Setup weight pointers for all events - BASIC WEIGHTS ONLY
  for (size_t i = 0; i < dunemcSamples.size(); ++i) {
    MCSamples[i].total_weight_pointers.clear();
    
    // Add only the standard weights
    MCSamples[i].total_weight_pointers.push_back(&(dunemcSamples[i].flux_w));   // Original flux weight
    MCSamples[i].total_weight_pointers.push_back(MCSamples[i].osc_w_pointer);    // Oscillation weight  
    MCSamples[i].total_weight_pointers.push_back(&(MCSamples[i].xsec_w));        // Cross-section weight
    MCSamples[i].total_weight_pointers.push_back(&(ExposureScaling));            // Exposure scaling
    
    // NO ATMOSPHERIC FLUX LOGIC HERE - that's handled by the functional parameter system
  }
  
  MACH3LOG_INFO("=== SampleHandlerAtm::SetupWeightPointers() finished ===");
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
    double RecoCZ = -RecoNuMomentumVector.Y();
    if (std::isnan(RecoCZ)) {
      continue;
    }
    if (std::isnan(RecoENu)) {
      continue;
    }
    dunemcSamples[iEvent].rw_erec = RecoENu;
    dunemcSamples[iEvent].rw_theta = RecoCZ;
    
    std::string CurrFileName = Chain->GetCurrentFile()->GetName();
    dunemcSamples[iEvent].nupdgUnosc = GetInitPDGFromFileName(CurrFileName);
    dunemcSamples[iEvent].nupdg = GetFinalPDGFromFileName(CurrFileName);
    dunemcSamples[iEvent].OscChannelIndex = static_cast<double>(GetOscChannel(OscChannels, dunemcSamples[iEvent].nupdgUnosc, dunemcSamples[iEvent].nupdg));
    
    int M3Mode = Modes->GetModeFromGenerator(std::abs(sr->mc.nu[0].mode));
    if (!sr->mc.nu[0].iscc) M3Mode += 14;
    if (M3Mode > 15) M3Mode -= 1;
    dunemcSamples[iEvent].mode = M3Mode;
    
    dunemcSamples[iEvent].rw_isCC = sr->mc.nu[0].iscc;
    dunemcSamples[iEvent].Target = kTarget_Ar;
    
    dunemcSamples[iEvent].rw_etru = static_cast<double>(sr->mc.nu[0].E);

    TVector3 TrueNuMomentumVector = (TVector3(sr->mc.nu[0].momentum.X(),sr->mc.nu[0].momentum.Y(),sr->mc.nu[0].momentum.Z())).Unit();
    dunemcSamples[iEvent].rw_truecz = -TrueNuMomentumVector.Y();

    dunemcSamples[iEvent].flux_w = sr->mc.nu[0].genweight;
  }
  
  // setup original weight
  original_flux_weights.resize(nEntries);
  for (int i = 0; i < nEntries; i++) {
    original_flux_weights[i] = dunemcSamples[i].flux_w;
  }

  // IMPORTANT: Cache spline responses AFTER loading all events
  CacheAllSplineResponses();

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