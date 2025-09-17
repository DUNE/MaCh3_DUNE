#include "SampleHandlerAtm.h"

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wfloat-conversion"
//Standard Record includes
#include "duneanaobj/StandardRecord/StandardRecord.h"
#pragma GCC diagnostic pop

#include <cstdlib>
#include <fstream>

// Muyuan's current approach
// LoadFluxSplines() → calls LoadSingleSpline() for each neutrino type
// SetupWeightPointers() → calls CalculateFluxWeight() for each event
// CalculateFluxWeight() → evaluates splines and returns combined weight
// Original flux weight gets permanently modified: dunemcSamples[i].flux_w = original_weight * spline_weight

// MUYUAN: LoadSingleSpline function
TSpline3* SampleHandlerAtm::LoadSingleSpline(const std::string& filename, const std::string& splinename) {
  MACH3LOG_INFO("Attempting to load spline '{}' from file '{}'", splinename, filename);
  
  TFile* file = TFile::Open(filename.c_str(), "READ");
  if (!file || file->IsZombie()) {
    MACH3LOG_ERROR("Failed to open spline file: {}", filename);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  MACH3LOG_INFO("Successfully opened file: {}", filename);

  TSpline3* spline = dynamic_cast<TSpline3*>(file->Get(splinename.c_str()));
  if (!spline) {
    MACH3LOG_ERROR("Failed to retrieve spline: {} from file: {}", splinename, filename);
    MACH3LOG_ERROR("Available objects in file:");
    file->ls(); // This will print all objects in the ROOT file
    file->Close();
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  MACH3LOG_INFO("Successfully retrieved spline: {}", splinename);
  
  // Print some spline info
  MACH3LOG_INFO("Spline '{}' has {} points, range: [{}, {}]", 
                splinename, spline->GetNp(), spline->GetXmin(), spline->GetXmax());
  
  TSpline3* cloned_spline = dynamic_cast<TSpline3*>(spline->Clone());
  if (!cloned_spline) {
    MACH3LOG_ERROR("Failed to clone spline: {}", splinename);
    file->Close();
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  MACH3LOG_INFO("Successfully cloned spline: {}", splinename);
  
  file->Close();
  MACH3LOG_INFO("Closed file and completed loading spline: {}", splinename);
  
  return cloned_spline;
}

// MUYUAN: LoadFluxSplines for energy and cosine splines
void SampleHandlerAtm::LoadFluxSplines() {
  MACH3LOG_INFO("=== TESTING: Loading numu energy and cosine splines ===");
  
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
    // Print what keys are actually available
    MACH3LOG_ERROR("Available keys under FluxSplines:");
    for (auto it = SampleManager->raw()["FluxSplines"].begin(); it != SampleManager->raw()["FluxSplines"].end(); ++it) {
      MACH3LOG_ERROR("  - {}", it->first.as<std::string>());
    }
    
    // Let's just try hardcoding the path from your YAML for now
    spline_file = "Inputs/DUNE_atmospheric_spline_files/atmflux_SolMin_splines.root";
    MACH3LOG_WARN("Using hardcoded spline file path: {}", spline_file);
  }
  
  MACH3LOG_INFO("Loading test splines from: {}", spline_file);
  
  // Check if file exists
  std::ifstream file_check(spline_file);
  if (!file_check.good()) {
    MACH3LOG_ERROR("Spline file does not exist or is not readable: {}", spline_file);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  file_check.close();
  MACH3LOG_INFO("Spline file exists and is readable");
  
  // Load all atmospheric flux splines
  MACH3LOG_INFO("Loading nue splines...");
  spline_nue_E = LoadSingleSpline(spline_file, "atmflux_NuE_E");
  spline_nue_Cos = LoadSingleSpline(spline_file, "atmflux_NuE_Cos");
  
  MACH3LOG_INFO("Loading anti-nue splines...");
  spline_anue_E = LoadSingleSpline(spline_file, "atmflux_ANuE_E");
  spline_anue_Cos = LoadSingleSpline(spline_file, "atmflux_ANuE_Cos");
  
  MACH3LOG_INFO("Loading numu splines...");
  spline_numu_E = LoadSingleSpline(spline_file, "atmflux_NuM_E");
  spline_numu_Cos = LoadSingleSpline(spline_file, "atmflux_NuM_Cos");
  
  MACH3LOG_INFO("Loading anti-numu splines...");
  spline_anumu_E = LoadSingleSpline(spline_file, "atmflux_ANuM_E");
  spline_anumu_Cos = LoadSingleSpline(spline_file, "atmflux_ANuM_Cos");
  
  MACH3LOG_INFO("=== All atmospheric flux splines loaded successfully! ===");
  MACH3LOG_INFO("Spline summary:");
  MACH3LOG_INFO("  nue_E:    {} points, range [{:.3f}, {:.3f}]", spline_nue_E->GetNp(), spline_nue_E->GetXmin(), spline_nue_E->GetXmax());
  MACH3LOG_INFO("  nue_Cos:  {} points, range [{:.3f}, {:.3f}]", spline_nue_Cos->GetNp(), spline_nue_Cos->GetXmin(), spline_nue_Cos->GetXmax());
  MACH3LOG_INFO("  anue_E:   {} points, range [{:.3f}, {:.3f}]", spline_anue_E->GetNp(), spline_anue_E->GetXmin(), spline_anue_E->GetXmax());
  MACH3LOG_INFO("  anue_Cos: {} points, range [{:.3f}, {:.3f}]", spline_anue_Cos->GetNp(), spline_anue_Cos->GetXmin(), spline_anue_Cos->GetXmax());
  MACH3LOG_INFO("  numu_E:   {} points, range [{:.3f}, {:.3f}]", spline_numu_E->GetNp(), spline_numu_E->GetXmin(), spline_numu_E->GetXmax());
  MACH3LOG_INFO("  numu_Cos: {} points, range [{:.3f}, {:.3f}]", spline_numu_Cos->GetNp(), spline_numu_Cos->GetXmin(), spline_numu_Cos->GetXmax());
  MACH3LOG_INFO("  anumu_E:  {} points, range [{:.3f}, {:.3f}]", spline_anumu_E->GetNp(), spline_anumu_E->GetXmin(), spline_anumu_E->GetXmax());
  MACH3LOG_INFO("  anumu_Cos:{} points, range [{:.3f}, {:.3f}]", spline_anumu_Cos->GetNp(), spline_anumu_Cos->GetXmin(), spline_anumu_Cos->GetXmax());
}

// MUYUAN: Flux weight calculation using both energy and cosine splines
// Complete CalculateFluxWeight function to handle all neutrino types
double SampleHandlerAtm::CalculateFluxWeight(int iEvent) {
  int pdg = dunemcSamples[iEvent].nupdgUnosc;
  double energy = dunemcSamples[iEvent].rw_etru;
  double cosZ = dunemcSamples[iEvent].rw_truecz;
  
  // Initialize weights
  double energy_weight = 1.0;
  double cosZ_weight = 1.0;
  
  // Select appropriate splines based on PDG code
  TSpline3* energy_spline = nullptr;
  TSpline3* cosZ_spline = nullptr;
  std::string neutrino_type = "";
  
  switch(pdg) {
    case 12:  // nue
      energy_spline = spline_nue_E;
      cosZ_spline = spline_nue_Cos;
      neutrino_type = "nue";
      break;
    case -12: // anti-nue
      energy_spline = spline_anue_E;
      cosZ_spline = spline_anue_Cos;
      neutrino_type = "anue";
      break;
    case 14:  // numu
      energy_spline = spline_numu_E;
      cosZ_spline = spline_numu_Cos;
      neutrino_type = "numu";
      break;
    case -14: // anti-numu
      energy_spline = spline_anumu_E;
      cosZ_spline = spline_anumu_Cos;
      neutrino_type = "anumu";
      break;
    default:
      MACH3LOG_DEBUG("Unknown/unsupported PDG code: {} for event {}, returning default weight", pdg, iEvent);
      return 1.0;
  }
  
  // Check energy spline
  if (energy_spline && energy >= energy_spline->GetXmin() && energy <= energy_spline->GetXmax()) {
    energy_weight = energy_spline->Eval(energy);
  } else if (energy_spline) {
    MACH3LOG_DEBUG("{} energy {:.3f} outside spline range [{:.3f}, {:.3f}] for event {}", 
                   neutrino_type, energy, energy_spline->GetXmin(), energy_spline->GetXmax(), iEvent);
  }
  
  // Check cosine spline
  if (cosZ_spline && cosZ >= cosZ_spline->GetXmin() && cosZ <= cosZ_spline->GetXmax()) {
    cosZ_weight = cosZ_spline->Eval(cosZ);
  } else if (cosZ_spline) {
    MACH3LOG_DEBUG("{} cosZ {:.3f} outside spline range [{:.3f}, {:.3f}] for event {}", 
                   neutrino_type, cosZ, cosZ_spline->GetXmin(), cosZ_spline->GetXmax(), iEvent);
  }
  
  // Combined weight: energy spline * cosine spline
  double final_weight = energy_weight * cosZ_weight;
  
  MACH3LOG_DEBUG("Event {}: {} E={:.3f}, cosZ={:.3f}, E_weight={:.6f}, cosZ_weight={:.6f}, final_weight={:.6f}", 
                 iEvent, neutrino_type, energy, cosZ, energy_weight, cosZ_weight, final_weight);
  
  // return the response (deviation from 1.0), not the absolute weight
  return final_weight - 1.0;

  // Can implement more complex logic if needed
  // Below, we scale the final weight by a parameter specific to the neutrino type
  // double param = 1.0;

  // switch (pdg) {
  //   case 12: param = param_atmflux_nue; break;
  //   case -12: param = param_atmflux_anue; break;
  //   case 14: param = param_atmflux_numu; break;
  //   case -14: param = param_atmflux_anumu; break;
  // }

  // double weight = 1.0 + (param - 1.0) * (final_weight - 1.0);
  // return weight;
  
}

// * correct way of applying weight 
// Muyuan: apply the flux weight
// following similar logic https://github.com/DUNE/MaCh3_DUNE/blob/e47b67efeccf0dec2305d320440a540c338506b9/Samples/SampleHandlerBeamFD.cpp#L57-L63
void SampleHandlerAtm::AtmFluxShift(const double* par, std::size_t iSample, std::size_t iEvent) {
  (void)iSample; // supres unused parameter warning
  //calculate the spline response for this event
  double spline_response = CalculateFluxWeight(static_cast<int>(iEvent)); // return respone response instead of abs weight
  // Apply the systematic by modifying the flux weight directly
  // Original flux weight * (1.0 + Dial_Value * Spline_Response)

  // Apply the systematic: Weight = Original * (1.0 + Dial_Value * Spline_Response)
  dunemcSamples[iEvent].flux_w = original_flux_weights[iEvent] * (1.0 + (*par) * spline_response);
  
  // Debug output for first few calls
  static int debug_calls = 0;
  debug_calls++;
  if (debug_calls <= 5) {
    int pdg = dunemcSamples[iEvent].nupdgUnosc;
    std::string neutrino_type;
    switch(pdg) {
      case 12: neutrino_type = "nue"; break;
      case -12: neutrino_type = "anue"; break;
      case 14: neutrino_type = "numu"; break;
      case -14: neutrino_type = "anumu"; break;
      default: neutrino_type = "unknown"; break;
    }
    
    // MACH3LOG_INFO("AtmFluxShift {} Event {}: dial={:.6f}, response={:.6f}, orig_flux={:.6f}, new_flux={:.6f}", 
                  // neutrino_type, iEvent, *par, spline_response, original_flux_weights[iEvent], dunemcSamples[iEvent].flux_w);
  }
}





// Enhanced constructor initialization
SampleHandlerAtm::SampleHandlerAtm(std::string mc_version_, ParameterHandlerGeneric* xsec_cov_, const std::shared_ptr<OscillationHandler>&  Oscillator_) : SampleHandlerFD(mc_version_, xsec_cov_, Oscillator_) {
  
  // Initialize all spline pointers to nullptr
  spline_nue_E = nullptr;
  spline_nue_Cos = nullptr;
  spline_anue_E = nullptr;
  spline_anue_Cos = nullptr;
  spline_numu_E = nullptr;
  spline_numu_Cos = nullptr;
  spline_anumu_E = nullptr;
  spline_anumu_Cos = nullptr;
  
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

// Enhanced destructor to clean up all splines
SampleHandlerAtm::~SampleHandlerAtm() {
  MACH3LOG_INFO("Cleaning up atmospheric flux splines");
  
  // Clean up all splines
  delete spline_nue_E;
  delete spline_nue_Cos;
  delete spline_anue_E;
  delete spline_anue_Cos;
  delete spline_numu_E;
  delete spline_numu_Cos;
  delete spline_anumu_E;
  delete spline_anumu_Cos;
  
  // Set all pointers to nullptr for safety
  spline_nue_E = nullptr;
  spline_nue_Cos = nullptr;
  spline_anue_E = nullptr;
  spline_anue_Cos = nullptr;
  spline_numu_E = nullptr;
  spline_numu_Cos = nullptr;
  spline_anumu_E = nullptr;
  spline_anumu_Cos = nullptr;
}

void SampleHandlerAtm::Init() {
  std::cout << "MUYUAN DEBUG: SampleHandlerAtm::Init() called!" << std::endl;
  // std::cout << __FUNC__ << std::endl;

  dunemcSamples.resize(nSamples,dunemc_base());
  
  IsELike = Get<bool>(SampleManager->raw()["SampleOptions"]["IsELike"],__FILE__,__LINE__);
  ExposureScaling = Get<double>(SampleManager->raw()["SampleOptions"]["ExposureScaling"],__FILE__,__LINE__);
}

void SampleHandlerAtm::SetupSplines() {
  SplineHandler = nullptr;
}

// Muyuan Enhanced RegisterFunctionalParameters to handle all flux parameters
void SampleHandlerAtm::RegisterFunctionalParameters() {
  MACH3LOG_INFO("Registering atmospheric flux functional parameters");
  
  RegisterIndividualFunctionalParameter("atmflux_0", 0,
    [this](const double* par, std::size_t iEvent) { 
      this->AtmFluxShift(par, 0, iEvent);  // Pass 0 for iSample
    }
  );
  
  MACH3LOG_INFO("Registered atmflux_0 functional parameter");
}

// Enhanced SetupWeightPointers with proper parameter handling
// void SampleHandlerAtm::SetupWeightPointers() { // no ability to be recalled after initialization
//   MACH3LOG_INFO("=== SampleHandlerAtm::SetupWeightPointers() called ===");
  
//   // Initialize flux parameters
//   param_atmflux_nue = 1.0;
//   param_atmflux_anue = 1.0;
//   param_atmflux_numu = 1.0;
//   param_atmflux_anumu = 1.0;
  
//   std::vector<NormParameter> norm_parameters = ParHandler->GetNormParsFromSampleName(GetSampleName());
//   MACH3LOG_INFO("Found {} norm parameters for sample {}", norm_parameters.size(), GetSampleName());

//   // Calculate flux weights for all events using splines
//   MACH3LOG_INFO("Calculating flux weights for {} events", dunemcSamples.size());
  
//   std::map<int, int> pdg_counts;
//   std::map<int, std::string> pdg_names = {{12, "nue"}, {-12, "anue"}, {14, "numu"}, {-14, "anumu"}};
  
//   for (size_t i = 0; i < dunemcSamples.size(); ++i) {
//     int pdg = dunemcSamples[i].nupdgUnosc;
//     pdg_counts[pdg]++;
    
//     // Calculate the new flux weight using splines
//     double original_weight = dunemcSamples[i].flux_w;
//     double spline_weight = CalculateFluxWeight(static_cast<int>(i));
    
//     // Multiply original flux weight by spline weight
//     dunemcSamples[i].flux_w = original_weight * spline_weight;
    
//     // Debug output for first few events of each type
//     if (pdg_counts[pdg] <= 3 && pdg_names.count(pdg)) {
//       MACH3LOG_INFO("{} Event {}: original_flux={:.6f}, spline={:.6f}, final_flux={:.6f}, energy={:.3f}, cosZ={:.3f}", 
//                     pdg_names[pdg], i, original_weight, spline_weight, dunemcSamples[i].flux_w, 
//                     dunemcSamples[i].rw_etru, dunemcSamples[i].rw_truecz);
//     }
//   }
  
//   // Log statistics
//   for (const auto& pair : pdg_counts) {
//     std::string name = pdg_names.count(pair.first) ? pdg_names[pair.first] : "unknown";
//     MACH3LOG_INFO("Found {} {} events (PDG: {})", pair.second, name, pair.first);
//   }

//   // Setup weight pointers for all events
//   for (size_t i = 0; i < dunemcSamples.size(); ++i) {
//     MCSamples[i].total_weight_pointers.clear(); // Clear any existing pointers
    
//     // Add basic weights
//     MCSamples[i].total_weight_pointers.push_back(&(dunemcSamples[i].flux_w));
//     MCSamples[i].total_weight_pointers.push_back(MCSamples[i].osc_w_pointer);
//     MCSamples[i].total_weight_pointers.push_back(&(MCSamples[i].xsec_w));
//     MCSamples[i].total_weight_pointers.push_back(&(ExposureScaling));

//     // Add appropriate flux parameter based on neutrino type
//     int pdg = dunemcSamples[i].nupdgUnosc;
//     switch(pdg) {
//       case 12:  // nue
//         MCSamples[i].total_weight_pointers.push_back(&param_atmflux_nue);
//         break;
//       case -12: // anti-nue
//         MCSamples[i].total_weight_pointers.push_back(&param_atmflux_anue);
//         break;
//       case 14:  // numu
//         MCSamples[i].total_weight_pointers.push_back(&param_atmflux_numu);
//         break;
//       case -14: // anti-numu
//         MCSamples[i].total_weight_pointers.push_back(&param_atmflux_anumu);
//         break;
//       default:
//         // For unknown PDG codes, add a dummy parameter (always 1.0)
//         static double dummy_param = 1.0;
//         MCSamples[i].total_weight_pointers.push_back(&dummy_param);
//         break;
//     }
    
//     // Add any additional norm parameters from the parameter handler
//     for (const auto& param : norm_parameters) {
//       std::string param_name = param.name;
//       if (param_name.find("atmflux") != std::string::npos) {
//         MCSamples[i].total_weight_pointers.push_back(ParHandler->RetPointer(param.index));
//         if (i == 0) { // Log only once
//           MACH3LOG_INFO("Adding norm parameter {} with index {}", param_name, param.index);
//         }
//       }
//     }
//   }
  
//   MACH3LOG_INFO("=== SampleHandlerAtm::SetupWeightPointers() finished ===");
// }

// Muyuan brought this back to original form. this method only sets up the pointers
// at the beginning of the fit. 
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

int SampleHandlerAtm::SetupExperimentMC() { // Heart of the code
  // Loads CAF files constraining DUNE atmospheric neutrino events
  // Loop through ~6.6M events in total
  // Extracts key info for each event: RecoENu, RecoCZ(cosine zenith), nupdg, nupdgUnosc, mode, flux_w
  int CurrErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  
  caf::StandardRecord* sr = new caf::StandardRecord();

  TChain* Chain = new TChain("cafTree"); // load CAF tree
  for (size_t iSample=0;iSample<mc_files.size();iSample++) {
    Chain->Add(mc_files[iSample].c_str());
  }
  
  Chain->SetBranchStatus("*", 1);
  Chain->SetBranchAddress("rec", &sr);

  int nEntries = static_cast<int>(Chain->GetEntries()); // Get number of entries in the chain
  dunemcSamples.resize(nEntries);
 
  for (int iEvent=0;iEvent<nEntries;iEvent++) { // Loop through all events
    Chain->GetEntry(iEvent);

    if ((iEvent % (nEntries/10))==0) {
      MACH3LOG_INFO("\tProcessing event: {}/{}",iEvent,nEntries);
    }
    
    if(sr->common.ixn.pandora.size() != 1) {
      MACH3LOG_WARN("Skipping event {}/{} -> Number of neutrino slices found in event: {}",iEvent,nEntries,sr->common.ixn.pandora.size());
      continue;
    }
    
    TVector3 RecoNuMomentumVector; // Reconstructed neutrino momentum vector
    double RecoENu;
    if (IsELike) { // ELike events have different reconstructed neutrino energy and direction
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
    dunemcSamples[iEvent].rw_erec = RecoENu; // Reconstructed neutrino energy
    dunemcSamples[iEvent].rw_theta = RecoCZ; // Reconstructed cosine zenith angle
    
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
  // setup original weight
  original_flux_weights.resize(nEntries);
  for (int i = 0; i < nEntries; i++) {
    original_flux_weights[i] = dunemcSamples[i].flux_w;
  }

  delete Chain;
  gErrorIgnoreLevel = CurrErrorLevel;


  return nEntries; // Return the number of events processed
}

void SampleHandlerAtm::SetupFDMC() { // Setup FD MC samples
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
  double* KinematicValue; // Pointer to the kinematic value for the event

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