#include "SampleHandlerBeamNDGAr.h"

SampleHandlerBeamNDGAr::SampleHandlerBeamNDGAr(std::string mc_version_, ParameterHandlerGeneric* ParHandler_) : SampleHandlerFD(mc_version_, ParHandler_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;

  KinematicVectors = &KinematicVectorsDUNE;
  ReversedKinematicVectors = &ReversedKinematicVectorsDUNE;

  Initialise();
}

SampleHandlerBeamNDGAr::~SampleHandlerBeamNDGAr() {
}

void SampleHandlerBeamNDGAr::Init() {
  energy_resolution_threshold = SampleManager->raw()["DetectorVariables"]["energy_resolution_threshold"].as<double>(); //NK energy resolution threshold, total as a fraction of momentum
  pixel_spacing = SampleManager->raw()["DetectorVariables"]["pixel_spacing"].as<double>(); //NK pixel spacing in mm to find num hits in y,z plane
  spatial_resolution = SampleManager->raw()["DetectorVariables"]["spatial_resolution"].as<double>(); //NK spatial resolution in mm to find  in y,z plane
  adc_sampling_frequency = SampleManager->raw()["DetectorVariables"]["adc_sampling_frequency"].as<double>(); //NK sampling frequency for ADC - needed to find timing resolution and spatial resolution in x dir in MHz
  drift_velocity = SampleManager->raw()["DetectorVariables"]["drift_velocity"].as<double>(); //NK drift velocity of electrons in gas - needed to find timing resolution and spatial resolution in x dir in cm/microsecond
  downsampling = SampleManager->raw()["DetectorVariables"]["downsampling"].as<double>(); //JM downsampling fraction
  do_geometric_correction = SampleManager->raw()["DetectorVariables"]["DoGeometricCorrection"].as<bool>(); //JM whether to apply geometric correction
  TPCFidLength = SampleManager->raw()["DetectorVariables"]["TPCFidLength"].as<double>();
  TPCFidRadius = SampleManager->raw()["DetectorVariables"]["TPCFidRadius"].as<double>();
  B_field = SampleManager->raw()["DetectorVariables"]["B_field"].as<double>(); //NK B field value in T
  TPCInstrumentedLength = SampleManager->raw()["DetectorVariables"]["TPCInstrumentedLength"].as<double>();
  TPCInstrumentedRadius = SampleManager->raw()["DetectorVariables"]["TPCInstrumentedRadius"].as<double>();
  ECALSciX0 = SampleManager->raw()["DetectorVariables"]["ECALSciX0"].as<double>();
  ECALBarrelForwardDepth = SampleManager->raw()["DetectorVariables"]["ECALBarrelForwardDepth"].as<double>();
  ECALBarrelBackwardDepth = SampleManager->raw()["DetectorVariables"]["ECALBarrelBackwardDepth"].as<double>();
  nECALBackSegments = SampleManager->raw()["DetectorVariables"]["nECALBackSegments"].as<int>();
  ECALEndCapDepth = SampleManager->raw()["DetectorVariables"]["ECALEndCapDepth"].as<double>();
  interaction_model = SampleManager->raw()["DetectorVariables"]["interaction_model"].as<std::string>();

  beamNDGArSampleDetails.resize(static_cast<size_t>(GetNsamples()));
  
  auto EnabledSamples = Get<std::vector<std::string>>(SampleManager->raw()["Samples"], __FILE__ , __LINE__);

  for (size_t i = 0; i < static_cast<size_t>(GetNsamples()); i++){
    const auto TempTitle = EnabledSamples[i];
    beamNDGArSampleDetails[i].isFHC = SampleManager->raw()[TempTitle]["DUNESampleBools"]["isFHC"].as<double>();
    beamNDGArSampleDetails[i].iselike = SampleManager->raw()[TempTitle]["DUNESampleBools"]["iselike"].as<bool>();
    beamNDGArSampleDetails[i].pot = SampleManager->raw()[TempTitle]["POT"].as<double>();

    beamNDGArSampleDetails[i].norm_s = 1;
    beamNDGArSampleDetails[i].pot_s = (beamNDGArSampleDetails[i].pot)/(downsampling*1e21);

    MACH3LOG_INFO("Setting up beam ND sample {}", GetSampleTitle(static_cast<int>(i)));
    MACH3LOG_INFO("- isFHC: {}", beamNDGArSampleDetails[i].isFHC);
    MACH3LOG_INFO("- iselike: {}", beamNDGArSampleDetails[i].iselike);
  }

  MACH3LOG_INFO("-------------------------------------------------------------------");
}

void SampleHandlerBeamNDGAr::SetupSplines() {
  ///@todo move all of the spline setup into core
  if(ParHandler->GetNumParamsFromSampleName(SampleHandlerName, kSpline) > 0){
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", ParHandler->GetNumParamsFromSampleName(SampleHandlerName, kSpline));
    SplineHandler = std::unique_ptr<BinnedSplineHandler>(new BinnedSplineHandlerDUNE(ParHandler,Modes.get()));
    InitialiseSplineObject();
  }
  else{
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or evaluate splines", ParHandler->GetNumParamsFromSampleName(SampleHandlerName, kSpline));
    SplineHandler = nullptr;
  }

  return;
}

void SampleHandlerBeamNDGAr::AddAdditionalWeightPointers() {
  for (size_t i = 0; i < dunendgarmcFitting.size(); ++i) {
    // MACH3LOG_INFO("pot: {}\nnorm: {}\nberpa: {}\nflux: {}\ngeom: {}", dunendgarmcFitting[i].pot_s, dunendgarmcFitting[i].norm_s, dunendgarmcFitting[i].rw_berpaacvwgt, dunendgarmcFitting[i].flux_w, dunendgarmcPlotting[i].geometric_correction);
    MCSamples[i].total_weight_pointers.push_back(&(dunendgarmcFitting[i].pot_s));
    MCSamples[i].total_weight_pointers.push_back(&(dunendgarmcFitting[i].norm_s));
    MCSamples[i].total_weight_pointers.push_back(&(dunendgarmcFitting[i].rw_berpaacvwgt));
    MCSamples[i].total_weight_pointers.push_back(&(dunendgarmcFitting[i].flux_w));
    MCSamples[i].total_weight_pointers.push_back(&(dunendgarmcPlotting[i].geometric_correction));
  }
}

void SampleHandlerBeamNDGAr::CleanMemoryBeforeFit() {
  CleanVector(dunendgarmcPlotting);
}

bool SampleHandlerBeamNDGAr::isCoordOnTrack(int charge, double ycoord, double zcoord, double centre_circle_y, double centre_circle_z, double theta_start, double theta_spanned) {
  double theta_coord = atan2(ycoord - centre_circle_y, zcoord - centre_circle_z);
  bool iscoordinTPC = (ycoord - TPC_centre_y)*(ycoord - TPC_centre_y)+(zcoord - TPC_centre_z)*(zcoord - TPC_centre_z) < TPCInstrumentedRadius*TPCInstrumentedRadius;
  bool iscoordinarc;

  if (charge == 1) {
    double theta_end = theta_start + theta_spanned;
    if (theta_coord < theta_start) theta_coord += 2*M_PI;
    iscoordinarc = theta_coord > theta_start && theta_coord < theta_end;
  }
  else {
    double theta_end = theta_start - theta_spanned;
    if (theta_coord > theta_start) theta_start -= 2*M_PI;
    iscoordinarc = theta_coord < theta_start && theta_coord > theta_end;
  }

  return iscoordinTPC && iscoordinarc;
}

double SampleHandlerBeamNDGAr::FindNHits(double pixel_spacing_cm, double centre_circle_y, double centre_circle_z, double rad_curvature, double theta_start, double theta_spanned, int charge) {

  int num_intersections = 0;
  int num_vertices = 0;

  //Define yz pixel grid which covers entire tpc cross section
  int numpixelboundaries = static_cast<int>(floor(TPCInstrumentedRadius*2/(pixel_spacing_cm)))+1;
  if (numpixelboundaries%2 == 0) numpixelboundaries += 1;
  int pixelmin = -(numpixelboundaries-1)/2;
  int pixelmax = (numpixelboundaries-1)/2;

  //Loop through all pixel boundaries
  for (int pixel=pixelmin; pixel<=pixelmax; pixel++) {

    //Check if corresponding y boundary is crossed
    double ycoord = TPC_centre_y + pixel*pixel_spacing_cm;
    double quadratic_ineq_y = rad_curvature*rad_curvature-(ycoord-centre_circle_y)*(ycoord-centre_circle_y);

    if (quadratic_ineq_y >= 0) {
      double zcoord1 = centre_circle_z + std::sqrt(quadratic_ineq_y);
      if (isCoordOnTrack(charge, ycoord, zcoord1, centre_circle_y, centre_circle_z, theta_start, theta_spanned)) num_intersections++;
      if (fmod(std::abs(zcoord1 - TPC_centre_z),pixel_spacing_cm)==0) num_vertices++;

      double zcoord2 = centre_circle_z - std::sqrt(quadratic_ineq_y);
      if (isCoordOnTrack(charge, ycoord, zcoord2, centre_circle_y, centre_circle_z, theta_start, theta_spanned)) num_intersections++;
      if (fmod(std::abs(zcoord2 - TPC_centre_z),pixel_spacing_cm)==0) num_vertices++;
    }

    //Check if corresponding z boundary is crossed
    double zcoord = TPC_centre_z + pixel*pixel_spacing_cm;
    double quadratic_ineq_z = rad_curvature*rad_curvature-(zcoord-centre_circle_z)*(zcoord-centre_circle_z);

    if (quadratic_ineq_z >= 0) {
      double ycoord1 = centre_circle_y + std::sqrt(quadratic_ineq_z);
      if (isCoordOnTrack(charge, ycoord1, zcoord, centre_circle_y, centre_circle_z, theta_start, theta_spanned)) num_intersections++;
      if (fmod(std::abs(ycoord1 - TPC_centre_y),pixel_spacing_cm)==0) num_vertices++;

      double ycoord2 = centre_circle_y - std::sqrt(quadratic_ineq_z);
      if (isCoordOnTrack(charge, ycoord2, zcoord, centre_circle_y, centre_circle_z, theta_start, theta_spanned)) num_intersections++;
      if (fmod(std::abs(ycoord2 - TPC_centre_y),pixel_spacing_cm)==0) num_vertices++;
    }
  }
  return num_intersections - num_vertices;
}

double SampleHandlerBeamNDGAr::CalcBeta(double p_mag, double& bg, double& gamma, double mass){ //calculate beta (v/c)
  bg = p_mag/mass; //beta*gamma
  gamma = std::sqrt(1+bg*bg); //gamma
  double beta = bg/gamma; //beta (velocity)
  return beta;
}

int SampleHandlerBeamNDGAr::GetChargeFromPDG(const int pdg) {
  if (pdg > 1000000000) return 0;
  switch(pdg) {
    case 11:    // e-
    case 13:    // mu-
    case 15:    // tau-
    case -211:  // pi-
    case -321:  // K-
    case -213:  // rho-
    case -411:  // D-
    case -431:  // Ds-
    case -2212: // p-
    case 3112:  // Sig-
    case -3222: // anti Sig+
    case 3312:  // Xi-
    case 3334:  // Omega-
      return -1;

    case -11:   // e+
    case -13:   // mu+
    case -15:   // tau+
    case 211:   // pi+
    case 321:   // K+
    case 213:   // rho+
    case 411:   // D+
    case 431:   // Ds+
    case 2212:  // p+
    case -3112: // anti Sig-
    case 3222:  // Sig+
    case -3312: // Xi+
    case -3334: // Omega+
    case 4122:  // Lambda_c+
    case 4212:  // Sigma_c+
      return 1;

    case 4222:  // Sigma_c++
      return 2;

    default:
      switch(std::abs(pdg)) {
        case 12:    // nue
        case 14:    // numu
        case 16:    // nutau
        case 22:    // photon
        case 111:   // pi0
        case 113:   // rho0
        case 221:   // eta
        case 331:   // eta'
        case 311:   // K0
        case 130:   // K0L
        case 310:   // K0S
        case 223:   // omega
        case 421:   // D0
        case 2112:  // n
        case 3122:  // Lambda
        case 3212:  // Sig0
        case 3322:  // Xi0
          return 0; 
        default:
          MACH3LOG_ERROR("No saved charge for PDG: {}", pdg);
          throw MaCh3Exception(__FILE__,__LINE__);
      }
  }
}

void SampleHandlerBeamNDGAr::EraseDescendants(int motherID, std::unordered_map<int, std::vector<int>>& mother_to_daughter_ID) {
  for (int daughterID : mother_to_daughter_ID[motherID]) {
    EraseDescendants(daughterID, mother_to_daughter_ID);
  }
  mother_to_daughter_ID.erase(motherID);
}

// Removes descendants of primID which can be reconstructed from curvature from mother_to_daughter_ID map
bool SampleHandlerBeamNDGAr::CurvatureResolutionFilter(int id, std::unordered_map<int, std::vector<int>>& mother_to_daughter_ID, 
                                                       const std::unordered_map<int, size_t>& ID_to_index, dunemc_plotting& plotting_vars, double pixel_spacing_cm) {
  if (IsResolvedFromCurvature(plotting_vars, ID_to_index.at(id), pixel_spacing_cm)) { // If mother can be reconstructed from curvature, remove her and her descendants
    size_t index = ID_to_index.at(id);
    int motherID = _MCPMotherTrkID->at(index);
    if (motherID != 0) {
      EraseDescendants(id, mother_to_daughter_ID);
    }
    return true;
  }
  auto& daughters = mother_to_daughter_ID.at(id); 
  std::vector<int> filtered_daughters = {};
  for (int daughterID : daughters) { // Check secondary tree for any others which can be reconstructed from curvature
    if (!CurvatureResolutionFilter(daughterID, mother_to_daughter_ID, ID_to_index, plotting_vars, pixel_spacing_cm)) {
      filtered_daughters.push_back(daughterID); 
    }
  }
  daughters = std::move(filtered_daughters); // Update with filtered daughters
  return false;
}

double SampleHandlerBeamNDGAr::GetDCalBoundary(const std::vector<double>& startpos, const std::vector<double>& dir, size_t& boundary_index) {
  double dmin = std::numeric_limits<double>::max();
  int best_plane = -1;
  int best_segment = -1;

  for (size_t iPlane=0; iPlane < outerECalP.size(); iPlane++) {
    double num = 0;
    double denom = 0;
    double norm = 0;
    for (size_t j=0; j<3; j++) {
      num += (outerECalP[iPlane][j] - startpos[j]) * outerECalA[iPlane][j];
      denom += dir[j] * outerECalA[iPlane][j];
      norm += dir[j] * dir[j];
    }
    if (denom == 0) continue;

    double d = num * std::sqrt(norm) / denom;
    // Must cross plane with d>0 (moving forwards)
    if (d < 0) continue;

    std::vector<double> intersection(3, 0.);
    for (size_t j=0; j<3; j++) {
      intersection[j] = startpos[j] + d*dir[j]/std::sqrt(norm);
    }

    if (static_cast<int>(iPlane) < _NumCalSides) { // Barrel plane
      int segment = GetCalSegment(intersection[1], intersection[2]);

      // If crossing the barrel in the right segment within the ecal length, must be correct segment
      if ((segment == static_cast<int>(iPlane)) && (std::abs(intersection[0]) <= ECALEndCapEnd)) {
        boundary_index = iPlane;
        return d / ECALSciX0;
      }

      if (d < dmin) { // Keep track of best barrel segment
        dmin = d;
        best_plane = static_cast<int>(iPlane);
        best_segment = segment;
      }
    }

    // After trying all barrel planes, must be end cap or dividing plane (if using shorter calorimeter), but dividing planes don't live in the vector
    // If end cap crossing is within the depth radially, must be the correct one
    double depth = GetCalDepth(0., intersection[1], intersection[2]);
    if (depth != 999.) {
      boundary_index = iPlane;
      return d / ECALSciX0;
    }
  }

  // If we reach here, the track must be pointing at a dividing plane (between the front and back calorimeter)
  // There are two, but we can tell which by comparing our best_segment with our best_plane (can see this from a sketch)
  if (best_plane == -1 || best_segment == -1) return -1;

  size_t dividing_plane_idx = 0;
  if (best_segment > best_plane) {
    dividing_plane_idx = 1;
  }

  // Calculate distance to dividing plane
  double num = 0;
  double denom = 0;
  double norm = 0;
  for (size_t j=0; j<3; j++) {
    num += (dividingPlaneP[dividing_plane_idx][j] - startpos[j]) * dividingPlaneA[dividing_plane_idx][j];
    denom += dir[j] * dividingPlaneA[dividing_plane_idx][j];
    norm += dir[j] * dir[j];
  }
  double d = num * std::sqrt(norm) / denom;
  boundary_index = static_cast<size_t>(_NumCalSides)+2;

  return d / ECALSciX0;
}

bool SampleHandlerBeamNDGAr::IsPrimContained(int id, const std::unordered_map<int, std::vector<int>>& mother_to_daughter_ID, const std::unordered_map<int, size_t>& ID_to_index,
                                             const std::unordered_map<int, std::vector<double>>& eID_to_showerstart, dunemc_plotting& plotting_vars) {
  size_t idx = ID_to_index.at(id);
  int pdg = _MCPPDG->at(idx);

  // Don't require containment of secondary neutrons/neutrinos/nuclei
  if ((std::abs(pdg) == 2112) || (std::abs(pdg) == 12) || (std::abs(pdg) ==14) || (std::abs(pdg) == 16) || (std::abs(pdg) > 1000000000)) return true;

  const std::vector<double> startpos = {_MCPStartX->at(idx), _MCPStartY->at(idx), _MCPStartZ->at(idx)};
  const std::vector<double> endpos = {_MCPEndX->at(idx), _MCPEndY->at(idx), _MCPEndZ->at(idx)};

  bool leavesTPC = (GetCalDepth(startpos[0], startpos[1], startpos[2]) == -999.) && (GetCalDepth(endpos[0], endpos[1], endpos[2]) != -999.);

  // Only check for containment particles which leave the TPC - otherwise move on to secondaries
  if (leavesTPC) {
    bool isContained = true;

    // First deal with EM showers. Have identified region of [ecalP, D_wall] space with resolution better than different thresholds in external studies.
    if ((std::abs(pdg) == 22) || (std::abs(pdg) == 11) || (std::abs(pdg) == 111)) {
      // Get momentum at the point of ecal entry
      const std::vector<double> ecalP = {_MCPCalPX->at(idx)/1000., _MCPCalPY->at(idx)/1000., _MCPCalPZ->at(idx)/1000.};
      if (std::isnan(ecalP[0])) {
        // These shouldn't exist in priciple - should only be particles which end in the ECal without making a step there - shower not recorded for plots but considered contained
        MACH3LOG_WARN("Particle (pdg {}) which leaves TPC has no stored ecal momentum. Depth = {} cm.", pdg, GetCalDepth(endpos[0], endpos[1], endpos[2]));
        return true;
      }
      double ecalE = std::sqrt(ecalP[0]*ecalP[0] + ecalP[1]*ecalP[1] + ecalP[2]*ecalP[2] + MaCh3Utils::GetMassFromPDG(pdg)*MaCh3Utils::GetMassFromPDG(pdg));

      // Now find the shower start point
      std::vector<double> showerstart;
      double a_par, b_par, c_par;

      double px = _MCPStartPX->at(idx)/1000.;
      double py = _MCPStartPY->at(idx)/1000.;
      double pz = _MCPStartPZ->at(idx)/1000.;
      double ptot = std::sqrt(px*px + py*py + pz*pz);

      plotting_vars.shower_energy.push_back(ecalE);
      plotting_vars.shower_bangle.push_back(acos(px/ptot)*180/M_PI);
      plotting_vars.shower_pdg.push_back(pdg);
      if (_MCPEndProcess->at(idx) != "conv") plotting_vars.shower_isconv.push_back(false);
      else plotting_vars.shower_isconv.push_back(true);

      // For electrons, start point is the location of the first hit above threshold. This is stored in eID_to_showerstart.
      if (std::abs(pdg) == 11) {
        if (eID_to_showerstart.find(id) != eID_to_showerstart.end()) showerstart = eID_to_showerstart.at(id);
        else {
          if (std::abs(GetCalDepth(endpos[0], endpos[1], endpos[2])) == 999)  {
            MACH3LOG_WARN("Electron of energy {} GeV has no hits above threshold and does not stop in the ECal - considered uncontained.", ecalE);
            plotting_vars.shower_dcalboundary.push_back(std::numeric_limits<double>::quiet_NaN());
            plotting_vars.shower_iscontained.push_back(false);
            plotting_vars.shower_cosnorm.push_back(std::numeric_limits<double>::quiet_NaN());
            return false;
          }
          else showerstart = endpos;
        }
      }
      else showerstart = endpos; // For photons and pi0s and enclosed electrons with no hits above threshold, start point is end of trajectory (pair produce/decay)

      if (GetCalDepth(showerstart[0], showerstart[1], showerstart[2]) == 999.) {
        plotting_vars.shower_dcalboundary.push_back(std::numeric_limits<double>::quiet_NaN());
        plotting_vars.shower_cosnorm.push_back(std::numeric_limits<double>::quiet_NaN());

        // End position only good approximator of shower start for photons which pair produce. This is dominant above 20 MeV. Below this energy, showers which 'escape' haven't really - they have just deposited energy by Compton scattering instead, but the photon leaves with some small remaining energy.
        if (ecalE > 0.020) {
          plotting_vars.shower_iscontained.push_back(false);
          // MACH3LOG_WARN("Shower of energy {}, pdg {} starts beyond calorimeter boundary - treating as uncontained - startpos ({}, {}, {}), showerstart ({}, {}, {})", ecalE, pdg, startpos[0], startpos[1], startpos[2], showerstart[0], showerstart[1], showerstart[2]);
          return false;
        }
        // Otherwise, isContained remains true
        plotting_vars.shower_iscontained.push_back(true);
      }
      else  {
        // Find the parameters governing the containment threshold curve
        int resolution_threshold_int = static_cast<int>(std::round(energy_resolution_threshold*100));
        if (threshold_to_containment_params.find(std::abs(pdg)) != threshold_to_containment_params.end() && threshold_to_containment_params.at(std::abs(pdg)).find(resolution_threshold_int) != threshold_to_containment_params.at(std::abs(pdg)).end()) {
          a_par = threshold_to_containment_params.at(std::abs(pdg)).at(resolution_threshold_int)[0];
          b_par = threshold_to_containment_params.at(std::abs(pdg)).at(resolution_threshold_int)[1];
          c_par = threshold_to_containment_params.at(std::abs(pdg)).at(resolution_threshold_int)[2];
        }
        else {
          MACH3LOG_ERROR("No parameters available for determining shower containment at energy resolution threshold of {}", energy_resolution_threshold);
          throw MaCh3Exception(__FILE__,__LINE__);
        }

        // Get D_wall
        size_t boundary_index = 999;
        double d_boundary = GetDCalBoundary(showerstart, ecalP, boundary_index);

        // Check if contained (above the containment threshold curve)
        isContained = d_boundary >= a_par * std::pow(ecalE, b_par) + c_par;

        // Fill shower-level kinpars
        if (boundary_index == 999) {
          MACH3LOG_WARN("No boundary found for shower of pdg {}", pdg);
          MACH3LOG_WARN("Start Position ({} {}, {})", showerstart[0], showerstart[1], showerstart[2]);
          MACH3LOG_WARN("Direction ({}, {}, {})", ecalP[0], ecalP[1], ecalP[2]);
          plotting_vars.shower_cosnorm.push_back(std::numeric_limits<double>::quiet_NaN());
        }
        else if (boundary_index < outerECalA.size()) {
          std::vector<double> normal = outerECalA[boundary_index];
          double dotProd = ecalP[0]*normal[0] + ecalP[1]*normal[1] + ecalP[2]*normal[2];
          plotting_vars.shower_cosnorm.push_back(dotProd/ptot);
        }
        else {
          plotting_vars.shower_cosnorm.push_back(std::numeric_limits<double>::quiet_NaN());
        }

        plotting_vars.shower_dcalboundary.push_back(d_boundary);
        plotting_vars.shower_iscontained.push_back(isContained);
      }
    }
    else { // For non-EM showers, need the particle to stop within the calorimeter and not re-interact
      isContained = GetCalDepth(endpos[0], endpos[1], endpos[2]) != 999.;
      std::string end_process = _MCPEndProcess->at(idx);

      // For protons, not re-interacting means the end process is 'hIoni'
      if (pdg == 2212) isContained = isContained && (end_process == "hIoni");
      // For everything else, not re-interacting means the end process is 'Decay'
      else if (std::abs(pdg) == 211 || std::abs(pdg) == 321) isContained = isContained && (end_process == "Decay");
    }
    // No need to check secondaries if uncontained - just return false
    if (!isContained) return false;
  }
  // Check secondaries - if any are uncontained, just return false
  auto& daughters = mother_to_daughter_ID.at(id); 
  for (int daughterID : daughters) {
    if (!IsPrimContained(daughterID, mother_to_daughter_ID, ID_to_index, eID_to_showerstart, plotting_vars)) return false;
  }

  return true; // No daughters are uncontained, so return true
}

int SampleHandlerBeamNDGAr::GetCalSegment(double y, double z) {
  double theta = atan2(y, z);
  // Convert to 0 to 2pi
  theta = (theta < 0) ? theta + 2*M_PI : theta;
  int segment = static_cast<int>(std::round(theta/(M_PI/6.))) % _NumCalSides;
  return segment;
}

// Returns depth of a coordinate position in the dodecoganol ecal, or +/- 999 if it is beyond/within the ecal boundary
double SampleHandlerBeamNDGAr::GetCalDepth(double x, double y, double z) {
  x = x - TPC_centre_x;
  y = y - TPC_centre_y;
  z = z - TPC_centre_z;
  double r = std::sqrt(y*y + z*z);

  if (std::abs(x) > ECALEndCapEnd) return 999.; // beyond ecal length-wise
  if ((r < ECALInnerRadius) && (std::abs(x) >= ECALEndCapStart)) return std::abs(x) - ECALEndCapStart; // give depth relative to end cap

  double theta = atan2(y, z);
  theta = (theta < 0) ? theta + 2*M_PI : theta;
  int segment = GetCalSegment(y, z);
  double theta_segment = segment*M_PI/6.;
  double projected_r = r*cos(theta_segment - theta);

  double outer_radius = ECALOuterFrontRadius;
  if (std::find(ECalBackSegments.begin(), ECalBackSegments.end(), segment) != ECalBackSegments.end()) {
    outer_radius = ECALOuterBackRadius;
  }

  if (projected_r > outer_radius) return 999.; // beyond ecal radially
  if (projected_r < ECALInnerRadius) return -999.; // in tpc

  return projected_r - ECALInnerRadius; // give depth relative to barrel
}

// Calculate the layer from the depth. Barrel: 8x0.673cm, 34x1.142cm. Endcap: 6x0.673cm, 36x1.142cm.
// double SampleHandlerBeamNDGAr::DepthToLayer(double depth, double r) {
//   int n_thin_layers;
//   if (r < ECALInnerRadius) n_thin_layers = _NEndCapHG;
//   else n_thin_layers = _NBarrelHG;

//   double hg_width = _HGAbsWidth+_HGSciWidth+_HGBoardWidth;
//   double lg_width = _LGAbsWidth+_LGSciWidth;

//   double layer;
//   if (depth < static_cast<double>(n_thin_layers)*hg_width) {
//     layer = static_cast<int>(depth/hg_width) + 1;
//   } else {
//     layer = static_cast<int>((depth-static_cast<double>(n_thin_layers)*hg_width)/lg_width) + 1 + n_thin_layers;
//   }

//   return layer;
// }

// double SampleHandlerBeamNDGAr::CalcEDepCal(int motherID, std::unordered_map<int, std::vector<int>>& mother_to_daughter_ID, const std::unordered_map<int, std::vector<double>>& ID_to_ECalDep, const int tot_layers) {
//   auto it = mother_to_daughter_ID.find(motherID);
//   double EDepCrit = 0.;
//   if (it != mother_to_daughter_ID.end()) {
//     for (int i_layer=tot_layers-crit_layers; i_layer<tot_layers; i_layer++) {
//       EDepCrit += ID_to_ECalDep.at(motherID)[static_cast<size_t>(i_layer)];
//     }
//     for (int daughterID : it->second) {
//       EDepCrit += CalcEDepCal(daughterID, mother_to_daughter_ID, ID_to_ECalDep, tot_layers);
//     }
//   }
//   return EDepCrit;
// }

bool SampleHandlerBeamNDGAr::IsResolvedFromCurvature(dunemc_plotting& plotting_vars, size_t i_particle, double pixel_spacing_cm){
  // Get particle properties from Anatree
  double xstart = _MCPStartX->at(i_particle);
  double ystart = _MCPStartY->at(i_particle);
  double zstart = _MCPStartZ->at(i_particle);
  double start_radius = std::sqrt((ystart-TPC_centre_y)*(ystart-TPC_centre_y) + (zstart-TPC_centre_z)*(zstart-TPC_centre_z));
  double start_length = xstart - TPC_centre_x;
  bool starts_in_tpc = std::abs(start_length)<=TPCInstrumentedLength && start_radius<=TPCInstrumentedRadius;
  if (!starts_in_tpc) return false;

  int pdg = _MCPPDG->at(i_particle);
  if (pdg > 1000000000) return false;
  int charge = GetChargeFromPDG(pdg);
  double mass = MaCh3Utils::GetMassFromPDG(pdg);
  double p_x = _MCPStartPX->at(i_particle)/1000.;
  double p_y = _MCPStartPY->at(i_particle)/1000.;
  double p_z = _MCPStartPZ->at(i_particle)/1000.;
  double p_beam = p_x*BeamDirection[0] + p_y*BeamDirection[1] + p_z*BeamDirection[2];

  double energy = std::sqrt(p_x*p_x + p_y*p_y + p_z*p_z + mass*mass);
  double transverse_mom = std::sqrt(p_y*p_y + p_z*p_z); //transverse to B-field
  double mom_tot = std::sqrt(p_x*p_x + transverse_mom*transverse_mom);

  double xend = _MCPEndX->at(i_particle);
  double yend = _MCPEndY->at(i_particle);
  double zend = _MCPEndZ->at(i_particle);
  double end_radius = std::sqrt((yend-TPC_centre_y)*(yend-TPC_centre_y) + (zend-TPC_centre_z)*(zend-TPC_centre_z)); 
  double end_length = xend-TPC_centre_x;
  double end_ecaldepth = GetCalDepth(xend, yend, zend);

  bool stops_in_tpc = std::abs(end_length)<=TPCInstrumentedLength && end_radius<=TPCInstrumentedRadius;
  bool stops_before_ecal = end_ecaldepth == -999.;
  bool stops_beyond_ecal = end_ecaldepth == 999.;
  bool stops_in_ecal = !(stops_before_ecal || stops_beyond_ecal);

  //Fill particle-level kinematic variables with default or actual (if possible at this stage) values
  if (_MCPMotherTrkID->at(i_particle) == 0) {
    double invis_mass = (pdg == 2212 || pdg == 2112) ? mass : 0.; 
    plotting_vars.prim_evis[i_particle] = std::sqrt(energy*energy - invis_mass*invis_mass); // don't include rest mass energy for p/n
    plotting_vars.prim_momentum[i_particle] = mom_tot;
    plotting_vars.prim_transversemomentum[i_particle] = transverse_mom; //momentum transverse to B-field
    plotting_vars.prim_bangle[i_particle] = acos(p_x/mom_tot)*180/M_PI; //angle to B-field
    plotting_vars.prim_beamangle[i_particle] = acos(p_beam/mom_tot)*180/M_PI;

    plotting_vars.prim_startx[i_particle] = start_length;
    plotting_vars.prim_startr2[i_particle] = start_radius*start_radius;
    plotting_vars.prim_endr[i_particle] = end_radius;
    plotting_vars.prim_enddepth[i_particle] = end_ecaldepth;
    plotting_vars.prim_endx[i_particle] = end_length;
    plotting_vars.prim_endy[i_particle] = yend-TPC_centre_y;
    plotting_vars.prim_endz[i_particle] = zend-TPC_centre_z;

    plotting_vars.prim_isstoppedingap[i_particle] = !stops_in_tpc && stops_before_ecal;
    plotting_vars.prim_isstoppedinbarrelgap[i_particle] = !stops_in_tpc && stops_before_ecal && std::abs(end_length)<=TPCInstrumentedLength; 
    plotting_vars.prim_isstoppedinendgap[i_particle] = !stops_in_tpc && stops_before_ecal && std::abs(end_length)>TPCInstrumentedLength; 
    plotting_vars.prim_isstoppedinecal[i_particle] = stops_in_ecal;
    plotting_vars.prim_isstoppedinbarrel[i_particle] = stops_in_ecal && end_radius>=ECALInnerRadius;
    plotting_vars.prim_isstoppedinendcap[i_particle] = stops_in_ecal && end_radius<ECALInnerRadius;
    plotting_vars.prim_isstoppedintpc[i_particle] = stops_in_tpc;
    plotting_vars.prim_isescaped[i_particle] = stops_beyond_ecal;
  }

  if (charge == 0) return false;

  double rad_curvature = 100*transverse_mom/(0.3*B_field); //p = 0.3*B*r where p in GeV/c, B in T, r in m (*100 to convert to cm)
  double theta_xT = atan2(p_x, transverse_mom); //helix pitch angle
  double pitch = std::abs(2*M_PI*rad_curvature*tan(theta_xT)); //distance between two turns of a helix in cm
  double tan_theta = tan(theta_xT);

  //Find centre of circular path
  double centre_circle_y;
  double centre_circle_z;
  double L_yz; //length of curved track in y-z plane
  double theta_spanned = 0;

  if (charge == 1) {
    centre_circle_y = ystart + (rad_curvature*p_z/transverse_mom); //Note plus sign here as cross product gives F in direction of ( pz j - py k) F= q v x B
    centre_circle_z = zstart - (rad_curvature*p_y/transverse_mom);
  }
  else {
    centre_circle_y = ystart - (rad_curvature*p_z/transverse_mom); //Note minus sign here as cross product gives F in direction of ( -pz j + py k)
    centre_circle_z = zstart + (rad_curvature*p_y/transverse_mom);
  }
  double theta_start = atan2(ystart - centre_circle_y, zstart - centre_circle_z);

  if (stops_in_tpc) { 
    theta_spanned = std::abs(xstart - _MCPEndX->at(i_particle))*2*M_PI/pitch; 
  }
  else { //Case where primary exits TPC

    // Find position where track leaves TPC. Intersection of two circles.
    double m_const = (TPC_centre_z - centre_circle_z)/(TPC_centre_y-centre_circle_y); //gradient of line between two intersection points
    double a_const = (TPCInstrumentedRadius*TPCInstrumentedRadius-rad_curvature*rad_curvature - (TPC_centre_y*TPC_centre_y - centre_circle_y*centre_circle_y)-(TPC_centre_z*TPC_centre_z - centre_circle_z*centre_circle_z))/(2*(centre_circle_y-TPC_centre_y));
    double quadraticformula_b = -(2*m_const*(a_const - TPC_centre_y) + 2*TPC_centre_z);
    double quadraticformula_a = m_const*m_const + 1;
    double quadraticformula_c = (a_const - TPC_centre_y)*(a_const - TPC_centre_y) + TPC_centre_z*TPC_centre_z - TPCInstrumentedRadius*TPCInstrumentedRadius;

    // If circles do intersect:
    if(quadraticformula_b*quadraticformula_b - 4*quadraticformula_a*quadraticformula_c > 0){
      double z_intersect_1 = (-quadraticformula_b + std::sqrt(quadraticformula_b*quadraticformula_b - 4*quadraticformula_a*quadraticformula_c))/(2*quadraticformula_a);
      double y_intersect_1 = -m_const*z_intersect_1 + a_const;
      double z_intersect_2 = (-quadraticformula_b - std::sqrt(quadraticformula_b*quadraticformula_b - 4*quadraticformula_a*quadraticformula_c))/(2*quadraticformula_a);
      double y_intersect_2 = -m_const*z_intersect_2 + a_const;

      //Find angle wrt y in yz plane where track starts and where it intersects TPC boundary
      double theta_1 = atan2(y_intersect_1 - centre_circle_y, z_intersect_1 - centre_circle_z);
      double theta_2 = atan2(y_intersect_2 - centre_circle_y, z_intersect_2 - centre_circle_z);

      double theta_diff_1, theta_diff_2;
      //Lorentz force - if positively charged, theta is increasing and vice versa
      if(charge == 1){ 
        theta_diff_1 = (theta_1 > theta_start) ? (theta_1 - theta_start) : (2*M_PI + theta_1 - theta_start);
        theta_diff_2 = (theta_2 > theta_start) ? (theta_2 - theta_start) : (2*M_PI + theta_2 - theta_start);
      }
      else {
        theta_diff_1 = (theta_1 < theta_start) ? (theta_start - theta_1) : (2*M_PI + theta_start - theta_1);
        theta_diff_2 = (theta_2 < theta_start) ? (theta_start - theta_2) : (2*M_PI + theta_start - theta_2);
      }

      if(theta_diff_1 < theta_diff_2){
        theta_spanned = theta_diff_1;
      }
      else {
        theta_spanned = theta_diff_2;
      }

      double x_end = (p_x > 0) ? (xstart + theta_spanned*pitch/(2*M_PI)) : (xstart - theta_spanned*pitch/(2*M_PI));

      //Check if escapes through end caps
      if (x_end > TPC_centre_x + TPCInstrumentedLength) {
        theta_spanned = (TPC_centre_x + TPCInstrumentedLength - xstart)*2*M_PI/pitch;
      }
      else if (x_end < TPC_centre_x - TPCInstrumentedLength) { 
        theta_spanned = (xstart - (TPC_centre_x - TPCInstrumentedLength))*2*M_PI/pitch;
      }
    }
    else { //Radius of curvature small enough to spiral within TPC (must leave through end caps)
      if (p_x > 0) theta_spanned = (TPC_centre_x + TPCInstrumentedLength - xstart)*2*M_PI/pitch;
      else theta_spanned = (xstart - (TPC_centre_x - TPCInstrumentedLength))*2*M_PI/pitch;
    }
  }
  double nturns = theta_spanned/(2*M_PI);
  L_yz = rad_curvature*theta_spanned;

  double length_track_x = (theta_spanned/(2*M_PI))*pitch;
  double bg = 0; 
  double gamma = 0;
  double beta = CalcBeta(mom_tot, bg, gamma, mass);

  double nhits = FindNHits(pixel_spacing_cm, centre_circle_y, centre_circle_z, rad_curvature, theta_start, theta_spanned, charge);
  if (nturns > 1) nhits *= nturns;

  double sigmax = (drift_velocity/100)*(1/(adc_sampling_frequency));
  double sigmax_frac = sigmax/(std::abs(length_track_x)/100);
  double sigmayz = (spatial_resolution/(1000)); //needs to be in m
  double momres_yz = transverse_mom*(std::sqrt(720/(nhits+4)) * (sigmayz*transverse_mom/(0.3*B_field*(L_yz/100)*(L_yz/100)))
    * std::sqrt(1.-(1./21.)*(L_yz/rad_curvature)*(L_yz/rad_curvature)));
  double momres_ms = transverse_mom*(0.016/(0.3*B_field*(L_yz/100)*cos(theta_xT)*beta))*std::sqrt(L_yz/X0);
  double momres_tottransverse = std::sqrt(momres_yz*momres_yz + momres_ms*momres_ms)/transverse_mom;
  double sigma_theta = (cos(theta_xT)*cos(theta_xT) * (pitch/(2*M_PI*rad_curvature)) *
    std::sqrt(sigmax_frac*sigmax_frac + momres_tottransverse*momres_tottransverse));
  double momres_frac = std::sqrt(momres_tottransverse*momres_tottransverse + (sigma_theta*tan_theta)*(sigma_theta*tan_theta));

  if (_MCPMotherTrkID->at(i_particle) == 0) {
    plotting_vars.prim_nturns[i_particle] = nturns*2*M_PI;
    plotting_vars.prim_nhits[i_particle] = nhits;
    plotting_vars.prim_tracklengthyz[i_particle] = L_yz;
    plotting_vars.prim_momresms[i_particle] = momres_ms/transverse_mom;
    plotting_vars.prim_momresyz[i_particle] = momres_yz/transverse_mom;
    plotting_vars.prim_momresx[i_particle] = std::abs(sigma_theta*tan_theta);
  }

  if(momres_frac*(mom_tot/energy)*(mom_tot/energy) > energy_resolution_threshold) return false;
  return true;
}

void SampleHandlerBeamNDGAr::clearBranchVectors() {
  _MCPStartX->clear();
  _MCPStartY->clear();
  _MCPStartZ->clear();
  _MCPEndX->clear();
  _MCPEndY->clear();
  _MCPEndZ->clear();
  _MCPStartPX->clear();
  _MCPStartPY->clear();
  _MCPStartPZ->clear();
  _MCPCalPX->clear();
  _MCPCalPY->clear();
  _MCPCalPZ->clear();
  _MCPPDG->clear();
  _MCPTrkID->clear();
  _MCPMotherTrkID->clear();
  _MCPEndProcess->clear();
  _TPCHitTrkID->clear();
  _TPCHitIsSec->clear();
  _TPCHitEnergy->clear();
  _TPCHitX->clear();
  _TPCHitY->clear();
  _TPCHitZ->clear();
  _CalHitTrkID->clear();
  _CalHitEnergy->clear();
  _CalHitIsSec->clear();
  _CalHitTime->clear();
  _CalHitX->clear();
  _CalHitY->clear();
  _CalHitZ->clear();
}

// FastGArSim output has beam x, Genie has beam z but left handed coordinate system. This code uses beam z but right handed coordinate system
void SampleHandlerBeamNDGAr::fixCoordinates() {
  auto switchCoords = [&](float& xcoord, float& zcoord) {
    float tmp = xcoord;
    xcoord = -zcoord;
    zcoord = tmp;
  };

  for (size_t i = 0; i < _MCPStartX->size(); i++) {
    switchCoords(_MCPStartX->at(i), _MCPStartZ->at(i));
    switchCoords(_MCPEndX->at(i), _MCPEndZ->at(i));
    switchCoords(_MCPStartPX->at(i), _MCPStartPZ->at(i));
    switchCoords(_MCPCalPX->at(i), _MCPCalPZ->at(i));
  }
  for (size_t i = 0; i < _TPCHitX->size(); i++) {
    switchCoords(_TPCHitX->at(i), _TPCHitZ->at(i));
  }
  for (size_t i = 0; i < _CalHitX->size(); i++) {
    switchCoords(_CalHitX->at(i), _CalHitZ->at(i));
  }
  _PXnu = -_PXnu;
  _PXlep = -_PXlep;
}

void SampleHandlerBeamNDGAr::FillGeoVars() {
  if (B_field != _BField) {
    MACH3LOG_ERROR("B-field specified in config ({} T) does not match the input file ({} T)", B_field, _BField);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (TPCInstrumentedRadius != _TPCRad) {
    MACH3LOG_ERROR("TPC radius specified in config does not match the input file");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (TPCInstrumentedLength != _TPCLen/2.) {
    MACH3LOG_ERROR("TPC length specified in config {} does not match the input file {}", TPCInstrumentedLength, _TPCLen/2.);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  ECALInnerRadius = _TPCRad + _BarrelGap;
  ECALOuterFrontRadius = ECALInnerRadius + ECALBarrelForwardDepth*ECALSciX0;
  ECALOuterBackRadius = ECALInnerRadius + ECALBarrelBackwardDepth*ECALSciX0;
  double maxECALOuterRadius = ECALInnerRadius + (_NBarrelHG*(_HGAbsWidth+_HGSciWidth+_HGBoardWidth) + _NBarrelLG*(_LGAbsWidth+_LGSciWidth));
  ECALEndCapStart = _TPCLen/2. + _EndCapGap + 0.5; // 0.5cm readout PCB
  ECALEndCapEnd = ECALEndCapStart + ECALEndCapDepth*ECALSciX0;
  double maxECALEndCapEnd   = ECALEndCapStart + (_NEndCapHG*(_HGAbsWidth+_HGSciWidth+_HGBoardWidth) + _NEndCapLG*(_LGAbsWidth+_LGSciWidth));

  // Check variable ecal parameters are all reasonable
  if ((nECALBackSegments%2 != 1 && nECALBackSegments != 0) || nECALBackSegments < 0 || nECALBackSegments >= _NumCalSides) {
    MACH3LOG_ERROR("Ensure that nECALBackSegments is set to an odd positive integer less than {}, or 0.", _NumCalSides);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (ECALBarrelBackwardDepth > ECALBarrelForwardDepth && nECALBackSegments > 0) {
    MACH3LOG_ERROR("Ensure that ECALBarrelBackwardDepth <= ECALBarrelForwardDepth.");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (ECALOuterFrontRadius > maxECALOuterRadius) {
    MACH3LOG_ERROR("ECAL radius in config ({} cm) cannot be larger than radius in sample ({} cm).", ECALOuterFrontRadius, maxECALOuterRadius);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (ECALEndCapEnd > maxECALEndCapEnd) {
    MACH3LOG_ERROR("ECAL end cap depth in config ({} cm) cannot be larger than depth in sample ({} cm).", ECALEndCapEnd, maxECALEndCapEnd);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  ECalBackSegments = {};
  int first_seg = _NumCalSides/2 - (nECALBackSegments-1)/2;
  for (int iSeg=0; iSeg<nECALBackSegments; iSeg++) {
    ECalBackSegments.push_back(first_seg + iSeg);
  }

  outerECalA = {};
  outerECalP = {};

  for (int iPlane=0; iPlane<_NumCalSides; iPlane++) {
    double outer_radius = ECALOuterFrontRadius;
    if (std::find(ECalBackSegments.begin(), ECalBackSegments.end(), iPlane) != ECalBackSegments.end()) {
      outer_radius = ECALOuterBackRadius;
    }
    double y = outer_radius*sin(2. * M_PI * iPlane / _NumCalSides);
    double z = outer_radius*cos(2. * M_PI * iPlane / _NumCalSides);
    outerECalP.push_back({0, y, z});
    outerECalA.push_back({0, y/outer_radius, z/outer_radius});
  }
  outerECalP.push_back({ECALEndCapEnd, 0., 0.});
  outerECalP.push_back({-ECALEndCapEnd, 0., 0.});
  outerECalA.push_back({1., 0., 0.});
  outerECalA.push_back({-1., 0., 0.});

  if (nECALBackSegments > 0) {
    dividingPlaneP = {{0., 0., 0.}, {0., 0., 0.}};
    double theta_back1 = (2*M_PI/_NumCalSides)*ECalBackSegments[0] - M_PI/_NumCalSides - M_PI/2;
    double theta_back2 = (2*M_PI/_NumCalSides)*ECalBackSegments.back() - M_PI/_NumCalSides - M_PI/2;
    dividingPlaneA = {{0., sin(theta_back1), cos(theta_back1)}, {0., sin(theta_back2), cos(theta_back2)}};
  }
}

int SampleHandlerBeamNDGAr::SetupExperimentMC() {

  MACH3LOG_INFO("-------------------------------------------------------------------");

  TChain* _data = new TChain("AnaTree");
  TChain* _geometry = new TChain("GeoTree");
  TChain* _genie = new TChain("gst");
  // Maps the file index within the TChain (GetTreeNumber()) to its sample index.
  std::vector<size_t> fileIndexToSample;
  for (size_t iSample=0;iSample<SampleDetails.size();iSample++) {
    for (const std::string& filename : SampleDetails[iSample].mc_files) {
      MACH3LOG_INFO("Adding FastGArSim file to TChain: {}", filename);
      // HH: Check whether the file exists, see https://root.cern/doc/master/classTChain.html#a78a896924ac6c7d3691b7e013bcbfb1c
      int _add_rtn = _data->Add(filename.c_str(), -1);
      if(_add_rtn == 0){
        MACH3LOG_ERROR("Could not add file {} to TChain, please check the file exists and is readable", filename);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      int _add_rtn_geo = _geometry->Add(filename.c_str(), -1);
      if(_add_rtn_geo == 0){
        MACH3LOG_ERROR("Could not add file {} to TChain, please check the file exists and is readable", filename);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      // Read Genie file
      std::string genieFileStr = "Inputs/DUNE_NDGAr_files/FastGArSim/genie_inputs/";
      if (interaction_model == "hA") genieFileStr += "numu_argon_G18_10a_gst.root";
      else if (interaction_model == "hN") genieFileStr += "numu_argon_G18_10b_gst.root";
      else if (interaction_model == "INCL") genieFileStr += "numu_argon_G18_10c_gst.root";
      else if (interaction_model == "G4BC") genieFileStr += "numu_argon_G18_10d_gst.root";
      else {
        MACH3LOG_ERROR("{} is not an availble interaction model.", interaction_model);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      MACH3LOG_INFO("Adding genie file to TChain: {}", genieFileStr);

      int _add_rtn_genie = _genie->Add(genieFileStr.c_str(), -1);
      if(_add_rtn_genie == 0){
        MACH3LOG_ERROR("Could not add file {} to TChain, please check the file exists and is readable", genieFileStr);
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      // Each file added (glob patterns may expand to multiple) gets the same sample index.
      // We query how many files are now in the chain to know how many entries to push.
      const int nFilesNow = _data->GetListOfFiles()->GetEntries();
      while(static_cast<int>(fileIndexToSample.size()) < nFilesNow){
        fileIndexToSample.push_back(iSample);
      }
    }
  }

  _data->SetBranchStatus("*", 0);
  _genie->SetBranchStatus("*", 0);

  auto readBranch = [&](TTree* tree, const char* name, void* addr) {
    tree->SetBranchStatus(name, 1);
    tree->SetBranchAddress(name, addr);
  };

  readBranch(_data, "eventID", &_EventID);
  readBranch(_data, "startX", &_MCPStartX);
  readBranch(_data, "startY", &_MCPStartY);
  readBranch(_data, "startZ", &_MCPStartZ);
  readBranch(_data, "endX", &_MCPEndX);
  readBranch(_data, "endY", &_MCPEndY);
  readBranch(_data, "endZ", &_MCPEndZ);
  readBranch(_data, "startPX", &_MCPStartPX);
  readBranch(_data, "startPY", &_MCPStartPY);
  readBranch(_data, "startPZ", &_MCPStartPZ);
  readBranch(_data, "ecalPX", &_MCPCalPX);
  readBranch(_data, "ecalPY", &_MCPCalPY);
  readBranch(_data, "ecalPZ", &_MCPCalPZ);
  readBranch(_data, "pdgCode", &_MCPPDG);
  readBranch(_data, "trackID", &_MCPTrkID);
  readBranch(_data, "motherID", &_MCPMotherTrkID);
  readBranch(_data, "endProcess", &_MCPEndProcess);
  readBranch(_data, "tpcHitTrackID", &_TPCHitTrkID);
  readBranch(_data, "tpcHitEdep", &_TPCHitEnergy);
  readBranch(_data, "tpcHitX", &_TPCHitX);
  readBranch(_data, "tpcHitY", &_TPCHitY);
  readBranch(_data, "tpcHitZ", &_TPCHitZ);
  readBranch(_data, "tpcHitIsSec", &_TPCHitIsSec);
  readBranch(_data, "ecalHitTrackID", &_CalHitTrkID);
  readBranch(_data, "ecalHitEdep", &_CalHitEnergy);
  readBranch(_data, "ecalHitIsSec", &_CalHitIsSec);
  readBranch(_data, "ecalHitTime", &_CalHitTime);
  readBranch(_data, "ecalHitX", &_CalHitX);
  readBranch(_data, "ecalHitY", &_CalHitY);
  readBranch(_data, "ecalHitZ", &_CalHitZ);

  readBranch(_geometry, "gar_tpc_radius", &_TPCRad);
  readBranch(_geometry, "gar_tpc_length", &_TPCLen);
  readBranch(_geometry, "gar_magnetic_field", &_BField);
  readBranch(_geometry, "ecal_num_sides", &_NumCalSides);
  readBranch(_geometry, "ecal_barrel_gap", &_BarrelGap);
  readBranch(_geometry, "ecal_endcap_gap", &_EndCapGap);
  readBranch(_geometry, "ecal_hg_absorber_thickness", &_HGAbsWidth);
  readBranch(_geometry, "ecal_lg_absorber_thickness", &_LGAbsWidth);
  readBranch(_geometry, "ecal_hg_scintillator_thickness", &_HGSciWidth);
  readBranch(_geometry, "ecal_lg_scintillator_thickness", &_LGSciWidth);
  readBranch(_geometry, "ecal_hg_board_thickness", &_HGBoardWidth);
  readBranch(_geometry, "ecal_barrel_hg_layers", &_NBarrelHG);
  readBranch(_geometry, "ecal_barrel_lg_layers", &_NBarrelLG);
  readBranch(_geometry, "ecal_endcap_hg_layers", &_NEndCapHG);
  readBranch(_geometry, "ecal_endcap_lg_layers", &_NEndCapLG);
  _geometry->GetEntry(0);
  FillGeoVars();

  readBranch(_genie, "Ev", &_Enu);
  readBranch(_genie, "pxv", &_PXnu);
  readBranch(_genie, "pyv", &_PYnu);
  readBranch(_genie, "pzv", &_PZnu);
  readBranch(_genie, "neu", &_nuPDG);
  readBranch(_genie, "El", &_Elep);
  readBranch(_genie, "pxl", &_PXlep);
  readBranch(_genie, "pyl", &_PYlep);
  readBranch(_genie, "pzl", &_PZlep);
  readBranch(_genie, "cc", &_isCC);
  readBranch(_genie, "nfpip", &_npip);
  readBranch(_genie, "nfpim", &_npim);
  readBranch(_genie, "nfpi0", &_npi0);
  readBranch(_genie, "neut_code", &_neut_code);

  size_t nEntries = static_cast<size_t>(downsampling*static_cast<double>(_data->GetEntries()));
  size_t countwidth = nEntries / 50;

  dunendgarmcFitting.resize(nEntries);
  dunendgarmcPlotting.resize(nEntries);
  _data->GetEntry(0);
  _genie->GetEntry(0);

  double pixel_spacing_cm = pixel_spacing/10; //convert to cm
  int numCC = 0;
  int num_in_fdv = 0;

  for (unsigned int i_event = 0; i_event < nEntries; ++i_event) { 
    if (i_event != 0) clearBranchVectors();
    _data->GetEntry(i_event);
    _genie->GetEntry(_EventID);

    // FastGArSim output has a different coordinate system to the genie file and this code. This is fixed here.
    fixCoordinates();

    if (i_event % countwidth == 0) {
      MaCh3Utils::PrintProgressBar(i_event, static_cast<Long64_t>(nEntries));
    }

    const Int_t treeNum = _data->GetTreeNumber();
    if(treeNum < 0 || static_cast<size_t>(treeNum) >= fileIndexToSample.size()){
      MACH3LOG_ERROR("GetTreeNumber() returned {} which is out of range [0, {})", treeNum, fileIndexToSample.size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    const size_t sample_index = fileIndexToSample[static_cast<size_t>(treeNum)];

    dunendgarmcFitting[i_event].SampleIndex = static_cast<unsigned int>(sample_index);
    dunendgarmcFitting[i_event].rw_etru = _Enu;
    dunendgarmcFitting[i_event].rw_isCC = _isCC;
    dunendgarmcFitting[i_event].nupdg = _nuPDG;
    dunendgarmcFitting[i_event].nupdgUnosc = _nuPDG;
    dunendgarmcFitting[i_event].OscChannelIndex = static_cast<double>(GetOscChannel(SampleDetails[sample_index].OscChannels, dunendgarmcFitting[i_event].nupdgUnosc, dunendgarmcFitting[i_event].nupdg));
    dunendgarmcFitting[i_event].rw_berpaacvwgt = _BeRPA_cvwgt;
    dunendgarmcFitting[i_event].Target = 40; // Assume everything is Argon
    dunendgarmcFitting[i_event].rw_Q0 = _Enu - _Elep;
    dunendgarmcFitting[i_event].rw_Q3 = std::sqrt((_PXnu-_PXlep)*(_PXnu-_PXlep) + (_PYnu-_PYlep)*(_PYnu-_PYlep) + (_PZnu-_PZlep)*(_PZnu-_PZlep));
    dunendgarmcFitting[i_event].norm_s = 1.;
    dunendgarmcFitting[i_event].pot_s = beamNDGArSampleDetails[sample_index].pot/(downsampling*1e21);
    dunendgarmcFitting[i_event].flux_w = 1.;

    int M3Mode = Modes->GetModeFromGenerator(std::abs(_neut_code));
    if (!_isCC) M3Mode += 14; //Account for no ability to distinguish CC/NC
    if (M3Mode > 15) M3Mode -= 1; //Account for no NCSingleKaon
    dunendgarmcFitting[i_event].mode = M3Mode;

    std::vector<double> vertex = {M3::_BAD_DOUBLE_, M3::_BAD_DOUBLE_, M3::_BAD_DOUBLE_};

    std::unordered_map<int, std::vector<int>> mother_to_daughter_ID;
    std::unordered_map<int, size_t> ID_to_index;
    // std::unordered_map<int, std::vector<double>> ID_to_ECalDep; // particle track ID -> total energy deposited in each ecal layer
    std::unordered_map<int, double> eID_to_smallest_t;
    std::unordered_map<int, std::vector<double>> eID_to_showerstart; // electron track ID -> shower start position (x, y, z)
    std::unordered_map<int, double> ID_to_TPCDep;
    size_t n_particles = _MCPTrkID->size();
    dunendgarmcPlotting[i_event].rw_ePi0 = 0.; 
    dunendgarmcPlotting[i_event].npi0 = 0; 

    // const int tot_ecal_layers = std::max(_NBarrelHG+_NBarrelHG, _NEndCapHG+_NEndCapLG);

    // Fill maps
    for (size_t i_particle=0; i_particle<n_particles; i_particle++) {
      int trkID = _MCPTrkID->at(i_particle);
      int motherID = _MCPMotherTrkID->at(i_particle);
      // ID_to_ECalDep[trkID] = std::vector<double>(static_cast<size_t>(tot_ecal_layers),0.);
      ID_to_index[trkID] = i_particle;
      mother_to_daughter_ID[motherID].push_back(trkID);
      mother_to_daughter_ID[trkID]; //Ensure all particles are added to the map (even if no secondaries)

      // Fill particle-level variables
      if (_MCPPDG->at(i_particle) == 22) {
        double photon_energy = std::sqrt(_MCPStartPX->at(i_particle)*_MCPStartPX->at(i_particle) + _MCPStartPY->at(i_particle)*_MCPStartPY->at(i_particle) + _MCPStartPZ->at(i_particle)*_MCPStartPZ->at(i_particle))/1000.;
        dunendgarmcPlotting[i_event].photon_energy.push_back(photon_energy);
      }
      else if (std::abs(_MCPPDG->at(i_particle)) == 11) {
        eID_to_smallest_t[trkID] = std::numeric_limits<double>::max();
      }
    }

    // Fill map from particle ID to ECal deposited energy
    for (size_t i_calhit=0; i_calhit<_CalHitTrkID->size(); i_calhit++) {
      int dep_trkid = _CalHitTrkID->at(i_calhit);
      if (dep_trkid <= 0) continue;
      double dep_energy = _CalHitEnergy->at(i_calhit)/1000.;
      // double dep_depth = GetCalDepth(_CalHitX->at(i_calhit), _CalHitY->at(i_calhit), _CalHitZ->at(i_calhit));
      // ID_to_ECalDep[dep_trkid][static_cast<size_t>(dep_layer)] += dep_energy;

      // Also store the position of first hit above energy threshold for electrons/positrons
      bool dep_issec = _CalHitIsSec->at(i_calhit);
      double energy_threshold = 0.0005;
      if ((dep_energy > energy_threshold) && (!dep_issec) && (std::abs(_MCPPDG->at(ID_to_index.at(dep_trkid))) == 11)) {
        double time = _CalHitTime->at(i_calhit);
        if (time < eID_to_smallest_t.at(dep_trkid)) {
          eID_to_smallest_t[dep_trkid] = time;
          eID_to_showerstart[dep_trkid] = {_CalHitX->at(i_calhit), _CalHitY->at(i_calhit), _CalHitZ->at(i_calhit)};
        }
      }
    }

    // Fill map from particle ID to TPC deposited energy
    for (size_t i_tpchit=0; i_tpchit<_TPCHitTrkID->size(); i_tpchit++) {
      int trkid = _TPCHitTrkID->at(i_tpchit);
      if (trkid <= 0) continue;
      double dep_energy = _TPCHitEnergy->at(i_tpchit)/1000.;
      if (!_TPCHitIsSec->at(i_tpchit)) ID_to_TPCDep[trkid] += dep_energy;
    }

    bool isEventAccepted = true;
    double pi0_p2 = 0.;
    // double mu_p2 = 0.;

    // Resize vectors for prim-level parameters
    size_t n_prim_in_event = mother_to_daughter_ID[0].size();
    dunendgarmcPlotting[i_event].prim_pdg.resize(n_prim_in_event, M3::_BAD_INT_);
    dunendgarmcPlotting[i_event].prim_evis.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_momentum.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_transversemomentum.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_bangle.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_beamangle.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_isaccepted.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].prim_iscurvatureresolved.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].prim_isstoppedintpc.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].prim_isstoppedinecal.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].prim_isstoppedingap.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].prim_isstoppedinbarrelgap.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].prim_isstoppedinendgap.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].prim_isstoppedinbarrel.resize(n_prim_in_event, M3::_BAD_INT_);
    dunendgarmcPlotting[i_event].prim_isstoppedinendcap.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].prim_isescaped.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].prim_startx.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_startr2.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_endr.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_enddepth.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_endx.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_endy.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_endz.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_nturns.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_nhits.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_tracklengthyz.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_momresms.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_momresyz.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_momresx.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_tpcedepfrac.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].prim_iscontained.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 

    // Loop through primaries
    for (int& primID : mother_to_daughter_ID[0]) {

      // Do not require the reconstruction of neutrons and neutrinos
      size_t prim_index = ID_to_index[primID];
      int pdg = _MCPPDG->at(prim_index);
      if (pdg == 2112 || std::abs(pdg) == 12 || std::abs(pdg) == 14 || std::abs(pdg) == 16) continue;

      dunendgarmcPlotting[i_event].prim_pdg[prim_index] = pdg; 

      // Check if primary is resolved from curvature
      bool isCurvatureResolved = false;
      // Remove all descendants of primID from mother_to_daughter_ID whose momentum we get from curvature (or whose parent we get from curvature)
      if (CurvatureResolutionFilter(primID, mother_to_daughter_ID, ID_to_index, dunendgarmcPlotting[i_event], pixel_spacing_cm)) {
        isCurvatureResolved = true;
      }
      dunendgarmcPlotting[i_event].prim_iscurvatureresolved[prim_index] = isCurvatureResolved;

      // Check for containment
      bool isContained = IsPrimContained(primID, mother_to_daughter_ID, ID_to_index, eID_to_showerstart, dunendgarmcPlotting[i_event]);

      // Primary is accepted if contained or curvature resolved
      bool isPrimAccepted = true;
      if (!(isContained || isCurvatureResolved)) {
        isPrimAccepted = false;
        isEventAccepted = false;
      }
      dunendgarmcPlotting[i_event].prim_isaccepted[prim_index] = isPrimAccepted;
      dunendgarmcPlotting[i_event].prim_iscontained[prim_index] = isContained;

      double p_x = _MCPStartPX->at(prim_index)/1000.;
      double p_y = _MCPStartPY->at(prim_index)/1000.;
      double p_z = _MCPStartPZ->at(prim_index)/1000.;
      double p2 = p_x*p_x + p_y*p_y + p_z*p_z;

      if (!isContained && p2 < 0.0001) {
        MACH3LOG_WARN("Particle of pdg {} with momentum {} is uncontained", pdg, std::sqrt(p2));
      }
      if (pdg < 1000000000) {
        double mass = MaCh3Utils::GetMassFromPDG(pdg);
        double energy = std::sqrt(p2+mass*mass);
        dunendgarmcPlotting[i_event].prim_tpcedepfrac[prim_index] = ID_to_TPCDep[primID]/(energy-mass);
      }

      // Get vertex from primary muon start position
      if (pdg == 13 && std::abs((p_x-_PXlep)/_PXlep) < 0.00001 && std::abs((p_y-_PYlep)/_PYlep) < 0.00001 && std::abs((p_z-_PZlep)/_PZlep) < 0.00001) {
        vertex = {_MCPStartX->at(prim_index), _MCPStartY->at(prim_index), _MCPStartZ->at(prim_index)};
        dunendgarmcPlotting[i_event].lep_tracklengthyz = dunendgarmcPlotting[i_event].prim_tracklengthyz[prim_index];
      }

      // Get pi0 information 
      if(pdg == 111) {
        dunendgarmcPlotting[i_event].npi0 ++;
        if (p2 > pi0_p2) {
          pi0_p2 = p2;
          double mass = MaCh3Utils::GetMassFromPDG(pdg);
          dunendgarmcPlotting[i_event].rw_ePi0 = std::sqrt(p2+mass*mass);
        }
      }
    }
    if (vertex[0] == M3::_BAD_DOUBLE_) MACH3LOG_ERROR("No vertex found for event {}.", i_event);

    dunendgarmcFitting[i_event].rw_lep_pX = _PXlep;
    dunendgarmcFitting[i_event].rw_lep_pY = _PYlep;
    dunendgarmcFitting[i_event].rw_lep_pZ = _PZlep;
    dunendgarmcFitting[i_event].rw_LepE = _Elep;
    dunendgarmcPlotting[i_event].is_accepted = isEventAccepted;
    dunendgarmcFitting[i_event].rw_vtx_x = vertex[0]-TPC_centre_x;
    dunendgarmcFitting[i_event].rw_vtx_y = vertex[1]-TPC_centre_y;
    dunendgarmcFitting[i_event].rw_vtx_z = vertex[2]-TPC_centre_z;

    // Find lepton kinematic variables
    double lep_momentum = std::sqrt(_PXlep*_PXlep + _PYlep*_PYlep + _PZlep*_PZlep);
    double lep_pBeam = (_PYlep*BeamDirection[1] + _PZlep*BeamDirection[2]);
    double lep_pB = _PXlep;
    double lep_pPerp = (_PYlep*BeamDirection[2] - _PZlep*BeamDirection[1]);
    double lep_beamangle = acos(lep_pBeam/lep_momentum)*180/M_PI; //Angle to beam 
    double lep_bangle = acos(lep_pB/lep_momentum)*180/M_PI; //Angle to B-field
    double lep_perpangle = acos(lep_pPerp/lep_momentum)*180/M_PI; //Angle to axis perpendicular to beam and B
    double lep_phi = atan2(lep_pPerp, lep_pB)*180/M_PI;

    dunendgarmcPlotting[i_event].rw_lep_theta = lep_beamangle;
    dunendgarmcPlotting[i_event].rw_lep_phi = lep_phi;
    dunendgarmcPlotting[i_event].rw_lep_bangle = lep_bangle;
    dunendgarmcPlotting[i_event].rw_lep_p = lep_momentum;
    dunendgarmcFitting[i_event].rw_lep_pT = std::sqrt(lep_momentum*lep_momentum - lep_pBeam*lep_pBeam); 

    double radius = std::sqrt((vertex[1]-TPC_centre_y)*(vertex[1]-TPC_centre_y) 
                              + (vertex[2]-TPC_centre_z)*(vertex[2]-TPC_centre_z)); //find radius of interaction vertex
    dunendgarmcFitting[i_event].rw_rad = radius;

    if(std::abs(vertex[0] - TPC_centre_x) <= TPCFidLength && radius<=TPCFidRadius){
      num_in_fdv++;
      dunendgarmcPlotting[i_event].in_fdv = 1;
    } else{
      dunendgarmcPlotting[i_event].in_fdv = 0;
    }
    if(_isCC) numCC++;

    // Perform 'geometric correction' if do_geometric_correction set to true
    dunendgarmcPlotting[i_event].geometric_correction = 1.;
    if (do_geometric_correction) {
      if ((lep_bangle < 45 || lep_bangle > 135) && lep_momentum > 0.3) dunendgarmcPlotting[i_event].geometric_correction = 0.;
      else if ((lep_perpangle < 45 || lep_perpangle > 135) && lep_momentum > 0.3) dunendgarmcPlotting[i_event].geometric_correction = 2.;
    }

  }
  MACH3LOG_INFO("nEntries = {}, numCC = {}, numFDV = {}", nEntries, numCC, num_in_fdv);

  _data->Reset();
  delete _data;
  _geometry->Reset();
  delete _geometry;
  _genie->Reset();
  delete _genie;

  return static_cast<int>(nEntries);
}

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wswitch-enum"
const double* SampleHandlerBeamNDGAr::GetPointerToKinematicParameter(KinematicTypes KinematicParameter, size_t iEvent) {
  switch(KinematicParameter) {
    case kTrueNeutrinoEnergy:
      return &dunendgarmcFitting[iEvent].rw_etru; 
    case kMode:
      return &dunendgarmcFitting[iEvent].mode;
    case kOscChannel:
      return &dunendgarmcFitting[iEvent].OscChannelIndex;
    case kTrueXPos:
      return &dunendgarmcFitting[iEvent].rw_vtx_x;
    case kTrueYPos:
      return &dunendgarmcFitting[iEvent].rw_vtx_y;
    case kTrueZPos:
      return &dunendgarmcFitting[iEvent].rw_vtx_z;
    case kTrueRad:
      return &dunendgarmcFitting[iEvent].rw_rad;
    case kTrueLepEnergy:
      return &dunendgarmcFitting[iEvent].rw_LepE;
    case kLepPT:
      return &dunendgarmcFitting[iEvent].rw_lep_pT;
    case kLepPZ:
      return &dunendgarmcFitting[iEvent].rw_lep_pZ;
    case kLepTheta:
      return &dunendgarmcPlotting[iEvent].rw_lep_theta;
    case kLepPhi:
      return &dunendgarmcPlotting[iEvent].rw_lep_phi;
    case kLepBAngle:
      return &dunendgarmcPlotting[iEvent].rw_lep_bangle;
    case kLepP:
      return &dunendgarmcPlotting[iEvent].rw_lep_p;
    case kLepTrackLengthYZ:
      return &dunendgarmcPlotting[iEvent].lep_tracklengthyz;
    case kTrueQ0:
      return &dunendgarmcFitting[iEvent].rw_Q0;
    case kTrueQ3:
      return &dunendgarmcFitting[iEvent].rw_Q3;
    case kEPi0:
      return &dunendgarmcPlotting[iEvent].rw_ePi0;
    default:
      MACH3LOG_ERROR("Did not recognise Kinematic Parameter {}", static_cast<int>(KinematicParameter));
      throw MaCh3Exception(__FILE__, __LINE__);
      return nullptr;
  }
}
#pragma GCC diagnostic pop

const double* SampleHandlerBeamNDGAr::GetPointerToKinematicParameter(double KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(std::round(KinematicVariable));
  return GetPointerToKinematicParameter(KinPar,static_cast<size_t>(iEvent));
}

const double* SampleHandlerBeamNDGAr::GetPointerToKinematicParameter(std::string KinematicParameter, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return GetPointerToKinematicParameter(KinPar,static_cast<size_t>(iEvent));
}

double SampleHandlerBeamNDGAr::ReturnKinematicParameter(int KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(std::round(KinematicVariable));
  return ReturnKinematicParameter(KinPar,static_cast<size_t>(iEvent));
}

double SampleHandlerBeamNDGAr::ReturnKinematicParameter(std::string KinematicParameter, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return ReturnKinematicParameter(KinPar,static_cast<size_t>(iEvent));
}

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wswitch-enum"
double SampleHandlerBeamNDGAr::ReturnKinematicParameter(KinematicTypes KinPar, size_t iEvent) {
  //HH: Special cases for dealing with non-doubles
  switch(KinPar) {
    case kEvent_IsAccepted:
      return static_cast<double>(dunendgarmcPlotting[iEvent].is_accepted);
    case kIsCC:
      return static_cast<double>(dunendgarmcFitting[iEvent].rw_isCC);
    case kInFDV:
      return static_cast<double>(dunendgarmcPlotting[iEvent].in_fdv);
    case kNPi0:
      return static_cast<double>(dunendgarmcPlotting[iEvent].npi0);
    default:
      return *GetPointerToKinematicParameter(KinPar, iEvent);
  }
}
#pragma GCC diagnostic pop

std::vector<double> SampleHandlerBeamNDGAr::ReturnKinematicVector(KinematicVecs KinVec, size_t iEvent) {
  switch(KinVec) {
    case kPrim_IsAccepted:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].prim_isaccepted.begin(),
        dunendgarmcPlotting[iEvent].prim_isaccepted.end());
    case kPrim_IsCurvatureResolved:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].prim_iscurvatureresolved.begin(),
        dunendgarmcPlotting[iEvent].prim_iscurvatureresolved.end());
    case kPrim_PDG:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].prim_pdg.begin(),
        dunendgarmcPlotting[iEvent].prim_pdg.end());
    case kPrim_IsStoppedInTPC:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].prim_isstoppedintpc.begin(),
        dunendgarmcPlotting[iEvent].prim_isstoppedintpc.end());
    case kPrim_IsStoppedInECal:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].prim_isstoppedinecal.begin(),
        dunendgarmcPlotting[iEvent].prim_isstoppedinecal.end());
    case kPrim_IsStoppedInBarrel:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].prim_isstoppedinbarrel.begin(),
        dunendgarmcPlotting[iEvent].prim_isstoppedinbarrel.end());
    case kPrim_IsStoppedInEndCap:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].prim_isstoppedinendcap.begin(),
        dunendgarmcPlotting[iEvent].prim_isstoppedinendcap.end());
    case kPrim_IsStoppedInGap:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].prim_isstoppedingap.begin(),
        dunendgarmcPlotting[iEvent].prim_isstoppedingap.end());
    case kPrim_IsStoppedInEndGap:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].prim_isstoppedinendgap.begin(),
        dunendgarmcPlotting[iEvent].prim_isstoppedinendgap.end());
    case kPrim_IsStoppedInBarrelGap:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].prim_isstoppedinbarrelgap.begin(),
        dunendgarmcPlotting[iEvent].prim_isstoppedinbarrelgap.end());
    case kPrim_IsEscaped:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].prim_isescaped.begin(),
        dunendgarmcPlotting[iEvent].prim_isescaped.end());
    case kPrim_EVis:
      return dunendgarmcPlotting[iEvent].prim_evis;
    case kPrim_Momentum:
      return dunendgarmcPlotting[iEvent].prim_momentum;
    case kPrim_TransverseMomentum:
      return dunendgarmcPlotting[iEvent].prim_transversemomentum;
    case kPrim_BAngle:
      return dunendgarmcPlotting[iEvent].prim_bangle;
    case kPrim_BeamAngle:
      return dunendgarmcPlotting[iEvent].prim_beamangle;
    case kPrim_NTurns:
      return dunendgarmcPlotting[iEvent].prim_nturns;
    case kPrim_NHits:
      return dunendgarmcPlotting[iEvent].prim_nhits;
    case kPrim_TrackLengthYZ:
      return dunendgarmcPlotting[iEvent].prim_tracklengthyz;
    case kPrim_MomResMS:
      return dunendgarmcPlotting[iEvent].prim_momresms;
    case kPrim_MomResYZ:
      return dunendgarmcPlotting[iEvent].prim_momresyz;
    case kPrim_MomResX:
      return dunendgarmcPlotting[iEvent].prim_momresx;
    case kPrim_StartR2:
      return dunendgarmcPlotting[iEvent].prim_startr2;
    case kPrim_EndR:
      return dunendgarmcPlotting[iEvent].prim_endr;
    case kPrim_EndDepth:
      return dunendgarmcPlotting[iEvent].prim_enddepth;
    case kPrim_EndX:
      return dunendgarmcPlotting[iEvent].prim_endx;
    case kPrim_EndY:
      return dunendgarmcPlotting[iEvent].prim_endy;
    case kPrim_EndZ:
      return dunendgarmcPlotting[iEvent].prim_endz;
    case kPrim_StartX:
      return dunendgarmcPlotting[iEvent].prim_startx;
    case kPrim_TPCEDepFrac:
      return dunendgarmcPlotting[iEvent].prim_tpcedepfrac;
    case kPrim_IsContained:
      return dunendgarmcPlotting[iEvent].prim_iscontained;
    case kShower_CosNorm:
      return dunendgarmcPlotting[iEvent].shower_cosnorm;
    case kShower_PDG:
      return dunendgarmcPlotting[iEvent].shower_pdg;
    case kShower_DCalBoundary:
      return dunendgarmcPlotting[iEvent].shower_dcalboundary;
    case kShower_Energy:
      return dunendgarmcPlotting[iEvent].shower_energy;
    case kShower_BAngle:
      return dunendgarmcPlotting[iEvent].shower_bangle;
    case kShower_IsContained:
      return dunendgarmcPlotting[iEvent].shower_iscontained;
    case kShower_IsConv:
      return dunendgarmcPlotting[iEvent].shower_isconv;
    case kPhoton_Energy:
      return dunendgarmcPlotting[iEvent].photon_energy;
    default:
      MACH3LOG_ERROR("Unrecognized Kinematic Vector: {}", static_cast<int>(KinVec));
      throw MaCh3Exception(__FILE__, __LINE__);
  }
}

std::vector<double> SampleHandlerBeamNDGAr::ReturnKinematicVector(int KinematicVector, int iEvent) {
  KinematicVecs KinVec = static_cast<KinematicVecs>(KinematicVector);
  return ReturnKinematicVector(KinVec, static_cast<size_t>(iEvent));
}

std::vector<double> SampleHandlerBeamNDGAr::ReturnKinematicVector(std::string KinematicVector, int iEvent) {
  KinematicVecs KinVec = static_cast<KinematicVecs>(ReturnKinematicVectorFromString(KinematicVector));
  return ReturnKinematicVector(KinVec, static_cast<size_t>(iEvent));
}

void SampleHandlerBeamNDGAr::SetupFDMC() {
  for(size_t iEvent = 0 ;iEvent < GetNEvents() ; ++iEvent){
    MCSamples[iEvent].rw_etru = &(dunendgarmcFitting[iEvent].rw_etru);
    MCSamples[iEvent].mode = &(dunendgarmcFitting[iEvent].mode);
    MCSamples[iEvent].Target = &(dunendgarmcFitting[iEvent].Target);    
    MCSamples[iEvent].isNC = !dunendgarmcFitting[iEvent].rw_isCC;
    MCSamples[iEvent].nupdg = &(dunendgarmcFitting[iEvent].nupdg);
    MCSamples[iEvent].nupdgUnosc = &(dunendgarmcFitting[iEvent].nupdgUnosc);
    MCSamples[iEvent].NominalSample = static_cast<int>(dunendgarmcFitting[iEvent].SampleIndex);
  }
}
