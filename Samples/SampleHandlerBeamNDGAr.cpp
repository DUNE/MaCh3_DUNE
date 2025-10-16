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
  pot = SampleManager->raw()["POT"].as<double>();
  B_field = SampleManager->raw()["DetectorVariables"]["B_field"].as<double>(); //NK B field value in T
  momentum_resolution_threshold = SampleManager->raw()["DetectorVariables"]["momentum_resolution_threshold"].as<double>(); //NK momentum_resolution threshold, total as a fraction of momentum
  pixel_spacing = SampleManager->raw()["DetectorVariables"]["pixel_spacing"].as<double>(); //NK pixel spacing in mm to find num hits in y,z plane
  spatial_resolution = SampleManager->raw()["DetectorVariables"]["spatial_resolution"].as<double>(); //NK spatial resolution in mm to find  in y,z plane
  adc_sampling_frequency = SampleManager->raw()["DetectorVariables"]["adc_sampling_frequency"].as<double>(); //NK sampling frequency for ADC - needed to find timing resolution and spatial resolution in x dir in MHz
  drift_velocity = SampleManager->raw()["DetectorVariables"]["drift_velocity"].as<double>(); //NK drift velocity of electrons in gas - needed to find timing resolution and spatial resolution in x dir in cm/microsecond
  downsampling = SampleManager->raw()["DetectorVariables"]["downsampling"].as<double>(); //JM downsampling fraction
  TPCFidLength = SampleManager->raw()["DetectorVariables"]["TPCFidLength"].as<double>();
  TPCFidRadius = SampleManager->raw()["DetectorVariables"]["TPCFidRadius"].as<double>();
  TPCInstrumentedLength = SampleManager->raw()["DetectorVariables"]["TPCInstrumentedLength"].as<double>();
  TPCInstrumentedRadius = SampleManager->raw()["DetectorVariables"]["TPCInstrumentedRadius"].as<double>();
  ECALInnerRadius = SampleManager->raw()["DetectorVariables"]["ECALInnerRadius"].as<double>();
  ECALOuterRadius = SampleManager->raw()["DetectorVariables"]["ECALOuterRadius"].as<double>();
  ECALEndCapStart = SampleManager->raw()["DetectorVariables"]["ECALEndCapStart"].as<double>();
  ECALEndCapEnd = SampleManager->raw()["DetectorVariables"]["ECALEndCapEnd"].as<double>();
}

void SampleHandlerBeamNDGAr::SetupSplines() {
  ///@todo move all of the spline setup into core
  if(ParHandler->GetNumParamsFromSampleName(SampleName, kSpline) > 0){
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", ParHandler->GetNumParamsFromSampleName(SampleName, kSpline));
    SplineHandler = std::unique_ptr<BinnedSplineHandler>(new BinnedSplineHandlerDUNE(ParHandler,Modes.get()));
    InitialiseSplineObject();
  }
  else{
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or evaluate splines", ParHandler->GetNumParamsFromSampleName(SampleName, kSpline));
    SplineHandler = nullptr;
  }

  return;
}

void SampleHandlerBeamNDGAr::SetupWeightPointers() {
  for (int i = 0; i < static_cast<int>(dunendgarmcFitting.size()); ++i) {
      MCSamples[i].total_weight_pointers.push_back(&(dunendgarmcFitting[i].pot_s));
      MCSamples[i].total_weight_pointers.push_back(&(dunendgarmcFitting[i].norm_s));
      MCSamples[i].total_weight_pointers.push_back(MCSamples[i].osc_w_pointer);
      MCSamples[i].total_weight_pointers.push_back(&(dunendgarmcFitting[i].rw_berpaacvwgt));
      MCSamples[i].total_weight_pointers.push_back(&(dunendgarmcFitting[i].flux_w));
      MCSamples[i].total_weight_pointers.push_back(&(MCSamples[i].xsec_w));
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
    case 2212:  // p+
    case -3112: // anti Sig-
    case 3222:  // Sig+
    case -3312: // Xi+
    case -3334: // Omega+
      return 1;

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
bool SampleHandlerBeamNDGAr::CurvatureResolutionFilter(int id, std::unordered_map<int, std::vector<int>>& mother_to_daughter_ID, const std::unordered_map<int, size_t>& ID_to_index, dunemc_plotting& plotting_vars, double pixel_spacing_cm) {
  if (IsResolvedFromCurvature(plotting_vars, static_cast<int>(ID_to_index.at(id)), pixel_spacing_cm)) { // If mother can be reconstructed from curvature, remove her and her descendants
    int index = static_cast<int>(ID_to_index.at(id));
    int motherID = _MotherTrkID->at(index);
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

// Returns depth of a coordinate position in the dodecoganol ecal, or +999 if it is beyond the ecal boundaries
double SampleHandlerBeamNDGAr::GetCalDepth(double x, double y, double z) {
  x = x - TPC_centre_x;
  y = y - TPC_centre_y;
  z = z - TPC_centre_z;
  double r = std::sqrt(y*y + z*z);
  if (x > ECALEndCapEnd) return 999.;
  if (r < ECALInnerRadius) return std::abs(x) - ECALEndCapStart; // give depth relative to end cap

  double theta = atan2(y, z);
  double theta_segment = std::round(theta/(M_PI/6.))*M_PI/6.;
  double projected_r = r*cos(theta_segment - theta);
  if (projected_r > ECALOuterRadius) return 999.;

  return projected_r - ECALInnerRadius; // give depth relative to barrel
}

// Calculate the layer from the depth. Barrel: 8x0.673cm, 34x1.142cm. Endcap: 6x0.673cm, 36x1.142cm.
double SampleHandlerBeamNDGAr::DepthToLayer(double depth, double r) {
  int n_thin_layers;
  if (r < ECALInnerRadius) n_thin_layers = 6;
  else n_thin_layers = 8;

  double layer;
  if (depth < static_cast<double>(n_thin_layers)*0.673) {
    layer = static_cast<int>(depth/0.673) + 1;
  } else {
    layer = static_cast<int>((depth-static_cast<double>(n_thin_layers)*0.673)/1.142) + 1 + n_thin_layers;
  }

  return layer;
}

double SampleHandlerBeamNDGAr::CalcEDepCal(int motherID, std::unordered_map<int, std::vector<int>>& mother_to_daughter_ID, const std::unordered_map<int, std::vector<double>>& ID_to_ECalDep, const int tot_layers, const int crit_layers) {
  auto it = mother_to_daughter_ID.find(motherID);
  double EDepCrit = 0.;
  if (it != mother_to_daughter_ID.end()) {
    for (int i_layer=tot_layers-crit_layers; i_layer<tot_layers; i_layer++) {
      EDepCrit += ID_to_ECalDep.at(motherID)[i_layer];
    }
    for (int daughterID : it->second) {
      EDepCrit += CalcEDepCal(daughterID, mother_to_daughter_ID, ID_to_ECalDep, tot_layers, crit_layers);
    }
  }
  return EDepCrit;
}

bool SampleHandlerBeamNDGAr::IsResolvedFromCurvature(dunemc_plotting& plotting_vars, int i_anapart, double pixel_spacing_cm){
  // Get particle properties from Anatree
  double xstart = _MCPStartX->at(i_anapart);
  double ystart = _MCPStartY->at(i_anapart);
  double zstart = _MCPStartZ->at(i_anapart);
  double start_radius = std::sqrt((ystart-TPC_centre_y)*(ystart-TPC_centre_y) + (zstart-TPC_centre_z)*(zstart-TPC_centre_z));
  double start_length = xstart - TPC_centre_x;
  bool starts_in_tpc = std::abs(start_length)<=TPCInstrumentedLength && start_radius<=TPCInstrumentedRadius;
  if (!starts_in_tpc) return false;

  int pdg = _PDG->at(i_anapart);
  if (pdg > 1000000000) return false;
  int charge = GetChargeFromPDG(pdg);
  double mass = MaCh3Utils::GetMassFromPDG(pdg);
  double p_x = _MCPStartPX->at(i_anapart);
  double p_y = _MCPStartPY->at(i_anapart);
  double p_z = _MCPStartPZ->at(i_anapart);
  double p_beam = p_y*BeamDirection[1] + p_z*BeamDirection[2];

  double energy = std::sqrt(p_x*p_x + p_y*p_y + p_z*p_z + mass*mass);
  double transverse_mom = std::sqrt(p_y*p_y + p_z*p_z); //transverse to B-field
  double mom_tot = std::sqrt(p_x*p_x + transverse_mom*transverse_mom);
  double end_mom = std::sqrt(_MCPEndPX->at(i_anapart)*_MCPEndPX->at(i_anapart) + _MCPEndPY->at(i_anapart)*_MCPEndPY->at(i_anapart) + _MCPEndPZ->at(i_anapart)*_MCPEndPZ->at(i_anapart)); 

  double xend = _MCPEndX->at(i_anapart);
  double yend = _MCPEndY->at(i_anapart);
  double zend = _MCPEndZ->at(i_anapart);
  double end_radius = std::sqrt((yend-TPC_centre_y)*(yend-TPC_centre_y) + (zend-TPC_centre_z)*(zend-TPC_centre_z)); 
  double end_length = xend-TPC_centre_x;
  double end_ecaldepth = GetCalDepth(xend, yend, zend);

  bool stops_in_tpc = std::abs(end_length)<=TPCInstrumentedLength && end_radius<=TPCInstrumentedRadius;
  bool stops_before_ecal = end_ecaldepth < 0.;
  bool stops_beyond_ecal = (end_ecaldepth == 999.);
  bool stops_in_ecal = !(stops_before_ecal || stops_beyond_ecal);

  //Fill particle-level kinematic variables with default or actual (if possible at this stage) values
  if (_MotherTrkID->at(i_anapart) == 0) {
    plotting_vars.particle_energy[i_anapart] = energy;
    plotting_vars.particle_momentum[i_anapart] = mom_tot;
    plotting_vars.particle_endmomentum[i_anapart] = end_mom;
    plotting_vars.particle_transversemomentum[i_anapart] = transverse_mom; //momentum transverse to B-field
    plotting_vars.particle_bangle[i_anapart] = acos(p_x/mom_tot)*180/M_PI; //angle to B-field
    plotting_vars.particle_beamangle[i_anapart] = acos(p_beam/mom_tot)*180/M_PI;

    plotting_vars.particle_startx[i_anapart] = start_length;
    plotting_vars.particle_startr2[i_anapart] = start_radius*start_radius;
    plotting_vars.particle_endr[i_anapart] = end_radius;
    plotting_vars.particle_enddepth[i_anapart] = end_ecaldepth;
    plotting_vars.particle_endx[i_anapart] = end_length;
    plotting_vars.particle_endy[i_anapart] = yend-TPC_centre_y;
    plotting_vars.particle_endz[i_anapart] = zend-TPC_centre_z;

    plotting_vars.particle_isdecayed[i_anapart] = _MCPEndProc->at(i_anapart) == "Decay";
    plotting_vars.particle_isstoppedingap[i_anapart] = !stops_in_tpc && stops_before_ecal;
    plotting_vars.particle_isstoppedinbarrelgap[i_anapart] = !stops_in_tpc && stops_before_ecal && std::abs(end_length)<=TPCInstrumentedLength; 
    plotting_vars.particle_isstoppedinendgap[i_anapart] = !stops_in_tpc && stops_before_ecal && std::abs(end_length)>TPCInstrumentedLength; 
    plotting_vars.particle_isstoppedinecal[i_anapart] = stops_in_ecal;
    plotting_vars.particle_isstoppedinbarrel[i_anapart] = stops_in_ecal && end_radius>=ECALInnerRadius;
    plotting_vars.particle_isstoppedinendcap[i_anapart] = stops_in_ecal && end_radius<ECALInnerRadius;
    plotting_vars.particle_isstoppedintpc[i_anapart] = stops_in_tpc;
    plotting_vars.particle_isescaped[i_anapart] = stops_beyond_ecal;
  }

  if (charge == 0) return false;

  double rad_curvature = 100*transverse_mom/(0.3*B_field); //p = 0.3*B*r where p in GeV/c, B in T, r in m (*100 to convert to cm)
  double theta_xT = atan(p_x/transverse_mom); //helix pitch angle
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
    theta_spanned = std::abs(xstart - _MCPEndX->at(i_anapart))*2*M_PI/pitch; 
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
    * std::sqrt(1-(1/21)*(L_yz/rad_curvature)*(L_yz/rad_curvature)));
  double momres_ms = transverse_mom*(0.016/(0.3*B_field*(L_yz/100)*cos(theta_xT)*beta))*std::sqrt(L_yz/X0);
  double momres_tottransverse = std::sqrt(momres_yz*momres_yz + momres_ms*momres_ms)/transverse_mom;
  double sigma_theta = (cos(theta_xT)*cos(theta_xT) * (pitch/(2*M_PI*rad_curvature)) *
    std::sqrt(sigmax_frac*sigmax_frac + momres_tottransverse*momres_tottransverse));
  double momres_frac = std::sqrt(momres_tottransverse*momres_tottransverse + (sigma_theta*tan_theta)*(sigma_theta*tan_theta));

  if (_MotherTrkID->at(i_anapart) == 0) {
    plotting_vars.particle_nturns[i_anapart] = nturns;
    plotting_vars.particle_nhits[i_anapart] = nhits;
    plotting_vars.particle_tracklengthyz[i_anapart] = L_yz;
    plotting_vars.particle_momresms[i_anapart] = momres_ms/transverse_mom;
    plotting_vars.particle_momresyz[i_anapart] = momres_yz/transverse_mom;
    plotting_vars.particle_momresx[i_anapart] = std::abs(sigma_theta*tan_theta);
  }

  if(momres_frac > momentum_resolution_threshold) return false;
  return true;
}

void SampleHandlerBeamNDGAr::clearBranchVectors() {
  _MCVertX->clear();
  _MCVertY->clear();
  _MCVertZ->clear();
  _MCNuPx->clear();
  _MCNuPy->clear();
  _MCNuPz->clear();
  _IsNC->clear();
  _MCMode->clear();
  _MCPStartX->clear();
  _MCPStartY->clear();
  _MCPStartZ->clear();
  _MCPEndX->clear();
  _MCPEndY->clear();
  _MCPEndZ->clear();
  _MCPStartPX->clear();
  _MCPStartPY->clear();
  _MCPStartPZ->clear();
  _MCPEndPX->clear();
  _MCPEndPY->clear();
  _MCPEndPZ->clear();
  _PDG ->clear();
  _MCPTrkID->clear();
  _MCPProc->clear();
  _MCPEndProc->clear();
  _MotherTrkID->clear();
  _SimHitTrkID->clear();
  _SimHitLayer->clear();
  _SimHitEnergy->clear();
  _SimHitX->clear();
  _SimHitY->clear();
  _SimHitZ->clear();
}

int SampleHandlerBeamNDGAr::SetupExperimentMC() {

  MACH3LOG_INFO("-------------------------------------------------------------------");

  TChain* _data = new TChain("GArAnaTree");
  for (size_t iSample=0;iSample<mc_files.size();iSample++) {
    MACH3LOG_INFO("Adding file to TChain: {}", mc_files[iSample]);

    int _add_rtn = _data->Add(mc_files[iSample].c_str(), -1);
    if(_add_rtn == 0){
      MACH3LOG_ERROR("Could not add file {} to TChain, please check the file exists and is readable", mc_files[iSample]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  _data->SetBranchStatus("*", 0);

  auto readBranch = [&](const char* name, void* addr) {
    _data->SetBranchStatus(name, 1);
    _data->SetBranchAddress(name, addr);
  };

  readBranch("MCVertX", &_MCVertX);
  readBranch("MCVertY", &_MCVertY);
  readBranch("MCVertZ", &_MCVertZ);
  readBranch("MCNuPx", &_MCNuPx);
  readBranch("MCNuPy", &_MCNuPy);
  readBranch("MCNuPz", &_MCNuPz);
  readBranch("CCNC", &_IsNC);
  readBranch("Mode", &_MCMode);
  readBranch("MCPStartX", &_MCPStartX);
  readBranch("MCPStartY", &_MCPStartY);
  readBranch("MCPStartZ", &_MCPStartZ);
  readBranch("MCPEndX", &_MCPEndX);
  readBranch("MCPEndY", &_MCPEndY);
  readBranch("MCPEndZ", &_MCPEndZ);
  readBranch("MCPStartPX", &_MCPStartPX);
  readBranch("MCPStartPY", &_MCPStartPY);
  readBranch("MCPStartPZ", &_MCPStartPZ);
  readBranch("MCPEndPX", &_MCPEndPX);
  readBranch("MCPEndPY", &_MCPEndPY);
  readBranch("MCPEndPZ", &_MCPEndPZ);
  readBranch("PDG", &_PDG);
  readBranch("MCPTrkID", &_MCPTrkID);
  readBranch("MCPProc", &_MCPProc);
  readBranch("MCPEndProc", &_MCPEndProc);
  readBranch("MotherTrkID", &_MotherTrkID);
  readBranch("SimHitTrkID", &_SimHitTrkID);
  readBranch("SimHitLayer", &_SimHitLayer);
  readBranch("SimHitEnergy", &_SimHitEnergy);
  readBranch("SimHitX", &_SimHitX);
  readBranch("SimHitY", &_SimHitY);
  readBranch("SimHitZ", &_SimHitZ);

  int nEntries = static_cast<int>(downsampling*static_cast<double>(_data->GetEntries()));

  dunendgarmcFitting.resize(nEntries);
  dunendgarmcPlotting.resize(nEntries);

  double pixel_spacing_cm = pixel_spacing/10; //convert to cm
  int numCC = 0;
  int num_in_fdv = 0;
  bool do_geometric_correction = false;

  for (int i_event = 0; i_event < nEntries; ++i_event) { 
    if (i_event != 0) clearBranchVectors();
    _data->GetEntry(i_event);

    if (i_event % (nEntries/100) == 0) {
      MACH3LOG_INFO("\tNow processing event: {}/{}",i_event,nEntries);
    }

    double radius = std::sqrt((_MCVertY->at(0)-TPC_centre_y)*(_MCVertY->at(0)-TPC_centre_y) + (_MCVertZ->at(0)-TPC_centre_z)*(_MCVertZ->at(0)-TPC_centre_z)); //find radius of interaction vertex

    if(std::abs(_MCVertX->at(0) - TPC_centre_x)<=TPCFidLength &&  radius<=TPCFidRadius){
      num_in_fdv++;
      dunendgarmcPlotting[i_event].in_fdv = 1;
    } else{
      dunendgarmcPlotting[i_event].in_fdv = 0;
    }

    dunendgarmcFitting[i_event].rw_etru = std::sqrt(_MCNuPx->at(0)*_MCNuPx->at(0) + _MCNuPy->at(0)*_MCNuPy->at(0) + _MCNuPz->at(0)*_MCNuPz->at(0));
    dunendgarmcFitting[i_event].rw_isCC = static_cast<int>(!_IsNC->at(0));
    dunendgarmcFitting[i_event].nupdg = 14;
    dunendgarmcFitting[i_event].nupdgUnosc = 14;
    dunendgarmcFitting[i_event].OscChannelIndex = static_cast<double>(GetOscChannel(OscChannels, dunendgarmcFitting[i_event].nupdgUnosc, dunendgarmcFitting[i_event].nupdg));
    dunendgarmcFitting[i_event].rw_berpaacvwgt = _BeRPA_cvwgt;

    dunendgarmcFitting[i_event].rw_vtx_x = _MCVertX->at(0);
    dunendgarmcFitting[i_event].rw_vtx_y = _MCVertY->at(0);
    dunendgarmcFitting[i_event].rw_vtx_z = _MCVertZ->at(0);

    std::unordered_map<int, std::vector<int>> mother_to_daughter_ID; // particle track ID -> vector of daughter IDs
    std::unordered_map<int, size_t> ID_to_index; // particle track ID -> index in anatree
    std::unordered_map<int, std::vector<double>> ID_to_ECalDep; // particle track ID -> total energy deposited in each ecal layer
    size_t n_particles_in_event = _MCPTrkID->size();
    const int tot_ecal_layers = 42;

    // Fill maps
    for (size_t i_anapart=0; i_anapart<n_particles_in_event; i_anapart++) {
      int secID = _MCPTrkID->at(i_anapart);
      int motherID = _MotherTrkID->at(i_anapart);
      ID_to_ECalDep[secID] = std::vector<double>(tot_ecal_layers,0.);
      ID_to_index[secID] = i_anapart;
      mother_to_daughter_ID[motherID].push_back(secID);
      mother_to_daughter_ID[secID]; // Ensure all particles are added to the map (even if no secondaries)
    }

    for (size_t i_ecaldep=0; i_ecaldep<_SimHitTrkID->size(); i_ecaldep++) {
      int simhit_trkid = _SimHitTrkID->at(i_ecaldep);
      int simhit_layer = _SimHitLayer->at(i_ecaldep);
      if (simhit_trkid <= 0 || simhit_layer <= 0) {
        continue;
      }
      ID_to_ECalDep[simhit_trkid][simhit_layer-1] += _SimHitEnergy->at(i_ecaldep);

      double hitx = _SimHitX->at(i_ecaldep);
      double hity = _SimHitY->at(i_ecaldep);
      double hitz = _SimHitZ->at(i_ecaldep);
      double hitr = std::sqrt((hity - TPC_centre_y)*(hity - TPC_centre_y) + (hitz - TPC_centre_z)*(hitz - TPC_centre_z));

      bool isEndCapHit = hitr < ECALInnerRadius; // Whether the hit is in the end cap radially (agnostic to x coordinate)
      double dep_depth = GetCalDepth(hitx, hity, hitz);
      double layer_calc = DepthToLayer(dep_depth, hitr);

      if (layer_calc < 0 || layer_calc > 42) MACH3LOG_INFO("Hit recorded outside ECAL dimensions. IsEndCap: {}. TrueLayer: {}. Calculated layer: {}. Depth: {}.", isEndCapHit, simhit_layer, layer_calc, dep_depth);
      if (layer_calc != simhit_layer) {
        MACH3LOG_INFO("Miscalculated layer for hit. IsEndCap: {}. Simhitlayer = {}, calculated layer = {}, depth = {} cm.", isEndCapHit, simhit_layer, layer_calc, dep_depth);
      }
    }

    double muon_p2 = 0.;
    bool isEventAccepted = true;
    int crit_layers = 2; // Number of outer layers of the calorimeter forming the 'critical' region

    // Resize vectors for particle-level parameters
    size_t n_prim_in_event = mother_to_daughter_ID[0].size();
    dunendgarmcPlotting[i_event].particle_pdg.resize(n_prim_in_event, M3::_BAD_INT_);
    dunendgarmcPlotting[i_event].particle_energy.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_momentum.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_endmomentum.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_transversemomentum.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_bangle.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_beamangle.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_isaccepted.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].particle_iscurvatureresolved.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].particle_isdecayed.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].particle_isstoppedintpc.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].particle_isstoppedinecal.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].particle_isstoppedingap.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].particle_isstoppedinbarrelgap.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].particle_isstoppedinendgap.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].particle_isstoppedinbarrel.resize(n_prim_in_event, M3::_BAD_INT_);
    dunendgarmcPlotting[i_event].particle_isstoppedinendcap.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].particle_isescaped.resize(n_prim_in_event, M3::_BAD_INT_); 
    dunendgarmcPlotting[i_event].particle_startx.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_startr2.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_endr.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_enddepth.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_endx.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_endy.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_endz.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_nturns.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_nhits.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_tracklengthyz.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_momresms.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_momresyz.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_momresx.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 
    dunendgarmcPlotting[i_event].particle_edepcrit.resize(n_prim_in_event, M3::_BAD_DOUBLE_); 

    // Loop through primaries
    for (int& primID : mother_to_daughter_ID[0]) {

      // Do not require the reconstruction of neutrons and neutrinos
      size_t i_anaprim = ID_to_index[primID];
      int pdg = _PDG->at(i_anaprim);
      if (pdg == 2112 || std::abs(pdg) == 12 || std::abs(pdg) == 14 || std::abs(pdg) == 16) continue;

      dunendgarmcPlotting[i_event].particle_edepcrit[i_anaprim] = CalcEDepCal(primID, mother_to_daughter_ID, ID_to_ECalDep, tot_ecal_layers, crit_layers);
      dunendgarmcPlotting[i_event].particle_pdg[i_anaprim] = pdg; 

      // Check if primary is resolved from curvature
      bool isCurvatureResolved = false;
      // Remove secondaries (and their descendants) from mother_to_daughter_ID who's momentum we get from curvature
      if (CurvatureResolutionFilter(primID, mother_to_daughter_ID, ID_to_index, dunendgarmcPlotting[i_event], pixel_spacing_cm)) {
        isCurvatureResolved = true;
      }
      dunendgarmcPlotting[i_event].particle_iscurvatureresolved[i_anaprim] = isCurvatureResolved;

      // Find energy deposited by by primary and non-curvature-resolved descendants in critical region of calorimeter
      double EDepCrit = CalcEDepCal(primID, mother_to_daughter_ID, ID_to_ECalDep, tot_ecal_layers, crit_layers);
      // Check for containment
      bool isContained = true;
      if (EDepCrit > 0.002) {
        isContained = false;
      }

      // Primary is accepted if contained or curvature resolved
      bool isParticleAccepted = true;
      if (!(isContained || isCurvatureResolved)) {
        isParticleAccepted = false;
        isEventAccepted = false;
      } 
      dunendgarmcPlotting[i_event].particle_isaccepted[i_anaprim] = isParticleAccepted;

      // Get primary muon information
      if(pdg == 13) {
        double p_x = _MCPStartPX->at(i_anaprim);
        double p_y = _MCPStartPY->at(i_anaprim);
        double p_z = _MCPStartPZ->at(i_anaprim);
        double p2 = p_x*p_x + p_y*p_y + p_z*p_z;
        if (p2 > muon_p2) {
          muon_p2 = p2;
          double mass = MaCh3Utils::GetMassFromPDG(pdg);
          dunendgarmcFitting[i_event].rw_lep_pX = p_x;
          dunendgarmcFitting[i_event].rw_lep_pY = p_y;
          dunendgarmcFitting[i_event].rw_lep_pZ = p_z;
          dunendgarmcFitting[i_event].rw_LepE = std::sqrt(p2 + mass*mass);
        }
      }
    }
    dunendgarmcPlotting[i_event].is_accepted = isEventAccepted;

    // Find lepton kinematic variables
    double lep_momentum = std::sqrt(dunendgarmcFitting[i_event].rw_lep_pX*dunendgarmcFitting[i_event].rw_lep_pX + dunendgarmcFitting[i_event].rw_lep_pY*dunendgarmcFitting[i_event].rw_lep_pY + dunendgarmcFitting[i_event].rw_lep_pZ*dunendgarmcFitting[i_event].rw_lep_pZ);
    double lep_pBeam = dunendgarmcFitting[i_event].rw_lep_pY*BeamDirection[1] + dunendgarmcFitting[i_event].rw_lep_pZ*BeamDirection[2];
    double lep_pB = dunendgarmcFitting[i_event].rw_lep_pX;
    double lep_pPerp = dunendgarmcFitting[i_event].rw_lep_pY*BeamDirection[2] - dunendgarmcFitting[i_event].rw_lep_pZ*BeamDirection[1];
    double lep_beamangle = acos(lep_pBeam/lep_momentum)*180/M_PI; //Angle to beam (beam direction: [0.0,-0.101,0.995])
    double lep_bangle = acos(lep_pB/lep_momentum)*180/M_PI; //Angle to B-field (b-field along x)
    double lep_perpangle = acos(lep_pPerp/lep_momentum)*180/M_PI; //Angle to axis perpendicular to beam and B
    double lep_phi = atan2(lep_pPerp, lep_pB)*180/M_PI;

    dunendgarmcFitting[i_event].rw_lep_theta = lep_beamangle;
    dunendgarmcFitting[i_event].rw_lep_phi = lep_phi;
    dunendgarmcFitting[i_event].rw_lep_bangle = lep_bangle;
    dunendgarmcFitting[i_event].rw_lep_p = lep_momentum;

    // Perform 'geometric correction' if do_geometric_correction set to true
    dunendgarmcPlotting[i_event].geometric_correction = 1.;
    if (do_geometric_correction) {
      if ((lep_bangle < 45 || lep_bangle > 135) && lep_momentum > 0.3) dunendgarmcPlotting[i_event].geometric_correction = 0.;
      else if ((lep_perpangle < 45 || lep_perpangle > 135) && lep_momentum > 0.3) dunendgarmcPlotting[i_event].geometric_correction = 2.;
    }

    // Fill remaining event-level kinematic parameters
    dunendgarmcFitting[i_event].rw_rad = radius;
    dunendgarmcFitting[i_event].Target = 40; // Assume everything is Argon
    dunendgarmcFitting[i_event].rw_Q0 = dunendgarmcFitting[i_event].rw_etru - dunendgarmcFitting[i_event].rw_LepE;
    dunendgarmcFitting[i_event].rw_Q3 = std::sqrt((_MCNuPx->at(0)-dunendgarmcFitting[i_event].rw_lep_pX)*(_MCNuPx->at(0)-dunendgarmcFitting[i_event].rw_lep_pX) + 
                                                  (_MCNuPy->at(0)-dunendgarmcFitting[i_event].rw_lep_pY)*(_MCNuPy->at(0)-dunendgarmcFitting[i_event].rw_lep_pY) + 
                                                  (_MCNuPz->at(0)-dunendgarmcFitting[i_event].rw_lep_pZ)*(_MCNuPz->at(0)-dunendgarmcFitting[i_event].rw_lep_pZ));
    dunendgarmcFitting[i_event].rw_lep_pT = std::sqrt(lep_momentum*lep_momentum - lep_pBeam*lep_pBeam); 
    dunendgarmcFitting[i_event].mode = _MCMode->at(0);
    dunendgarmcFitting[i_event].norm_s = 1.;
    dunendgarmcFitting[i_event].pot_s = pot/(downsampling*1e21);
    dunendgarmcFitting[i_event].flux_w = 1.0;
    if(dunendgarmcFitting[i_event].rw_isCC == 1) numCC++;
  }
  MACH3LOG_INFO("nEntries = {}, numCC = {}, numFDV = {}", nEntries, numCC, num_in_fdv);

  _data->Reset();
  delete _data;
  return nEntries;
}

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wswitch-enum"
const double* SampleHandlerBeamNDGAr::GetPointerToKinematicParameter(KinematicTypes KinematicParameter, int iEvent) {
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
      return &dunendgarmcFitting[iEvent].rw_lep_theta;
    case kLepPhi:
      return &dunendgarmcFitting[iEvent].rw_lep_phi;
    case kLepBAngle:
      return &dunendgarmcFitting[iEvent].rw_lep_bangle;
    case kLepP:
      return &dunendgarmcFitting[iEvent].rw_lep_p;
    case kTrueQ0:
      return &dunendgarmcFitting[iEvent].rw_Q0;
    case kTrueQ3:
      return &dunendgarmcFitting[iEvent].rw_Q3;
    default:
      MACH3LOG_ERROR("Did not recognise Kinematic Parameter {}", static_cast<int>(KinematicParameter));
      throw MaCh3Exception(__FILE__, __LINE__);
      return nullptr;
  }
}
#pragma GCC diagnostic pop

const double* SampleHandlerBeamNDGAr::GetPointerToKinematicParameter(double KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(std::round(KinematicVariable));
  return GetPointerToKinematicParameter(KinPar,iEvent);
}

const double* SampleHandlerBeamNDGAr::GetPointerToKinematicParameter(std::string KinematicParameter, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return GetPointerToKinematicParameter(KinPar,iEvent);
}

double SampleHandlerBeamNDGAr::ReturnKinematicParameter(int KinematicVariable, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(std::round(KinematicVariable));
  return ReturnKinematicParameter(KinPar,iEvent);
}

double SampleHandlerBeamNDGAr::ReturnKinematicParameter(std::string KinematicParameter, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return ReturnKinematicParameter(KinPar,iEvent);
}

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wswitch-enum"
double SampleHandlerBeamNDGAr::ReturnKinematicParameter(KinematicTypes KinPar, int iEvent) {
  //HH: Special cases for dealing with non-doubles
  switch(KinPar) {
    case kEvent_IsAccepted:
      return static_cast<double>(dunendgarmcPlotting[iEvent].is_accepted);
    case kIsCC:
      return static_cast<double>(dunendgarmcFitting[iEvent].rw_isCC);
    case kInFDV:
      return static_cast<double>(dunendgarmcPlotting[iEvent].in_fdv);
    default:
      return *GetPointerToKinematicParameter(KinPar, iEvent);
  }
}
#pragma GCC diagnostic pop

std::vector<double> SampleHandlerBeamNDGAr::ReturnKinematicVector(KinematicVecs KinVec, int iEvent) {
  switch(KinVec) {
    case kParticle_IsAccepted:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].particle_isaccepted.begin(),
        dunendgarmcPlotting[iEvent].particle_isaccepted.end());
    case kParticle_IsCurvatureResolved:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].particle_iscurvatureresolved.begin(),
        dunendgarmcPlotting[iEvent].particle_iscurvatureresolved.end());
    case kParticle_IsDecayed:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].particle_isdecayed.begin(),
        dunendgarmcPlotting[iEvent].particle_isdecayed.end());
    case kParticle_PDG:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].particle_pdg.begin(),
        dunendgarmcPlotting[iEvent].particle_pdg.end());
    case kParticle_IsStoppedInTPC:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].particle_isstoppedintpc.begin(),
        dunendgarmcPlotting[iEvent].particle_isstoppedintpc.end());
    case kParticle_IsStoppedInECal:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].particle_isstoppedinecal.begin(),
        dunendgarmcPlotting[iEvent].particle_isstoppedinecal.end());
    case kParticle_IsStoppedInBarrel:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].particle_isstoppedinbarrel.begin(),
        dunendgarmcPlotting[iEvent].particle_isstoppedinbarrel.end());
    case kParticle_IsStoppedInEndCap:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].particle_isstoppedinendcap.begin(),
        dunendgarmcPlotting[iEvent].particle_isstoppedinendcap.end());
    case kParticle_IsStoppedInGap:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].particle_isstoppedingap.begin(),
        dunendgarmcPlotting[iEvent].particle_isstoppedingap.end());
    case kParticle_IsStoppedInEndGap:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].particle_isstoppedinendgap.begin(),
        dunendgarmcPlotting[iEvent].particle_isstoppedinendgap.end());
    case kParticle_IsStoppedInBarrelGap:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].particle_isstoppedinbarrelgap.begin(),
        dunendgarmcPlotting[iEvent].particle_isstoppedinbarrelgap.end());
    case kParticle_IsEscaped:
      return std::vector<double>(
        dunendgarmcPlotting[iEvent].particle_isescaped.begin(),
        dunendgarmcPlotting[iEvent].particle_isescaped.end());
    case kParticle_Energy:
      return dunendgarmcPlotting[iEvent].particle_energy;
    case kParticle_Momentum:
      return dunendgarmcPlotting[iEvent].particle_momentum;
    case kParticle_EndMomentum:
      return dunendgarmcPlotting[iEvent].particle_endmomentum;
    case kParticle_TransverseMomentum:
      return dunendgarmcPlotting[iEvent].particle_transversemomentum;
    case kParticle_BAngle:
      return dunendgarmcPlotting[iEvent].particle_bangle;
    case kParticle_BeamAngle:
      return dunendgarmcPlotting[iEvent].particle_beamangle;
    case kParticle_NTurns:
      return dunendgarmcPlotting[iEvent].particle_nturns;
    case kParticle_NHits:
      return dunendgarmcPlotting[iEvent].particle_nhits;
    case kParticle_TrackLengthYZ:
      return dunendgarmcPlotting[iEvent].particle_tracklengthyz;
    case kParticle_MomResMS:
      return dunendgarmcPlotting[iEvent].particle_momresms;
    case kParticle_MomResYZ:
      return dunendgarmcPlotting[iEvent].particle_momresyz;
    case kParticle_MomResX:
      return dunendgarmcPlotting[iEvent].particle_momresx;
    case kParticle_StartR2:
      return dunendgarmcPlotting[iEvent].particle_startr2;
    case kParticle_EndR:
      return dunendgarmcPlotting[iEvent].particle_endr;
    case kParticle_EndDepth:
      return dunendgarmcPlotting[iEvent].particle_enddepth;
    case kParticle_EndX:
      return dunendgarmcPlotting[iEvent].particle_endx;
    case kParticle_EndY:
      return dunendgarmcPlotting[iEvent].particle_endy;
    case kParticle_EndZ:
      return dunendgarmcPlotting[iEvent].particle_endz;
    case kParticle_StartX:
      return dunendgarmcPlotting[iEvent].particle_startx;
    case kParticle_EDepCrit:
      return dunendgarmcPlotting[iEvent].particle_edepcrit;
    default:
      MACH3LOG_ERROR("Unrecognized Kinematic Vector: {}", static_cast<int>(KinVec));
      throw MaCh3Exception(__FILE__, __LINE__);
  }
}

std::vector<double> SampleHandlerBeamNDGAr::ReturnKinematicVector(int KinematicVector, int iEvent) {
  KinematicVecs KinVec = static_cast<KinematicVecs>(KinematicVector);
  return ReturnKinematicVector(KinVec, iEvent);
}

std::vector<double> SampleHandlerBeamNDGAr::ReturnKinematicVector(std::string KinematicVector, int iEvent) {
  KinematicVecs KinVec = static_cast<KinematicVecs>(ReturnKinematicVectorFromString(KinematicVector));
  return ReturnKinematicVector(KinVec, iEvent);
}

void SampleHandlerBeamNDGAr::SetupFDMC() {
  for(int iEvent = 0 ;iEvent < int(GetNEvents()) ; ++iEvent){
    MCSamples[iEvent].rw_etru = &(dunendgarmcFitting[iEvent].rw_etru);
    MCSamples[iEvent].mode = &(dunendgarmcFitting[iEvent].mode);
    MCSamples[iEvent].Target = &(dunendgarmcFitting[iEvent].Target);    
    MCSamples[iEvent].isNC = !dunendgarmcFitting[iEvent].rw_isCC;
    MCSamples[iEvent].nupdg = &(dunendgarmcFitting[iEvent].nupdg);
    MCSamples[iEvent].nupdgUnosc = &(dunendgarmcFitting[iEvent].nupdgUnosc);
  }
}
