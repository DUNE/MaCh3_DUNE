#include "SampleHandlerBeamNDGAr.h"

SampleHandlerBeamNDGAr::SampleHandlerBeamNDGAr(std::string mc_version_, ParameterHandlerGeneric* ParHandler_) : SampleHandlerFD(mc_version_, ParHandler_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;

  Initialise();
}

SampleHandlerBeamNDGAr::~SampleHandlerBeamNDGAr() {
}

void SampleHandlerBeamNDGAr::Init() {
  dunendgarmcSamples.resize(nSamples,dunemc_base());
  nparticlesinsample = new int[nSamples]();
  pot = SampleManager->raw()["POT"].as<double>();
  B_field = SampleManager->raw()["SampleCuts"]["B_field"].as<double>(); //NK B field value in T
  momentum_resolution_threshold = SampleManager->raw()["SampleCuts"]["momentum_resolution_threshold"].as<double>(); //NK momentum_resolution threshold, total as a fraction of momentum
  pixel_spacing = SampleManager->raw()["SampleCuts"]["pixel_spacing"].as<double>(); //NK pixel spacing in mm to find num hits in y,z plane
  spatial_resolution = SampleManager->raw()["SampleCuts"]["spatial_resolution"].as<double>(); //NK spatial resolution in mm to find  in y,z plane
  adc_sampling_frequency = SampleManager->raw()["SampleCuts"]["adc_sampling_frequency"].as<double>(); //NK sampling frequency for ADC - needed to find timing resolution and spatial resolution in x dir in MHz
  drift_velocity = SampleManager->raw()["SampleCuts"]["drift_velocity"].as<double>(); //NK drift velocity of electrons in gas - needed to find timing resolution and spatial resolution in x dir in cm/microsecond
  TPCFidLength = SampleManager->raw()["SampleCuts"]["TPCFidLength"].as<double>();
  TPCFidRadius = SampleManager->raw()["SampleCuts"]["TPCFidRadius"].as<double>();
  TPCInstrumentedLength = SampleManager->raw()["SampleCuts"]["TPCInstrumentedLength"].as<double>();
  TPCInstrumentedRadius = SampleManager->raw()["SampleCuts"]["TPCInstrumentedRadius"].as<double>();
  ECALInnerRadius = SampleManager->raw()["SampleCuts"]["ECALInnerRadius"].as<double>();
  ECALOuterRadius = SampleManager->raw()["SampleCuts"]["ECALOuterRadius"].as<double>();
  ECALEndCapStart = SampleManager->raw()["SampleCuts"]["ECALEndCapStart"].as<double>();
  ECALEndCapEnd = SampleManager->raw()["SampleCuts"]["ECALEndCapEnd"].as<double>();
}

void SampleHandlerBeamNDGAr::SetupSplines() {
  ///@todo move all of the spline setup into core
  if(ParHandler->GetNumParamsFromSampleName(SampleName, kSpline) > 0){
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", ParHandler->GetNumParamsFromSampleName(SampleName, kSpline));
    SplineHandler = std::unique_ptr<BinnedSplineHandler>(new BinnedSplineHandlerDUNE(ParHandler,Modes));
    InitialiseSplineObject();
  }
  else{
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or evaluate splines", ParHandler->GetNumParamsFromSampleName(SampleName, kSpline));
    SplineHandler = nullptr;
  }

  return;
}

void SampleHandlerBeamNDGAr::SetupWeightPointers() {
  for (int i = 0; i < static_cast<int>(dunendgarmcSamples.size()); ++i) {
    for (int j = 0; j < dunendgarmcSamples[i].nEvents; ++j) {
      MCSamples[i].ntotal_weight_pointers[j] = 7;
      //MCSamples[i].total_weight_pointers[j] = new const double*[MCSamples[i].ntotal_weight_pointers[j]];
      MCSamples[i].total_weight_pointers[j].resize(MCSamples[i].ntotal_weight_pointers[j]);
      MCSamples[i].total_weight_pointers[j][0] = &(dunendgarmcSamples[i].pot_s);
      MCSamples[i].total_weight_pointers[j][1] = &(dunendgarmcSamples[i].norm_s);
      MCSamples[i].total_weight_pointers[j][2] = MCSamples[i].osc_w_pointer[j];
      MCSamples[i].total_weight_pointers[j][3] = &(dunendgarmcSamples[i].rw_berpaacvwgt[j]);
      MCSamples[i].total_weight_pointers[j][4] = &(dunendgarmcSamples[i].flux_w[j]);
      MCSamples[i].total_weight_pointers[j][5] = &(MCSamples[i].xsec_w[j]);
      MCSamples[i].total_weight_pointers[j][6] = &(dunendgarmcSamples[i].geometric_correction[j]);

      double totalweight=1.;
      for (int weighting=0; weighting<MCSamples[i].ntotal_weight_pointers[j]; weighting++) {
        totalweight *= *(MCSamples[i].total_weight_pointers[j][weighting]);
      }
      if (totalweight != 1) MACH3LOG_INFO("Non-unitary weighting! Event {}, Weighting {}", j, totalweight);
    }
  }
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
bool SampleHandlerBeamNDGAr::CurvatureResolutionFilter(int id, std::unordered_map<int, std::vector<int>>& mother_to_daughter_ID, const std::unordered_map<int, size_t>& ID_to_index, dunemc_base *duneobj, double pixel_spacing_cm) {
  if (IsResolvedFromCurvature(duneobj, static_cast<int>(ID_to_index.at(id)), pixel_spacing_cm)) { // If mother can be reconstructed from curvature, remove her and her descendants
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
    if (!CurvatureResolutionFilter(daughterID, mother_to_daughter_ID, ID_to_index, duneobj, pixel_spacing_cm)) {
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

bool SampleHandlerBeamNDGAr::IsResolvedFromCurvature(dunemc_base *duneobj, int i_anapart, double pixel_spacing_cm){
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
    duneobj->particle_pdg->back() = pdg;
    duneobj->particle_energy->back() = energy;
    duneobj->particle_momentum->back() = mom_tot;
    duneobj->particle_endmomentum->back() = end_mom;
    duneobj->particle_transversemomentum->back() = transverse_mom; //momentum transverse to B-field
    duneobj->particle_bangle->back() = acos(p_x/mom_tot)*180/M_PI; //angle to B-field
    duneobj->particle_beamangle->back() = acos(p_beam/mom_tot)*180/M_PI;

    duneobj->particle_startx->back() = start_length;
    duneobj->particle_startr2->back() = start_radius*start_radius;
    duneobj->particle_endr->back() = end_radius;
    duneobj->particle_enddepth->back() = end_ecaldepth;
    duneobj->particle_endx->back() = end_length;
    duneobj->particle_endy->back() = yend-TPC_centre_y;
    duneobj->particle_endz->back() = zend-TPC_centre_z;

    duneobj->particle_isdecayed->back() = _MCPEndProc->at(i_anapart) == "Decay";
    duneobj->particle_isstoppedingap->back() = !stops_in_tpc && stops_before_ecal;
    duneobj->particle_isstoppedinbarrelgap->back() = !stops_in_tpc && stops_before_ecal && std::abs(end_length)<=TPCInstrumentedLength; 
    duneobj->particle_isstoppedinendgap->back() = !stops_in_tpc && stops_before_ecal && std::abs(end_length)>TPCInstrumentedLength; 
    duneobj->particle_isstoppedinecal->back() = stops_in_ecal;
    duneobj->particle_isstoppedinbarrel->back() = stops_in_ecal && end_radius>=ECALInnerRadius;
    duneobj->particle_isstoppedinendcap->back() = stops_in_ecal && end_radius<ECALInnerRadius;
    duneobj->particle_isstoppedintpc->back() = stops_in_tpc;
    duneobj->particle_isescaped->back() = stops_beyond_ecal;
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
    duneobj->particle_nturns->back() = nturns;
    duneobj->particle_nhits->back() = nhits;
    duneobj->particle_tracklengthyz->back() = L_yz;
    duneobj->particle_momresms->back() = momres_ms/transverse_mom;
    duneobj->particle_momresyz->back() = momres_yz/transverse_mom;
    duneobj->particle_momresx->back() = std::abs(sigma_theta*tan_theta);
  }

  if(momres_frac > momentum_resolution_threshold) return false;
  return true;
}

int SampleHandlerBeamNDGAr::SetupExperimentMC(int iSample) {

  dunemc_base *duneobj = &(dunendgarmcSamples[iSample]);
  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("Input File: {}", mc_files.at(iSample));

  //Read anatree file
  _sampleFile = TFile::Open(mc_files.at(iSample).c_str(), "READ");
  _data = static_cast<TTree*>(_sampleFile->Get("GArAnaTree"));
  if(_data){
    MACH3LOG_INFO("Found \"anatree\" tree in {}", mc_files.at(iSample));
    MACH3LOG_INFO("With number of entries: {}", _data->GetEntries());
  }
  else{
    MACH3LOG_ERROR("Could not find \"anatree\" tree in {}", mc_files.at(iSample));
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // _data->SetBranchStatus("*", 1);
  _data->SetBranchAddress("MCVertX", &_MCVertX);
  _data->SetBranchAddress("MCVertY", &_MCVertY);
  _data->SetBranchAddress("MCVertZ", &_MCVertZ);
  _data->SetBranchAddress("MCNuPx", &_MCNuPx);
  _data->SetBranchAddress("MCNuPy", &_MCNuPy);
  _data->SetBranchAddress("MCNuPz", &_MCNuPz);
  _data->SetBranchAddress("CCNC", &_IsNC);
  _data->SetBranchAddress("Mode", &_MCMode);
  _data->SetBranchAddress("MCPStartX", &_MCPStartX);
  _data->SetBranchAddress("MCPStartY", &_MCPStartY);
  _data->SetBranchAddress("MCPStartZ", &_MCPStartZ);
  _data->SetBranchAddress("MCPEndX", &_MCPEndX);
  _data->SetBranchAddress("MCPEndY", &_MCPEndY);
  _data->SetBranchAddress("MCPEndZ", &_MCPEndZ);
  _data->SetBranchAddress("MCPStartPX", &_MCPStartPX);
  _data->SetBranchAddress("MCPStartPY", &_MCPStartPY);
  _data->SetBranchAddress("MCPStartPZ", &_MCPStartPZ);
  _data->SetBranchAddress("MCPEndPX", &_MCPEndPX);
  _data->SetBranchAddress("MCPEndPY", &_MCPEndPY);
  _data->SetBranchAddress("MCPEndPZ", &_MCPEndPZ);
  _data->SetBranchAddress("PDG", &_PDG);
  _data->SetBranchAddress("MCPTrkID", &_MCPTrkID);
  _data->SetBranchAddress("MCPProc", &_MCPProc);
  _data->SetBranchAddress("MCPEndProc", &_MCPEndProc);
  _data->SetBranchAddress("MotherTrkID", &_MotherTrkID);
  _data->SetBranchAddress("SimHitTrkID", &_SimHitTrkID);
  _data->SetBranchAddress("SimHitLayer", &_SimHitLayer);
  _data->SetBranchAddress("SimHitEnergy", &_SimHitEnergy);
  _data->SetBranchAddress("SimHitX", &_SimHitX);
  _data->SetBranchAddress("SimHitY", &_SimHitY);
  _data->SetBranchAddress("SimHitZ", &_SimHitZ);

  duneobj->norm_s = 1.;
  double downsampling = 1; //default 1, set to eg. 0.01 for quick testing
  bool do_geometric_correction = false;
  duneobj->pot_s = pot/(downsampling*1e21);
  duneobj->nEvents = static_cast<int>(std::round(downsampling*static_cast<double>(_data->GetEntries())));

  // allocate memory for dunendgarmc variables
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 
  duneobj->geometric_correction = new double[duneobj->nEvents];

  duneobj->rw_Q0 = new double[duneobj->nEvents];
  duneobj->rw_Q3 = new double[duneobj->nEvents];
  duneobj->mode = new double[duneobj->nEvents];

  duneobj->ntruemuon = new int[duneobj->nEvents];
  duneobj->ntruemuonprim = new int[duneobj->nEvents];

  duneobj->in_fdv = new bool[duneobj->nEvents];
  duneobj->rw_elep_true = new double[duneobj->nEvents];

  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];
  duneobj->rw_rad = new double[duneobj->nEvents];

  duneobj->rw_lep_pT = new double[duneobj->nEvents];
  duneobj->rw_lep_pX = new double[duneobj->nEvents];
  duneobj->rw_lep_pY = new double[duneobj->nEvents];
  duneobj->rw_lep_pZ = new double[duneobj->nEvents];
  duneobj->rw_lep_p = new double[duneobj->nEvents];
  duneobj->rw_lep_theta = new double[duneobj->nEvents]; //angle to z (beam)
  duneobj->rw_lep_phi = new double[duneobj->nEvents]; //angle in xy plane to x (B-field)
  duneobj->rw_lep_bangle = new double[duneobj->nEvents];

  duneobj->Target = new int[duneobj->nEvents];

  duneobj->is_accepted = new bool[duneobj->nEvents];

  //Particle-level kinematic variables
  duneobj->particle_event = new std::vector<int>; duneobj->particle_event->reserve(7*duneobj->nEvents);
  duneobj->particle_pdg = new std::vector<int>; duneobj->particle_pdg->reserve(7*duneobj->nEvents); 
  duneobj->particle_energy = new std::vector<double>; duneobj->particle_energy->reserve(7*duneobj->nEvents); 
  duneobj->particle_momentum = new std::vector<double>; duneobj->particle_momentum->reserve(7*duneobj->nEvents); 
  duneobj->particle_endmomentum = new std::vector<double>; duneobj->particle_endmomentum->reserve(7*duneobj->nEvents); 
  duneobj->particle_transversemomentum = new std::vector<double>; duneobj->particle_transversemomentum->reserve(7*duneobj->nEvents); 
  duneobj->particle_bangle = new std::vector<double>; duneobj->particle_bangle->reserve(7*duneobj->nEvents);
  duneobj->particle_beamangle = new std::vector<double>; duneobj->particle_beamangle->reserve(7*duneobj->nEvents);
  duneobj->particle_isaccepted = new std::vector<bool>; duneobj->particle_isaccepted->reserve(7*duneobj->nEvents);
  duneobj->particle_iscurvatureresolved = new std::vector<bool>; duneobj->particle_iscurvatureresolved->reserve(7*duneobj->nEvents);
  duneobj->particle_isdecayed = new std::vector<bool>; duneobj->particle_isdecayed->reserve(7*duneobj->nEvents);
  duneobj->particle_isstoppedintpc = new std::vector<bool>; duneobj->particle_isstoppedintpc->reserve(7*duneobj->nEvents); 
  duneobj->particle_isstoppedinecal = new std::vector<bool>; duneobj->particle_isstoppedinecal->reserve(7*duneobj->nEvents); 
  duneobj->particle_isstoppedingap = new std::vector<bool>; duneobj->particle_isstoppedingap->reserve(7*duneobj->nEvents); 
  duneobj->particle_isstoppedinbarrelgap = new std::vector<bool>; duneobj->particle_isstoppedinbarrelgap->reserve(7*duneobj->nEvents); 
  duneobj->particle_isstoppedinendgap = new std::vector<bool>; duneobj->particle_isstoppedinendgap->reserve(7*duneobj->nEvents); 
  duneobj->particle_isstoppedinbarrel= new std::vector<bool>; duneobj->particle_isstoppedinbarrel->reserve(7*duneobj->nEvents);
  duneobj->particle_isstoppedinendcap = new std::vector<bool>; duneobj->particle_isstoppedinendcap->reserve(7*duneobj->nEvents); 
  duneobj->particle_isescaped = new std::vector<bool>; duneobj->particle_isescaped->reserve(7*duneobj->nEvents); 
  duneobj->particle_startx = new std::vector<double>; duneobj->particle_startx->reserve(7*duneobj->nEvents); 
  duneobj->particle_startr2 = new std::vector<double>; duneobj->particle_startr2->reserve(7*duneobj->nEvents); 
  duneobj->particle_endr = new std::vector<double>; duneobj->particle_endr->reserve(7*duneobj->nEvents); 
  duneobj->particle_enddepth = new std::vector<double>; duneobj->particle_enddepth->reserve(7*duneobj->nEvents); 
  duneobj->particle_endx = new std::vector<double>; duneobj->particle_endx->reserve(7*duneobj->nEvents); 
  duneobj->particle_endy = new std::vector<double>; duneobj->particle_endy->reserve(7*duneobj->nEvents); 
  duneobj->particle_endz = new std::vector<double>; duneobj->particle_endz->reserve(7*duneobj->nEvents); 
  duneobj->particle_nturns = new std::vector<double>; duneobj->particle_nturns->reserve(7*duneobj->nEvents); 
  duneobj->particle_nhits = new std::vector<double>; duneobj->particle_nhits->reserve(7*duneobj->nEvents); 
  duneobj->particle_tracklengthyz = new std::vector<double>; duneobj->particle_tracklengthyz->reserve(7*duneobj->nEvents); 
  duneobj->particle_momresms = new std::vector<double>; duneobj->particle_momresms->reserve(7*duneobj->nEvents); 
  duneobj->particle_momresyz = new std::vector<double>; duneobj->particle_momresyz->reserve(7*duneobj->nEvents);
  duneobj->particle_momresx = new std::vector<double>; duneobj->particle_momresx->reserve(7*duneobj->nEvents); 
  duneobj->particle_edepcrit = new std::vector<double>; duneobj->particle_edepcrit->reserve(7*duneobj->nEvents); 

  int numCC = 0;
  int num_in_fdv = 0;

  double pixel_spacing_cm = pixel_spacing/10; //convert to cm

  for (int i_event = 0; i_event < (duneobj->nEvents); ++i_event) { 
    _data->GetEntry(i_event);

    if (i_event % (duneobj->nEvents/100) == 0) {
      MACH3LOG_INFO("\tNow processing event: {}/{}",i_event,duneobj->nEvents);
    }
    double radius = std::sqrt((_MCVertY->at(0)-TPC_centre_y)*(_MCVertY->at(0)-TPC_centre_y) + (_MCVertZ->at(0)-TPC_centre_z)*(_MCVertZ->at(0)-TPC_centre_z)); //find radius of interaction vertex
    if(std::abs(_MCVertX->at(0) - TPC_centre_x)<=TPCFidLength &&  radius<=TPCFidRadius){
      num_in_fdv++;
      duneobj->in_fdv[i_event] = 1;
    } else{
      duneobj->in_fdv[i_event] = 0;
    }

    duneobj->rw_etru[i_event] = std::sqrt(_MCNuPx->at(0)*_MCNuPx->at(0) + _MCNuPy->at(0)*_MCNuPy->at(0) + _MCNuPz->at(0)*_MCNuPz->at(0)); // in GeV
    duneobj->rw_isCC[i_event] = static_cast<int>(!_IsNC->at(0));
    duneobj->rw_berpaacvwgt[i_event] = _BeRPA_cvwgt;

    duneobj->rw_vtx_x[i_event] = _MCVertX->at(0);
    duneobj->rw_vtx_y[i_event] = _MCVertY->at(0);
    duneobj->rw_vtx_z[i_event] = _MCVertZ->at(0);

    // Deal with truth-level information 
    std::unordered_map<int, std::vector<int>> mother_to_daughter_ID;
    std::unordered_map<int, size_t> ID_to_index;
    std::unordered_map<int, std::vector<double>> ID_to_ECalDep;
    const int tot_ecal_layers = 42;
    size_t n_ana_particles = _MCPTrkID->size();

    // Fill maps
    for (size_t i_anapart=0; i_anapart<n_ana_particles; i_anapart++) {
      int pdg = _PDG->at(i_anapart);
      if (pdg == 13) duneobj->ntruemuon[i_event]++;

      int secID = _MCPTrkID->at(i_anapart);
      int motherID = _MotherTrkID->at(i_anapart);
      ID_to_ECalDep[secID] = std::vector<double>(tot_ecal_layers,0.);
      ID_to_index[secID] = i_anapart;
      mother_to_daughter_ID[motherID].push_back(secID);
      mother_to_daughter_ID[secID]; // Ensure all particles are added to the map (even if no secondaries)
    }

    // Fill map from particle ID to ECal deposited energy
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

    for (int& primID : mother_to_daughter_ID[0]) {
      // Do not require the reconstruction of neutrons and neutrinos
      size_t i_anaprim = ID_to_index[primID];
      int pdg = _PDG->at(i_anaprim);
      if (pdg == 2112 || std::abs(pdg) == 12 || std::abs(pdg) == 14 || std::abs(pdg) == 16) continue;

      // Fill default values for particle-level parameters:
      duneobj->particle_event->push_back(i_event);
      duneobj->particle_pdg->push_back(M3::_BAD_INT_);
      duneobj->particle_energy->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_momentum->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_endmomentum->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_transversemomentum->push_back(M3::_BAD_DOUBLE_); // momentum transverse to B-field
      duneobj->particle_bangle->push_back(M3::_BAD_DOUBLE_); // angle to B-field
      duneobj->particle_beamangle->push_back(M3::_BAD_DOUBLE_);

      duneobj->particle_startx->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_startr2->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_endr->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_enddepth->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_endx->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_endy->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_endz->push_back(M3::_BAD_DOUBLE_);

      duneobj->particle_isaccepted->push_back(true); 
      duneobj->particle_iscurvatureresolved->push_back(false);
      duneobj->particle_isdecayed->push_back(false);
      duneobj->particle_isstoppedingap->push_back(false);
      duneobj->particle_isstoppedinbarrelgap->push_back(false);
      duneobj->particle_isstoppedinendgap->push_back(false);
      duneobj->particle_isstoppedinecal->push_back(false);
      duneobj->particle_isstoppedinbarrel->push_back(false);
      duneobj->particle_isstoppedinendcap->push_back(false);
      duneobj->particle_isstoppedintpc->push_back(false);
      duneobj->particle_isescaped->push_back(false);

      duneobj->particle_momresms->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_momresyz->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_momresx->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_nturns->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_nhits->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_tracklengthyz->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_edepcrit->push_back(CalcEDepCal(primID, mother_to_daughter_ID, ID_to_ECalDep, tot_ecal_layers, crit_layers));

      nparticlesinsample[iSample]++;

      bool isCurvatureResolved = false;
      // Remove descendants (and their descendants) from mother_to_daughter_ID who's momentum we get from curvature
      if (CurvatureResolutionFilter(primID, mother_to_daughter_ID, ID_to_index, duneobj, pixel_spacing_cm)) {
        isCurvatureResolved = true;
        duneobj->particle_iscurvatureresolved->back() = true;
      }

      // Find energy deposited by by primary and non-curvature-resolved descendants in critical region of calorimeter
      double EDepCrit = CalcEDepCal(primID, mother_to_daughter_ID, ID_to_ECalDep, tot_ecal_layers, crit_layers);
      double EDepTot = CalcEDepCal(primID, mother_to_daughter_ID, ID_to_ECalDep, tot_ecal_layers, tot_ecal_layers);
      bool isContained = true;
      if (EDepCrit > 0.002 && EDepCrit/EDepTot > 0.02) {
        isContained = false;
      }

      if (!(isContained || isCurvatureResolved)) {
        isEventAccepted = false;
        duneobj->particle_isaccepted->back() = false;
      }

      // Get primary muon information
      if(pdg == 13) {
        duneobj->ntruemuonprim[i_event]++;
        double p_x = _MCPStartPX->at(i_anaprim);
        double p_y = _MCPStartPY->at(i_anaprim);
        double p_z = _MCPStartPZ->at(i_anaprim);
        double p2 = p_x*p_x + p_y*p_y + p_z*p_z;
        if (p2 > muon_p2) {
          muon_p2 = p2;
          double mass = MaCh3Utils::GetMassFromPDG(pdg);
          duneobj->rw_lep_pX[i_event] = p_x;
          duneobj->rw_lep_pY[i_event] = p_y;
          duneobj->rw_lep_pZ[i_event] = p_z;
          duneobj->rw_elep_true[i_event] = std::sqrt(p2 + mass*mass);
        }
      }
    }

    if(isEventAccepted) {duneobj->is_accepted[i_event]=1;}
    else {duneobj->is_accepted[i_event]=0;} 

    double lep_momentum = std::sqrt(duneobj->rw_lep_pX[i_event]*duneobj->rw_lep_pX[i_event] + duneobj->rw_lep_pY[i_event]*duneobj->rw_lep_pY[i_event] + duneobj->rw_lep_pZ[i_event]*duneobj->rw_lep_pZ[i_event]);
    double lep_pBeam = duneobj->rw_lep_pY[i_event]*BeamDirection[1] + duneobj->rw_lep_pZ[i_event]*BeamDirection[2];
    double lep_pB = duneobj->rw_lep_pX[i_event];
    double lep_pPerp = duneobj->rw_lep_pY[i_event]*BeamDirection[2] - duneobj->rw_lep_pZ[i_event]*BeamDirection[1];

    double lep_beamangle = acos(lep_pBeam/lep_momentum)*180/M_PI; //Angle to beam (beam direction: [0.0,-0.101,0.995])
    double lep_bangle = acos(lep_pB/lep_momentum)*180/M_PI; //Angle to B-field (b-field along x)
    double lep_perpangle = acos(lep_pPerp/lep_momentum)*180/M_PI; //Angle to axis perpendicular to beam and B
    double lep_phi = atan2(lep_pPerp, lep_pB)*180/M_PI;

    duneobj->rw_lep_theta[i_event] = lep_beamangle;
    duneobj->rw_lep_phi[i_event] = lep_phi;
    duneobj->rw_lep_bangle[i_event] = lep_bangle;
    duneobj->rw_lep_p[i_event] = lep_momentum;

    duneobj->geometric_correction[i_event] = 1.;
    if (do_geometric_correction) {
      if ((lep_bangle < 45 || lep_bangle > 135) && lep_momentum > 0.3) duneobj->geometric_correction[i_event] = 0.;
      else if ((lep_perpangle < 45 || lep_perpangle > 135) && lep_momentum > 0.3) duneobj->geometric_correction[i_event] = 2.;
    }

    duneobj->rw_rad[i_event] = radius;

    //Assume everything is on Argon for now....
    duneobj->Target[i_event] = 40;

    duneobj->rw_Q0[i_event] = duneobj->rw_etru[i_event] - duneobj->rw_elep_true[i_event];
    duneobj->rw_Q3[i_event] = std::sqrt((_MCNuPx->at(0)-duneobj->rw_lep_pX[i_event])*(_MCNuPx->at(0)-duneobj->rw_lep_pX[i_event]) + 
                                        (_MCNuPy->at(0)-duneobj->rw_lep_pY[i_event])*(_MCNuPy->at(0)-duneobj->rw_lep_pY[i_event]) + 
                                        (_MCNuPz->at(0)-duneobj->rw_lep_pZ[i_event])*(_MCNuPz->at(0)-duneobj->rw_lep_pZ[i_event]));
    duneobj->rw_lep_pT[i_event] = std::sqrt(lep_momentum*lep_momentum - lep_pBeam*lep_pBeam); 

    duneobj->mode[i_event] = _MCMode->at(0);

    duneobj->flux_w[i_event] = 1.0;
    if(duneobj->rw_isCC[i_event] == 1) numCC++;
  }
  MACH3LOG_INFO("nEvents = {}, numCC = {}, numFDV = {}", duneobj->nEvents, numCC, num_in_fdv);

  _sampleFile->Close();
  return duneobj->nEvents;
}

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wswitch-enum"
const double* SampleHandlerBeamNDGAr::GetPointerToKinematicParameter(KinematicTypes KinematicParameter, int iSample, int iEvent) {
  switch(KinematicParameter) {
    case kTrueNeutrinoEnergy:
      return &dunendgarmcSamples[iSample].rw_etru[iEvent]; 
    case kMode:
      return &dunendgarmcSamples[iSample].mode[iEvent];
    case kTrueXPos:
      return &dunendgarmcSamples[iSample].rw_vtx_x[iEvent];
    case kTrueYPos:
      return &dunendgarmcSamples[iSample].rw_vtx_y[iEvent];
    case kTrueZPos:
      return &dunendgarmcSamples[iSample].rw_vtx_z[iEvent];
    case kTrueRad:
      return &dunendgarmcSamples[iSample].rw_rad[iEvent];
    case kTrueLepEnergy:
      return &dunendgarmcSamples[iSample].rw_elep_true[iEvent];
    case kLepPT:
      return &dunendgarmcSamples[iSample].rw_lep_pT[iEvent];
    case kLepPZ:
      return &dunendgarmcSamples[iSample].rw_lep_pZ[iEvent];
    case kLepTheta:
      return &dunendgarmcSamples[iSample].rw_lep_theta[iEvent];
    case kLepPhi:
      return &dunendgarmcSamples[iSample].rw_lep_phi[iEvent];
    case kLepBAngle:
      return &dunendgarmcSamples[iSample].rw_lep_bangle[iEvent];
    case kLepP:
      return &dunendgarmcSamples[iSample].rw_lep_p[iEvent];
    case kTrueQ0:
      return &dunendgarmcSamples[iSample].rw_Q0[iEvent];
    case kTrueQ3:
      return &dunendgarmcSamples[iSample].rw_Q3[iEvent];
    case kParticle_Momentum:
      return &dunendgarmcSamples[iSample].particle_momentum->at(iEvent);
    case kParticle_EndMomentum:
      return &dunendgarmcSamples[iSample].particle_endmomentum->at(iEvent);
    case kParticle_TransverseMomentum:
      return &dunendgarmcSamples[iSample].particle_transversemomentum->at(iEvent);
    case kParticle_BAngle:
      return &dunendgarmcSamples[iSample].particle_bangle->at(iEvent);
    case kParticle_BeamAngle:
      return &dunendgarmcSamples[iSample].particle_beamangle->at(iEvent);
    case kParticle_NTurns:
      return &dunendgarmcSamples[iSample].particle_nturns->at(iEvent);
    case kParticle_NHits:
      return &dunendgarmcSamples[iSample].particle_nhits->at(iEvent);
    case kParticle_TrackLengthYZ:
      return &dunendgarmcSamples[iSample].particle_tracklengthyz->at(iEvent);
    case kParticle_MomResMS:
      return &dunendgarmcSamples[iSample].particle_momresms->at(iEvent);
    case kParticle_MomResYZ:
      return &dunendgarmcSamples[iSample].particle_momresyz->at(iEvent);
    case kParticle_MomResX:
      return &dunendgarmcSamples[iSample].particle_momresx->at(iEvent);
    case kParticle_StartR2:
      return &dunendgarmcSamples[iSample].particle_startr2->at(iEvent);
    case kParticle_EndR:
      return &dunendgarmcSamples[iSample].particle_endr->at(iEvent);
    case kParticle_EndDepth:
      return &dunendgarmcSamples[iSample].particle_enddepth->at(iEvent);
    case kParticle_EndX:
      return &dunendgarmcSamples[iSample].particle_endx->at(iEvent);
    case kParticle_EndY:
      return &dunendgarmcSamples[iSample].particle_endy->at(iEvent);
    case kParticle_EndZ:
      return &dunendgarmcSamples[iSample].particle_endz->at(iEvent);
    case kParticle_StartX:
      return &dunendgarmcSamples[iSample].particle_startx->at(iEvent);
    case kParticle_EDepCrit:
      return &dunendgarmcSamples[iSample].particle_edepcrit->at(iEvent);
    default:
      MACH3LOG_ERROR("Did not recognise Kinematic Parameter {}", static_cast<int>(KinematicParameter));
      throw MaCh3Exception(__FILE__, __LINE__);
      return nullptr;
  }
}
#pragma GCC diagnostic pop

const double* SampleHandlerBeamNDGAr::GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(std::round(KinematicVariable));
  return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

const double* SampleHandlerBeamNDGAr::GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

double SampleHandlerBeamNDGAr::ReturnKinematicParameter(int KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(std::round(KinematicVariable));
  return ReturnKinematicParameter(KinPar,iSample,iEvent);
}

double SampleHandlerBeamNDGAr::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return ReturnKinematicParameter(KinPar,iSample,iEvent);
}

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wswitch-enum"
double SampleHandlerBeamNDGAr::ReturnKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent) {
  //HH: Special cases for dealing with non-doubles
  switch(KinPar) {
    case kEvent_IsAccepted:
      return static_cast<double>(dunendgarmcSamples[iSample].is_accepted[iEvent]);
    case kIsCC:
      return static_cast<double>(dunendgarmcSamples[iSample].rw_isCC[iEvent]);
    case kInFDV:
      return static_cast<double>(dunendgarmcSamples[iSample].in_fdv[iEvent]);
    case kParticle_Event:
      return static_cast<double>(dunendgarmcSamples[iSample].particle_event->at(iEvent));
    case kParticle_IsAccepted:
      return static_cast<double>(dunendgarmcSamples[iSample].particle_isaccepted->at(iEvent));
    case kParticle_IsCurvatureResolved:
      return static_cast<double>(dunendgarmcSamples[iSample].particle_iscurvatureresolved->at(iEvent));
    case kParticle_IsDecayed:
      return static_cast<double>(dunendgarmcSamples[iSample].particle_isdecayed->at(iEvent));
    case kParticle_PDG:
      return static_cast<double>(dunendgarmcSamples[iSample].particle_pdg->at(iEvent));
    case kParticle_IsStoppedInTPC:
      return static_cast<double>(dunendgarmcSamples[iSample].particle_isstoppedintpc->at(iEvent));
    case kParticle_IsStoppedInECal:
      return static_cast<double>(dunendgarmcSamples[iSample].particle_isstoppedinecal->at(iEvent));
    case kParticle_IsStoppedInGap:
      return static_cast<double>(dunendgarmcSamples[iSample].particle_isstoppedingap->at(iEvent));
    case kParticle_IsStoppedInEndGap:
      return static_cast<double>(dunendgarmcSamples[iSample].particle_isstoppedinendgap->at(iEvent));
    case kParticle_IsStoppedInBarrelGap:
      return static_cast<double>(dunendgarmcSamples[iSample].particle_isstoppedinbarrelgap->at(iEvent));
    case kParticle_IsEscaped:
      return static_cast<double>(dunendgarmcSamples[iSample].particle_isescaped->at(iEvent));
    default:
      return *GetPointerToKinematicParameter(KinPar, iSample, iEvent);
  }
}
#pragma GCC diagnostic pop

void SampleHandlerBeamNDGAr::SetupFDMC(int iSample) {
  dunemc_base *duneobj = &(dunendgarmcSamples[iSample]);
  FarDetectorCoreInfo *fdobj = &(MCSamples[iSample]);

  //fdobj->nutype = duneobj->nutype;
  //fdobj->oscnutype = duneobj->oscnutype;
  //fdobj->signal = duneobj->signal;
  //fdobj->SampleName = SampleName;

  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent){
    fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
    fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
    fdobj->Target[iEvent] = &(duneobj->Target[iEvent]);
  }
}

std::vector<double> SampleHandlerBeamNDGAr::ReturnKinematicParameterBinning(std::string KinematicParameterStr) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameterStr));
  return ReturnKinematicParameterBinning(KinPar);
}

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wswitch-enum"
std::vector<double> SampleHandlerBeamNDGAr::ReturnKinematicParameterBinning(KinematicTypes KinPar) {
  std::vector<double> binningVector;
  switch(KinPar){
    case kTrueNeutrinoEnergy:
      for(double ibins =0; ibins<10*10; ibins++){
        double binval = ibins/10;
        binningVector.push_back(binval);
      }
      break;
    case kTrueXPos:
      for(double ibins =0; ibins<259*2; ibins++){
        binningVector.push_back(ibins-259);
      }
      break;
    case kTrueYPos:
      for(double ibins =0; ibins<277*2; ibins++){
        binningVector.push_back(ibins-277-150);
      }
      break;
    case kTrueZPos:
      for(double ibins =0; ibins<277*2; ibins++){
        binningVector.push_back(ibins-277+1486);
      }
      break;
    case kTrueLepEnergy:
      for(double ibins =0; ibins<10*10; ibins++){
        binningVector.push_back(ibins/10);
      } 
      break;
    case kTrueRad:
      for(double ibins =0; ibins<300; ibins++){
        binningVector.push_back(ibins);
      }
      break;
    case kLepPT:
    case kLepPZ:
      for(double ibins =0; ibins<10*10; ibins++){
        binningVector.push_back(ibins/10);
      }
      break;
    case kTrueQ0:
    case kTrueQ3:
      for(double ibins =0; ibins<3*50+1; ibins++){
        binningVector.push_back(ibins/50.);
      }
      break;
    case kMode:
      for (int ibins=0; ibins<Modes->GetNModes(); ibins++) {
        binningVector.push_back(ibins);
      }
      break;
    default:
      for(double ibins =0; ibins<10*100; ibins++){
        binningVector.push_back(ibins/100);
      }
      break;
  }
  return binningVector;
}
#pragma GCC diagnostic pop

TH1* SampleHandlerBeamNDGAr::Get1DParticleVarHist(std::string ProjectionVar_Str, std::vector< KinematicCut > SelectionVec, int WeightStyle, TAxis* Axis) {
  //DB Grab the associated enum with the argument string
  int ProjectionVar_Int = ReturnKinematicParameterFromString(ProjectionVar_Str);

  //DB Need to overwrite the Selection member variable so that IsEventSelected function operates correctly.
  //   Consequently, store the selection cuts already saved in the sample, overwrite the Selection variable, then reset
  std::vector< KinematicCut > tmp_Selection = Selection;
  std::vector< KinematicCut > SelectionVecToApply;

  //DB Add all the predefined selections to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<Selection.size();iSelec++) {
    SelectionVecToApply.emplace_back(Selection[iSelec]);
  }

  //DB Add all requested cuts from the argument to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<SelectionVec.size();iSelec++) {
    SelectionVecToApply.emplace_back(SelectionVec[iSelec]);
  }

  //DB Set the member variable to be the cuts to apply
  Selection = SelectionVecToApply;

  //DB Define the histogram which will be returned
  TH1D* _h1DVar;
  if (Axis) {
    _h1DVar = new TH1D("","",Axis->GetNbins(),Axis->GetXbins()->GetArray());
  } else {
    std::vector<double> xBinEdges = ReturnKinematicParameterBinning(ProjectionVar_Str);
    _h1DVar = new TH1D("", "", int(xBinEdges.size())-1, xBinEdges.data());
  }

  //Loop over all primary particles
  for (int iSample=0;iSample<GetNMCSamples();iSample++) {
    for (int iParticle=0;iParticle<nparticlesinsample[iSample];iParticle++) {
      int event = static_cast<int>(ReturnKinematicParameter("Particle_Event", iSample, iParticle));
      if (IsParticleSelected(iSample,event,iParticle)) {
        double Weight = GetEventWeight(iSample,event);
        if (WeightStyle==1) {
          Weight = 1.;
        }
        double Var = ReturnKinematicParameter(ProjectionVar_Int,iSample,iParticle);
        _h1DVar->Fill(Var,Weight);
      }
    }
  }

  //DB Reset the saved selection
  Selection = tmp_Selection;

  return _h1DVar;
}

TH2* SampleHandlerBeamNDGAr::Get2DParticleVarHist(std::string ProjectionVar_StrX, std::string ProjectionVar_StrY, std::vector< KinematicCut > SelectionVec, int WeightStyle, TAxis* AxisX, TAxis* AxisY) {
  //DB Grab the associated enum with the argument string
  int ProjectionVar_IntX = ReturnKinematicParameterFromString(ProjectionVar_StrX);
  int ProjectionVar_IntY = ReturnKinematicParameterFromString(ProjectionVar_StrY);

  //DB Need to overwrite the Selection member variable so that IsEventSelected function operates correctly.
  //   Consequently, store the selection cuts already saved in the sample, overwrite the Selection variable, then reset
  std::vector< KinematicCut > tmp_Selection = Selection;
  std::vector< KinematicCut > SelectionVecToApply;

  //DB Add all the predefined selections to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<Selection.size();iSelec++) {
    SelectionVecToApply.emplace_back(Selection[iSelec]);
  }

  //DB Add all requested cuts from the argument to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<SelectionVec.size();iSelec++) {
    SelectionVecToApply.emplace_back(SelectionVec[iSelec]);
  }

  //DB Set the member variable to be the cuts to apply
  Selection = SelectionVecToApply;

  //DB Define the histogram which will be returned
  TH2D* _h2DVar;
  if (AxisX && AxisY) {
    _h2DVar = new TH2D("","",AxisX->GetNbins(),AxisX->GetXbins()->GetArray(),AxisY->GetNbins(),AxisY->GetXbins()->GetArray());
  } else {
    std::vector<double> xBinEdges = ReturnKinematicParameterBinning(ProjectionVar_StrX);
    std::vector<double> yBinEdges = ReturnKinematicParameterBinning(ProjectionVar_StrY);
    _h2DVar = new TH2D("", "", int(xBinEdges.size())-1, xBinEdges.data(), int(yBinEdges.size())-1, yBinEdges.data());
  }

  //JM Loop over all particles
  for (int iSample=0;iSample<GetNMCSamples();iSample++) {
    for (int iParticle=0;iParticle<nparticlesinsample[iSample];iParticle++) {
      int event = static_cast<int>(ReturnKinematicParameter("Particle_Event", iSample, iParticle));
      if (IsParticleSelected(iSample,event,iParticle)) {
        double Weight = GetEventWeight(iSample,event);
        if (WeightStyle==1) {
          Weight = 1.;
        }
        double VarX = ReturnKinematicParameter(ProjectionVar_IntX,iSample,iParticle);
        double VarY = ReturnKinematicParameter(ProjectionVar_IntY,iSample,iParticle);
        _h2DVar->Fill(VarX,VarY,Weight);
      }
    }
  }

  //DB Reset the saved selection
  Selection = tmp_Selection;

  return _h2DVar;
}

bool SampleHandlerBeamNDGAr::IsParticleSelected(const int iSample, const int iEvent, const int iParticle) {
  double Val;
  for (unsigned int iSelection=0;iSelection < Selection.size() ;iSelection++) {
    std::string SelecVarStr = ReturnStringFromKinematicParameter(Selection[iSelection].ParamToCutOnIt);
    int SelectionVarElement;

    if (SelecVarStr.find("Particle_") != std::string::npos) {SelectionVarElement = iParticle;}
    else {SelectionVarElement = iEvent;}

    Val = ReturnKinematicParameter(Selection[iSelection].ParamToCutOnIt, iSample, SelectionVarElement);
    if ((Val<Selection[iSelection].LowerBound)||(Val>=Selection[iSelection].UpperBound)) {
      return false;
    }
  }
  return true;
}
