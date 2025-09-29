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
          return 0;
          // throw MaCh3Exception(__FILE__,__LINE__);
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
    int motherID = _MCPMotherTrkID->at(index);
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

// Returns depth of a coordinate position in the dodecoganol ecal, or +/- 999 if it is beyond or within the ecal boundaries
double SampleHandlerBeamNDGAr::GetCalDepth(double x, double y, double z) {
  x = x - TPC_centre_x;
  y = y - TPC_centre_y;
  z = z - TPC_centre_z;
  if (std::abs(x) > ECALEndCapEnd) return 999.; // beyond ecal lengthways

  double theta = atan2(y, z) + M_PI + M_PI/12.;
  double remainder = fmod(theta, M_PI/6.);
  double theta_segment;
  if (remainder < M_PI/12.) theta_segment = theta - remainder;
  else theta_segment = theta + M_PI/6. - remainder;

  double r = std::sqrt(y*y + z*z);
  double projected_r = r*cos(theta_segment - theta);

  if (projected_r > ECALOuterRadius) return 999.; // beyond ecal radially
  else if (projected_r >= ECALInnerRadius) return projected_r - ECALInnerRadius; // in barrel
  else if (std::abs(x) >= ECALEndCapStart && r <= ECALInnerRadius) return std::abs(x) - ECALEndCapStart; // in end cap
  else return -999.; // before ecal
}

double SampleHandlerBeamNDGAr::CalcEDepCal(int trkID, const std::unordered_map<int, std::vector<int>>& mother_to_daughter_ID, const std::unordered_map<int, std::vector<std::vector<double>>>& ID_to_ECalDep, double crit_reg, std::unordered_map<int, double>& ecal_layer_deposits, bool store_deposits = false) {
  double EDepCrit = 0.;
  auto it = mother_to_daughter_ID.find(trkID);
  if (it != mother_to_daughter_ID.end()) {
    if(ID_to_ECalDep.find(trkID) != ID_to_ECalDep.end()) {
      for (const auto& calhit : ID_to_ECalDep.at(trkID)) {
        if (ECALOuterRadius - ECALInnerRadius - calhit[1] < crit_reg) EDepCrit += calhit[0];
        if (store_deposits) {
          int layer = static_cast<int>(calhit[1]/1.5);
          ecal_layer_deposits[layer] += calhit[0];
        }
      }
    }
    for (int daughterID : it->second) {
      EDepCrit += CalcEDepCal(daughterID, mother_to_daughter_ID, ID_to_ECalDep, crit_reg, ecal_layer_deposits, store_deposits);
    }
  }
  return EDepCrit;
}

bool SampleHandlerBeamNDGAr::IsResolvedFromCurvature(dunemc_base *duneobj, int i_particle, double pixel_spacing_cm){
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
  double end_mom = std::sqrt(_MCPEndPX->at(i_particle)*_MCPEndPX->at(i_particle) + _MCPEndPY->at(i_particle)*_MCPEndPY->at(i_particle) + _MCPEndPZ->at(i_particle)*_MCPEndPZ->at(i_particle))/1000.; 

  double yend_centre = _MCPEndY->at(i_particle)-TPC_centre_y;
  double zend_centre = _MCPEndZ->at(i_particle)-TPC_centre_z;
  double end_radius = std::sqrt(yend_centre*yend_centre + zend_centre*zend_centre); 
  double end_length = _MCPEndX->at(i_particle)-TPC_centre_x;
  double ecal_depth = GetCalDepth(_MCPEndX->at(i_particle), _MCPEndY->at(i_particle), _MCPEndZ->at(i_particle));

  bool stops_in_tpc = std::abs(end_length)<=TPCInstrumentedLength && end_radius<=TPCInstrumentedRadius;
  bool stops_before_ecal = ecal_depth == -999.;
  bool stops_beyond_ecal = ecal_depth == 999.;
  bool stops_in_ecal = !(stops_before_ecal || stops_beyond_ecal);

  //Fill particle-level kinematic variables with default or actual (if possible at this stage) values
  if (_MCPMotherTrkID->at(i_particle) == 0) {
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
    duneobj->particle_endx->back() = end_length;
    duneobj->particle_endy->back() = yend_centre;
    duneobj->particle_endz->back() = zend_centre;
    duneobj->particle_ecaldepth->back() = ecal_depth;

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
  double theta_xT = atan2(p_x,transverse_mom); //helix pitch angle
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

    //Find Position where track leaves TPC. Intersection of two circles.
    double m_const = (TPC_centre_z - centre_circle_z)/(TPC_centre_y-centre_circle_y); //gradient of line between two intersection points
    double a_const = (TPCInstrumentedRadius*TPCInstrumentedRadius-rad_curvature*rad_curvature - (TPC_centre_y*TPC_centre_y - centre_circle_y*centre_circle_y)-(TPC_centre_z*TPC_centre_z - centre_circle_z*centre_circle_z))/(2*(centre_circle_y-TPC_centre_y));
    double quadraticformula_b = -(2*m_const*(a_const - TPC_centre_y) + 2*TPC_centre_z);
    double quadraticformula_a = m_const*m_const + 1;
    double quadraticformula_c = (a_const - TPC_centre_y)*(a_const - TPC_centre_y) + TPC_centre_z*TPC_centre_z - TPCInstrumentedRadius*TPCInstrumentedRadius;


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

  if (_MCPMotherTrkID->at(i_particle) == 0) {
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
  ecal_deposits.open("CalorimeterDeposits.txt");

  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("Input File: {}", mc_files.at(iSample));

  // Read fastgarsim output file
  std::string simFileStr = mc_files.at(iSample);
  simFile = TFile::Open(simFileStr.c_str(), "READ");
  simTree = static_cast<TTree*>(simFile->Get("AnaTree"));

  if(simTree){
    MACH3LOG_INFO("Found anatree in {}", mc_files[iSample]);
    MACH3LOG_INFO("With number of entries: {}", simTree->GetEntries());
  }
  else{
    MACH3LOG_ERROR("Could not find anatree in {}", mc_files[iSample]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // Read Genie input file
  std::string genieFileStr = "";
  if (simFileStr.find("hA") != std::string::npos) genieFileStr = "Inputs/DUNE_NDGAr_FastGarSim/GenieTrees/numu_argon_G18_10a_gst.root";
  else if (simFileStr.find("hN") != std::string::npos) genieFileStr = "Inputs/DUNE_NDGAr_FastGarSim/GenieTrees/numu_argon_G18_10b_gst.root";
  genieFile = TFile::Open(genieFileStr.c_str(), "READ");
  genieTree = static_cast<TTree*>(genieFile->Get("gst"));

  if(genieTree){
    MACH3LOG_INFO("Found Genie tree in {}", genieFileStr);
    MACH3LOG_INFO("With number of entries: {}", genieTree->GetEntries());
  }
  else{
    MACH3LOG_ERROR("Could not find Genie tree in {}", genieFileStr);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // Switching of z and x to ensure x is B-field and z is (approx) beam
  simTree->SetBranchAddress("eventID", &_EventID);
  simTree->SetBranchAddress("startZ", &_MCPStartX);
  simTree->SetBranchAddress("startY", &_MCPStartY);
  simTree->SetBranchAddress("startX", &_MCPStartZ);
  simTree->SetBranchAddress("endZ", &_MCPEndX);
  simTree->SetBranchAddress("endY", &_MCPEndY);
  simTree->SetBranchAddress("endX", &_MCPEndZ);
  simTree->SetBranchAddress("startPZ", &_MCPStartPX);
  simTree->SetBranchAddress("startPY", &_MCPStartPY);
  simTree->SetBranchAddress("startPX", &_MCPStartPZ);
  simTree->SetBranchAddress("endPZ", &_MCPEndPX);
  simTree->SetBranchAddress("endPY", &_MCPEndPY);
  simTree->SetBranchAddress("endPX", &_MCPEndPZ);
  simTree->SetBranchAddress("pdgCode", &_MCPPDG);
  simTree->SetBranchAddress("trackID", &_MCPTrkID);
  simTree->SetBranchAddress("motherID", &_MCPMotherTrkID);
  simTree->SetBranchAddress("tpcHitTrackID", &_TPCHitTrkID);
  simTree->SetBranchAddress("tpcHitEdep", &_TPCHitEnergy);
  simTree->SetBranchAddress("tpcHitZ", &_TPCHitX);
  simTree->SetBranchAddress("tpcHitY", &_TPCHitY);
  simTree->SetBranchAddress("tpcHitX", &_TPCHitZ);
  simTree->SetBranchAddress("ecalHitTrackID", &_CalHitTrkID);
  simTree->SetBranchAddress("ecalHitEdep", &_CalHitEnergy);
  simTree->SetBranchAddress("ecalHitZ", &_CalHitX);
  simTree->SetBranchAddress("ecalHitY", &_CalHitY);
  simTree->SetBranchAddress("ecalHitX", &_CalHitZ);
  simTree->SetBranchAddress("muidHitTrackID", &_MuIDHitTrkID);
  simTree->SetBranchAddress("muidHitEdep", &_MuIDHitEnergy);
  simTree->SetBranchAddress("muidHitZ", &_MuIDHitX);
  simTree->SetBranchAddress("muidHitY", &_MuIDHitY);
  simTree->SetBranchAddress("muidHitX", &_MuIDHitZ);

  genieTree->SetBranchAddress("Ev", &_Enu);
  genieTree->SetBranchAddress("pxv", &_PXnu);
  genieTree->SetBranchAddress("pyv", &_PYnu);
  genieTree->SetBranchAddress("pzv", &_PZnu);
  genieTree->SetBranchAddress("neu", &_nuPDG);
  genieTree->SetBranchAddress("El", &_Elep);
  genieTree->SetBranchAddress("pxl", &_PXlep);
  genieTree->SetBranchAddress("pyl", &_PYlep);
  genieTree->SetBranchAddress("pzl", &_PZlep);
  genieTree->SetBranchAddress("cc", &_isCC);
  genieTree->SetBranchAddress("nfpip", &_npip);
  genieTree->SetBranchAddress("nfpim", &_npim);
  genieTree->SetBranchAddress("nfpi0", &_npi0);

  duneobj->norm_s = 1.;
  double downsampling = 1.; //default 1, set to eg. 0.01 for quick testing
  bool do_geometric_correction = false;
  duneobj->pot_s = pot/(downsampling*1e21);
  duneobj->nEvents = static_cast<int>(std::round(downsampling*static_cast<double>(simTree->GetEntries())));

  // allocate memory for dunendgarmc variables
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_nuPDG = new int[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 
  duneobj->geometric_correction = new double[duneobj->nEvents];

  duneobj->rw_Q0 = new double[duneobj->nEvents];
  duneobj->rw_Q3 = new double[duneobj->nEvents];

  duneobj->nproton = new int[duneobj->nEvents];
  duneobj->nneutron = new int[duneobj->nEvents];
  duneobj->npip = new int[duneobj->nEvents];
  duneobj->npim = new int[duneobj->nEvents];
  duneobj->npi0 = new int[duneobj->nEvents];

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
  int avg_particles_per_event = 30;
  duneobj->particle_event = new std::vector<int>; duneobj->particle_event->reserve(avg_particles_per_event*duneobj->nEvents);
  duneobj->particle_pdg = new std::vector<int>; duneobj->particle_pdg->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_energy = new std::vector<double>; duneobj->particle_energy->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_momentum = new std::vector<double>; duneobj->particle_momentum->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_endmomentum = new std::vector<double>; duneobj->particle_endmomentum->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_transversemomentum = new std::vector<double>; duneobj->particle_transversemomentum->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_bangle = new std::vector<double>; duneobj->particle_bangle->reserve(avg_particles_per_event*duneobj->nEvents);
  duneobj->particle_beamangle = new std::vector<double>; duneobj->particle_beamangle->reserve(avg_particles_per_event*duneobj->nEvents);
  duneobj->particle_isaccepted = new std::vector<bool>; duneobj->particle_isaccepted->reserve(avg_particles_per_event*duneobj->nEvents);
  duneobj->particle_iscurvatureresolved = new std::vector<bool>; duneobj->particle_iscurvatureresolved->reserve(avg_particles_per_event*duneobj->nEvents);
  duneobj->particle_isstoppedintpc = new std::vector<bool>; duneobj->particle_isstoppedintpc->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_isstoppedinecal = new std::vector<bool>; duneobj->particle_isstoppedinecal->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_isstoppedingap = new std::vector<bool>; duneobj->particle_isstoppedingap->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_isstoppedinbarrelgap = new std::vector<bool>; duneobj->particle_isstoppedinbarrelgap->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_isstoppedinendgap = new std::vector<bool>; duneobj->particle_isstoppedinendgap->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_isstoppedinbarrel= new std::vector<bool>; duneobj->particle_isstoppedinbarrel->reserve(avg_particles_per_event*duneobj->nEvents);
  duneobj->particle_isstoppedinendcap = new std::vector<bool>; duneobj->particle_isstoppedinendcap->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_isescaped = new std::vector<bool>; duneobj->particle_isescaped->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_startx = new std::vector<double>; duneobj->particle_startx->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_startr2 = new std::vector<double>; duneobj->particle_startr2->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_endr = new std::vector<double>; duneobj->particle_endr->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_endx = new std::vector<double>; duneobj->particle_endx->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_endy = new std::vector<double>; duneobj->particle_endy->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_endz = new std::vector<double>; duneobj->particle_endz->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_ecaldepth = new std::vector<double>; duneobj->particle_ecaldepth->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_nturns = new std::vector<double>; duneobj->particle_nturns->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_nhits = new std::vector<double>; duneobj->particle_nhits->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_tracklengthyz = new std::vector<double>; duneobj->particle_tracklengthyz->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_momresms = new std::vector<double>; duneobj->particle_momresms->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_momresyz = new std::vector<double>; duneobj->particle_momresyz->reserve(avg_particles_per_event*duneobj->nEvents);
  duneobj->particle_momresx = new std::vector<double>; duneobj->particle_momresx->reserve(avg_particles_per_event*duneobj->nEvents); 
  duneobj->particle_edepcrit = new std::vector<double>; duneobj->particle_edepcrit->reserve(avg_particles_per_event*duneobj->nEvents); 

  int numCC = 0;
  int num_in_fdv = 0;

  double pixel_spacing_cm = pixel_spacing/10; //convert to cm

  for (int i_event = 0; i_event < (duneobj->nEvents); ++i_event) { 
    simTree->GetEntry(i_event);
    genieTree->GetEntry(_EventID);
    if (i_event % (duneobj->nEvents/100) == 0) {
      MACH3LOG_INFO("\tNow processing event: {}/{}",i_event,duneobj->nEvents);
    }
    std::vector<double> vertex = {M3::_BAD_DOUBLE_, M3::_BAD_DOUBLE_, M3::_BAD_DOUBLE_};

    duneobj->rw_etru[i_event] = _Enu; // in GeV
    duneobj->rw_isCC[i_event] = _isCC;
    duneobj->rw_nuPDG[i_event] = _nuPDG;
    duneobj->rw_berpaacvwgt[i_event] = _BeRPA_cvwgt;

    // Deal with truth-level information 
    std::unordered_map<int, std::vector<int>> mother_to_daughter_ID;
    std::unordered_map<int, size_t> ID_to_index;
    std::unordered_map<int, std::vector<std::vector<double>>> ID_to_ECalDep; // Each deposition is vector of size 3: [energy, radius, x]
    std::unordered_map<int, double> ID_to_TPCDep;  // Each deposition is vector of size 3: [energy, radius, x]
    std::unordered_map<int, double> ID_to_MuIDDep; // Each deposition is vector of size 3: [energy, radius, x]
    size_t n_particles = _MCPTrkID->size();

    // Fill maps
    for (size_t i_particle=0; i_particle<n_particles; i_particle++) {
      int pdg = _MCPPDG->at(i_particle);
      if (pdg == 13) duneobj->ntruemuon[i_event]++;

      int trkID = _MCPTrkID->at(i_particle);
      int motherID = _MCPMotherTrkID->at(i_particle);
      ID_to_index[trkID] = i_particle;
      mother_to_daughter_ID[motherID].push_back(trkID);
      mother_to_daughter_ID[trkID]; //Ensure all particles are added to the map (even if no secondaries)
    }

    // Fill map from particle ID to ECal deposited energy
    for (size_t i_calhit=0; i_calhit<_CalHitTrkID->size(); i_calhit++) {
      int trkid = _CalHitTrkID->at(i_calhit);
      if (trkid <= 0) continue;
      double dep_energy = _CalHitEnergy->at(i_calhit)/1000.;
      double dep_depth = GetCalDepth(_CalHitX->at(i_calhit), _CalHitY->at(i_calhit), _CalHitZ->at(i_calhit));
      if (std::abs(dep_depth) == 999.) MACH3LOG_ERROR("Recorded calorimeter hit outside the calorimeter");
      ID_to_ECalDep[trkid].push_back({dep_energy, dep_depth});
    }

    // Fill map from particle ID to TPC deposited energy
    for (size_t i_tpchit=0; i_tpchit<_TPCHitTrkID->size(); i_tpchit++) {
      int trkid = _TPCHitTrkID->at(i_tpchit);
      if (trkid <= 0) continue;
      double dep_energy = _TPCHitEnergy->at(i_tpchit)/1000.;
      ID_to_TPCDep[trkid] = dep_energy;
    }

    // Fill map from particle ID to MuID deposited energy
    for (size_t i_muidhit=0; i_muidhit<_MuIDHitTrkID->size(); i_muidhit++) {
      int trkid = _MuIDHitTrkID->at(i_muidhit);
      if (trkid <= 0) continue;
      double dep_energy = _MuIDHitEnergy->at(i_muidhit)/1000.;
      ID_to_MuIDDep[trkid] = dep_energy;
    }

    double muon_p2 = 0.;
    bool isEventAccepted = true;
    double crit_reg = 3.; // Outer number of cm defining the calorimeter's 'critical region'

    for (int& primID : mother_to_daughter_ID[0]) {
      // Do not require the reconstruction of neutrons and neutrinos
      size_t prim_index = ID_to_index[primID];
      int pdg = _MCPPDG->at(prim_index);
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
      duneobj->particle_endx->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_endy->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_endz->push_back(M3::_BAD_DOUBLE_);
      duneobj->particle_ecaldepth->push_back(M3::_BAD_DOUBLE_);

      duneobj->particle_isaccepted->push_back(true); 
      duneobj->particle_iscurvatureresolved->push_back(false);
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

      ecal_deposits << "\n" << primID << ")" << std::endl;
      std::unordered_map<int, double> ecal_layer_deposits; 
      duneobj->particle_edepcrit->push_back(CalcEDepCal(primID, mother_to_daughter_ID, ID_to_ECalDep, crit_reg, ecal_layer_deposits, true));
      for (int i_layer=0; i_layer<42; i_layer++) {
        ecal_deposits << "Energy: " << ecal_layer_deposits[i_layer] << " Depth: " << i_layer << std::endl;
      }

      nparticlesinsample[iSample]++;

      bool isCurvatureResolved = false;
      // Remove descendants (and their descendants) from mother_to_daughter_ID who's momentum we get from curvature
      if (CurvatureResolutionFilter(primID, mother_to_daughter_ID, ID_to_index, duneobj, pixel_spacing_cm)) {
        isCurvatureResolved = true;
      }
      duneobj->particle_iscurvatureresolved->back() = isCurvatureResolved;

      ecal_deposits << "Primary: PDG = " << pdg << " Energy = " << duneobj->particle_energy->back() << " TPC = " << duneobj->particle_isstoppedintpc->back() <<
        " Gap = " << duneobj->particle_isstoppedingap->back() << " ECal = " << duneobj->particle_isstoppedinecal->back() << std::endl;

      // Find energy deposited by by primary and non-curvature-resolved descendants in critical region of calorimeter
      double EDepCrit = CalcEDepCal(primID, mother_to_daughter_ID, ID_to_ECalDep, crit_reg, ecal_layer_deposits);
      double EDepTot = CalcEDepCal(primID, mother_to_daughter_ID, ID_to_ECalDep, ECALOuterRadius, ecal_layer_deposits);
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
        double p_x = _MCPStartPX->at(prim_index)/1000.;
        double p_y = _MCPStartPY->at(prim_index)/1000.;
        double p_z = _MCPStartPZ->at(prim_index)/1000.;
        double p2 = p_x*p_x + p_y*p_y + p_z*p_z;
        if (p2 > muon_p2) {
          muon_p2 = p2;
          double mass = MaCh3Utils::GetMassFromPDG(pdg);
          duneobj->rw_lep_pX[i_event] = p_x;
          duneobj->rw_lep_pY[i_event] = p_y;
          duneobj->rw_lep_pZ[i_event] = p_z;
          duneobj->rw_elep_true[i_event] = std::sqrt(p2 + mass*mass);
          vertex = {_MCPStartX->at(prim_index), _MCPStartY->at(prim_index), _MCPStartZ->at(prim_index)};
        }
      }
    }
    if (vertex[0] == M3::_BAD_DOUBLE_) MACH3LOG_INFO("No vertex found in event {}", i_event);
    if (std::abs((duneobj->rw_lep_pX[i_event]-_PXlep)/_PXlep) > 0.0001 || std::abs((duneobj->rw_lep_pY[i_event]-_PYlep)/_PYlep) > 0.0001 || std::abs((duneobj->rw_lep_pZ[i_event]-_PZlep)/_PZlep) > 0.0001) {
      MACH3LOG_INFO("Genie p_lep = [{}, {}, {}] and FastGArSim p_lep = [{}, {}, {}]", _PXlep, _PYlep, _PZlep, duneobj->rw_lep_pX[i_event], duneobj->rw_lep_pY[i_event], duneobj->rw_lep_pZ[i_event]);
    }
    if (std::abs((_Elep-duneobj->rw_elep_true[i_event])/_Elep) > 0.001) MACH3LOG_INFO("Found primary muon of energy {} MeV. Should be {} MeV.", duneobj->rw_elep_true[i_event], _Elep);

    if (isEventAccepted) duneobj->is_accepted[i_event]=1;
    else duneobj->is_accepted[i_event]=0;

    double lep_momentum = std::sqrt(duneobj->rw_lep_pX[i_event]*duneobj->rw_lep_pX[i_event] + duneobj->rw_lep_pY[i_event]*duneobj->rw_lep_pY[i_event] + duneobj->rw_lep_pZ[i_event]*duneobj->rw_lep_pZ[i_event]);
    double lep_pBeam = duneobj->rw_lep_pX[i_event]*BeamDirection[0] + duneobj->rw_lep_pY[i_event]*BeamDirection[1] + duneobj->rw_lep_pZ[i_event]*BeamDirection[2];
    double lep_pB = duneobj->rw_lep_pX[i_event];
    double lep_pPerp = duneobj->rw_lep_pY[i_event]*BeamDirection[2] - duneobj->rw_lep_pZ[i_event]*BeamDirection[1];

    double lep_beamangle = acos(lep_pBeam/lep_momentum)*180/M_PI; //Angle to beam (beam direction: [0.0,0.101,0.995])
    double lep_bangle = acos(lep_pB/lep_momentum)*180/M_PI; //Angle to B-field (b-field along x)
    double lep_perpangle = acos(lep_pPerp/lep_momentum)*180/M_PI; //Angle to axis perpendicular to beam and B
    double lep_phi = atan2(lep_pPerp, lep_pB)*180/M_PI;

    duneobj->rw_lep_theta[i_event] = lep_beamangle;
    duneobj->rw_lep_phi[i_event] = lep_phi;
    duneobj->rw_lep_bangle[i_event] = lep_bangle;
    duneobj->rw_lep_p[i_event] = lep_momentum;

    double radius = std::sqrt((vertex[1]-TPC_centre_y)*(vertex[1]-TPC_centre_y) 
                              + (vertex[2]-TPC_centre_z)*(vertex[2]-TPC_centre_z)); //find radius of interaction vertex
    duneobj->rw_vtx_x[i_event] = vertex[0]-TPC_centre_x;
    duneobj->rw_vtx_y[i_event] = vertex[1]-TPC_centre_y;
    duneobj->rw_vtx_z[i_event] = vertex[2]-TPC_centre_z;

    if(std::abs(vertex[0] - TPC_centre_x) <= TPCFidLength && radius<=TPCFidRadius){
      num_in_fdv++;
      duneobj->in_fdv[i_event] = 1;
    } else{
      duneobj->in_fdv[i_event] = 0;
    }

    duneobj->geometric_correction[i_event] = 1.;
    if (do_geometric_correction) {
      if ((lep_bangle < 45 || lep_bangle > 135) && lep_momentum > 0.3) duneobj->geometric_correction[i_event] = 0.;
      else if ((lep_perpangle < 45 || lep_perpangle > 135) && lep_momentum > 0.3) duneobj->geometric_correction[i_event] = 2.;
    }

    duneobj->rw_rad[i_event] = radius;

    //Assume everything is on Argon for now....
    duneobj->Target[i_event] = 40;

    duneobj->rw_Q0[i_event] = _Enu - _Elep;
    duneobj->rw_Q3[i_event] = std::sqrt((_PXnu-_PXlep)*(_PXnu-_PXlep) + (_PYnu-_PYlep)*(_PYnu-_PYlep) + (_PZnu-_PZlep)*(_PZnu-_PZlep));
    duneobj->rw_lep_pT[i_event] = std::sqrt(lep_momentum*lep_momentum - lep_pBeam*lep_pBeam); 

    duneobj->flux_w[i_event] = 1.0;
    if(duneobj->rw_isCC[i_event] == 1) numCC++;
  }
  MACH3LOG_INFO("nEvents = {}, numCC = {}, numFDV = {}", duneobj->nEvents, numCC, num_in_fdv);

  simFile->Close();
  genieFile->Close();
  ecal_deposits.close();
  return duneobj->nEvents;
}

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wswitch-enum"
const double* SampleHandlerBeamNDGAr::GetPointerToKinematicParameter(KinematicTypes KinematicParameter, int iSample, int iEvent) {
  switch(KinematicParameter) {
    case kTrueNeutrinoEnergy:
      return &dunendgarmcSamples[iSample].rw_etru[iEvent]; 
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
    case kParticle_EndX:
      return &dunendgarmcSamples[iSample].particle_endx->at(iEvent);
    case kParticle_EndY:
      return &dunendgarmcSamples[iSample].particle_endy->at(iEvent);
    case kParticle_EndZ:
      return &dunendgarmcSamples[iSample].particle_endz->at(iEvent);
    case kParticle_ECALDepth:
      return &dunendgarmcSamples[iSample].particle_ecaldepth->at(iEvent);
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

  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent){
    fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
    // fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
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
