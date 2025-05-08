#include "samplePDFDUNEBeamNDGAr.h"

samplePDFDUNEBeamNDGAr::samplePDFDUNEBeamNDGAr(std::string mc_version_, covarianceXsec* XsecCov_) : samplePDFFDBase(mc_version_, XsecCov_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;

  Initialise();
}

samplePDFDUNEBeamNDGAr::~samplePDFDUNEBeamNDGAr() {
}

void samplePDFDUNEBeamNDGAr::Init() {
  dunendgarmcSamples.resize(nSamples,dunemc_base());
  nparticlesinsample = new int[nSamples]();
  pot = SampleManager->raw()["POT"].as<double>();
  iscalo_reco = SampleManager->raw()["SampleBools"]["iscalo_reco"].as<bool>(); //NK determine what reco used
  iselike = SampleManager->raw()["SampleBools"]["iselike"].as<bool>();
  incl_geant = SampleManager->raw()["SampleBools"]["incl_geant"].as<bool>(); //NK whether to read GEANT files
  ecal_containment = SampleManager->raw()["SampleBools"]["ecal_containment"].as<bool>(); //NK do we count containment if its stopped in ECAL
  muonscore_threshold = SampleManager->raw()["SampleCuts"]["muonscore_threshold"].as<double>(); //NK determine what muon score threshold to use
  protondEdxscore = SampleManager->raw()["SampleCuts"]["protondEdxscore_threshold"].as<double>(); //NK determine what proton score threshold to use
  protontofscore = SampleManager->raw()["SampleCuts"]["protontofscore_threshold"].as<double>();  //NK determine what muon score threshold to use
  recovertexradiusthreshold =  SampleManager->raw()["SampleCuts"]["recovertexradius_threshold"].as<double>();  //NK determine what radius threshold to use
  pionenergy_threshold = (SampleManager->raw()["SampleCuts"]["pionenergy_threshold"].as<double>())/1000; //NK determine what muon score threshold to use
  B_field = SampleManager->raw()["SampleCuts"]["B_field"].as<double>(); //NK B field value in T
  momentum_resolution_threshold = SampleManager->raw()["SampleCuts"]["momentum_resolution_threshold"].as<double>(); //NK momentum_resolution threshold, total as a fraction of momentum
  pixel_spacing = SampleManager->raw()["SampleCuts"]["pixel_spacing"].as<double>(); //NK pixel spacing in mm to find num hits in y,z plane
  spatial_resolution = SampleManager->raw()["SampleCuts"]["spatial_resolution"].as<double>(); //NK spatial resolution in mm to find  in y,z plane
  adc_sampling_frequency = SampleManager->raw()["SampleCuts"]["adc_sampling_frequency"].as<double>(); //NK sampling frequency for ADC - needed to find timing resolution and spatial resolution in x dir in MHz
  drift_velocity = SampleManager->raw()["SampleCuts"]["drift_velocity"].as<double>(); //NK drift velocity of electrons in gas - needed to find timing resolution and spatial resolution in x dir in cm/microsecond
  //  average_gain = SampleManager->raw()["SampleCuts"]["average_gain"].as<double>();
  pi0_reco_efficiency = SampleManager->raw()["SampleCuts"]["pi0_reco_efficiency"].as<double>(); //efficiency for pi0 reco in ECAL
  gamma_reco_efficiency = SampleManager->raw()["SampleCuts"]["gamma_reco_efficiency"].as<double>(); //efficiency for gamma reco in ECAL
  TPCFidLength = SampleManager->raw()["SampleCuts"]["TPCFidLength"].as<double>();
  TPCFidRadius = SampleManager->raw()["SampleCuts"]["TPCFidRadius"].as<double>();
  TPCInstrumentedLength = SampleManager->raw()["SampleCuts"]["TPCInstrumentedLength"].as<double>();
  TPCInstrumentedRadius = SampleManager->raw()["SampleCuts"]["TPCInstrumentedRadius"].as<double>();
  ECALInnerRadius = SampleManager->raw()["SampleCuts"]["ECALInnerRadius"].as<double>();
  ECALOuterRadius = SampleManager->raw()["SampleCuts"]["ECALOuterRadius"].as<double>();
  ECALEndCapStart = SampleManager->raw()["SampleCuts"]["ECALEndCapStart"].as<double>();
  ECALEndCapEnd = SampleManager->raw()["SampleCuts"]["ECALEndCapEnd"].as<double>();
}

void samplePDFDUNEBeamNDGAr::SetupSplines() {
  ///@todo move all of the spline setup into core
  if(XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline) > 0){
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    //splinesDUNE* DUNESplines = new splinesDUNE(XsecCov);
    SplineHandler = std::unique_ptr<splineFDBase>(new splinesDUNE(XsecCov,Modes));
    //splineFile = (splineFDBase*)DUNESplines;
    InitialiseSplineObject();
  }
  else{
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or evaluate splines", XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    //splineFile = nullptr;
    SplineHandler = nullptr;
  }

  return;
}

void samplePDFDUNEBeamNDGAr::SetupWeightPointers() {
  for (int i = 0; i < static_cast<int>(dunendgarmcSamples.size()); ++i) {
    for (int j = 0; j < dunendgarmcSamples[i].nEvents; ++j) {
      MCSamples[i].ntotal_weight_pointers[j] = 6;
      //MCSamples[i].total_weight_pointers[j] = new const double*[MCSamples[i].ntotal_weight_pointers[j]];
      MCSamples[i].total_weight_pointers[j].resize(MCSamples[i].ntotal_weight_pointers[j]);
      MCSamples[i].total_weight_pointers[j][0] = &(dunendgarmcSamples[i].pot_s);
      MCSamples[i].total_weight_pointers[j][1] = &(dunendgarmcSamples[i].norm_s);
      MCSamples[i].total_weight_pointers[j][2] = MCSamples[i].osc_w_pointer[j];
      MCSamples[i].total_weight_pointers[j][3] = &(dunendgarmcSamples[i].rw_berpaacvwgt[j]);
      MCSamples[i].total_weight_pointers[j][4] = &(dunendgarmcSamples[i].flux_w[j]);
      MCSamples[i].total_weight_pointers[j][5] = &(MCSamples[i].xsec_w[j]);
    }
  }
}

void samplePDFDUNEBeamNDGAr::makePixelGrid(double pixel_spacing_cm){ //make a square pixel grid with spacing defined in yaml file. Spacing must be input in the yaml file in mm, and then is converted to cm later in the code
  int numpixelrows = static_cast<int>(floor(TPCInstrumentedRadius*2/(pixel_spacing_cm))); //find number of pixels along y and z axis.
  double centre_yboundary, centre_zboundary;
  if(numpixelrows % 2 == 0){centre_yboundary=TPC_centre_y; centre_zboundary=TPC_centre_z;}
  else{centre_yboundary=TPC_centre_y-(pixel_spacing_cm/2); centre_zboundary=TPC_centre_z-(pixel_spacing_cm/2);}
  pixelymin = centre_yboundary - floor(numpixelrows/2)*pixel_spacing_cm;
  pixelymax = centre_yboundary + floor(numpixelrows/2)*pixel_spacing_cm;
  pixelzmin = centre_zboundary - floor(numpixelrows/2)*pixel_spacing_cm;
  pixelzmax = centre_zboundary + floor(numpixelrows/2)*pixel_spacing_cm;

  for(int i_pixel = 0; i_pixel<numpixelrows; i_pixel++){
    yboundarypositions.push_back(pixelymin+i_pixel*pixel_spacing_cm);
    zboundarypositions.push_back(pixelzmin+i_pixel*pixel_spacing_cm);
  }
}

double samplePDFDUNEBeamNDGAr::FindNHits(bool positivecharged, double pixel_spacing_cm, double centre_circle_y, double centre_circle_z, double rad_curvature, double theta_start, double theta_end){
  //use the pixel grid method to find number of pixels hit in a track.
  int num_vertices =0; // number of vertices hit. Counting to avoid duplicating
  int num_intersections =0;

  //equation for circle = (y-y0)^2 + (z-z0)^2 = r^2

  for(unsigned int i_intersect = 0; i_intersect<yboundarypositions.size(); i_intersect++){ //check every boundary line to see if it crossed within the TPC instrumented region
    double quadratic_ineq_y= (rad_curvature*100)*(rad_curvature*100)-(yboundarypositions[i_intersect]-centre_circle_y)*(yboundarypositions[i_intersect]-centre_circle_y);
    if (quadratic_ineq_y < 0) continue; //circle does not intersect this boundary
    
    double z_coord1 = centre_circle_z + std::sqrt(quadratic_ineq_y);
    double y_coord = yboundarypositions[i_intersect];
    if(isCoordOnTrack(positivecharged, y_coord, z_coord1, centre_circle_y, centre_circle_z, theta_start, theta_end)){ // check point is in tpc volume and on arc of track
      num_intersections++;
      if(fmod((z_coord1 - pixelzmin), pixel_spacing_cm) == 0){ // this is the case when a vertex is crossed so to avoid double counting pixels
        num_vertices++;
      }
    }

    if (quadratic_ineq_y == 0) continue; //pixel boundary tangent so one z solution

    double z_coord2 = centre_circle_z - std::sqrt(quadratic_ineq_y);
    if(isCoordOnTrack(positivecharged, y_coord, z_coord2, centre_circle_y, centre_circle_z, theta_start, theta_end)){ // check point is in tpc volume and on arc of track
      num_intersections++;
      if(fmod((z_coord2 - pixelzmin), pixel_spacing_cm) == 0){ // this is the case when a vertex is crossed so to avoid double counting pixels
        num_vertices++;
      }
    }
  }
  for(unsigned int i_intersect = 0; i_intersect<zboundarypositions.size(); i_intersect++){
    double quadratic_ineq_z = (rad_curvature*100)*(rad_curvature*100)-(zboundarypositions[i_intersect]-centre_circle_z)*(zboundarypositions[i_intersect]-centre_circle_z);
    if(quadratic_ineq_z < 0) continue;
    
    double y_coord1 = centre_circle_y + std::sqrt(quadratic_ineq_z);
    double z_coord = zboundarypositions[i_intersect];
    if(isCoordOnTrack(positivecharged, y_coord1, z_coord, centre_circle_y, centre_circle_z, theta_start, theta_end)){
      num_intersections++;
    }

    if (quadratic_ineq_z == 0) continue; //pixel boundary tangent so one y solution

    double y_coord2 = centre_circle_y - std::sqrt(quadratic_ineq_z)-TPC_centre_y;
    if(isCoordOnTrack(positivecharged, y_coord2, z_coord, centre_circle_y, centre_circle_z, theta_start, theta_end)){ 
      num_intersections++;
    }
    // already checked all vertices for duplicates before so no need to repeat that
  }
  double N_hits = num_intersections - num_vertices + 1; //Add one for the pixel that it starts on
  return N_hits;
}

bool samplePDFDUNEBeamNDGAr::isCoordOnTrack(bool positivecharged, double ycoord, double zcoord, double centre_circle_y, double centre_circle_z, double theta_start, double theta_end) {
  double theta_coord = atan2(ycoord - centre_circle_y, zcoord - centre_circle_z);
  if (theta_coord<0) theta_coord += 2*M_PI;

  bool iscoordinTPC = (ycoord - TPC_centre_y)*(ycoord - TPC_centre_y)+(zcoord - TPC_centre_z)*(zcoord - TPC_centre_z) < TPCInstrumentedRadius*TPCInstrumentedRadius;
  
  bool fullrotation = false;
  if (theta_end > 2*M_PI) {
    if (theta_end < 2*M_PI + theta_start) theta_end -= 2*M_PI;
    else fullrotation = true;
  }

  bool iscoordinarc;
  if (!positivecharged) {
    if (theta_start < theta_end) iscoordinarc = theta_coord < theta_start || theta_coord > theta_end;
    else iscoordinarc = theta_coord > theta_end && theta_coord < theta_start;
  }
  else {
    if (theta_start < theta_end) iscoordinarc = theta_coord > theta_start && theta_coord < theta_end;
    else iscoordinarc = theta_coord > theta_start || theta_coord < theta_end;
  }

  return iscoordinTPC && (iscoordinarc || fullrotation);
}
/*
double samplePDFDUNEBeamNDGAr::CalcMomResYZ(double pixel_spacing_cm, double centre_circle_y, double centre_circle_z, double rad_curvature, double theta_start, double theta_end, bool positivecharged) {

  //summations for resolution calculation
  double sumy=0, sumz=0, sumr2=0, sumyz=0, sumyr2=0, sumzr2=0;
  int nhits = 0;

  //define yz pixel grid which covers entire tpc cross section
  int numpixelrows = static_cast<int>(floor(TPCInstrumentedRadius*2/(pixel_spacing_cm)))+1;
  if (numpixelrows%2 == 0) numpixelrows += 1;
  int pixelmin = -(numpixelrows-1)/2;
  int pixelmax = (numpixelrows+1)/2;

  //Loop through all pixels in grid
  for (int ypixel=pixelmin; ypixel<=pixelmax; ypixel++) {
    for (int zpixel=pixelmin; zpixel<=pixelmax; zpixel++) {
      double pixel_centre_y = TPC_centre_y + ypixel*pixel_spacing_cm;
      double pixel_centre_z = TPC_centre_z + zpixel*pixel_spacing_cm;

      //Check if track passes through pixel: nearest corner is contained in the circle the furthest corner is not
      double near_corner_y, near_corner_z;
      double far_corner_y, far_corner_z;
      if (centre_circle_y > pixel_centre_y) near_corner_y = pixel_centre_y + 0.5*pixel_spacing_cm;
      else near_corner_y = pixel_centre_y - 0.5*pixel_spacing_cm;
      if (centre_circle_z > pixel_centre_z) near_corner_z = pixel_centre_z + 0.5*pixel_spacing_cm;
      else near_corner_z = pixel_centre_z - 0.5*pixel_spacing_cm;
      
      double near_corner_r2 = ((near_corner_y - centre_circle_y)*(near_corner_y - centre_circle_y) 
          + (near_corner_z - centre_circle_z)*(near_corner_z - centre_circle_z));
      double far_corner_r2 = ((far_corner_y - centre_circle_y)*(far_corner_y - centre_circle_y) 
          + (far_corner_z - centre_circle_z)*(far_corner_z - centre_circle_z));

      if (near_corner_r2 < rad_curvature*rad_curvature && far_corner_r2 > rad_curvature*rad_curvature) {
        if (isCoordOnTrack(positivecharged, pixel_centre_y, pixel_centre_z, centre_circle_y, centre_circle_z, theta_start, theta_end)) {
          double r2 = pixel_centre_y*pixel_centre_y + pixel_centre_z*pixel_centre_z;
          sumy += pixel_centre_y;
          sumz += pixel_centre_z;
          sumr2 += r2;
          sumyy += pixel_centre_y*pixel_centre_y;
          sumzz += pixel_centre_z*pixel_centre_z;
          sumr2r2 += r2*r2;
          sumyz += pixel_centre_y*pixel_centre_z;
          sumyr2 += pixel_centre_y*r2;
          sumzr2 += pixel_centre_z*r2;
          nhits++;
        }
      }
    }
  }
  //Calculation of momres follows V. Karimaki CMS note 
  //'Explicit Covariance Matrix for Particle Measurement Precision'
  double sigyy = sumyy/nhits - sumy*sumy/(nhits*nhits);
  double sigzz = sumzz/nhits - sumz*sumz/(nhits*nhits);
  double sigr2r2 = sumr2r2/nhits - sumr2*sumr2/(nhits*nhits);
  double sigyz = sumyz/nhits - sumy*sumz/(nhits*nhits);
  double sigyr2 = sumyr2/nhits - sumy*sumr2/(nhits*nhits);
  double sigzr2 = sumzr2/nhits - sumz*sumr2/(nhits*nhits;

  double phi = 0.5*atan(2*(sigr2r2*sigyz - sigyr2*sigzr2)/(sigr2r2*(sigyy-sigzz) - sigyr2*sigyr2 + sigzr2*sigzr2));
  double rho = 2*(sin(phi)*sigyr2 - cos(phi)*sigzr2)/sigr2r2;
  double rad_curv_est = 1/rho;
  double momres_yz = 0.3*B_field*(std::abs(rad_curv_est) - rad_curvature);
  return momres_yz;
}
*/
double samplePDFDUNEBeamNDGAr::CalcBeta(double p_mag, double& bg, double& gamma, double pdgmass){ //calculate beta (v/c)
  bg = p_mag/pdgmass; //beta*gamma
  gamma = std::sqrt(1+bg*bg); //gamma
  double beta = bg/gamma; //beta (velocity)
  return beta;
}

bool samplePDFDUNEBeamNDGAr::IsParticleAccepted(dunemc_base *duneobj, int i_sample, int i_event, int i_truepart, double pixel_spacing_cm, bool *isgoodcafparticle, double pdgmass){
  static std::random_device rd;
  static std::mt19937 gen(rd());
  static std::uniform_real_distribution<double> dist(0.0, 1.0);

  bool isAccepted = true;
  int nummatched = 0;
  int pdgcaf = sr->mc.nu[0].prim[i_truepart].pdg;

  //Loop through particles in anatree
  for(unsigned int i_anapart =0; i_anapart<_MCPStartPX->size(); i_anapart++){

    //Select particle in the anatree that is the same as in the caf (same pdg and momentum)
    if((_PDG->at(i_anapart) != pdgcaf) || 
        (std::abs(_MCPStartPX->at(i_anapart)-static_cast<double>(sr->mc.nu[0].prim[i_truepart].p.px))>0.001*std::abs(_MCPStartPX->at(i_anapart))) ||
        (std::abs(_MCPStartPY->at(i_anapart)-static_cast<double>(sr->mc.nu[0].prim[i_truepart].p.py))>0.001*std::abs(_MCPStartPY->at(i_anapart))) ||
        (std::abs(_MCPStartPZ->at(i_anapart)-static_cast<double>(sr->mc.nu[0].prim[i_truepart].p.pz))>0.001*std::abs(_MCPStartPZ->at(i_anapart)))){
      continue;
    }
    nummatched++;

    auto it = std::find(_MCPTrkID->begin(), _MCPTrkID->end(), sr->mc.nu[0].prim[i_truepart].G4ID+1);
    if (i_anapart == it - _MCPTrkID->begin()) *isgoodcafparticle = true;

    //Ignore neutrons and neutrinos (they will be accepted by default for now but will not appear on particle-level plots)
    if(std::abs(pdgcaf) == 2112 || std::abs(pdgcaf) == 14 || std::abs(pdgcaf) == 12 
        /*||(std::abs(_MCPStartX->at(i_anapart))-TPC_centre_x)>=TPCFidLength || start_radius>=TPCFidRadius*/){
      continue;
    }

    //Determine if stopped in ecal length/radius
    double ystart = _MCPStartY->at(i_anapart)-TPC_centre_y;
    double zstart = _MCPStartZ->at(i_anapart)-TPC_centre_z;
    double start_radius = std::sqrt(ystart*ystart + zstart*zstart);

    double yend = _MCPEndY->at(i_anapart)-TPC_centre_y;
    double zend = _MCPEndZ->at(i_anapart)-TPC_centre_z;
    double end_radius = std::sqrt(yend*yend + zend*zend); 
    double end_length = _MCPEndX->at(i_anapart)-TPC_centre_x;

    bool stops_in_tpc = std::abs(end_length)<=TPCInstrumentedLength && end_radius<=TPCInstrumentedRadius;
    bool stops_before_ecal = std::abs(end_length)<ECALEndCapStart && end_radius<ECALInnerRadius;
    bool stops_beyond_ecal = std::abs(end_length)>ECALEndCapEnd || end_radius>ECALOuterRadius;
    bool stops_in_ecal = !(stops_before_ecal || stops_beyond_ecal);

    double transverse_mom = std::sqrt(_MCPStartPY->at(i_anapart)*_MCPStartPY->at(i_anapart)
        + _MCPStartPZ->at(i_anapart)*_MCPStartPZ->at(i_anapart));
    double mom_tot = std::sqrt(_MCPStartPX->at(i_anapart)*_MCPStartPX->at(i_anapart) + transverse_mom*transverse_mom);

    //Fill particle-level kinematic variables with default or actual (if possible at this stage) values
    duneobj->particle_event->push_back(i_event);
    duneobj->particle_pdg->push_back(pdgcaf);
    duneobj->particle_energy->push_back(static_cast<double>(sr->mc.nu[0].prim[i_truepart].p.E));
    duneobj->particle_momentum->push_back(mom_tot);
    duneobj->particle_transversemomentum->push_back(transverse_mom);
    duneobj->particle_bangle->push_back(acos(_MCPStartPX->at(i_anapart)/mom_tot)*180/M_PI);//Angle to B-field
    duneobj->particle_startx->push_back(_MCPStartX->at(i_anapart)-TPC_centre_x);
    duneobj->particle_startr2->push_back(start_radius*start_radius);
    duneobj->particle_isaccepted->push_back(true); //default (updates later)
    duneobj->particle_isstoppedingap->push_back(!stops_in_tpc && stops_before_ecal);
    duneobj->particle_isstoppedinbarrelgap->push_back(!stops_in_tpc && stops_before_ecal && std::abs(end_length)<=TPCInstrumentedLength); 
    duneobj->particle_isstoppedinendgap->push_back(!stops_in_tpc && stops_before_ecal && std::abs(end_length)>TPCInstrumentedLength); 
    duneobj->particle_isstoppedinecal->push_back(stops_in_ecal);
    duneobj->particle_isstoppedinbarrel->push_back(stops_in_ecal && end_radius>=ECALInnerRadius);
    duneobj->particle_isstoppedinendcap->push_back(stops_in_ecal && end_radius<ECALInnerRadius);
    duneobj->particle_isstoppedintpc->push_back(stops_in_tpc);
    duneobj->particle_momresms->push_back(std::numeric_limits<double>::quiet_NaN()); //default (updates later)
    duneobj->particle_momrestransfrac->push_back(std::numeric_limits<double>::quiet_NaN()); //default (updates later)
    duneobj->particle_momrestrans->push_back(std::numeric_limits<double>::quiet_NaN()); //default (updates later)
    duneobj->particle_ecaldepositfraction->push_back(std::numeric_limits<double>::quiet_NaN()); //default (updates later)
    duneobj->particle_nhits->push_back(std::numeric_limits<double>::quiet_NaN()); //default (updates later)
    duneobj->particle_nturns->push_back(std::numeric_limits<double>::quiet_NaN()); //default (updates later)
    nparticlesinsample[i_sample]++;

    //If particle is not stopped in the tpc or ecal 
    if(!stops_in_tpc && !stops_in_ecal){
      //Check if charged (p +/- , pi +/- , mu +/- , e +/- , K +/-)
      //JM why not more (eg. sig +/-)
      if(std::abs(pdgcaf) == 2212 || std::abs(pdgcaf) == 211 || std::abs(pdgcaf) == 13 || std::abs(pdgcaf) == 11 || std::abs(pdgcaf) == 321) {

        double length_track_x;
        if(std::abs(end_length)>TPCInstrumentedLength){
          if((end_length)>=0){length_track_x = TPCInstrumentedLength - (_MCPStartX->at(i_anapart)-TPC_centre_x);} //in cm
          else{length_track_x = -TPCInstrumentedLength - (_MCPStartX->at(i_anapart)-TPC_centre_x);} //in cm
        }
        else{length_track_x = _MCPEndX->at(i_anapart) - _MCPStartX->at(i_anapart);} //in cm

        double rad_curvature = transverse_mom/(0.3*B_field); //p = 0.3*B*r where p in GeV/c, B in T, r in m
        double theta_xT = atan(_MCPStartPX->at(i_anapart)/transverse_mom); //helix pitch angle
        double pitch = std::abs(2*M_PI*rad_curvature*tan(theta_xT)); //distance between two turns of a helix in m
        double tan_theta = tan(theta_xT);

        //Find centre of circular path
        bool positivecharged = 0;
        double centre_circle_y;
        double centre_circle_z;
        double L_yz, L_yz_chord; //length of curved track in y-z plane
        if(pdgcaf == 2212 || pdgcaf == 211 || pdgcaf == -13 || pdgcaf == -11 || pdgcaf == 321){positivecharged = 1;}
        if(positivecharged){
          centre_circle_y = _MCPStartY->at(i_anapart) + (rad_curvature*100*_MCPStartPZ->at(i_anapart)/transverse_mom); //Note plus sign here as cross product gives F in direction of ( pz j - py k) F= q v x B
          centre_circle_z = _MCPStartZ->at(i_anapart) - (rad_curvature*100*_MCPStartPY->at(i_anapart)/transverse_mom);
        }
        else if(!positivecharged){
          centre_circle_y = _MCPStartY->at(i_anapart) - (rad_curvature*100*_MCPStartPZ->at(i_anapart)/transverse_mom); //Note minus sign here as cross product gives F in direction of ( -pz j + py k)
          centre_circle_z = _MCPStartZ->at(i_anapart) + (rad_curvature*100*_MCPStartPY->at(i_anapart)/transverse_mom);
        }
        //Find Position where track leaves TPC. Intersection of two circles.
        //JM do you consider the case where the particle leaves out the end cap?
        double m_const = (TPC_centre_z - centre_circle_z)/(TPC_centre_y-centre_circle_y); //gradient of line between two intersection points
        double a_const = (TPCInstrumentedRadius*TPCInstrumentedRadius-(rad_curvature*100)*(rad_curvature*100) - (TPC_centre_y*TPC_centre_y - centre_circle_y*centre_circle_y)-(TPC_centre_z*TPC_centre_z - centre_circle_z*centre_circle_z))/(2*(centre_circle_y-TPC_centre_y));
        double quadraticformula_b = -(2*m_const*(a_const - TPC_centre_y) + 2*TPC_centre_z);
        double quadraticformula_a = m_const*m_const + 1;
        double quadraticformula_c = (a_const - TPC_centre_y)*(a_const - TPC_centre_y) + TPC_centre_z*TPC_centre_z - TPCInstrumentedRadius*TPCInstrumentedRadius;

        double z_intersect_1, y_intersect_1, z_intersect_2, y_intersect_2, theta_1, theta_2, theta_chosen, theta_diff_1, theta_diff_2, theta_end;
        double nturns = 0;

        double theta_start = atan2(_MCPStartY->at(i_anapart) - centre_circle_y, _MCPStartZ->at(i_anapart) - centre_circle_z);

        if(quadraticformula_b*quadraticformula_b - 4*quadraticformula_a*quadraticformula_c > 0){
          z_intersect_1 = (-quadraticformula_b + std::sqrt(quadraticformula_b*quadraticformula_b - 4*quadraticformula_a*quadraticformula_c))/(2*quadraticformula_a);
          y_intersect_1 = -m_const*z_intersect_1 + a_const;
          z_intersect_2 = (-quadraticformula_b - std::sqrt(quadraticformula_b*quadraticformula_b - 4*quadraticformula_a*quadraticformula_c))/(2*quadraticformula_a);
          y_intersect_2 = -m_const*z_intersect_2 + a_const;

          //Find angle wrt y in yz plane where track starts and where it intersects TPC boundary
          theta_1 = atan2(y_intersect_1 - centre_circle_y, z_intersect_1 - centre_circle_z);
          theta_2 = atan2(y_intersect_2 - centre_circle_y, z_intersect_2 - centre_circle_z);

          //Ensure they are in the range [0,2pi)
          if (theta_1 < 0) theta_1 += 2*M_PI;
          if (theta_2 < 0) theta_2 += 2*M_PI;
          if (theta_start < 0) theta_start += 2*M_PI;

          //Lorentz force law, if positively charged, counter clockwise, if negative charged, clockwise (looking down the B-field)
          if(!positivecharged){
            theta_diff_1 = (theta_1 < theta_start) ? (theta_start - theta_1) : (2*M_PI - (theta_1 - theta_start));
            theta_diff_2 = (theta_2 < theta_start) ? (theta_start - theta_2) : (2*M_PI - (theta_2 - theta_start));
          }
          else if(positivecharged){ 
            theta_diff_1 = (theta_1 > theta_start) ? (theta_1 - theta_start) : (2*M_PI - (theta_start - theta_1));
            theta_diff_2 = (theta_2 > theta_start) ? (theta_2 - theta_start) : (2*M_PI - (theta_start - theta_2));
          }

          //JM Why is the second check necessary?
          if(theta_diff_1<theta_diff_2 && (rad_curvature*100*theta_diff_1 > (TPCInstrumentedRadius-start_radius))){
            theta_chosen = theta_diff_1;
            theta_end = theta_1;
          }
          else {
            theta_chosen = theta_diff_2;
            theta_end = theta_2;
          }

          L_yz = rad_curvature*100*theta_chosen;
          L_yz_chord = std::abs(2*rad_curvature*100*sin(theta_chosen/2));
          nturns = std::abs(length_track_x/100)/pitch;

          //If track ends before it would intersect, adjust
          if(std::abs(L_yz*tan_theta) > std::abs(length_track_x)){
            theta_chosen = fmod(nturns, 1)*2*M_PI; //JM why do we need the fmod? nturns<1 in this case surely?
            L_yz = rad_curvature*100*theta_chosen;
            L_yz_chord = std::abs(2*rad_curvature*100*sin(theta_chosen/2));
          } 
          else {
            length_track_x = std::abs(L_yz*tan_theta);
          }
        }
        else { //Radius of curvature small enough to spiral within TPC
          nturns = std::abs(length_track_x/100)/pitch;
          double theta_intersect = nturns*2*M_PI;
          L_yz = theta_intersect*rad_curvature*100;
          theta_end = theta_start + theta_intersect;
          L_yz_chord = L_yz;
        }
        double N_hits = FindNHits(positivecharged, pixel_spacing_cm, centre_circle_y, centre_circle_z, rad_curvature, theta_start, theta_end);
        double p_mag = sr->mc.nu[0].prim[i_truepart].p.Mag();
        double bg = 0; 
        double gamma = 0;
        double beta = CalcBeta(p_mag, bg, gamma, pdgmass);

        //JM removed avg_betapT calc (seems to be unused)

        double sigmax = (drift_velocity/100)*(1/(adc_sampling_frequency));
        double sigmax_frac = sigmax/(std::abs(length_track_x)/100);
        double sigmayz = (spatial_resolution/(1000)); //needs to be in m              
        double momres_yz = transverse_mom*(std::sqrt(720/(N_hits+4)) * (sigmayz*transverse_mom/(0.3*B_field*(L_yz_chord/100)*(L_yz_chord/100))) 
            * std::sqrt(1-(1/21)*(L_yz_chord/(rad_curvature*100))*(L_yz_chord/(rad_curvature*100))));
        double momres_ms = transverse_mom*(0.016/(0.3*B_field*(L_yz/100)*cos(theta_xT)*beta))*std::sqrt(L_yz/X0);
        double momres_tottransverse = std::sqrt(momres_yz*momres_yz + momres_ms*momres_ms)/transverse_mom;
        double sigma_theta = (cos(theta_xT)*cos(theta_xT) * (pitch/(2*M_PI*rad_curvature)) *
            std::sqrt(sigmax_frac*sigmax_frac + momres_tottransverse*momres_tottransverse));
        double momres_frac = std::sqrt(momres_tottransverse*momres_tottransverse + (sigma_theta*tan_theta)*(sigma_theta*tan_theta));

        duneobj->particle_nhits->back() = N_hits;
        duneobj->particle_nturns->back() = nturns;
        duneobj->particle_momresms->back() = momres_ms;
        duneobj->particle_momrestransfrac->back() = momres_tottransverse/momres_frac;
        duneobj->particle_momrestrans->back() = momres_tottransverse;

        if(momres_frac > momentum_resolution_threshold){
          isAccepted = false;
        }
      }
      else{ //Neutral particles: accept pi0 and gamma with some efficiency
        if((std::abs(pdgcaf)) == 111){
          double pi0recoprob = dist(gen);
          if(pi0recoprob>pi0_reco_efficiency){
            isAccepted = false;
          }
        }
        else if((std::abs(pdgcaf)) == 22){
          double gammarecoprob = dist(gen);
          if(gammarecoprob>gamma_reco_efficiency){
            isAccepted = false;
          }
        }
        else{
          isAccepted = false;
        }
      }
    }
    else if (stops_in_ecal){
      double energydepsum = 0;
      for(unsigned int i_ecaldep=0; i_ecaldep<_SimHitEnergy->size(); i_ecaldep++){
        if(_SimHitTrkID->at(i_ecaldep) == _MCPTrkID->at(i_anapart)){
          energydepsum = energydepsum + _SimHitEnergy->at(i_ecaldep);
        }
      }
      duneobj->particle_ecaldepositfraction->back() = energydepsum/(sr->mc.nu[0].prim[i_truepart].p.E);
    }
    if (isAccepted == false) {duneobj->particle_isaccepted->back() = false;}
    break;
  }
  if(nummatched != 1){
    MACH3LOG_INFO("Found {} matching particles in anatree", nummatched);
    MACH3LOG_INFO("PDG: {}, momentum: ({},{},{})",pdgcaf,sr->mc.nu[0].prim[i_truepart].p.px,sr->mc.nu[0].prim[i_truepart].p.py,sr->mc.nu[0].prim[i_truepart].p.pz);
  }
  return isAccepted;
}

int samplePDFDUNEBeamNDGAr::setupExperimentMC(int iSample) {

  dunemc_base *duneobj = &(dunendgarmcSamples[iSample]);
  //int nutype = sample_nutype[iSample];
  //int oscnutype = sample_oscnutype[iSample];
  //bool signal = sample_signal[iSample];

  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("Input File: {}", mc_files.at(iSample));

  //Read caf file
  _sampleFile = TFile::Open(mc_files.at(iSample).c_str(), "READ");
  _data = static_cast<TTree*>(_sampleFile->Get("cafTree"));

  if(_data){
    MACH3LOG_INFO("Found \"caf\" tree in {}", mc_files[iSample]);
    MACH3LOG_INFO("With number of entries: {}", _data->GetEntries());
  }
  else{
    MACH3LOG_ERROR("Could not find \"caf\" tree in {}", mc_files[iSample]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  //Read geant file
  if(incl_geant){
    std::string geantfile(mc_files.at(iSample));

    std::string prefix = "inputs/DUNE_NDGAr_CAF_files/";
    std::string::size_type i_prefix = geantfile.find(prefix);
    if(i_prefix != std::string::npos){
      geantfile.erase(i_prefix, prefix.length());
    }

    std::string suffix = ".root";
    std::string::size_type i_suffix = geantfile.find(suffix);
    if(i_suffix != std::string::npos){
      geantfile.erase(i_suffix, suffix.length());
    }

    std::string geantfilename = "inputs/DUNE_NDGAr_AnaTrees/" + geantfile + "_geant.root"; 
    _sampleFile_geant = TFile::Open(geantfilename.c_str(),"READ");
    _data_geant = static_cast<TTree*>(_sampleFile_geant->Get("GArAnaTree"));

    if(_data_geant){
      MACH3LOG_INFO("Found \"anatree\" tree in {}", geantfilename);
      MACH3LOG_INFO("With number of entries: {}", _data_geant->GetEntries());
    }
    else{
      MACH3LOG_ERROR("Could not find \"anatree\" tree in {}", geantfilename);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    _data->AddFriend(_data_geant);
  }

  _data->SetBranchStatus("*", 1);
  _data->SetBranchAddress("rec", &sr);

  if(incl_geant){
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
    _data->SetBranchAddress("SimHitTrkID", &_SimHitTrkID);
    _data->SetBranchAddress("SimHitEnergy", &_SimHitEnergy);
  }
  duneobj->norm_s = 1.;
  double downsampling = 1; //default 1, set to eg. 0.01 for quick testing
  duneobj->pot_s = pot/(downsampling*1e21);
  duneobj->nEvents = static_cast<int>(std::round(downsampling*static_cast<double>(_data->GetEntries())));

  // allocate memory for dunendgarmc variables
  duneobj->rw_yrec = new double[duneobj->nEvents];
  duneobj->rw_elep_reco = new double[duneobj->nEvents];
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_erec = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_nuPDGunosc = new int[duneobj->nEvents];
  duneobj->rw_nuPDG = new int[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 

  duneobj->rw_Q0 = new double[duneobj->nEvents];
  duneobj->rw_Q3 = new double[duneobj->nEvents];
  duneobj->mode = new double[duneobj->nEvents];

  duneobj->nproton = new int[duneobj->nEvents];
  duneobj->nneutron = new int[duneobj->nEvents];
  duneobj->npip = new int[duneobj->nEvents];
  duneobj->npim = new int[duneobj->nEvents];
  duneobj->npi0 = new int[duneobj->nEvents];

  duneobj->nrecomuon = new int[duneobj->nEvents];
  duneobj->ntruemuon = new int[duneobj->nEvents];
  duneobj->nmuonsratio = new double[duneobj->nEvents];
  duneobj->ntruemuonprim = new int[duneobj->nEvents];

  duneobj->nrecoparticles = new int[duneobj->nEvents];
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

  duneobj->rw_reco_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_reco_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_reco_vtx_z = new double[duneobj->nEvents];
  duneobj->rw_reco_rad = new double[duneobj->nEvents];

  duneobj->Target = new int[duneobj->nEvents];

  duneobj->is_accepted = new bool[duneobj->nEvents];
  duneobj->is_good_caf_event = new bool[duneobj->nEvents];

  //Particle-level kinematic variables
  duneobj->particle_event = new std::vector<int>; duneobj->particle_event->reserve(7*duneobj->nEvents);
  duneobj->particle_ecaldepositfraction = new std::vector<double>; duneobj->particle_ecaldepositfraction->reserve(7*duneobj->nEvents);
  duneobj->particle_pdg = new std::vector<int>; duneobj->particle_pdg->reserve(7*duneobj->nEvents); 
  duneobj->particle_energy = new std::vector<double>; duneobj->particle_energy->reserve(7*duneobj->nEvents); 
  duneobj->particle_momentum = new std::vector<double>; duneobj->particle_momentum->reserve(7*duneobj->nEvents); 
  duneobj->particle_transversemomentum = new std::vector<double>; duneobj->particle_transversemomentum->reserve(7*duneobj->nEvents); 
  duneobj->particle_bangle = new std::vector<double>; duneobj->particle_bangle->reserve(7*duneobj->nEvents);
  duneobj->particle_isaccepted = new std::vector<bool>; duneobj->particle_isaccepted->reserve(7*duneobj->nEvents);
  duneobj->particle_isstoppedintpc = new std::vector<bool>; duneobj->particle_isstoppedintpc->reserve(7*duneobj->nEvents); 
  duneobj->particle_isstoppedinecal = new std::vector<bool>; duneobj->particle_isstoppedinecal->reserve(7*duneobj->nEvents); 
  duneobj->particle_isstoppedingap = new std::vector<bool>; duneobj->particle_isstoppedingap->reserve(7*duneobj->nEvents); 
  duneobj->particle_isstoppedinbarrelgap = new std::vector<bool>; duneobj->particle_isstoppedinbarrelgap->reserve(7*duneobj->nEvents); 
  duneobj->particle_isstoppedinendgap = new std::vector<bool>; duneobj->particle_isstoppedinendgap->reserve(7*duneobj->nEvents); 
  duneobj->particle_isstoppedinbarrel= new std::vector<bool>; duneobj->particle_isstoppedinbarrel->reserve(7*duneobj->nEvents);
  duneobj->particle_isstoppedinendcap = new std::vector<bool>; duneobj->particle_isstoppedinendcap->reserve(7*duneobj->nEvents); 
  duneobj->particle_startx = new std::vector<double>; duneobj->particle_startx->reserve(7*duneobj->nEvents); 
  duneobj->particle_startr2 = new std::vector<double>; duneobj->particle_startr2->reserve(7*duneobj->nEvents); 
  duneobj->particle_nhits = new std::vector<double>; duneobj->particle_nhits->reserve(7*duneobj->nEvents); 
  duneobj->particle_nturns = new std::vector<double>; duneobj->particle_nturns->reserve(7*duneobj->nEvents); 
  duneobj->particle_momresms = new std::vector<double>; duneobj->particle_momresms->reserve(7*duneobj->nEvents); 
  duneobj->particle_momrestransfrac = new std::vector<double>; duneobj->particle_momrestransfrac->reserve(7*duneobj->nEvents); 
  duneobj->particle_momrestrans= new std::vector<double>; duneobj->particle_momrestrans->reserve(7*duneobj->nEvents);

  int numCC = 0;
  int num_no_ixns = 0;
  int num_no_recparticles = 0;
  int num_in_fdv = 0;
  int num_in_fdv_noreco = 0;
  int num_notin_fdv = 0;
  int num_nanenergy = 0;
  int num_nanparticles = 0;
  int ngoodcafevents = 0;
  int nbadcafevents = 0;

  double pixel_spacing_cm = pixel_spacing/10; //convert to cm
  makePixelGrid(pixel_spacing_cm); //make square pixel grid and fill vectors with y and z pixel boundaries	

  for (int i_event = 1; i_event < (duneobj->nEvents); ++i_event) { 
    _data->GetEntry(i_event);

    if ((i_event % (duneobj->nEvents/100))==0) {
      MACH3LOG_INFO("\tProcessing event: {}/{}",i_event,duneobj->nEvents);
    }
    double radius = std::sqrt((sr->mc.nu[0].vtx.y-TPC_centre_y)*(sr->mc.nu[0].vtx.y-TPC_centre_y) + (sr->mc.nu[0].vtx.z-TPC_centre_z)*(sr->mc.nu[0].vtx.z-TPC_centre_z)); //find radius of interaction vertex
    if(std::abs(sr->mc.nu[0].vtx.x)<=TPCFidLength &&  radius<=TPCFidRadius){
      num_in_fdv++;
      duneobj->in_fdv[i_event] = 1;
    } else{
      num_notin_fdv++;
      duneobj->in_fdv[i_event] = 0;
    }

    if(sr->common.ixn.ngsft == 0){ //if no reconstructed interaction, fill reco variables with 0
      //duneobj->rw_erec[i_event] = static_cast<double>(0);
      duneobj->rw_elep_reco[i_event] = double(0);
      duneobj->rw_yrec[i_event] = static_cast<double>(0);
      num_no_ixns++;
      duneobj->nrecoparticles[i_event] = 0;
    } else{
      duneobj->nrecoparticles[i_event] = 0;
      double erec_total =0;
      double elep_reco =0;
      double muonscore = muonscore_threshold;
      int nixns = static_cast<int>(sr->common.ixn.ngsft);
      for(int i_ixn =0; i_ixn<nixns; i_ixn++) {
        int nrecoparticles = sr->common.ixn.gsft[i_ixn].part.ngsft;
        duneobj->nrecoparticles[i_event] += sr->common.ixn.gsft[i_ixn].part.ngsft;
        int nanparticles = 0;
        if(nrecoparticles ==0){
          if(duneobj->in_fdv[i_event]) {
            num_in_fdv_noreco++;
          }   
          num_no_recparticles++;
        }

        for(int i_part =0; i_part<nrecoparticles; i_part++) {
          double erec_part = static_cast<double>(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].E);
          if(std::isnan(erec_part)){nanparticles++;}
          erec_total+=erec_part;
          if(static_cast<double>(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].score.gsft_pid.muon_score>muonscore)){
            if(erec_part>elep_reco){
              elep_reco = erec_part;
              duneobj->rw_reco_vtx_x[i_event] = static_cast<double>(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.x);
              duneobj->rw_reco_vtx_y[i_event] = static_cast<double>(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.y);
              duneobj->rw_reco_vtx_z[i_event] = static_cast<double>(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.z);
            }
            duneobj->nrecomuon[i_event]++; 
          }
        }
        num_nanparticles = num_nanparticles + (nanparticles/nrecoparticles);
      } //ADD PRIMARY LEPTON ENERGY ELEP_RECO
      if(std::isnan(erec_total)){
        num_nanenergy++; 
        erec_total = static_cast<double>(sr->common.ixn.gsft[0].Enu.lep_calo);
      }
      if(iscalo_reco){
        duneobj->rw_erec[i_event]=static_cast<double>(sr->common.ixn.gsft[0].Enu.lep_calo);
      }
      else {
        duneobj->rw_erec[i_event]=erec_total;
      }
      duneobj->rw_elep_reco[i_event] = elep_reco;
    }

    if(duneobj->rw_erec[i_event] != 0.){
      duneobj->rw_yrec[i_event] = ((duneobj->rw_erec[i_event])-(duneobj->rw_elep_reco[i_event]))/(duneobj->rw_erec[i_event]);
    }
    else {
      duneobj->rw_yrec[i_event] = static_cast<double>(0);
    }
    duneobj->rw_etru[i_event] = static_cast<double>(sr->mc.nu[0].E); // in GeV
    duneobj->rw_isCC[i_event] = static_cast<int>(sr->mc.nu[0].iscc);
    duneobj->rw_nuPDGunosc[i_event] = sr->mc.nu[0].pdgorig;
    duneobj->rw_nuPDG[i_event] = sr->mc.nu[0].pdg;
    duneobj->rw_berpaacvwgt[i_event] = _BeRPA_cvwgt;

    duneobj->nproton[i_event] = sr->mc.nu[0].nproton;
    duneobj->nneutron[i_event] = sr->mc.nu[0].nneutron;
    duneobj->npip[i_event] = sr->mc.nu[0].npip;
    duneobj->npim[i_event] = sr->mc.nu[0].npim;
    duneobj->npi0[i_event] = sr->mc.nu[0].npi0;
    duneobj->rw_vtx_x[i_event] = static_cast<double>(sr->mc.nu[0].vtx.x);
    duneobj->rw_vtx_y[i_event] = static_cast<double>(sr->mc.nu[0].vtx.y);
    duneobj->rw_vtx_z[i_event] = static_cast<double>(sr->mc.nu[0].vtx.z);

    double muonenergy = 0.;
    bool isEventAccepted = true;
    bool isgoodcafevent = true;

    int ntrueparticles = sr->mc.nu[0].nprim;

    for(int i_truepart =0; i_truepart<ntrueparticles; i_truepart++){
      int partpdg = sr->mc.nu[0].prim[i_truepart].pdg;
      double pdgmass = MaCh3Utils::GetMassFromPDG(partpdg);
      if(std::abs(partpdg) == 13) {
        duneobj->ntruemuon[i_event]++;
        duneobj->ntruemuonprim[i_event]++;
        if(static_cast<double>(sr->mc.nu[0].prim[i_truepart].p.E)>muonenergy){
          muonenergy = static_cast<double>(sr->mc.nu[0].prim[i_truepart].p.E);
          duneobj->rw_elep_true[i_event] = sr->mc.nu[0].prim[i_truepart].p.E - pdgmass;
          duneobj->rw_lep_pX[i_event] = static_cast<double>(sr->mc.nu[0].prim[i_truepart].p.px);
          duneobj->rw_lep_pY[i_event] = static_cast<double>(sr->mc.nu[0].prim[i_truepart].p.py);
          duneobj->rw_lep_pZ[i_event] = static_cast<double>(sr->mc.nu[0].prim[i_truepart].p.pz);
        }
      }

      bool isgoodcafparticle = false;
      if (!IsParticleAccepted(duneobj, iSample, i_event, i_truepart, pixel_spacing_cm, &isgoodcafparticle, pdgmass)) {
        isEventAccepted = false;
      }

      //Check if particle is stored properly in CAF
      if (isgoodcafparticle == false && isgoodcafevent == true) {
        isgoodcafevent = false;
        //MACH3LOG_INFO("Bad CAF event. Bad CAF particle has pdg = {}, start pos = ({},{},{}), start mom = ({},{},{})",
        //    partpdg, sr->mc.nu[0].prim[i_truepart].start_pos.x, sr->mc.nu[0].prim[i_truepart].start_pos.y, sr->mc.nu[0].prim[i_truepart].start_pos.z,
        //    sr->mc.nu[0].prim[i_truepart].p.px, sr->mc.nu[0].prim[i_truepart].p.py, sr->mc.nu[0].prim[i_truepart].p.pz);
      }
      else if (isgoodcafparticle == true && isgoodcafevent == true && i_truepart == ntrueparticles-1) {
        //MACH3LOG_INFO("Good CAF event. Good CAF particle has pdg = {}, start pos = ({},{},{}), start mom = ({},{},{})",
        //    partpdg, sr->mc.nu[0].prim[i_truepart].start_pos.x, sr->mc.nu[0].prim[i_truepart].start_pos.y, sr->mc.nu[0].prim[i_truepart].start_pos.z,
        //    sr->mc.nu[0].prim[i_truepart].p.px, sr->mc.nu[0].prim[i_truepart].p.py, sr->mc.nu[0].prim[i_truepart].p.pz);
      }
    }

    if(isEventAccepted) {duneobj->is_accepted[i_event]=1;}
    else {duneobj->is_accepted[i_event]=0;} 

    if (!isgoodcafevent) {
      nbadcafevents++;
      duneobj->is_good_caf_event[i_event]=0;
    }
    else {
      ngoodcafevents++;
      duneobj->is_good_caf_event[i_event]=1;
    }

    int ntruesecparticles = sr->mc.nu[0].nsec;
    for(int i_truepart =0; i_truepart<ntruesecparticles; i_truepart++){
      if(std::abs(sr->mc.nu[0].sec[i_truepart].pdg) == 13){
        duneobj->ntruemuon[i_event]++;
      }
    }

    duneobj->nmuonsratio[i_event] = static_cast<double>(duneobj->nrecomuon[i_event])/static_cast<double>(duneobj->ntruemuonprim[i_event]);
    duneobj->rw_rad[i_event] = radius;
    double reco_y = duneobj->rw_reco_vtx_y[i_event]-TPC_centre_y;
    double reco_z = duneobj->rw_reco_vtx_z[i_event]-TPC_centre_z;
    duneobj->rw_reco_rad[i_event] = std::sqrt(reco_y*reco_y + reco_z*reco_z);

    //Assume everything is on Argon for now....
    duneobj->Target[i_event] = 40;

    duneobj->rw_Q0[i_event] = static_cast<double>(sr->mc.nu[0].q0);
    duneobj->rw_Q3[i_event] = static_cast<double>(sr->mc.nu[0].modq);
    duneobj->rw_lep_pT[i_event] = std::sqrt((duneobj->rw_lep_pX[i_event])*(duneobj->rw_lep_pX[i_event]) + (duneobj->rw_lep_pY[i_event])*(duneobj->rw_lep_pY[i_event])); 

    //int M3Mode = Modes->GetModeFromGenerator(std::abs(sr->mc.nu[0].mode));
    duneobj->mode[i_event] = sr->mc.nu[0].mode;

    duneobj->flux_w[i_event] = 1.0;
    if(duneobj->rw_isCC[i_event] == 1) numCC++;
  }
  MACH3LOG_INFO("nEvents = {}, numCC = {}, numFDV = {}", duneobj->nEvents, numCC, num_in_fdv);
  MACH3LOG_INFO("ngoodcafevents = {}, nbadcafevents = {}", ngoodcafevents, nbadcafevents);
  _sampleFile->Close();
  _sampleFile_geant->Close();
  return duneobj->nEvents;
}

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wswitch-enum"
const double* samplePDFDUNEBeamNDGAr::GetPointerToKinematicParameter(KinematicTypes KinematicParameter, int iSample, int iEvent) {
  switch(KinematicParameter) {
    case kTrueNeutrinoEnergy:
      return &dunendgarmcSamples[iSample].rw_etru[iEvent]; 
    case kRecoNeutrinoEnergy:
      return &dunendgarmcSamples[iSample].rw_erec[iEvent];
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
    case kNMuonsRecoOverTruth:
      return &dunendgarmcSamples[iSample].nmuonsratio[iEvent];
    case kRecoLepEnergy:
      return &dunendgarmcSamples[iSample].rw_elep_reco[iEvent];
    case kTrueLepEnergy:
      return &dunendgarmcSamples[iSample].rw_elep_true[iEvent];
    case kRecoXPos:
      return &dunendgarmcSamples[iSample].rw_reco_vtx_x[iEvent];
    case kRecoYPos:
      return &dunendgarmcSamples[iSample].rw_reco_vtx_y[iEvent];
    case kRecoZPos:
      return &dunendgarmcSamples[iSample].rw_reco_vtx_z[iEvent];
    case kRecoRad:
      return &dunendgarmcSamples[iSample].rw_reco_rad[iEvent];
    case kLepPT:
      return &dunendgarmcSamples[iSample].rw_lep_pT[iEvent];
    case kLepPZ:
      return &dunendgarmcSamples[iSample].rw_lep_pZ[iEvent];
    case kTrueQ0:
      return &dunendgarmcSamples[iSample].rw_Q0[iEvent];
    case kTrueQ3:
      return &dunendgarmcSamples[iSample].rw_Q3[iEvent];
    case kParticle_Momentum:
      return &dunendgarmcSamples[iSample].particle_momentum->at(iEvent);
    case kParticle_TransverseMomentum:
      return &dunendgarmcSamples[iSample].particle_transversemomentum->at(iEvent);
    case kParticle_BAngle:
      return &dunendgarmcSamples[iSample].particle_bangle->at(iEvent);
    case kParticle_NHits:
      return &dunendgarmcSamples[iSample].particle_nhits->at(iEvent);
    case kParticle_NTurns:
      return &dunendgarmcSamples[iSample].particle_nturns->at(iEvent);
    case kParticle_MomResMS:
      return &dunendgarmcSamples[iSample].particle_momresms->at(iEvent);
    case kParticle_MomResTrans:
      return &dunendgarmcSamples[iSample].particle_momrestrans->at(iEvent);
    default:
      MACH3LOG_ERROR("Did not recognise Kinematic Parameter {}", KinematicParameter);
      throw MaCh3Exception(__FILE__, __LINE__);
      return nullptr;
  }
}
#pragma GCC diagnostic pop

const double* samplePDFDUNEBeamNDGAr::GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(std::round(KinematicVariable));
  return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

const double* samplePDFDUNEBeamNDGAr::GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

double samplePDFDUNEBeamNDGAr::ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(std::round(KinematicVariable));
  return ReturnKinematicParameter(KinPar,iSample,iEvent);
}

double samplePDFDUNEBeamNDGAr::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return ReturnKinematicParameter(KinPar,iSample,iEvent);
}

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wswitch-enum"
double samplePDFDUNEBeamNDGAr::ReturnKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent) {
  //HH: Special cases for dealing with non-doubles
  switch(KinPar) {
    case kEvent_IsAccepted:
      return static_cast<double>(dunendgarmcSamples[iSample].is_accepted[iEvent]);
    case kIsGoodCAFEvent:
      return static_cast<double>(dunendgarmcSamples[iSample].is_good_caf_event[iEvent]);
    case kIsCC:
      return static_cast<double>(dunendgarmcSamples[iSample].rw_isCC[iEvent]);
    case kInFDV:
      return static_cast<double>(dunendgarmcSamples[iSample].in_fdv[iEvent]);
    case kParticle_Event:
      return static_cast<double>(dunendgarmcSamples[iSample].particle_event->at(iEvent));
    case kParticle_IsAccepted:
      return static_cast<double>(dunendgarmcSamples[iSample].particle_isaccepted->at(iEvent));
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
    default:
      return *GetPointerToKinematicParameter(KinPar, iSample, iEvent);
  }
}
#pragma GCC diagnostic pop

void samplePDFDUNEBeamNDGAr::setupFDMC(int iSample) {
  dunemc_base *duneobj = &(dunendgarmcSamples[iSample]);
  FarDetectorCoreInfo *fdobj = &(MCSamples[iSample]);

  //fdobj->nutype = duneobj->nutype;
  //fdobj->oscnutype = duneobj->oscnutype;
  //fdobj->signal = duneobj->signal;
  //fdobj->SampleDetID = SampleDetID;

  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent){
    fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
    fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
    fdobj->Target[iEvent] = &(duneobj->Target[iEvent]);
  }
}

std::vector<double> samplePDFDUNEBeamNDGAr::ReturnKinematicParameterBinning(std::string KinematicParameterStr) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameterStr));
  return ReturnKinematicParameterBinning(KinPar);
}

#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wswitch-enum"
std::vector<double> samplePDFDUNEBeamNDGAr::ReturnKinematicParameterBinning(KinematicTypes KinPar) {
  std::vector<double> binningVector;
  switch(KinPar){
    case kTrueNeutrinoEnergy:
      for(double ibins =0; ibins<10*10; ibins++){
        double binval = ibins/10;
        binningVector.push_back(binval);
      }
      break;
    case kRecoNeutrinoEnergy:
      for(double ibins =0; ibins<10*10; ibins++){
        double binval = ibins/10;
        binningVector.push_back(binval);
      } 
      break;
    case kRecoXPos:
    case kTrueXPos:
      for(double ibins =0; ibins<259*2; ibins++){
        binningVector.push_back(ibins-259);
      }
      break;
    case kRecoYPos:
    case kTrueYPos:
      for(double ibins =0; ibins<277*2; ibins++){
        binningVector.push_back(ibins-277-150);
      }
      break;
    case kRecoZPos:
    case kTrueZPos:
      for(double ibins =0; ibins<277*2; ibins++){
        binningVector.push_back(ibins-277+1486);
      }
      break;
    case kNMuonsRecoOverTruth:
      for(double ibins =0; ibins<20*10; ibins++){
        binningVector.push_back(-10+ibins/10);
      }
      break;
    case kRecoLepEnergy:
      for(double ibins =0; ibins<10*10; ibins++){
        binningVector.push_back(ibins/10);
      } 
      break;
    case kTrueLepEnergy:
      for(double ibins =0; ibins<10*10; ibins++){
        binningVector.push_back(ibins/10);
      } 
      break;
    case kTrueRad:
    case kRecoRad:
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

TH1* samplePDFDUNEBeamNDGAr::get1DParticleVarHist(std::string ProjectionVar_Str, std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* Axis) {
  //DB Grab the associated enum with the argument string
  int ProjectionVar_Int = ReturnKinematicParameterFromString(ProjectionVar_Str);

  //DB Need to overwrite the Selection member variable so that IsEventSelected function operates correctly.
  //   Consequently, store the selection cuts already saved in the sample, overwrite the Selection variable, then reset
  std::vector< std::vector<double> > tmp_Selection = Selection;
  std::vector< std::vector<double> > SelectionVecToApply;

  //DB Add all the predefined selections to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<Selection.size();iSelec++) {
    SelectionVecToApply.emplace_back(Selection[iSelec]);
  }

  //DB Add all requested cuts from the argument to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<SelectionVec.size();iSelec++) {
    SelectionVecToApply.emplace_back(SelectionVec[iSelec]);
  }

  //DB Check the formatting of all requested cuts, should be [cutPar,lBound,uBound]
  for (size_t iSelec=0;iSelec<SelectionVecToApply.size();iSelec++) {
    if (SelectionVecToApply[iSelec].size()!=3) {
      MACH3LOG_ERROR("Selection Vector[{}] is not formed correctly. Expect size == 3, given: {}",iSelec,SelectionVecToApply[iSelec].size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
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

  //Loop over all particles
  for (int iSample=0;iSample<getNMCSamples();iSample++) {
    for (int iParticle=0;iParticle<nparticlesinsample[iSample];iParticle++) {
      int event = static_cast<int>(std::round(ReturnKinematicParameter("Particle_Event", iSample, iParticle)));
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

TH2* samplePDFDUNEBeamNDGAr::get2DParticleVarHist(std::string ProjectionVar_StrX, std::string ProjectionVar_StrY, std::vector< std::vector<double> > SelectionVec, int WeightStyle, TAxis* AxisX, TAxis* AxisY) {
  //DB Grab the associated enum with the argument string
  int ProjectionVar_IntX = ReturnKinematicParameterFromString(ProjectionVar_StrX);
  int ProjectionVar_IntY = ReturnKinematicParameterFromString(ProjectionVar_StrY);

  //DB Need to overwrite the Selection member variable so that IsEventSelected function operates correctly.
  //   Consequently, store the selection cuts already saved in the sample, overwrite the Selection variable, then reset
  std::vector< std::vector<double> > tmp_Selection = Selection;
  std::vector< std::vector<double> > SelectionVecToApply;

  //DB Add all the predefined selections to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<Selection.size();iSelec++) {
    SelectionVecToApply.emplace_back(Selection[iSelec]);
  }

  //DB Add all requested cuts from the argument to the selection vector which will be applied
  for (size_t iSelec=0;iSelec<SelectionVec.size();iSelec++) {
    SelectionVecToApply.emplace_back(SelectionVec[iSelec]);
  }

  //DB Check the formatting of all requested cuts, should be [cutPar,lBound,uBound]
  for (size_t iSelec=0;iSelec<SelectionVecToApply.size();iSelec++) {
    if (SelectionVecToApply[iSelec].size()!=3) {
      MACH3LOG_ERROR("Selection Vector[{}] is not formed correctly. Expect size == 3, given: {}",iSelec,SelectionVecToApply[iSelec].size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
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
  for (int iSample=0;iSample<getNMCSamples();iSample++) {
    for (int iParticle=0;iParticle<nparticlesinsample[iSample];iParticle++) {
      int event = static_cast<int>(std::round(ReturnKinematicParameter("Particle_Event", iSample, iParticle)));
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

bool samplePDFDUNEBeamNDGAr::IsParticleSelected(const int iSample, const int iEvent, const int iParticle) {
  double Val;
  for (unsigned int iSelection=0;iSelection < Selection.size() ;iSelection++) {
    std::string SelecVarStr = ReturnStringFromKinematicParameter(static_cast<int>(std::round(Selection[iSelection][0])));
    int SelectionVarElement;

    if (SelecVarStr.find("Particle_") != std::string::npos) {SelectionVarElement = iParticle;}
    else {SelectionVarElement = iEvent;}

    Val = ReturnKinematicParameter(Selection[iSelection][0], iSample, SelectionVarElement);
    if ((Val<Selection[iSelection][1])||(Val>=Selection[iSelection][2])) {
      return false;
    }
  }
  return true;
}
