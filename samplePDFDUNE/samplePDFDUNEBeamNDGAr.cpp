#include "samplePDFDUNEBeamNDGAr.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

samplePDFDUNEBeamNDGAr::samplePDFDUNEBeamNDGAr(std::string mc_version_, covarianceXsec* XsecCov_) : samplePDFFDBase(mc_version_, XsecCov_) {
	KinematicParameters = &KinematicParametersDUNE;
	ReversedKinematicParameters = &ReversedKinematicParametersDUNE;

	Initialise();
}

samplePDFDUNEBeamNDGAr::~samplePDFDUNEBeamNDGAr() {
}

void samplePDFDUNEBeamNDGAr::Init() {
	dunendgarmcSamples.resize(nSamples,dunemc_base());
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

double samplePDFDUNEBeamNDGAr::FindNHits(double pixel_spacing_cm, double centre_circle_y, double centre_circle_z, double rad_curvature){
	//use the pixel grid method to find number of pixels hit in a track.
	int num_vertices =0; // number of vertices hit. Counting to avoid duplicating
	int num_intersections =0;

	//equation for circle = (y-y0)^2 + (z-z0)^2 = r^2

	for(unsigned int i_intersect = 0; i_intersect<yboundarypositions.size(); i_intersect++){ //check every boundary line to see if it crossed within the TPC instrumented region
		double quadratic_ineq_y= pow(rad_curvature*100, 2)-pow((yboundarypositions[i_intersect]-centre_circle_y), 2);
		if(quadratic_ineq_y > 0){
			if((pow((pow((centre_circle_z + pow(quadratic_ineq_y, 0.5)-TPC_centre_z), 2) + pow((yboundarypositions[i_intersect]-TPC_centre_y), 2)), 0.5) <=TPCInstrumentedRadius) && (pixelzmin<(centre_circle_z + pow(quadratic_ineq_y, 0.5))) && (centre_circle_z + pow(quadratic_ineq_y, 0.5)<=pixelzmax)){ //check that the z coord is also on pixel plane
				num_intersections++;
				if(fmod((centre_circle_z + pow(quadratic_ineq_y, 0.5) - pixelzmin), pixel_spacing_cm) == 0){ // this is the case when a vertex is crossed so to avoid double counting pixels
					num_vertices++;
				}
			}
			if((pow((pow((centre_circle_z - pow(quadratic_ineq_y, 0.5)-TPC_centre_z), 2) + pow((yboundarypositions[i_intersect]-TPC_centre_y), 2)), 0.5) <TPCInstrumentedRadius) && pixelzmin<(centre_circle_z - pow(quadratic_ineq_y, 0.5)) && (centre_circle_z - pow(quadratic_ineq_y, 0.5))<=pixelzmax){ //check that the z coord is also on pixel plane
				num_intersections++;
				if(fmod((centre_circle_z - pow(quadratic_ineq_y, 0.5) - pixelzmin), pixel_spacing_cm) == 0){ // this is the case when a vertex is crossed so to avoid double counting pixels
					num_vertices++;
				}
			}
		}
		else if(quadratic_ineq_y == 0){ //when the pixel boundary is a tangent to the circle z = z0
			num_intersections++;
		}
	}
	for(unsigned int i_intersect = 0; i_intersect<zboundarypositions.size(); i_intersect++){
		double quadratic_ineq_z = pow(rad_curvature*100, 2)-pow((zboundarypositions[i_intersect]-centre_circle_z), 2);
		if(quadratic_ineq_z > 0){
			if( (pow((pow((centre_circle_y + pow(quadratic_ineq_z, 0.5)-TPC_centre_y), 2) + pow((zboundarypositions[i_intersect]-TPC_centre_z), 2)), 0.5) <=TPCInstrumentedRadius) && pixelymin<(centre_circle_y + pow(quadratic_ineq_z, 0.5)) && centre_circle_y + pow(quadratic_ineq_z, 0.5)<=pixelymax){ //check that the z coord is also on pixel plane
				num_intersections++;
			}
			if((pow((pow((centre_circle_y - pow(quadratic_ineq_z, 0.5)-TPC_centre_y), 2) + pow((zboundarypositions[i_intersect]-TPC_centre_z), 2)), 0.5) <=TPCInstrumentedRadius) && pixelymin<(centre_circle_y - pow(quadratic_ineq_z, 0.5)) && centre_circle_y - pow(quadratic_ineq_z, 0.5)<=pixelymax){ //check that the z coord is also on pixel plane
				num_intersections++;
			}
			// already checked all vertices for duplicates before so no need to repeat that
		}
		else if(quadratic_ineq_z == 0){ //when the pixel boundary is a tangent to the circle y = y0
			num_intersections++;
		}
	}
	double N_hits = num_intersections - num_vertices + 1; //Add one for the pixel that it starts on
	return N_hits;
}

double samplePDFDUNEBeamNDGAr::CalcBeta(double p_mag, double& bg, double& gamma){ //calculate beta (v/c)
	bg = p_mag/pdgmass; //beta*gamma
	gamma = pow((1+bg*bg), 0.5); //gamma
	double beta = bg/gamma; //beta (velocity)
	return beta;
}

double samplePDFDUNEBeamNDGAr::GetMass(int partpdg){
	double mass=0;
	switch(std::abs(partpdg)) {
		case 13: 
			mass = m_mu; break;
		case 211: 
			mass = m_chargedpi; break;
		case 111: 
			mass = m_pi0; break;
		case 2212: 
			mass = m_p; break;
		case 2112: 
			mass = m_n; break;
		case 11: 
			mass = m_e; break;
		case 321: 
			mass = m_chargedk; break;
		case 311:
			mass = m_k0; break;
		case 3122:
			mass = m_lambda; break;
		case 3222: //sig+
			mass = 1.118937; break;
		case 3112: //sig-
			mass = 1.197449; break;
		case 3212: //sig0
			mass = 1.192642; break;
		case 1000180400: //Ar40 
		  mass = 37.27272; break;
		case 1000060120: //C12
			mass = 11.17793; break;
		case 1000080160: //O16
			mass = 14.89589; break;
		case 1000190390: //K39
			mass = 36.32719; break;
		case 1000110230: //Na23
			mass = 21.45502; break;
		case 1000300640: //Zn64
			mass = 59.60401; break;
		case 1000140280: //Si28
			mass = 26.19938; break;
		case 1000070140: //N14
			mass = 13043.78; break;
		case 130: //K_l
    case 310: //K_s
			mass = 0.497611; break;
		case 12: //nue
		case 14: //numu
		case 16: //nutau
		case 22: //gamma
			mass = 0.; break;
		default:
			MACH3LOG_INFO("Particle mass not assigned: PDG = {}", partpdg);
			break;
	}
	return mass;
}

bool samplePDFDUNEBeamNDGAr::IsParticleAccepted(dunemc_base *duneobj, int i_event, int i_truepart, double pixel_spacing_cm){
	static std::random_device rd;
  static std::mt19937 gen(rd());
	static std::uniform_real_distribution<double> dist(0.0, 1.0);
	
	bool isAccepted = true;
	int nummatched = 0;

	for(unsigned int i_anapart =0; i_anapart<_MCPStartPX->size(); i_anapart++){

		double mom_tot = pow(pow(_MCPStartPX->at(i_anapart), 2)+pow(_MCPStartPY->at(i_anapart), 2)+pow(_MCPStartPZ->at(i_anapart), 2), 0.5);
		//Now select particle in the anatree that is the same as in the caf (same pdg and momentum)
		//if((_PDG->at(i_anapart) == sr->mc.nu[0].prim[i_truepart].pdg) && (mom_tot >= 0.999*static_cast<double>(sr->mc.nu[0].prim[i_truepart].p.Mag())) && (mom_tot <= 1.001*static_cast<double>(sr->mc.nu[0].prim[i_truepart].p.Mag()))){			
		if((_PDG->at(i_anapart) == sr->mc.nu[0].prim[i_truepart].pdg) && 
			(std::abs(_MCPStartPX->at(i_anapart)-static_cast<double>(sr->mc.nu[0].prim[i_truepart].p.px))<=0.001*std::abs(_MCPStartPX->at(i_anapart))) &&
			(std::abs(_MCPStartPY->at(i_anapart)-static_cast<double>(sr->mc.nu[0].prim[i_truepart].p.py))<=0.001*std::abs(_MCPStartPY->at(i_anapart))) &&
			(std::abs(_MCPStartPZ->at(i_anapart)-static_cast<double>(sr->mc.nu[0].prim[i_truepart].p.pz))<=0.001*std::abs(_MCPStartPZ->at(i_anapart)))){			
			nummatched++;
			//Ignore neutrons and neutrinos
			if((std::abs(sr->mc.nu[0].prim[i_truepart].pdg))!= 2112 && (std::abs(sr->mc.nu[0].prim[i_truepart].pdg))!= 14 && (std::abs(sr->mc.nu[0].prim[i_truepart].pdg))!= 12){
				bool stopsinecal_radius = false;
				bool stopsinecal_length = false;
				double end_radius = pow((pow(_MCPEndY->at(i_anapart)-TPC_centre_y, 2)+pow(_MCPEndZ->at(i_anapart)-TPC_centre_z, 2)), 0.5);
				double end_length = _MCPEndX->at(i_anapart)-TPC_centre_x;
				if(ecal_containment && end_radius>ECALInnerRadius && end_radius<ECALOuterRadius){stopsinecal_radius = true;}
				if(ecal_containment && std::abs(end_length)>ECALEndCapStart && std::abs(end_length)<ECALEndCapEnd){stopsinecal_length = true;}

				//Momentum transverse to B field
				double transverse_mom = pow(pow(_MCPStartPY->at(i_anapart), 2)+pow(_MCPStartPZ->at(i_anapart), 2), 0.5);
				double start_radius =pow((pow(_MCPStartY->at(i_anapart)-TPC_centre_y, 2)+pow(_MCPStartZ->at(i_anapart)-TPC_centre_z, 2)), 0.5);

				//Fill particle-level kinematic variables
				duneobj->particle[i_event].pdg[i_truepart]=_PDG->at(i_anapart);
				duneobj->particle[i_event].energy[i_truepart]=static_cast<double>(sr->mc.nu[0].prim[i_truepart].p.E);
				duneobj->particle[i_event].momentum[i_truepart]=mom_tot;
				//particleBAngle defined as angle wrt B-field (x)
				duneobj->particle[i_event].bangle[i_truepart]=acos(_MCPStartPX->at(i_anapart)/mom_tot)*180/M_PI;
				//If particle is not stopped in the tpc or ecal 
				if((std::abs(end_length)>TPCInstrumentedLength && !stopsinecal_length) || (end_radius>TPCInstrumentedRadius && !stopsinecal_radius)){
					//If stopped between tpc and ecal, determine where
					if(end_radius<ECALOuterRadius && std::abs(end_length)<ECALEndCapEnd) {
						duneobj->particle[i_event].isstoppedingap[i_truepart] = true;
						if (std::abs(end_length)>TPCInstrumentedLength) {duneobj->particle[i_event].isstoppedinendgap[i_truepart] = true;}
						else {duneobj->particle[i_event].isstoppedinbarrelgap[i_truepart] = true;}
					}
					//Check if originated inside fiducial volume
					if((std::abs(_MCPStartX->at(i_anapart))-TPC_centre_x)<=TPCFidLength && start_radius<=TPCFidRadius){
						//Check if charged (p +/- , pi +/- , mu +/- , e +/- , K +/-)
						if(std::abs(_PDG->at(i_anapart)) == 2212 || std::abs(_PDG->at(i_anapart)) == 211 || std::abs(_PDG->at(i_anapart)) == 13 || std::abs(_PDG->at(i_anapart)) == 11 || std::abs(_PDG->at(i_anapart)) == 321) {

							double length_track_x;
							if(std::abs(end_length)>TPCInstrumentedLength){
								if((end_length)>=0){ length_track_x = TPCInstrumentedLength - (_MCPStartX->at(i_anapart)-TPC_centre_x);} //in cm
								else{ length_track_x = -TPCInstrumentedLength - (_MCPStartX->at(i_anapart)-TPC_centre_x);} //in cm
							}
							else{length_track_x = _MCPEndX->at(i_anapart) - _MCPStartX->at(i_anapart);} //in cm

							double rad_curvature = transverse_mom/(0.3*B_field); //p = 0.3*B*r where p in GeV/c, B in T, r in m
							double theta_xT = atan(_MCPStartPX->at(i_anapart)/transverse_mom); //helix pitch angle
							double pitch = std::abs(2*M_PI*rad_curvature*tan(theta_xT)); //distance between two turns of a helix in m
							double tan_theta = tan(theta_xT);

							//Find centre of circular path
							bool positivecharged =0;
							double centre_circle_y;
							double centre_circle_z;
							double L_yz, L_yz_chord; //length of curved track in y-z plane
							if(_PDG->at(i_anapart) == 2212 || _PDG->at(i_anapart) == 211 || _PDG->at(i_anapart) == -13 || _PDG->at(i_anapart) == -11 || _PDG->at(i_anapart) == 321){ positivecharged = 1;}
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
							double a_const = (pow(TPCInstrumentedRadius, 2)-pow(rad_curvature*100, 2) - (pow(TPC_centre_y, 2)-pow(centre_circle_y, 2))-(pow(TPC_centre_z, 2)-pow(centre_circle_z, 2)))/(2*(centre_circle_y-TPC_centre_y));
							double quadraticformula_b = -(2*m_const*(a_const -TPC_centre_y)+2*TPC_centre_z);
							double quadraticformula_a = pow(m_const, 2)+1;
							double quadraticformula_c = pow((a_const - TPC_centre_y), 2) +pow(TPC_centre_z, 2) - pow(TPCInstrumentedRadius,2);

							double z_intersect_1, y_intersect_1, z_intersect_2, y_intersect_2, theta_1, theta_2, theta_start, theta_chosen, theta_diff_1, theta_diff_2;
							if(pow(quadraticformula_b, 2)-4*quadraticformula_a*quadraticformula_c>0){
								z_intersect_1 = (-quadraticformula_b+pow((pow(quadraticformula_b, 2)-4*quadraticformula_a*quadraticformula_c), 0.5))/(2*quadraticformula_a);
								y_intersect_1 = -m_const*z_intersect_1+a_const;
								z_intersect_2 = (-quadraticformula_b-pow((pow(quadraticformula_b, 2)-4*quadraticformula_a*quadraticformula_c), 0.5))/(2*quadraticformula_a);
								y_intersect_2 = -m_const*z_intersect_2+a_const;
								theta_1 = std::abs(atan((y_intersect_1-centre_circle_y)/(z_intersect_1-centre_circle_z)));
								theta_2 = std::abs(atan((y_intersect_2-centre_circle_y)/(z_intersect_2-centre_circle_z)));
								theta_start = std::abs(atan((_MCPStartY->at(i_anapart)-centre_circle_y)/(_MCPStartZ->at(i_anapart)-centre_circle_z)));
								//Adjust angles based on quadrants
								if((z_intersect_2-centre_circle_z)<0 && (y_intersect_2-centre_circle_y)<0){theta_2 = M_PI+theta_2;}
								else if((y_intersect_2-centre_circle_y)<0){theta_2 = 2*M_PI-theta_2;}
								else if((z_intersect_2-centre_circle_z)<0){theta_2 = M_PI - theta_2;}

								if((z_intersect_1-centre_circle_z)<0 && (y_intersect_1-centre_circle_y)<0){theta_1 = M_PI+theta_1;}
								else if((y_intersect_1-centre_circle_y)<0){theta_1 = 2*M_PI-theta_1;}
								//JM corrected typo below (theta_2->theta_1)
								else if((z_intersect_1-centre_circle_z)<0){theta_1 = M_PI - theta_1;}

								if((_MCPStartZ->at(i_anapart)-centre_circle_z)<0 && (_MCPStartY->at(i_anapart)-centre_circle_y)<0){theta_start = M_PI+theta_start;}
								else if((_MCPStartY->at(i_anapart)-centre_circle_y)<0){theta_start = 2*M_PI-theta_start;}
								else if((_MCPStartZ->at(i_anapart)-centre_circle_z)<0){theta_start = M_PI - theta_start;}

								//Lorentz force law, if positively charged, counter clockwise, if negative charged, clockwise
								if(!positivecharged){
									theta_diff_1 = (theta_1 < theta_start) ? (theta_start - theta_1) : (2*M_PI - (theta_1 - theta_start));
									theta_diff_2 = (theta_2 < theta_start) ? (theta_start - theta_2) : (2*M_PI - (theta_2 - theta_start));
								}
								else if(positivecharged){ 
									theta_diff_1 = (theta_1 > theta_start) ? (theta_1 - theta_start) : (2*M_PI - (theta_start - theta_1));
									theta_diff_2 = (theta_2 > theta_start) ? (theta_2 - theta_start) : (2*M_PI - (theta_start - theta_2));
								}

								if(theta_diff_1<theta_diff_2 && (rad_curvature*100*theta_diff_1 > (TPCInstrumentedRadius-start_radius))){
									theta_chosen = theta_diff_1; 
								} 
								else {
									theta_chosen = theta_diff_2; 
								}
								L_yz = rad_curvature*100*theta_chosen;
								L_yz_chord = std::abs(2*rad_curvature*100*sin(theta_chosen/2));
								if(std::abs(L_yz*tan_theta) > std::abs(length_track_x)){
									double nturns = std::abs(length_track_x/100)/pitch;
									double theta_intersect = fmod(nturns, 1)*2*M_PI;
									L_yz_chord = std::abs(2*rad_curvature*100*sin(theta_intersect/2));
								} 
								else {
									length_track_x = std::abs(L_yz*tan_theta);
								}
							}
							else {
								double nturns = (length_track_x/100)/pitch;
								double theta_intersect;
								if (nturns>=1) {theta_intersect = 2*M_PI;}
								else {theta_intersect = nturns*2*M_PI;}
								L_yz_chord = std::abs(2*rad_curvature*100*sin(theta_intersect/2));
								L_yz = rad_curvature*100*theta_intersect; //JM this was previously uninitialised in this case
								theta_chosen = theta_intersect;
							}

							double N_hits = FindNHits(pixel_spacing_cm, centre_circle_y, centre_circle_z, rad_curvature);
							double p_mag = sr->mc.nu[0].prim[i_truepart].p.Mag();
							double bg = 0; double gamma = 0;
							double beta = CalcBeta(p_mag, bg, gamma);

							//JM removed avg_betapT calc (seems to be unused)

							double sigmax = (drift_velocity/100)*(1/(adc_sampling_frequency));
							double sigmayz = (spatial_resolution/(1000)); //needs to be in m              
							double momres_yz = transverse_mom*(pow(720/(N_hits+4), 0.5)*(sigmayz*transverse_mom/(0.3*B_field*pow(L_yz_chord/100, 2)))*pow((1-(1/21)*pow((L_yz_chord/(rad_curvature*100)), 2)), 0.5));
							double momres_ms = transverse_mom*(0.016/(0.3*B_field*(L_yz/100)*cos(theta_xT)*beta))*pow(L_yz/X0, 0.5);
							double momres_tottransverse = pow(pow(momres_yz, 2) + pow(momres_ms, 2), 0.5);
							double sigma_theta = (pow(cos(theta_xT), 2))*(pitch/(2*M_PI*rad_curvature))*pow((pow(sigmax/(std::abs(length_track_x)/100),2) +pow(momres_tottransverse/transverse_mom, 2)), 0.5);
							double momres_frac = pow(pow((momres_tottransverse/transverse_mom), 2)+pow(sigma_theta*tan_theta, 2), 0.5);
							duneobj->particle[i_event].momresms[i_truepart] = momres_ms;
							duneobj->particle[i_event].momrestransfrac[i_truepart] = (momres_tottransverse/transverse_mom)/momres_frac;
							duneobj->particle[i_event].momrestrans[i_truepart] = momres_tottransverse/transverse_mom;
							if(momres_frac > momentum_resolution_threshold){
								isAccepted = false;
							}
						}
						else{ //Accept pi0 and gamma with some efficiency
							if((std::abs(sr->mc.nu[0].prim[i_truepart].pdg)) == 111){
								double pi0recoprob = dist(gen);
								if(pi0recoprob>pi0_reco_efficiency){
									isAccepted = false;
								}
							}
							else if((std::abs(sr->mc.nu[0].prim[i_truepart].pdg)) == 22){
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
					else { //Originates outside fiducial volume
						isAccepted = false; 
					}
				}
				else if ((std::abs(end_length)>TPCInstrumentedLength) || (end_radius>TPCInstrumentedRadius)){//Stops in ECal
					double energydepsum = 0;
					for(unsigned int i_ecaldep=0; i_ecaldep<_SimHitEnergy->size(); i_ecaldep++){
						if(_SimHitTrkID->at(i_ecaldep) == _MCPTrkID->at(i_anapart)){
							energydepsum = energydepsum + _SimHitEnergy->at(i_ecaldep);
						}
					}
					duneobj->particle[i_event].ecaldepositfraction[i_truepart] = energydepsum/(sr->mc.nu[0].prim[i_truepart].p.E);
					duneobj->particle[i_event].isstoppedinecal[i_truepart] = true;
					if(end_radius>TPCInstrumentedRadius) {duneobj->particle[i_event].isstoppedinbarrel[i_truepart] = true;}
					else {duneobj->particle[i_event].isstoppedinendcap[i_truepart] = true;}
				}
				else {duneobj->particle[i_event].isstoppedintpc[i_truepart] = true;}
				if(isAccepted) {
					//fill particle-level accepted kinematic parameters
					duneobj->particle[i_event].accmomentum[i_truepart]=mom_tot;
					duneobj->particle[i_event].accbangle[i_truepart]=acos(_MCPStartPX->at(i_anapart)/mom_tot)*180/M_PI;
					duneobj->particle[i_event].accpdg[i_truepart]=_PDG->at(i_anapart);
				}
				else {
					//fill particle-level rejected kinematic parameters
					duneobj->particle[i_event].rejx[i_truepart] = _MCPStartX->at(i_anapart)-TPC_centre_x;
					duneobj->particle[i_event].rejr2[i_truepart] = start_radius*start_radius;
				}	
			}
			//break;
		}
	}
	if(nummatched != 1){
		MACH3LOG_INFO("Found {} matching particles in anatree", nummatched);
		MACH3LOG_INFO("PDG: {}, momentum: ({},{},{})",sr->mc.nu[0].prim[i_truepart].pdg,sr->mc.nu[0].prim[i_truepart].p.px,sr->mc.nu[0].prim[i_truepart].p.py,sr->mc.nu[0].prim[i_truepart].p.pz);
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
	double downsampling = 0.01;
	duneobj->pot_s = pot/(downsampling*1e21);
	duneobj->nEvents = static_cast<int>(downsampling*static_cast<double>(_data->GetEntries()));

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

	//Particle-level kinematic variables
	duneobj->particle = new particlevar[duneobj->nEvents];
	//int ntrueparticles = sr->mc.nu[0].nprim;
	/*int ntrueparticles = 7*duneobj->nEvents;
		duneobj->particleecaldepositfraction = new double[ntrueparticles];
		duneobj->particleevent = new int[ntrueparticles];
		duneobj->particlepdg = new int[ntrueparticles];
		duneobj->accparticlepdg = new int[ntrueparticles];
		duneobj->particleenergy = new double[ntrueparticles];
		duneobj->particlemomentum = new double[ntrueparticles];
		duneobj->particleBAngle = new double[ntrueparticles];
		duneobj->accparticleBAngle = new double[ntrueparticles];
		duneobj->accparticlemomentum = new double[ntrueparticles];
		duneobj->isparticlestoppedintpc = new bool[ntrueparticles];
		duneobj->isparticlestoppedinecal = new bool[ntrueparticles];
		duneobj->isparticlestoppedingap = new bool[ntrueparticles];
		duneobj->isparticlestoppedinbarrelgap = new bool[ntrueparticles];
		duneobj->isparticlestoppedinendgap = new bool[ntrueparticles];
		duneobj->isparticlestoppedinbarrel= new bool[ntrueparticles];
		duneobj->isparticlestoppedinendcap = new bool[ntrueparticles];
		duneobj->rejparticlex = new double[ntrueparticles];
		duneobj->rejparticler2 = new double[ntrueparticles];
		duneobj->particlemomresms = new double[ntrueparticles];
		duneobj->particlemomrestransfrac = new double[ntrueparticles];
		duneobj->particlemomrestrans= new double[ntrueparticles];*/

	int num_no_ixns =0;
	int num_no_recparticles = 0;
	int num_in_fdv = 0;
	int num_in_fdv_noreco = 0;
	int num_notin_fdv =0;
	int num_nanenergy =0;
	int num_nanparticles =0;

	double pixel_spacing_cm = pixel_spacing/10; //convert to cm
	makePixelGrid(pixel_spacing_cm); //make square pixel grid and fill vectors with y and z pixel boundaries	

	for (int i_event = 0; i_event < (duneobj->nEvents); ++i_event) { 
		_data->GetEntry(i_event);

		if ((i_event % (duneobj->nEvents/100))==0) {
			MACH3LOG_INFO("\tProcessing event: {}/{}",i_event,duneobj->nEvents);
		}
		double radius = pow((pow((sr->mc.nu[0].vtx.y-TPC_centre_y),2) + pow((sr->mc.nu[0].vtx.z-TPC_centre_z),2)),0.5); //find radius of interaction vertex
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

		double muonenergy = 0.;
		bool isEventAccepted = true;

		int ntrueparticles = sr->mc.nu[0].nprim;

		//allocate memory to particle-level kinematic variables
		duneobj->particle[i_event].nparticles = ntrueparticles;
		duneobj->particle[i_event].ecaldepositfraction = new double[ntrueparticles];
		duneobj->particle[i_event].pdg = new int[ntrueparticles];
		duneobj->particle[i_event].accpdg = new int[ntrueparticles];
		duneobj->particle[i_event].energy = new double[ntrueparticles];
		duneobj->particle[i_event].momentum = new double[ntrueparticles];
		duneobj->particle[i_event].bangle = new double[ntrueparticles];
		duneobj->particle[i_event].accbangle = new double[ntrueparticles];
		duneobj->particle[i_event].accmomentum = new double[ntrueparticles];
		duneobj->particle[i_event].isstoppedintpc = new bool[ntrueparticles];
		duneobj->particle[i_event].isstoppedinecal = new bool[ntrueparticles];
		duneobj->particle[i_event].isstoppedingap = new bool[ntrueparticles];
		duneobj->particle[i_event].isstoppedinbarrelgap = new bool[ntrueparticles];
		duneobj->particle[i_event].isstoppedinendgap = new bool[ntrueparticles];
		duneobj->particle[i_event].isstoppedinbarrel= new bool[ntrueparticles];
		duneobj->particle[i_event].isstoppedinendcap = new bool[ntrueparticles];
		duneobj->particle[i_event].rejx = new double[ntrueparticles];
		duneobj->particle[i_event].rejr2 = new double[ntrueparticles];
		duneobj->particle[i_event].momresms = new double[ntrueparticles];
		duneobj->particle[i_event].momrestransfrac = new double[ntrueparticles];
		duneobj->particle[i_event].momrestrans= new double[ntrueparticles];


		for(int i_truepart =0; i_truepart<ntrueparticles; i_truepart++){
			int partpdg = sr->mc.nu[0].prim[i_truepart].pdg;
			pdgmass = GetMass(partpdg);
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
			if(!IsParticleAccepted(duneobj, i_event, i_truepart, pixel_spacing_cm)) {
				isEventAccepted = false;
			}
		}

		if(isEventAccepted) {duneobj->is_accepted[i_event]=1;}
		else {duneobj->is_accepted[i_event]=0;} 

		int ntruesecparticles = sr->mc.nu[0].nsec;
		for(int i_truepart =0; i_truepart<ntruesecparticles; i_truepart++){
			if(std::abs(sr->mc.nu[0].sec[i_truepart].pdg) == 13){
				duneobj->ntruemuon[i_event]++;
			}
		}

		duneobj->nmuonsratio[i_event] = static_cast<double>(duneobj->nrecomuon[i_event])/static_cast<double>(duneobj->ntruemuonprim[i_event]);
		duneobj->rw_vtx_x[i_event] = static_cast<double>(sr->mc.nu[0].vtx.x);
		duneobj->rw_vtx_y[i_event] = static_cast<double>(sr->mc.nu[0].vtx.y);
		duneobj->rw_vtx_z[i_event] = static_cast<double>(sr->mc.nu[0].vtx.z);

		duneobj->rw_rad[i_event] = pow((pow((duneobj->rw_vtx_y[i_event]+150),2) + pow((duneobj->rw_vtx_z[i_event]-1486),2)),0.5);
		duneobj->rw_reco_rad[i_event] = pow(pow((duneobj->rw_reco_vtx_y[i_event]+150),2) + pow((duneobj->rw_reco_vtx_z[i_event]-1486), 2), 0.5);

		//Assume everything is on Argon for now....
		duneobj->Target[i_event] = 40;

		duneobj->rw_Q0[i_event] = static_cast<double>(sr->mc.nu[0].q0);
		duneobj->rw_Q3[i_event] = static_cast<double>(sr->mc.nu[0].modq);
		duneobj->rw_lep_pT[i_event] = (pow(pow(duneobj->rw_lep_pX[i_event], 2) + pow(duneobj->rw_lep_pY[i_event], 2), 0.5)); 

		_mode = sr->mc.nu[0].mode;
		_isCC = static_cast<int>(sr->mc.nu[0].iscc);

		int M3Mode = Modes->GetModeFromGenerator(std::abs(sr->mc.nu[0].mode));
		duneobj->mode[i_event] = M3Mode;
		//int mode= TMath::Abs(_mode);       
		//duneobj->mode[i_event]=static_cast<double>GENIEMode_ToMaCh3Mode(mode, _isCC);

		duneobj->flux_w[i_event] = 1.0;
	}

	_sampleFile->Close();
	_sampleFile_geant->Close();
	return duneobj->nEvents;
}

const double* samplePDFDUNEBeamNDGAr::GetPointerToKinematicParameter(KinematicTypes KinematicParameter, int iSample, int iEvent) {
	double* KinematicValue;

	switch(KinematicParameter) {
		case kTrueNeutrinoEnergy:
			KinematicValue = &dunendgarmcSamples[iSample].rw_etru[iEvent]; 
			break;
		case kRecoNeutrinoEnergy:
			KinematicValue = &dunendgarmcSamples[iSample].rw_erec[iEvent];
			break;
		case kMode:
			KinematicValue = &dunendgarmcSamples[iSample].mode[iEvent];
			break;
		case kTrueXPos:
			KinematicValue = &dunendgarmcSamples[iSample].rw_vtx_x[iEvent];
			break;
		case kTrueYPos:
			KinematicValue = &dunendgarmcSamples[iSample].rw_vtx_y[iEvent];
			break;
		case kTrueZPos:
			KinematicValue = &dunendgarmcSamples[iSample].rw_vtx_z[iEvent];
			break;
		case kTrueRad:
			KinematicValue = &dunendgarmcSamples[iSample].rw_rad[iEvent];
			break;
		case kNMuonsRecoOverTruth:
			KinematicValue = &dunendgarmcSamples[iSample].nmuonsratio[iEvent];
			break;
		case kRecoLepEnergy:
			KinematicValue = &dunendgarmcSamples[iSample].rw_elep_reco[iEvent];
			break;
		case kTrueLepEnergy:
			KinematicValue = &dunendgarmcSamples[iSample].rw_elep_true[iEvent];
			break;
		case kRecoXPos:
			KinematicValue = &dunendgarmcSamples[iSample].rw_reco_vtx_x[iEvent];
			break;
		case kRecoYPos:
			KinematicValue = &dunendgarmcSamples[iSample].rw_reco_vtx_y[iEvent];
			break;
		case kRecoZPos:
			KinematicValue = &dunendgarmcSamples[iSample].rw_reco_vtx_z[iEvent];
			break;
		case kRecoRad:
			KinematicValue = &dunendgarmcSamples[iSample].rw_reco_rad[iEvent];
			break;
		case kLepPT:
			KinematicValue = &dunendgarmcSamples[iSample].rw_lep_pT[iEvent];
			break;
		case kLepPZ:
			KinematicValue = &dunendgarmcSamples[iSample].rw_lep_pZ[iEvent];
			break;
		case kTrueQ0:
			KinematicValue = &dunendgarmcSamples[iSample].rw_Q0[iEvent];
			break;
		case kTrueQ3:
			KinematicValue = &dunendgarmcSamples[iSample].rw_Q3[iEvent];
			break;
		default:
			MACH3LOG_ERROR("Did not recognise Kinematic Parameter {}", KinematicParameter);
			throw MaCh3Exception(__FILE__, __LINE__);
	}

	return KinematicValue;
}

const double* samplePDFDUNEBeamNDGAr::GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
	KinematicTypes KinPar = static_cast<KinematicTypes>(std::round(KinematicVariable));
	return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

const double* samplePDFDUNEBeamNDGAr::GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
	KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
	return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

double samplePDFDUNEBeamNDGAr::ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
	return *GetPointerToKinematicParameter(KinematicVariable, iSample, iEvent);
}

double samplePDFDUNEBeamNDGAr::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
	return *GetPointerToKinematicParameter(KinematicParameter, iSample, iEvent);
}

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

