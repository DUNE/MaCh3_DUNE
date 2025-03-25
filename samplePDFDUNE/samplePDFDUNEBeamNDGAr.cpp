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
	muonscore_threshold = SampleManager->raw()["SampleCuts"]["muonscore_threshold"].as<float>(); //NK determine what muon score threshold to use
	protondEdxscore = SampleManager->raw()["SampleCuts"]["protondEdxscore_threshold"].as<float>(); //NK determine what proton score threshold to use
	protontofscore = SampleManager->raw()["SampleCuts"]["protontofscore_threshold"].as<float>();  //NK determine what muon score threshold to use
	recovertexradiusthreshold =  SampleManager->raw()["SampleCuts"]["recovertexradius_threshold"].as<float>();  //NK determine what radius threshold to use
	pionenergy_threshold = (SampleManager->raw()["SampleCuts"]["pionenergy_threshold"].as<float>())/1000; //NK determine what muon score threshold to use
	B_field = SampleManager->raw()["SampleCuts"]["B_field"].as<float>(); //NK B field value in T
	momentum_resolution_threshold = SampleManager->raw()["SampleCuts"]["momentum_resolution_threshold"].as<float>(); //NK momentum_resolution threshold, total as a fraction of momentum
	pixel_spacing = SampleManager->raw()["SampleCuts"]["pixel_spacing"].as<float>(); //NK pixel spacing in mm to find num hits in y,z plane
	spatial_resolution = SampleManager->raw()["SampleCuts"]["spatial_resolution"].as<float>(); //NK spatial resolution in mm to find  in y,z plane
	adc_sampling_frequency = SampleManager->raw()["SampleCuts"]["adc_sampling_frequency"].as<float>(); //NK sampling frequency for ADC - needed to find timing resolution and spatial resolution in x dir in MHz
	drift_velocity = SampleManager->raw()["SampleCuts"]["drift_velocity"].as<float>(); //NK drift velocity of electrons in gas - needed to find timing resolution and spatial resolution in x dir in cm/microsecond
	//  average_gain = SampleManager->raw()["SampleCuts"]["average_gain"].as<float>();
	pi0_reco_efficiency = SampleManager->raw()["SampleCuts"]["pi0_reco_efficiency"].as<float>(); //efficiency for pi0 reco in ECAL
	gamma_reco_efficiency = SampleManager->raw()["SampleCuts"]["gamma_reco_efficiency"].as<float>(); //efficiency for gamma reco in ECAL
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
	for (int i = 0; i < (int)dunendgarmcSamples.size(); ++i) {
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

void samplePDFDUNEBeamNDGAr::makePixelGrid(float pixel_spacing_cm){ //make a square pixel grid with spacing defined in yaml file. Spacing must be input in the yaml file in mm, and then is converted to cm later in the code
	int numpixelrows = floor(TPCInstrumentedRadius*2/(pixel_spacing_cm)); //find number of pixels along y and z axis.
	float centre_yboundary, centre_zboundary;
	if(numpixelrows % 2 == 0){centre_yboundary=TPC_centre_y; centre_zboundary=TPC_centre_z;}
	else{centre_yboundary=TPC_centre_y-(pixel_spacing_cm/2); centre_zboundary=TPC_centre_z-(pixel_spacing_cm/2);}
	std::cout<<"num pixel rows: "<<numpixelrows<<std::endl;
	pixelymin = centre_yboundary - floor(numpixelrows/2)*pixel_spacing_cm;
	pixelymax = centre_yboundary + floor(numpixelrows/2)*pixel_spacing_cm;
	pixelzmin = centre_zboundary - floor(numpixelrows/2)*pixel_spacing_cm;
	pixelzmax = centre_zboundary + floor(numpixelrows/2)*pixel_spacing_cm;

	for(int i_pixel = 0; i_pixel<numpixelrows; i_pixel++){
		yboundarypositions.push_back(pixelymin+i_pixel*pixel_spacing_cm);
		zboundarypositions.push_back(pixelzmin+i_pixel*pixel_spacing_cm);
	}
}

double samplePDFDUNEBeamNDGAr::FindNHits(float pixel_spacing_cm, float centre_circle_y, float centre_circle_z, double rad_curvature){
	//use the pixel grid method to find number of pixels hit in a track.

	int num_vertices =0; // number of vertices hit. Counting to avoid duplicating
	int num_intersections =0;

	//equation for circle = (y-y0)^2 + (z-z0)^2 = r^2

	for(int i_intersect = 1; i_intersect<=yboundarypositions.size(); i_intersect++){ //check every boundary line to see if it crossed within the TPC instrumented region
		float quadratic_ineq_y= pow(rad_curvature*100, 2)-pow((yboundarypositions[i_intersect]-centre_circle_y), 2);
		if(quadratic_ineq_y > 0){
			if((pow((pow((centre_circle_z + pow(quadratic_ineq_y, 0.5)-TPC_centre_z), 2) + pow((yboundarypositions[i_intersect]-TPC_centre_y), 2)), 0.5) <=TPCInstrumentedRadius) && (pixelzmin<(centre_circle_z + pow(quadratic_ineq_y, 0.5))<=pixelzmax)){ //check that the z coord is also on pixel plane
				num_intersections++;
				if((double)(fmod((centre_circle_z + (float)(pow(quadratic_ineq_y, 0.5)) - pixelzmin), pixel_spacing_cm)) == 0){ // this is the case when a vertex is crossed so to avoid double counting pixels
					num_vertices++;
				}
			}
			if((pow((pow((centre_circle_z - pow(quadratic_ineq_y, 0.5)-TPC_centre_z), 2) + pow((yboundarypositions[i_intersect]-TPC_centre_y), 2)), 0.5) <TPCInstrumentedRadius) && pixelzmin<(centre_circle_z - pow(quadratic_ineq_y, 0.5))<=pixelzmax){ //check that the z coord is also on pixel plane
				num_intersections++;
				if((double)(fmod((centre_circle_z - (float)(pow(quadratic_ineq_y, 0.5)) - pixelzmin), pixel_spacing_cm)) == 0){ // this is the case when a vertex is crossed so to avoid double counting pixels
					num_vertices++;
				}
			}
		}
		else if(quadratic_ineq_y == 0){ //when the pixel boundary is a tangent to the circle z = z0
			num_intersections++;
		}
	}
	for(int i_intersect = 1; i_intersect<=zboundarypositions.size(); i_intersect++){
		float quadratic_ineq_z = pow(rad_curvature*100, 2)-pow((zboundarypositions[i_intersect]-centre_circle_z), 2);
		if(quadratic_ineq_z > 0){
			if( (pow((pow((centre_circle_y + pow(quadratic_ineq_z, 0.5)-TPC_centre_y), 2) + pow((zboundarypositions[i_intersect]-TPC_centre_z), 2)), 0.5) <=TPCInstrumentedRadius) && pixelymin<(centre_circle_y + pow(quadratic_ineq_z, 0.5))<=pixelymax){ //check that the z coord is also on pixel plane
				num_intersections++;
			}
			if((pow((pow((centre_circle_y - pow(quadratic_ineq_z, 0.5)-TPC_centre_y), 2) + pow((zboundarypositions[i_intersect]-TPC_centre_z), 2)), 0.5) <=TPCInstrumentedRadius) && pixelymin<(centre_circle_y - pow(quadratic_ineq_z, 0.5))<=pixelymax){ //check that the z coord is also on pixel plane
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
	bg = (double)(p_mag/pdgmass); //beta*gamma
	gamma = pow((1+bg*bg), 0.5); //gamma
	double beta = bg/gamma; //beta (velocity)
	return beta;
}

double samplePDFDUNEBeamNDGAr::CalcDeDx(double beta, double bg, double gamma){ //calc de/dx in 10bar argon
	double mer = m_e/pdgmass; //electron mass/mass of particle
	double tmax = 2.*1000*m_e* bg*bg / (1. + 2.*gamma*mer + mer*mer);  // Maximum delta ray energy (MeV).
	//Assume tcut = tmax
	double tcut = tmax;
	//Find density effect correction (delta). Sternheimer values set in header file
	double log_bg = std::log10(bg);
	double delta = 0;
	if(log_bg >= sternheimer_X0){
		delta = 2. * std::log(10.)*log_bg - sternheimer_Cbar;
		if(log_bg < sternheimer_X1){
			delta += sternheimer_A * std::pow(sternheimer_X1 - log_bg, sternheimer_K);
		}
	}

	//Calculate Stopping number, B
	double B = 0.5 * std::log(2.*1000*m_e*bg*bg*tcut / (1.e-12 * pow(excitationenergy, 0.5)))- 0.5*beta*beta*(1. + tcut / tmax) - 0.5*delta;
	if(B<1.){B=1.;} //Don't let B become negative

	//Calculate dE/dX 
	double dedx = density*K_const*18*B / (39.981 * beta*beta); //18 is atomic number and 39.981 is atomic mass in g/mol

	return dedx;
}

bool samplePDFDUNEBaseNDGAr::IsParticleAccepted(dunendgarmc_base *duneobj, int& i_truepart, int& i, double& highestpT, float pixel_spacing_cm, int& tot_particles){
	
	bool isAccepted = true;
	int nummatched = 0;
	
	for(int i_anapart =0; i_anapart<_MCPStartPX->size(); i_anapart++){

		float mom_tot = pow(pow(_MCPStartPX->at(i_anapart), 2)+pow(_MCPStartPY->at(i_anapart), 2)+pow(_MCPStartPZ->at(i_anapart), 2), 0.5);
		//Now select particle in the anatree that is the same as in the caf (same pdg and momentum)
		if((_PDG->at(i_anapart) == sr->mc.nu[0].prim[i_truepart].pdg) && ((double)(mom_tot) >= 0.999*(double)(sr->mc.nu[0].prim[i_truepart].p.Mag())) && ((double)(mom_tot) <= 1.001*(double)(sr->mc.nu[0].prim[i_truepart].p.Mag()))){			
			nummatched++;
			//Ignore neutrons and neutrinos
			if((std::abs(sr->mc.nu[0].prim[i_truepart].pdg))!= 2112 && (std::abs(sr->mc.nu[0].prim[i_truepart].pdg))!= 14 && (std::abs(sr->mc.nu[0].prim[i_truepart].pdg))!= 12){
				bool stopsinecal_radius = false;
				bool stopsinecal_length = false;
				float end_radius = pow((pow(_MCPEndY->at(i_anapart)-TPC_centre_y, 2)+pow(_MCPEndZ->at(i_anapart)-TPC_centre_z, 2)), 0.5);
				float end_length = _MCPEndX->at(i_anapart)-TPC_centre_x;
				if(ecal_containment && end_radius>ECALInnerRadius && end_radius<ECALOuterRadius){stopsinecal_radius = true;}
				if(ecal_containment && std::abs(end_length)>ECALEndCapStart && std::abs(end_length)<ECALEndCapEnd){stopsinecal_length = true;}

				//Momentum transverse to B field
				double transverse_mom = pow(pow(_MCPStartPY->at(i_anapart), 2)+pow(_MCPStartPZ->at(i_anapart), 2), 0.5);
				float start_radius =pow((pow(_MCPStartY->at(i_anapart)-TPC_centre_y, 2)+pow(_MCPStartZ->at(i_anapart)-TPC_centre_z, 2)), 0.5);

				//Fill particle-level kinematic variables
				duneobj->particlepdg[tot_particles]=_PDG->at(i_anapart);
				duneobj->particleenergy[tot_particles]=(double)(sr->mc.nu[0].prim[i_truepart].p.E);
				duneobj->particlemomentum[tot_particles]=(double)(mom_tot);
				//particleBAngle defined as angle wrt B-field (x)
				duneobj->particleBAngle[tot_particles]=(double)(acos(_MCPStartPX->at(i_anapart)/mom_tot)*180/M_PI);

				//If particle is not stopped in the tpc or ecal 
				if((std::abs(end_length)>TPCInstrumentedLength && !stopsinecal_length) || (end_radius>TPCInstrumentedRadius && !stopsinecal_radius)){
					//If stopped between tpc and ecal, determine where
					if(end_radius<ECALOuterRadius && std::abs(end_length)<ECALEndCapEnd) {
						duneobj->isparticlestoppedingap[tot_particles] = true;
						if (std::abs(end_length)>TPCInstrumentedLength) {duneobj->isparticlestoppedinendgap[tot_particles] = true;}
						else {duneobj->isparticlestoppedinbarrelgap[tot_particles] = true;}
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
							double helixlength = ((std::abs(length_track_x)/100)/pitch)*pow((pow(M_PI*2*rad_curvature, 2) + pow(pitch, 2)), 0.5); //L = height/pitch*sqrt((pi*diameter)**2 + pitch**2)
							double tan_theta = tan(theta_xT);

							//Find centre of circular path
							bool positivecharged =0;
							float centre_circle_y;
							float centre_circle_z;
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
							float m_const = (TPC_centre_z - centre_circle_z)/(TPC_centre_y-centre_circle_y); //gradient of line between two intersection points
							float a_const = (pow(TPCInstrumentedRadius, 2)-pow(rad_curvature*100, 2) - (pow(TPC_centre_y, 2)-pow(centre_circle_y, 2))-(pow(TPC_centre_z, 2)-pow(centre_circle_z, 2)))/(2*(centre_circle_y-TPC_centre_y));
							float quadraticformula_b = -(2*m_const*(a_const -TPC_centre_y)+2*TPC_centre_z);
							float quadraticformula_a = pow(m_const, 2)+1;
							float quadraticformula_c = pow((a_const - TPC_centre_y), 2) +pow(TPC_centre_z, 2) - pow(TPCInstrumentedRadius,2);

							double z_intersect_1, y_intersect_1, z_intersect_2, y_intersect_2, z_intersect_chosen, y_intersect_chosen, theta_1, theta_2, theta_start, theta_chosen, theta_diff_1, theta_diff_2;
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
									y_intersect_chosen = y_intersect_1; 
									z_intersect_chosen = z_intersect_1;
								} 
								else {
									theta_chosen = theta_diff_2; 
									y_intersect_chosen = y_intersect_2; 
									z_intersect_chosen = z_intersect_2;
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
								theta_chosen = theta_intersect;
							}

							double N_hits = FindNHits(pixel_spacing_cm, centre_circle_y, centre_circle_z, rad_curvature);

							double p_mag = sr->mc.nu[0].prim[i_truepart].p.Mag();
							double bg = 0; double gamma = 0;
							double beta = CalcBeta(p_mag, bg, gamma);
							double dedx = CalcDeDx(beta, bg, gamma);

							double sigmax = (drift_velocity/100)*(1/(adc_sampling_frequency));
							double sigmayz = (spatial_resolution/(1000)); //needs to be in m              
							double momres_yz = transverse_mom*(pow(720/(N_hits+4), 0.5)*(sigmayz*transverse_mom/(0.3*B_field*pow(L_yz_chord/100, 2)))*pow((1-(1/21)*pow((L_yz_chord/(rad_curvature*100)), 2)), 0.5));
							double momres_ms = transverse_mom*(0.016/(0.3*B_field*(L_yz/100)*cos(theta_xT)*beta))*pow(L_yz/X0, 0.5);
							double momres_tottransverse = pow(pow(momres_yz, 2) + pow(momres_ms, 2), 0.5);
							double momres_yz_old = transverse_mom*(pow(720/(N_hits+4), 0.5)*(sigmayz*transverse_mom/(0.3*B_field*pow(L_yz/100, 2))));
							double sigma_theta = (pow(cos(theta_xT), 2))*(pitch/(2*M_PI*rad_curvature))*pow((pow(sigmax/(std::abs(length_track_x)/100),2) +pow(momres_tottransverse/transverse_mom, 2)), 0.5);
							double momres_frac = pow(pow((momres_tottransverse/transverse_mom), 2)+pow(sigma_theta*tan_theta, 2), 0.5);

							duneobj->particlemomresms[tot_particles] = momres_ms;
							duneobj->particlemomrestransfrac[tot_particles] = (momres_tottransverse/transverse_mom)/momres_frac;
							duneobj->particlemomrestrans[tot_particles] = momres_tottransverse/transverse_mom;

							if(momres_frac > momentum_resolution_threshold){
								isAccepted = false;
							}
						}
						else{ //Accept pi0 and gamma with some efficiency
							if((std::abs(sr->mc.nu[0].prim[i_truepart].pdg)) == 111){
								double pi0recoprob =(double)(std::rand())/RAND_MAX;
								if(pi0recoprob>pi0_reco_efficiency){
									isAccepted = false;
								}
							}
							else if((std::abs(sr->mc.nu[0].prim[i_truepart].pdg)) == 22){
								double gammarecoprob =(double)(std::rand())/RAND_MAX;
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
					for(int i_ecaldep=0; i_ecaldep<_SimHitEnergy->size(); i_ecaldep++){
						if(_SimHitTrkID->at(i_ecaldep) == _MCPTrkID->at(i_anapart)){
							energydepsum = energydepsum + _SimHitEnergy->at(i_ecaldep);
						}
					}
					duneobj->particleecaldepositfraction[tot_particles] = energydepsum/(sr->mc.nu[0].prim[i_truepart].p.E);
					duneobj->isparticlestoppedinecal[tot_particles] = true;
					if(end_radius>TPCInstrumentedRadius) {duneobj->isparticlestoppedinbarrel[tot_particles] = true;}
					else {duneobj->isparticlestoppedinendcap[tot_particles] = true;}
				}
				else {duneobj->isparticlestoppedintpc[tot_particles] = true;}
				
				if(isAccepted) {
					//fill particle-level accepted kinematic parameters
					duneobj->accparticlemomentum[tot_particles]=(double)(mom_tot);
					duneobj->accparticleBAngle[tot_particles]=(double)(acos(_MCPStartPX->at(i_anapart)/mom_tot)*180/M_PI);
					duneobj->accparticlepdg[tot_particles]=_PDG->at(i_anapart);
				}
				else {
					//fill particle-level rejected kinematic parameters
					duneobj->rejparticlex[tot_particles] = _MCPStartX->at(i_anapart)-TPC_centre_x;
					duneobj->rejparticler2[tot_particles] = start_radius*start_radius;
				}	
			}
			//break;
		}
	}
	if(nummatched != 1){MACH3LOG_INFO("Found {} matching particles in anatree");}
	return isAccepted;
}

int samplePDFDUNEBeamNDGAr::setupExperimentMC(int iSample) {

	dunemc_base *duneobj = &(dunendgarmcSamples[iSample]);
	//int nutype = sample_nutype[iSample];
	//int oscnutype = sample_oscnutype[iSample];
	//bool signal = sample_signal[iSample];

	MACH3LOG_INFO("-------------------------------------------------------------------");
	MACH3LOG_INFO("Input File: {}", mc_files.at(iSample));

	_sampleFile = TFile::Open(mc_files.at(iSample).c_str(), "READ");
	_data = (TTree*)_sampleFile->Get("cafTree");

	if(_data){
		MACH3LOG_INFO("Found \"caf\" tree in {}", mc_files[iSample]);
		MACH3LOG_INFO("With number of entries: {}", _data->GetEntries());
	}
	else{
		MACH3LOG_ERROR("Could not find \"caf\" tree in {}", mc_files[iSample]);
		throw MaCh3Exception(__FILE__, __LINE__);
	}

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
		_data_geant = (TTree*)_sampleFile_geant->Get("GArAnaTree");

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

	duneobj->norm_s = 1.;
	duneobj->pot_s = pot/1e21;
	duneobj->osc_channel = (double)iSample;
	duneobj->nEvents = _data->GetEntries();
	//duneobj->nutype = nutype;
	//duneobj->oscnutype = oscnutype;
	//duneobj->signal = signal;

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
	duneobj->rw_lep_pZ = new double[duneobj->nEvents];

	duneobj->rw_reco_vtx_x = new double[duneobj->nEvents];
	duneobj->rw_reco_vtx_y = new double[duneobj->nEvents];
	duneobj->rw_reco_vtx_z = new double[duneobj->nEvents];
	duneobj->rw_reco_rad = new double[duneobj->nEvents];

	duneobj->Target = new int[duneobj->nEvents];

	int num_no_ixns =0;
	int num_no_recparticles = 0;
	int num_in_fdv = 0;
	int num_in_fdv_noreco = 0;
	int num_notin_fdv =0;
	int num_nanenergy =0;
	int num_nanparticles =0;

	//FILL DUNE STRUCT
	for (int i = 0; i < (duneobj->nEvents); ++i) { // Loop through tree
		_data->GetEntry(i);
		double radius = pow((pow((sr->mc.nu[0].vtx.y+150),2) + pow((sr->mc.nu[0].vtx.z-1486),2)),0.5);
		if(std::abs(sr->mc.nu[0].vtx.x)<=209.0 &&  radius<=227.02){
			num_in_fdv++;
			duneobj->in_fdv[i] = 1;
		} else{
			num_notin_fdv++;
			duneobj->in_fdv[i] = 0;
		}

		if(sr->common.ixn.ngsft == 0){
			//duneobj->rw_erec[i] = (double)(0);
			float erec_total =0;
			duneobj->rw_elep_reco[i] = double(0);
			duneobj->rw_yrec[i] = (double)(0);
			num_no_ixns++;
			duneobj->nrecoparticles[i] = (int)(0);
		} else{
			duneobj->nrecoparticles[i] = (int)(0);
			float erec_total =0;
			float elep_reco =0;
			float muonscore = muonscore_threshold;
			int nixns = (int)(sr->common.ixn.ngsft);
			for(int i_ixn =0; i_ixn<nixns; i_ixn++) {
				int nrecoparticles = (int)(sr->common.ixn.gsft[i_ixn].part.ngsft);
				duneobj->nrecoparticles[i] += (int)(sr->common.ixn.gsft[i_ixn].part.ngsft);
				int nanparticles = 0;
				if(nrecoparticles ==0){
					double radius = pow((pow((sr->mc.nu[0].vtx.y+150.),2) + pow((sr->mc.nu[0].vtx.z-1486.),2)),0.5);
					if(std::abs(sr->mc.nu[0].vtx.x)<=209.0 || radius<=227.02) {
						num_in_fdv_noreco++;
					}   
					num_no_recparticles++;}
				for(int i_part =0; i_part<nrecoparticles; i_part++) {
					float erec_part = (float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].E);
					if(std::isnan(erec_part)){nanparticles++;}
					erec_total+=erec_part;
					if((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].score.gsft_pid.muon_score>muonscore)){
						if(erec_part>elep_reco){
							elep_reco = erec_part;
							duneobj->rw_reco_vtx_x[i] = (double)((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.x));
							duneobj->rw_reco_vtx_y[i] = (double)((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.y));
							duneobj->rw_reco_vtx_z[i] = (double)((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].start.z));
							duneobj->rw_lep_pT[i] = (double)(pow(pow((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.x), 2) + pow((float)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.y), 2), 0.5));
							duneobj->rw_lep_pZ[i] = (double)(sr->common.ixn.gsft[i_ixn].part.gsft[i_part].p.z);
						}
						duneobj->nrecomuon[i]++; 
					}
				}
				num_nanparticles = num_nanparticles + (nanparticles/nrecoparticles);
			} //ADD PRIMARY LEPTON ENERGY ELEP_RECO
			if(std::isnan(erec_total)){num_nanenergy++; erec_total = (float)(sr->common.ixn.gsft[0].Enu.lep_calo);}
			if(iscalo_reco){duneobj->rw_erec[i]=(double)(sr->common.ixn.gsft[0].Enu.lep_calo);}
			else{duneobj->rw_erec[i]=(double)(erec_total);}
			duneobj->rw_elep_reco[i] = (double)(elep_reco);
		}

		if(duneobj->rw_erec[i] != 0){duneobj->rw_yrec[i] = (double)(((duneobj->rw_erec[i])-(duneobj->rw_elep_reco[i]))/(duneobj->rw_erec[i]));}
		else{duneobj->rw_yrec[i] = (double)(0);}
		duneobj->rw_etru[i] = (double)(sr->mc.nu[0].E); // in GeV
		duneobj->rw_isCC[i] = (int)(sr->mc.nu[0].iscc);
		duneobj->rw_nuPDGunosc[i] = sr->mc.nu[0].pdgorig;
		duneobj->rw_nuPDG[i] = sr->mc.nu[0].pdg;
		duneobj->rw_berpaacvwgt[i] = _BeRPA_cvwgt;

		int ntrueparticles = (int)(sr->mc.nu[0].nprim);
		for(int i_truepart =0; i_truepart<ntrueparticles; i_truepart++){
			if(std::abs(sr->mc.nu[0].prim[i_truepart].pdg) == 13){
				duneobj->ntruemuon[i]++;
				duneobj->ntruemuonprim[i]++;
			}
		}
		int ntruesecparticles = (int)(sr->mc.nu[0].nsec);
		for(int i_truepart =0; i_truepart<ntruesecparticles; i_truepart++){
			if(std::abs(sr->mc.nu[0].sec[i_truepart].pdg) == 13){
				duneobj->ntruemuon[i]++;
			}
		}

		duneobj->nproton[i] = sr->mc.nu[0].nproton;
		duneobj->nneutron[i] = sr->mc.nu[0].nneutron;
		duneobj->npip[i] = sr->mc.nu[0].npip;
		duneobj->npim[i] = sr->mc.nu[0].npim;
		duneobj->npi0[i] = sr->mc.nu[0].npi0;

		duneobj->nmuonsratio[i] = (double)(duneobj->nrecomuon[i])/(double)(duneobj->ntruemuonprim[i]);
		duneobj->rw_vtx_x[i] = (double)(sr->mc.nu[0].vtx.x);
		duneobj->rw_vtx_y[i] = (double)(sr->mc.nu[0].vtx.y);
		duneobj->rw_vtx_z[i] = (double)(sr->mc.nu[0].vtx.z);

		duneobj->rw_rad[i] = (double)(pow((pow((duneobj->rw_vtx_y[i]+150),2) + pow((duneobj->rw_vtx_z[i]-1486),2)),0.5)); 
		duneobj->rw_reco_rad[i] = (double)(pow(pow((duneobj->rw_reco_vtx_y[i]+150),2) + pow((duneobj->rw_reco_vtx_z[i]-1486), 2), 0.5));
		duneobj->rw_elep_true[i] = (double)(sr->mc.nu[0].prim[0].p.E);

		//Assume everything is on Argon for now....
		duneobj->Target[i] = 40;

		duneobj->rw_Q0[i] = (double)(sr->mc.nu[0].q0);
		duneobj->rw_Q3[i] = (double)(sr->mc.nu[0].modq);

		_mode = sr->mc.nu[0].mode;
		_isCC = (int)(sr->mc.nu[0].iscc);

		int M3Mode = Modes->GetModeFromGenerator(std::abs(sr->mc.nu[0].mode));
		duneobj->mode[i] = M3Mode;
		//int mode= TMath::Abs(_mode);       
		//duneobj->mode[i]=(double)GENIEMode_ToMaCh3Mode(mode, _isCC);

		duneobj->flux_w[i] = 1.0;
	}

	_sampleFile->Close();
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
		case kOscChannel:
			KinematicValue = &dunendgarmcSamples[iSample].osc_channel;	
			break;
		case kTrueQ0:
			KinematicValue = &dunendgarmcSamples[iSample].rw_Q0[iEvent];
			break;
		case kTrueQ3:
			KinematicValue = &dunendgarmcSamples[iSample].rw_Q3[iEvent];
			break;
		default:
			std::cout << KinematicParameter << std::endl;
			MACH3LOG_ERROR("Did not recognise Kinematic Parameter type...");
			throw MaCh3Exception(__FILE__, __LINE__);
	}

	return KinematicValue;
}

const double* samplePDFDUNEBeamNDGAr::GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
	KinematicTypes KinPar = (KinematicTypes) std::round(KinematicVariable);
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
		case kPionMultiplicity:
			for(double ibins =0; ibins<10; ibins++){
				binningVector.push_back(ibins);
			}
			break;
		case kNRecoParticles:
			for(double ibins =0; ibins<50; ibins++){
				binningVector.push_back(ibins);
			}
			break; 
		case kInFDV:
			for(double ibins =0; ibins<3; ibins++){
				binningVector.push_back(ibins);
			}
			break;
		case kNMuonsRecoOverTruth:
		case kTrueMinusRecoEnergyRatio:
			for(double ibins =0; ibins<20*10; ibins++){
				binningVector.push_back(-10+(double)(ibins)/10);
			}
			break;
		case kTrueMinusRecoEnergy:
			for(double ibins =0; ibins<20*10; ibins++){
				binningVector.push_back(-10+(double)(ibins)/10);
			}
			break;
		case kNTrueMuons:
		case kNRecoMuons:
			for(double ibins =0; ibins<10; ibins++){
				binningVector.push_back(ibins);
			}
			break;
		case kRecoLepEnergy:
			for(double ibins =0; ibins<10*10; ibins++){
				binningVector.push_back((double)(ibins)/10);
			} 
			break;
		case kTrueLepEnergy:
			for(double ibins =0; ibins<10*10; ibins++){
				binningVector.push_back((double)(ibins)/10);
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
				binningVector.push_back((double)(ibins)/10);
			}
			break;
		case kOscChannel:
			for (int ibins=0; ibins<GetNsamples(); ibins++){
				binningVector.push_back(ibins);
			}
		case kTrueQ0:
		case kTrueQ3:
			for(double ibins =0; ibins<3*50+1; ibins++){
				double binval = ibins/50;
				binningVector.push_back(binval);
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

