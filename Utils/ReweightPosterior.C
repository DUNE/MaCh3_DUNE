#include <math.h> 

void ReweightPosterior(std::string ReducedChain, int prior = -1){
	if (prior == -1 ){
		std::cout << "Plase select which alternative prior you would like to reweight for: " << std::endl;
		std::cout << "Prior 0: Usual Reactor Constraint" << std::endl;
		std::cout << "Prior 1: Reactor constraint in sin^2 2th13" << std::endl;
		std::cout << "Prior 2: Flat in th13" << std::endl;
		std::cout << "Prior 3: Flat in sin^2 2th13" << std::endl;
		std::cout << "Prior 4: Flat in th23" << std::endl;
		std::cout << "Prior 5: Flat in sin^2 2th23" << std::endl;
		std::cout << "Prior 6: Flat in sindcp" << std::endl;
		std::cout << "Prior 7: Flat in cosdcp" << std::endl;

		throw;
	}
	
	// Load in file
	TFile *fin = new TFile(ReducedChain.c_str(),"READ");
	TTree *treeold, *tree;
	bool reduced; 
	if (fin->Get("osc_posteriors")){
		std::cout << "Loading a reduced chain" << std::endl;
		treeold = (TTree*)fin->Get("osc_posteriors");
		reduced = true; // Using a reduced chain
	}
	else if (fin->Get("posteriors")) {
		// TH: highly recommend reducing chain before reweighting if possible
		std::cout << "Loading a non-reduced chain (posteriors)" << std::endl;
		treeold = (TTree*)fin->Get("posteriors");
		reduced = false; // using a full chain
	}
	else {
		std::cerr << "Make sure you're loading a chain from MaCh3!" << std::endl;
		std::exit(1);
	}
	
	// First get in/out name
	std::size_t pos;
	pos = ReducedChain.find(".root");
	std::string sout = ReducedChain.substr(0,pos) + "_reweighted_sin.root";

	// Make a new file and clone the old TTree
	TFile *fout = new TFile(sout.c_str(),"RECREATE");
	tree = treeold->CloneTree();
	fin->Close();
	delete fin;
	
	// Set brach name to store weights
	std::string branchname, param;
	if (prior == 0) {branchname = "RCreweight";}
	if (prior == 1) {branchname = "RC_sin2_2th13_weights";}
	if (prior == 2) {branchname = "flat_th13_weights";}
	if (prior == 3) {branchname = "flat_sin2_2th13_weights";}
	if (prior == 4) {branchname = "flat_th23_weights";}
	if (prior == 5) {branchname = "flat_sin2_2th23_weights";}
	if (prior == 6) {branchname = "flat_sindcp_weights";}
	if (prior == 7) {branchname = "flat_cosdcp_weights";}

	if (reduced == true){
		if (prior == 0) {param = "theta13";}
		if (prior == 1) {param = "theta13";}
		if (prior == 2) {param = "theta13";}
		if (prior == 3) {param = "theta13";}
		if (prior == 4) {param = "theta23";}
		if (prior == 5) {param = "theta23";}
		if (prior == 6) {param = "dcp";}
		if (prior == 7) {param = "dcp";}
	}
	else if (reduced == false){
		if (prior == 0) {param = "sin2th_13";}
		if (prior == 1) {param = "sin2th_13";}
		if (prior == 2) {param = "sin2th_13";}
		if (prior == 3) {param = "sin2th_13";}
		if (prior == 4) {param = "sin2th_23";}
		if (prior == 5) {param = "sin2th_23";}
		if (prior == 6) {param = "delta_cp";}
		if (prior == 7) {param = "delta_cp";}
	}

	// Variables needed for RC reweight in sin^2 theta_13 (NuFit 4.0)
	double new_central = 0.02241;
	double new_error = 0.00065;
	double old_central = -999;
	double old_error = -999;
	double old_chi, old_prior, new_chi, new_prior;

    // Gaussian in sin^2 2th13 (NuFit4)	
	double new_central_2th13 = 0.088;
	double new_error_2th13 = 0.003;

	// Weights and parameter value
	Double_t reweight;
	Double_t parameter;

	// Make brach to store weights and set parameter of interest
	TBranch *brout = tree->Branch(branchname.c_str(), &reweight);
	tree->SetBranchAddress(param.c_str(), &parameter);
	// turn off branches that aren't needed in the reweighting, helps speed things up
	tree->SetBranchStatus("*",1);
	tree->SetBranchStatus(branchname.c_str(),1);
	tree->SetBranchStatus(param.c_str(), 1);

	Long64_t nEntries = tree->GetEntries();
	std::cout << "Reweighting chain with " << nEntries << " steps" << std::endl;

	for (Long64_t i = 0; i < nEntries; ++i){
		tree->GetEntry(i);

		switch(prior) {
			case 0: {
				// PDG 2020
				new_chi = (parameter - new_central)/new_error;
				new_prior = std::exp(-0.5 * new_chi * new_chi);
				old_prior = 1.0;
				reweight = new_prior/old_prior;
				break;
			}
			case 1: {
				// RC reweight in sin^2 2theta_13 (NuFit4)
				double current = 4. * parameter * (1. - parameter);
				new_chi = (current - new_central_2th13)/new_error_2th13;
				new_prior = std::exp(-0.5 * new_chi * new_chi);
				old_prior = TMath::Abs(4. * (1. - parameter));
				reweight = new_prior/old_prior;
				break;
			}
			case 2:
			case 4: { // flat in theta_ij
				reweight = 1./(2.*sqrt(parameter)*sqrt(1.-parameter));
				break;
			}
			case 3:
			case 5:{ // flat in sin^2 2theta_ij
				reweight = TMath::Abs( 4. * (1. - parameter));
				break;
			}
			case 6:{ // flat in sin(dcp)
				reweight = TMath::Abs((TMath::Cos(parameter)));
				break;
			}
			case 7: { // flat in cos(dcp)
				reweight = TMath::Abs((TMath::Sin(parameter)));
				break;
			}
	
	}
	brout->Fill();
  }
		
	
	tree->Write();
	fout->Close();

}
