#include <iostream>
#include <vector>
/*#---------------------------------- Some common variables ----------------------------------------#*/
//YSP: The sample names are slightly different between OA2023 and OA2024
static const int nSamples = 6; 

TString samplenames[nSamples] = {"numu","numubar","numucc1pi","nue","nuebar","nue1pi"}; //this is the generic year independent name of the samples

TString *samplenameuser;

TFile* outfile = NULL;

TString samplenameOA2023[nSamples] = {"FHC1Rmu_2023","RHC1Rmu_2023","FHCnumuCC1pi_2023","FHC1Re_2023","RHC1Re_2023","FHC1Re1de_2023"};
TString samplenameOA2024[nSamples] = {"FHC1Rmu_2024","RHC1Rmu_2024","FHCnumuCC1pi_2024","FHC1Re_2024","RHC1Re_2024","FHCnueCC1pi_2024"};

//YSP: This sets the max y axis range for the posterior distributions. WIP to remove this altogether
int maxysample[nSamples] = {150,150,150,25,25,5};

//Added this array to reduce all the if statements and string comparisons look cleaner now
TString M3Mode[16] = {"CCQE","CC1pipm","CCcoh","CCMpi","CCDIS","NC1pi0","NC1pipm","NCcoh","NCoth","2p2h","NC1gam","CCMisc","NCMpi","NCDIS","CC1pi0","UNKNOWN_BAD"};
TString SplineMode[11] = {"ccqe","cccoh","ccmpi","ccdis","nc1pi0","nc1pipm","nccoh","ncoth","cc2p2h","nc1gam","ccmisc"};

/*#---------------------------------- Some common variables ----------------------------------------#*/



//NOTE use run_makeSKspectra.C to process many files
void makeSKspectra_sample(const char *infile = "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/Chain_testing/newconfig_eventrate.root", TString name ="ND_FHC_CCnumu", TString outname = "output.root", int maxy = 1, int nthrows = 1, bool do_var = false, TString var =""){
  TFile* file = new TFile(infile, "open");
  // ETA -Set the number of z-axis bins to be something not crazy
  //this makes the plots easier on the eyes
  std::cout << "infile is " << infile << std::endl;
  int nbins = maxy*10;
  //std::cout << "Setting nbinsy to be " << nbins << std::endl;
  //std::cout << "Max bins is " << maxy << std::endl;
  TString tempname = name+"_1";
  std::cout << "Got temp name as " << tempname << std::endl;
  std::cout << "Got name as " << name << std::endl;
  TH1D* temp = (TH1D*)file->Get(tempname);
  TH1D* nominal = (TH1D*)file->Get(name);

  //Check to see if reweight tree exists
  TTree *reweight = 0;
  double weight = 1;
  file->GetObject("reweight", reweight);

  TFile* outfile = new TFile("output.root", "RECREATE");


  // ETA - as of OA2020 we no longer run chains for both woRC
  // and wRC, instead we reweight the woRC. To make posterior predictions
  // for wRC we need to fill the TH2D with appropriate weights
  // depending on the th13 value. If the reweight tree exists in the file
  // then weights will be applied. You add the reweight tree to the
  // reweight_prior.C script
  //
  if(reweight){
		std::cout << "FOUND REWEIGHT TREE" << std::endl;
		std::cout << "I WILL USE THE WEIGHTS WHEN FILLING THE SPECTRA!" << std::endl;
		reweight->SetBranchAddress("weight", &weight);
  }
  else{
  	std::cout << "DIDN'T FIND REWEIGHT TREE" << std::endl;
  }

  int nbinsx = temp->GetNbinsX();
  //double binedges[100];

  std::cout << "BLAH" << std::endl;
 
  /*
  std::vector<double> binedges(nbinsx + 1);

  for (int i = 0; i < nbinsx; i++) {
	binedges[i] = temp->GetXaxis()->GetBinLowEdge(i + 1);
  }
  binedges[nbinsx] = temp->GetXaxis()->GetBinUpEdge(nbinsx);
  */

  std::vector<double> binedges(nbinsx + 1);
  for (int i = 0; i <= nbinsx; i++) {
	  binedges[i] = temp->GetXaxis()->GetBinLowEdge(i + 1);
	}

  std::cout << "Made Th2D spectra" << std::endl;
  TH2D* spectra = new TH2D("test","test",nbinsx,binedges.data(),nbins+1,0,maxy); //400,0,/*20*/24/*100*/); //500,0,10
  spectra->Sumw2(true);
  std::cout << "Made Th2D spectra" << std::endl;
  printf("%i \n",nbinsx);
  double average=0;

  // I don't remember what these histograms were for or why they were useful,
  // and I never use them anyway, so commented out
  // ETA - they're useful for error on the event rate!!!!!
  TH1D* hAvg = new TH1D("hAvg","hAvg",1000,0,500);
  hAvg->Sumw2(true);

  int nsteps = nthrows;
  std::cout << "Number of throws is " << nthrows << std::endl;

  TRandom3* rnd = new TRandom3(0);
  for(int i = 0; i < nsteps; i++){
	std::cout << "Getting histogram for throw " << i << std::endl;
    if(reweight){
		std::cout << "In reweight" << i << std::endl;
        reweight->GetEntry(i);
	    std::cout << "Done reweight" << i << std::endl;
	 	int print_step = std::max(1, nsteps / 10);  // avoid division by 0

		if (i % print_step == 0) {
			std::cout << "Weight is " << weight << std::endl;
		}
		std::cout << "nsteps = " << nsteps << ", print_step = " << print_step << std::endl;
	}

	TString name1 = "";
	name1.Form("%s_%i", name.Data(), i);
	std::cout << "Looking for histogram with name " << name1 << std::endl;
	temp = (TH1D*)file->Get(name1);
	average+=temp->Integral();
	hAvg->Fill(temp->Integral(), weight);
	if (!temp || !spectra) {
    std::cerr << "Either 'temp' or 'spectra' is null!" << std::endl;
    return;
	}
	if (spectra->GetDimension() != 2) {
		std::cerr << "Spectra is not a 2D histogram!" << std::endl;
		return;
	}

	for (int j = 1; j <= temp->GetNbinsX(); ++j) {
		double x_val = temp->GetXaxis()->GetBinCenter(j);
		double y_val = temp->GetBinContent(j);
		spectra->Fill(x_val, y_val, weight);
	}
	std::cout << "Filled spectra " << std::endl;
	
	average/=(double)nsteps;
	printf("%f %f %f\n",average,hAvg->GetMean(),hAvg->GetRMS());
	//printf("%f \n",average);

	TCanvas* c = new TCanvas("c","c",0,0,700,700);
	spectra->Draw("colz");
	c->Update();
  	outfile->cd();

	////////////////////////////////////////
	//
	//  Now make 1D projections
	//  Nice to look at these to see how
	//  funcky the prediction is.
	//  (Hopefully it's nice and guassian :p)
	//
	////////////////////////////////////////

  std::cout << "Now looping through to make projections" << std::endl;
  double x[100], y[100], ex[100], ey[100];
  c->Clear();
  TH1D* proj=0;

  for (int i = 1; i <= nbinsx; i++) {
	std::cout << "on bin " << i << std::endl;

	if (!spectra) {
		std::cerr << "Error: 'spectra' histogram is null!" << std::endl;
		break;
	}
	if (i < 1 || i > spectra->GetNbinsX()) {
		std::cerr << "Invalid bin index i = " << i << std::endl;
		continue;
	}

	TString namepy = name + "_py" + i;
	TH1D* proj = spectra->ProjectionY(namepy, i, i);
	if (!proj) {
		std::cerr << "Projection failed on bin " << i << std::endl;
		continue;
	}

	if (proj->GetEntries() > 0) {
		int first_bin = proj->GetNbinsX() + 1;
		int last_bin = 0;

		for (int bin_i = 1; bin_i <= proj->GetNbinsX(); ++bin_i) {
			double val = proj->GetBinContent(bin_i);
			if (val > 0) {
				if (bin_i < first_bin) first_bin = bin_i;
				if (bin_i > last_bin) last_bin = bin_i;
			}
		}

		if (first_bin > last_bin) continue;

		double diff = proj->GetBinCenter(last_bin) - proj->GetBinCenter(first_bin);
		double opt_bin_width = diff / 10.0;
		double current_bin_width = proj->GetBinWidth(last_bin);
		int rebin = 1;

		if (diff > 0) {
			rebin = int(opt_bin_width / current_bin_width);
			if (rebin == 0) rebin = 1;
			while ((proj->GetNbinsX() % rebin) != 0 && rebin > 1) {
				rebin--;
			}
		}
		if (rebin > 1) proj->Rebin(rebin);

		// Set the graph arrays
		x[i-1]  = spectra->GetXaxis()->GetBinCenter(i);
		ex[i-1] = spectra->GetXaxis()->GetBinWidth(i) / 2.0;
		y[i-1]  = proj->GetMean();
		ey[i-1] = proj->GetRMS();
	}
	proj->Write(namepy);
}


  //////////////////////////////////
  //
  //  Time to write things to file!!
  //
  /////////////////////////////////

  TGraphErrors* errorbars = new TGraphErrors(nbinsx,x,y,ex,ey);
  errorbars->Draw("A E2");

  TString dummyname="";

  for(int i=0; i<nSamples; ++i){
  	if(do_var){
  		if(name==samplenameuser[i]+var){dummyname = samplenames[i];}
  	}
  	else{
  		if(name==samplenameuser[i]){dummyname = samplenames[i];}
	}
  }

  // make TH1D with y-errors only
  TString histname = "hist_1d_"+dummyname;
  TH1D *hist1d = (TH1D*)temp->Clone(histname);
  hist1d->Reset();
  for (int i=0; i<hist1d->GetXaxis()->GetNbins(); i++)
  {
		hist1d->SetBinContent(i+1,y[i]);
		hist1d->SetBinError(i+1,ey[i]);
  }

  errorbars->Write(name);
  histname="spectra2d_"+dummyname;
  spectra->Write(histname);
  histname="hist_1d_"+dummyname;
  hist1d->Write(histname);
  histname=name+"_nominal";
  nominal->Write(histname);
  histname="hist_1d_"+dummyname+"_events";
  hAvg->Write(histname);

  file->Close();
  //delete file;
  //delete outfile;
  //delete spectra;
  //delete errorbars;
  //delete rnd;
  //delete c;
}
// -------------------------------------------------------------------------------- //
// Same function as above really but has to deal with MaCh3 modes
// this is made by setting do_by_mode=true in posteriorPredictiveSK2020
void makeSKspectra_sample_by_mode(const char *infile = "sk_nominal.root", TString name = "hist", TString outname = "output.root", int maxy = 1, int nthrows = 1){

  TFile* file = new TFile(infile, "open");
  std::cout << "infile is " << infile << std::endl;
  int nbins = maxy*1000;
  //std::cout << "Setting nbinsy to be " << nbins << std::endl;
  //std::cout << "Max bins is " << maxy << std::endl;

  TString tempname = name+"_ccqe0";
  std::cout << "Got temp name as " << tempname << std::endl;
  std::cout << "Got name as " << name << std::endl;
  TH1D* temp = (TH1D*)file->Get(tempname);
  //TH1D* nominal = (TH1D*)file->Get(name);

  int nbinsx = temp->GetNbinsX();
  double binedges[100];

  std::cout << "BLAH" << std::endl;
  printf("%i \n",nbinsx);

  for(int i=0; i<=nbinsx; i++){
  	binedges[i]=temp->GetXaxis()->GetBinLowEdge(i+1);
  }

  std::cout << "Made Th2D spectra" << std::endl;
  std::vector<TH2D*> spectra;
  std::vector<TH1D*> hAvg;
  TString mode = "";
  //Find the correct mode string
  //Reduced hardcoding, see details in Structs.h
  for(int mode_i = 0 ; mode_i < 11 ; mode_i++){
  	mode = M3Mode[mode_i];
    //Create histograms to store stuff in
    spectra.push_back(new TH2D("spectra_"+mode,"",nbinsx,binedges,nbins+1,0,maxy)); //400,0,/*20*/24/*100*/); //500,0,10
    hAvg.push_back(new TH1D("hAvg"+mode,"",1000,0,500));
  }


  double average=0;

  int nsteps=nthrows;
  std::cout << "Number of throws is " << nthrows << std::endl;

  TRandom3* rnd = new TRandom3(0);
  for(int mode_i = 0; mode_i < spectra.size() ; mode_i++){
		mode = M3Mode[mode_i];

	for(int i=0; i<nsteps; i++)
	{
	  TString name1=name+"_"+mode;
	  name1+=i;

	  temp = (TH1D*)file->Get(name1);
	  average+=temp->Integral();
	  hAvg[mode_i]->Fill(temp->Integral());

	  for(int j=1; j<=temp->GetNbinsX(); j++)
		spectra[mode_i]->Fill(temp->GetXaxis()->GetBinCenter(j),temp->GetBinContent(j));
	}
	std::cout << "Filled spectra " << std::endl;
  }
  outfile->cd();

  ///////////////////////////////
  //
  //  Now make projections
  //
  ///////////////////////////////

  std::cout << "Now looping through to make projections" << std::endl;
  double x[100], y[100], ex[100], ey[100];
  TH1D* proj=0;
  //vector to store all the TGraphs for each mode
  std::vector<TGraphErrors*> errorbars;
  TString dummyname;
  TString histname;
  std::vector<TH1D*> hist1d;

  for(int mode_i = 0 ; mode_i < spectra.size() ; mode_i++){	
	mode = SplineMode[mode_i];
	if (name=="sk_nue"){dummyname = "nue_"+mode;}
	if (name=="sk_nuebar"){dummyname = "nuebar_"+mode;}
	if (name=="sk_nue1pi"){dummyname = "nue1pi_"+mode;}
	if (name=="sk_numu"){dummyname = "numu_"+mode;}
	if (name=="sk_numubar"){dummyname = "numubar_"+mode;}

    std::cout << "Starting mode " << mode << " " << mode_i << " / " << spectra.size() << std::endl;
	// make TH1D with y-errors only
	histname = "hist_1d_"+dummyname;
	hist1d.push_back((TH1D*)temp->Clone(histname));
	hist1d[mode_i]->Reset();
	hist1d[mode_i]->GetYaxis()->SetRange(0,20);
		//Some hard coded rebinning based on mode and sample
		//Some modes don't have many events in so we need to
		//rebin
		if(name == "sk_numu" || name == "sk_numubar"){
	  	if(mode_i != 0 && mode_i != 4){
				spectra[mode_i]->RebinX(4);
				hist1d[mode_i]->RebinX(4);
	  	}
	  	else if(mode_i == 4){
        spectra[mode_i]->RebinX(8);
				hist1d[mode_i]->RebinX(8);
	  	}
		}
		else{
	  	if(mode_i>3){
      	spectra[mode_i]->RebinX(2);
				hist1d[mode_i]->RebinX(2);
	  	}
		}

		for(int i=1; i<=spectra[mode_i]->GetNbinsX() ; i++){

	  	//std::cout << "on bin " << i << std::endl;
	  	TString namepy=name+"_"+mode+"_py";
	  	namepy+=i;
	  	proj=spectra[mode_i]->ProjectionY(namepy,i,i);

	  	int first_bin = 99999;
	  	int last_bin = 0;

	  	for(int bin_i = 1 ; bin_i < proj->GetNbinsX() - 1 ; bin_i++){
				double val = proj->GetBinContent(bin_i);
				if(val > 0 && bin_i < first_bin){
		  		first_bin = bin_i;
				}

				if(val > 0){
		  		last_bin = bin_i;
				}
	  	}

		  //Get the range the filled bins cover
		  double diff = proj->GetBinCenter(last_bin) - proj->GetBinCenter(first_bin);
		  double opt_bin_width = diff/25;
		  double current_bin_width = proj->GetBinWidth(last_bin);
		  int rebin = 1;
		  if(diff > 0){rebin= (opt_bin_width/current_bin_width);}
	  	//std::cout << "rebin is " << rebin << std::endl;
	  	//std::cout << "diff is " << diff << std::endl;
	  	if(rebin==0){rebin+=1;}
	  	while((proj->GetNbinsX() % rebin) != 0){
				rebin--;
	  	}

	  	if(rebin > 0){
				proj->Rebin(rebin);
	  	}

		  //Get the bin centre array
		  // i - 1 is just because TGraphs bins are different to TH1
		  x[i-1]=spectra[mode_i]->GetXaxis()->GetBinCenter(i);
		  //Set the error on as being half the bin width
		  ex[i-1]=spectra[mode_i]->GetXaxis()->GetBinWidth(i)/2.0;
		  //Set the centre of the bin to be
		  //y[i-1]=proj->GetBinCenter(proj->GetMaximumBin());
		  y[i-1]=proj->GetMean();

		  //Set the error on y to be RMS of the bin
		  ey[i-1]=proj->GetRMS();

		  //TFitResultPtr r = proj->Fit("gaus","QS","goff",y[i-1]-ey[i-1],y[i-1]+ey[i-1]);
		  //if(y[i-1]>0.1 && r->Parameter(1)>0){y[i-1]=r->Parameter(1);}
		  //      c->Update();

	  	if(y[i-1]> 1E-4){
				if(mode_i==0){std::cout << "Filling bin " << i << " with " << y[i-1] << std::endl;}
				hist1d[mode_i]->SetBinContent(i,y[i-1]);
				hist1d[mode_i]->SetBinError(i,ey[i-1]);
	  	}
	 		else{
				hist1d[mode_i]->SetBinContent(i,0);
				hist1d[mode_i]->SetBinError(i,0);
	  	}

	  	proj->Write(namepy);
		}

		errorbars.push_back(new TGraphErrors(nbinsx,x,y,ex,ey));

		errorbars[mode_i]->Write("graph_"+dummyname);
		histname="spectra2d_"+dummyname;
		spectra[mode_i]->Write(histname);
		histname="hist_1d_"+dummyname;
		hist1d[mode_i]->Write(histname);
		histname=name+"_nominal";
		//nominal[mode_i]->Write(histname);
		histname="hist_1d_"+dummyname+"_events";
		hAvg[mode_i]->Write(histname);
  }

  file->Close();
  //delete file;
  //delete outfile;
  //delete spectra;
  //delete errorbars;
  //delete rnd;
  //delete c;
}

// -------------------------------------------------------------------------------- //
// ETA - added several knew options to this function for OA 2020. Can do parse in whether
// you have done posteriorPredictions by mode (do_by_mode = true in posteriorPredictionsSK2020.cpp)
// and also make spectra in different variables ( if PLOT1DVARHISTS is set in config then
// the posterior predictive code will make throws in different kinematic variables, see
// samplePDF/Structs.h ND280KinematicTypes for info. Despite these options you can just run
// the code as before (omitting the last two variables) and all will be well.
// You will probably need to change the binning range in the makeSKspectra_sample
// calls (assuming we have more POT!!!).

void makeSKspectra(const char *infile="sk_nominal.root", TString outname = "output.root", int nthrows = 1, int do_var = -1, bool do_by_mode = false, int which_oa = 2023)
{
  if(!do_by_mode){
  	bool do_var_bool = false;
  	TString var;
	  switch(do_var){
		  case 0:
				do_var_bool = true;
				var = "plep";
				break;
		  case 1:
				do_var_bool = true;
				var = "costheta_lep";
				break;
		  case 2:
				do_var_bool = true;
				var = "trueEnu";
				break;
		  case 3:
				do_var_bool = true;
				var = "trueQ2";
				break;
		  case 4:
				do_var_bool = true;
				var = "ErecQE";
				break;
		  case 5:
				do_var_bool = true;
				var = "Q2QE";
				break;
		  default:
				break;
	  }
	  switch(which_oa){
	  	case 2023:
	  		std::cout << "Generating OA2023 spectra" << std::endl;
	  		samplenameuser = samplenameOA2023;
			break;
	  	case 2024:
	  		std::cout << "Generating OA2024 spectra" << std::endl;
	  		samplenameuser = samplenameOA2024;
	  		break;
	  	default:
			std::cerr <<"Inserted " << which_oa<<" OA " << std::endl;
	  		std::cerr <<"Wrong OA year, please correct" << std::endl;
	  	 	break;
	  }

	  outfile = new TFile(outname+var+".root","recreate");


	  //Show samples in different kinematic variable
	  for(int i=0; i<nSamples; ++i){
	  	
	  	std::cout << std::endl << " ---------------- " << samplenameuser[i] << " ---------------- " << std::endl << std::endl;
	  	
	  	if(do_var_bool){
	  		makeSKspectra_sample(infile, samplenameuser[i]+"_"+var, outname, maxysample[i], nthrows, do_var, var);
	  	}	
	  	else{
	  		makeSKspectra_sample(infile, samplenameuser[i], outname, maxysample[i], nthrows);
	  	}
	  }
  }

  else{ //Samples split by mode
  	
  	for(int i=0; i<nSamples; ++i){
  		std::cout << std::endl << " ---------------- " << samplenameuser[i] << " ---------------- " << std::endl << std::endl;
  		makeSKspectra_sample_by_mode(infile, samplenameuser[i], outname, maxysample[i], nthrows);
  	}
  }		
  outfile->Close();
}

// -------------------------------------------------------------------------------- //
void makeSKspectra_beta_sample(const char *infile="sk_nominal_beta.root", TString name = "hist", int nbinsy = 1, int maxy = 1)
{

  TFile* file = new TFile(infile, "open");

  TString tmpname = name+"0";
  TString tmpname_beta0 = name+"_beta0_0";
  TString tmpname_beta1 = name+"_beta1_1";
  TH1D* temp = (TH1D*)file->Get(tmpname);
  TH1D* temp_beta0 = (TH1D*)file->Get(tmpname_beta0);
  TH1D* temp_beta1 = (TH1D*)file->Get(tmpname_beta1);
  int nbinsx = temp->GetNbinsX(); // Assume binning is the same for all three
  const int n = nbinsx;
  //double binedges[n+1];
  std::vector<double> binedges(nbinsx + 1);

  printf("%i \n",nbinsx);

  for(int i=0; i<=nbinsx; i++)
  {
	binedges[i]=temp->GetXaxis()->GetBinLowEdge(i+1);
  }
  //binedges[nbinsx+1] = temp->GetXaxis()->GetBinUpEdge(nbinsx);

  //TH2D* spectra = new TH2D("spectra","spectra",nbinsx,binedges.data(),nbinsy,0,maxy);
  // Assuming 'nbins' is the number of bins on Y and 'maxy' is the upper bound
  TH2D* spectra = new TH2D("test", "test", nbinsx, binedges.data(), nbinsx + 1, 0, maxy);
  TH2D* spectra_beta0 = new TH2D("spectra_beta0","spectra_beta0",nbinsx,binedges.data(),nbinsy,0,maxy);
  TH2D* spectra_beta1 = new TH2D("spectra_beta1","spectra_beta1",nbinsx,binedges.data(),nbinsy,0,maxy);

  double average=0, average_beta0=0, average_beta1=0;

  // I don't remember what these histograms were for or why they were useful,
  // and I never use them anyway, so commented out
  //TH1D* hAvg = new TH1D("hAvg","hAvg",200,0,10);
  //TH1D* hAvg_beta0 = new TH1D("hAvg_beta0","hAvg_beta0",200,0,10);
  //TH1D* hAvg_beta1 = new TH1D("hAvg_beta1","hAvg_beta1",200,0,10);
  TRandom3* rnd = new TRandom3(0);

  for(int i=0; i<1000; i++)
  {
	TString name1=name;
	TString name_beta0=name+"_beta0_";
	TString name_beta1=name+"_beta1_";
	name1+=i;
	name_beta0+=i;
	name_beta1+=i;

	temp = (TH1D*)file->Get(name1);
	average+=temp->Integral();
	//hAvg->Fill(temp->Integral());

	temp_beta0 = (TH1D*)file->Get(name_beta0);
	average_beta0+=temp_beta0->Integral();
	//hAvg_beta0->Fill(temp_beta0->Integral());

	temp_beta1 = (TH1D*)file->Get(name_beta1);
	average_beta1+=temp_beta1->Integral();
	//hAvg_beta1->Fill(temp_beta1->Integral());

	for(int j=1; j<=temp->GetNbinsX(); j++)
	{
	  spectra->Fill(temp->GetXaxis()->GetBinCenter(j),temp->GetBinContent(j));
	  spectra_beta0->Fill(temp_beta0->GetXaxis()->GetBinCenter(j),temp_beta0->GetBinContent(j));
	  spectra_beta1->Fill(temp_beta1->GetXaxis()->GetBinCenter(j),temp_beta1->GetBinContent(j));
	}
  }

  average/=2000.0;
  average_beta0/=2000.0;
  average_beta1/=2000.0;
  //   printf("combined: %f %f %f\n",average,hAvg->GetMean(),hAvg->GetRMS());
  //   printf("beta=0: %f %f %f\n",average_beta0,hAvg_beta0->GetMean(),hAvg_beta0->GetRMS());
  //   printf("beta=1: %f %f %f\n",average_beta1,hAvg_beta1->GetMean(),hAvg_beta1->GetRMS());
  printf("combined average: %f\n",average);
  printf("beta=0 average: %f\n",average_beta0);
  printf("beta=1 average: %f\n",average_beta1);


  TCanvas* c = new TCanvas("c","c",0,0,700,700);
  spectra->Draw("colz");
  c->Update();
  //   hAvg->Draw();
  //c->Update();
  spectra_beta0->Draw("colz");
  c->Update();
  //   hAvg_beta0->Draw();
  //c->Update();
  spectra_beta1->Draw("colz");
  c->Update();
  //   hAvg_beta1->Draw();
  //c->Update();

  outfile->cd();

  double x[100], y[100], ex[100], ey[100], x_beta0[100], y_beta0[100], ex_beta0[100], ey_beta0[100], x_beta1[100], y_beta1[100], ex_beta1[100], ey_beta1[100];
  c->Clear();

  TH1D* proj=0;
  for(int i=1; i<=nbinsx; i++)
  {
	// Combined
	std::cout << "-------------------- combined (" << i << ")-------------------" << std::endl;
	TString namepy_beta=name+"_py";
	namepy_beta+=i;
	proj=spectra->ProjectionY(namepy_beta,i,i);

	x[i-1]=spectra->GetXaxis()->GetBinCenter(i);
	ex[i-1]=spectra->GetXaxis()->GetBinWidth(i)/2.0;
	y[i-1]=proj->GetMean();
	ey[i-1]=proj->GetRMS();
	TFitResultPtr r = proj->Fit("gaus","QS","goff",y[i-1]-ey[i-1],y[i-1]+ey[i-1]);
	if(y[i-1]>0.02 && r->Parameter(1)>0)
	  y[i-1]=r->Parameter(1);
	proj->Write(namepy_beta);

	// beta=0
	std::cout << "-------------------- beta = 0 (" << i << ")-------------------" << std::endl;
	TString namepy_beta0=name+"_py_beta0";
	namepy_beta0+=i;
	proj=spectra_beta0->ProjectionY(namepy_beta0,i,i);

	x_beta0[i-1]=spectra_beta0->GetXaxis()->GetBinCenter(i);
	ex_beta0[i-1]=spectra_beta0->GetXaxis()->GetBinWidth(i)/2.0;
	y_beta0[i-1]=proj->GetMean();
	ey_beta0[i-1]=proj->GetRMS();
	TFitResultPtr r_beta0 = proj->Fit("gaus","QS","goff",y_beta0[i-1]-ey_beta0[i-1],y_beta0[i-1]+ey_beta0[i-1]);
	if(y_beta0[i-1]>0.02 && r_beta0->Parameter(1)>0)
	  y_beta0[i-1]=r_beta0->Parameter(1);
	proj->Write(namepy_beta0);

	// beta=1
	std::cout << "-------------------- beta = 1 (" << i << ")-------------------" << std::endl;
	TString namepy_beta1=name+"_py_beta1";
	namepy_beta1+=i;
	proj=spectra_beta1->ProjectionY(namepy_beta1,i,i);


	x_beta1[i-1]=spectra_beta1->GetXaxis()->GetBinCenter(i);
	ex_beta1[i-1]=spectra_beta1->GetXaxis()->GetBinWidth(i)/2.0;
	y_beta1[i-1]=proj->GetMean();
	ey_beta1[i-1]=proj->GetRMS();
	TFitResultPtr r_beta1 = proj->Fit("gaus","QS","goff",y_beta1[i-1]-ey_beta1[i-1],y_beta1[i-1]+ey_beta1[i-1]);
	if(y_beta1[i-1]>0.02 && r_beta1->Parameter(1)>0)
	  y_beta1[i-1]=r_beta1->Parameter(1);
	proj->Write(namepy_beta1);


  }

  TGraphErrors* errorbars = new TGraphErrors(nbinsx,x,y,ex,ey);
  errorbars->Draw("A E2");

  TGraphErrors* errorbars_beta0 = new TGraphErrors(nbinsx,x_beta0,y_beta0,ex_beta0,ey_beta0);
  errorbars_beta0->Draw("A E2");

  TGraphErrors* errorbars_beta1 = new TGraphErrors(nbinsx,x_beta1,y_beta1,ex_beta1,ey_beta1);
  errorbars_beta1->Draw("A E2");

  // make TH1D with y-errors only
  TString histname = name+"_hist1d";
  TH1D *hist1d = (TH1D*)temp->Clone(histname);
  hist1d->Reset();
  for (int i=0; i<hist1d->GetXaxis()->GetNbins(); i++)
  {
	hist1d->SetBinContent(i+1,y[i]);
	hist1d->SetBinError(i+1,ey[i]);
  }

  histname = name+"_hist1d_beta0";
  TH1D *hist1d_beta0 = (TH1D*)temp_beta0->Clone(histname);
  hist1d_beta0->Reset();
  for (int i=0; i<hist1d_beta0->GetXaxis()->GetNbins(); i++)
  {
	hist1d_beta0->SetBinContent(i+1,y_beta0[i]);
	hist1d_beta0->SetBinError(i+1,ey_beta0[i]);
  }

  histname = name+"_hist1d_beta1";
  TH1D *hist1d_beta1 = (TH1D*)temp_beta1->Clone(histname);
  hist1d_beta1->Reset();
  for (int i=0; i<hist1d_beta1->GetXaxis()->GetNbins(); i++)
  {
	hist1d_beta1->SetBinContent(i+1,y_beta1[i]);
	hist1d_beta1->SetBinError(i+1,ey_beta1[i]);
  }

  TString dummyname;
  if (name=="sk_nue"){dummyname = "nue";}
  if (name=="sk_nuebar"){dummyname = "nuebar";}
  if (name=="sk_nue1pi"){dummyname = "nue1pi";}
  if (name=="sk_numu"){dummyname = "numu";}
  if (name=="sk_numubar"){dummyname = "numubar";}

  outfile->cd();
  errorbars->Write(name);
  histname = "spectra2d_"+dummyname;
  spectra->Write(histname);
  histname = "hist_1d_"+dummyname;
  hist1d->Write(histname);
  histname = name+"_beta0";
  errorbars_beta0->Write(histname);
  histname = "spectra2d_"+dummyname+"_beta0";
  spectra_beta0->Write(histname);
  histname = "hist_1d_"+dummyname+"_beta0";
  hist1d_beta0->Write(histname);
  histname = name+"_beta1";
  errorbars_beta1->Write(histname);
  histname = "spectra2d_"+dummyname+"_beta1";
  spectra_beta1->Write(histname);
  histname = "hist_1d_"+dummyname+"_beta1";
  hist1d_beta1->Write(histname);

  // Finally, calculate likelihood ratio between beta=1 and beta=0 (NOTE: THIS WILL ONLY WORK IF MC IS BINNED THE SAME AS THE DATA)
  // Commented out for now
  /*double negLogL_beta0=0, negLogL_beta1=0;
	TFile *f_dat = new TFile("inputs/Run5-6Data_June2015.root","open");
	TH1D *dathist = (TH1D*)f_dat->Get("nue");

	for (int i=0; i<hist1d_beta0->GetXaxis()->GetNbins(); i++)
	{
	double mc = hist1d_beta0->GetBinContent(i+1);
	double dat = dathist->GetBinContent(i+1);
	std::cout << mc << "  " << dat << std::endl;
	if (dat == 0)
	negLogL_beta0 += (mc - dat);
	else
	negLogL_beta0 += (mc - dat + dat*TMath::Log(dat/mc));

	std::cout << negLogL_beta0 << std::endl;
	}

	for (int i=0; i<hist1d_beta1->GetXaxis()->GetNbins(); i++)
	{
	double mc = hist1d_beta1->GetBinContent(i+1);
	double dat = dathist->GetBinContent(i+1);
	std::cout << mc << "  " << dat << std::endl;
	if (dat ==0)
	negLogL_beta1 += (mc - dat);
	else
	negLogL_beta1 += (mc - dat + dat*TMath::Log(dat/mc));

	std::cout << negLogL_beta1 << std::endl;
	}

	std::cout << "negLogL_beta0 = " << negLogL_beta0 << ", negLogL_beta1 = " << negLogL_beta1 << std::endl;
	std::cout << "Likelihood ratio = " << TMath::Exp(negLogL_beta1-negLogL_beta0) << std::endl;
	*/
}

void makeSKspectra_beta(const char *infile="sk_nominal.root")
{
  outfile = new TFile("spectra_beta.root","recreate");
  std::cout << std::endl << " ---------------- numu (beta) ----------------" <<std::endl << std::endl;
  makeSKspectra_beta_sample(infile, "sk_numu", 400, 80);
  std::cout << std::endl << " ---------------- nue (beta) ----------------" <<std::endl << std::endl;
  makeSKspectra_beta_sample(infile, "sk_nue", 400, 14);
  std::cout << std::endl << " ---------------- nue1pi (beta) ----------------" <<std::endl << std::endl;
  makeSKspectra_beta_sample(infile, "sk_nue1pi", 400, 14);
  std::cout << std::endl << " ---------------- numubar (beta) ----------------" <<std::endl << std::endl;
  makeSKspectra_beta_sample(infile, "sk_numubar", 400, 10);
  std::cout << std::endl << " ---------------- nuebar (beta) ----------------" <<std::endl << std::endl;
  makeSKspectra_beta_sample(infile, "sk_nuebar", 100, 6);
  outfile->Close();
}

// -------------------------------------------------------------------------------- //



