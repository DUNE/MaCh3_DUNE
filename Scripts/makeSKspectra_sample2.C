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
void makeSKspectra_sample(const char *infile = "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/Chain_testing/Enurec_Elep_test_largerxsecscale4_Process.root", TString name = "hist", TString outname = "output.root", int maxy = 1, int nthrows = 1, bool do_var = false, TString var ="")
{

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
  double binedges[100];

  std::cout << "BLAH" << std::endl;
  printf("%i \n",nbinsx);

  for(int i=0; i<=nbinsx; i++){
  	binedges[i]=temp->GetXaxis()->GetBinLowEdge(i+1);
  }

  std::cout << "Made Th2D spectra" << std::endl;
  TH2D* spectra = new TH2D("test","test",nbinsx,binedges,nbins+1,0,maxy); //400,0,/*20*/24/*100*/); //500,0,10
  spectra->Sumw2(true);
  std::cout << "Made Th2D spectra" << std::endl;

  double average=0;

  // I don't remember what these histograms were for or why they were useful,
  // and I never use them anyway, so commented out
  // ETA - they're useful for error on the event rate!!!!!
  TH1D* hAvg = new TH1D("hAvg","hAvg",1000,0,500);
  hAvg->Sumw2(true);

  int nsteps = nthrows;
  std::cout << "Number of throws is " << nthrows << std::endl;

  TRandom3* rnd = new TRandom3(0);
  for(int i = 0; i < nsteps; i++)
  {
		std::cout << "Getting histogram for throw " << nsteps << std::endl;
    if(reweight){
      reweight->GetEntry(i);
	  if( i % (nsteps/10) == 0){
	  	std::cout << "Weight is " << weight << std::endl;
	  }
	}

	TString name1 = "";
	name1.Form("%s_%i", name.Data(), i);
	//name1+=i;

	std::cout << "Looking for histogram with name " << name1 << std::endl;
	temp = (TH1D*)file->Get(name1);
	average+=temp->Integral();
	hAvg->Fill(temp->Integral(), weight);

	for(int j=1; j<=temp->GetNbinsX(); j++)
	  spectra->Fill(temp->GetXaxis()->GetBinCenter(j),temp->GetBinContent(j), weight);
  }
  
  std::cout << "Filled spectra " << std::endl;

  average/=(double)nsteps;
  printf("%f %f %f\n",average,hAvg->GetMean(),hAvg->GetRMS());
  //printf("%f \n",average);

  TCanvas* c = new TCanvas("c","c",0,0,700,700);
  spectra->Draw("colz");
  c->Update();
  //hAvg->Draw();
  //c->Update();

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

  for(int i=1; i<=nbinsx; i++)
  {
		std::cout << "on bin " << i << std::endl;
		TString namepy=name+"_py";
		namepy+=i;
		proj=spectra->ProjectionY(namepy,i,i);
		//proj->Sumw2(true);

		if(proj->GetEntries() > 0){

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

	  	// ETA - this is an attempt to automate
	  	// the bin width in the 1D projections.
	  	// Often we fit 1D gaussians to these projections
	  	// and for this the number of bins should really be
	  	// at least 10 across the range of the projection.

	  	//Get the range the filled bins cover
	  	double diff = proj->GetBinCenter(last_bin) - proj->GetBinCenter(first_bin);
	  	double opt_bin_width = diff/10;
	  	double current_bin_width = proj->GetBinWidth(last_bin);
	  	int rebin = 1;
	  	
	  	if(diff > 0){rebin= (opt_bin_width/current_bin_width);}
	  	//std::cout << "rebin is " << rebin << std::endl;
	  	//std::cout << "diff is " << diff << std::endl;
	  	if(rebin==0){rebin+=1;}
	  	//If the number of bins in the rebinning isn't a factor
	  	//of the total number of bins keep reducing it until it is
	  	while((proj->GetNbinsX() % rebin) != 0){
				rebin--;
	  	}

      //Just in case the above loop does a bad job
	  	if(rebin > 0){
				proj->Rebin(rebin);
	  	}

	  	if(diff / proj->GetBinWidth(last_bin) < 50){
	   	/*	
	   	std::cout << "This hasn't worked!! " << std::endl;
			std::cout << "Current width is " << current_bin_width << std::endl;
			std::cout << "Omptimal width is " << opt_bin_width << std::endl;
			std::cout << "Bin width now is " << proj->GetBinWidth(last_bin) << std::endl;
			std::cout << "Rebin is " << rebin << std::endl;
			std::cout << "First bin is " << first_bin << std::endl;
			std::cout << "Last bin is " << last_bin << std::endl;
			*/
	  	}

	  	//Get the bin centre array
	  	x[i-1]=spectra->GetXaxis()->GetBinCenter(i);
	  	//Set the error on as being half the bin width
	  	ex[i-1]=spectra->GetXaxis()->GetBinWidth(i)/2.0;
		  //Set the centre of the bin to be
		  //y[i-1]=proj->GetBinCenter(proj->GetMaximumBin());
		  y[i-1]=proj->GetMean();
		  //Set the error on y to be RMS of the bin
		  ey[i-1]=proj->GetRMS();
	    //Option to fit a gaussian
		  //TFitResultPtr r = proj->Fit("gaus","QS","goff",y[i-1]-ey[i-1],y[i-1]+ey[i-1]);
		  //if(y[i-1]>0.1 && r->Parameter(1)>0){y[i-1]=r->Parameter(1);}
		  //c->Update();
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
void makeSKspectra_sample_by_mode(const char *infile = "sk_nominal.root", TString name = "hist", TString outname = "output.root", int maxy = 1, int nthrows = 1)
{

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
    spectra.push_back(new TH2D("spectra_"+mode,"",nbinsx,binedges,nbins+1,0,maxy)); //400,0,/*20*/24/*100*/); //500,0,1...
  }
}