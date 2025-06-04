// ROOT script to convert the DUNE flux systematics ROOT file 'total_covariance_DUNE_opt.root'
// into a format used by MaCh3.
// It will be input to a python script Flux_RootCovToXML_DUNE.py that outputs an XML file format
// for the covariances

#include "TMatrixD.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"

#include <iostream>

 // The covariance matrix is a 208x208 matrix with bins corresponding to:
// 1-19: Near detector, neutrino mode, nu_mu 19
// 20-38: Near detector, neutrino mode, nu_mu_bar 19
// 39-45: Near detector, neutrino mode, nu_e 7
// 46-52: Near detector, neutrino mode, nu_e_bar 7
// 53-71: Near detector, antineutrino mode, nu_mu 19
// 72-90: Near detector, antineutrino mode, nu_mu_bar 19
// 91-97: Near detector, antineutrino mode, nu_e 7
// 98-104: Near detector, antineutrino mode, nu_e_bar 7
// 105-123: Far detector, neutrino mode, nu_mu 19
// 124-142: Far detector, neutrino mode, nu_mu_bar 19
// 143-149: Far detector, neutrino mode, nu_e 7
// 150-156: Far detector, neutrino mode, nu_e_bar 7
// 157-175: Far detector, antineutrino mode, nu_mu 19
// 176-194: Far detector, antineutrino mode, nu_mu_bar 19
// 195-201: Far detector, antineutrino mode, nu_e 7
// 202-208: Far detector, antineutrino mode, nu_e_bar 7


const int Nbins_numu = 19;
const double edges_numu[Nbins_numu+1] = {0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 12, 16, 20, 40, 100};

const int Nbins_nue = 7;
const double edges_nue[Nbins_nue+1] = {0, 2, 4, 6, 8, 10, 20, 100};

const int Nspectra = 16;

const int Nbins = (Nspectra/2) * Nbins_numu + (Nspectra/2) * Nbins_nue;


TMatrixD mat(TH2* h)
{
  TMatrixD m(Nbins, Nbins);
  for(int i = 1; i <= Nbins; ++i){
    for(int j = 1; j <= Nbins; ++j){
      m(i-1, j-1) = h->GetBinContent(i, j);
    }
  }
  return m;
}

void make_flux_systs_DUNE(){

  // Output file
  TFile * fout = new TFile("DUNE_flux_systs.root", "RECREATE");

  // Pick up the DUNE  total covariance input: note this is fractional, so will
  // have to multiply each entry by the corresponding absolute flux (see below)
  TFile* fcov = new TFile("total_covariance_DUNE_opt.root", "READ");
  TH2* hcov = (TH2*)fcov->Get("total_covariance");

  TH1D * hfluxes = new TH1D("hfluxes", "hfluxes", Nbins, 0, Nbins);

  std::string det, hc, flav;
  
  // Pick up all the 16 nominal flux predictions, only used for the axis bins
  for(std::string hcStr: {"neutrino", "antineutrino"}){
    if(hcStr == "neutrino")
      hc = std::string("numode");
    else
      hc = std::string("anumode");
    for(std::string detStr: {"ND", "FD"}){
      TFile* f = new TFile(("histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_"+hcStr+"_LBNE"+detStr+"_fastmc.root").c_str(), "READ");
      if(detStr == "ND")
	det = std::string("nd5");
      else
	det = std::string("sk");
      
      for(std::string flavStr: {"numu", "numubar", "nue", "nuebar"}){
        TH1* h = (TH1*)f->Get((flavStr+"_flux").c_str());
	if(flavStr == "numu")
	  flav = "numu";
	if(flavStr == "numubar")
	  flav = "numub";
	if(flavStr == "nue")
	  flav = "nue";
	if(flavStr == "nuebar")
	  flav = "nueb";
	

	if(flavStr == std::string("numu") || flavStr == std::string("numubar")){ 
	  //	  TH1D * hrebin = (TH1D*)h->Rebin(Nbins_numu, (detStr + "_" + hcStr + "_" + flavStr + "_flux").c_str(), edges_numu);
	  TH1D * hrebin = (TH1D*)h->Rebin(Nbins_numu, (det + "_" + hc + "_" + flav + "_bins").c_str(), edges_numu);
	  hrebin->SetTitle((detStr + "_" + hcStr + "_" + flavStr + "_flux").c_str());
	  TAxis * hax = hrebin->GetXaxis();
	  hax->SetName((det + "_" + hc + "_" + flav + "_bins").c_str());
	  fout->cd();
	  //hrebin->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
	  //hrebin->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
	  //	  hrebin->Write();
	  hax->Write();
	  f->cd();
	  //hrebin->Scale(1/hrebin->Integral(0, -1));// Why???

	  // This is just filling the bins of the stitched histo hfluxes...not very elegant
	  if(detStr == std::string("ND")){
	    if(hcStr == std::string("neutrino")){
	      if(flavStr == std::string("numu")){
		for(int i = 1; i <= Nbins_numu; i++)
		  hfluxes->SetBinContent(i, hrebin->GetBinContent(i));
	      }else{//numubar
		for(int i = 1; i <= Nbins_numu; i++)
		  hfluxes->SetBinContent(i+Nbins_numu, hrebin->GetBinContent(i));
	      }
	    }else{// antineutrino
	      if(flavStr == std::string("numu")){
		for(int i = 1; i <= Nbins_numu; i++)
		  hfluxes->SetBinContent(i+2*Nbins_numu+2*Nbins_nue, hrebin->GetBinContent(i));
	      }else{//numubar
		for(int i = 1; i <= Nbins_numu; i++)
		  hfluxes->SetBinContent(i+2*Nbins_numu+2*Nbins_nue+Nbins_numu, hrebin->GetBinContent(i));
	      }
	    }
	  }else{// FD
	    if(hcStr == std::string("neutrino")){
	      if(flavStr == std::string("numu")){
		for(int i = 1; i <= Nbins_numu; i++)
		  hfluxes->SetBinContent(i+4*Nbins_numu+4*Nbins_nue, hrebin->GetBinContent(i));
	      }else{//numubar
		for(int i = 1; i <= Nbins_numu; i++)
		  hfluxes->SetBinContent(i+4*Nbins_numu+4*Nbins_nue+Nbins_numu, hrebin->GetBinContent(i));
	      }
	    }else{// antineutrino
	      if(flavStr == std::string("numu")){
		for(int i = 1; i <= Nbins_numu; i++)
		  hfluxes->SetBinContent(i+6*Nbins_numu+6*Nbins_nue, hrebin->GetBinContent(i));
	      }else{// numubar
		for(int i = 1; i <= Nbins_numu; i++)
		  hfluxes->SetBinContent(i+6*Nbins_numu+6*Nbins_nue+Nbins_numu, hrebin->GetBinContent(i));
	      }
	    }
	  }

	      

	}else{// nue, nuebar
	  //	  TH1D * hrebin = (TH1D*)h->Rebin(Nbins_nue, (detStr + "_" + hcStr + "_" + flavStr + "_flux").c_str(), edges_nue);
	  TH1D * hrebin = (TH1D*)h->Rebin(Nbins_nue, (det + "_" + hc + "_" + flav + "_bins").c_str(), edges_nue);
	  hrebin->SetTitle((detStr + "_" + hcStr + "_" + flavStr + "_flux").c_str());
	  TAxis * hax = hrebin->GetXaxis();
	  hax->SetName((det + "_" + hc + "_" + flav + "_bins").c_str());
	  fout->cd();
	  //hrebin->Write();
	  hax->Write();
	  f->cd();


	  //	  hrebin->Scale(1/hrebin->Integral(0, -1));// Why??
	  // This is just filling the bins of the stitched histo hfluxes...not very elegant
	  if(detStr == std::string("ND")){
	    if(hcStr == std::string("neutrino")){
	      if(flavStr == std::string("nue")){
		for(int i = 1; i <= Nbins_nue; i++)
		  hfluxes->SetBinContent(i+2*Nbins_numu, hrebin->GetBinContent(i));
	      }else{//nuebar
		for(int i = 1; i <= Nbins_nue; i++)
		  hfluxes->SetBinContent(i+2*Nbins_numu+Nbins_nue, hrebin->GetBinContent(i));
	      }
	    }else{// antineutrino
	      if(flavStr == std::string("nue")){
		for(int i = 1; i <= Nbins_nue; i++)
		  hfluxes->SetBinContent(i+2*Nbins_numu+2*Nbins_nue+2*Nbins_numu, hrebin->GetBinContent(i));
	      }else{//nuebar
		for(int i = 1; i <= Nbins_nue; i++)
		  hfluxes->SetBinContent(i+2*Nbins_numu+2*Nbins_nue+2*Nbins_numu+Nbins_nue, hrebin->GetBinContent(i));
	      }
	    }
	  }else{// FD
	    if(hcStr == std::string("neutrino")){
	      if(flavStr == std::string("nue")){
		for(int i = 1; i <= Nbins_nue; i++)
		  hfluxes->SetBinContent(i+6*Nbins_numu+4*Nbins_nue, hrebin->GetBinContent(i));
	      }else{//nuebar
		for(int i = 1; i <= Nbins_nue; i++)
		  hfluxes->SetBinContent(i+6*Nbins_numu+4*Nbins_nue+Nbins_nue, hrebin->GetBinContent(i));
	      }
	    }else{// antineutrino
	      if(flavStr == std::string("nue")){
		for(int i = 1; i <= Nbins_nue; i++)
		  hfluxes->SetBinContent(i+8*Nbins_numu+6*Nbins_nue, hrebin->GetBinContent(i));
	      }else{// nuebar
		for(int i = 1; i <= Nbins_nue; i++)
		  hfluxes->SetBinContent(i+8*Nbins_numu+6*Nbins_nue+Nbins_nue, hrebin->GetBinContent(i));
	      }
	    }
	  }

	  
	}
      }
    }
  }

  fout->cd();
  
  hfluxes->Write("hfluxes");

  // Normalize to absolute flux uncertainty
  for(int i = 1; i <= Nbins; ++i){
    for(int j = 1; j <= Nbins; ++j){
      const double z = hcov->GetBinContent(i, j);
            hcov->SetBinContent(i, j, z * hfluxes->GetBinContent(i) * hfluxes->GetBinContent(j));
    }
  }

  TMatrixD mcov = mat(hcov);
  
  mcov.Write("total_flux_cov");
  
  fout->Close();


    

}
