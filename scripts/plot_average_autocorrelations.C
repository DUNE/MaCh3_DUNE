#include "TH1D.h"
#include "TString.h"
#include "TFile.h"
#include "THStack.h"
#include "TTree.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TPad.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include <iostream>
#include <string>
//Plot autocorrelations and traces and save to png

void plot_average_ac(TString diagfile, TString output)
{
  TFile *fin = new TFile(diagfile, "open");

  if(fin->IsZombie()){
   std::cerr<<"Error couldnt open "<<diagfile<<std::endl;
   throw;
  }

  TDirectoryFile* autocor = (TDirectoryFile*)fin->Get("Auto_corr");

  if(!autocor){
    std::cerr<<"Error couldn't find Auto_Corr in "<<diagfile<<std::endl;
    throw;
  }


  TKey* key;
  TIter next(autocor->GetListOfKeys());
  
  std::cout<<"loaded in files"<<std::endl;
  
  TCanvas* c = new TCanvas("c", "c", 1200, 600);
  c->Draw();
  c->cd();
  c->SetGrid();

  std::cout<<"Beginning loop!"<<std::endl;

  gStyle->SetPalette(kRainbow);

  TH1D* auto_hist=nullptr;

  int npars = 0;

  THStack* s = new THStack();

  TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  while( (key = (TKey*) next()))
  { 
      TString auto_syst(key->GetName());
      std::cout<<"Getting AC info"<<auto_syst<<std::endl;

      auto split_string = auto_syst.Tokenize("_");
      auto param_index = split_string->At(1);

      //needs a TH2D for everything to go in number of params
      TH2D* h_q0q3 = new TH2D("h_q0q3", "q_{0} vs q_{3}; q_{3} [GeV]; q_{0} [GeV]",
                          q3_edges.size() - 1, &q3_edges[0],
                          q0_edges.size() - 1, &q0_edges[0]);
      h_q0q3->SetDirectory(nullptr);

  

      TH1D* plotting_hist = dynamic_cast<TH1D*>(autocor->Get(auto_syst.c_str())->Clone());
      plotting_hist->SetLineColorAlpha(kAzure-2, 0.1);
      double integral_auto_corr = plotting_hist->Integral();

        int nbins_q3 = q3_edges.size() - 1;
    int nbins_q0 = q0_edges.size() - 1;
    for (int i = 0; i < nbins_q0 * nbins_q3; ++i) {
        int bin_q0 = i / nbins_q3 + 1;
        int bin_q3 = i % nbins_q3 + 1;
        double q0 = h_q0q3->GetYaxis()->GetBinCenter(bin_q0);
        double q3 = h_q0q3->GetXaxis()->GetBinCenter(bin_q3);
        if (q0 <= q3) h_q0q3->SetBinContent(bin_q3, bin_q0, h_flat->GetBinContent(i + 1));
    }
      std::cout << "Integral of autocorrelation = " << integral_auto_corr << std::endl;
  


}
