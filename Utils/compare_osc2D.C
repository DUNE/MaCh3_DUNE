////////////////////////////////////////////////////////////////////////////
//
// Quick plotting script to show the affect of marginalising 2D distributions 
// Give this exe a reduced chain and it'll draw th23 vs. dm32 and show the 1D
// and 2D best-fit points and make panels either side of the main canvas
// showing the 1D distributions 
//  
/////////////////////////////////////////////////////////////////////////////

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "TFile.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TMarker.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TTree.h"
#include "TROOT.h"
#include "TColor.h"
#include "TImage.h"

void compare_osc2D(const char * infile1, const char * infile2){

  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPalette(51,0);
  gStyle->SetNumberContours(104);

  std::vector<std::string> title = {"#nu_{#mu} #rightarrow #nu_{#mu}","#nu_{#mu} #rightarrow #nu_{e}","#nu_{#mu} #rightarrow #nu_{#tau}", "#nu_{e} #rightarrow #nu_{#mu}","#nu_{e} #rightarrow #nu_{e}","#nu_{e} #rightarrow #nu_{#tau}", "#nu_{#tau} #rightarrow #nu_{#mu}","#nu_{tau} #rightarrow #nu_{e}","#nu_{#tau} #rightarrow #nu_{#tau}"};
  std::vector<std::string> title_nu = {"#bar{#nu}_{#mu}","#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e}","#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{#tau}"};
  std::vector<std::string> filename = {"mm", "me", "mt", "em", "ee", "et", "tm", "te", "tt"};

  TFile *f1 = new TFile(infile1, "READ");
  TFile *f2 = new TFile(infile2, "READ");


  TH2D* hOsc1 = 0;
  TH2D* hOsc2 = 0;
  TH2D* hOscDiff = 0;

  //TKey *key1;
  //TKey *key2;

  //TIter next1( f1->GetListOfKeys());
  //TIter next2( f2->GetListOfKeys());
  
  TList * list = f1->GetListOfKeys();
 
 for (int i=0; i<list->GetEntries(); i++) {

 // while ((key1 = (TKey *) next1())) {
    //while ((key2 = (TKey *) next2())) {
    
      std::cout << "hello" << std::endl;
  
      TString keyname = list->At(i)->GetName();
      std::cout << keyname << std::endl;
      hOsc1 = (TH2D*)f1->Get(keyname);
      hOsc2 = (TH2D*)f2->Get(keyname);
      hOscDiff = (TH2D*)hOsc1->Clone(filename[i].c_str());
      hOscDiff->Add(hOsc2, -1);

      
      hOsc1->GetXaxis()->SetTitle("True Energy (GeV)");
      hOsc1->GetXaxis()->SetLabelSize(0.05);
      hOsc1->GetYaxis()->SetLabelSize(0.05);
      hOsc1->GetXaxis()->SetTitleSize(0.06);
      hOsc1->GetXaxis()->SetTitleOffset(0.9);
      hOsc1->GetYaxis()->SetTitle("cosine zenith (radians)");
      hOsc1->GetYaxis()->SetTitleSize(0.06);
      hOsc1->GetYaxis()->SetTitleOffset(0.8);
      hOsc1->GetYaxis()->SetTickLength(0);
      hOsc1->GetXaxis()->SetTickLength(0);
      hOsc1->GetYaxis()->SetMaxDigits(3);
      hOsc1->GetZaxis()->SetTitle("Osc Prob");
      hOsc1->GetZaxis()->SetTitleSize(0.06);
      hOsc1->GetZaxis()->SetLabelSize(0.04);
      hOsc1->GetZaxis()->SetTitleOffset(0.8);

      hOsc2->GetXaxis()->SetTitle("True Energy (GeV)");
      hOsc2->GetXaxis()->SetLabelSize(0.05);
      hOsc2->GetYaxis()->SetLabelSize(0.05);
      hOsc2->GetXaxis()->SetTitleSize(0.06);
      hOsc2->GetXaxis()->SetTitleOffset(0.9);
      hOsc2->GetYaxis()->SetTitle("cosine zenith (radians)");
      hOsc2->GetYaxis()->SetTitleSize(0.06);
      hOsc2->GetYaxis()->SetTitleOffset(0.8);
      hOsc2->GetYaxis()->SetTickLength(0);
      hOsc2->GetXaxis()->SetTickLength(0);
      hOsc2->GetYaxis()->SetMaxDigits(3);
      hOsc1->GetZaxis()->SetTitle("Osc Prob");
      hOsc1->GetZaxis()->SetTitleSize(0.06);
      hOsc1->GetZaxis()->SetLabelSize(0.04);
      hOsc1->GetZaxis()->SetTitleOffset(0.8);
      
      hOscDiff->GetXaxis()->SetTitle("True Energy (GeV)");
      hOscDiff->GetXaxis()->SetLabelSize(0.05);
      hOscDiff->GetYaxis()->SetLabelSize(0.04);
      hOscDiff->GetXaxis()->SetTitleSize(0.04);
      hOscDiff->GetXaxis()->SetTitleOffset(0.9);
      hOscDiff->GetYaxis()->SetTitle("cosine zenith (radians)");
      hOscDiff->GetYaxis()->SetTitleSize(0.06);
      hOscDiff->GetYaxis()->SetTitleOffset(0.8);
      hOscDiff->GetYaxis()->SetTickLength(0);
      hOscDiff->GetXaxis()->SetTickLength(0);
      hOscDiff->GetYaxis()->SetMaxDigits(3);
      hOscDiff->GetZaxis()->SetTitle("Osc Prob difference");
      hOscDiff->GetZaxis()->SetTitleSize(0.06);
      hOscDiff->GetZaxis()->SetLabelSize(0.04);
      hOscDiff->GetZaxis()->SetTitleOffset(0.8);

      hOsc1->SetTitle(title[i].c_str());
      hOsc2->SetTitle(title[i].c_str());
      hOscDiff->SetTitle(title[i].c_str());

      TCanvas C1 ("C", "C", 2000, 2000);
      TCanvas C2 ("C", "C", 2000, 2000);
      TCanvas Cdiff ("C", "C", 2000, 2000);

      C1.cd();
      hOsc1->Draw("colz");
      C1.Update();
      C1.SaveAs(("CUDAProb_liban_atm_poly_" + filename[i] + ".png").c_str());
      
      C2.cd();
      hOsc2->Draw("colz");
      C2.Update();
      C2.SaveAs(("CUDAProb_dan_atm_poly_" + filename[i] + ".png").c_str());
      C2.cd();

      Cdiff.cd();
      hOscDiff->Draw("colz");
      Cdiff.Update();
      Cdiff.SaveAs(("CUDAProb_validation_atm_poly_" + filename[i] + ".png").c_str());

      /* TLegend *leg = new TLegend(0.45,0.7,0.9,0.9);
      leg->AddEntry(plot, "P_{CUDAProb}", "l");
      leg->AddEntry(plot2, "(P_{Prob3++} - P_{CUDAProb}) / P_{Prob3++}", "l");
      leg->AddEntry(h1, "P_{Prob3++} #leq 0.02", "f");
      leg->SetTextSize(0.03);
      leg->AddEntry(plot3, "#nu_{#mu} #rightarrow #nu_{e}" "Matter Oscillations", "l");
      leg->AddEntry(plot4, "#bar{#nu}_{#mu} #rightarrow #bar{#nu}_{e} Vacuum Oscillations", "l");
      leg->Draw(); */
      

  
  }

  return;

}
