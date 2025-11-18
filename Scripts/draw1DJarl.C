#include "TH2D.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TMath.h"
#include <iostream>

// sindcp is a boolean indicating whether we want to reweight to put a flat prior on sin(dcp)
// instead of the usual flat in dcp
// hieararchy = 1-->NH, 0-->both, -1-->IH
void draw1DJarl(int hierarchy, bool sindcp = false)
{

  TFile* f = new TFile("/vols/t2k/users/ljw20/software/MaCh3_DUNE_clean/MaCh3_DUNE/jarlskog_DUNE_FD_nufit.root");
  TH1D* j_hist;
  if(!sindcp) {
    if(hierarchy == 1) {
      j_hist = (TH1D*)f->Get("jarlskog_NH");
      j_hist->SetTitle("Jarlskog Invariant, Normal Hierarchy;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;Probability Density");
    }
    else if(hierarchy == 0) {
      j_hist = (TH1D*)f->Get("jarlskog_both");
      //j_hist->SetTitle("Jarlskog Invariant, Both Hierarchies;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;Probability Density");
      j_hist->SetTitle(";J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;Probability Density");
    }
    else if(hierarchy == -1) {
      j_hist = (TH1D*)f->Get("jarlskog_IH");
      j_hist->SetTitle("Jarlskog Invariant, Inverted Hierarchy;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;Probability Density");
    }
    else {
     std::cout << "Error: invalid hierarchy option. 1 for NH, 0 for both, -1 for IH" <<std::endl;
     throw;
    }
  }
  else {
    if(hierarchy == 1) {
      j_hist = (TH1D*)f->Get("jarlskog_NH_flatsindcp");
      j_hist->SetTitle("Jarlskog Invariant, Normal Hierarchy;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;Probability Density");
    }
    else if(hierarchy == 0) {
      j_hist = (TH1D*)f->Get("jarlskog_both_flatsindcp");
      //j_hist->SetTitle("Jarlskog Invariant, Both Hierarchies;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;Probability Density");
      j_hist->SetTitle(";J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;Probability Density");
    }
    else if(hierarchy == -1) {
      j_hist = (TH1D*)f->Get("jarlskog_IH_flatsindcp");
      j_hist->SetTitle("Jarlskog Invariant, Inverted Hierarchy;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;Probability Density");
    }
    else {
     std::cout << "Error: invalid hierarchy option. 1 for NH, 0 for both, -1 for IH" <<std::endl;
     throw;
    }
  }

  j_hist->GetYaxis()->SetLabelSize(0.);

  j_hist->Rebin(5);
  TH1D* j_copy = (TH1D*)j_hist->Clone("j_copy");

  TH1D* j_cred_1sig = (TH1D*)j_hist->Clone("j_cred_1sig");
  TH1D* j_cred_90p = (TH1D*)j_hist->Clone("j_cred_90p");
  TH1D* j_cred_2sig = (TH1D*)j_hist->Clone("j_cred_2sig");
  TH1D* j_cred_99p = (TH1D*)j_hist->Clone("j_cred_99p");
  TH1D* j_cred_3sig = (TH1D*)j_hist->Clone("j_cred_3sig");

  double contlevel1=0.68;
  double contlevel2=0.90;
  double contlevel3=0.954;
  double contlevel4=0.99;
  double contlevel5=0.9973;

  double integral, tsum=0.;

  integral = j_copy->Integral();

  while((tsum/integral)<contlevel5) {

    double tmax = j_copy->GetMaximum();
    int bin = j_copy->GetMaximumBin();
    if((tsum/integral)<contlevel1) {
      j_copy->SetBinContent(bin,-1.0);
      j_cred_1sig->SetBinContent(bin,0.);
      j_cred_90p->SetBinContent(bin,0.);
      j_cred_2sig->SetBinContent(bin,0.);
      j_cred_99p->SetBinContent(bin,0.);
      j_cred_3sig->SetBinContent(bin,0.);
    }
    if((tsum/integral)<contlevel2  && (tsum / integral > contlevel1) ) {
      j_copy->SetBinContent(bin,-3.0);
      j_cred_90p->SetBinContent(bin,0.);
      j_cred_2sig->SetBinContent(bin,0.);
      j_cred_99p->SetBinContent(bin,0.);
      j_cred_3sig->SetBinContent(bin,0.);
    }
    if((tsum/integral)<contlevel3  && (tsum / integral > contlevel2) ) {
      j_copy->SetBinContent(bin,-5.0);
      j_cred_2sig->SetBinContent(bin,0.);
      j_cred_99p->SetBinContent(bin,0.);
      j_cred_3sig->SetBinContent(bin,0.);
    }
    if((tsum/integral)<contlevel4  && (tsum / integral > contlevel3) ) {
      j_copy->SetBinContent(bin,-7.0);
      j_cred_99p->SetBinContent(bin,0.);
      j_cred_3sig->SetBinContent(bin,0.);
    }
    if((tsum/integral)<contlevel5  && (tsum / integral > contlevel4) ) {
      j_copy->SetBinContent(bin,-9.0);
      j_cred_3sig->SetBinContent(bin,0.);
    }
    tsum+=tmax;
  }

  j_hist->SetLineColor(kBlack);
  j_cred_1sig->SetLineColor(kBlack);
  j_cred_90p->SetLineColor(kBlack);
  j_cred_2sig->SetLineColor(kBlack);
  j_cred_99p->SetLineColor(kBlack);
  j_cred_3sig->SetLineColor(kBlack);

  j_hist->SetFillColor(17);
  j_cred_1sig->SetFillColor(15);
  j_cred_90p->SetFillColor(15);
  j_cred_2sig->SetFillColor(13);
  j_cred_99p->SetFillColor(10);
  j_cred_3sig->SetFillColor(10);

  TCanvas* c = new TCanvas("c","c",600,600);
  //c->Print("jarl.pdf]");

  c->Draw();
  c->cd();
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetLeftMargin(0.13);

  j_hist->GetYaxis()->SetRangeUser(1000.,j_hist->GetMaximum()*5.);
  c->SetLogy();

  TLegend* leg;
  if(sindcp) leg = new TLegend(0.5,0.6,1.0,0.85);
  else leg = new TLegend(0.5,0.65,1.15,0.88);
  //TLegend* leg = new TLegend(0.7,0.35,1.2,0.65);
  leg->SetFillStyle(0);
  //leg->SetLineColor(0);
  leg->SetLineColorAlpha(0,0);
  leg->SetTextFont(132);
  gStyle->SetLegendBorderSize(0);
  if(sindcp) leg->SetHeader("#splitline{T2K 2022 Credible Intervals}{flat prior on sin#delta_{CP}}");
  else leg->SetHeader("DUNE FD TDR NuFit 4.0");
  leg->AddEntry(j_hist,"1#sigma","f");
  leg->AddEntry(j_cred_1sig,"2#sigma","f");
  leg->AddEntry(j_cred_2sig,"3#sigma","f");

  j_hist->Draw();
  j_cred_1sig->Draw("same");
  j_cred_2sig->Draw("same");
  j_cred_3sig->Draw("same");
  leg->Draw("same");

  gPad->RedrawAxis();
  c->Update();
  gPad->Update();

  if(!sindcp) {
    if(hierarchy == 1) c->Print("jarl1D_NH.pdf");
    else if(hierarchy == 0) c->Print("jarl1D_both.pdf");
    else if(hierarchy == -1) c->Print("jarl1D_IH.pdf");
  }
  else {
    if(hierarchy == 1) c->Print("jarl1D_NH_flatsindcp.pdf");
    else if(hierarchy == 0) c->Print("jarl1D_both_flatsindcp.pdf");
    else if(hierarchy == -1) c->Print("jarl1D_IH_flatsindcp.pdf");
  } 

}
