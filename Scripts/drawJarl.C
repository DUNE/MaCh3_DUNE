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

double* getInterval2D(TH2D *hist, double &p68, double &p90, double &p95, double &p99 ,double &p3sig)//, TH2D &h68, TH2D &h90)
{ 
  std::cout << "getting interval" << std::endl;
  TH2D *hCopy = (TH2D*)hist->Clone("hCopy");
  TH2D *hCont68 = (TH2D*)hist->Clone("hCont68");
  TH2D *hCont90 = (TH2D*)hist->Clone("hCont90");
  TH2D *hCont95 = (TH2D*)hist->Clone("hCont95");
  TH2D *hCont99 = (TH2D*)hist->Clone("hCont99");
  TH2D *hCont3sig = (TH2D*)hist->Clone("hCont3sig");
  
  double integral = hCopy->Integral();
  double tsum = 0; 
  double cont68lev = 0;// array('d', [0.0])
  double cont90lev = 0;//array('d', [0.0])
  double cont95lev = 0;//array('d', [0.0])
  double cont99lev = 0;//array('d', [0.0])
  double cont3siglev = 0;//array('d', [0.0])
  
  std::cout << integral << " " << tsum << std::endl;
  
  while ((tsum / float(integral)) < 0.9973)
    { 
      double tmax = hCopy->GetMaximum();
      tsum = tsum + tmax;
      int bin = hCopy->GetMaximumBin();
      if (tsum / float(integral) < 0.68)
        { 
          cont68lev = tmax;
          hCopy->SetBinContent(bin, -1.0);
        }
      if ((tsum / float(integral) < 0.9) && (tsum / float(integral) > 0.68))
        { 
          cont90lev = tmax;
          hCopy->SetBinContent(bin, -3.0);
        }
      if ((tsum / float(integral) < 0.955) && (tsum / float(integral) > 0.9))
        { 
          cont95lev = tmax;
          hCopy->SetBinContent(bin, -5.0);
        }
      if ((tsum / float(integral) < 0.99) && (tsum / float(integral) > 0.955))
        { 
          cont99lev = tmax;
          hCopy->SetBinContent(bin, -7.0);
        }       
      if ((tsum / float(integral) < 0.9973) && (tsum / float(integral) > 0.99))
        { 
          cont3siglev = tmax;
          hCopy->SetBinContent(bin, -9.0);
        }
    
    }
  
  double quant[5];
  quant[0] = cont90lev;
  quant[1] = cont68lev;
  
  quant[2] = cont95lev;
  quant[3] = cont99lev;
  quant[4] = cont3siglev;
  
  p90 = cont90lev;
  p68 = cont68lev;
  
  p95 = cont95lev;
  p99 = cont99lev;
  p3sig = cont3siglev;

  
  std::cout << "p68 = " << p68 << ", p90 = " << p90 << std::endl;
  std::cout << "p95 = " << p95 << ", p99 = " << p99 << std::endl;
  std::cout << "p3sig = " << p3sig << std::endl;
  
  return quant;//, hCont68, hCont90;
}

// sindcp is a boolean indicating whether we want to reweight to put a flat prior on sin(dcp)
// instead of the usual flat in dcp
// hieararchy = 1-->NH, 0-->both, -1-->IH
void drawJarl(int hierarchy, bool sindcp = false)
{

  TFile* f = new TFile("/vols/t2k/users/ljw20/software/MaCh3_DUNE_clean/MaCh3_DUNE/jarlskog_DUNE_FD_nufit.root");
  TH2D* j_th23;
  TH2D* j_dcp;
  if(!sindcp) {
    if(hierarchy == 1) {
      j_th23 = (TH2D*)f->Get("jarlskog_th23_NH");
      j_dcp  = (TH2D*)f->Get("jarlskog_dcp_NH");
    }
    else if(hierarchy == 0) {
      j_th23 = (TH2D*)f->Get("jarlskog_th23_both");
      j_dcp  = (TH2D*)f->Get("jarlskog_dcp_both");
    }
    else if(hierarchy == -1) {
      j_th23 = (TH2D*)f->Get("jarlskog_th23_IH");
      j_dcp  = (TH2D*)f->Get("jarlskog_dcp_IH");
    }
    else {
     std::cout << "Error: invalid hierarchy option. 1 for NH, 0 for both, -1 for IH" <<std::endl;
     throw;
    }
  }
  else {
    if(hierarchy == 1) j_th23 = (TH2D*)f->Get("jarlskog_th23_NH_flatsindcp");
    else if(hierarchy == 0) j_th23 = (TH2D*)f->Get("jarlskog_th23_both_flatsindcp");
    else if(hierarchy == -1) j_th23 = (TH2D*)f->Get("jarlskog_th23_IH_flatsindcp");
    else {
     std::cout << "Error: invalid hierarchy option. 1 for NH, 0 for both, -1 for IH" <<std::endl;
     throw;
    }
  }

  j_th23->SetTitle("Jarlskog Invariant;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;sin^{2}#theta_{23}");
  j_dcp->SetTitle("Jarlskog Invariant;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;#delta_{CP}");

  //if(hierarchy == -1)j_th23->Rebin2D(10,11);
  //else j_th23->Rebin2D(12,12);
 
  j_th23->Rebin2D(8,8);
  j_dcp->Rebin2D(8,8);

  TH2D* j_th23_cont_1sig = (TH2D*)j_th23->Clone("j_th23_cont_1sig");
  TH2D* j_th23_cont_2sig = (TH2D*)j_th23->Clone("j_th23_cont_2sig");
  TH2D* j_th23_cont_3sig = (TH2D*)j_th23->Clone("j_th23_cont_3sig");
  double min;
  double p68, p90, p95, p99, p3sig;
  double *p = getInterval2D(j_th23, p68, p90, p95, p99, p3sig);

  double tpp68[1];
  double tpp90[1];
  double tpp95[1];
  double tpp99[1];
  double tpp3sig[1];

  tpp68[0] = p68;
  tpp90[0] = p90;
  tpp95[0] = p95;
  tpp99[0] = p99;
  tpp3sig[0] = p3sig;

  j_th23->SetContour(255);

  j_th23_cont_1sig->Smooth(1);
  j_th23_cont_2sig->Smooth(1);
  j_th23_cont_3sig->Smooth(1);
  j_th23->Smooth(1);

  j_th23_cont_1sig->SetContour(1,tpp68);
  j_th23_cont_2sig->SetContour(1,tpp95);
  j_th23_cont_3sig->SetContour(1,tpp3sig);

  j_th23_cont_1sig->SetLineColor(kRed);
  j_th23_cont_2sig->SetLineColor(kBlue);
  j_th23_cont_3sig->SetLineColor(kGreen+3);
  
  j_th23_cont_1sig->SetLineWidth(2);
  j_th23_cont_2sig->SetLineWidth(2);
  j_th23_cont_3sig->SetLineWidth(2);

  // make J-dcp contours
  TH2D* j_dcp_cont_1sig = (TH2D*)j_dcp->Clone("j_dcp_cont_1sig");
  TH2D* j_dcp_cont_2sig = (TH2D*)j_dcp->Clone("j_dcp_cont_2sig");
  TH2D* j_dcp_cont_3sig = (TH2D*)j_dcp->Clone("j_dcp_cont_3sig");
  double min_dcp;
  double p68_dcp, p90_dcp, p95_dcp, p99_dcp, p3sig_dcp;
  double *p_dcp = getInterval2D(j_dcp, p68_dcp, p90_dcp, p95_dcp, p99_dcp, p3sig_dcp);

  double tpp68_dcp[1];
  double tpp90_dcp[1];
  double tpp95_dcp[1];
  double tpp99_dcp[1];
  double tpp3sig_dcp[1];

  tpp68_dcp[0] = p68_dcp;
  tpp90_dcp[0] = p90_dcp;
  tpp95_dcp[0] = p95_dcp;
  tpp99_dcp[0] = p99_dcp;
  tpp3sig_dcp[0] = p3sig_dcp;

  j_dcp->SetContour(255);

  j_dcp_cont_1sig->Smooth(1);
  j_dcp_cont_2sig->Smooth(1);
  j_dcp_cont_3sig->Smooth(1);
  j_dcp->Smooth(1);

  j_dcp_cont_1sig->SetContour(1,tpp68_dcp);
  j_dcp_cont_2sig->SetContour(1,tpp95_dcp);
  j_dcp_cont_3sig->SetContour(1,tpp3sig_dcp);
  j_dcp_cont_1sig->SetLineColor(kRed);
  j_dcp_cont_2sig->SetLineColor(kBlue);
  j_dcp_cont_3sig->SetLineColor(kGreen+3);
  j_dcp_cont_1sig->SetLineWidth(2);
  j_dcp_cont_2sig->SetLineWidth(2);
  j_dcp_cont_3sig->SetLineWidth(2);

  TCanvas* c = new TCanvas("c","c",600,600);
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  //c->Print("jarl.pdf]");

  c->Draw();
  c->cd();
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.13);
  gPad->SetTickx();
  gPad->SetTicky();
  TLegend* leg;

  if(sindcp) leg = new TLegend(0.6,0.15,0.95,0.4);
  else leg = new TLegend(0.55,0.15,0.97,0.4);
  //TLegend* leg = new TLegend(0.7,0.35,1.2,0.65);
  leg->SetFillStyle(0);
  //leg->SetLineColor(0);
  leg->SetLineColorAlpha(0,0);
  leg->SetTextFont(132);
  gStyle->SetLegendBorderSize(0);
  if(sindcp) leg->SetHeader("#splitline{T2K 2022, Bayesian}{flat prior on sin#delta_{CP}}");
  else leg->SetHeader("T2K 2022, Bayesian");
  leg->AddEntry(j_th23_cont_1sig,"1#sigma","l");
  leg->AddEntry(j_th23_cont_2sig,"2#sigma","l");
  leg->AddEntry(j_th23_cont_3sig,"3#sigma","l");

  // dummy th2 to extend the y-axis
  TH2D* dummy = new TH2D("dummy","dummy",10,-0.05,0.05,10,0.3,0.7);
  if(hierarchy == 1)dummy->SetTitle("Jarlskog Invariant, Normal Hierarchy;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;sin^{2}#theta_{23}");
  //else if(hierarchy == 0)dummy->SetTitle("Jarlskog Invariant, Both Hierarchies;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;sin^{2}#theta_{23}");
  else if(hierarchy == 0)dummy->SetTitle(";J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;sin^{2}#theta_{23}");
  else if(hierarchy == -1)dummy->SetTitle("Jarlskog Invariant, Inverted Hierarchy;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;sin^{2}#theta_{23}");
  dummy->GetXaxis()->SetNdivisions(508);

  //j_th23_cont_3sig->GetYaxis()->SetRangeUser(0.3,0.8);
  dummy->Draw();
  //j_th23->Draw("colz");
  j_th23_cont_3sig->Draw("same cont3");
  j_th23_cont_2sig->Draw("same cont3");
  j_th23_cont_1sig->Draw("same cont3");
  leg->Draw("same");

  c1->Draw();
  c1->cd();
  gStyle->SetOptStat(0);
  gPad->SetLeftMargin(0.13);
  gPad->SetTickx();
  gPad->SetTicky();

  TH2D* dummy_dcp = new TH2D("dummy_dcp","dummy_dcp",10,-0.05,0.05,10,-1.*TMath::Pi(),TMath::Pi());
  if(hierarchy == 1)dummy_dcp->SetTitle("Jarlskog Invariant, Normal Hierarchy;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;#delta_{CP}");
  //else if(hierarchy == 0)dummy_dcp->SetTitle("Jarlskog Invariant, Both Hierarchies;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;#delta_{CP}");
  else if(hierarchy == 0)dummy_dcp->SetTitle(";J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;#delta_{CP}");
  else if(hierarchy == -1)dummy_dcp->SetTitle("Jarlskog Invariant, Inverted Hierarchy;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;#delta_{CP}");
  dummy_dcp->GetXaxis()->SetNdivisions(508);

  dummy_dcp->Draw();
  //j_th23->Draw("colz");
  j_dcp_cont_3sig->Draw("same cont3");
  j_dcp_cont_2sig->Draw("same cont3");
  j_dcp_cont_1sig->Draw("same cont3");
  leg->Draw("same");

  if(!sindcp) {
    if(hierarchy == 1) {
      c->Print("jarl_NH.pdf");
      c1->Print("jarldcp_NH.pdf");
    }
    else if(hierarchy == 0) {
      c->Print("jarl_both.pdf");
      c1->Print("jarldcp_both.pdf");
    }
    else if(hierarchy == -1) {
      c->Print("jarl_IH.pdf");
      c1->Print("jarldcp_IH.pdf");
    }
  }
  else {
    if(hierarchy == 1) c->Print("jarl_NH_flatsindcp.pdf");
    else if(hierarchy == 0) c->Print("jarl_both_flatsindcp.pdf");
    else if(hierarchy == -1) c->Print("jarl_IH_flatsindcp.pdf");
  } 

}
