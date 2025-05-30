//////////////////////// Macro by Kevin Wood /////////////////////////////////////
//  Need to make directory 'Contours2D' for script to write output to //////////////
//  Makes all contours (woRC and wRC) and percent/sigma intervals   //////////////
//////////////////////////////////////////////////////////////////////////////////

#include "TPaveText.h"
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
#include "TPad.h"
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


// hieararchy = 1-->NH, 0-->both, -1-->IH
// app = true for appearance contours, false for disappearance
// levOpt = 0 --> 1,2,3sigma; levOpt = 1 --> 68%,90%,99%
void makePlot(TString ReducedChain, int hierarchy, int *burnin, int N_burnin, bool app=true, bool th23dcp=false, bool RC=false, int levOpt=0) {

  TCanvas* c = new TCanvas("c","c",600,500);
  c->Draw();
  c->cd();
  

  if(levOpt==1) {
    if(RC) {
      if(app && th23dcp) {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth23_NH_wRC_perc.gif[");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth23_both_wRC_perc.gif[");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth23_IH_wRC_perc.gif[");
      }
      else if (app == true && th23dcp == false){
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth13_NH_wRC_perc.gif[");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth13_both_wRC_perc.gif[");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth13_IH_wRC_perc.gif[");
      }
      else {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dm2th23_NH_wRC_perc.gif[");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dm2th23_both_wRC_perc.gif[");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dm2th23_IH_wRC_perc.gif[");
      }
    }
    else {
      if(app && th23dcp) {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth23_NH_perc.gif[");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth23_both_perc.gif[");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth23_IH_perc.gif[");
      }
      else if (app == true && th23dcp == false ){
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth13_NH_perc.gif[");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth13_both_perc.gif[");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth13_IH_perc.gif[");
      }
      else {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dm2th23_NH_perc.gif[");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dm2th23_both_perc.gif[");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dm2th23_IH_perc.gif[");
      }
    }
  }
  
  else {
    if(RC) {
      if(app && th23dcp) {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth23_NH_wRC_sigs.gif[");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth23_both_wRC_sigs.gif[");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth23_IH_wRC_sigs.gif[");
      }
      else if (app == true && th23dcp == false){
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth13_NH_wRC_sigs.gif[");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth13_both_wRC_sigs.gif[");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth13_IH_wRC_sigs.gif[");
      }
      else {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dm2th23_NH_wRC_sigs.gif[");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dm2th23_both_wRC_sigs.gif[");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dm2th23_IH_wRC_sigs.gif[");
      }
    }
    else {
      if(app && th23dcp) {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth23_NH_sigs.gif[");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth23_both_sigs.gif[");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth23_IH_sigs.gif[");
      }
      else if (app == true && th23dcp == false){
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth13_NH_sigs.gif[");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth13_both_sigs.gif[");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth13_IH_sigs.gif[");
      }
      else {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dm2th23_NH_sigs.gif[");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dm2th23_both_sigs.gif[");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dm2th23_IH_sigs.gif[");
      }
    }
  }

  for(int burn=0; burn < N_burnin; burn++) {
  
  
    bool drawAsimovPoint = true;
    double asim_dcp = -1.601;
    double asim_th13 = 0.0218;
    double asim_th23 = 0.528;
    double asim_dm2 = 2.509e-3;
  
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTextFont(42);
    gStyle->SetTextSize(0.08);
  
    TGraph* g_asim_dcp_th13 = new TGraph(1);
    g_asim_dcp_th13->SetPoint(0,asim_th13,asim_dcp);
    TGraph* g_asim_dm2_th23 = new TGraph(1);
    g_asim_dm2_th23->SetPoint(0,asim_th23,asim_dm2);
    TGraph* g_asim_dcp_th23 = new TGraph(1);
    g_asim_dcp_th23->SetPoint(0,asim_th23,asim_dcp);
    g_asim_dcp_th13->SetMarkerStyle(47);
    g_asim_dm2_th23->SetMarkerStyle(47);
    g_asim_dcp_th23->SetMarkerStyle(47);
  
    TFile* f_DUNE = new TFile(ReducedChain);
    TTree* t_DUNE = (TTree*)f_DUNE->Get("osc_posteriors");
  
    TH2D* h_dcp_th13_DUNE;
    TH2D* h_dcp_th23_DUNE;
    TH2D* h_dm32_th23_DUNE;
  
    h_dcp_th13_DUNE = new TH2D("h_dcp_th13_DUNE","DUNE FD;sin^{2}#theta_{13};#delta_{CP}",80,0.01,0.05,100,-1.*TMath::Pi(),TMath::Pi());
    h_dcp_th23_DUNE = new TH2D("h_dcp_th23_DUNE","DUNE FD;sin^{2}#theta_{23};#delta_{CP}",50,0.4,0.6,100,-1.*TMath::Pi(),TMath::Pi());
    h_dm32_th23_DUNE = new TH2D("h_dm32_th23_DUNE","DUNE FD;sin^{2}#theta_{23};#Delta m^{2}_{32}",75,0.4,0.6,400,-0.003,0.003);
  
    h_dcp_th13_DUNE->GetXaxis()->SetTitleFont(132);
    h_dcp_th13_DUNE->GetXaxis()->SetTitleOffset(0.9);
    h_dcp_th13_DUNE->GetXaxis()->SetTitleSize(0.07);
    h_dcp_th13_DUNE->GetXaxis()->SetLabelFont(132);
    h_dcp_th13_DUNE->GetXaxis()->SetLabelSize(0.07);
    TGaxis::SetMaxDigits(3) ;
    h_dcp_th13_DUNE->GetYaxis()->SetTitleFont(132);
    h_dcp_th13_DUNE->GetYaxis()->SetTitleOffset(0.9);
    h_dcp_th13_DUNE->GetYaxis()->SetTitleSize(0.07);
    h_dcp_th13_DUNE->GetYaxis()->SetLabelFont(132);
    h_dcp_th13_DUNE->GetYaxis()->SetLabelSize(0.07); 
  
    h_dcp_th23_DUNE->GetXaxis()->SetTitleFont(132);
    h_dcp_th23_DUNE->GetXaxis()->SetTitleOffset(0.9);
    h_dcp_th23_DUNE->GetXaxis()->SetTitleSize(0.07);
    h_dcp_th23_DUNE->GetXaxis()->SetLabelFont(132);
    h_dcp_th23_DUNE->GetXaxis()->SetLabelSize(0.07);
    TGaxis::SetMaxDigits(3) ;
    h_dcp_th23_DUNE->GetYaxis()->SetTitleFont(132);
    h_dcp_th23_DUNE->GetYaxis()->SetTitleOffset(0.9);
    h_dcp_th23_DUNE->GetYaxis()->SetTitleSize(0.07);
    h_dcp_th23_DUNE->GetYaxis()->SetLabelFont(132);
    h_dcp_th23_DUNE->GetYaxis()->SetLabelSize(0.07); 
  
    h_dm32_th23_DUNE->GetXaxis()->SetTitleFont(132);
    h_dm32_th23_DUNE->GetXaxis()->SetTitleOffset(0.9);
    h_dm32_th23_DUNE->GetXaxis()->SetTitleSize(0.07);
    h_dm32_th23_DUNE->GetXaxis()->SetLabelFont(132);
    h_dm32_th23_DUNE->GetXaxis()->SetLabelSize(0.07);
    TGaxis::SetMaxDigits(3) ;
    h_dm32_th23_DUNE->GetYaxis()->SetTitleFont(132);
    h_dm32_th23_DUNE->GetYaxis()->SetTitleOffset(0.95);
    h_dm32_th23_DUNE->GetYaxis()->SetTitleSize(0.07);
    h_dm32_th23_DUNE->GetYaxis()->SetLabelFont(132);
    h_dm32_th23_DUNE->GetYaxis()->SetLabelSize(0.07); 
  
    TPaveText* hlab = new TPaveText(0.72,0.83,0.855,0.88,"NDC");
  
      if(hierarchy==1) { // normal hiearchy --> dm32 > 0.
        t_DUNE->Draw(("dcp:theta13>>h_dcp_th13_DUNE","(dm23>0.)*(step>"+std::to_string(burnin[burn])+")").c_str());
        t_DUNE->Draw(("dcp:theta23>>h_dcp_th23_DUNE","(dm23>0.)*(step>"+std::to_string(burnin[burn])+")").c_str());
        t_DUNE->Draw(("dm23:theta23>>h_dm32_th23_DUNE","(dm23>0.)*(step>"+std::to_string(burnin[burn])+")").c_str());
        hlab->AddText("NH only");
      }
      else if(hierarchy==0) {
        t_DUNE->Draw(("dcp:theta13>>h_dcp_th13_DUNE","step>"+std::to_string(burnin[burn])).c_str());
        t_DUNE->Draw(("dcp:theta23>>h_dcp_th23_DUNE","step>"+std::to_string(burnin[burn])).c_str());
        t_DUNE->Draw(("dm23:theta23>>h_dm32_th23_DUNE","step>"+std::to_string(burnin[burn])).c_str());
        //hlab->AddText("NH+IH");
      }
      else if(hierarchy==-1) { // inverted hiearchy --> dm32 < 0.
        t_DUNE->Draw(("dcp:theta13>>h_dcp_th13_DUNE","(dm23<0.)*(step>"+std::to_string(burnin[burn])+")").c_str());
        t_DUNE->Draw(("dcp:theta23>>h_dcp_th23_DUNE","(dm23<0.)*(step>"+std::to_string(burnin[burn])+")").c_str());
        t_DUNE->Draw(("dm23:theta23>>h_dm32_th23_DUNE","(dm23<0.)*(step>"+std::to_string(burnin[burn])+")").c_str());
        hlab->AddText("IH only");
      }
  
    //h_dm32_th23_DUNE->GetYaxis()->SetMaxDigits(2);
  
    // DUNE-only fit
    TH2D* h_dcp_th13_cont_1sig_DUNE;
    TH2D* h_dcp_th13_cont_2sig_DUNE;
    TH2D* h_dcp_th13_cont_3sig_DUNE;
  
    if(app && th23dcp) {
      h_dcp_th13_cont_1sig_DUNE = (TH2D*)h_dcp_th23_DUNE->Clone("h_dcp_th23_cont_1sig_DUNE");
      h_dcp_th13_cont_2sig_DUNE = (TH2D*)h_dcp_th23_DUNE->Clone("h_dcp_th23_cont_2sig_DUNE");
      h_dcp_th13_cont_3sig_DUNE = (TH2D*)h_dcp_th23_DUNE->Clone("h_dcp_th23_cont_3sig_DUNE");
    }
    else if (app == true && th23dcp == false){
      h_dcp_th13_cont_1sig_DUNE = (TH2D*)h_dcp_th13_DUNE->Clone("h_dcp_th23_cont_1sig_DUNE");
      h_dcp_th13_cont_2sig_DUNE = (TH2D*)h_dcp_th13_DUNE->Clone("h_dcp_th23_cont_2sig_DUNE");
      h_dcp_th13_cont_3sig_DUNE = (TH2D*)h_dcp_th13_DUNE->Clone("h_dcp_th23_cont_3sig_DUNE");
    }
    else {
      h_dcp_th13_cont_1sig_DUNE = (TH2D*)h_dm32_th23_DUNE->Clone("h_dm32_th23_cont_1sig_DUNE");
      h_dcp_th13_cont_2sig_DUNE = (TH2D*)h_dm32_th23_DUNE->Clone("h_dm32_th23_cont_2sig_DUNE");
      h_dcp_th13_cont_3sig_DUNE = (TH2D*)h_dm32_th23_DUNE->Clone("h_dm32_th23_cont_3sig_DUNE");
    }
  
    double min_DUNE;
    double p68_DUNE, p90_DUNE, p95_DUNE, p99_DUNE, p3sig_DUNE;
    double *p_DUNE;
    if(app && th23dcp) p_DUNE = getInterval2D(h_dcp_th23_DUNE, p68_DUNE, p90_DUNE, p95_DUNE, p99_DUNE, p3sig_DUNE);
    else if(app == true && th23dcp == false) p_DUNE = getInterval2D(h_dcp_th13_DUNE, p68_DUNE, p90_DUNE, p95_DUNE, p99_DUNE, p3sig_DUNE);
    else p_DUNE = getInterval2D(h_dm32_th23_DUNE, p68_DUNE, p90_DUNE, p95_DUNE, p99_DUNE, p3sig_DUNE);
  
    double tpp68_DUNE[1];
    double tpp90_DUNE[1];
    double tpp95_DUNE[1];
    double tpp99_DUNE[1];
    double tpp3sig_DUNE[1];
  
    tpp68_DUNE[0] = p68_DUNE;
    tpp90_DUNE[0] = p90_DUNE;
    tpp95_DUNE[0] = p95_DUNE;
    tpp99_DUNE[0] = p99_DUNE;
    tpp3sig_DUNE[0] = p3sig_DUNE;
  
    h_dcp_th13_DUNE->SetContour(255);
    h_dcp_th13_cont_1sig_DUNE->Smooth(1);
    h_dcp_th13_cont_2sig_DUNE->Smooth(1);
    h_dcp_th13_cont_3sig_DUNE->Smooth(1);
    h_dcp_th13_DUNE->Smooth(1);
  
    if(levOpt==1) {
      h_dcp_th13_cont_1sig_DUNE->SetContour(1,tpp68_DUNE);
      h_dcp_th13_cont_2sig_DUNE->SetContour(1,tpp90_DUNE);
      h_dcp_th13_cont_3sig_DUNE->SetContour(1,tpp99_DUNE);
    }
    else {
      h_dcp_th13_cont_1sig_DUNE->SetContour(1,tpp68_DUNE);
      h_dcp_th13_cont_2sig_DUNE->SetContour(1,tpp95_DUNE);
      h_dcp_th13_cont_3sig_DUNE->SetContour(1,tpp3sig_DUNE);
    }
    h_dcp_th13_cont_1sig_DUNE->SetLineColor(kRed);
    h_dcp_th13_cont_2sig_DUNE->SetLineColor(kRed);
    h_dcp_th13_cont_3sig_DUNE->SetLineColor(kRed);
    h_dcp_th13_cont_1sig_DUNE->SetLineWidth(2);
    h_dcp_th13_cont_2sig_DUNE->SetLineWidth(2);
    h_dcp_th13_cont_3sig_DUNE->SetLineWidth(2);
    h_dcp_th13_cont_1sig_DUNE->SetLineStyle(2);
    h_dcp_th13_cont_3sig_DUNE->SetLineStyle(3);
  
    // find best fit point
    int Mbx, Mby, Mbz, MaCh3bin;
    double Mx, My; 
    TGraph *bestfitM = new TGraph(1);
  
    if(app && th23dcp) {
      MaCh3bin = h_dcp_th23_DUNE->GetMaximumBin();
      h_dcp_th23_DUNE->GetBinXYZ(MaCh3bin,Mbx,Mby,Mbz);
      Mx = h_dcp_th23_DUNE->GetXaxis()->GetBinCenter(Mbx);
      My = h_dcp_th23_DUNE->GetYaxis()->GetBinCenter(Mby);
    }
    else if (app == true && th23dcp == false ){
      MaCh3bin = h_dcp_th13_DUNE->GetMaximumBin();
      h_dcp_th13_DUNE->GetBinXYZ(MaCh3bin,Mbx,Mby,Mbz);
      Mx = h_dcp_th13_DUNE->GetXaxis()->GetBinCenter(Mbx);
      My = h_dcp_th13_DUNE->GetYaxis()->GetBinCenter(Mby);
    }
    else {
      MaCh3bin = h_dm32_th23_DUNE->GetMaximumBin();
      h_dm32_th23_DUNE->GetBinXYZ(MaCh3bin,Mbx,Mby,Mbz);
      Mx = h_dm32_th23_DUNE->GetXaxis()->GetBinCenter(Mbx);
      My = h_dm32_th23_DUNE->GetYaxis()->GetBinCenter(Mby);
    }
  
    bestfitM->SetPoint(0,Mx,My);
    bestfitM->SetMarkerStyle(26);
    bestfitM->SetMarkerSize(1);
    bestfitM->SetMarkerColor(kBlue);
  
    if (app && th23dcp) {
  	  std::cout << "App = " << app << "   Hierarchy = " << hierarchy << "   RC = " << RC << std::endl; 
  	  std::cout << "Best fit dCP = " << My << " --- " << "Best fit th23: " << Mx << std::endl;
    }
    else if (app == true && th23dcp == false){
      std::cout << "App = " << app << "   Hierarchy = " << hierarchy << "   RC = " << RC << std::endl; 
      std::cout << "Best fit dCP = " << My << " --- " << "Best fit th13: " << Mx << std::endl;
    }
    else{
  	  std::cout << "App = " << app << "   Hierarchy = " << hierarchy << "   RC = " << RC << std::endl; 
  	  std::cout << "Best fit dm32 = " << My << " --- " << "Best fit th23: " << Mx << std::endl;
  
    }
  
    //c->Print("jarl.pdf]");
  
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.13);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.08);
    TGaxis::SetExponentOffset(0.0, 0.01, "y");
  
    TLegend* leg;
  
    hlab->SetFillStyle(0);
    hlab->SetFillStyle(0);
    hlab->SetTextFont(132);
    hlab->SetBorderSize(0);
  
    if(!app && hierarchy==0) leg = new TLegend(0.65,0.53,0.95,0.93);
    else if(!app && hierarchy!=0) leg = new TLegend(0.50,0.6,0.94,0.8);
    else leg = new TLegend(0.50,0.55,0.95,0.75);
    //TLegend* leg = new TLegend(0.7,0.35,1.2,0.65);
    leg->SetFillStyle(0);
    //leg->SetLineColor(0);
    leg->SetLineColorAlpha(0,0);
    leg->SetTextFont(132);
    leg->SetTextSize(0.04);
    gStyle->SetLegendBorderSize(0);
    leg->SetHeader("DUNE Credible Regions");
    TH2D* solid = new TH2D("solid","solid",10,-0.05,0.05,10,0.3,0.7);
    TH2D* dashed = new TH2D("dashed","dashed",10,-0.05,0.05,10,0.3,0.7);
    solid->SetLineColor(kBlack);
    dashed->SetLineColor(kBlack);
    solid->SetLineWidth(2);
    dashed->SetLineWidth(2);
    dashed->SetLineStyle(2);
    //leg->AddEntry(dashed,"1#sigma");
    //leg->AddEntry(solid,"2#sigma");
    leg->AddEntry(bestfitM,"Best Fit","p");
    if(levOpt==1) {
      leg->AddEntry(h_dcp_th13_cont_1sig_DUNE,"68%","l");
      leg->AddEntry(h_dcp_th13_cont_2sig_DUNE,"90%","l");
      leg->AddEntry(h_dcp_th13_cont_3sig_DUNE,"99%","l");
    }
    else {
      leg->AddEntry(h_dcp_th13_cont_1sig_DUNE,"1#sigma","l");
      leg->AddEntry(h_dcp_th13_cont_2sig_DUNE,"2#sigma","l");
      leg->AddEntry(h_dcp_th13_cont_3sig_DUNE,"3#sigma","l");
    }
    if(drawAsimovPoint) leg->AddEntry(g_asim_dcp_th13,"Asimov Point (NH)","p");
  
    h_dcp_th13_cont_1sig_DUNE->GetXaxis()->SetNdivisions(507);
  
    if(!app && hierarchy == 0) {
  
      TPad *p_NH = new TPad("p_NH","p_NH",0.0,0.5,1.0,1.0);
      TPad *p_IH = new TPad("p_IH","p_IH",0.0,0.0,1.0,0.5);
      p_NH->SetBottomMargin(0.01);
      p_NH->SetTopMargin(0.22);
      p_IH->SetBottomMargin(0.22);
      p_IH->SetTopMargin(0.01);
      p_NH->SetLeftMargin(0.15);
      p_IH->SetLeftMargin(0.15);
      p_NH->Draw(); p_IH->Draw();
      TH2D* dummyNH = (TH2D*)h_dm32_th23_DUNE->Clone("dummyNH");
      TH2D* dummyIH = (TH2D*)h_dm32_th23_DUNE->Clone("dummyNH");
      dummyNH->Reset(); dummyIH->Reset();
      dummyNH->GetYaxis()->SetRangeUser(0.00225,0.00275);
      dummyIH->GetYaxis()->SetRangeUser(-0.00275,-0.00225);
      dummyNH->GetXaxis()->SetLabelSize(0.00);
      dummyNH->GetYaxis()->SetLabelSize(0.011);
      dummyNH->GetXaxis()->SetTitle("");
      dummyNH->GetYaxis()->SetTitleSize(0.1);
      dummyNH->GetYaxis()->SetTitleOffset(0.55);
      dummyIH->GetXaxis()->SetTitleSize(0.1);
      dummyIH->GetXaxis()->SetTitleOffset(0.9);
      dummyIH->GetXaxis()->SetLabelSize(0.011);
      dummyIH->GetYaxis()->SetLabelSize(0.011);
      dummyIH->GetYaxis()->SetTitle("");
      dummyIH->SetTitle("");
  
      p_NH->cd();
      dummyNH->Draw();
      h_dcp_th13_cont_1sig_DUNE->Draw("cont3 same");
      h_dcp_th13_cont_2sig_DUNE->Draw("cont3 same");
      h_dcp_th13_cont_3sig_DUNE->Draw("cont3 same");
      bestfitM->Draw("SAME.P");
      if(drawAsimovPoint) g_asim_dm2_th23->Draw("same p");
       
      p_IH->cd();
      dummyIH->Draw();
      h_dcp_th13_cont_1sig_DUNE->Draw("cont3 same");
      h_dcp_th13_cont_2sig_DUNE->Draw("cont3 same");
      h_dcp_th13_cont_3sig_DUNE->Draw("cont3 same");
      bestfitM->Draw("SAME.P");
      leg->Draw("same");
    }
    else {
      if(!app && hierarchy==1) h_dcp_th13_cont_1sig_DUNE->GetYaxis()->SetRangeUser(0.0023,0.0028);
      else if(!app && hierarchy==-1) h_dcp_th13_cont_1sig_DUNE->GetYaxis()->SetRangeUser(-0.0028,-0.0023);
      h_dcp_th13_cont_1sig_DUNE->Draw("cont3 same");
      h_dcp_th13_cont_2sig_DUNE->Draw("cont3 same");
      h_dcp_th13_cont_3sig_DUNE->Draw("cont3 same");
      bestfitM->Draw("SAME.P");
      if(drawAsimovPoint && app) g_asim_dcp_th13->Draw("same p");
      if(drawAsimovPoint && th23dcp) g_asim_dcp_th23->Draw("same p");
      else if(drawAsimovPoint && !app) g_asim_dm2_th23->Draw("same p");
      //h_dcp_th13_cont_2sig_joint->Draw("cont3 same");
      //h_dcp_th13_cont_2sig_DUNE->Draw("cont3 same");
      //h_dcp_th13_cont_2sig_nova->Draw("cont3 same");
      leg->Draw("same");
      if(hierarchy!=0) hlab->Draw();
    }
  
    if(levOpt==1) {
      if(RC) {
        if(app && th23dcp) {
          if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth23_NH_wRC_perc.gif");
          else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth23_both_wRC_perc.gif");
          else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth23_IH_wRC_perc.gif");
        }
        else if (app == true && th23dcp == false){
          if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth13_NH_wRC_perc.gif");
          else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth13_both_wRC_perc.gif");
          else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth13_IH_wRC_perc.gif");
        }
        else {
          if(hierarchy == 1) c->Print("Contours2D/DUNE_dm2th23_NH_wRC_perc.gif");
          else if(hierarchy == 0) c->Print("Contours2D/DUNE_dm2th23_both_wRC_perc.gif");
          else if(hierarchy == -1) c->Print("Contours2D/DUNE_dm2th23_IH_wRC_perc.gif");
        }
      }
      else {
        if(app && th23dcp) {
          if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth23_NH_perc.gif");
          else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth23_both_perc.gif");
          else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth23_IH_perc.gif");
        }
        else if (app == true && th23dcp == false ){
          if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth13_NH_perc.gif");
          else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth13_both_perc.gif");
          else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth13_IH_perc.gif");
        }
        else {
          if(hierarchy == 1) c->Print("Contours2D/DUNE_dm2th23_NH_perc.gif");
          else if(hierarchy == 0) c->Print("Contours2D/DUNE_dm2th23_both_perc.gif");
          else if(hierarchy == -1) c->Print("Contours2D/DUNE_dm2th23_IH_perc.gif");
        }
      }
    }
  
    else {
      if(RC) {
        if(app && th23dcp) {
          if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth23_NH_wRC_sigs.gif");
          else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth23_both_wRC_sigs.gif");
          else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth23_IH_wRC_sigs.gif");
        }
        else if (app == true && th23dcp == false){
          if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth13_NH_wRC_sigs.gif");
          else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth13_both_wRC_sigs.gif");
          else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth13_IH_wRC_sigs.gif");
        }
        else {
          if(hierarchy == 1) c->Print("Contours2D/DUNE_dm2th23_NH_wRC_sigs.gif");
          else if(hierarchy == 0) c->Print("Contours2D/DUNE_dm2th23_both_wRC_sigs.gif");
          else if(hierarchy == -1) c->Print("Contours2D/DUNE_dm2th23_IH_wRC_sigs.gif");
        }
      }
      else {
        if(app && th23dcp) {
          if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth23_NH_sigs.gif");
          else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth23_both_sigs.gif");
          else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth23_IH_sigs.gif");
        }
        else if (app == true && th23dcp == false){
          if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth13_NH_sigs.gif");
          else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth13_both_sigs.gif");
          else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth13_IH_sigs.gif");
        }
        else {
          if(hierarchy == 1) c->Print("Contours2D/DUNE_dm2th23_NH_sigs.gif");
          else if(hierarchy == 0) c->Print("Contours2D/DUNE_dm2th23_both_sigs.gif");
          else if(hierarchy == -1) c->Print("Contours2D/DUNE_dm2th23_IH_sigs.gif");
        }
      }
    }
  
    delete t_DUNE;
    delete f_DUNE;

  }

  if(levOpt==1) {
    if(RC) {
      if(app && th23dcp) {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth23_NH_wRC_perc.gif]");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth23_both_wRC_perc.gif]");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth23_IH_wRC_perc.gif]");
      }
      else if (app == true && th23dcp == false){
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth13_NH_wRC_perc.gif]");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth13_both_wRC_perc.gif]");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth13_IH_wRC_perc.gif]");
      }
      else {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dm2th23_NH_wRC_perc.gif]");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dm2th23_both_wRC_perc.gif]");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dm2th23_IH_wRC_perc.gif]");
      }
    }
    else {
      if(app && th23dcp) {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth23_NH_perc.gif]");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth23_both_perc.gif]");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth23_IH_perc.gif]");
      }
      else if (app == true && th23dcp == false ){
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth13_NH_perc.gif]");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth13_both_perc.gif]");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth13_IH_perc.gif]");
      }
      else {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dm2th23_NH_perc.gif]");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dm2th23_both_perc.gif]");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dm2th23_IH_perc.gif]");
      }
    }
  }
  
  else {
    if(RC) {
      if(app && th23dcp) {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth23_NH_wRC_sigs.gif]");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth23_both_wRC_sigs.gif]");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth23_IH_wRC_sigs.gif]");
      }
      else if (app == true && th23dcp == false){
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth13_NH_wRC_sigs.gif]");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth13_both_wRC_sigs.gif]");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth13_IH_wRC_sigs.gif]");
      }
      else {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dm2th23_NH_wRC_sigs.gif]");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dm2th23_both_wRC_sigs.gif]");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dm2th23_IH_wRC_sigs.gif]");
      }
    }
    else {
      if(app && th23dcp) {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth23_NH_sigs.gif]");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth23_both_sigs.gif]");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth23_IH_sigs.gif]");
      }
      else if (app == true && th23dcp == false){
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dcpth13_NH_sigs.gif]");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dcpth13_both_sigs.gif]");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dcpth13_IH_sigs.gif]");
      }
      else {
        if(hierarchy == 1) c->Print("Contours2D/DUNE_dm2th23_NH_sigs.gif]");
        else if(hierarchy == 0) c->Print("Contours2D/DUNE_dm2th23_both_sigs.gif]");
        else if(hierarchy == -1) c->Print("Contours2D/DUNE_dm2th23_IH_sigs.gif]");
      }
    }
  }

  delete c;
}

void make2DContoursGIF(TString ReducedChain) {
//void makePlot(TString ReducedChain, int hierarchy,bool app=true, bool th23dcp=false, bool RC=false, int levOpt=0, int *burnin=80000, int N_burnin=1) {
  //Appearance = th13_cp
  //makePlot(ReducedChain,0,true,false,1);
  int N_burnin = 7;
  int burninList[7] = {100, 1000, 10000, 50000, 80000, 100000, 500000};
  makePlot(ReducedChain,0, burninList, N_burnin, true, false);

  //Dissappearance = th23_dm32
  //makePlot(ReducedChain,0,false,false,1);
  makePlot(ReducedChain,1, burninList, N_burnin, false, false);
  
  // th23dCP plots
  //makePlot(ReducedChain,0,true,true,false,1);
  makePlot(ReducedChain,0, burninList, N_burnin, true, true);
}
