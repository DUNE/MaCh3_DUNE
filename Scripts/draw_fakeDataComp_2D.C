/*
 * ETA - this script is adapted from krw_2Dcreds.C but is now meant to be used for fake data studies.
 * This will basically do the same as the 2D script but also do it on a stat-only fit
 */

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

TString nominal_file;
TString stat_file;

TString nominal_title;
TString stat_title;

TString outdir;

TString plot_title;// = "T2K+NOvA FDS comparison";

void saveContours(TString file_name, bool RC, int app, int levOpt, int hierarchy, TH2D* h_dcp_th13_cont_1sig_asimov, TH2D* h_dcp_th13_cont_2sig_asimov, TH2D* h_dcp_th13_cont_3sig_asimov){

  TFile *out_file_asimov = new TFile(file_name, "UPDATE");
  out_file_asimov->cd();
 
  if(RC){
	if(app == 0){
	  if(levOpt == 1){
		if(hierarchy == 1){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth13_NH_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth13_NH_wRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth13_NH_wRC_99per");
		}
		else if(hierarchy == 0){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth13_both_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth13_both_wRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth13_both_wRC_99per");
		}
		else if(hierarchy == -1 ){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth13_IH_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth13_IH_wRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth13_IH_wRC_99per");
		}
	  }
	  else{
		if(hierarchy == 1){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth13_NH_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth13_NH_wRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth13_NH_wRC_3sig");
		}
		else if(hierarchy == 0){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth13_both_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth13_both_wRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth13_both_wRC_3sig");
		}
		else if(hierarchy == -1 ){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth13_IH_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth13_IH_wRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth13_IH_wRC_3sig");
		}
	  }
	}
	else if(app == 1){
	  if(levOpt == 1){
		if(hierarchy == 1){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dm2th23_NH_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dm2th23_NH_wRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dm2th23_NH_wRC_99per");
		}
		else if(hierarchy == 0){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dm2th23_both_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dm2th23_both_wRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dm2th23_both_wRC_99per");
		}
		else if(hierarchy == -1 ){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dm2th23_IH_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dm2th23_IH_wRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dm2th23_IH_wRC_99per");
		}
	  }	
	  else{
		if(hierarchy == 1){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dm2th23_NH_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dm2th23_NH_wRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dm2th23_NH_wRC_3sig");
		}
		else if(hierarchy == 0){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dm2th23_both_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dm2th23_both_wRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dm2th23_both_wRC_3sig");
		}
		else if(hierarchy == -1 ){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dm2th23_IH_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dm2th23_IH_wRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dm2th23_IH_wRC_3sig");
		}
	  }
	}
	else if(app == 2){
	  if(levOpt == 1){
		if(hierarchy == 1){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth23_NH_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth23_NH_wRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth23_NH_wRC_99per");
		}
		else if(hierarchy == 0){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth23_both_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth23_both_wRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth23_both_wRC_99per");
		}
		else if(hierarchy == -1 ){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth23_IH_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth23_IH_wRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth23_IH_wRC_99per");
		}
	  }	
	  else{
		if(hierarchy == 1){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth23_NH_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth23_NH_wRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth23_NH_wRC_3sig");
		}
		else if(hierarchy == 0){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth23_both_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth23_both_wRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth23_both_wRC_3sig");
		}
		else if(hierarchy == -1 ){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth23_IH_wRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth23_IH_wRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth23_IH_wRC_3sig");
		}
	  }
	}
  }
  else{
	if(app == 0){
	  if(levOpt == 1){
		if(hierarchy == 1){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth13_NH_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth13_NH_woRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth13_NH_woRC_99per");
		}
		else if(hierarchy == 0){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth13_both_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth13_both_woRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth13_both_woRC_99per");
		}
		else if(hierarchy == -1 ){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth13_IH_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth13_IH_woRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth13_IH_woRC_99per");
		}
	  }
	  else{
		if(hierarchy == 1){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth13_NH_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth13_NH_woRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth13_NH_woRC_3sig");
		}
		else if(hierarchy == 0){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth13_both_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth13_both_woRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth13_both_woRC_3sig");
		}
		else if(hierarchy == -1 ){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth13_IH_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth13_IH_woRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth13_IH_woRC_3sig");
		}
	  }
	}
	else if(app == 1){
	  if(levOpt == 1){
		if(hierarchy == 1){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dm2th23_NH_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dm2th23_NH_woRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dm2th23_NH_woRC_99per");
		}
		else if(hierarchy == 0){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dm2th23_both_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dm2th23_both_woRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dm2th23_both_woRC_99per");
		}
		else if(hierarchy == -1 ){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dm2th23_IH_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dm2th23_IH_woRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dm2th23_IH_woRC_99per");
		}
	  }	
	  else{
		if(hierarchy == 1){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dm2th23_NH_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dm2th23_NH_woRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dm2th23_NH_woRC_3sig");
		}
		else if(hierarchy == 0){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dm2th23_both_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dm2th23_both_woRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dm2th23_both_woRC_3sig");
		}
		else if(hierarchy == -1 ){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dm2th23_IH_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dm2th23_IH_woRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dm2th23_IH_woRC_3sig");
		}
	  }
	}
	else if(app == 2){
	  if(levOpt == 1){
		if(hierarchy == 1){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth23_NH_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth23_NH_woRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth23_NH_woRC_99per");
		}
		else if(hierarchy == 0){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth23_both_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth23_both_woRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth23_both_woRC_99per");
		}
		else if(hierarchy == -1 ){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth23_IH_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth23_IH_woRC_90per");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth23_IH_woRC_99per");
		}
	  }	
	  else{
		if(hierarchy == 1){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth23_NH_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth23_NH_woRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth23_NH_woRC_3sig");
		}
		else if(hierarchy == 0){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth23_both_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth23_both_woRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth23_both_woRC_3sig");
		}
		else if(hierarchy == -1 ){
		  h_dcp_th13_cont_1sig_asimov->Write("cont_dcpth23_IH_woRC_1sig");
		  h_dcp_th13_cont_2sig_asimov->Write("cont_dcpth23_IH_woRC_2sig");
		  h_dcp_th13_cont_3sig_asimov->Write("cont_dcpth23_IH_woRC_3sig");
		}
	  }
	}
  }

  out_file_asimov->Write();
  out_file_asimov->Close();
  delete out_file_asimov;
}

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
// app = 0 --> dcp vs. sin2th13; app = 1 --> dm32 vs. sin2th23; app = 2 --> dcp vs. sin2th23
//void makePlot(int hierarchy, int app=0 ,bool RC=false, int octant=0, int levOpt=0, bool statRC = true, bool asimovRC = true, int burn_in=150000, int burn_in_stat=200000) {

void makePlot(int hierarchy, int app=0, bool LogY = false, int octant=0, bool levSigs = true, int burn_in=80000, int burn_in_stat=80000, TString title_1="w/ Syst", TString title_2="stat-only", bool RC = false, bool asimovRC = false, bool statRC = true, int asimov_point=1){


  bool saveConts = false;

  bool drawAsimovPoint = true;

  //NuFit 4.0 NH
  double asim_dcp = -2.498;
  double asim_th13 = 0.0224;
  double asim_th23 = 0.582;
  double asim_dm2 = 2.525e-3;

  //double asim_dcp = -1.601;
  //double asim_th13 = 0.0218;
  //double asim_th23 = 0.528;
  //double asim_dm2 = 2.509e-3;

  if(asimov_point == 4){
	asim_dcp = -1.601;
	asim_th13 = 0.0218;
	asim_th23 = 0.55;
	asim_dm2 = -2.45e-3;
  }


  //Asimov 4 points
  /*

  */

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

  TFile* f_asimov = new TFile(nominal_file, "READ");
  TTree* t_asimov = 0;
  t_asimov = (TTree*)f_asimov->Get("osc_posteriors");
  if(!t_asimov){
	std::cout << "[WARNING] ~~ Couldn't find \"osc_posteriors\" tree, this could be because the chain isn't a reduced chain so will look for \"posteriors\" instead" << std::endl;
	t_asimov = (TTree*)f_asimov->Get("posteriors");

	//If the TTree still can't be found then this is actually a problem
	if(!t_asimov){
      std::cout << "[ERROR] ~~ Couldn't find \"posteriors\" TTree either, this is actually a problem now!" << std::endl;
	  throw;
	}
  }

  //Now do the same for the stat only file
  TFile* f_stat = new TFile(stat_file, "READ");
  TTree* t_stat = 0;
  t_stat = (TTree*)f_stat->Get("osc_posteriors");
  if(!t_stat){
	std::cout << "[WARNING] ~~ Couldn't find \"osc_posteriors\" tree, this could be because the chain isn't a reduced chain so will look for \"posteriors\" instead" << std::endl;
	t_stat = (TTree*)f_stat->Get("posteriors");

	//If the TTree still can't be found then this is actually a problem
	if(!t_stat){
      std::cout << "[ERROR] ~~ Couldn't find \"posteriors\" TTree either, this is actually a problem now!" << std::endl;
	  throw;
	}
  }


  //Histograms to store posterior
  TH2D* h_dcp_th13_asimov;
  TH2D* h_dm32_th23_asimov;
  TH2D* h_dcp_th23_asimov;

  // histogram from stat only
  TH2D* h_dcp_th13_stat;
  TH2D* h_dm32_th23_stat;
  TH2D* h_dcp_th23_stat;

  int nbins_dm32 = 300;//375;
  int nbins_th23 = 80;//40
  int nbins_th13 = 60;//48
  int nbins_dcp = 40;//90
/*

  if(RC) {
    h_dcp_th13_asimov = new TH2D("h_dcp_th13_asimov", Form("%s, With RC;sin^{2}#theta_{13};#delta_{CP}", plot_title.Data()), nbins_th13,0.018,0.03, nbins_dcp,-1.*TMath::Pi(),TMath::Pi());
    h_dm32_th23_asimov = new TH2D("h_dm32_th23_asimov", Form("%s, With RC;sin^{2}#theta_{23};#Delta m^{2}_{32}", plot_title.Data()),nbins_th23, 0.4, 0.8, nbins_dm32,-0.003,0.003);
    h_dcp_th23_asimov = new TH2D("h_dcp_th23_asimov", Form("%s, With RC;sin^{2}#theta_{23};#delta_{CP}", plot_title.Data()), nbins_th23, 0.4, 0.8, nbins_dcp, -1.*TMath::Pi(), TMath::Pi());

    h_dcp_th13_stat = new TH2D("h_dcp_th13_stat", Form("%s, With RC;sin^{2}#theta_{13};#delta_{CP}", plot_title.Data()), nbins_th13, 0.018,0.03, nbins_dcp, -1.*TMath::Pi(),TMath::Pi());
    h_dm32_th23_stat = new TH2D("h_dm32_th23_stat", Form("%s, With RC;sin^{2}#theta_{23};#Delta m^{2}_{32}", plot_title.Data()), nbins_th23, 0.4, 0.8, nbins_dm32,-0.003,0.003);
    h_dcp_th23_stat = new TH2D("h_dcp_th23_stat", Form("%s, With RC;sin^{2}#theta_{23};#delta_{CP}", plot_title.Data()), nbins_th23,0.4,0.8, nbins_dcp, -1.*TMath::Pi(), TMath::Pi());
  }
  else {
    h_dcp_th13_asimov = new TH2D("h_dcp_th13_asimov", Form("%s;sin^{2}#theta_{13};#delta_{CP}", plot_title.Data()), nbins_th13, 0.01,0.08,75,-1.*TMath::Pi(),TMath::Pi());
    h_dm32_th23_asimov = new TH2D("h_dm32_th23_asimov", Form("%s;sin^{2}#theta_{23};#Delta m^{2}_{32}", plot_title.Data()),75,0.4,0.72, nbins_dm32,-0.003,0.003);
    h_dcp_th23_asimov = new TH2D("h_dcp_th23_asimov", Form("%s;sin^{2}#theta_{23};#delta_{CP}", plot_title.Data()),75,0.4,0.8, nbins_dcp, -1.*TMath::Pi(), TMath::Pi());

    h_dcp_th13_stat = new TH2D("h_dcp_th13_stat", Form("%s;sin^{2}#theta_{13};#delta_{CP}", plot_title.Data()), nbins_th13,0.01,0.08,75,-1.*TMath::Pi(),TMath::Pi());
    h_dm32_th23_stat = new TH2D("h_dm32_th23_stat", Form("%s;sin^{2}#theta_{23};#Delta m^{2}_{32}", plot_title.Data()),75,0.4,0.72,nbins_dm32,-0.003,0.003);
    h_dcp_th23_stat = new TH2D("h_dcp_th23_stat", Form("%s;sin^{2}#theta_{23};#delta_{CP}", plot_title.Data()),75,0.4,0.8, nbins_dcp, -1.*TMath::Pi(), TMath::Pi());
  }
*/

// Changed Binning for 1sig only

  if(RC) {
    h_dcp_th13_asimov = new TH2D("h_dcp_th13_asimov", Form("%s, With RC;sin^{2}#theta_{13};#delta_{CP}", plot_title.Data()), nbins_th13,0.018,0.03, nbins_dcp,-1.*TMath::Pi(),0);
    h_dm32_th23_asimov = new TH2D("h_dm32_th23_asimov", Form("%s, With RC;sin^{2}#theta_{23};#Delta m^{2}_{32}", plot_title.Data()),nbins_th23, 0.4, 0.8, nbins_dm32,-0.003,0.003);
    h_dcp_th23_asimov = new TH2D("h_dcp_th23_asimov", Form("%s, With RC;sin^{2}#theta_{23};#delta_{CP}", plot_title.Data()), nbins_th23, 0.4, 0.6, nbins_dcp, -1.*TMath::Pi(), 0);

    h_dcp_th13_stat = new TH2D("h_dcp_th13_stat", Form("%s, With RC;sin^{2}#theta_{13};#delta_{CP}", plot_title.Data()), nbins_th13, 0.018,0.03, nbins_dcp, -1.*TMath::Pi(),TMath::Pi());
    h_dm32_th23_stat = new TH2D("h_dm32_th23_stat", Form("%s, With RC;sin^{2}#theta_{23};#Delta m^{2}_{32}", plot_title.Data()), nbins_th23, 0.4, 0.8, nbins_dm32,-0.003,0.003);
    h_dcp_th23_stat = new TH2D("h_dcp_th23_stat", Form("%s, With RC;sin^{2}#theta_{23};#delta_{CP}", plot_title.Data()), nbins_th23,0.4,0.8, nbins_dcp, -1.*TMath::Pi(), TMath::Pi());
  }
  else {
    h_dcp_th13_asimov = new TH2D("h_dcp_th13_asimov", Form("%s;sin^{2}#theta_{13};#delta_{CP}", plot_title.Data()), nbins_th13, 0.01,0.05,nbins_dcp,-1.*TMath::Pi(),TMath::Pi());
    h_dm32_th23_asimov = new TH2D("h_dm32_th23_asimov", Form("%s;sin^{2}#theta_{23};#Delta m^{2}_{32}", plot_title.Data()),nbins_th23,0.42,0.7, nbins_dm32,-0.003,0.003);
    h_dcp_th23_asimov = new TH2D("h_dcp_th23_asimov", Form("%s;sin^{2}#theta_{23};#delta_{CP}", plot_title.Data()),nbins_th23,0.4,0.7, nbins_dcp, -1.*TMath::Pi(),TMath::Pi() );


    h_dcp_th13_stat = new TH2D("h_dcp_th13_stat", Form("%s;sin^{2}#theta_{13};#delta_{CP}", plot_title.Data()), nbins_th13, 0.01,0.05,nbins_dcp,-1.*TMath::Pi(),TMath::Pi());
    h_dm32_th23_stat = new TH2D("h_dm32_th23_stat", Form("%s;sin^{2}#theta_{23};#Delta m^{2}_{32}", plot_title.Data()),nbins_th23,0.42,0.7, nbins_dm32,-0.003,0.003);
    h_dcp_th23_stat = new TH2D("h_dcp_th23_stat", Form("%s;sin^{2}#theta_{23};#delta_{CP}", plot_title.Data()),nbins_th23,0.4,0.7, nbins_dcp, -1.*TMath::Pi(),TMath::Pi() );

  }

  h_dcp_th13_asimov->GetXaxis()->SetTitleFont(132);
  h_dcp_th13_asimov->GetXaxis()->SetTitleOffset(0.9);
  h_dcp_th13_asimov->GetXaxis()->SetTitleSize(0.07);
  h_dcp_th13_asimov->GetXaxis()->SetLabelFont(132);
  h_dcp_th13_asimov->GetXaxis()->SetLabelSize(0.07);
  TGaxis::SetMaxDigits(3) ;
  h_dcp_th13_asimov->GetYaxis()->SetTitleFont(132);
  h_dcp_th13_asimov->GetYaxis()->SetTitleOffset(0.9);
  h_dcp_th13_asimov->GetYaxis()->SetTitleSize(0.07);
  h_dcp_th13_asimov->GetYaxis()->SetLabelFont(132);
  h_dcp_th13_asimov->GetYaxis()->SetLabelSize(0.07); 

  h_dm32_th23_asimov->GetXaxis()->SetTitleFont(132);
  h_dm32_th23_asimov->GetXaxis()->SetTitleOffset(0.9);
  h_dm32_th23_asimov->GetXaxis()->SetTitleSize(0.07);
  h_dm32_th23_asimov->GetXaxis()->SetLabelFont(132);
  h_dm32_th23_asimov->GetXaxis()->SetLabelSize(0.07);
  TGaxis::SetMaxDigits(3) ;
  h_dm32_th23_asimov->GetYaxis()->SetTitleFont(132);
  h_dm32_th23_asimov->GetYaxis()->SetTitleOffset(0.95);
  h_dm32_th23_asimov->GetYaxis()->SetTitleSize(0.07);
  h_dm32_th23_asimov->GetYaxis()->SetLabelFont(132);
  h_dm32_th23_asimov->GetYaxis()->SetLabelSize(0.05); 

  h_dcp_th23_asimov->GetXaxis()->SetTitleFont(132);
  h_dcp_th23_asimov->GetXaxis()->SetTitleOffset(0.9);
  h_dcp_th23_asimov->GetXaxis()->SetTitleSize(0.07);
  h_dcp_th23_asimov->GetXaxis()->SetLabelFont(132);
  h_dcp_th23_asimov->GetXaxis()->SetLabelSize(0.07);
  TGaxis::SetMaxDigits(3) ;
  h_dcp_th23_asimov->GetYaxis()->SetTitleFont(132);
  h_dcp_th23_asimov->GetYaxis()->SetTitleOffset(0.95);
  h_dcp_th23_asimov->GetYaxis()->SetTitleSize(0.07);
  h_dcp_th23_asimov->GetYaxis()->SetLabelFont(132);
  h_dcp_th23_asimov->GetYaxis()->SetLabelSize(0.07); 

  TPaveText* hlab = new TPaveText(0.72,0.83,0.855,0.88,"NDC");

  TString octant_cut;

  if(octant == -1){octant_cut="*(th23<0.5)";}
  if(octant == 1){octant_cut="*(th23>0.5)";}

  if(RC){
    if(hierarchy==1) { // normal hiearchy --> dm32 > 0.

	  if(asimovRC){
		t_asimov->Draw("dcp:theta13>>h_dcp_th13_asimov", Form("RCreweight*(dm23>0.)*(step>%d)%s", burn_in, octant_cut.Data()));
		t_asimov->Draw("dm23:theta23>>h_dm32_th23_asimov", Form("RCreweight*(dm23>0.)*(step>%d)%s", burn_in, octant_cut.Data()));
		t_asimov->Draw("dcp:theta23>>h_dcp_th23_asimov", Form("RCreweight*(dm23>0.)*(step>%d)%s", burn_in, octant_cut.Data()));
	  }
	  else{
		t_asimov->Draw("dcp:theta13>>h_dcp_th13_asimov", Form("(dm23>0.)*(step>%d)%s", burn_in, octant_cut.Data()));
		t_asimov->Draw("dm23:theta23>>h_dm32_th23_asimov", Form("(dm23>0.)*(step>%d)%s", burn_in, octant_cut.Data()));
		t_asimov->Draw("dcp:theta23>>h_dcp_th23_asimov", Form("(dm23>0.)*(step>%d)%s", burn_in, octant_cut.Data()));
	  }

	  //if the stat-only fit has been ran with RC then you don't want to apply a RC again!
	  if(statRC){
		t_stat->Draw("dcp:theta13>>h_dcp_th13_stat", Form("RCreweight*(dm23>0.)*(step>%d)%s", burn_in_stat, octant_cut.Data()));
		t_stat->Draw("dm23:theta23>>h_dm32_th23_stat", Form("RCreweight*(dm23>0.)*(step>%d)%s", burn_in_stat, octant_cut.Data()));
		t_stat->Draw("dcp:theta23>>h_dcp_th23_stat", Form("RCreweight*(dm23>0.)*(step>%d)%s", burn_in_stat, octant_cut.Data()));
	  }
	  else{
		t_stat->Draw("dcp:theta13>>h_dcp_th13_stat", Form("(dm23>0.)*(step>%d)", burn_in_stat));
		t_stat->Draw("dm23:theta23>>h_dm32_th23_stat", Form("(dm23>0.)*(step>%d)", burn_in_stat));
		t_stat->Draw("dcp:theta23>>h_dcp_th23_stat", Form("(dm23>0.)*(step>%d)", burn_in_stat));
	  }

      //hlab->AddText("NH only");
    }
    else if(hierarchy==0) {

	  if(asimovRC){
		t_asimov->Draw("dcp:theta13>>h_dcp_th13_asimov", Form("RCreweight*(step>%d)", burn_in));
		t_asimov->Draw("dm23:theta23>>h_dm32_th23_asimov", Form("RCreweight*(step>%d)", burn_in));
		t_asimov->Draw("dcp:theta23>>h_dcp_th23_asimov", Form("RCreweight*(step>%d)", burn_in));
	  }
	  else{
		t_asimov->Draw("dcp:theta13>>h_dcp_th13_asimov", Form("(step>%d)", burn_in));
		t_asimov->Draw("dm23:theta23>>h_dm32_th23_asimov", Form("(step>%d)", burn_in));
		t_asimov->Draw("dcp:theta23>>h_dcp_th23_asimov", Form("(step>%d)", burn_in));
	  }

	  //if the stat-only fit has been ran with RC then you don't want to apply a RC again!
	  if(statRC){
		t_stat->Draw("dcp:theta13>>h_dcp_th13_stat", Form("RCreweight*(step>%d)", burn_in_stat));
		t_stat->Draw("dm23:theta23>>h_dm32_th23_stat", Form("RCreweight*(step>%d)", burn_in_stat));
		t_stat->Draw("dcp:theta23>>h_dcp_th23_stat", Form("RCreweight*(step>%d)", burn_in_stat));
	  }
	  else{
		t_stat->Draw("dcp:theta13>>h_dcp_th13_stat", Form("(step>%d)", burn_in_stat));
		t_stat->Draw("dm23:theta23>>h_dm32_th23_stat", Form("(step>%d)", burn_in_stat));
		t_stat->Draw("dcp:theta23>>h_dcp_th23_stat", Form("(step>%d)", burn_in_stat));
	  }

      //hlab->AddText("NH+IH");
    }
    else if(hierarchy==-1) { // inverted hiearchy --> dm32 < 0.

	  if(asimovRC){
		t_asimov->Draw("dcp:theta13>>h_dcp_th13_asimov", Form("RCreweight*(dm23<0.)*(step>%d)", burn_in));
		t_asimov->Draw("dm23:theta23>>h_dm32_th23_asimov", Form("RCreweight*(dm23<0.)*(step>%d)", burn_in));
		t_asimov->Draw("dcp:theta23>>h_dcp_th23_asimov", Form("RCreweight*(dm23<0.)*(step>%d)", burn_in));
	  }
	  else{
		t_asimov->Draw("dcp:theta13>>h_dcp_th13_asimov", Form("(dm23<0.)*(step>%d)", burn_in));
		t_asimov->Draw("dm23:theta23>>h_dm32_th23_asimov", Form("(dm23<0.)*(step>%d)", burn_in));
		t_asimov->Draw("dcp:theta23>>h_dcp_th23_asimov", Form("(dm23<0.)*(step>%d)", burn_in));
	  }

	  //if the stat-only fit has been ran with RC then you don't want to apply a RC again!
	  if(statRC){
		t_stat->Draw("dcp:theta13>>h_dcp_th13_stat", Form("RCreweight*(dm23<0.)*(step>%d)", burn_in_stat));
		t_stat->Draw("dm23:theta23>>h_dm32_th23_stat", Form("RCreweight*(dm23<0.)*(step>%d)", burn_in_stat));
		t_stat->Draw("dcp:theta23>>h_dcp_th23_stat", Form("RCreweight*(dm23<0.)*(step>%d)", burn_in_stat));
	  }
	  else{
		t_stat->Draw("dcp:theta13>>h_dcp_th13_stat", Form("(dm23<0.)*(step>%d)", burn_in_stat));
		t_stat->Draw("dm23:theta23>>h_dm32_th23_stat", Form("(dm23<0.)*(step>%d)", burn_in_stat));
		t_stat->Draw("dcp:theta23>>h_dcp_th23_stat", Form("(dm23<0.)*(step>%d)", burn_in_stat));
	  }
      hlab->AddText("IH only");
    }
  }
  else {
	if(hierarchy==1) { // normal hiearchy --> dm32 > 0.

	  t_asimov->Draw("dcp:theta13>>h_dcp_th13_asimov", Form("(dm23>0.)*(step>%d)%s", burn_in, octant_cut.Data()));
	  t_asimov->Draw("dm23:theta23>>h_dm32_th23_asimov", Form("(dm23>0.)*(step>%d)%s", burn_in, octant_cut.Data()));
	  t_asimov->Draw("dcp:theta23>>h_dcp_th23_asimov", Form("(dm23>0.)*(step>%d)%s", burn_in, octant_cut.Data()));

	  t_stat->Draw("dcp:theta13>>h_dcp_th13_stat", Form("(dm23>0.)*(step>%d)", burn_in_stat));
	  t_stat->Draw("dm23:theta23>>h_dm32_th23_stat", Form("(dm23>0.)*(step>%d)", burn_in_stat));
	  t_stat->Draw("dcp:theta23>>h_dcp_th23_stat", Form("(dm23>0.)*(step>%d)", burn_in_stat)); 

	  //hlab->AddText("NH only");
	}
	else if(hierarchy==0) {

	  t_asimov->Draw("dcp:theta13>>h_dcp_th13_asimov", Form("(step>%d)", burn_in));
	  t_asimov->Draw("dm23:theta23>>h_dm32_th23_asimov", Form("(step>%d)", burn_in));
	  t_asimov->Draw("dcp:theta23>>h_dcp_th23_asimov", Form("(step>%d)", burn_in));

	  t_stat->Draw("dcp:theta13>>h_dcp_th13_stat", Form("(step>%d)", burn_in_stat));
	  t_stat->Draw("dm23:theta23>>h_dm32_th23_stat", Form("(step>%d)", burn_in_stat));
	  t_stat->Draw("dcp:theta23>>h_dcp_th23_stat", Form("(step>%d)", burn_in_stat));

	  hlab->AddText("NH+IH");
	}
	else if(hierarchy==-1) { // inverted hiearchy --> dm32 < 0.

	  t_asimov->Draw("dcp:theta13>>h_dcp_th13_asimov", Form("(dm23<0.)*(step>%d)", burn_in));
	  t_asimov->Draw("dm23:theta23>>h_dm32_th23_asimov", Form("(dm23<0.)*(step>%d)", burn_in));
	  t_asimov->Draw("dcp:theta23>>h_dcp_th23_asimov", Form("(dm23<0.)*(step>%d)", burn_in));

	  t_stat->Draw("dcp:theta13>>h_dcp_th13_stat", Form("(dm23<0.)*(step>%d)", burn_in_stat));
	  t_stat->Draw("dm23:theta23>>h_dm32_th23_stat", Form("(dm23<0.)*(step>%d)", burn_in_stat));
	  t_stat->Draw("dcp:theta23>>h_dcp_th23_stat", Form("(dm23<0.)*(step>%d)", burn_in_stat));

	  hlab->AddText("IH only");
	}

  }
   //j_th23->SetTitle("Jarlskog Invariant;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;sin^{2}#theta_{23}");

  //h_dm32_th23_asimov->GetYaxis()->SetMaxDigits(2);

  TH2D* h_dcp_th13_cont_1sig_asimov;
  TH2D* h_dcp_th13_cont_2sig_asimov;
  TH2D* h_dcp_th13_cont_3sig_asimov;
 
  //stat-only fit 
  TH2D* h_dcp_th13_cont_1sig_stat;
  TH2D* h_dcp_th13_cont_2sig_stat;
  TH2D* h_dcp_th13_cont_3sig_stat;


  if(app == 0) {
    h_dcp_th13_cont_1sig_asimov = (TH2D*)h_dcp_th13_asimov->Clone("h_dcp_th13_cont_1sig_asimov");
    h_dcp_th13_cont_2sig_asimov = (TH2D*)h_dcp_th13_asimov->Clone("h_dcp_th13_cont_2sig_asimov");
    h_dcp_th13_cont_3sig_asimov = (TH2D*)h_dcp_th13_asimov->Clone("h_dcp_th13_cont_3sig_asimov");

	//stat-only
    h_dcp_th13_cont_1sig_stat = (TH2D*)h_dcp_th13_stat->Clone("h_dcp_th13_cont_1sig_stat");
    h_dcp_th13_cont_2sig_stat = (TH2D*)h_dcp_th13_stat->Clone("h_dcp_th13_cont_2sig_stat");
    h_dcp_th13_cont_3sig_stat = (TH2D*)h_dcp_th13_stat->Clone("h_dcp_th13_cont_3sig_stat");

  }
  else if(app == 1){
    h_dcp_th13_cont_1sig_asimov = (TH2D*)h_dm32_th23_asimov->Clone("h_dm32_th23_cont_1sig_asimov");
    h_dcp_th13_cont_2sig_asimov = (TH2D*)h_dm32_th23_asimov->Clone("h_dm32_th23_cont_2sig_asimov");
    h_dcp_th13_cont_3sig_asimov = (TH2D*)h_dm32_th23_asimov->Clone("h_dm32_th23_cont_3sig_asimov");

	//stat-only 
	h_dcp_th13_cont_1sig_stat = (TH2D*)h_dm32_th23_stat->Clone("h_dm32_th23_cont_1sig_stat");
    h_dcp_th13_cont_2sig_stat = (TH2D*)h_dm32_th23_stat->Clone("h_dm32_th23_cont_2sig_stat");
    h_dcp_th13_cont_3sig_stat = (TH2D*)h_dm32_th23_stat->Clone("h_dm32_th23_cont_3sig_stat");	
  }
  else if(app == 2){
    h_dcp_th13_cont_1sig_asimov = (TH2D*)h_dcp_th23_asimov->Clone("h_dcp_th23_cont_1sig_asimov");
    h_dcp_th13_cont_2sig_asimov = (TH2D*)h_dcp_th23_asimov->Clone("h_dcp_th23_cont_2sig_asimov");
    h_dcp_th13_cont_3sig_asimov = (TH2D*)h_dcp_th23_asimov->Clone("h_dcp_th23_cont_3sig_asimov");

	//stat-only 
	h_dcp_th13_cont_1sig_stat = (TH2D*)h_dcp_th23_stat->Clone("h_dcp_th23_cont_1sig_stat");
    h_dcp_th13_cont_2sig_stat = (TH2D*)h_dcp_th23_stat->Clone("h_dcp_th23_cont_2sig_stat");
    h_dcp_th13_cont_3sig_stat = (TH2D*)h_dcp_th23_stat->Clone("h_dcp_th23_cont_3sig_stat");	
  }

  h_dcp_th13_cont_1sig_asimov->SetDirectory(0);
  h_dcp_th13_cont_2sig_asimov->SetDirectory(0);
  h_dcp_th13_cont_3sig_asimov->SetDirectory(0);

  h_dcp_th13_cont_1sig_stat->SetDirectory(0);
  h_dcp_th13_cont_2sig_stat->SetDirectory(0);
  h_dcp_th13_cont_3sig_stat->SetDirectory(0);

  double min_asimov;
  double p68_asimov, p90_asimov, p95_asimov, p99_asimov, p3sig_asimov;
  double *p_asimov;
  std::cout << "Now calling getInterval for fake data" << std::endl;
  if(app == 0){p_asimov = getInterval2D(h_dcp_th13_asimov, p68_asimov, p90_asimov, p95_asimov, p99_asimov, p3sig_asimov);}
  else if(app == 1){p_asimov = getInterval2D(h_dm32_th23_asimov, p68_asimov, p90_asimov, p95_asimov, p99_asimov, p3sig_asimov);}
  else if(app == 2){p_asimov = getInterval2D(h_dcp_th23_asimov, p68_asimov, p90_asimov, p95_asimov, p99_asimov, p3sig_asimov);}

  double tpp68_asimov[1];
  double tpp90_asimov[1];
  double tpp95_asimov[1];
  double tpp99_asimov[1];
  double tpp3sig_asimov[1];

  tpp68_asimov[0] = p68_asimov;
  tpp90_asimov[0] = p90_asimov;
  tpp95_asimov[0] = p95_asimov;
  tpp99_asimov[0] = p99_asimov;
  tpp3sig_asimov[0] = p3sig_asimov;


  //////////////////////////////////////
  //Same but for stat-only contours
  //////////////////////////////////////
  double min_stat;
  double p68_stat, p90_stat, p95_stat, p99_stat, p3sig_stat;
  double *p_stat;
  std::cout << "Now calling getInterval for stat-only fit" << std::endl;
  if(app == 0){p_stat = getInterval2D(h_dcp_th13_stat, p68_stat, p90_stat, p95_stat, p99_stat, p3sig_stat);}
  else if(app == 1){p_stat = getInterval2D(h_dm32_th23_stat, p68_stat, p90_stat, p95_stat, p99_stat, p3sig_stat);}
  else if(app == 2){p_stat = getInterval2D(h_dcp_th23_stat, p68_stat, p90_stat, p95_stat, p99_stat, p3sig_stat);}

  double tpp68_stat[1];
  double tpp90_stat[1];
  double tpp95_stat[1];
  double tpp99_stat[1];
  double tpp3sig_stat[1];

  tpp68_stat[0] = p68_stat;
  tpp90_stat[0] = p90_stat;
  tpp95_stat[0] = p95_stat;
  tpp99_stat[0] = p99_stat;
  tpp3sig_stat[0] = p3sig_stat;

  h_dcp_th13_stat->SetContour(255);
  h_dcp_th13_cont_1sig_stat->Smooth(1);
  h_dcp_th13_cont_2sig_stat->Smooth(1);
  h_dcp_th13_cont_3sig_stat->Smooth(1);
  h_dcp_th13_stat->Smooth(1);

  //Now set the contours using the appropriate significance level
  if(!levSigs) {
    h_dcp_th13_cont_1sig_asimov->SetContour(1,tpp68_asimov);
    h_dcp_th13_cont_2sig_asimov->SetContour(1,tpp90_asimov);
    h_dcp_th13_cont_3sig_asimov->SetContour(1,tpp99_asimov);

	//stat-only
    h_dcp_th13_cont_1sig_stat->SetContour(1,tpp68_stat);
    h_dcp_th13_cont_2sig_stat->SetContour(1,tpp90_stat);
    h_dcp_th13_cont_3sig_stat->SetContour(1,tpp99_stat);

  }
  else {
    h_dcp_th13_cont_1sig_asimov->SetContour(1,tpp68_asimov);
    h_dcp_th13_cont_2sig_asimov->SetContour(1,tpp95_asimov);
    h_dcp_th13_cont_3sig_asimov->SetContour(1,tpp3sig_asimov);

	//stat-only
    h_dcp_th13_cont_1sig_stat->SetContour(1,tpp68_stat);
    h_dcp_th13_cont_2sig_stat->SetContour(1,tpp95_stat);
    h_dcp_th13_cont_3sig_stat->SetContour(1,tpp3sig_stat);	
  }

  h_dcp_th13_cont_1sig_asimov->SetLineColor(kGreen+2);
  h_dcp_th13_cont_2sig_asimov->SetLineColor(kGreen+2);
  h_dcp_th13_cont_3sig_asimov->SetLineColor(kGreen+2);
  h_dcp_th13_cont_1sig_asimov->SetLineWidth(2);
  h_dcp_th13_cont_2sig_asimov->SetLineWidth(2);
  h_dcp_th13_cont_3sig_asimov->SetLineWidth(2);
  h_dcp_th13_cont_1sig_asimov->SetLineStyle(2);
  h_dcp_th13_cont_3sig_asimov->SetLineStyle(3);

  h_dcp_th13_cont_1sig_stat->SetLineColor(kRed);
  h_dcp_th13_cont_2sig_stat->SetLineColor(kRed);
  h_dcp_th13_cont_3sig_stat->SetLineColor(kRed);
  h_dcp_th13_cont_1sig_stat->SetLineWidth(2);
  h_dcp_th13_cont_2sig_stat->SetLineWidth(2);
  h_dcp_th13_cont_3sig_stat->SetLineWidth(2);
  h_dcp_th13_cont_1sig_stat->SetLineStyle(2);
  h_dcp_th13_cont_3sig_stat->SetLineStyle(3);

  // find best fit point
  int Mbx, Mby, Mbz, MaCh3bin;
  double Mx, My; 
  TGraph *bestfitM = new TGraph(1);

  // find best fit point for stat-only
  int Mbx_stat, Mby_stat, Mbz_stat, MaCh3bin_stat;
  double Mx_stat, My_stat; 
  TGraph *bestfitM_stat = new TGraph(1);


  if(app == 0) {
    MaCh3bin = h_dcp_th13_asimov->GetMaximumBin();
    h_dcp_th13_asimov->GetBinXYZ(MaCh3bin,Mbx,Mby,Mbz);
    Mx = h_dcp_th13_asimov->GetXaxis()->GetBinCenter(Mbx);
    My = h_dcp_th13_asimov->GetYaxis()->GetBinCenter(Mby);

	//stat-only
    MaCh3bin_stat = h_dcp_th13_stat->GetMaximumBin();
    h_dcp_th13_stat->GetBinXYZ(MaCh3bin_stat,Mbx_stat, Mby_stat, Mbz_stat);
    Mx_stat = h_dcp_th13_stat->GetXaxis()->GetBinCenter(Mbx_stat);
    My_stat = h_dcp_th13_stat->GetYaxis()->GetBinCenter(Mby_stat);
  }
  else if(app == 1){
    MaCh3bin = h_dm32_th23_asimov->GetMaximumBin();
    h_dm32_th23_asimov->GetBinXYZ(MaCh3bin,Mbx,Mby,Mbz);
    Mx = h_dm32_th23_asimov->GetXaxis()->GetBinCenter(Mbx);
    My = h_dm32_th23_asimov->GetYaxis()->GetBinCenter(Mby);

	//stat-only
    MaCh3bin_stat = h_dm32_th23_stat->GetMaximumBin();
    h_dm32_th23_stat->GetBinXYZ(MaCh3bin_stat, Mbx_stat, Mby_stat, Mbz_stat);
    Mx_stat = h_dm32_th23_stat->GetXaxis()->GetBinCenter(Mbx_stat);
    My_stat = h_dm32_th23_stat->GetYaxis()->GetBinCenter(Mby_stat);
  }
  else if(app == 2){
    MaCh3bin = h_dcp_th23_asimov->GetMaximumBin();
    h_dcp_th23_asimov->GetBinXYZ(MaCh3bin,Mbx,Mby,Mbz);
    Mx = h_dcp_th23_asimov->GetXaxis()->GetBinCenter(Mbx);
    My = h_dcp_th23_asimov->GetYaxis()->GetBinCenter(Mby);

	//stat-only
    MaCh3bin_stat = h_dcp_th23_stat->GetMaximumBin();
    h_dcp_th23_stat->GetBinXYZ(MaCh3bin_stat, Mbx_stat, Mby_stat, Mbz_stat);
    Mx_stat = h_dcp_th23_stat->GetXaxis()->GetBinCenter(Mbx_stat);
    My_stat = h_dcp_th23_stat->GetYaxis()->GetBinCenter(Mby_stat);
  }

  bestfitM->SetPoint(0,Mx,My);
  bestfitM->SetMarkerStyle(22);
  bestfitM->SetMarkerSize(1);
  bestfitM->SetMarkerColor(kGreen+2);

  //stat-only best-fit
  bestfitM_stat->SetPoint(0,Mx_stat, My_stat);
  bestfitM_stat->SetMarkerStyle(34);
  bestfitM_stat->SetMarkerSize(1);
  bestfitM_stat->SetMarkerColor(kRed);

  TCanvas* c = new TCanvas("c","c",600,500);
  //c->Print("jarl.pdf]");

  c->Draw();
  c->cd();

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

  if(app == 1 && hierarchy==0) leg = new TLegend(0.6,0.53,0.85,0.88);
  else if(app == 1 && hierarchy!=0) leg = new TLegend(0.65,0.53,0.90,0.88);
  else leg = new TLegend(0.55,0.53,0.85,0.88);
  //TLegend* leg = new TLegend(0.7,0.35,1.2,0.65);
  leg->SetFillStyle(0);
  //leg->SetLineColor(0);
  leg->SetLineColorAlpha(0,0);
  leg->SetTextFont(132);
  leg->SetTextSize(0.025);
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
  leg->AddEntry(h_dcp_th13_cont_2sig_asimov, nominal_title, "l");
  leg->AddEntry(bestfitM,"Best Fit","p");
  if(!levSigs) {
    leg->AddEntry(h_dcp_th13_cont_1sig_asimov,"68%","l");
    leg->AddEntry(h_dcp_th13_cont_2sig_asimov,"90%","l");
    leg->AddEntry(h_dcp_th13_cont_3sig_asimov,"99%","l");

    leg->AddEntry(h_dcp_th13_cont_2sig_stat, stat_title, "l");
    leg->AddEntry(bestfitM_stat,"Best Fit","p");

	//stat-only
    leg->AddEntry(h_dcp_th13_cont_1sig_stat,"68%","l");
    leg->AddEntry(h_dcp_th13_cont_2sig_stat,"90%","l");
    leg->AddEntry(h_dcp_th13_cont_3sig_stat,"99%","l");
  }
  else {
    leg->AddEntry(h_dcp_th13_cont_1sig_asimov,"1#sigma","l");
    leg->AddEntry(h_dcp_th13_cont_2sig_asimov,"2#sigma","l");
    leg->AddEntry(h_dcp_th13_cont_3sig_asimov,"3#sigma","l");


    leg->AddEntry(h_dcp_th13_cont_2sig_stat, stat_title, "l");
    leg->AddEntry(bestfitM_stat,"Best Fit","p");
	//stat-only
    leg->AddEntry(h_dcp_th13_cont_1sig_stat,"1#sigma","l");
    leg->AddEntry(h_dcp_th13_cont_2sig_stat,"2#sigma","l");
    leg->AddEntry(h_dcp_th13_cont_3sig_stat,"3#sigma","l");
  }
  if(drawAsimovPoint) leg->AddEntry(g_asim_dcp_th13,"Asimov Point (NH)","p");

/*
  // dummy th2 to extend the y-axis
  TH2D* dummy = new TH2D("dummy","dummy",10,-0.05,0.05,10,0.3,0.7);
  if(hierarchy == 1)dummy->SetTitle("Jarlskog Invariant, Normal Hierarchy;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;sin^{2}#theta_{23}");
  else if(hierarchy == 0)dummy->SetTitle("Jarlskog Invariant, Both Hierarchies;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;sin^{2}#theta_{23}");
  else if(hierarchy == -1)dummy->SetTitle("Jarlskog Invariant, Inverted Hierarchy;J #equiv s_{13}c^{2}_{13}s_{12}c_{12}s_{23}c_{23}sin#delta;sin^{2}#theta_{23}");
  dummy->GetXaxis()->SetNdivisions(508);
*/

  //j_th23_cont_3sig->GetYaxis()->SetRangeUser(0.3,0.8);
  //dummy->Draw();
  //j_th23->Draw("colz");
  //h_dcp_th13->Draw("colz");
  h_dcp_th13_cont_1sig_asimov->GetXaxis()->SetNdivisions(507);

  if(app == 1 && hierarchy == 0) {

    TPad *p_NH = new TPad("p_NH","p_NH",0.0,0.5,1.0,1.0);
    TPad *p_IH = new TPad("p_IH","p_IH",0.0,0.0,1.0,0.5);
    p_NH->SetBottomMargin(0.01);
    p_NH->SetTopMargin(0.22);
    p_IH->SetBottomMargin(0.22);
    p_IH->SetTopMargin(0.01);
    p_NH->SetLeftMargin(0.15);
    p_IH->SetLeftMargin(0.15);
    p_NH->Draw(); p_IH->Draw();
    TH2D* dummyNH = (TH2D*)h_dm32_th23_asimov->Clone("dummyNH");
    TH2D* dummyIH = (TH2D*)h_dm32_th23_asimov->Clone("dummyNH");
    dummyNH->Reset(); dummyIH->Reset();
    dummyNH->GetYaxis()->SetRangeUser(0.00225,0.00275);
    dummyIH->GetYaxis()->SetRangeUser(-0.00275,-0.00225);
    dummyNH->GetXaxis()->SetLabelSize(0.00);
    dummyNH->GetYaxis()->SetLabelSize(0.09);
    dummyNH->GetXaxis()->SetTitle("");
    dummyNH->GetYaxis()->SetTitleSize(0.1);
    dummyNH->GetYaxis()->SetTitleOffset(0.55);
    dummyIH->GetXaxis()->SetTitleSize(0.1);
    dummyIH->GetXaxis()->SetTitleOffset(0.9);
    dummyIH->GetXaxis()->SetLabelSize(0.09);
    dummyIH->GetYaxis()->SetLabelSize(0.09);
    dummyIH->GetYaxis()->SetTitle("");
    dummyIH->SetTitle("");

    p_NH->cd();
    dummyNH->Draw();
    h_dcp_th13_cont_1sig_asimov->Draw("cont3 same");
    h_dcp_th13_cont_2sig_asimov->Draw("cont3 same");
    h_dcp_th13_cont_3sig_asimov->Draw("cont3 same");

	//now draw the stat as well!
    h_dcp_th13_cont_1sig_stat->Draw("cont3 same");
    h_dcp_th13_cont_2sig_stat->Draw("cont3 same");
    h_dcp_th13_cont_3sig_stat->Draw("cont3 same");

    bestfitM->Draw("SAME.P");
    bestfitM_stat->Draw("SAME.P");

    if(drawAsimovPoint) g_asim_dm2_th23->Draw("same p");
     
    p_IH->cd();
    dummyIH->Draw();
    h_dcp_th13_cont_1sig_asimov->Draw("cont3 same");
    h_dcp_th13_cont_2sig_asimov->Draw("cont3 same");
    h_dcp_th13_cont_3sig_asimov->Draw("cont3 same");

	//now draw the stat as well
    h_dcp_th13_cont_1sig_stat->Draw("cont3 same");
    h_dcp_th13_cont_2sig_stat->Draw("cont3 same");
    h_dcp_th13_cont_3sig_stat->Draw("cont3 same");

    bestfitM->Draw("SAME.P");
    //bestfitM_fake->Draw("SAME.P");
    bestfitM_stat->Draw("SAME.P");
    leg->Draw("same");
  }
  else {
    if(app == 1 && hierarchy==1){h_dcp_th13_cont_1sig_asimov->GetYaxis()->SetRangeUser(0.0024,0.0027);}
    //if(app == 1 && hierarchy==1){h_dcp_th13_cont_1sig_asimov->GetYaxis()->SetRangeUser(0.0023,0.0028);}
    else if(app == 1 && hierarchy==-1){h_dcp_th13_cont_1sig_asimov->GetYaxis()->SetRangeUser(-0.0028,-0.0023);}

    h_dcp_th13_cont_1sig_asimov->Draw("cont3 same");
    h_dcp_th13_cont_2sig_asimov->Draw("cont3 same");
    h_dcp_th13_cont_3sig_asimov->Draw("cont3 same");

    h_dcp_th13_cont_1sig_stat->Draw("cont3 same");
    h_dcp_th13_cont_2sig_stat->Draw("cont3 same");
    h_dcp_th13_cont_3sig_stat->Draw("cont3 same");

    bestfitM->Draw("SAME.P");
    bestfitM_stat->Draw("SAME.P");

    if(drawAsimovPoint && app == 0){g_asim_dcp_th13->Draw("same p");}
    else if(drawAsimovPoint && app == 1){g_asim_dm2_th23->Draw("same p");}
    else if(drawAsimovPoint && app == 2){g_asim_dcp_th23->Draw("same p");}
    //h_dcp_th13_cont_2sig_joint->Draw("cont3 same");
    //h_dcp_th13_cont_2sig_asimov->Draw("cont3 same");
    //h_dcp_th13_cont_2sig_nova->Draw("cont3 same");
    leg->Draw("same");
    if(hierarchy!=0) hlab->Draw();
  }

  if(!levSigs) {
    if(RC) {
      if(app == 0) {
        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth13_NH_wRC_perc.pdf");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth13_both_wRC_perc.pdf");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth13_IH_wRC_perc.pdf");

        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth13_NH_wRC_perc.png");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth13_both_wRC_perc.png");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth13_IH_wRC_perc.png");

      }
      else if(app == 1){
        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dm2th23_NH_wRC_perc.pdf");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dm2th23_both_wRC_perc.pdf");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dm2th23_IH_wRC_perc.pdf");

        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dm2th23_NH_wRC_perc.png");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dm2th23_both_wRC_perc.png");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dm2th23_IH_wRC_perc.png");
      }
	  else if(app == 2){
        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth23_NH_wRC_perc.pdf");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth23_both_wRC_perc.pdf");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth23_IH_wRC_perc.pdf");

        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth23_NH_wRC_perc.png");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth23_both_wRC_perc.png");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth23_IH_wRC_perc.png");
	  }
    }
    else {
      if(app == 0) {
        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth13_NH_woRC_perc.pdf");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth13_both_woRC_perc.pdf");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth13_IH_woRC_perc.pdf");

        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth13_NH_woRC_perc.png");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth13_both_woRC_perc.png");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth13_IH_woRC_perc.png");
      }
      else if(app == 1){
        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dm2th23_NH_woRC_perc.pdf");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dm2th23_both_woRC_perc.pdf");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dm2th23_IH_woRC_perc.pdf");

        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dm2th23_NH_woRC_perc.png");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dm2th23_both_woRC_perc.png");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dm2th23_IH_woRC_perc.png");
      }
      else if(app == 2){
        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth23_NH_woRC_perc.pdf");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth23_both_woRC_perc.pdf");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth23_IH_woRC_perc.pdf");

        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth23_NH_woRC_perc.png");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth23_both_woRC_perc.png");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth23_IH_woRC_perc.png");
      }
    }
  }
  else {
    if(RC) {
      if(app == 0) {
        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth13_NH_wRC_sigs.pdf");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth13_both_wRC_sigs.pdf");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth13_IH_wRC_sigs.pdf");

        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth13_NH_wRC_sigs.png");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth13_both_wRC_sigs.png");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth13_IH_wRC_sigs.png");
      }
      else if(app == 1){
        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dm2th23_NH_wRC_sigs.pdf");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dm2th23_both_wRC_sigs.pdf");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dm2th23_IH_wRC_sigs.pdf");

        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dm2th23_NH_wRC_sigs.png");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dm2th23_both_wRC_sigs.png");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dm2th23_IH_wRC_sigs.png");
      }
      else if(app == 2){
        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth23_NH_wRC_sigs.pdf");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth23_both_wRC_sigs.pdf");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth23_IH_wRC_sigs.pdf");

        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth23_NH_wRC_sigs.png");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth23_both_wRC_sigs.png");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth23_IH_wRC_sigs.png");
      }
    }
    else {
      if(app == 0) {
        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth13_NH_woRC_sigs.pdf");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth13_both_woRC_sigs.pdf");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth13_IH_woRC_sigs.pdf");

        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth13_NH_woRC_sigs.png");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth13_both_woRC_sigs.png");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth13_IH_woRC_sigs.png");
      }
      else if(app == 1){
        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dm2th23_NH_woRC_sigs.pdf");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dm2th23_both_woRC_sigs.pdf");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dm2th23_IH_woRC_sigs.pdf");

        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dm2th23_NH_woRC_sigs.png");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dm2th23_both_woRC_sigs.png");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dm2th23_IH_woRC_sigs.png");
      }
      else if(app == 2){
        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth23_NH_woRC_sigs.pdf");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth23_both_woRC_sigs.pdf");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth23_IH_woRC_sigs.pdf");

        if(hierarchy == 1) c->Print(outdir+"/"+outdir+"_dcpth23_NH_woRC_sigs.png");
        else if(hierarchy == 0) c->Print(outdir+"/"+outdir+"_dcpth23_both_woRC_sigs.png");
        else if(hierarchy == -1) c->Print(outdir+"/"+outdir+"_dcpth23_IH_woRC_sigs.png");
      }
    }
  }

  delete t_asimov;
  delete f_asimov;

  delete t_stat;
  delete f_stat;

  delete c;

  ///////////////
  // Now write the contours out to file
  ///////////////
  
  if(saveConts){
	//std::cout << "Now saving the contours in file " << outfile_asimov_name << std::endl;
	//saveContours(outfile_asimov_name, RC, app, levOpt, hierarchy, h_dcp_th13_cont_1sig_asimov, h_dcp_th13_cont_2sig_asimov, h_dcp_th13_cont_3sig_asimov); 
	//std::cout << "Now saving the contours in file " << outfile_stat_name << std::endl;
	//saveContours(outfile_stat_name, RC, app, levOpt, hierarchy, h_dcp_th13_cont_1sig_stat, h_dcp_th13_cont_2sig_stat, h_dcp_th13_cont_3sig_stat); 
  }

}

//(TString infile_1_name, TString infile_2_name, TString infile_3_name, TString title, TString title_1, TString title_2, TString title_3, TString out_directory, int hierarchy, bool LogY=false, int octant=0, int levOpt=0, bool doRC = false, bool apply_RC_to_asimov = false, bool apply_RC_to_stat = false, bool apply_RC_to_fds = false, int burn_in=150000, int burn_in_stat=20000, int burn_in_fake = 100000) {




/*
void draw_fakeDataComp_2D(TString infile_1_name, TString infile_2_name, TString outfile_1_name, TString outfile_2_name, TString title, TString title_1, TString title_2, TString out_directory, int hierarchy, int app=0 ,bool RC=false, bool stat_apply_RC = true,  bool asimov_apply_RC = true, int octant=0, int levOpt=0, int burn_in=150000, int burn_in_stat=50000) {

  fake_data_file = infile_1_name;
  outfile_asimov_name = outfile_1_name;
  title_asimov = title_1;

  stat_only_file = infile_2_name;
  outfile_stat_name = outfile_2_name;
  title_stat = title_2;


  plot_title = title;

  outdir = out_directory;
  makePlot(hierarchy, app, RC, levOpt, stat_apply_RC, asimov_apply_RC, burn_in, burn_in_stat);

}
*/

void draw_fakeDataComp_2D(TString infile_1_name, TString infile_2_name, TString title, TString title_1, TString title_2, TString out_directory, int hierarchy, bool LogY=false, int octant=0, int levOpt=0, bool doRC = false, bool apply_RC_to_asimov = false, bool apply_RC_to_stat = false, int burn_in=80000, int burn_in_stat=80000) {


  nominal_file = infile_1_name;
  nominal_title = title_1;

  stat_file = infile_2_name;
  stat_title = title_2;

  plot_title = title;

  outdir = out_directory;
  std::cout << "Using " << nominal_file << " for the nominal file" << std::endl;
  std::cout << "Using " << stat_file << " for the stat only fit" << std::endl;

  std::cout << "Using titles " << nominal_title << ", " << stat_title << std::endl;
/*
  makePlot(hierarchy, 0, LogY, octant, levOpt, burn_in, burn_in_stat, nominal_title, stat_title, doRC, apply_RC_to_asimov, apply_RC_to_stat);
  makePlot(hierarchy, 1, LogY, octant, levOpt, burn_in, burn_in_stat, nominal_title, stat_title, doRC, apply_RC_to_asimov, apply_RC_to_stat);
  makePlot(hierarchy, 2, LogY, octant, levOpt, burn_in, burn_in_stat, nominal_title, stat_title, doRC, apply_RC_to_asimov, apply_RC_to_stat);
*/

  makePlot(hierarchy, 0);
  makePlot(1, 1);
  makePlot(hierarchy, 2);

  return;
}

/*
void krw_drawContours() {


  ///////
  //wRC
  ///////
  //app
  std::cout << "~~~~~~~" << std::endl;
  std::cout << "Making ploys for dcp vs. sin2th13 for both hierarchies" << std::endl;
  makePlot(0, 0, true);//both Hierarchy
  std::cout << "Making ploys for dcp vs. sin2th13 for NH" << std::endl;
  makePlot(1, 0, true);// NH
  std::cout << "Making ploys for dcp vs. sin2th13 for IH" << std::endl;
  makePlot(-1, 0, true);// IH

  //disapp
  makePlot(0, 1, true);//both Hierarchy
  makePlot(1, 1, true);// NH
  makePlot(-1, 1, true);// IH

  //dcp vs. sin2th23
  makePlot(0, 2, true);//both Hierarchy
  makePlot(1, 2, true);// NH
  makePlot(-1, 2, true);// IH
 
  //levOpt == 1 -> print out 68, 90% and 99% credible intervals
  //app
  makePlot(0, 0, true, 1);//both Hierarchy
  makePlot(1, 0, true, 1);// NH
  makePlot(-1, 0, true, 1);// IH

  //disapp
  makePlot(0, 1, true, 1);//both Hierarchy
  makePlot(1, 1, true, 1);// NH
  makePlot(-1, 1, true, 1);// IH

  //dcp vs. sin2th23
  makePlot(0, 2, true, 1);//both Hierarchy
  makePlot(1, 2, true, 1);// NH
  makePlot(-1, 2, true, 1);// IH
  

  //throw;

}
*/
