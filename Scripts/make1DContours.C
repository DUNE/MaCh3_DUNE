//////////////////////// Macro by Kevin Wood /////////////////////////////////////
//  Need to make directory 'Contours1D' for script to write output to //////////////
//  Makes all contours (woRC and wRC) and percent/sigma intervals   //////////////
//////////////////////////////////////////////////////////////////////////////////
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

// WARNING: ASSUMPTIONS EXIST IN THIS FUNCTION!!!!!
std::string get1Dcred3sig_delta_cp(TH1D* hist, int actuallyNsig=3) {
  // return range as string:
  std::string range;

  double level = 0.9973;
  if(actuallyNsig==2) level = 0.9545;
  if(actuallyNsig==1) level = 0.68;

  TH1D* hist_copy = (TH1D*)hist->Clone("hist_copy");
  double integral, tsum=0.;
  integral = hist_copy->Integral();
  double xval;
  double xwidth;

  // get *worst* fit
  double wf = hist_copy->GetXaxis()->GetBinCenter(hist_copy->GetMinimumBin());
  double up_pos = 9999999.;
  double low_pos = -9999999.;

  while((tsum/integral)<level) {
    double tmax = hist_copy->GetMaximum();
    int bin = hist_copy->GetMaximumBin();
    xval = hist_copy->GetXaxis()->GetBinCenter(bin);
    xwidth = hist_copy->GetXaxis()->GetBinWidth(bin);
    if((tsum/integral)<level) {
      hist_copy->SetBinContent(bin,-1.0);
      if(xval>0. && xval<wf && xval>low_pos) low_pos = xval + xwidth/2.;
      else if(xval>0. && xval>wf && xval<up_pos) up_pos = xval - xwidth/2.;
    tsum+=tmax;
  }

    range = "[-$\\pi$,"+std::to_string(low_pos)+"] $\\cup$ ["+std::to_string(up_pos)+",$\\pi$]";
  }
  return range;

}

//std::string get1Dcred_dm2(TH1d* hist, double level, double &bf, double &up_IH, double &low_IH, double &up_NH, double &low_IH) {
std::string get1Dcred_dm2(TH1D* hist, double level) {

  // return range as string:
  std::string range;

  if(level>1.) {
    std::cout << "credible interval must be < 1." << std::endl;
    throw;
  }

  TH1D* hist_copy = (TH1D*)hist->Clone("hist_copy");
  double integral, tsum=0.;
  integral = hist_copy->Integral();
  double xval;
  double xwidth;

  // get best fit
  double bf = hist_copy->GetXaxis()->GetBinCenter(hist_copy->GetMaximumBin());
  double up_IH = -9999999.;
  double low_IH = 9999999.;
  double up_NH = -9999999.;
  double low_NH = 9999999.;

  while((tsum/integral)<level) {
    double tmax = hist_copy->GetMaximum();
    int bin = hist_copy->GetMaximumBin();
    xval = hist_copy->GetXaxis()->GetBinCenter(bin);
    xwidth = hist_copy->GetXaxis()->GetBinWidth(bin);
    if((tsum/integral)<level) {
      hist_copy->SetBinContent(bin,-1.0);
      if(xval>0.) {
        if(xval<low_NH && xval<bf) low_NH = xval - xwidth/2.;
        if(xval>up_NH && xval>bf) up_NH = xval + xwidth/2.;
      }
      else if(xval<0.) {
        if(xval<low_IH) low_IH = xval - xwidth/2.;
        if(xval>up_IH) up_IH = xval + xwidth/2.;
      }
    }
    tsum+=tmax;
  }

  //unit conversion and rounding:
  low_NH *= 1.e3;
  low_IH *= 1.e3;
  up_NH *= 1.e3;
  up_IH *= 1.e3;

  // WARNING: hardcoding the NH preference here, because I'm lazy
  if(low_IH>10. || up_IH<-10.) range = "["+std::to_string(low_NH)+","+std::to_string(up_NH)+"]";
  else range = "["+std::to_string(low_IH)+","+std::to_string(up_IH)+"] $\\cup$ ["+std::to_string(low_NH)+","+std::to_string(up_NH)+"]";

  return range;

}

TH1D* get1Dcred(TH1D* hist, double level, double &bf, double &up, double &low, int &nbins) {

  if(level>1.) {
    std::cout << "credible interval must be < 1." << std::endl;
    throw;
  }

  TH1D* hist_copy = (TH1D*)hist->Clone("hist_copy");
  TH1D* hist_ret = (TH1D*)hist->Clone("hist_ret");
  double integral, tsum=0.;
  integral = hist_copy->Integral();
  double xval;
  double xwidth;

  // get best fit
  bf = hist_copy->GetXaxis()->GetBinCenter(hist_copy->GetMaximumBin());
  up = -9999999.;
  low = 9999999.;
  nbins = 0;

  while((tsum/integral)<level) {
    double tmax = hist_copy->GetMaximum();
    int bin = hist_copy->GetMaximumBin();
    xval = hist_copy->GetXaxis()->GetBinCenter(bin);
    xwidth = hist_copy->GetXaxis()->GetBinWidth(bin);
    if((tsum/integral)<level) {
      hist_copy->SetBinContent(bin,-1.0);
      hist_ret->SetBinContent(bin,0.);
      if(xval<low && xval<bf) low = xval - xwidth/2.;
      if(xval>up && xval>bf) up = xval + xwidth/2.;
      nbins++;
    }
    tsum+=tmax;
  }

  return hist_ret;

}

// hieararchy = 1-->NH, 0-->both, -1-->IH
void makePlot(TString ReducedChain, int hierarchy, bool doRC = false, bool LogY = false, bool levSigs = true)
{
  TFile* f = new TFile(ReducedChain);
  TTree* t = (TTree*)f->Get("osc_posteriors");

  TH1D* delta_cp_hist = new TH1D("delta_cp_hist",";#delta_{CP};posterior probability",180,-1.*TMath::Pi(),TMath::Pi());
  TH1D* th23_hist = new TH1D("th23_hist",";sin^{2}#theta_{23};posterior probability",150,0.35,0.65);
  TH1D* th13_hist;
  if(doRC) th13_hist = new TH1D("th13_hist",";sin^{2}#theta_{13};posterior probability",150,0.019,0.025);
  else th13_hist = new TH1D("th13_hist",";sin^{2}#theta_{13};posterior probability",150,0.012,0.05);
  TH1D* dm2_hist = new TH1D("dm2_hist",";#Delta m^{2}_{32};posterior probability",1000,-0.003,0.003);

  //t->Draw("delta_cp>>j_hist_fsd","abs(cos(delta_cp))*RCreweight*(step>80000)");
  TH1D* hist_arr[4] = {delta_cp_hist,th23_hist,th13_hist,dm2_hist};
  for(unsigned int i=0; i<4;i++) {
    hist_arr[i]->GetXaxis()->SetTitleFont(132) ;
    hist_arr[i]->GetXaxis()->SetTitleOffset(0.9) ;
    hist_arr[i]->GetXaxis()->SetTitleSize(0.07) ;
    hist_arr[i]->GetXaxis()->SetLabelFont(132) ;
    hist_arr[i]->GetXaxis()->SetLabelSize(0.05) ;
    hist_arr[i]->GetYaxis()->SetTitleFont(132) ;
    hist_arr[i]->GetYaxis()->SetTitleOffset(0.6) ;
    hist_arr[i]->GetYaxis()->SetTitleSize(0.07) ;
    hist_arr[i]->GetYaxis()->SetLabelFont(132) ;
    hist_arr[i]->GetYaxis()->SetLabelSize(0.05) ;
  }

  if(hierarchy == 1) {
    t->Draw("dcp>>delta_cp_hist","(step>80000)*(dm23>0.)");
    t->Draw("theta23>>th23_hist","(step>80000)*(dm23>0.)");
    t->Draw("theta13>>th13_hist","(step>80000)*(dm23>0.)");
    t->Draw("dm23>>dm2_hist","(step>80000)*(dm23>0.)");
  }
  else if(hierarchy == 0) {
    t->Draw("dcp>>delta_cp_hist","(step>80000)");
    t->Draw("theta23>>th23_hist","(step>80000)");
    t->Draw("theta13>>th13_hist","(step>80000)");
    t->Draw("dm23>>dm2_hist","(step>80000)");
  }
  else if(hierarchy == -1) {
    t->Draw("dcp>>delta_cp_hist","(step>80000)*(dm23<0.)");
    t->Draw("theta23>>th23_hist","(step>80000)*(dm23<0.)");
    t->Draw("theta13>>th13_hist","(step>80000)*(dm23<0.)");
    t->Draw("dm23>>dm2_hist","(step>80000)*(dm23<0.)");
  }
  else {
   std::cout << "Error: invalid hierarchy option. 1 for NH, 0 for both, -1 for IH" <<std::endl;
   throw;
  }

  delta_cp_hist->GetYaxis()->SetLabelSize(0.);
  th23_hist->GetYaxis()->SetLabelSize(0.);
  th13_hist->GetYaxis()->SetLabelSize(0.);
  dm2_hist->GetYaxis()->SetLabelSize(0.);

  //int delta_cp_maxbin, th23_maxbin, th13_maxbin, dm2_maxbin;
  double delta_cp_bf, th23_bf, th13_bf, dm2_bf;
  double delta_cp_up, th23_up, th13_up, dm2_up;
  double delta_cp_lo, th23_lo, th13_lo, dm2_lo;
  int delta_cp_nbins, th23_nbins, th13_nbins, dm2_nbins;

  double delta_cp_up2sig, th23_up2sig, th13_up2sig, dm2_up2sig;
  double delta_cp_lo2sig, th23_lo2sig, th13_lo2sig, dm2_lo2sig;
  int delta_cp_nbins2sig, th23_nbins2sig, th13_nbins2sig, dm2_nbins2sig;

  double delta_cp_up3sig, th23_up3sig, th13_up3sig, dm2_up3sig;
  double delta_cp_lo3sig, th23_lo3sig, th13_lo3sig, dm2_lo3sig;
  int delta_cp_nbins3sig, th23_nbins3sig, th13_nbins3sig, dm2_nbins3sig;

  double dummy_bf;

  double level_1 = 0.68;
  double level_2 = 0.9545;
  double level_3 = 0.9973;
  if(!levSigs) {
    level_1 = 0.68;
    level_2 = 0.90;
    level_3 = 0.99;
  }

  TH1D* delta_cp_cred = get1Dcred(delta_cp_hist,level_1,delta_cp_bf,delta_cp_up,delta_cp_lo,delta_cp_nbins);
  TH1D* th23_cred = get1Dcred(th23_hist,level_1,th23_bf,th23_up,th23_lo,th23_nbins);
  TH1D* th13_cred = get1Dcred(th13_hist,level_1,th13_bf,th13_up,th13_lo,th13_nbins);
  TH1D* dm2_cred = get1Dcred(dm2_hist,level_1,dm2_bf,dm2_up,dm2_lo,dm2_nbins);

  TH1D* delta_cp_cred2sig = get1Dcred(delta_cp_hist,level_2,dummy_bf,delta_cp_up2sig,delta_cp_lo2sig,delta_cp_nbins2sig);
  TH1D* th23_cred2sig = get1Dcred(th23_hist,level_2,dummy_bf,th23_up2sig,th23_lo2sig,th23_nbins2sig);
  TH1D* th13_cred2sig = get1Dcred(th13_hist,level_2,dummy_bf,th13_up2sig,th13_lo2sig,th13_nbins2sig);
  TH1D* dm2_cred2sig = get1Dcred(dm2_hist,level_2,dummy_bf,dm2_up2sig,dm2_lo2sig,dm2_nbins2sig);

  TH1D* delta_cp_cred3sig = get1Dcred(delta_cp_hist,level_3,dummy_bf,delta_cp_up3sig,delta_cp_lo3sig,delta_cp_nbins3sig);
  TH1D* th23_cred3sig = get1Dcred(th23_hist,level_3,dummy_bf,th23_up3sig,th23_lo3sig,th23_nbins3sig);
  TH1D* th13_cred3sig = get1Dcred(th13_hist,level_3,dummy_bf,th13_up3sig,th13_lo3sig,th13_nbins3sig);
  TH1D* dm2_cred3sig = get1Dcred(dm2_hist,level_3,dummy_bf,dm2_up3sig,dm2_lo3sig,dm2_nbins3sig);

  // get1Dcred doesn't handle dm2 well, so recalculate all of the intervals with this function:
  std::string dm2_cred1sig_str = get1Dcred_dm2(dm2_hist,level_1);
  std::string dm2_cred2sig_str = get1Dcred_dm2(dm2_hist,level_2);
  std::string dm2_cred3sig_str = get1Dcred_dm2(dm2_hist,level_3);;

  std::string delta_cp_cred3sig_str = get1Dcred3sig_delta_cp(delta_cp_hist);
  std::string delta_cp_cred2sig_str = "dummy";
  std::string delta_cp_cred1sig_str = "dummy";
  if(!doRC&&hierarchy==1) delta_cp_cred1sig_str = get1Dcred3sig_delta_cp(delta_cp_hist,1);
  if(!doRC) delta_cp_cred2sig_str = get1Dcred3sig_delta_cp(delta_cp_hist,2);

  std::cout << "Printing 68\% and 90\% credible intervals: " << std::endl;

  if (!levSigs){
    std::cout << "hierarchy = " << hierarchy << "doRC = " << doRC << std::endl << std::endl;
    std::cout << "sin2th23   " << th23_bf << " & [" << th23_lo << "," << th23_up << "]" << " & [" << th23_lo2sig << "," << th23_up2sig << "]"  << " & [" << th23_lo3sig << "," << th23_up3sig << "]     \\\\" << std::endl;
    std::cout << "dm32:      " << dm2_bf*1.e3 << " & " << dm2_cred1sig_str << " & " << dm2_cred2sig_str  << " & " << dm2_cred3sig_str << "\\\\" << std::endl;
    std::cout << "sin2th13:  " << th13_bf << " & [" << th13_lo << "," << th13_up << "]" << " & [" << th13_lo2sig << "," << th13_up2sig << "]"  << " & [" << th13_lo3sig << "," << th13_up3sig << "]     \\\\" << std::endl;
    std::cout << "delta_cp     :  " << delta_cp_bf << " & [" << delta_cp_lo << "," << delta_cp_up << "]" << " & [" << delta_cp_lo2sig << "," << delta_cp_up2sig << "]"  << " & [" << delta_cp_lo3sig << "," << delta_cp_up3sig << "]     \\\\" << std::endl;
  }

  if(levSigs) {
    std::cout << "hierarchy = " << hierarchy << "doRC = " << doRC << std::endl << std::endl;
    std::cout << "delta_cp credible interval (" << delta_cp_nbins << "): " << delta_cp_bf << " + " << delta_cp_up - delta_cp_bf << " - " << delta_cp_bf - delta_cp_lo << std::endl;
    std::cout << "                       [" << delta_cp_lo << "," << delta_cp_up << "]" << std::endl;

    std::cout << "th23 credible interval (" << th23_nbins << "): " << th23_bf << " + " << th23_up - th23_bf << " - " << th23_bf - th23_lo << std::endl;
    std::cout << "                       [" << th23_lo << "," << th23_up << "]" << std::endl;

    std::cout << "th13 credible interval (" << th13_nbins << "): " << th13_bf << " + " << th13_up - th13_bf << " - " << th13_bf - th13_lo << std::endl;
    std::cout << "                       [" << th13_lo << "," << th13_up << "]" << std::endl;

    std::cout << "dm2 credible interval (" << dm2_nbins << "): "  << dm2_bf << " + " << dm2_up - dm2_bf << " - " << dm2_bf - dm2_lo << std::endl;
    std::cout << "                       [" << dm2_lo << "," << dm2_up << "]" << std::endl;

    std::cout << std::endl <<  std::endl << std::endl << std::endl;
    std::cout << "\\begin{table}" << std::endl;
    std::cout << "\\centering" << std::endl;
    std::cout << "\\begin{tabular}{|c|c|c|c|c|}" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "& Central Value &1$\\sigma$ C.I. & 2$\\sigma$ C.I.& 3$\\sigma$ C.I. \\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    if(!doRC&&hierarchy==1) std::cout << "$\\deltacp$ [rad.] &" << delta_cp_bf << " & " << delta_cp_cred1sig_str  << " & " << delta_cp_cred2sig_str << " & " << delta_cp_cred3sig_str << "\\\\" << std::endl;
    else if(!doRC) std::cout << "$\\deltacp$ [rad.] &" << delta_cp_bf << " & [" << delta_cp_lo << "," << delta_cp_up << "]" << " & " << delta_cp_cred2sig_str << " & " << delta_cp_cred3sig_str << "\\\\" << std::endl;
    else std::cout << "$\\deltacp$ [rad.] &" << delta_cp_bf << " & [" << delta_cp_lo << "," << delta_cp_up << "]" << " & [" << delta_cp_lo2sig << "," << delta_cp_up2sig << "]"  << " & " << delta_cp_cred3sig_str << "\\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$\\sin^{2}\\theta_{23}$ &" << th23_bf << " & [" << th23_lo << "," << th23_up << "]" << " & [" << th23_lo2sig << "," << th23_up2sig << "]"  << " & [" << th23_lo3sig << "," << th23_up3sig << "]     \\\\" << std::endl;
//  std::cout << "$\\sin^{2}\\theta_{23}$ &" << th23_bf << " & [" << th23_lo << "," << th23_up << "] \\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "$\\sin^{2}\\theta_{13}$ &" << th13_bf << " & [" << th13_lo << "," << th13_up << "]" << " & [" << th13_lo2sig << "," << th13_up2sig << "]"  << " & [" << th13_lo3sig << "," << th13_up3sig << "]     \\\\" << std::endl;

//  std::cout << "$\\sin^{2}\\theta_{13}$ &" << th13_bf << " & [" << th13_lo << "," << th13_up << "] \\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    //std::cout << "$\\Delta m^{2}_{32}$ [$\\times10^{-3}$ eV$^{2}$] &" << dm2_bf << " & [" << dm2_lo << "," << dm2_up << "]" << " & [" << dm2_lo2sig << "," << dm2_up2sig << "]"  << " & [" << dm2_lo3sig << "," << dm2_up3sig << "]     \\\\" << std::endl;
    std::cout << "$\\Delta m^{2}_{32}$ [$\\times10^{-3}$ eV$^{2}$] &" << dm2_bf*1.e3 << " & " << dm2_cred1sig_str << " & " << dm2_cred2sig_str  << " & " << dm2_cred3sig_str << "\\\\" << std::endl;
//  std::cout << "$\\Delta m^{2}_{32}$ [$\\times10^{-3}$ eV$^{2}$] &" << dm2_bf*1000. << " & [" << dm2_lo*1000. << "," << dm2_up*1000. << "] \\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
    std::cout << "\\end{tabular}" << std::endl;
    std::cout << "\\caption{Central values and $1\\sigma$ credible intervals from 1D marginalized posterior distributions, with a Gaussian $\\sin^{2}\\theta_{13}$ prior, for the specified oscillation parameters.}" << std::endl;
    std::cout << "\\label{table:credIntervals}" << std::endl;
    std::cout << "\\end{table}" << std::endl;
    std::cout << std::endl <<  std::endl << std::endl << std::endl;
  }

  TCanvas* c = new TCanvas("c","c",1200,1200);
  c->cd(); c->Draw();
  TPad* p1 = new TPad("p1","p1",0.0,0.0,0.5,0.5);
  TPad* p2 = new TPad("p2","p2",0.5,0.0,1.0,0.5);
  TPad* p3 = new TPad("p3","p3",0.0,0.5,0.5,1.0);
  TPad* p4 = new TPad("p4","p4",0.5,0.5,1.0,1.0);
  p1->Draw(); p2->Draw(); p3->Draw(); p4->Draw();

  // cosmetics
  gStyle->SetOptStat(0);
  gPad->SetTickx();
  gPad->SetTicky();
  int onesigcolor = 17;
  int twosigcolor = 15;
  int threesigcolor = 13;
  int outcolor = 10;

  delta_cp_hist->SetFillColor(onesigcolor);
  th23_hist->SetFillColor(onesigcolor);
  th13_hist->SetFillColor(onesigcolor);
  dm2_hist->SetFillColor(onesigcolor);

  delta_cp_cred->SetFillColor(twosigcolor);
  th23_cred->SetFillColor(twosigcolor);
  th13_cred->SetFillColor(twosigcolor);
  dm2_cred->SetFillColor(twosigcolor);

  delta_cp_cred2sig->SetFillColor(threesigcolor);
  th23_cred2sig->SetFillColor(threesigcolor);
  th13_cred2sig->SetFillColor(threesigcolor);
  dm2_cred2sig->SetFillColor(threesigcolor);

  delta_cp_cred3sig->SetFillColor(outcolor);
  th23_cred3sig->SetFillColor(outcolor);
  th13_cred3sig->SetFillColor(outcolor);
  dm2_cred3sig->SetFillColor(outcolor);

  delta_cp_cred->SetLineColor(kBlack);
  th23_cred->SetLineColor(kBlack);
  th13_cred->SetLineColor(kBlack);
  dm2_cred->SetLineColor(kBlack);

  delta_cp_cred2sig->SetLineColor(kBlack);
  th23_cred2sig->SetLineColor(kBlack);
  th13_cred2sig->SetLineColor(kBlack);
  dm2_cred2sig->SetLineColor(kBlack);

  delta_cp_cred3sig->SetLineColor(kBlack);
  th23_cred3sig->SetLineColor(kBlack);
  th13_cred3sig->SetLineColor(kBlack);
  dm2_cred3sig->SetLineColor(kBlack);

  delta_cp_hist->SetLineColor(kBlack);
  th23_hist->SetLineColor(kBlack);
  th13_hist->SetLineColor(kBlack);
  dm2_hist->SetLineColor(kBlack);
  dm2_hist->GetXaxis()->SetMaxDigits(2);

  TLegend* leg = new TLegend(0.25,0.7,0.9,0.87);
  TLegend* leg1 = new TLegend(0.05,0.7,0.9,0.87);
  if(levSigs) {
    leg->AddEntry(delta_cp_hist,"1#sigma credible interval","f");
    leg->AddEntry(delta_cp_cred,"2#sigma credible interval","f");
    leg->AddEntry(delta_cp_cred2sig,"3#sigma credible interval","f");
    leg->SetFillStyle(0);                                                                       
    leg1->AddEntry(delta_cp_hist,"1#sigma credible interval","f");
    leg1->AddEntry(delta_cp_cred,"2#sigma credible interval","f");
    leg1->AddEntry(delta_cp_cred2sig,"3#sigma credible interval","f");
    leg1->SetFillStyle(0);
  }
  else {
    leg->AddEntry(delta_cp_hist,"68% credible interval","f");
    leg->AddEntry(delta_cp_cred,"90% credible interval","f");
    leg->AddEntry(delta_cp_cred2sig,"99% credible interval","f");
    leg->SetFillStyle(0);                                                                       
    leg1->AddEntry(delta_cp_hist,"68% credible interval","f");
    leg1->AddEntry(delta_cp_cred,"90% credible interval","f");
    leg1->AddEntry(delta_cp_cred2sig,"99% credible interval","f");
    leg1->SetFillStyle(0);
  }

  leg->SetLineColorAlpha(0,0);
  leg->SetTextFont(132);
  leg->SetTextSize(0.05);

  leg1->SetLineColorAlpha(0,0);
  leg1->SetTextFont(132);
  leg1->SetTextSize(0.07);
  gStyle->SetLegendBorderSize(0);

  delta_cp_hist->GetYaxis()->SetRangeUser(0.,delta_cp_hist->GetMaximum()*1.5);
  th13_hist->GetYaxis()->SetRangeUser(0.,th13_hist->GetMaximum()*1.5);
  th23_hist->GetYaxis()->SetRangeUser(0.,th23_hist->GetMaximum()*1.5);
  if(hierarchy != 0) dm2_hist->GetYaxis()->SetRangeUser(0.,dm2_hist->GetMaximum()*1.5);
  TPaveText* hlab = new TPaveText(0.75,0.8,0.87,0.9,"NDC");
  if(hierarchy == 1) {
    dm2_hist->GetXaxis()->SetRangeUser(0.00225,0.00275);
    hlab->AddText("NH only");
  }

  if(hierarchy == -1) {
    dm2_hist->GetXaxis()->SetRangeUser(-0.00275,-0.00225);
    hlab->AddText("IH only");
  }

  hlab->SetFillStyle(0);
  hlab->SetFillStyle(0);
  hlab->SetTextFont(132);
  hlab->SetBorderSize(0);

  p1->cd(); delta_cp_hist->Draw("hist"); delta_cp_cred->Draw("hist same"); leg->Draw();
  p2->cd(); th23_hist->Draw("hist"); th23_cred->Draw("hist same"); leg->Draw(); 
  p3->cd(); th13_hist->Draw("hist"); th13_cred->Draw("hist same"); leg->Draw(); 
  p4->cd(); dm2_hist->Draw("hist"); dm2_cred->Draw("hist same"); leg->Draw(); 
  //dm2_cred->GetXaxis()->SetRangeUser(0.002,0.003);

  //if(hierarchy == 1)  c->Print("credInt_debugHists_NH.pdf");
  //else if(hierarchy == 0)  c->Print("credInt_debugHists_both.pdf");
  //else if(hierarchy == -1) c->Print("credInt_debugHists_IH.pdf");


  TCanvas* delta_cp_c = new TCanvas("delta_cp_c","delta_cp_c",600,500);
  TCanvas* th23_c = new TCanvas("th23_c","th23_c",600,500);
  TCanvas* th13_c = new TCanvas("th13_c","th13_c",600,500);
  TCanvas* dm2_c = new TCanvas("dm2_c","dm2_c",600,500);

  delta_cp_c->SetBottomMargin(0.15);
  th23_c->SetBottomMargin(0.15);
  th13_c->SetBottomMargin(0.15);
  dm2_c->SetBottomMargin(0.15);

  if(LogY) {
   delta_cp_c->SetLogy();
   th23_c->SetLogy();
   th13_c->SetLogy();
   dm2_c->SetLogy();
  }

  delta_cp_c->Draw(); delta_cp_c->cd();
  delta_cp_hist->Draw("hist");
  delta_cp_cred->Draw("hist same");
  delta_cp_cred2sig->Draw("hist same");
  delta_cp_cred3sig->Draw("hist same");
  leg->Draw();
  hlab->Draw();
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->RedrawAxis();

  th23_c->Draw(); th23_c->cd();
  th23_hist->Draw("hist");
  th23_cred->Draw("hist same");
  th23_cred2sig->Draw("hist same");
  th23_cred3sig->Draw("hist same");
  leg->Draw();
  hlab->Draw();
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->RedrawAxis();

  th13_c->Draw(); th13_c->cd();
  th13_hist->Draw("hist");
  th13_cred->Draw("hist same");
  th13_cred2sig->Draw("hist same");
  th13_cred3sig->Draw("hist same");
  leg->Draw();
  hlab->Draw();
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->RedrawAxis();

  if(hierarchy == 0) {
    dm2_c->Draw(); dm2_c->cd();
    TPad *p_NH = new TPad("p_NH","p_NH",0.5,0.0,1.0,1.0);
    TPad *p_IH = new TPad("p_IH","p_IH",0.0,0.0,0.5,1.0);
    p_NH->SetLeftMargin(0.01);
    p_NH->SetRightMargin(0.18);
    p_IH->SetRightMargin(0.01);
    p_IH->SetLeftMargin(0.18);
    p_NH->SetBottomMargin(0.15);
    p_IH->SetBottomMargin(0.15);
    p_NH->Draw(); p_IH->Draw();
    TH1D* dummyNH = (TH1D*)dm2_hist->Clone("dummyNH");
    TH1D* dummyIH = (TH1D*)dm2_hist->Clone("dummyIH");
    dummyNH->Reset(); dummyIH->Reset();
    dummyNH->GetXaxis()->SetRangeUser(0.00225,0.00275);
    dummyIH->GetXaxis()->SetRangeUser(-0.00275,-0.00225);
    dummyNH->GetYaxis()->SetRangeUser(0.,dm2_hist->GetMaximum()*1.5);
    dummyIH->GetYaxis()->SetRangeUser(0.,dm2_hist->GetMaximum()*1.5);

    //dummyNH->GetYaxis()->SetLabelSize(0.00);
    dummyNH->GetXaxis()->SetLabelSize(0.1);
    dummyIH->GetXaxis()->SetLabelSize(0.1);
    //dummyNH->GetYaxis()->SetTitle("");
    dummyNH->GetXaxis()->SetTitleSize(0.1);
    dummyIH->GetXaxis()->SetTitleSize(0.1);
    dummyIH->GetYaxis()->SetTitleSize(0.1);
    dummyNH->GetXaxis()->SetTitleOffset(0.6);
    dummyIH->GetXaxis()->SetTitleOffset(1.5);
    dummyIH->GetYaxis()->SetTitleOffset(0.8);
    dummyNH->GetXaxis()->SetLabelOffset(-0.01);
    dummyIH->GetXaxis()->SetLabelOffset(-0.01);
    dummyIH->GetXaxis()->SetNdivisions(512);
    dummyNH->GetXaxis()->SetNdivisions(512);
    TGaxis::SetExponentOffset(0.01, 0.0, "x"); 
    //dummyIH->GetXaxis()->SetLabelSize(0.06);
    //dummyIH->GetYaxis()->SetLabelSize(0.06);
    //dummyIH->GetYaxis()->SetTitle("");
    //dummyIH->SetTitle("");

    if(LogY) {
      //dm2_hist->GetYaxis()->SetRangeUser(100.,dm2_hist->GetMaximum()*100.);
      dummyIH->GetYaxis()->SetRangeUser(100.,dm2_hist->GetMaximum()*100.);
      dummyNH->GetYaxis()->SetRangeUser(100.,dm2_hist->GetMaximum()*100.);
      p_NH->SetLogy();
      p_IH->SetLogy();
    }

    p_NH->cd();
    dummyNH->Draw();
    dm2_hist->Draw("hist same");
    dm2_cred->Draw("hist same");
    dm2_cred2sig->Draw("hist same");
    dm2_cred3sig->Draw("hist same");
    leg1->Draw();
    gPad->SetTickx();
    gPad->SetTicky();
    gPad->RedrawAxis();
     
    p_IH->cd();
    dummyIH->Draw();
    dm2_hist->Draw("hist same");
    dm2_cred->Draw("hist same");
    dm2_cred2sig->Draw("hist same");
    dm2_cred3sig->Draw("hist same");
    gPad->SetTickx();
    gPad->SetTicky();
    gPad->RedrawAxis();

  }
  else {
    dm2_c->Draw(); dm2_c->cd();
    dm2_hist->Draw("hist");
    dm2_cred->Draw("hist same");
    dm2_cred2sig->Draw("hist same");
    dm2_cred3sig->Draw("hist same");
    leg->Draw();
    hlab->Draw();
    gPad->SetTickx();
    gPad->SetTicky();
    gPad->RedrawAxis();
  }

  if(LogY) {
    if(doRC) {
      delta_cp_hist->GetYaxis()->SetRangeUser(500.,delta_cp_hist->GetMaximum()*10.);
      th23_hist->GetYaxis()->SetRangeUser(100.,th23_hist->GetMaximum()*10.);
      th13_hist->GetYaxis()->SetRangeUser(100.,th13_hist->GetMaximum()*10.);
    }
    else {
      delta_cp_hist->GetYaxis()->SetRangeUser(50000.,delta_cp_hist->GetMaximum()*10.);
      th23_hist->GetYaxis()->SetRangeUser(10000.,th23_hist->GetMaximum()*10.);
      th13_hist->GetYaxis()->SetRangeUser(10000.,th13_hist->GetMaximum()*10.);
    }
    delta_cp_c->cd(); gPad->RedrawAxis();
    th23_c->cd(); gPad->RedrawAxis();
    th13_c->cd(); gPad->RedrawAxis();
    if(hierarchy != 0) {
      dm2_hist->GetYaxis()->SetRangeUser(100.,dm2_hist->GetMaximum()*100.);
      dm2_c->cd(); gPad->RedrawAxis();
    }
  }

  if(levSigs) {
      if(hierarchy == 1) {
        if(LogY) {
          delta_cp_c->Print("ContoursJoint1D/creds_1D_delta_cp_NH_LogY.pdf");
          th23_c->Print("ContoursJoint1D/creds_1D_th23_NH_LogY.pdf");
          th13_c->Print("ContoursJoint1D/creds_1D_th13_NH_LogY.pdf");
          dm2_c->Print("ContoursJoint1D/creds_1D_dm2_NH_LogY.pdf"); 
        }
        else {
          delta_cp_c->Print("ContoursJoint1D/creds_1D_delta_cp_NH.pdf");
          th23_c->Print("ContoursJoint1D/creds_1D_th23_NH.pdf");
          th13_c->Print("ContoursJoint1D/creds_1D_th13_NH.pdf");
          dm2_c->Print("ContoursJoint1D/creds_1D_dm2_NH.pdf");
        }
      }
      if(hierarchy == 0) {
        if(LogY) {
          delta_cp_c->Print("ContoursJoint1D/creds_1D_delta_cp_both_LogY.pdf");
          th23_c->Print("ContoursJoint1D/creds_1D_th23_both_LogY.pdf");
          th13_c->Print("ContoursJoint1D/creds_1D_th13_both_LogY.pdf");
          dm2_c->Print("ContoursJoint1D/creds_1D_dm2_both_LogY.pdf"); 
        }
        else {
          delta_cp_c->Print("ContoursJoint1D/creds_1D_delta_cp_both.pdf");
          th23_c->Print("ContoursJoint1D/creds_1D_th23_both.pdf");
          th13_c->Print("ContoursJoint1D/creds_1D_th13_both.pdf");
          dm2_c->Print("ContoursJoint1D/creds_1D_dm2_both.pdf");
        }
      }
      if(hierarchy == -1) {
        if(LogY) {
          delta_cp_c->Print("ContoursJoint1D/creds_1D_delta_cp_IH_LogY.pdf");
          th23_c->Print("ContoursJoint1D/creds_1D_th23_IH_LogY.pdf");
          th13_c->Print("ContoursJoint1D/creds_1D_th13_IH_LogY.pdf");
          dm2_c->Print("ContoursJoint1D/creds_1D_dm2_IH_LogY.pdf"); 
        }
        else {
          delta_cp_c->Print("ContoursJoint1D/creds_1D_delta_cp_IH.pdf");
          th23_c->Print("ContoursJoint1D/creds_1D_th23_IH.pdf");
          th13_c->Print("ContoursJoint1D/creds_1D_th13_IH.pdf");
          dm2_c->Print("ContoursJoint1D/creds_1D_dm2_IH.pdf");
        }
      }
  }
  else {
      if(hierarchy == 1) {
        if(LogY) {
          delta_cp_c->Print("ContoursJoint1D/creds_1D_delta_cp_NH_LogY_perc.pdf");
          th23_c->Print("ContoursJoint1D/creds_1D_th23_NH_LogY_perc.pdf");
          th13_c->Print("ContoursJoint1D/creds_1D_th13_NH_LogY_perc.pdf");
          dm2_c->Print("ContoursJoint1D/creds_1D_dm2_NH_LogY_perc.pdf"); 
        }
        else {
          delta_cp_c->Print("ContoursJoint1D/creds_1D_delta_cp_NH_perc.pdf");
          th23_c->Print("ContoursJoint1D/creds_1D_th23_NH_perc.pdf");
          th13_c->Print("ContoursJoint1D/creds_1D_th13_NH_perc.pdf");
          dm2_c->Print("ContoursJoint1D/creds_1D_dm2_NH_perc.pdf");
        }
      }
      if(hierarchy == 0) {
        if(LogY) {
          delta_cp_c->Print("ContoursJoint1D/creds_1D_delta_cp_both_LogY_perc.pdf");
          th23_c->Print("ContoursJoint1D/creds_1D_th23_both_LogY_perc.pdf");
          th13_c->Print("ContoursJoint1D/creds_1D_th13_both_LogY_perc.pdf");
          dm2_c->Print("ContoursJoint1D/creds_1D_dm2_both_LogY_perc.pdf"); 
        }
        else {
          delta_cp_c->Print("ContoursJoint1D/creds_1D_delta_cp_both_perc.pdf");
          th23_c->Print("ContoursJoint1D/creds_1D_th23_both_perc.pdf");
          th13_c->Print("ContoursJoint1D/creds_1D_th13_both_perc.pdf");
          dm2_c->Print("ContoursJoint1D/creds_1D_dm2_both_perc.pdf");
          //dm2_c->Print("thContours/creds_1D_dm2_both.root");
        }
      }
      if(hierarchy == -1) {
        if(LogY) {
          delta_cp_c->Print("ContoursJoint1D/creds_1D_delta_cp_IH_LogY_perc.pdf");
          th23_c->Print("ContoursJoint1D/creds_1D_th23_IH_LogY_perc.pdf");
          th13_c->Print("ContoursJoint1D/creds_1D_th13_IH_LogY_perc.pdf");
          dm2_c->Print("ContoursJoint1D/creds_1D_dm2_IH_LogY_perc.pdf"); 
        }
        else {
          delta_cp_c->Print("ContoursJoint1D/creds_1D_delta_cp_IH_perc.pdf");
          th23_c->Print("ContoursJoint1D/creds_1D_th23_IH_perc.pdf");
          th13_c->Print("ContoursJoint1D/creds_1D_th13_IH_perc.pdf");
          dm2_c->Print("ContoursJoint1D/creds_1D_dm2_IH_perc.pdf");
        }
      }
  }

  delete t;
  delete f;
  delete delta_cp_c;
  delete th23_c;
  delete th13_c;
  delete dm2_c;

}

void make1DContours(TString ReducedChain) {
//void makePlot(int hierarchy, bool doRC = false, bool LogY = false, bool levSigs = true)
  // 1,2,3 sigma
  //makePlot(ReducedChain,0,true,true);
  makePlot(ReducedChain,1,false,false);
  //makePlot(ReducedChain,0,true,false);
  //makePlot(ReducedChain,0,false,true);

  // 68%, 90%, 99%
  //makePlot(ReducedChain,0,false,false,false);
  //makePlot(ReducedChain,0,true,false,false);

  //makePlot(ReducedChain,0,false,true,false);
  //makePlot(ReducedChain,0,true,true,false);
}
