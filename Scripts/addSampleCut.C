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

//DUNE ND CC inclusive Selection Cut
inline bool IsCCInclusive(int reco_numu, int muon_contained, int muon_tracker, double Ehad_veto) 
{
  return (reco_numu && (muon_contained || muon_tracker) && Ehad_veto < 30);
}



//DUNE ND FV cut
inline bool IsInNDFV(double pos_x_cm, double pos_y_cm, double pos_z_cm) 
{
  bool inDeadRegion = false;
  for (int i = -3; i <= 3; ++i) {
    // 0.5cm cathode in the middle of each module, plus 0.5cm buffer
    double cathode_center = i * 102.1;
    if (pos_x_cm > cathode_center - 0.75 && pos_x_cm < cathode_center + 0.75) {
      inDeadRegion = true;
	}

    // 1.6cm dead region between modules (0.5cm module wall and 0.3cm pixel
    // plane, x2) don't worry about outer boundary because events are only
    // generated in active Ar + insides
    double module_boundary = i * 102.1 + 51.05;
    if (i <= 2 && pos_x_cm > module_boundary - 1.3 && pos_x_cm < module_boundary + 1.3) {
       inDeadRegion = true; 
	}
  }
  for (int i = 1; i <= 4; ++i) {
    // module boundaries in z are 1.8cm (0.4cm ArCLight plane + 0.5cm module
    // wall, x2) module is 102.1cm wide, but only 101.8cm long due to cathode
    // (0.5cm) being absent in length but ArCLight is 0.1cm thicker than pixel
    // plane so it's 0.3cm difference positions are off-set by 0.6 because I
    // defined 0 to be the upstream edge based on the active volume by
    // inspecting a plot, and aparently missed by 3 mm, but whatever add 8mm =
    // 2 pad buffer due to worse position resolution in spatial dimension z
    // compared to timing direction x so total FV gap will be 1.8 + 2*0.8
    // = 3.4cm
    double module_boundary = i * 101.8 - 0.6;
    if (pos_z_cm > module_boundary - 1.7 && pos_z_cm < module_boundary + 1.7) {
	  inDeadRegion = true;
    }
  }
    
  return (abs(pos_x_cm) < 200 && abs(pos_y_cm) < 100 && pos_z_cm > 50 &&
		              pos_z_cm < 350 && !inDeadRegion);
}


void addSampleCut(std::string infile){

  gROOT->ProcessLine( "gErrorIgnoreLevel = 6001;");
  gROOT->ProcessLine( "gPrintViaErrorHandler = kTRUE;");

  int mu_contained;
  int  mu_tracker;
  double ehad_veto;
  int reco_numu;
  double vtx_x;
  double vtx_y;
  double vtx_z;
  int fv_cc_cut;

  TFile *f = new TFile(infile.c_str(), "update");

  TTree * t = (TTree*)f->Get("caf");
  
  t->SetBranchAddress("muon_contained", &mu_contained);
  t->SetBranchAddress("muon_tracker", &mu_tracker);
  t->SetBranchAddress("Ehad_veto", &ehad_veto);
  t->SetBranchAddress("reco_numu", &reco_numu);
  t->SetBranchAddress("vtx_x", &vtx_x);
  t->SetBranchAddress("vtx_y", &vtx_y);
  t->SetBranchAddress("vtx_z", &vtx_z);

  int N = t->GetEntries();

  TBranch *fv_ccinc_cut = t->Branch("fv_ccinc_cut", &fv_cc_cut, "fv_ccinc_cut/I");
  int j = 0;

  for (int i=0; i < N; i++) {
	t->GetEntry(i);

    if( i % 100000 == 0 ) printf( "Event %d of %d...\n", i, N );

    if(IsInNDFV(vtx_x, vtx_y, vtx_z) && IsCCInclusive(reco_numu, mu_contained, mu_tracker, ehad_veto))
	{
	  fv_cc_cut = 1;
	  j++;
	}
	else
	{
	  fv_cc_cut = 0;
	}

    fv_ccinc_cut->Fill();

  }

  std::cout << "Finished writing to " << infile << std::endl;
  std::cout << "Events that passed the cut: " << j << "/" << N << std::endl;
  std::cout << "Works out to:" << j*100/N << "%" << std::endl;
  f->cd();
  t->Write("", TObject::kOverwrite);
  f->Close();
  delete f;

  return;
}
