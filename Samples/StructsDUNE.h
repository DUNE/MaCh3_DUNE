#ifndef _StructsDUNE_h_
#define _StructsDUNE_h_

#include "TMatrixD.h"

// HH: A little struct to store the ND covariance matrices together, to avoid having to pass multiple TMatrixD* around in the code. 
// Also includes a bool to specify whether to use a combined ND covariance matrix or separate ones for FHC and RHC
struct BeamNDCov{
  TMatrixD* NDCov_FHC;
  TMatrixD* NDCov_RHC;
  TMatrixD* NDCov_all;
  bool useCombinedNDCov;
};

struct BeamSampleInfo {
  bool iselike;
  double isFHC;
  double pot;
};

struct BeamFDSampleInfo : public BeamSampleInfo {
};

struct BeamNDSampleInfo : public BeamSampleInfo {
  double pot_s;
  double norm_s;
};

struct dunemc_base { // Store variables used in fitting

  double pot_s;
  double norm_s;
  
  double Target; //Target the interaction was on
  
  int nupdg;
  int nupdgUnosc;
  int rw_isCC;
  double OscChannelIndex;
  
  double rw_erec;
  double enu_true;

  double flux_w;
  double mode;

  // Analysis sample
  unsigned int SampleIndex = M3::_BAD_INT_;
};

struct dunemc_atm : public dunemc_base { // Store variables used by SampleHandlerAtm
  double rw_theta;
  double coszenith_true;
};

struct dunemc_beamfd : public dunemc_base { // Store variables used by SampleHandlerBeamFD
  double rw_erec_shifted;
  double rw_erec_had;
  double rw_erec_lep;
  double rw_eRecoP;
  double rw_eRecoPip;
  double rw_eRecoPim;
  double rw_eRecoPi0;
  double rw_eRecoN;
  double rw_LepE;
  double rw_eP;
  double rw_ePip;
  double rw_ePim;
  double rw_ePi0;
  double rw_eN;

  double rw_berpaacvwgt;

  double rw_cvnnumu;
  double rw_cvnnue;
  double rw_cvnnumu_shifted;
  double rw_cvnnue_shifted;

  double rw_vtx_x;
  double rw_vtx_y;
  double rw_vtx_z;

  double rw_erec_had_sqrt;
  double rw_erec_lep_sqrt;
  double rw_eRecoPi0_sqrt;
  double rw_eRecoN_sqrt;
  double rw_sum_ehad;
  double rw_sum_ehad_sqrt;
  double rw_trueccnue;
  double rw_trueccnumu;
};

struct dunemc_beamnd : public dunemc_base { // Store variables used by SampleHandlerBeamND
  double rw_erec_shifted;
  double rw_erec_had;
  double rw_erec_lep;
  double rw_yrec;
  double rw_theta;
  double rw_berpaacvwgt;
  double rw_reco_q;
};

struct dunemc_beamndgar : public dunemc_base { // Store variables used by SampleHandlerBeamNDGAr
  double rw_LepE;
  double rw_berpaacvwgt;
  double rw_vtx_x;
  double rw_vtx_y;
  double rw_vtx_z;
  double rw_Q0;
  double rw_Q3;
  double rw_lep_pT;  //transverse lepton momentum
  double rw_lep_pX;
  double rw_lep_pY;
  double rw_lep_pZ; //parallel lepton momentum
  double rw_rad;
};

struct dunemc_plotting { // Store variables just used in plotting (cleared from memory before a fit)
  bool in_fdv;
  bool is_accepted;
  double geometric_correction;
  double rw_lep_p;
  double rw_lep_phi;
  double rw_lep_theta;
  double rw_lep_bangle;
  double lep_tracklengthyz;
  int npi0;
  double rw_ePi0;
 
  std::vector<int> prim_event = {};
  std::vector<int> prim_trkid = {};
  std::vector<int> prim_pdg = {};
  std::vector<double> prim_evis = {};
  std::vector<double> prim_theta = {};
  std::vector<double> prim_bangle = {};
  std::vector<double> prim_beamangle = {};
  std::vector<double> prim_dedx = {};
  std::vector<double> prim_momentum = {};
  std::vector<double> prim_endmomentum = {};
  std::vector<double> prim_transversemomentum = {};
  std::vector<int> prim_isaccepted = {};
  std::vector<int> prim_iscurvatureresolved = {};
  std::vector<int> prim_isdecayed = {};
  std::vector<int> prim_isstoppedintpc = {};
  std::vector<int> prim_isstoppedinecal = {};
  std::vector<int> prim_isstoppedingap = {};
  std::vector<int> prim_isstoppedinbarrelgap = {};
  std::vector<int> prim_isstoppedinendgap = {};
  std::vector<int> prim_isstoppedinbarrel = {};
  std::vector<int> prim_isstoppedinendcap = {};
  std::vector<int> prim_isescaped = {};
  std::vector<double> prim_startx = {};
  std::vector<double> prim_startr2 = {};
  std::vector<double> prim_endr = {};
  std::vector<double> prim_enddepth = {};
  std::vector<double> prim_endx = {};
  std::vector<double> prim_endy = {};
  std::vector<double> prim_endz = {};
  std::vector<double> prim_nturns = {};
  std::vector<double> prim_nhits = {};
  std::vector<double> prim_tracklengthyz = {};
  std::vector<double> prim_momresms = {};
  std::vector<double> prim_momresyz = {};
  std::vector<double> prim_momresx = {};
  std::vector<double> prim_edepcrit = {};
  std::vector<double> prim_tpcedepfrac = {};
  std::vector<double> prim_iscontained = {};

  std::vector<double> shower_dcalboundary = {};
  std::vector<double> shower_pdg = {};
  std::vector<double> shower_energy = {};
  std::vector<double> shower_bangle = {};
  std::vector<double> shower_iscontained = {};
  std::vector<double> shower_cosnorm = {};

  std::vector<double> photon_energy = {};
  std::vector<double> photon_endx = {};
  std::vector<double> photon_endy = {};
  std::vector<double> photon_endz = {};
};

#endif
