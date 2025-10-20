#ifndef _StructsDUNE_h_
#define _StructsDUNE_h_

struct dunemc_base { // Store variables used in fitting

  int nEvents; // how many MC events are there
  double osc_channel;
  double pot_s;
  double norm_s;
  
  int Target; //Target the interaction was on
  
  int nupdg;
  int nupdgUnosc;
  double OscChannelIndex;
  
  double rw_erec;
  double rw_erec_shifted;
  double rw_erec_had;
  double rw_erec_lep;
  double rw_yrec;
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
  
  double rw_etru;
  double rw_mom;
  double rw_theta;
  //double rw_Q2;

  double rw_cvnnumu;
  double rw_cvnnue;
  double rw_cvnnumu_shifted;
  double rw_cvnnue_shifted;
  int rw_reco_nue;
  int rw_reco_numu;
  double rw_berpaacvwgt;
  int rw_isCC;

  double rw_vtx_x;
  double rw_vtx_y;
  double rw_vtx_z;
  double rw_reco_q;
  double reco_numu;
  double rw_Q0;
  double rw_Q3;

  //double beam_w;
  double flux_w;

  double mode;
  //int isbound;

  double rw_truecz;

  /*
  int nproton; ///< number of (post-FSI) primary protons
  int nneutron; ///< number of (post-FSI) primary neutrons
  int npip; ///< number of (post-FSI) primary pi+
  int npim; ///< number of (post-FSI) primary pi-
  int npi0; ///< number of (post-FSI) primary pi0
  int ntruemuon; //number of true muons
  int ntruemuonprim; //number of true primary muons
  int nrecomuon; //number of reconstructed muons
  double nmuonsratio; //number of reco muons divided by number of true muons
  */

  double rw_lep_pT;  //transverse lepton momentum
  double rw_lep_pX;
  double rw_lep_pY;
  double rw_lep_pZ; //parallel lepton momentum
  double rw_rad;

  /*
  double rw_reco_vtx_x;
  double rw_reco_vtx_y;
  double rw_reco_vtx_z;
  double rw_reco_rad;
  
  double rw_elep_reco;
  double rw_elep_true;
*/

  double rw_erec_had_sqrt;
  double rw_erec_lep_sqrt;
  double rw_eRecoPi0_sqrt;
  double rw_eRecoN_sqrt;
  double rw_sum_ehad;
  double rw_sum_ehad_sqrt;
  double rw_trueccnue;
  double rw_trueccnumu;
};

struct dunemc_plotting { // Store variables just used in plotting (cleared from memory before a fit)

  bool in_fdv;
  bool is_accepted;
  double geometric_correction;
  double rw_lep_p;
  double rw_lep_phi;
  double rw_lep_theta;
  double rw_lep_bangle;

  std::vector<int> particle_event = {};
  std::vector<int> particle_trkid = {};
  std::vector<int> particle_pdg = {};
  std::vector<double> particle_energy = {};
  std::vector<double> particle_theta = {};
  std::vector<double> particle_bangle = {};
  std::vector<double> particle_beamangle = {};
  std::vector<double> particle_dedx = {};
  std::vector<double> particle_momentum = {};
  std::vector<double> particle_endmomentum = {};
  std::vector<double> particle_transversemomentum = {};
  std::vector<int> particle_isaccepted = {};
  std::vector<int> particle_iscurvatureresolved = {};
  std::vector<int> particle_isdecayed = {};
  std::vector<int> particle_isstoppedintpc = {};
  std::vector<int> particle_isstoppedinecal = {};
  std::vector<int> particle_isstoppedingap = {};
  std::vector<int> particle_isstoppedinbarrelgap = {};
  std::vector<int> particle_isstoppedinendgap = {};
  std::vector<int> particle_isstoppedinbarrel = {};
  std::vector<int> particle_isstoppedinendcap = {};
  std::vector<int> particle_isescaped = {};
  std::vector<double> particle_startx = {};
  std::vector<double> particle_startr2 = {};
  std::vector<double> particle_endr = {};
  std::vector<double> particle_enddepth = {};
  std::vector<double> particle_endx = {};
  std::vector<double> particle_endy = {};
  std::vector<double> particle_endz = {};
  std::vector<double> particle_nturns = {};
  std::vector<double> particle_nhits = {};
  std::vector<double> particle_tracklengthyz = {};
  std::vector<double> particle_momresms = {};
  std::vector<double> particle_momresyz = {};
  std::vector<double> particle_momresx = {};
  std::vector<double> particle_edepcrit = {};
};

#endif
