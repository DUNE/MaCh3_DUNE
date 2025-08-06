#ifndef _StructsDUNE_h_
#define _StructsDUNE_h_

struct dunemc_base {

  /*
  int nEvents; // how many MC events are there
  double osc_channel;
  */
  double pot_s;
  double norm_s;
  
  int Target; //Target the interaction was on
  
  int nupdg;
  int nupdgUnosc;
  double OscChannelIndex;
  
  double rw_lovere;
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
  /*
  int rw_nuPDGunosc;
  int rw_nuPDG;
  int rw_run;
  bool    *rw_isFHC;
  */
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

  double rw_lep_pT;  //transverse lepton momentum
  double rw_lep_pX;
  double rw_lep_pY;
  double rw_lep_pZ; //parallel lepton momentum
  double rw_reco_vtx_x;
  double rw_reco_vtx_y;
  double rw_reco_vtx_z;
  double rw_reco_rad;
  double rw_rad;

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
  /*

  int nrecoparticles;
  bool *in_fdv;
  bool *is_accepted;
  bool *is_good_caf_event;

  //Particle-level kinematic parameters (JM for NDGAr)
  double *particle_ecaldepositfraction;
  int *particle_event;
  int *particle_pdg;
  double *particle_energy;
  double *particle_theta;
  double *particle_bangle;
  double *particle_dedx;
  double *particle_momentum;
  double *particle_transversemomentum;
  bool *particle_isaccepted;
  bool *particle_isstoppedintpc;
  bool *particle_isstoppedinecal;
  bool *particle_isstoppedingap;
  bool *particle_isstoppedinbarrelgap;
  bool *particle_isstoppedinendgap;
  bool *particle_isstoppedinbarrel;
  bool *particle_isstoppedinendcap;
  double *particle_startx;
  double *particle_startr2;
  double *particle_nhits;
  double *particle_nturns;
  double *particle_momresms;
  double *particle_momrestransfrac;
  double *particle_momrestrans;
  */
};

#endif
