#ifndef _StructsDUNE_h_
#define _StructsDUNE_h_

struct dunemc_base {

  // General Event Information
  int nEvents; // how many MC events are there
  int *Target; // Target the interaction was on

  // Neutrino Information
  int *nupdg;
  int *nupdgUnosc;

  // Reconstructed Energy Variables
  double *rw_erec;
  double *rw_erec_shifted;
  double *rw_erec_had;
  double *rw_erec_lep;
  double *rw_yrec;

  // Reconstructed Particle Energies
  double *rw_eRecoP;
  double *rw_eRecoPip;
  double *rw_eRecoPim;
  double *rw_eRecoPi0;
  double *rw_eRecoN;

  // True Particle Energies
  double *rw_LepE;
  double *rw_eP;
  double *rw_ePip;
  double *rw_ePim;
  double *rw_ePi0;
  double *rw_eN;

  // Kinematic Variables
  double *rw_etru;
  double *rw_mom;
  double *rw_theta;
  double *rw_Q2;

  double *rw_reco_theta;

  // CVN Scores
  double *rw_cvnnumu;
  double *rw_cvnnue;
  double *rw_cvnnumu_shifted;
  double *rw_cvnnue_shifted;

  // Reconstructed Event Information
  int *rw_reco_nue;
  int *rw_reco_numu;
  double *rw_berpaacvwgt;
  int *rw_isCC;
  int *rw_nuPDGunosc;
  int *rw_nuPDG;
  int *rw_run;
  bool *rw_isFHC;

  // Vertex Information
  double *rw_vtx_x;
  double *rw_vtx_y;
  double *rw_vtx_z;
  double dummy_y;
  double *rw_reco_q;

  // Muon Information
  double *reco_numu;
  double *rw_Q0;
  double *rw_Q3;

  // Scaling and Normalization
  double pot_s;
  double norm_s;
  double osc_channel;
  double *beam_w;
  double *flux_w;

  // Interaction Mode and Binding
  double *mode;
  int *isbound;

  // True Cosine of Zenith Angle
  double *rw_truecz;

  // Particle Counts
  int *nproton; ///< number of (post-FSI) primary protons
  int *nneutron; ///< number of (post-FSI) primary neutrons
  int *npip; ///< number of (post-FSI) primary pi+
  int *npim; ///< number of (post-FSI) primary pi-
  int *npi0; ///< number of (post-FSI) primary pi0

  // Muon Counts and Ratios
  int *ntruemuon; // number of true muons
  int *ntruemuonprim; // number of true primary muons
  int *nrecomuon; // number of reconstructed muons
  double *nmuonsratio; // number of reco muons divided by number of true muons

  // Lepton Momentum
  double *rw_lep_pT;  // transverse lepton momentum
  double *rw_lep_pX;
  double *rw_lep_pY;
  double *rw_lep_pZ; // parallel lepton momentum

  // Reconstructed Vertex Information
  double *rw_reco_vtx_x;
  double *rw_reco_vtx_y;
  double *rw_reco_vtx_z;
  double *rw_reco_rad;
  double *rw_rad;

  // Electron Energy
  double *rw_elep_reco;
  double *rw_elep_true;

  // ********************************
  // ND GAr Variables
  // ********************************
  int *nrecoparticles;
  bool *in_fdv;
  bool *is_accepted;
  bool *is_good_caf_event;

  // Particle-level Kinematic Parameters (ND GAr)
  std::vector<double> *particle_ecaldepositfraction;
  std::vector<int> *particle_event;
  std::vector<int> *particle_pdg;
  std::vector<double> *particle_energy;
  std::vector<double> *particle_theta;
  std::vector<double> *particle_bangle;
  std::vector<double> *particle_dedx;
  std::vector<double> *particle_momentum;
  std::vector<double> *particle_transversemomentum;
  std::vector<bool> *particle_isaccepted;
  std::vector<bool> *particle_isstoppedintpc;
  std::vector<bool> *particle_isstoppedinecal;
  std::vector<bool> *particle_isstoppedingap;
  std::vector<bool> *particle_isstoppedinbarrelgap;
  std::vector<bool> *particle_isstoppedinendgap;
  std::vector<bool> *particle_isstoppedinbarrel;
  std::vector<bool> *particle_isstoppedinendcap;
  std::vector<double> *particle_startx;
  std::vector<double> *particle_startr2;
  std::vector<double> *particle_nhits;
  std::vector<double> *particle_nturns;
  std::vector<double> *particle_momresms;
  std::vector<double> *particle_momrestransfrac;
  std::vector<double> *particle_momrestrans;
  int *nparticlesinsample;

  // ********************************
  // ND LAr Particle Variables
  // ********************************
  
  // ND LAr Variables
  std::vector<double> *rw_ndLAr_particle_eRecoMuon;
  std::vector<double> *rw_ndLAr_particle_eRecoP;
  std::vector<double> *rw_ndLAr_particle_eRecoPip;
  std::vector<double> *rw_ndLAr_particle_eRecoPim;
  std::vector<double> *rw_ndLAr_particle_eRecoPi0;
  std::vector<double> *rw_ndLAr_particle_eRecoN;
  std::vector<double> *rw_ndLAr_particle_eRecoElectron;
  std::vector<double> *rw_ndLAr_particle_eRecoGamma;

  std::vector<double> *rw_ndLAr_particle_eMuon;
  std::vector<double> *rw_ndLAr_particle_eP;
  std::vector<double> *rw_ndLAr_particle_ePip;
  std::vector<double> *rw_ndLAr_particle_ePim;
  std::vector<double> *rw_ndLAr_particle_ePi0;
  std::vector<double> *rw_ndLAr_particle_eN;
  std::vector<double> *rw_ndLAr_particle_eElectron;
  std::vector<double> *rw_ndLAr_particle_eGamma;

  std::vector<double> *rw_ndLAr_particle_MomRecoMu;
  std::vector<double> *rw_ndLAr_particle_MomRecoPip;
  std::vector<double> *rw_ndLAr_particle_MomRecoPim;
  std::vector<double> *rw_ndLAr_particle_MomRecoPi0;
  std::vector<double> *rw_ndLAr_particle_MomRecoP;
  std::vector<double> *rw_ndLAr_particle_MomRecoN;
  std::vector<double> *rw_ndLAr_particle_MomRecoElectron;
  std::vector<double> *rw_ndLAr_particle_MomRecoGamma;

  std::vector<double> *rw_ndLAr_particle_MomMu;
  std::vector<double> *rw_ndLAr_particle_MomPip;
  std::vector<double> *rw_ndLAr_particle_MomPim;
  std::vector<double> *rw_ndLAr_particle_MomPi0;
  std::vector<double> *rw_ndLAr_particle_MomP;
  std::vector<double> *rw_ndLAr_particle_MomN;
  std::vector<double> *rw_ndLAr_particle_MomElectron;
  std::vector<double> *rw_ndLAr_particle_MomGamma;

  // Reconstructed Vertex Start Positions
  std::vector<double> *rw_ndLAr_particle_StartXvtxRecoMu;
  std::vector<double> *rw_ndLAr_particle_StartYvtxRecoMu;
  std::vector<double> *rw_ndLAr_particle_StartZvtxRecoMu;
  std::vector<double> *rw_ndLAr_particle_StartXvtxRecoPip;
  std::vector<double> *rw_ndLAr_particle_StartYvtxRecoPip;
  std::vector<double> *rw_ndLAr_particle_StartZvtxRecoPip;
  std::vector<double> *rw_ndLAr_particle_StartXvtxRecoPim;
  std::vector<double> *rw_ndLAr_particle_StartYvtxRecoPim;
  std::vector<double> *rw_ndLAr_particle_StartZvtxRecoPim;
  std::vector<double> *rw_ndLAr_particle_StartXvtxRecoPi0;
  std::vector<double> *rw_ndLAr_particle_StartYvtxRecoPi0;
  std::vector<double> *rw_ndLAr_particle_StartZvtxRecoPi0;
  std::vector<double> *rw_ndLAr_particle_StartXvtxRecoP;
  std::vector<double> *rw_ndLAr_particle_StartYvtxRecoP;
  std::vector<double> *rw_ndLAr_particle_StartZvtxRecoP;
  std::vector<double> *rw_ndLAr_particle_StartXvtxRecoN;
  std::vector<double> *rw_ndLAr_particle_StartYvtxRecoN;
  std::vector<double> *rw_ndLAr_particle_StartZvtxRecoN;
  std::vector<double> *rw_ndLAr_particle_StartXvtxRecoElectron;
  std::vector<double> *rw_ndLAr_particle_StartYvtxRecoElectron;
  std::vector<double> *rw_ndLAr_particle_StartZvtxRecoElectron;
  std::vector<double> *rw_ndLAr_particle_StartXvtxRecoGamma;
  std::vector<double> *rw_ndLAr_particle_StartYvtxRecoGamma;
  std::vector<double> *rw_ndLAr_particle_StartZvtxRecoGamma;

  // True Vertex Start Positions
  std::vector<double> *rw_ndLAr_particle_StartXvtxMu;
  std::vector<double> *rw_ndLAr_particle_StartYvtxMu;
  std::vector<double> *rw_ndLAr_particle_StartZvtxMu;
  std::vector<double> *rw_ndLAr_particle_StartXvtxPip;
  std::vector<double> *rw_ndLAr_particle_StartYvtxPip;
  std::vector<double> *rw_ndLAr_particle_StartZvtxPip;
  std::vector<double> *rw_ndLAr_particle_StartXvtxPim;
  std::vector<double> *rw_ndLAr_particle_StartYvtxPim;
  std::vector<double> *rw_ndLAr_particle_StartZvtxPim;
  std::vector<double> *rw_ndLAr_particle_StartXvtxPi0;
  std::vector<double> *rw_ndLAr_particle_StartYvtxPi0;
  std::vector<double> *rw_ndLAr_particle_StartZvtxPi0;
  std::vector<double> *rw_ndLAr_particle_StartXvtxP;
  std::vector<double> *rw_ndLAr_particle_StartYvtxP;
  std::vector<double> *rw_ndLAr_particle_StartZvtxP;
  std::vector<double> *rw_ndLAr_particle_StartXvtxN;
  std::vector<double> *rw_ndLAr_particle_StartYvtxN;
  std::vector<double> *rw_ndLAr_particle_StartZvtxN;
  std::vector<double> *rw_ndLAr_particle_StartXvtxElectron;
  std::vector<double> *rw_ndLAr_particle_StartYvtxElectron;
  std::vector<double> *rw_ndLAr_particle_StartZvtxElectron;
  std::vector<double> *rw_ndLAr_particle_StartXvtxGamma;
  std::vector<double> *rw_ndLAr_particle_StartYvtxGamma;
  std::vector<double> *rw_ndLAr_particle_StartZvtxGamma;

  // Reconstructed Vertex End Positions
  std::vector<double> *rw_ndLAr_particle_EndXvtxRecoMu;
  std::vector<double> *rw_ndLAr_particle_EndYvtxRecoMu;
  std::vector<double> *rw_ndLAr_particle_EndZvtxRecoMu;
  std::vector<double> *rw_ndLAr_particle_EndXvtxRecoPip;
  std::vector<double> *rw_ndLAr_particle_EndYvtxRecoPip;
  std::vector<double> *rw_ndLAr_particle_EndZvtxRecoPip;
  std::vector<double> *rw_ndLAr_particle_EndXvtxRecoPim;
  std::vector<double> *rw_ndLAr_particle_EndYvtxRecoPim;
  std::vector<double> *rw_ndLAr_particle_EndZvtxRecoPim;
  std::vector<double> *rw_ndLAr_particle_EndXvtxRecoPi0;
  std::vector<double> *rw_ndLAr_particle_EndYvtxRecoPi0;
  std::vector<double> *rw_ndLAr_particle_EndZvtxRecoPi0;
  std::vector<double> *rw_ndLAr_particle_EndXvtxRecoP;
  std::vector<double> *rw_ndLAr_particle_EndYvtxRecoP;
  std::vector<double> *rw_ndLAr_particle_EndZvtxRecoP;
  std::vector<double> *rw_ndLAr_particle_EndXvtxRecoN;
  std::vector<double> *rw_ndLAr_particle_EndYvtxRecoN;
  std::vector<double> *rw_ndLAr_particle_EndZvtxRecoN;
  std::vector<double> *rw_ndLAr_particle_EndXvtxRecoElectron;
  std::vector<double> *rw_ndLAr_particle_EndYvtxRecoElectron;
  std::vector<double> *rw_ndLAr_particle_EndZvtxRecoElectron;
  std::vector<double> *rw_ndLAr_particle_EndXvtxRecoGamma;
  std::vector<double> *rw_ndLAr_particle_EndYvtxRecoGamma;
  std::vector<double> *rw_ndLAr_particle_EndZvtxRecoGamma;

  // True Vertex End Positions
  std::vector<double> *rw_ndLAr_particle_EndXvtxMu;
  std::vector<double> *rw_ndLAr_particle_EndYvtxMu;
  std::vector<double> *rw_ndLAr_particle_EndZvtxMu;
  std::vector<double> *rw_ndLAr_particle_EndXvtxPip;
  std::vector<double> *rw_ndLAr_particle_EndYvtxPip;
  std::vector<double> *rw_ndLAr_particle_EndZvtxPip;
  std::vector<double> *rw_ndLAr_particle_EndXvtxPim;
  std::vector<double> *rw_ndLAr_particle_EndYvtxPim;
  std::vector<double> *rw_ndLAr_particle_EndZvtxPim;
  std::vector<double> *rw_ndLAr_particle_EndXvtxPi0;
  std::vector<double> *rw_ndLAr_particle_EndYvtxPi0;
  std::vector<double> *rw_ndLAr_particle_EndZvtxPi0;
  std::vector<double> *rw_ndLAr_particle_EndXvtxP;
  std::vector<double> *rw_ndLAr_particle_EndYvtxP;
  std::vector<double> *rw_ndLAr_particle_EndZvtxP;
  std::vector<double> *rw_ndLAr_particle_EndXvtxN;
  std::vector<double> *rw_ndLAr_particle_EndYvtxN;
  std::vector<double> *rw_ndLAr_particle_EndZvtxN;
  std::vector<double> *rw_ndLAr_particle_EndXvtxElectron;
  std::vector<double> *rw_ndLAr_particle_EndYvtxElectron;
  std::vector<double> *rw_ndLAr_particle_EndZvtxElectron;
  std::vector<double> *rw_ndLAr_particle_EndXvtxGamma;
  std::vector<double> *rw_ndLAr_particle_EndYvtxGamma;
  std::vector<double> *rw_ndLAr_particle_EndZvtxGamma;

  // Reconstructed Particle Angles
  std::vector<double> *rw_ndLAr_particle_ThetaRecoMu;
  std::vector<double> *rw_ndLAr_particle_ThetaRecoPip;
  std::vector<double> *rw_ndLAr_particle_ThetaRecoPim;
  std::vector<double> *rw_ndLAr_particle_ThetaRecoPi0;
  std::vector<double> *rw_ndLAr_particle_ThetaRecoP;
  std::vector<double> *rw_ndLAr_particle_ThetaRecoN;
  std::vector<double> *rw_ndLAr_particle_ThetaRecoElectron;
  std::vector<double> *rw_ndLAr_particle_ThetaRecoGamma;

  // True Particle Angles
  std::vector<double> *rw_ndLAr_particle_ThetaMu;
  std::vector<double> *rw_ndLAr_particle_ThetaPip;
  std::vector<double> *rw_ndLAr_particle_ThetaPim;
  std::vector<double> *rw_ndLAr_particle_ThetaPi0;
  std::vector<double> *rw_ndLAr_particle_ThetaP;
  std::vector<double> *rw_ndLAr_particle_ThetaN;
  std::vector<double> *rw_ndLAr_particle_ThetaElectron;
  std::vector<double> *rw_ndLAr_particle_ThetaGamma;

  // Reconstructed Track Length
  std::vector<double> *rw_ndLAr_particle_TrackLengthRecoMu;
  std::vector<double> *rw_ndLAr_particle_TrackLengthRecoPip;
  std::vector<double> *rw_ndLAr_particle_TrackLengthRecoPim;
  std::vector<double> *rw_ndLAr_particle_TrackLengthRecoPi0;
  std::vector<double> *rw_ndLAr_particle_TrackLengthRecoP;
  std::vector<double> *rw_ndLAr_particle_TrackLengthRecoN;
  std::vector<double> *rw_ndLAr_particle_TrackLengthRecoElectron;
  std::vector<double> *rw_ndLAr_particle_TrackLengthRecoGamma;
  
  // True Track Length
  std::vector<double> *rw_ndLAr_particle_TrackLengthMu;
  std::vector<double> *rw_ndLAr_particle_TrackLengthPip;
  std::vector<double> *rw_ndLAr_particle_TrackLengthPim;
  std::vector<double> *rw_ndLAr_particle_TrackLengthPi0;
  std::vector<double> *rw_ndLAr_particle_TrackLengthP;
  std::vector<double> *rw_ndLAr_particle_TrackLengthN;
  std::vector<double> *rw_ndLAr_particle_TrackLengthElectron;
  std::vector<double> *rw_ndLAr_particle_TrackLengthGamma;

  // True Contained
  std::vector<int> *rw_ndLAr_particle_ContainedMu;
  std::vector<int> *rw_ndLAr_particle_ContainedPip;
  std::vector<int> *rw_ndLAr_particle_ContainedPim;
  std::vector<int> *rw_ndLAr_particle_ContainedPi0;
  std::vector<int> *rw_ndLAr_particle_ContainedP;
  std::vector<int> *rw_ndLAr_particle_ContainedN;
  std::vector<int> *rw_ndLAr_particle_ContainedElectron;
  std::vector<int> *rw_ndLAr_particle_ContainedGamma;

};

// ********************************
// ND Detector Systematic Functions
// ********************************

// -------------------------------------------------------------------------
// Global ND Energy Scale Systematics - Essentially Calibration Uncertainty
// Don't shift muon energy since that isn't calculated calorimetrically
// -------------------------------------------------------------------------


// Total Energy Scale
inline void TotalEScaleND(const double * par, double * erec, double erecHad, double erecLep, bool NotCCnumu) {

  (*erec) += (*par) * erecHad;

  //if not true CC numu event AND reco nue event
  if (NotCCnumu)
  {
    (*erec) += (*par) * erecLep;
  }

}

// Total Energy Scale Sqrt
inline void TotalEScaleSqrtND(const double * par, double * erec, double erecHad, double erecLep, double sqrtErecHad, double sqrtErecLep, bool NotCCnumu) {

  (*erec) += (*par) * sqrtErecHad * erecHad ;

  //if not true CC numu AND reco nue event
  if (NotCCnumu)
  {
    (*erec) += (*par) * sqrtErecLep * erecLep;
  }

}

// Total Energy Scale Inverse Sqrt
inline void TotalEScaleInvSqrtND(const double * par, double * erec, double erecHad, double erecLep, double invSqrtErecHad, double invSqrtErecLep, bool NotCCnumu) {

  (*erec) += (*par) * invSqrtErecHad * erecHad ;

  //if not true CC numu AND reco nue event
  if (NotCCnumu)
  {
    (*erec) += (*par) * invSqrtErecLep * erecLep;
  }

}

// ---------------------------------------------------------------
// Particle-Specific Uncertainties - Essentially Particle Response
// ---------------------------------------------------------------


// ---------------------------------------------------------------
// CHARGED HADRONS
// ---------------------------------------------------------------


// Charged Hadron Energy Scale 
inline void HadEScaleND(const double * par, double * erec, double sumEhad) {

  // Protons + Positive Pions + Negative Pions
  (*erec) += (*par) * sumEhad;

}

// Charged Hadron Energy Scale Sqrt
inline void HadEScaleSqrtND(const double * par, double * erec, double sumEhad, double sqrtSumEhad) {

  // Protons + Positive Pions + Negative Pions
  (*erec) += (*par) * sqrtSumEhad * sumEhad;

}

// Charged Hadron Energy Scale Inv Sqrt
inline void HadEScaleInvSqrtND(const double * par, double * erec, double sumEhad, double invSqrtSumEhad) {

  // Protons + Positive Pions + Negative Pions
  (*erec) += (*par) * invSqrtSumEhad * sumEhad;

}


// ---------------------------------------------------------------
// Muons
// ---------------------------------------------------------------


// Muon Energy Scale
inline void MuEScaleND(const double * par, double * erec, double erecLep, bool CCnumu) {

  //if true CC numu AND reco numu event
  if (CCnumu)
  {
    (*erec) += (*par) * erecLep;
  }

}

// Muon Energy Scale Sqrt
inline void MuEScaleSqrtND(const double * par, double * erec, double erecLep, double sqrtErecLep, bool CCnumu) {

  //if true CC numu AND reco numu event
  if (CCnumu)
  {
    (*erec) += (*par) * sqrtErecLep * erecLep;
  }

}

// Muon Energy Scale Inverse Sqrt
inline void MuEScaleInvSqrtND(const double * par, double * erec, double erecLep, double invSqrtErecLep, bool CCnumu) {

  //if true CC numu AND reco numu event
  if (CCnumu)
  {
    (*erec) += (*par) * invSqrtErecLep * erecLep;
  }

}

// ---------------------------------------------------------------
// Neutrons
// ---------------------------------------------------------------


// Neutron Energy Scale
inline void NEScaleND(const double * par, double * erec, double eRecoN) {

  (*erec) += (*par) * eRecoN;

}

// Neutron Energy Scale Sqrt
inline void NEScaleSqrtND(const double * par, double * erec, double eRecoN, double sqrteRecoN) {

  (*erec) += (*par) * sqrteRecoN * eRecoN;

}

// Neutron Energy Scale Inverse Sqrt
inline void NEScaleInvSqrtND(const double * par, double * erec, double eRecoN, double invSqrteRecoN) {

  (*erec) += (*par) * invSqrteRecoN * eRecoN;

}

// ---------------------------------------------------------------
// Electromagnetic Shower
// ---------------------------------------------------------------


// Electromagnetic Shower Energy Scale
inline void EMEScaleND(const double * par, double * erec, double eRecoPi0, double erecLep, bool CCnue) {

  (*erec) += (*par) * eRecoPi0;

  //if true CC nue AND reco nue event
  if (CCnue)
  {
    (*erec) += (*par) * erecLep;
  }

}

// Electromagnetic Shower Energy Scale Sqrt
inline void EMEScaleSqrtND(const double * par, double * erec, double eRecoPi0, double erecLep, double sqrtErecLep, double sqrteRecoPi0, bool CCnue) {

  (*erec) += (*par) * sqrteRecoPi0 * eRecoPi0;

  //if true CC nue AND reco nue event
  if (CCnue)
  {
    (*erec) += (*par) * sqrtErecLep * erecLep;
  }

}

// Electromagnetic Shower Energy Scale Inverse Sqrt
inline void EMEScaleInvSqrtND(const double * par, double * erec, double eRecoPi0, double erecLep, double invSqrtErecLep, double invSqrteRecoPi0, bool CCnue) {

  (*erec) += (*par) * invSqrteRecoPi0 * eRecoPi0;

  //if true CC nue AND reco nue event
  if (CCnue)
  {
    (*erec) += (*par) * invSqrtErecLep * erecLep;
  }

}

// ---------------------------------------------------------------
// Resolution Uncertainties
// ---------------------------------------------------------------

// ---------------------------------------------------------------
// CHARGED HADRONS
// ---------------------------------------------------------------
inline void HadResND(const double * par, double * erec, double eRecoP, double eRecoPip, double eRecoPim, double eP, double ePip, double ePim) {

  // Reco Sum: Protons + Positive Pions + Negative Pions
  double recoSum = eRecoP + eRecoPip + eRecoPim;

  // True Sum: Protons + Positive Pions + Negative Pions
  double trueSum = eP + ePip + ePim;

  (*erec) += (*par) * (trueSum - recoSum);

}

// ---------------------------------------------------------------
// Muons
// ---------------------------------------------------------------
inline void MuResND(const double * par, double * erec, double erecLep, double LepE, bool CCnumu) {

  //if true CC numu AND reco numu event
  if (CCnumu)
  {
    (*erec) += (*par) * (LepE - erecLep);
  }

}


// ---------------------------------------------------------------
// Neutron
// ---------------------------------------------------------------
inline void NResND(const double * par, double * erec, double eRecoN, double eN) {

  (*erec) += (*par) * (eN - eRecoN);

}

// ---------------------------------------------------------------
// Electromagnetic Shower
// ---------------------------------------------------------------

inline void EMResND(const double * par, double * erec, double eRecoPi0, double ePi0, double erecLep, double LepE, bool CCnue) {

  (*erec) += (*par) * (ePi0 - eRecoPi0);

  //if true CC nue AND reco nue event
  if (CCnue)
  {
    (*erec) += (*par) * (LepE - erecLep);
  }

}

// ********************************
// FD Detector Systematic Functions
// ********************************

// -------------------------------------------------------------------------
// Global FD Energy Scale Systematics - Essentially Calibration Uncertainty
// Don't shift muon energy since that isn't calculated calorimetrically
// -------------------------------------------------------------------------


// Total Energy Scale
inline void TotalEScaleFD(const double * par, double * erec, double erecHad, double erecLep, bool NotCCnumu) {

  (*erec) += (*par) * erecHad;

  //if not true CC numu event AND reco nue event
  if (NotCCnumu)
  {
    (*erec) += (*par) * erecLep;
  }

}

// Total Energy Scale Sqrt
inline void TotalEScaleSqrtFD(const double * par, double * erec, double erecHad, double erecLep, double sqrtErecHad, double sqrtErecLep, bool NotCCnumu) {

  (*erec) += (*par) * sqrtErecHad * erecHad ;

  //if not true CC numu AND reco nue event
  if (NotCCnumu)
  {
    (*erec) += (*par) * sqrtErecLep * erecLep;
  }

}

// Total Energy Scale Inverse Sqrt
inline void TotalEScaleInvSqrtFD(const double * par, double * erec, double erecHad, double erecLep, double invSqrtErecHad, double invSqrtErecLep, bool NotCCnumu) {

  (*erec) += (*par) * invSqrtErecHad * erecHad ;

  //if not true CC numu AND reco nue event
  if (NotCCnumu)
  {
    (*erec) += (*par) * invSqrtErecLep * erecLep;
  }

}

// ---------------------------------------------------------------
// Particle-Specific Uncertainties - Essentially Particle Response
// ---------------------------------------------------------------


// ---------------------------------------------------------------
// CHARGED HADRONS
// ---------------------------------------------------------------


// Charged Hadron Energy Scale 
inline void HadEScaleFD(const double * par, double * erec, double sumEhad) {

  // Protons + Positive Pions + Negative Pions
  (*erec) += (*par) * sumEhad;

}

// Charged Hadron Energy Scale Sqrt
inline void HadEScaleSqrtFD(const double * par, double * erec, double sumEhad, double sqrtSumEhad) {

  // Protons + Positive Pions + Negative Pions
  (*erec) += (*par) * sqrtSumEhad * sumEhad;

}

// Charged Hadron Energy Scale Inv Sqrt
inline void HadEScaleInvSqrtFD(const double * par, double * erec, double sumEhad, double invSqrtSumEhad) {

  // Protons + Positive Pions + Negative Pions
  (*erec) += (*par) * invSqrtSumEhad * sumEhad;

}


// ---------------------------------------------------------------
// Muons
// ---------------------------------------------------------------


// Muon Energy Scale
inline void MuEScaleFD(const double * par, double * erec, double erecLep, bool CCnumu) {

  //if true CC numu AND reco numu event
  if (CCnumu)
  {
    (*erec) += (*par) * erecLep;
  }

}

// Muon Energy Scale Sqrt
inline void MuEScaleSqrtFD(const double * par, double * erec, double erecLep, double sqrtErecLep, bool CCnumu) {

  //if true CC numu AND reco numu event
  if (CCnumu)
  {
    (*erec) += (*par) * sqrtErecLep * erecLep;
  }

}

// Muon Energy Scale Inverse Sqrt
inline void MuEScaleInvSqrtFD(const double * par, double * erec, double erecLep, double invSqrtErecLep, bool CCnumu) {

  //if true CC numu AND reco numu event
  if (CCnumu)
  {
    (*erec) += (*par) * invSqrtErecLep * erecLep;
  }

}

// ---------------------------------------------------------------
// Neutrons
// ---------------------------------------------------------------


// Neutron Energy Scale
inline void NEScaleFD(const double * par, double * erec, double eRecoN) {

  (*erec) += (*par) * eRecoN;

}

// Neutron Energy Scale Sqrt
inline void NEScaleSqrtFD(const double * par, double * erec, double eRecoN, double sqrteRecoN) {

  (*erec) += (*par) * sqrteRecoN * eRecoN;

}

// Neutron Energy Scale Inverse Sqrt
inline void NEScaleInvSqrtFD(const double * par, double * erec, double eRecoN, double invSqrteRecoN) {

  (*erec) += (*par) * invSqrteRecoN * eRecoN;

}

// ---------------------------------------------------------------
// Electromagnetic Shower
// ---------------------------------------------------------------


// Electromagnetic Shower Energy Scale
inline void EMEScaleFD(const double * par, double * erec, double eRecoPi0, double erecLep, bool CCnue) {

  (*erec) += (*par) * eRecoPi0;

  //if true CC nue AND reco nue event
  if (CCnue)
  {
    (*erec) += (*par) * erecLep;
  }

}

// Electromagnetic Shower Energy Scale Sqrt
inline void EMEScaleSqrtFD(const double * par, double * erec, double eRecoPi0, double erecLep, double sqrtErecLep, double sqrteRecoPi0, bool CCnue) {

  (*erec) += (*par) * sqrteRecoPi0 * eRecoPi0;

  //if true CC nue AND reco nue event
  if (CCnue)
  {
    (*erec) += (*par) * sqrtErecLep * erecLep;
  }

}

// Electromagnetic Shower Energy Scale Inverse Sqrt
inline void EMEScaleInvSqrtFD(const double * par, double * erec, double eRecoPi0, double erecLep, double invSqrtErecLep, double invSqrteRecoPi0, bool CCnue) {

  (*erec) += (*par) * invSqrteRecoPi0 * eRecoPi0;

  //if true CC nue AND reco nue event
  if (CCnue)
  {
    (*erec) += (*par) * invSqrtErecLep * erecLep;
  }

}

// ---------------------------------------------------------------
// Resolution Uncertainties
// ---------------------------------------------------------------

// ---------------------------------------------------------------
// CHARGED HADRONS
// ---------------------------------------------------------------
inline void HadResFD(const double * par, double * erec, double eRecoP, double eRecoPip, double eRecoPim, double eP, double ePip, double ePim) {

  // Reco Sum: Protons + Positive Pions + Negative Pions
  double recoSum = eRecoP + eRecoPip + eRecoPim;

  // True Sum: Protons + Positive Pions + Negative Pions
  double trueSum = eP + ePip + ePim;

  (*erec) += (*par) * (trueSum - recoSum);

}

// ---------------------------------------------------------------
// Muons
// ---------------------------------------------------------------
inline void MuResFD(const double * par, double * erec, double erecLep, double LepE, bool CCnumu) {

  //if true CC numu AND reco numu event
  if (CCnumu)
  {
    (*erec) += (*par) * (LepE - erecLep);
  }

}


// ---------------------------------------------------------------
// Neutron
// ---------------------------------------------------------------
inline void NResFD(const double * par, double * erec, double eRecoN, double eN) {

  (*erec) += (*par) * (eN - eRecoN);

}

// ---------------------------------------------------------------
// Electromagnetic Shower
// ---------------------------------------------------------------

inline void EMResFD(const double * par, double * erec, double eRecoPi0, double ePi0, double erecLep, double LepE, bool CCnue) {

  (*erec) += (*par) * (ePi0 - eRecoPi0);

  //if true CC nue AND reco nue event
  if (CCnue)
  {
    (*erec) += (*par) * (LepE - erecLep);
  }

}

// ---------------------------------------------------------------
// FD Reconstruction Uncertainties - Shift on CVN Scores
// ---------------------------------------------------------------

//CVN Numu
inline void CVNNumuFD(const double * par, double * cvnnumu) {

  (*cvnnumu) += (*par);

}

//CVN Nue
inline void CVNNueFD(const double * par, double * cvnnue) {

  (*cvnnue) += (*par);

}

#endif
