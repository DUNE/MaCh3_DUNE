#ifndef _StructsDUNE_h_
#define _StructsDUNE_h_

struct dunemc_base {

  int nEvents; // how many MC events are there
  int *Target; //Target the interaction was on

  int *nupdg;
  int *nupdgUnosc;
  
  double *rw_erec;
  double *rw_erec_shifted;
  double *rw_erec_had;
  double *rw_erec_lep;
  double *rw_yrec;
  double *rw_eRecoP;
  double *rw_eRecoPip;
  double *rw_eRecoPim;
  double *rw_eRecoPi0;
  double *rw_eRecoN;
  double *rw_eRecoMuon;

  std::vector<double> *rw_particle_PMom;
  std::vector<double> *rw_particle_PipMom;
  std::vector<double> *rw_particle_PimMom;
  std::vector<double> *rw_particle_Pi0Mom;
  std::vector<double> *rw_particle_NMom;
  std::vector<double> *rw_particle_MuMom;

  std::vector<double> *rw_particle_PTheta;
  std::vector<double> *rw_particle_PipTheta;
  std::vector<double> *rw_particle_PimTheta;
  std::vector<double> *rw_particle_Pi0Theta;
  std::vector<double> *rw_particle_NTheta;
  std::vector<double> *rw_particle_MuTheta;

  std::vector<double> *rw_particle_PStartX;
  std::vector<double> *rw_particle_PStartY;
  std::vector<double> *rw_particle_PStartZ;
  std::vector<double> *rw_particle_PEndX;
  std::vector<double> *rw_particle_PEndY;
  std::vector<double> *rw_particle_PEndZ;

  std::vector<double> *rw_particle_PipStartX;
  std::vector<double> *rw_particle_PipStartY;
  std::vector<double> *rw_particle_PipStartZ;
  std::vector<double> *rw_particle_PipEndX;
  std::vector<double> *rw_particle_PipEndY;
  std::vector<double> *rw_particle_PipEndZ;

  std::vector<double> *rw_particle_PimStartX;
  std::vector<double> *rw_particle_PimStartY;
  std::vector<double> *rw_particle_PimStartZ;
  std::vector<double> *rw_particle_PimEndX;
  std::vector<double> *rw_particle_PimEndY;
  std::vector<double> *rw_particle_PimEndZ;

  std::vector<double> *rw_particle_Pi0StartX;
  std::vector<double> *rw_particle_Pi0StartY;
  std::vector<double> *rw_particle_Pi0StartZ;
  std::vector<double> *rw_particle_Pi0EndX;
  std::vector<double> *rw_particle_Pi0EndY;
  std::vector<double> *rw_particle_Pi0EndZ;

  std::vector<double> *rw_particle_NStartX;
  std::vector<double> *rw_particle_NStartY;
  std::vector<double> *rw_particle_NStartZ;
  std::vector<double> *rw_particle_NEndX;
  std::vector<double> *rw_particle_NEndY;
  std::vector<double> *rw_particle_NEndZ;

  std::vector<double> *rw_particle_MuStartX;
  std::vector<double> *rw_particle_MuStartY;
  std::vector<double> *rw_particle_MuStartZ;
  std::vector<double> *rw_particle_MuEndX;
  std::vector<double> *rw_particle_MuEndY;
  std::vector<double> *rw_particle_MuEndZ;

  double *rw_PMomReco;
  double *rw_PipMomReco;
  double *rw_PimMomReco;
  double *rw_Pi0MomReco;
  double *rw_NMomReco;
  double *rw_MuMomReco;

  double *rw_PThetaReco;
  double *rw_PipThetaReco;
  double *rw_PimThetaReco;
  double *rw_Pi0ThetaReco;
  double *rw_NThetaReco;
  double *rw_MuThetaReco;

  double *rw_E_diff;
  double *rw_E_diff_Muon;
  double *rw_E_diff_Pip;
  double *rw_E_diff_Pim;
  double *rw_E_diff_Pi0;
  double *rw_E_diff_P;
  double *rw_E_diff_N;

  double *rw_LepE;
  std::vector<double> *rw_particle_eP;
  std::vector<double> *rw_particle_ePip;
  std::vector<double> *rw_particle_ePim;
  std::vector<double> *rw_particle_ePi0;
  std::vector<double> *rw_particle_eN;
  std::vector<double> *rw_particle_eMuon;
  double *rw_eP;
  double *rw_ePip;
  double *rw_ePim;
  double *rw_ePi0;
  double *rw_eN;
  double *rw_eMuon;

  double *rw_etru;
  double *rw_etru_vis;
  double *rw_mom;
  double *rw_theta;
  double *rw_Q2;

  double *rw_cvnnumu;
  double *rw_cvnnue;
  double *rw_cvnnumu_shifted;
  double *rw_cvnnue_shifted;
  int *rw_reco_nue;
  int *rw_reco_numu;
  double *rw_berpaacvwgt;
  int    *rw_isCC;
  int    *rw_nuPDGunosc;
  int    *rw_nuPDG;
  int    *rw_run;
  bool    *rw_isFHC;
  double *rw_vtx_x;
  double *rw_vtx_y;
  double *rw_vtx_z;
  double *rw_vtx_end_x;
  double *rw_vtx_end_y;
  double *rw_vtx_end_z;
  double dummy_y;
  double *rw_reco_q;
  double *reco_numu;

  double pot_s;
  double norm_s;

  double *beam_w;
  double *flux_w;

  double *mode;
  int *isbound;

  double *rw_truecz;

  int *nproton; ///< number of (post-FSI) primary protons
  int *nneutron; ///< number of (post-FSI) primary neutrons
  int *npip; ///< number of (post-FSI) primary pi+
  int *npim; ///< number of (post-FSI) primary pi-
  int *npi0; ///< number of (post-FSI) primary pi0

  int *ntruemuon; //number of true muons
  int *ntruemuonprim; //number of true primary muons
  int *nrecomuon; //number of reconstructed muons
  double *nmuonsratio; //number of reco muons divided by number of true muons

  double *rw_lep_pT;  //transverse lepton momentum
  double *rw_lep_pZ; //parallel lepton momentum
  double *rw_reco_vtx_x;
  double *rw_reco_vtx_y;
  double *rw_reco_vtx_z;
  double *rw_reco_vtx_end_x;
  double *rw_reco_vtx_end_y;
  double *rw_reco_vtx_end_z;
  double *rw_reco_px;
  double *rw_reco_py;
  double *rw_reco_pz;
  double *rw_reco_rad;
  double *rw_rad;

  double *rw_px;
  double *rw_py;
  double *rw_pz;

  double *rw_elep_reco;
  double *rw_elep_true;

  int *nrecoparticles;
  bool *in_fdv;

  double *rw_reco_pid;
  int *rw_iscontained;
  // double* *rw_mu_track_length_truth_true;
  // double* *rw_mu_track_length_truth_reco;
  // double* *rw_mu_track_length_reco;
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
