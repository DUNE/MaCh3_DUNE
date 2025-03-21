#ifndef _StructsDUNE_h_
#define _StructsDUNE_h_

#include <iostream>

#include "TVector3.h"
#include "TLorentzVector.h"

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

  double *rw_LepE;
  double *rw_eP;
  double *rw_ePip;
  double *rw_ePim;
  double *rw_ePi0;
  double *rw_eN;

  double *rw_etru;
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
  double    *rw_isFHC;
  double *rw_vtx_x;
  double *rw_vtx_y;
  double *rw_vtx_z;
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
  double *rw_reco_rad;
  double *rw_rad;

  double *rw_elep_reco;
  double *rw_elep_true;

  double *rw_erec_had_sqrt;
  double *rw_erec_lep_sqrt;
  double *rw_eRecoPi0_sqrt;
  double *rw_eRecoN_sqrt;

  int *nrecoparticles;
  bool *in_fdv;
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

// ***************************
// Struct to describe the GENIE interaction modes that will be used during the fit
// Note: these modes can be found in https://github.com/GENIE-MC/Generator/blob/R-2_12_10/src/Interaction/ScatteringType.h
enum MaCh3_Mode {
  // ***************************

  // CCQE
  kMaCh3_CCQE              = 0,
  // 1 Kaon
  kMaCh3_CC_Single_Kaon     = 1,
  // CC DIS
  kMaCh3_CC_DIS            = 2,
  // CC RES
  kMaCh3_CC_RES             = 3,
  // CC COH
  kMaCh3_CC_COH             = 4,
  // CC Diffractive
  kMaCh3_CC_Diffractive     = 5,
  // CC Electron EL
  kMaCh3_CC_Nue_EL          = 6,
  // CC Inverse Muon Decay
  kMaCh3_CC_IMD             = 7,
  // CC Atmospheric Neutrino Gamma
  kMaCh3_CC_AMnuGamma       = 8,
  // CC MEC (aka 2p2h)
  kMaCh3_CC_MEC             = 9,
  // CC Coherent Elastic
  kMaCh3_CC_COHEL           = 10,
  // CC Inverse Beta Decay
  kMaCh3_CC_IBD             = 11,
  // CC GlashowRES 
  kMaCh3_CC_GlashowRES      = 12,
  // CC IMD Annihilation
  kMaCh3_CC_IMDAnnihalation = 13,
  
  // NCQE
  kMaCh3_NCQE              = 14,
  // NC DIS
  kMaCh3_NC_DIS            = 15,
  // NC RES
  kMaCh3_NC_RES             = 16,
  // NC COH
  kMaCh3_NC_COH             = 17,
  // CC Diffractive
  kMaCh3_NC_Diffractive     = 18,
  // NC Electron EL
  kMaCh3_NC_Nue_EL          = 19,
  // NC Inverse Muon Decay
  kMaCh3_NC_IMD             = 20,
  // NC Atmospheric Neutrino Gamma
  kMaCh3_NC_AMnuGamma       = 21,
  // NC MEC (aka 2p2h)
  kMaCh3_NC_MEC             = 22,
  // NC Coherent Elastic
  kMaCh3_NC_COHEL           = 23,
  // NC Inverse Beta Decay
  kMaCh3_NC_IBD             = 24,
  // NC GlashowRES 
  kMaCh3_NC_GlashowRES      = 25,
  // NC IMD Annihilation
  kMaCh3_NC_IMDAnnihalation = 26,

  // Keep a counter of the number of MaCh3 modes we have
  kMaCh3_nModes
};


// ***************************
// Struct to describe the SIMB interaction modes which are found inside the FD CAF files
// Note: this modes can be found in https://nusoft.fnal.gov/larsoft/doxsvn/html/namespacesimb.html
enum SIMB_Mode {
  // ***************************

  // Unknown
  kUnknownInteraction        =   -1,
  // QE
  kQE                        =    0,
  // Resonant
  kRes                       =    1,
  // DIS
  kDIS                       =    2,
  // Coherent
  kCoh                       =    3,
  // Coherent Elastic
  kCohElastic                =    4,
  // Electron Scattering
  kElectronScattering        =    5,
  // Inverse Muon Decay Annihliation
  kIMDAnnihilation           =    6,
  // Inverse Beta Decay
  kInverseBetaDecay          =    7, 
  // Glasgow Resonance
  kGlashowResonance          =    8,
  // Atmospheric Muon Nu Gamma
  kAMNuGamma                 =    9,
  // MEC aka 2p2h
  kMEC                       =   10,
  // Diffractive
  kDiffractive               =   11,
  //
  kEM                        =   12,
  //
  kWeakMix                   =   13,
  // Just keep a counter of the number of modes
  kSIMB_nModes
};

// ***************************
// Struct to describe the GENIE interaction modes which are found inside the ND CAF files
// Note: this modes can be found in https://github.com/GENIE-MC/Generator/blob/R-2_12_10/src/Interaction/ScatteringType.N
enum GENIE_Mode {
  // ***************************

  // Unknown
  gUnknownInteraction        =   -1,
  // QE
  gQE                        =    1,
  // Single Kaon
  gSingleKaon                =    2,
  // DIS
  gDIS                       =    3,
  // RES
  gRes                       =    4,
  // Coherent
  gCoh                       =    5,
  // Diffractive
  gDiffractive                =    6,
  // Nu-e Elastic
  gElectronScattering        =    7, 
  // Inverse Muon Decay
  gIMD                       =    8,
  // Atmospheric Muon Nu Gamma
  gAMNuGamma                 =    9,
  // MEC aka 2p2h
  gMEC                       =   10,
  // Coherent Elastic
  gCohElastic                =   11,
  // Inverse Beta Decay
  gInverseBetaDecay          =   12, 
  // Glasgow Resonance
  gGlashowResonance          =   13,
  // Inverse Muon Decay Annihliation
  gIMDAnnihilation           =   14,
  // Just keep a counter of the number of modes
  gGENIE_nModes
};

// **********************
// Convert an input SIMB mode to the MaCh3 (GENIE) modes
inline int SIMBMode_ToMaCh3Mode(int SIMB_mode, int isCC) {
  // **********************

  int ReturnMode = kMaCh3_nModes;
  
  if (isCC == 1) {
    switch (SIMB_mode) {
      //Unknown
    case kUnknownInteraction:
      std::cerr << "Unknown Interaction mode:" << SIMB_mode << std::endl;
      throw;
      //QE
    case kQE:
      ReturnMode = kMaCh3_CCQE; // CC QE in MaCh3 DUNE
      break;
      // DIS
    case kDIS:
      ReturnMode = kMaCh3_CC_DIS; //CC DIS in MaCh3 DUNE
      break;
      // RES
    case kRes:
      ReturnMode = kMaCh3_CC_RES; // CC RES in MaCh3 DUNE
      break;
      // Coherent
    case kCoh:
      ReturnMode = kMaCh3_CC_COH; // CC Coh in MaCh3 DUNE
      break;
      // Diffractive
    case kDiffractive:
      ReturnMode = kMaCh3_CC_Diffractive; // CC multi-pi in MaCh3
      break;
      // Nue Elastic
    case kElectronScattering:
      ReturnMode = kMaCh3_CC_Nue_EL; // CC Nue scattering in MaCh3 DUNE
      break;
      // Atmospheric Mu Gamma
    case kAMNuGamma:
      ReturnMode = kMaCh3_CC_AMnuGamma; // CC Am Nu Mu in MaCh3 DUNE
      break;
      // MEC
    case kMEC:
      ReturnMode = kMaCh3_CC_MEC; // CC MEC in MaCh3 DUNE
      break;
      // Coherent Elastic
    case kCohElastic:
      ReturnMode = kMaCh3_CC_COHEL; // CC Coherent Elastic in MaCh3 DUNE
      break;
      // Inverse Beta Decay
    case kInverseBetaDecay:
      ReturnMode = kMaCh3_CC_IBD; // CC Inverse Beta Decay in MaCh3 DUNE
      break;
      // Glashow Resonance
    case kGlashowResonance:
      ReturnMode = kMaCh3_CC_GlashowRES; // CC Glashow Reaonance in DUNE
      break;
      // IMD Annihilation
    case kIMDAnnihilation:
      ReturnMode = kMaCh3_CC_IMDAnnihalation; // CC Inverse Muon Decay Annihalation in DUNE
      break;
    default:
      std::cerr << "Invalid mode:" << SIMB_mode << std::endl;
      throw;
    }
    
  }
  
  
  if (isCC == 0) {
    switch (SIMB_mode) {
      //Unknown
    case kUnknownInteraction:
      std::cerr << "Unknown Interaction mode:" << SIMB_mode << std::endl;
      throw;
      //QE
    case kQE:
      ReturnMode = kMaCh3_NCQE; // NC QE in MaCh3 DUNE
      break;
      // DIS
    case kDIS:
      ReturnMode = kMaCh3_NC_DIS; // NC DIS in MaCh3 DUNE
      break;
      // RES
    case kRes:
      ReturnMode = kMaCh3_NC_RES; // NC RES in MaCh3 DUNE
      break;
      // Coherent
    case kCoh:
      ReturnMode = kMaCh3_NC_COH; // NC Coh in MaCh3 DUNE
      break;
      // Diffractive
    case kDiffractive:
      ReturnMode = kMaCh3_NC_Diffractive; // CC multi-pi in MaCh3
      break;
      // Nue Elastic
    case kElectronScattering:
      ReturnMode = kMaCh3_NC_Nue_EL; // NC Nue scattering in MaCh3 DUNE
      break;
      // Atmospheric Mu Gamma
    case kAMNuGamma:
      ReturnMode = kMaCh3_NC_AMnuGamma; // NC Am Nu Mu in MaCh3 DUNE
      break;
      // MEC
    case kMEC:
      ReturnMode = kMaCh3_NC_MEC; // NC MEC in MaCh3 DUNE
      break;
      // Coherent Elastic
    case kCohElastic:
      ReturnMode = kMaCh3_NC_COHEL; // NC Coherent Elastic in MaCh3 DUNE
      break;
      // Inverse Beta Decay
    case kInverseBetaDecay:
      ReturnMode = kMaCh3_NC_IBD; // Inverse Beta Decay in MaCh3 DUNE
      break;
      // Glashow Resonance
    case kGlashowResonance:
      ReturnMode = kMaCh3_NC_GlashowRES; 
      break;
      // IMD Annihilation
    case kIMDAnnihilation:
      ReturnMode = kMaCh3_NC_IMDAnnihalation; 
      break;
    default:
      std::cerr << "Invalid mode:" << SIMB_mode << std::endl;
      throw;
    }
    
  }

  return ReturnMode;
};

// **********************
// Convert an input SIMB mode to the MaCh3 (GENIE) modes
inline int GENIEMode_ToMaCh3Mode(int GENIE_mode, int isCC) {
  // **********************

  int ReturnMode = kMaCh3_nModes;
  
  if (isCC == 1) {
    switch (GENIE_mode) {
        //Unknown
      case gUnknownInteraction:
        ReturnMode = kMaCh3_nModes;
      	    break;
        //QE
      case gQE:
        ReturnMode = kMaCh3_CCQE; // CC QE in MaCh3 DUNE
        break;
        // DIS
      case gDIS:
        ReturnMode = kMaCh3_CC_DIS; //CC DIS in MaCh3 DUNE
        break;
        // RES
      case gRes:
        ReturnMode = kMaCh3_CC_RES; // CC RES in MaCh3 DUNE
        break;
        // Coherent
      case gCoh:
        ReturnMode = kMaCh3_CC_COH; // CC Coh in MaCh3 DUNE
        break;
        // Diffractive
      case gDiffractive:
        ReturnMode = kMaCh3_CC_Diffractive; // CC multi-pi in MaCh3
        break;
        // Nue Elastic
      case gElectronScattering:
        ReturnMode = kMaCh3_CC_Nue_EL; // CC Nue scattering in MaCh3 DUNE
        break;
        // Atmospheric Mu Gamma
      case gAMNuGamma:
        ReturnMode = kMaCh3_CC_AMnuGamma; // CC Am Nu Mu in MaCh3 DUNE
        break;
        // MEC
      case gMEC:
        ReturnMode = kMaCh3_CC_MEC; // CC MEC in MaCh3 DUNE
        break;
        // Coherent Elastic
      case gCohElastic:
        ReturnMode = kMaCh3_CC_COHEL; // CC Coherent Elastic in MaCh3 DUNE
        break;
        // Inverse Beta Decay
      case gInverseBetaDecay:
        ReturnMode = kMaCh3_CC_IBD; // CC Inverse Beta Decay in MaCh3 DUNE
        break;
        // Glashow Resonance
      case gGlashowResonance:
        ReturnMode = kMaCh3_CC_GlashowRES; // NC 1 gamma
        break;
        // IMD Annihilation
      case gIMDAnnihilation:
        ReturnMode = kMaCh3_CC_IMDAnnihalation; // NC other in MaCh3
        break;
      case gIMD:
        ReturnMode = kMaCh3_CC_IMDAnnihalation; // Stick Inverse Muon Decay into above
        break;
      case gSingleKaon:
        ReturnMode = kMaCh3_CC_IMDAnnihalation; // Stick Single Kaon into above
        break;
      default:
        ReturnMode = kMaCh3_nModes; // Something else in MaCh3 (sand?)
        break;
    }

  }


  if (isCC == 0) {
    switch (GENIE_mode) {
        //Unknown
      case kUnknownInteraction:
        ReturnMode = kMaCh3_nModes;
      	    break;
        //QE
      case gQE:
        ReturnMode = kMaCh3_NCQE; // NC QE in MaCh3 DUNE
        break;
        // DIS
      case gDIS:
        ReturnMode = kMaCh3_NC_DIS; // NC DIS in MaCh3 DUNE
        break;
        // RES
      case gRes:
        ReturnMode = kMaCh3_NC_RES; // NC RES in MaCh3 DUNE
        break;
        // Coherent
      case gCoh:
        ReturnMode = kMaCh3_NC_COH; // NC Coh in MaCh3 DUNE
        break;
        // Diffractive
      case gDiffractive:
        ReturnMode = kMaCh3_NC_Diffractive; // CC multi-pi in MaCh3
        break;
        // Nue Elastic
      case gElectronScattering:
        ReturnMode = kMaCh3_NC_Nue_EL; // NC Nue scattering in MaCh3 DUNE
        break;
        // Atmospheric Mu Gamma
      case gAMNuGamma:
        ReturnMode = kMaCh3_NC_AMnuGamma; // NC Am Nu Mu in MaCh3 DUNE
        break;
        // MEC
      case gMEC:
        ReturnMode = kMaCh3_NC_MEC; // NC MEC in MaCh3 DUNE
        break;
        // Coherent Elastic
      case gCohElastic:
        ReturnMode = kMaCh3_NC_COHEL; // NC Coherent Elastic in MaCh3 DUNE
        break;
        // Inverse Beta Decay
      case gInverseBetaDecay:
        ReturnMode = kMaCh3_NC_IBD; // Inverse Beta Decay in MaCh3 DUNE
        break;
        // Glashow Resonance
      case gGlashowResonance:
        ReturnMode = kMaCh3_NC_GlashowRES;
        break;
        // IMD Annihilation
      case gIMDAnnihilation:
        ReturnMode = kMaCh3_NC_IMDAnnihalation;
        break;
      case gIMD:
        ReturnMode = kMaCh3_NC_IMDAnnihalation; 
        break;
      case gSingleKaon:
        ReturnMode = kMaCh3_NC_IMDAnnihalation; 
        break;
      default:
        ReturnMode = kMaCh3_nModes; // Something else in MaCh3 (sand?)
        break;
    }

  }

  return ReturnMode;
};

// *******************************
// Convert MaCh3 mode to SK name
inline std::string MaCh3mode_ToDUNEString(MaCh3_Mode i) {
  // *******************************
  std::string name;

  switch(i) {
  case kMaCh3_CCQE:
    name = "ccqe";
    break;
  case kMaCh3_CC_Single_Kaon:
    name = "unknown";
    break;
  case kMaCh3_CC_DIS:
    name = "ccdis";
    break;
  case kMaCh3_CC_RES:
    name = "ccres";
    break;
  case kMaCh3_CC_COH:
    name = "cccoh";
    break;
  case kMaCh3_CC_COHEL:
    name = "cccohel";
    break;
  case kMaCh3_CC_Diffractive:
    name = "ccdiff";
    break;
  case kMaCh3_CC_Nue_EL:
    name = "ccnueel";
    break;
  case kMaCh3_CC_MEC:
    name = "ccmec";
    break;
  case kMaCh3_CC_IMD:
    name = "ccIMD";
    break;
  case kMaCh3_CC_GlashowRES:
    name = "ccglasres";
    break;
  case kMaCh3_CC_IMDAnnihalation:
    name = "ccimdannihilation";
    break;
  case kMaCh3_CC_AMnuGamma:
    name = "ccamnugamma";
    break;
  case kMaCh3_CC_IBD:
    name = "ccibd";
    break;
  case kMaCh3_NCQE:
    name = "ncqe";
    break;
  case kMaCh3_NC_DIS:
    name = "ncdis";
    break;
  case kMaCh3_NC_RES:
    name = "ncres";
    break;
  case kMaCh3_NC_COH:
    name = "nccoh";
    break;
  case kMaCh3_NC_COHEL:
    name = "nccohel";
    break;
  case kMaCh3_NC_Diffractive:
    name = "ncdiff";
    break;
  case kMaCh3_NC_Nue_EL:
    name = "ncnueel";
    break;
  case kMaCh3_NC_MEC:
    name = "ncmec";
    break;
  case kMaCh3_NC_IMD:
    name = "ncIMD";
    break;
  case kMaCh3_NC_GlashowRES:
    name = "ncglasres";
    break;
  case kMaCh3_NC_IMDAnnihalation:
    name = "ncimdannihilation";
    break;
  case kMaCh3_NC_AMnuGamma:
    name = "ncamnugamma";
    break;
  case kMaCh3_NC_IBD:
    name = "ncibd";
    break;
  case kMaCh3_nModes:
    name = "unknown";
    break;
  default:
    std::cerr << "Did not find the MaCh3 mode you specified " << i << std::endl;
    name = "UNKNOWN_BAD";
  }
  
  return name;
}

enum MaCh3_Spline_Modes {

  // ***************************
  // CCQE
  kMaCh3_Spline_CCQE              = 0,
  // 1 Kaon
  kMaCh3_Spline_CC_Single_Kaon     = 1,
  // CC DIS
  kMaCh3_Spline_CC_DIS            = 2,
  // CC RES
  kMaCh3_Spline_CC_RES             = 3,
  // CC COH
  kMaCh3_Spline_CC_COH             = 4,
  // CC Diffractive
  kMaCh3_Spline_CC_Diffractive     = 5,
  // CC Electron EL
  kMaCh3_Spline_CC_Nue_EL          = 6,
  // CC Inverse Muon Decay
  kMaCh3_Spline_CC_IMD             = 7,
  // CC Atmospheric Neutrino Gamma
  kMaCh3_Spline_CC_AMnuGamma       = 8,
  // CC MEC (aka 2p2h)
  kMaCh3_Spline_CC_MEC             = 9,
  // CC Coherent Elastic
  kMaCh3_Spline_CC_COHEL           = 10,
  // CC Inverse Beta Decay
  kMaCh3_Spline_CC_IBD             = 11,
  // CC GlashowRES 
  kMaCh3_Spline_CC_GlashowRES      = 12,
  // CC IMD Annihilation
  kMaCh3_Spline_CC_IMDAnnihalation = 13,
  
  // NCQE
  kMaCh3_Spline_NCQE              = 14,
  // NC DIS
  kMaCh3_Spline_NC_DIS            = 15,
  // NC RES
  kMaCh3_Spline_NC_RES             = 16,
  // NC COH
  kMaCh3_Spline_NC_COH             = 17,
  // CC Diffractive
  kMaCh3_Spline_NC_Diffractive     = 18,
  // NC Electron EL
  kMaCh3_Spline_NC_Nue_EL          = 19,
  // NC Inverse Muon Decay
  kMaCh3_Spline_NC_IMD             = 20,
  // NC Atmospheric Neutrino Gamma
  kMaCh3_Spline_NC_AMnuGamma       = 21,
  // NC MEC (aka 2p2h)
  kMaCh3_Spline_NC_MEC             = 22,
  // NC Coherent Elastic
  kMaCh3_Spline_NC_COHEL           = 23,
  // NC Inverse Beta Decay
  kMaCh3_Spline_NC_IBD             = 24,
  // NC GlashowRES 
  kMaCh3_Spline_NC_GlashowRES      = 25,
  // NC IMD Annihilation
  kMaCh3_Spline_NC_IMDAnnihalation = 26,

  // Keep a counter of the number of MaCh3 modes we have
  kMaCh3_Spline_nModes

};

inline int MaCh3Mode_to_SplineMode(int iMode){
  //No grouping of modes in MaCh3
  return iMode;
}

#endif
