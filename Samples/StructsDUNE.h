#ifndef _StructsDUNE_h_
#define _StructsDUNE_h_

#include "TMatrixD.h"
#include <stdexcept>
#include <cstddef>

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

/// @brief Fixed-capacity inline vector used for the few per-event systematic factors of the
///        atmospheric MC. Same interface as the std::vector it replaces (size, resize, [], begin/end), but no heap
///        block per event (3 per event before: ~2 GB and 45M allocations for 15M events).
template <typename T, size_t N>
struct SmallVecDUNE {
  T fData[N] = {};
  unsigned char fSize = 0;
  size_t size() const { return fSize; }
  bool empty() const { return fSize == 0; }
  void resize(const size_t n) {
    if (n > N) throw std::length_error("SmallVecDUNE: more per-event systematic factors than its capacity (see StructsDUNE.h)");
    for (size_t i = fSize; i < n; ++i) fData[i] = T();
    fSize = static_cast<unsigned char>(n);
  }
  void push_back(const T& v) { resize(fSize + 1u); fData[fSize - 1] = v; }
  void clear() { fSize = 0; }
  T& operator[](const size_t i) { return fData[i]; }
  const T& operator[](const size_t i) const { return fData[i]; }
  T* begin() { return fData; }
  T* end() { return fData + fSize; }
  const T* begin() const { return fData; }
  const T* end() const { return fData + fSize; }
};

struct dunemc_base { // Store variables used in fitting

  double pot_s;
  double norm_s;
  
  double Target; //Target the interaction was on
  
  int neutrinoflavour_interacting;
  int neutrinoflavour_generated;
  int IsCC;
  double OscChannelIndex;
  
  double neutrinoenergy_reco;
  double neutrinoenergy_true;

  double weights_flux;
  double mode;

  // Analysis sample
  unsigned int SampleIndex = M3::_BAD_INT_;
};

struct dunemc_atm : public dunemc_base { // Store variables used by SampleHandlerAtm
  double coszenith_reco;  
  double coszenith_true;
  uint eventID; //Event ID in the original CAF file (for per-event spline loading)
  double MinDistToWall;

  double hadronicenergy_reco;
  double leptonicenergy_reco;
  double inelasticity_reco;

  double TotalDetectorSystematicWeight;
  SmallVecDUNE<double, 4> DetectorSystematicRatios;  ///< one per detector systematic (inline, max 4)

  double TotalAtmFluxSystematicWeight;
  SmallVecDUNE<double, 2> AtmFluxSystFactors;         ///< inline, max 2
  SmallVecDUNE<const M3::float_t*, 2> AtmFluxSystParVals;
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
  double Target;
 
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
