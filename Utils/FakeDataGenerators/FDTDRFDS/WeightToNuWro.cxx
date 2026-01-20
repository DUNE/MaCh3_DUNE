#include "WeightToNuWro.h"

#include "BDTFeature.h"

#include "WeightToNuWro_nue_FD_FHC_BDT.h"
#include "WeightToNuWro_nue_FD_RHC_BDT.h"
#include "WeightToNuWro_nuebar_FD_RHC_BDT.h"
#include "WeightToNuWro_numu_FD_FHC_BDT.h"
#include "WeightToNuWro_numu_FD_RHC_BDT.h"
#include "WeightToNuWro_numu_ND_FHC_BDT.h"
#include "WeightToNuWro_numu_ND_RHC_BDT.h"
#include "WeightToNuWro_numubar_FD_RHC_BDT.h"
#include "WeightToNuWro_numubar_ND_RHC_BDT.h"

#include <cmath>

namespace WeightToNuWro {

float Weight(features const &fts) {
  if (!fts.isCC)
    return 1; // only reweight CC's

  // BDT does not know about taus.
  if (std::abs(fts.nuPDG) == 16)
    return 1;

  // Reweighting only applied up to 8 GeV
  if (fts.Ev > 8.)
    return 1;

  union BDTFeature features[22];

  features[0].fvalue = fts.Ev;
  features[1].fvalue = fts.LepE;
  features[2].fvalue = fts.LepNuAngle;
  features[3].fvalue = fts.Q2;
  features[4].fvalue = fts.W;
  features[5].fvalue = fts.X;
  features[6].fvalue = fts.Y;

  features[7].fvalue = fts.nP;
  features[8].fvalue = fts.nN;
  features[9].fvalue = fts.nipip;
  features[10].fvalue = fts.nipim;
  features[11].fvalue = fts.nipi0;
  features[12].fvalue = fts.niem;

  features[13].fvalue = fts.eP;
  features[14].fvalue = fts.eN;
  features[15].fvalue = fts.ePip;
  features[16].fvalue = fts.ePim;
  features[17].fvalue = fts.ePi0;

  // These last 4 are legacy of an initial attempt and I don't think the BDT
  //   is using them but let's fill them properly to be on the safe side.
  features[18].fvalue = (fts.isFD == 0 ? 1 : 0); // FD files have isFD = 999
  features[19].fvalue = fts.isFHC;
  features[20].fvalue =
      (fts.nuPDG > 0 ? 1 : 0); // Is a neutrino, as opposed to antineutrino
  features[21].fvalue =
      (std::abs(fts.nuPDG) == 14 ? 1
                                 : 0); // Is (anti)numu, as opposed to (anti)nue

  // BDT_norm_weight normalizes weights predicted by the BDT to Genie absolute
  //   normalization
  double BDT_norm_weight = 1.;
  // absolute_NuWro_over_Genie_weight gives the absolute normalization
  //   difference between Genie and NuWro, i.e., the inclusive xsec ratio
  //   integrated over the respective fluxes.
  double absolute_NuWro_over_Genie_weight = 1.;
  // Calibration parameters for BDT output (Platt scaling)
  double plattA = -1.;
  double plattB = 0.;

  double BDT_result = 0;

  if (fts.isFD == 0 && fts.isFHC && std::abs(fts.nuPDG) == 14 &&
      fts.nuPDG > 0) {
    BDT_norm_weight = 1.035;
    absolute_NuWro_over_Genie_weight = 1.013;
    plattA = -1.24585102;
    plattB = -0.01783353;
    BDT_result = WeightToNuWro_numu_ND_FHC(features);
  } else if (fts.isFD == 0 && !fts.isFHC && std::abs(fts.nuPDG) == 14 &&
             fts.nuPDG > 0) {
    BDT_norm_weight = 1.046;
    absolute_NuWro_over_Genie_weight = 0.917;
    plattA = -1.22806734;
    plattB = -0.01548346;
    BDT_result = WeightToNuWro_numu_ND_RHC(features);
  } else if (fts.isFD == 0 && !fts.isFHC && std::abs(fts.nuPDG) == 14 &&
             fts.nuPDG <= 0) {
    BDT_norm_weight = 1.019;
    absolute_NuWro_over_Genie_weight = 0.917;
    plattA = -1.3029002;
    plattB = -0.0314188;
    BDT_result = WeightToNuWro_numubar_ND_RHC(features);
  } else if (fts.isFD != 0 && fts.isFHC && std::abs(fts.nuPDG) == 14 &&
             fts.nuPDG > 0) {
    BDT_norm_weight = 1.042;
    absolute_NuWro_over_Genie_weight = 1.011;
    plattA = -1.21073292;
    plattB = -0.02231727;
    BDT_result = WeightToNuWro_numu_FD_FHC(features);
  } else if (fts.isFD != 0 && !fts.isFHC && std::abs(fts.nuPDG) == 14 &&
             fts.nuPDG > 0) {
    BDT_norm_weight = 1.040;
    absolute_NuWro_over_Genie_weight = 0.981;
    plattA = -1.21420542;
    plattB = -0.02423824;
    BDT_result = WeightToNuWro_numu_FD_RHC(features);
  } else if (fts.isFD != 0 && !fts.isFHC && std::abs(fts.nuPDG) == 14 &&
             fts.nuPDG <= 0) {
    BDT_norm_weight = 1.036;
    absolute_NuWro_over_Genie_weight = 0.904;
    plattA = -1.24604201;
    plattB = -0.03690337;
    BDT_result = WeightToNuWro_numubar_FD_RHC(features);
  } else if (fts.isFD != 0 && fts.isFHC && std::abs(fts.nuPDG) != 14 &&
             fts.nuPDG > 0) {
    BDT_norm_weight = 1.041;
    absolute_NuWro_over_Genie_weight = 1.038;
    plattA = -1.22058287;
    plattB = -0.01992826;
    BDT_result = WeightToNuWro_nue_FD_FHC(features);
  } else if (fts.isFD != 0 && !fts.isFHC && std::abs(fts.nuPDG) != 14 &&
             fts.nuPDG > 0) {
    BDT_norm_weight = 1.056;
    absolute_NuWro_over_Genie_weight = 0.988;
    plattA = -1.17760748;
    plattB = -0.02129378;
    BDT_result = WeightToNuWro_nue_FD_RHC(features);
  } else if (fts.isFD != 0 && !fts.isFHC && std::abs(fts.nuPDG) != 14 &&
             fts.nuPDG <= 0) {
    BDT_norm_weight = 1.029;
    absolute_NuWro_over_Genie_weight = 0.912;
    plattA = -1.27055928;
    plattB = -0.0294727;
    BDT_result = WeightToNuWro_nuebar_FD_RHC(features);
  }

  if (BDT_result) {
    return absolute_NuWro_over_Genie_weight * BDT_norm_weight *
           std::exp(plattA * BDT_result + plattB);
  }
  return 1;
}
}
