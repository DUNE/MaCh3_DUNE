#include <iostream>
#include "Fitters/MaCh3Factory.h"
#include "Samples/MaCh3DUNEFactory.h"
#include "Parameters/ParameterHandlerRegularised.h"
#include "Samples/SampleHandlerBeamOffAxis.h"

using dune::beamoffaxis::SampleHandlerBeamOffAxis;

int main(int argc, char * argv[]) {

  auto FitManager = MaCh3ManagerFactory(argc, argv);

  ParameterHandlerGeneric* xsec = nullptr;
  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  // Wire up regularisation
  auto* xsecReg = dynamic_cast<ParameterHandlerRegularised*>(xsec);
  if (!xsecReg) {
    MACH3LOG_ERROR("xsec is not a ParameterHandlerRegularised");
    return 1;
  }

  if (FitManager->raw()["General"]["BeamOffAxisRegularisation"]) {
    double lambda = FitManager->raw()["General"]["BeamOffAxisRegularisation"]["Lambda"].as<double>();
    for (auto* handler : DUNEPdfs) {
      if (auto* boaHandler = dynamic_cast<SampleHandlerBeamOffAxis*>(handler)) {
        boaHandler->BuildRegularisationMatrix(xsecReg, lambda);
      }
    }
  }

   //Start from prior (all 1.0) 
   int nPars = xsec->GetNumParams();
   std::vector<double> propVal(nPars);
   for (int i = 0; i < nPars; ++i) {
     propVal[i] = xsec->GetParInit(i); // all 1.0
   }
    // Find all Xsec Norm params and their indices
std::vector<int> xsecNormIndices;
for (int i = 0; i < nPars; ++i) {
    if (xsec->GetParamType(i) != kNorm) continue;
    if (!xsec->IsParFromGroup(i, "Xsec")) continue;
    xsecNormIndices.push_back(i);
}

 MACH3LOG_INFO("All Xsec Norm parameter names:");
auto normPars = xsec->GetNormParsFromSampleName("OffAxisND");
for (int i = 0; i < (int)xsecNormIndices.size(); ++i) {
    int gi = xsecNormIndices[i];
    auto& np = normPars[i];
    double enu_lo = -1, enubias_lo = -1;
    for (size_t iVar = 0; iVar < np.KinematicVarStr.size(); ++iVar) {
        if (np.KinematicVarStr[iVar] == "TrueNeutrinoEnergy") enu_lo = np.Selection[iVar][0][0];
        if (np.KinematicVarStr[iVar] == "Enubias") enubias_lo = np.Selection[iVar][0][0];
    }
    MACH3LOG_INFO("  [{}] {} : Enu_lo={:.4f} Enubias_lo={:.4f}", i, xsec->GetParName(gi), enu_lo, enubias_lo);
}

int nEnubias = 20; // known from the matrix build log

// Test several Enubias slices and Eν pairs within each
MACH3LOG_INFO("=== Testing adjacent Eν pairs (delta = 0.5) ===");
double delta = 0.5;

for (int iEnubias = 0; iEnubias < 20; ++iEnubias) { // test first 3 Enubias slices
    for (int iEnu = 0; iEnu < 30; ++iEnu) { // test first 3 Eν pairs in each slice

        int idxLo = iEnubias + iEnu * nEnubias;       // param in Eν bin iEnu
        int idxHi = iEnubias + (iEnu + 1) * nEnubias; // param in Eν bin iEnu+1

        if (idxLo >= int(xsecNormIndices.size()) || 
            idxHi >= int(xsecNormIndices.size())) continue;

        int gi = xsecNormIndices[idxLo];
        int gj = xsecNormIndices[idxHi];

        MACH3LOG_INFO("--- Enubias slice {}, Eν pair ({},{}) : {} vs {} ---",
                      iEnubias, iEnu, iEnu+1,
                      xsec->GetParName(gi), xsec->GetParName(gj));

        // Coherent
        std::fill(propVal.begin(), propVal.end(), 1.0);
        propVal[gi] = 1.0 + delta;
        propVal[gj] = 1.0 + delta;
        xsecReg->SetParameters(propVal);
        double pen_coherent = xsecReg->GetPenalty();

        // Anti-correlated
        std::fill(propVal.begin(), propVal.end(), 1.0);
        propVal[gi] = 1.0 + delta;
        propVal[gj] = 1.0 - delta;
        xsecReg->SetParameters(propVal);
        double pen_anti = xsecReg->GetPenalty();

        // Single shift
        std::fill(propVal.begin(), propVal.end(), 1.0);
        propVal[gi] = 1.0 + delta;
        xsecReg->SetParameters(propVal);
        double pen_single = xsecReg->GetPenalty();

        MACH3LOG_INFO("  Single shift:      {:.4f}", pen_single);
        MACH3LOG_INFO("  Coherent shift:    {:.4f} (expect ~= single if neighbours connected)", pen_coherent);
        MACH3LOG_INFO("  Anti-correlated:   {:.4f} (expect ~4x single if neighbours connected)", pen_anti);
        MACH3LOG_INFO("  Anti/Single ratio: {:.2f} (expect 4.0 if connected, 2.0 if not)", pen_anti / pen_single);
    }
}


  // // Get index of a parameter to perturb — first Xsec Norm param
  // int testPar = -1;
  // int testParNeighbour = -1;
  // for (int i = 0; i < nPars; ++i) {
  //   if (xsec->GetParamType(i) != kNorm) continue;
  //   if (!xsec->IsParFromGroup(i, "Xsec")) continue;
  //   if (testPar == -1) { testPar = i; continue; }
  //   if (testParNeighbour == -1) { testParNeighbour = i; break; }
  // }

  // if (testPar == -1 || testParNeighbour == -1) {
  //   MACH3LOG_ERROR("Could not find two Xsec Norm parameters to test");
  //   return 1;
  // }

  // MACH3LOG_INFO("Testing with parameters {} and {}", 
  //               xsec->GetParName(testPar), 
  //               xsec->GetParName(testParNeighbour));

  // // --- Scenario 1: all at prior, penalty should be 0 ---
  // xsecReg->SetParameters(propVal);
  // double pen0 = xsecReg->GetPenalty() - 
  //               dynamic_cast<ParameterHandlerGeneric*>(xsecReg)->ParameterHandlerGeneric::GetPenalty();
  // MACH3LOG_INFO("Penalty at prior (all 1.0): {:.6f}", pen0);

  // // --- Scenario 2: move one bin up by delta ---
  // std::vector<double> shifts = {0.1, 0.2, 0.5, 1.0};
  // for (double delta : shifts) {
  //   std::fill(propVal.begin(), propVal.end(), 1.0);
  //   propVal[testPar] = 1.0 + delta;
  //   xsecReg->SetParameters(propVal);
  //   MACH3LOG_INFO("Shift par {} by +{:.1f}, neighbour unchanged: penalty = {:.6f}",
  //                 xsec->GetParName(testPar), delta,
  //                 xsecReg->GetPenalty());
  // }

  // // --- Scenario 3: move two adjacent bins in opposite directions ---
  // MACH3LOG_INFO("--- Anti-correlated shifts ---");
  // for (double delta : shifts) {
  //   std::fill(propVal.begin(), propVal.end(), 1.0);
  //   propVal[testPar]          = 1.0 + delta;
  //   propVal[testParNeighbour] = 1.0 - delta;
  //   xsecReg->SetParameters(propVal);
  //   MACH3LOG_INFO("Par {} = +{:.1f}, par {} = -{:.1f}: penalty = {:.6f}",
  //                 xsec->GetParName(testPar), delta,
  //                 xsec->GetParName(testParNeighbour), delta,
  //                 xsecReg->GetPenalty());
  // }

  // // --- Scenario 4: move two adjacent bins coherently ---
  // MACH3LOG_INFO("--- Coherent shifts (both same direction) ---");
  // for (double delta : shifts) {
  //   std::fill(propVal.begin(), propVal.end(), 1.0);
  //   propVal[testPar]          = 1.0 + delta;
  //   propVal[testParNeighbour] = 1.0 + delta;
  //   xsecReg->SetParameters(propVal);
  //   MACH3LOG_INFO("Par {} = +{:.1f}, par {} = +{:.1f}: penalty = {:.6f}",
  //                 xsec->GetParName(testPar), delta,
  //                 xsec->GetParName(testParNeighbour), delta,
  //                 xsecReg->GetPenalty());
  // }

 
  return 0;
}