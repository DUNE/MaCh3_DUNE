#include "BinnedSplineHandlerDUNE.h"
#include "TROOT.h"

BinnedSplineHandlerDUNE::BinnedSplineHandlerDUNE(ParameterHandlerGeneric* xsec_cov, MaCh3Modes* Modes_) : BinnedSplineHandler(xsec_cov, Modes_) {
}

BinnedSplineHandlerDUNE::~BinnedSplineHandlerDUNE() {
}

std::vector<std::string> BinnedSplineHandlerDUNE::GetTokensFromSplineName(const std::string& FullSplineName) {
  std::vector<std::string> ReturnVec(5);

  TObjArray *tokens = TString(FullSplineName).Tokenize("_");
  
  ReturnVec[TokenOrdering::kSystToken] = ((TObjString *)(tokens->At(1)))->GetString();
  ReturnVec[TokenOrdering::kModeToken] = ((TObjString *)(tokens->At(2)))->GetString();
  // Skip 3 because it's "sp"
  ReturnVec[TokenOrdering::kVarBinToken + 0] = ((TObjString *)(tokens->At(4)))->GetString();
  ReturnVec[TokenOrdering::kVarBinToken + 1] = ((TObjString *)(tokens->At(5)))->GetString();
  ReturnVec[TokenOrdering::kVarBinToken + 2] = "0";
  
  if (tokens->GetEntries() == 7) {
    ReturnVec[TokenOrdering::kVarBinToken + 2] = ((TObjString *)(tokens->At(6)))->GetString();
  }

  return ReturnVec;
}
