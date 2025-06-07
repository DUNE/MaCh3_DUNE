#include "BinnedSplineHandlerDUNE.h"
#include "TROOT.h"

BinnedSplineHandlerDUNE::BinnedSplineHandlerDUNE(ParameterHandlerGeneric* xsec_cov, MaCh3Modes* Modes_) : BinnedSplineHandler(xsec_cov, Modes_) {
}

BinnedSplineHandlerDUNE::~BinnedSplineHandlerDUNE() {
}

std::vector<std::string> BinnedSplineHandlerDUNE::GetTokensFromSplineName(std::string FullSplineName) {
  std::vector<std::string> ReturnVec(TokenOrdering::kNTokens);

  TObjArray *tokens = TString(FullSplineName).Tokenize("_");
  
  ReturnVec[TokenOrdering::kSystToken] = ((TObjString *)(tokens->At(1)))->GetString();
  ReturnVec[TokenOrdering::kModeToken] = ((TObjString *)(tokens->At(2)))->GetString();
  // Skip 3 because it's "sp"
  ReturnVec[TokenOrdering::kVar1BinToken] = ((TObjString *)(tokens->At(4)))->GetString();
  ReturnVec[TokenOrdering::kVar2BinToken] = ((TObjString *)(tokens->At(5)))->GetString();
  ReturnVec[TokenOrdering::kVar3BinToken] = "0";
  
  if (tokens->GetEntries() == 7) {
    ReturnVec[TokenOrdering::kVar3BinToken] = ((TObjString *)(tokens->At(6)))->GetString();
  }

  return ReturnVec;
}
