#include "BinnedSplineHandlerDUNE.h"
#include "TROOT.h"

BinnedSplineHandlerDUNE::BinnedSplineHandlerDUNE(
    ParameterHandlerGeneric *xsec_cov, MaCh3Modes *Modes_)
    : BinnedSplineHandler(xsec_cov, Modes_) {}

BinnedSplineHandlerDUNE::~BinnedSplineHandlerDUNE() {}

inline std::vector<std::string> split_by(std::string s, char d = ',') {
  std::vector<std::string> splits;

  std::string s_left = s;
  size_t d_pos = s_left.find_first_of(d);
  while (d_pos != std::string::npos) {
    splits.push_back(s_left.substr(0, d_pos));
    s_left = s_left.substr(d_pos + 1);
    d_pos = s_left.find_first_of(d);
  }
  splits.push_back(s_left);

  return splits;
}

std::vector<std::string> BinnedSplineHandlerDUNE::GetTokensFromSplineName(
    const std::string &FullSplineName) {
  std::vector<std::string> ReturnVec;

  auto tokens = split_by(FullSplineName, '_');

  ReturnVec.push_back(tokens[1]);
  ReturnVec.push_back(tokens[2]);
  ReturnVec.push_back(tokens[4]);
  ReturnVec.push_back(tokens[5]);

  if (tokens.size() == 7) {
    ReturnVec.push_back(tokens[6]);
  }

  return ReturnVec;
}
