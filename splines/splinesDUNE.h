#ifndef _splinesDUNE_h_
#define _splinesDUNE_h_

#include "splines/splineFDBase.h"

class splinesDUNE : virtual public splineFDBase {
  public:
  splinesDUNE(covarianceXsec* xsec_cov, MaCh3Modes *Modes_);
  virtual ~splinesDUNE();
  
  std::vector<std::string> GetTokensFromSplineName(std::string FullSplineName);
  void FillSampleArray(std::string SampleName, std::vector<std::string> OscChanFileNames) override;
};

#endif
