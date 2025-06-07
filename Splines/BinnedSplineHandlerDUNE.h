#ifndef _splinesDUNE_h_
#define _splinesDUNE_h_

#include "Splines/BinnedSplineHandler.h"

class BinnedSplineHandlerDUNE : virtual public BinnedSplineHandler {
  public:
  BinnedSplineHandlerDUNE(ParameterHandlerGeneric* xsec_cov, MaCh3Modes *Modes_);
  virtual ~BinnedSplineHandlerDUNE();
  
  std::vector<std::string> GetTokensFromSplineName(std::string FullSplineName);
};

#endif
