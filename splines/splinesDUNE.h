#ifndef _splinesDUNE_h_
#define _splinesDUNE_h_

#include "splines/splineFDBase.h"

class splinesDUNE : virtual public splineFDBase
{
  public:
	splinesDUNE(covarianceXsec* xsec_cov);
	virtual ~splinesDUNE();
	void SetupSplines(int BinningOpt);
	virtual void FillSampleArray(std::string SampleName, std::vector<std::string> OscChanFileNames) override;
	virtual std::vector< std::vector<int> > StripDuplicatedModes(std::vector< std::vector<int> > InputVector) override;
	virtual std::vector< std::vector<int> > GetEventSplines(std::string SampleName, int iOscChan, int EventMode, double Var1Val, double Var2Val, double Var3Val) override;
};

#endif
