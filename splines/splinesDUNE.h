#ifndef _splinesDUNE_h_
#define _splinesDUNE_h_

#include "splines/splineFDBase.h"

/// @brief Specialisation of FD (binned) spline class
class splinesDUNE : virtual public splineFDBase
{
  public:
	/// @brief Constructor
	/// @param xsec_cov cross-section matrix
	splinesDUNE(covarianceXsec* xsec_cov);

	/// @brief Destructor
	virtual ~splinesDUNE();

	// HW: Doesn't appear to be implemented in DUNE splines, commented for now...
	// void SetupSplines(int BinningOpt);

	/// @brief Fills indexing for each sample and generates a large spline vector
	/// @param SampleName Name of sample
	/// @param OscChanFileNames names of oscillation channels in the sample
	virtual void FillSampleArray(std::string SampleName, std::vector<std::string> OscChanFileNames) override;

	/// @brief Find sand strips dubplicate modes to ensure everything corresponds to MaCh3 modes
	/// @param InputVector Input vector of smodes
	/// @return Processed vector with merged modes
	virtual std::vector< std::vector<int> > StripDuplicatedModes(std::vector< std::vector<int> > InputVector) override;

	/// @brief Getter method for each spline
	/// @param SampleName Name of sample
	/// @param iOscChan Oscillation channel of sample
	/// @param EventMode Mode of event
	/// @param Var1Val Bin 1 value
	/// @param Var2Val Bin 2 value
	/// @param Var3Val Bin 3 value
	/// @return Value of spline for given bin
	virtual std::vector< std::vector<int> > GetEventSplines(std::string SampleName, int iOscChan, int EventMode, double Var1Val, double Var2Val, double Var3Val) override;
};

#endif
