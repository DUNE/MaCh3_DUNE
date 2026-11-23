#pragma once

#include "BinnedSplineHandlerDUNE.h"
#include "MonolithSplineHandlerDUNE.h"
#include "FastAtmSplineHandler.h"
#include "Manager/MaCh3Modes.h"
#include "Samples/StructsDUNE.h"

enum SplineType {kUndef, kBinned, kMonolith};

class SplineHandlerFactoryDUNE  {  
public:
  /// @param useFastEngine  use FastAtmSplineHandler for PerEvent splines (false = original MonolithSplineHandlerDUNE)
  SplineHandlerFactoryDUNE(ParameterHandlerGeneric* xsec_params, MaCh3Modes* Modes_, const std::vector<dunemc_atm>& dunemcSamples,
			   const std::vector<std::string>& inputSplinesFile, const std::string& sampleName,
			   bool useFastEngine = false, const FastAtmSplineOptions& fastOptions = FastAtmSplineOptions());
  
  virtual ~SplineHandlerFactoryDUNE();
  std::unique_ptr<SplineBase> GetSplineHandler();
  SplineType GetSplineType() const { return fSplineType; }
  
private:
  std::unique_ptr<SplineBase> fSplineHandler;
  SplineType fSplineType;
  ParameterHandlerGeneric* fXsecParams;
};
