#pragma once

#include "BinnedSplineHandlerDUNE.h"
#include "MonolithSplineHandlerDUNE.h"
#include "Manager/MaCh3Modes.h"
#include "Samples/StructsDUNE.h"

enum SplineType {kUndef, kBinned, kMonolith};

class SplineHandlerFactoryDUNE  {
  
 public:
  SplineHandlerFactoryDUNE(ParameterHandlerGeneric* xsec_params, MaCh3Modes* Modes_,
    const std::vector<dunemc_base>& dunemcSamples, const std::string& inputSplinesFile,
    const std::string& sampleName);
  virtual ~SplineHandlerFactoryDUNE();
  std::unique_ptr<SplineBase> GetSplineHandler();
  SplineType GetSplineType() const { return fSplineType; }

  private:
    std::unique_ptr<SplineBase> fSplineHandler;
    SplineType fSplineType;
    ParameterHandlerGeneric* fXsecParams;
  };