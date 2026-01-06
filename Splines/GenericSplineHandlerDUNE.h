#pragma once

#include "Splines/SplineBase.h"
#include "Manager/MaCh3Modes.h"

class GenericSplineHandlerDUNE : public SplineBase {
 public:
  GenericSplineHandlerDUNE(ParameterHandlerGeneric* xsec, MaCh3Modes* Modes_);
  virtual ~GenericSplineHandlerDUNE();
};