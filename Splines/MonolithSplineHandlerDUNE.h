#pragma once

#include "Splines/SplineMonolith.h"

class MonolithSplineHandlerDUNE : virtual public SMonolith {
 public:
  MonolithSplineHandlerDUNE(ParameterHandlerGeneric* xsec);
  virtual ~MonolithSplineHandlerDUNE();

 private:
  void InitFromFile(std::string &spline_filename);
  MonolithSplineHandlerDUNE(std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>> initParams);
  static std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>> GetInitParamsFromConfig(ParameterHandlerGeneric* xsec);
};