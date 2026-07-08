#pragma once

#include "Parameters/ParameterHandlerGeneric.h"


class ParameterHandlerRegularised : public ParameterHandlerGeneric {

  std::function<double(std::vector<double> const &)> penalty;

public:
  template<typename... Args>
  ParameterHandlerRegularised(Args &&... args) : ParameterHandlerGeneric(args...) {
    penalty = [](std::vector<double> const &){ return 0; };
  }

  double GetLikelihood() override {
    double lhpenalty = penalty(_fPropVal);

    return ParameterHandlerGeneric::GetLikelihood() + lhpenalty;
  }
  void SetPenalty(std::function<double(std::vector<double> const &)> f) {
  penalty = std::move(f);
}
double GetPenalty() {
  return penalty(_fPropVal);
}

};
