#pragma once

#include "Parameters/ParameterHandlerRegularised.h"
#include "Samples/BeamOffAxis/EventInfo.h"
#include "Samples/BeamOffAxis/Projections.h"
#include "Samples/BeamOffAxis/Utility.h"
#include "Samples/BinningHandler.h"
#include "Samples/SampleHandlerFD.h"
#include "Splines/BinnedSplineHandlerDUNE.h"

#include "Fitters/FitterBase.h"
#include "Manager/MaCh3Exception.h"
#include "Manager/MaCh3Logger.h"
#include "Manager/Manager.h"
#include "Splines/SplineMonolith.h"
#include <iostream>
#include <string>

_MaCh3_Safe_Include_Start_ //{
#include "Eigen/Dense"
_MaCh3_Safe_Include_End_ //}

namespace dune::beamoffaxis {

class SampleHandlerBeamOffAxis : public SampleHandlerFD {

  Eigen::VectorXd unweighted_mc;

public:
  TH2D* MakeFancyPlot();
  SampleHandlerBeamOffAxis(
      std::string mc_version, ParameterHandlerGeneric *xsec_cov,
      const std::shared_ptr<OscillationHandler> &Oscillator);

  ~SampleHandlerBeamOffAxis() {}

  void BuildRegularisationMatrix(ParameterHandlerRegularised *ParHandler,
                                 double lambda);

  std::vector<double> ReturnKinematicParameterBinning(
      const int iSubSample,
      const std::string &KinematicParameter) const override {
    return SampleHandlerFD::ReturnKinematicParameterBinning(
        iSubSample, KinematicParameter);
  }
  enum class ECovLikelihoodType { kAbsolute, kFractional };

  void SetCovLikelihoodType(ECovLikelihoodType type) { fCovLikelihoodType = type; }

  friend std::vector<
      std::vector<std::vector<std::vector<std::unique_ptr<TH1>>>>>
  GetBinnedWeights(SampleHandlerBeamOffAxis &sample, int iSubSample,
                   std::vector<std::string> ParamNames,
                   std::vector<std::vector<int>> ParamModes,
                   std::vector<double> TrueEBins);

  std::vector<dune::beamoffaxis::EventInfo> DUNEMCEvents;
  std::vector<double> subsample_analysispot;
  std::vector<bool> subsample_is_numode;

  // Selects which covariance-matrix convention GetLikelihood() uses:
  //  - kAbsolute:   cvmx is in bin-content (absolute) units
  //  - kFractional: relative_cvmx holds fractional uncertainties, rescaled
  //                 to absolute units using the current MC prediction
  

  mutable Eigen::MatrixXd cvmx;
  mutable Eigen::MatrixXd relative_cvmx;
  mutable Eigen::MatrixXd icvmx;

  void CleanMemoryBeforeFit() {}

  ECovLikelihoodType GetCovLikelihoodType() const { return fCovLikelihoodType; }
  std::string GetCovLikelihoodTypeString() const;

private:
  //double GetLikelihood_AbsoluteCov() const;
  //double GetLikelihood_FractionalCov() const;

  ECovLikelihoodType fCovLikelihoodType = ECovLikelihoodType::kAbsolute;

protected:
  void Init() override;
  int SetupExperimentMC() override;
  void SetupFDMC() override;
  void AddAdditionalWeightPointers() override;
  void SetupSplines() override;
  void RegisterFunctionalParameters() override;
  void ResetShifts(int iEvent) override;
  void FinaliseShifts(int iEvent) override;

  double ReturnKinematicParameter(dune::beamoffaxis::KinematicTypes KinPar,
                                  int iEvent);
  double ReturnKinematicParameter(int KinematicVariable, int iEvent);
  double ReturnKinematicParameter(std::string KinematicParameter, int iEvent);

  const double *
  GetPointerToKinematicParameter(dune::beamoffaxis::KinematicTypes KinPar,
                                 int iEvent);
  const double *GetPointerToKinematicParameter(std::string KinematicParameter,
                                               int iEvent);
  const double *GetPointerToKinematicParameter(double KinematicVariable,
                                               int iEvent);

  Eigen::VectorXd GetUnweightedMCRate();

  double GetLikelihood() const override {
    switch (fCovLikelihoodType) {
      case ECovLikelihoodType::kAbsolute:
        return GetLikelihood_AbsoluteCov();
      case ECovLikelihoodType::kFractional:
        return GetLikelihood_FractionalCov();
    }
    MACH3LOG_ERROR("Unhandled ECovLikelihoodType");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  double GetLikelihood_AbsoluteCov() const {
    if (!cvmx.size()) {
      cvmx = Eigen::MatrixXd::Zero(SampleHandlerFD_data.size(),
                                    SampleHandlerFD_data.size());
    }

    Eigen::Map<Eigen::VectorXd const> data(SampleHandlerFD_data.data(),
                                            SampleHandlerFD_data.size());
    Eigen::Map<Eigen::VectorXd const> mc(SampleHandlerFD_array.data(),
                                          SampleHandlerFD_array.size());

    if (cvmx.rows() != Binning->GetNBins()) {
      MACH3LOG_ERROR("Covariance matrix ({}x{}) is not correct for "
                      "total number of bins: {}",
                      cvmx.rows(), cvmx.cols(), Binning->GetNBins());
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    if (!icvmx.size()) {
      cvmx.diagonal() += unweighted_mc + data;
      // Guard against structurally dead bins (kinematically forbidden,
      // always zero MC and zero data) which would otherwise leave a zero on
      // the diagonal and make cvmx singular.
      for (int i = 0; i < cvmx.rows(); ++i) {
        if (cvmx(i, i) <= 0.0) {
          cvmx(i, i) = 1.0;
        }
      }
      icvmx = cvmx.inverse();
    }

    Eigen::VectorXd residual = data - mc;
    double chi2 = residual.transpose() * icvmx * residual;
    //std::cout << "[GetLikelihood_AbsoluteCov] nbins=" << mc.size()   << " chi2=" << chi2 << " 0.5*chi2=" << 0.5 * chi2 << std::endl;
    return 0.5 * chi2;
  }

  double GetLikelihood_FractionalCov() const {
    if (!relative_cvmx.size()) {
      relative_cvmx = Eigen::MatrixXd::Zero(SampleHandlerFD_data.size(),
                                             SampleHandlerFD_data.size());
    }

    Eigen::Map<Eigen::VectorXd const> data(SampleHandlerFD_data.data(),
                                            SampleHandlerFD_data.size());
    Eigen::Map<Eigen::VectorXd const> mc(SampleHandlerFD_array.data(),
                                          SampleHandlerFD_array.size());

    if (relative_cvmx.rows() != Binning->GetNBins() ||
        relative_cvmx.cols() != Binning->GetNBins()) {
      MACH3LOG_ERROR("Covariance matrix ({}x{}) wrong dimensions, "
                      "total number of bins: {}",
                      relative_cvmx.rows(), relative_cvmx.cols(),
                      Binning->GetNBins());
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    const int n = static_cast<int>(mc.size());
    constexpr double mc_floor = 1e-3; // treat mc(i) below this as a dead bin

    Eigen::MatrixXd cvmx_local = relative_cvmx;
    Eigen::VectorXd residual = Eigen::VectorXd::Zero(n);

    for (int i = 0; i < n; ++i) {
      if (mc(i) > mc_floor) {
        cvmx_local(i, i) += (unweighted_mc(i) + data(i)) / (mc(i) * mc(i));
        residual(i) = (data(i) - mc(i)) / mc(i);
      }
    }

    for (int i = 0; i < n; ++i) {
      if (cvmx_local(i, i) <= 0.0 || mc(i) <= mc_floor) {
        cvmx_local(i, i) = 1.0;
        residual(i) = 0.0;
      }
    }

    Eigen::LLT<Eigen::MatrixXd> llt(cvmx_local);
    if (llt.info() != Eigen::Success) {
      MACH3LOG_ERROR("Covariance matrix is not positive-definite at this step.");
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    Eigen::VectorXd solved = llt.solve(residual);
    double chi2 = residual.dot(solved);
   // std::cout << "[GetLikelihood_FractionalCov] nbins=" << n<< " chi2=" << chi2 << " 0.5*chi2=" << 0.5 * chi2 << std::endl;

    return 0.5 * chi2;
  }
};

inline std::string ToString(SampleHandlerBeamOffAxis::ECovLikelihoodType type) {
  switch (type) {
    case SampleHandlerBeamOffAxis::ECovLikelihoodType::kAbsolute:
      return "Absolute";
    case SampleHandlerBeamOffAxis::ECovLikelihoodType::kFractional:
      return "Fractional";
  }
  return "Unknown";
}
inline SampleHandlerBeamOffAxis::ECovLikelihoodType
CovLikelihoodTypeFromString(const std::string &s) {
  if (s == "Absolute") return SampleHandlerBeamOffAxis::ECovLikelihoodType::kAbsolute;
  if (s == "Fractional") return SampleHandlerBeamOffAxis::ECovLikelihoodType::kFractional;
  MACH3LOG_ERROR("Unknown CovLikelihoodType '{}', expected 'Absolute' or 'Fractional'", s);
  throw MaCh3Exception(__FILE__, __LINE__);
}

inline std::string SampleHandlerBeamOffAxis::GetCovLikelihoodTypeString() const {
  return ToString(fCovLikelihoodType);
}



} // namespace dune::beamoffaxis