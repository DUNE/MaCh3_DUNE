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

_MaCh3_Safe_Include_Start_ //{
#include "Eigen/Dense"
_MaCh3_Safe_Include_End_ //}

namespace dune::beamoffaxis {

class SampleHandlerBeamOffAxis : public SampleHandlerFD {

  Eigen::VectorXd unweighted_mc;

public:
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

  friend std::vector<
      std::vector<std::vector<std::vector<std::unique_ptr<TH1>>>>>
  GetBinnedWeights(SampleHandlerBeamOffAxis &sample, int iSubSample,
                   std::vector<std::string> ParamNames,
                   std::vector<std::vector<int>> ParamModes,
                   std::vector<double> TrueEBins);

  std::vector<dune::beamoffaxis::EventInfo> DUNEMCEvents;
  std::vector<double> subsample_analysispot;
  std::vector<bool> subsample_is_numode;

  mutable Eigen::MatrixXd cvmx;
  mutable Eigen::MatrixXd icvmx;

  void CleanMemoryBeforeFit() {}

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
    if (!cvmx.size()) {
      // return SampleHandlerFD::GetLikelihood();
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
      icvmx = cvmx.inverse();
    }

    Eigen::VectorXd residual = data - mc;
    double chi2 = residual.transpose() * icvmx * residual;
    return chi2;
  }
};
}
