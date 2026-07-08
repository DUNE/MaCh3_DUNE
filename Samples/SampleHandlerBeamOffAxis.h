#pragma once

#include "Splines/BinnedSplineHandlerDUNE.h"
#include "Samples/SampleHandlerFD.h"
#include "Samples/BeamOffAxis/EventInfo.h"
#include "Samples/BeamOffAxis/Projections.h"
#include "Samples/BeamOffAxis/Utility.h"
#include "Parameters/ParameterHandlerRegularised.h"

#include "Manager/MaCh3Exception.h"
#include "Manager/MaCh3Logger.h"
#include "Splines/SplineMonolith.h"
#include "Fitters/FitterBase.h"
#include "Manager/Manager.h"

_MaCh3_Safe_Include_Start_ //{
#include "Eigen/Dense"
_MaCh3_Safe_Include_End_ //}

namespace dune::beamoffaxis {

class SampleHandlerBeamOffAxis : public SampleHandlerFD {
public:
  SampleHandlerBeamOffAxis(
      std::string mc_version, ParameterHandlerGeneric *xsec_cov,
      const std::shared_ptr<OscillationHandler> &Oscillator);

  ~SampleHandlerBeamOffAxis() {}

  void BuildRegularisationMatrix(ParameterHandlerRegularised *ParHandler, double lambda);

  std::vector<double> ReturnKinematicParameterBinning(
      const int iSubSample,
      const std::string &KinematicParameter) const override {
    return SampleHandlerFD::ReturnKinematicParameterBinning(
        iSubSample, KinematicParameter);
  }

  friend std::vector<std::vector<std::vector<std::vector<std::unique_ptr<TH1>>>>>
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

  double ReturnKinematicParameter(dune::beamoffaxis::KinematicTypes KinPar, int iEvent);
  double ReturnKinematicParameter(int KinematicVariable, int iEvent);
  double ReturnKinematicParameter(std::string KinematicParameter, int iEvent);

  const double *GetPointerToKinematicParameter(dune::beamoffaxis::KinematicTypes KinPar, int iEvent);
  const double *GetPointerToKinematicParameter(std::string KinematicParameter, int iEvent);
  const double *GetPointerToKinematicParameter(double KinematicVariable, int iEvent);

// ////////////////////////////////////////////ORIGINAL!!
// double GetLikelihood() const override {
//     // If no covariance matrix loaded, delegate to base class
//     if (!cvmx.size()) {
//         //std::cout << "[GetLikelihood] No covariance matrix loaded, delegating to base class" << std::endl;
//         return SampleHandlerFD::GetLikelihood();
//     }

//     Eigen::Map<Eigen::VectorXd const> data(SampleHandlerFD_data.data(),
//                                            SampleHandlerFD_data.size());
//     Eigen::Map<Eigen::VectorXd const> mc(SampleHandlerFD_array.data(),
//                                          SampleHandlerFD_array.size());

//     if (cvmx.rows() != static_cast<int> (data.size())) {
//         MACH3LOG_ERROR("Covariance matrix ({}x{}) is not correct for "
//                        "data array size: {}",
//                        cvmx.rows(), cvmx.cols(), data.rows());
//         throw MaCh3Exception(__FILE__, __LINE__);
//     }

    
//     Eigen::MatrixXd cvmx_with_stats = cvmx; // Build a fresh working copy each step so the inversion uses the current MC
//     cvmx_with_stats.diagonal() += data;


//     Eigen::MatrixXd icvmx_current = cvmx_with_stats.inverse();

//     Eigen::VectorXd residual = data - mc;
//     double chi2 = residual.transpose() * icvmx_current * residual;
//     return chi2;
// }

//////////

////////////////////////MCHist NEW Verssion!
Eigen::VectorXd GetUnweightedMCRate() const {
  Eigen::VectorXd unweighted_mc(GetNBins());
  unweighted_mc.setZero();

  for (int sample_i = 0; sample_i < NSamples; ++sample_i) {
    const int ndim = GetNDim(sample_i);

    if (ndim == 1) {
      auto mc_hist = const_cast<SampleHandlerBeamOffAxis*>(this)
                         ->Get1DVarHist(sample_i, GetXBinVarName(sample_i), SelectionVector, 1);
      if (!mc_hist) throw MaCh3Exception(__FILE__, __LINE__);

      for (int i = 0; i < mc_hist->GetNbinsX(); ++i) {
        const int global_bin = GetGlobalBinSafe(sample_i, {i});
        unweighted_mc(global_bin) = mc_hist->GetBinContent(i + 1);
      }
    } else if (ndim == 2) {
      auto mc_hist = const_cast<SampleHandlerBeamOffAxis*>(this)
                         ->Get2DVarHist(sample_i, GetXBinVarName(sample_i), GetYBinVarName(sample_i),
                                        SelectionVector, 1);
      if (!mc_hist) throw MaCh3Exception(__FILE__, __LINE__);

      for (int j = 0; j < mc_hist->GetNbinsY(); ++j) {
        for (int i = 0; i < mc_hist->GetNbinsX(); ++i) {
          const int global_bin = GetGlobalBinSafe(sample_i, {i, j});
          unweighted_mc(global_bin) = mc_hist->GetBinContent(i + 1, j + 1);
        }
      }
    } else {
      MACH3LOG_ERROR("GetUnweightedMCRate: sample {} has {} dimensions, only 1D or 2D supported",
                     sample_i, ndim);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  return unweighted_mc;
}

double GetLikelihood() const override {
  if (!cvmx.size()) {
    return SampleHandlerFD::GetLikelihood();
  }

  Eigen::Map<Eigen::VectorXd const> data(SampleHandlerFD_data.data(),
                                         SampleHandlerFD_data.size());
  Eigen::Map<Eigen::VectorXd const> mc(SampleHandlerFD_array.data(),
                                       SampleHandlerFD_array.size());

  if (cvmx.rows() != GetNBins()) {
    MACH3LOG_ERROR("Covariance matrix ({}x{}) is not correct for "
                   "total number of bins: {}",
                   cvmx.rows(), cvmx.cols(), GetNBins());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  Eigen::VectorXd unweighted_mc = GetUnweightedMCRate();

  Eigen::MatrixXd cvmx_with_stats = cvmx; // fresh working copy each step
  cvmx_with_stats.diagonal() += unweighted_mc;

  Eigen::MatrixXd icvmx_current = cvmx_with_stats.inverse();

  Eigen::VectorXd residual = data - mc;
  double chi2 = residual.transpose() * icvmx_current * residual;
  return chi2;
}
}

};
