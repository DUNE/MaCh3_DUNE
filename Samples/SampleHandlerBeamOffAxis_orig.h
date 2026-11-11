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
  mutable Eigen::MatrixXd relative_cvmx;
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

  // ////////////////////////////////////////////////Original LLH!
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
    
    // for (size_t i = 0; i < SampleHandlerFD_data.size() ; ++i){
    //   std::cout << " DATA"  << SampleHandlerFD_data[i] << " " << ""  <<  std::endl;
    //   std::cout << " Pred " << SampleHandlerFD_array[i] << " " << "" <<  std::endl;
    // }

    if (!icvmx.size()) {
      cvmx.diagonal() += unweighted_mc + data;
      // Guard against structurally dead bins (kinematically forbidden, always
      // zero MC and zero data) which would otherwise leave a zero on the
      // diagonal and make cvmx singular.
      for (int i = 0; i < cvmx.rows(); ++i) {
        if (cvmx(i, i) <= 0.0) {
          cvmx(i, i) = 1.0;
        }
      }
      icvmx = cvmx.inverse();
    }

    
        
    Eigen::VectorXd residual = data - mc;
    double chi2 = residual.transpose() * icvmx * residual; //added 0.5?
    
    // std::cout << "Bin | Data | MC | CovDiag | BinChi2" << std::endl;
    // double sum_diag_chi2 = 0.0;
    // for (int i = 0; i < data.size(); ++i) {
    //   double bin_residual = data(i) - mc(i);
    //   // Diagonal-only approximation: ignores off-diagonal covariance terms
    //   double bin_chi2 = bin_residual * bin_residual * icvmx(i, i);
    //   sum_diag_chi2 += bin_chi2;

    //   std::cout << i << " | " << data(i) << " | " << mc(i)
    //             << " | " << cvmx(i, i)
    //             << " | " << bin_chi2 << std::endl;
    // }
    // std::cout << "Total chi2 (full matrix) = " << chi2 << std::endl;
    // std::cout << "Sum of diagonal-only bin chi2 = " << sum_diag_chi2 << std::endl;

    return 0.5* chi2;
  }



// ///////////New version
//double GetLikelihood() const override {
//   if (!relative_cvmx.size()) {
//     relative_cvmx = Eigen::MatrixXd::Zero(SampleHandlerFD_data.size(),
//                                            SampleHandlerFD_data.size());
//   }

//   Eigen::Map<Eigen::VectorXd const> data(SampleHandlerFD_data.data(),
//                                           SampleHandlerFD_data.size());
//   Eigen::Map<Eigen::VectorXd const> mc(SampleHandlerFD_array.data(),
//                                         SampleHandlerFD_array.size());

//   if (relative_cvmx.rows() != Binning->GetNBins()) {
//     MACH3LOG_ERROR("Covariance matrix ({}x{}) is not correct for "
//                     "total number of bins: {}",
//                     relative_cvmx.rows(), relative_cvmx.cols(), Binning->GetNBins());
//     throw MaCh3Exception(__FILE__, __LINE__);
//   }

//   // Rescale fractional cov mat to absolute using the current prediction: 
//   // Cov_abs(i,j) = Cov_rel(i,j) * mc(i) * mc(j)
//   cvmx = relative_cvmx.array() * (mc * mc.transpose()).array();

//   // Statistical term (BB). unweighted_mc is fixed at(nominal), which is exact here since all systematics in
//   cvmx.diagonal() += unweighted_mc + data;

//   // Guard against structurally dead bins (with zero MC and zero data) which would otherwise leave a zero on the
//   // diagonal and make cvmx singular.
//   for (int i = 0; i < cvmx.rows(); ++i) {
//     if (cvmx(i, i) <= 0.0) {
//       cvmx(i, i) = 1.0;
//     }
//   }

//   //Cholesky decom
//   Eigen::LLT<Eigen::MatrixXd> llt(cvmx);
//   if (llt.info() != Eigen::Success) {
//     MACH3LOG_ERROR("Covariance matrix is not positive-definite at this step.");
//     throw MaCh3Exception(__FILE__, __LINE__);
//   }

//   Eigen::VectorXd residual = data - mc;
//   double chi2 = residual.dot(llt.solve(residual));

//   return chi2;
// }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////New 2.0 RELATIVE!!
//   double GetLikelihood() const override {
//   if (!relative_cvmx.size()) {
//     relative_cvmx = Eigen::MatrixXd::Zero(SampleHandlerFD_data.size(),
//                                            SampleHandlerFD_data.size());
//   }

//   Eigen::Map<Eigen::VectorXd const> data(SampleHandlerFD_data.data(),
//                                           SampleHandlerFD_data.size());
//   Eigen::Map<Eigen::VectorXd const> mc(SampleHandlerFD_array.data(),
//                                         SampleHandlerFD_array.size());

//   if (relative_cvmx.rows() != Binning->GetNBins() ||
//       relative_cvmx.cols() != Binning->GetNBins()) {
//     MACH3LOG_ERROR("Covariance matrix ({}x{}) wrong dimensions, "
//                     "total number of bins: {}",
//                     relative_cvmx.rows(), relative_cvmx.cols(), Binning->GetNBins());
//     throw MaCh3Exception(__FILE__, __LINE__);
//   }

//   const int n = static_cast<int>(mc.size());
//   constexpr double mc_floor = 1e-3; // treat mc(i) below this as a dead bin

//   cvmx = relative_cvmx;
//   Eigen::VectorXd residual = Eigen::VectorXd::Zero(n);

//   for (int i = 0; i < n; ++i) {
//     if (mc(i) > mc_floor) {
//       cvmx(i, i) += (unweighted_mc(i) + data(i)) / (mc(i) * mc(i));
//       residual(i) = (data(i) - mc(i)) / mc(i);
//     }
//   }

//   for (int i = 0; i < n; ++i) {
//     if (cvmx(i, i) <= 0.0 || mc(i) <= mc_floor) {
//       cvmx(i, i) = 1.0;
//       residual(i) = 0.0;
//     }
//   }

//   Eigen::LLT<Eigen::MatrixXd> llt(cvmx);
//   if (llt.info() != Eigen::Success) {
//     MACH3LOG_ERROR("Covariance matrix is not positive-definite at this step.");
//     throw MaCh3Exception(__FILE__, __LINE__);
//   }

//   Eigen::VectorXd solved = llt.solve(residual); // == cvmx^{-1} * residual
//   double chi2 = residual.dot(solved);
  
//   // ---------------- DEBUG OUTPUT ----------------
//   Eigen::MatrixXd cvmx_inv = llt.solve(Eigen::MatrixXd::Identity(n, n));

//   std::cout << "===== GetLikelihood() debug (relative-space) =====\n";
//   std::cout << "cvmx (relative, regularized):\n" << cvmx << "\n\n";
//   std::cout << "cvmx inverse:\n" << cvmx_inv << "\n\n";

//   std::cout << std::left
//             << std::setw(6)  << "bin"
//             << std::setw(14) << "mc"
//             << std::setw(14) << "data"
//             << std::setw(16) << "residual(rel)"
//             << std::setw(20) << "(cov^-1 * r)_i"
//             << std::setw(16) << "bin_chi2"
//             << "\n";
//   double running_total = 0.0;
//   for (int i = 0; i < n; ++i) {
//     double bin_chi2 = residual(i) * solved(i);
//     running_total += bin_chi2;
//     std::cout << std::left
//               << std::setw(6)  << i
//               << std::setw(14) << mc(i)
//               << std::setw(14) << data(i)
//               << std::setw(16) << residual(i)
//               << std::setw(20) << solved(i)
//               << std::setw(16) << bin_chi2
//               << "\n";
//   }
//   std::cout << "sum of per-bin contributions = " << running_total << "\n";
//   std::cout << "total chi2 (residual.dot(solved)) = " << chi2 << "\n";
//   std::cout << "====================================================\n";
//   // ------------------------------------------------

//   return 0.5 * chi2;
// }








};
}


