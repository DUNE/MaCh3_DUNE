#pragma once

#include "Splines/BinnedSplineHandlerDUNE.h"

#include "Samples/SampleHandlerFD.h"

#include "Samples/BeamOffAxis/EventInfo.h"
#include "Samples/BeamOffAxis/Projections.h"
#include "Samples/BeamOffAxis/Utility.h"

_MaCh3_Safe_Include_Start_ //{
#include "Eigen/Dense"
    _MaCh3_Safe_Include_End_ //}

    namespace dune::beamoffaxis {

  /// @brief Base class for handling FD Beam samples
  class SampleHandlerBeamOffAxis : public SampleHandlerFD {
  public:
    /// @brief SampleHandler FD beam Constructor
    /// @param mc_version Config Name
    /// @param xsec_cov Cross-section covariance matrix
    /// @param osc_cov Oscillation covariance matrix
    /// @param Oscillator_ Shared Oscillation Handler object
    SampleHandlerBeamOffAxis(
        std::string mc_version, ParameterHandlerGeneric *xsec_cov,
        const std::shared_ptr<OscillationHandler> &Oscillator);

    /// @brief destructor
    ~SampleHandlerBeamOffAxis() {}

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

  protected:
    /// @brief Initialises object
    void Init() override;

    /// @brief Function to setup MC from file
    /// @return Total number of events
    int SetupExperimentMC() override;

    /// @brief Tells FD base which variables to point to/be set to
    void SetupFDMC() override;

    void AddAdditionalWeightPointers() override;
    void SetupSplines() override;
    void RegisterFunctionalParameters() override;
    void ResetShifts(int iEvent) override;
    void FinaliseShifts(int iEvent) override;

    /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
    /// @param KinematicVariable Kinematic parameter Type
    /// @param iEvent Event ID
    /// @return Value of kinematic parameter corresponding for a given event
    double ReturnKinematicParameter(dune::beamoffaxis::KinematicTypes KinPar,
                                    int iEvent);

    /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
    /// @param KinematicVariable Kinematic parameter ID as int
    /// @param iEvent Event ID
    /// @return Value of kinematic parameter corresponding for a given event
    double ReturnKinematicParameter(int KinematicVariable, int iEvent);

    /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
    /// @param KinematicParameter Kinematic parameter name as string (gets cast
    /// -> int)
    /// @param iEvent Event ID
    /// @return Value of kinematic parameter corresponding for a given event
    double ReturnKinematicParameter(std::string KinematicParameter, int iEvent);

    /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
    /// @param KinPar Kinematic Parameter Type
    /// @param iEvent Event ID
    /// @return Pointer to KinPar for a given event
    const double *
    GetPointerToKinematicParameter(dune::beamoffaxis::KinematicTypes KinPar,
                                   int iEvent);

    /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
    /// @param KinematicParameter Kinematic parameter name as string (gets cast
    /// -> int)
    /// @param iEvent Event ID
    /// @return Pointer to KinPar for a given event
    const double *GetPointerToKinematicParameter(std::string KinematicParameter,
                                                 int iEvent);

    /// @brief Returns pointer to kinemtatic parameter for event in Structs DUNE
    /// @param KinematicVariable Kinematic parameter as double (gets cast ->
    /// int)
    /// @param iEvent Event ID
    /// @return Pointer to KinPar for a given event
    const double *GetPointerToKinematicParameter(double KinematicVariable,
                                                 int iEvent);

    double GetLikelihood() const override {

      if (!icvmx.size()) {
        return SampleHandlerFD::GetLikelihood();
      }

      Eigen::Map<Eigen::VectorXd const> data(SampleHandlerFD_data.data(),
                                             SampleHandlerFD_data.size());
      Eigen::Map<Eigen::VectorXd const> mc(SampleHandlerFD_array.data(),
                                           SampleHandlerFD_array.size());

      if (icvmx.rows() != data.size()) {
        MACH3LOG_ERROR("Inverse covariance matrix ({}x{}) is not correct for "
                       "data array size: {}",
                       icvmx.rows(), icvmx.cols(), data.rows());
        throw MaCh3Exception(__FILE__, __LINE__);
      }

      double lh = (data - mc).transpose() * icvmx * (data - mc);

      return lh;
    };

    std::vector<dune::beamoffaxis::EventInfo> DUNEMCEvents;

    std::vector<double> subsample_analysispot;
    std::vector<bool> subsample_is_numode;

    Eigen::MatrixXd cvmx, icvmx;

    void CleanMemoryBeforeFit() {}
  };

} // namespace dune::beamoffaxis
