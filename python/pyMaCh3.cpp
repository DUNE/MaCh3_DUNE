#include <memory>
#include <vector>
#include "python/pyMaCh3.h"

#include <pybind11/stl.h> // needed for std::vector<SampleHandlerBase*> <-> python list conversion

#include "Samples/SampleHandlerAtm.h"
#include "Samples/SampleHandlerBeamFD.h"
#include "Samples/SampleHandlerBeamND.h"
#include "Samples/SampleHandlerBeamNDGAr.h"
#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"


namespace py = pybind11;

class MaCh3DunePyBinder : public MaCh3PyBinder {

  public:

    void initSamplesExperiment(py::module &m_samples){

        std::cout << "Initializing SampleHandlerTutorial bindings... " << std::endl;

        // ####################################################
        // Beam FD
        py::class_<SampleHandlerBeamFD, SampleHandlerBase, SampleHandlerInterface>(m_samples, "SampleHandlerBeamFD")
        // ####################################################
        // Constructor with 2 arguments (no oscillation handler)
            .def(py::init([](const std::string& mc_version, ParameterHandlerGeneric* xsec_cov) {
                return new SampleHandlerBeamFD(mc_version, xsec_cov, nullptr);
            }),
                "Create SampleHandlerTutorial without oscillation handler",
                py::arg("mc_version"),
                py::arg("xsec_cov")
            )
            // Constructor with 3 arguments (with oscillation handler)
            .def(py::init([](const std::string& mc_version, 
                            ParameterHandlerGeneric* xsec_cov,
                            OscillationHandler* osc_cov) {
                std::shared_ptr<OscillationHandler> osc_ptr;
                if (osc_cov != nullptr) {
                    osc_ptr = std::shared_ptr<OscillationHandler>(osc_cov, [](OscillationHandler*){});
                }
                return new SampleHandlerBeamFD(mc_version, xsec_cov, osc_ptr);
            }),
                "Create SampleHandlerTutorial with oscillation handler",
                py::arg("mc_version"),
                py::arg("xsec_cov"),
                py::arg("osc_cov") = nullptr
            );

        // ####################################################
        // Beam ND
        py::class_<SampleHandlerBeamND, SampleHandlerBase, SampleHandlerInterface>(m_samples, "SampleHandlerBeamND")
        // ####################################################
            // Constructor with 3 arguments (with beam handler)
            .def(py::init([](const std::string& mc_version, 
                            ParameterHandlerGeneric* xsec_cov,
                            BeamNDCov beam_nd_cov) {

                return new SampleHandlerBeamND(mc_version, xsec_cov, beam_nd_cov);
            }),
                "Create SampleHandlerBeamND with oscillation handler",
                py::arg("mc_version"),
                py::arg("xsec_cov"),
                py::arg("beam_nd_cov")
            );

        // ####################################################
        // Atmospheric
        py::class_<SampleHandlerAtm, SampleHandlerBase, SampleHandlerInterface>(m_samples, "SampleHandlerAtm")
        // ####################################################
        // Constructor with 2 arguments (no oscillation handler)
            .def(py::init([](const std::string& mc_version, ParameterHandlerGeneric* xsec_cov) {
                return new SampleHandlerAtm(mc_version, xsec_cov, nullptr);
            }),
                "Create SampleHandlerAtm without oscillation handler",
                py::arg("mc_version"),
                py::arg("xsec_cov")
            )
            // Constructor with 3 arguments (with oscillation handler)
            .def(py::init([](const std::string& mc_version, 
                            ParameterHandlerGeneric* xsec_cov,
                            OscillationHandler* osc_cov) {
                std::shared_ptr<OscillationHandler> osc_ptr;
                if (osc_cov != nullptr) {
                    osc_ptr = std::shared_ptr<OscillationHandler>(osc_cov, [](OscillationHandler*){});
                }
                return new SampleHandlerAtm(mc_version, xsec_cov, osc_ptr);
            }),
                "Create SampleHandlerTutorial with oscillation handler",
                py::arg("mc_version"),
                py::arg("xsec_cov"),
                py::arg("osc_cov") = nullptr
            );

        // ####################################################
        // Atmospheric
        py::class_<SampleHandlerBeamNDGAr, SampleHandlerBase, SampleHandlerInterface>(m_samples, "SampleHandlerBeamNDGAr")
        // ####################################################
        // Constructor with 2 arguments (no oscillation handler)
            .def(py::init([](const std::string& mc_version, ParameterHandlerGeneric* xsec_cov) {
                return new SampleHandlerBeamNDGAr(mc_version, xsec_cov);
            }),
                "Create SampleHandlerBeamND gar without oscillation handler",
                py::arg("mc_version"),
                py::arg("xsec_cov")
            );



        // ####################################################
        // Beam ND Cov struct
        py::class_<BeamNDCov>(m_samples, "BeamNDCov")
            .def(py::init<>([](const std::string& nd_cov_file_name, bool use_combined ){
                auto nd_cov_file = M3::Open(nd_cov_file_name, "READ", __FILE__, __LINE__);

                // Grab TMatrices
                auto nd_cov_fhc = nd_cov_file->Get<TMatrixD>("nd_fhc_frac_cov");
                auto nd_cov_rhc = nd_cov_file->Get<TMatrixD>("nd_rhc_frac_cov");
                auto nd_cov_all = nd_cov_file->Get<TMatrixD>("nd_all_frac_cov");

                // Make sure they exists
                if (!(nd_cov_fhc && nd_cov_rhc && nd_cov_all))
                {
                    MACH3LOG_ERROR("Could not find NDCov objects from file: {}", nd_cov_file_name);
                    throw MaCh3Exception(__FILE__, __LINE__);
                }

                // Voila
                BeamNDCov beam_nd_cov;
                beam_nd_cov.NDCov_FHC = nd_cov_fhc;
                beam_nd_cov.NDCov_RHC = nd_cov_rhc;
                beam_nd_cov.NDCov_all = nd_cov_all;
                beam_nd_cov.useCombinedNDCov = use_combined;
                nd_cov_file->Close();
                return beam_nd_cov;
            })

        );

        // ####################################################
        // Factory function: GetMaCh3DuneInstance
        // ####################################################
        m_samples.def("GetMaCh3DuneInstance",
            [](const std::string& SampleType,
               const std::string& SampleConfig,
               ParameterHandlerGeneric* param_handler,
               OscillationHandler* BeamOscillator_,
               OscillationHandler* AtmOscillator_,
               BeamNDCov beamNDCov) -> SampleHandlerBase* {

                // GetMaCh3DuneInstance takes std::unique_ptr<ParameterHandlerGeneric>&,
                // so wrap the raw (Python-owned) pointer in a real unique_ptr, call
                // the function, then release it so we don't double-delete the
                // object when this local goes out of scope.
                std::unique_ptr<ParameterHandlerGeneric> param_handler_ptr(param_handler);

                std::shared_ptr<OscillationHandler> beam_osc_ptr;
                if (BeamOscillator_ != nullptr) {
                    beam_osc_ptr = std::shared_ptr<OscillationHandler>(BeamOscillator_, [](OscillationHandler*){});
                }

                std::shared_ptr<OscillationHandler> atm_osc_ptr;
                if (AtmOscillator_ != nullptr) {
                    atm_osc_ptr = std::shared_ptr<OscillationHandler>(AtmOscillator_, [](OscillationHandler*){});
                }

                SampleHandlerBase* Sample;
                try {
                    Sample = GetMaCh3DuneInstance(
                        SampleType, SampleConfig, param_handler_ptr,
                        beam_osc_ptr, atm_osc_ptr, beamNDCov);
                } catch (...) {
                    param_handler_ptr.release(); // ownership stays with the caller/Python side
                    throw;
                }

                param_handler_ptr.release(); // ownership stays with the caller/Python side

                return Sample;
            },
            "Create a MaCh3 DUNE SampleHandler instance based on SampleType "
            "(one of \"BeamFD\", \"BeamND\", \"Atm\", \"BeamNDGAr\")",
            py::arg("SampleType"),
            py::arg("SampleConfig"),
            py::arg("param_handler"),
            py::arg("BeamOscillator") = nullptr,
            py::arg("AtmOscillator") = nullptr,
            py::arg("beamNDCov") = BeamNDCov(),
            py::return_value_policy::take_ownership
        );

        // ####################################################
        // Factory function: MaCh3DuneSampleFactory
        // ####################################################
        m_samples.def("MaCh3DuneSampleFactory",
            [](Manager* fit_manager, ParameterHandlerGeneric* param_handler) -> std::vector<SampleHandlerBase*> {

                // Same ownership trick as GetMaCh3DuneInstance: MaCh3DuneSampleFactory
                // takes std::unique_ptr<Manager>& and std::unique_ptr<ParameterHandlerGeneric>&,
                // but Python owns these objects as raw pointers. Wrap them temporarily,
                // call the factory, then release before returning so Python retains
                // ownership and nothing gets double-freed.
                std::unique_ptr<Manager> fit_manager_ptr(fit_manager);
                std::unique_ptr<ParameterHandlerGeneric> param_handler_ptr(param_handler);

                std::vector<SampleHandlerBase*> result;
                try {
                    result = MaCh3DuneSampleFactory(fit_manager_ptr, param_handler_ptr);
                } catch (...) {
                    fit_manager_ptr.release();
                    param_handler_ptr.release();
                    throw;
                }

                fit_manager_ptr.release();
                param_handler_ptr.release();

                return result;
            },
            "Build the vector of DUNE SampleHandlers described by General:DUNESamples "
            "in the fit Manager config, using the given parameter handler for systematics "
            "and oscillation parameters",
            py::arg("fit_manager"),
            py::arg("param_handler"),
            py::return_value_policy::take_ownership
        );
    }
};

MAKE_PYMACH3_MDULE( MaCh3DunePyBinder )