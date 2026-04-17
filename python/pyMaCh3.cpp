#include <memory>
#include "python/pyMaCh3.h"

#include "SamplesTutorial/SampleHandlerAtm.h"
#include "SamplesTutorial/SampleHandlerBeamFD.h"
#include "SamplesTutorial/SampleHandlerBeamND.h"
#include "SamplesTutorial/SampleHandlerBeamNDGar.h"
#include "Samples/MaCh3DUNEFactory.h"

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
        // Beam ND GAr
        py::class_<SampleHandlerBeamNDGAr, SampleHandlerBase, SampleHandlerInterface>(m_samples, "SampleHandlerBeamNDGAr")
        // ####################################################
            .def(py::init([](const std::string& mc_version, 
                            ParameterHandlerGeneric* xsec_cov) {
                return new SampleHandlerBeamNDGAr(mc_version, xsec_cov);
            }),
                "Create SampleHandlerBeamNDGAr",
                py::arg("mc_version"),
                py::arg("xsec_cov")
            );

        // ####################################################
        // Top-level sample factory
        // ####################################################
        m_samples.def("MaCh3DuneSampleFactory",
            [](Manager* fit_manager, ParameterHandlerGeneric* param_handler)
                -> std::vector<SampleHandlerBase*> {

                auto fit_manager_unique = std::unique_ptr<Manager>(fit_manager);
                auto param_handler_unique = std::unique_ptr<ParameterHandlerGeneric>(param_handler);

                auto result = MaCh3DuneSampleFactory(fit_manager_unique, param_handler_unique);

                fit_manager_unique.release();
                param_handler_unique.release();
                return result;
            },
            "Build all DUNE SampleHandlerBase objects declared in the FitManager config",
            py::arg("fit_manager"),
            py::arg("param_handler"),
            py::return_value_policy::take_ownership
        );
    }
};

MAKE_PYMACH3_MDULE( MaCh3DunePyBinder )