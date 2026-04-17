#include <memory>
#include "python/pyMaCh3.h"

#include "Samples/SampleHandlerAtm.h"
#include "Samples/SampleHandlerBeamFD.h"
#include "Samples/SampleHandlerBeamND.h"
#include "Samples/SampleHandlerBeamNDGAr.h"
#include "Samples/MaCh3DUNEFactory.h"
#include "Fitters/MaCh3Factory.h"

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
        m_samples.def("MaCh3DuneFactory",
            [](const std::string& config_file, py::args extra_args)
                -> std::pair<
                    std::shared_ptr<ParameterHandlerGeneric>,
                    std::vector<std::shared_ptr<SampleHandlerBase>>
                >
            {
                // Create manager

                // Because of how MaCh3's setup we need to "spoof" it into thinking we're in the command line
                // For now I'm not bothering with override
                std::vector<std::string> args;
                args.emplace_back("mach3");        // (dummy arg)
                args.emplace_back(config_file);    // config file

                for (auto item : extra_args) {
                    args.emplace_back(py::cast<std::string>(item));
                }

                int mock_argc = static_cast<int>(args.size());

                std::vector<char*> mock_argv;
                mock_argv.reserve(args.size());
                for (auto& s : args) {
                    mock_argv.push_back(const_cast<char*>(s.c_str()));
                }

             
                auto FitManager = MaCh3ManagerFactory(mock_argc, mock_argv.data());

                // Create parameter handler
                auto xsec = MaCh3CovarianceFactory<ParameterHandlerGeneric>(
                    FitManager.get(), "Xsec"
                );

                if (CheckNodeExists(FitManager->raw(), "General", "OscillationParameters"))
                {
                    auto oscpars = Get<std::vector<double>>(
                        FitManager->raw()["General"]["OscillationParameters"],
                        __FILE__, __LINE__
                    );
                    xsec->SetGroupOnlyParameters("Osc", oscpars);
                }

                // Call existing factory 
                auto raw_samples = MaCh3DuneSampleFactory(FitManager, xsec);

                // Convert to shared_ptr for Python safety
                std::vector<std::shared_ptr<SampleHandlerBase>> samples;
                samples.reserve(raw_samples.size());

                for (auto* s : raw_samples) {
                    samples.emplace_back(std::shared_ptr<SampleHandlerBase>(s));
                }

                // Move ownership of xsec to shared_ptr
                std::shared_ptr<ParameterHandlerGeneric> xsec_shared = std::move(xsec);

                return {xsec_shared, samples};
            },
            py::arg("config_file"),
            "Create ParameterHandlerGeneric and DUNE SampleHandlers safely"
        );   
    }
};

MAKE_PYMACH3_MDULE( MaCh3DunePyBinder )