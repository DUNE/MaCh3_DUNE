#include <memory>
#include "python/pyMaCh3.h"

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
    }
};

MAKE_PYMACH3_MDULE( MaCh3DunePyBinder )