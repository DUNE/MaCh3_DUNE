#include "SplineHandlerFactoryDUNE.h"

SplineHandlerFactoryDUNE::SplineHandlerFactoryDUNE(ParameterHandlerGeneric* xsec_params,
    MaCh3Modes* Modes_, const std::vector<dunemc_base>& dunemcSamples,
    const std::string& inputSplinesFile, const std::string& sampleName)
: fXsecParams(xsec_params)
{
    YAML::Node ConfigCurrent = fXsecParams->GetConfig();
    fSplineType = kUndef;  // Default spline type
    //Get the systematics spline type
    for (const auto& systNode : ConfigCurrent["Systematics"]) {
        YAML::Node systDetails = systNode["Systematic"];
        if (systDetails && systDetails["Names"]) {
            std::string systName = systDetails["Names"]["ParameterName"].as<std::string>();
        }
 
        if (systDetails["SplineInformation"] && systDetails["SplineInformation"]["Type"]) {
            std::string splineTypeCurrent = systDetails["SplineInformation"]["Type"].as<std::string>();
            if (splineTypeCurrent == "Binned") {
                if (fSplineType != kUndef && fSplineType != kBinned) {
                    MACH3LOG_ERROR("Multiple spline types found in configuration. Unsupported yet.");
                    throw MaCh3Exception(__FILE__, __LINE__);
                }
                fSplineType = kBinned;
            } else if (splineTypeCurrent == "PerEvent") {
                if (fSplineType != kUndef && fSplineType != kMonolith) {
                    MACH3LOG_ERROR("Multiple spline types found in configuration. Unsupported yet.");
                    throw MaCh3Exception(__FILE__, __LINE__);
                }
                fSplineType = kMonolith;
            } else {
                MACH3LOG_ERROR("Unknown spline type: {}", splineTypeCurrent);
                throw MaCh3Exception(__FILE__, __LINE__);
            }
        }
    }

    if (fSplineType == kUndef) {
        MACH3LOG_WARN("No spline type defined in configuration. Defaulting to Binned.");
        fSplineType = kBinned;
    }
    switch (fSplineType) {
        case kBinned:
            MACH3LOG_INFO("Using BinnedSplineHandlerDUNE");
            if (!inputSplinesFile.empty()) {
                MACH3LOG_WARN("Input spline file provided but will not be used for BinnedSplineHandlerDUNE");
            }
            fSplineHandler = std::make_unique<BinnedSplineHandlerDUNE>(xsec_params, Modes_);
            break;
        case kMonolith:
            {
                MACH3LOG_INFO("Using MonolithSplineHandlerDUNE");
                if (inputSplinesFile.empty()) {
                    MACH3LOG_ERROR("Input spline file must be provided for MonolithSplineHandlerDUNE under SampleOptions/InputSplines");
                    throw MaCh3Exception(__FILE__, __LINE__);
                }
                std::vector<uint> eventIndices(dunemcSamples.size());
                for (size_t i = 0; i < dunemcSamples.size(); ++i) {
                    eventIndices[i] = dunemcSamples[i].eid;
                }
                const std::vector<SplineParameter> samplePars = xsec_params->GetSplineParsFromSampleName(sampleName);
                std::unique_ptr<MonolithSplineHandlerDUNE> monolith = std::make_unique<MonolithSplineHandlerDUNE>(samplePars, eventIndices, inputSplinesFile);
                std::vector<const double*> splineWeightPtrs(samplePars.size());
                for(uint i = 0; i < samplePars.size(); ++i) {
                    splineWeightPtrs[i] = xsec_params->RetPointer(samplePars[i].index);
                }
                monolith->setSplinePointers(splineWeightPtrs);
                fSplineHandler = std::move(monolith);
                break;
            }
            
        default:
            MACH3LOG_ERROR("Unsupported spline type encountered.");
            throw MaCh3Exception(__FILE__, __LINE__);
    }

}

std::unique_ptr<SplineBase> SplineHandlerFactoryDUNE::GetSplineHandler() {
    return std::move(fSplineHandler);
}

SplineHandlerFactoryDUNE::~SplineHandlerFactoryDUNE() {
}
