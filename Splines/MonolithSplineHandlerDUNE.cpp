#include "MonolithSplineHandlerDUNE.h"

#include <map>

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

// Delegating constructor: extracts spline data from file before initializing base class
MonolithSplineHandlerDUNE::MonolithSplineHandlerDUNE(
    const std::vector<SplineParameter>& splinePars, 
    const std::vector<uint>& eventIndices, 
    const std::string& spline_filename)
: MonolithSplineHandlerDUNE(GetInitParamsFromConfig(splinePars, eventIndices, spline_filename))
{
}

// Internal constructor that passes processed response functions to UnbinnedSplineHandler
MonolithSplineHandlerDUNE::MonolithSplineHandlerDUNE(
    std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>> initParams)
: UnbinnedSplineHandler(initParams.first, initParams.second)
{
}

// Main initialization method: loads event-by-event spline response curves from a monolith ROOT file
std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>> 
MonolithSplineHandlerDUNE::GetInitParamsFromConfig(
    const std::vector<SplineParameter>& splinePars,
    const std::vector<uint>& eventIndices,
    const std::string& spline_filename) {

    TFile *f = TFile::Open(spline_filename.c_str(), "READ");
    if (!f || f->IsZombie()) {
        MACH3LOG_ERROR("Could not open spline file: {}", spline_filename);
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    // Retrieve systematic parameter headers (knot positions, central values, etc.)
    std::vector<struct SplineHeader> fileSplinePars = 
        MonolithSplineHandlerDUNE::GetSplineParametersFromFile(f, "systsHeader");

    // Print summary table of available systematic splines in the input file
    MACH3LOG_INFO("Available splines in file: {}", spline_filename);
    MACH3LOG_INFO("╔══════════════════════════════════════════════════════════════════════════════╗");
    MACH3LOG_INFO("║ {:^76} ║", "Spline Parameters from File");
    MACH3LOG_INFO("╠══════════════════════════════════════════════════════════════════════════════╣");
    MACH3LOG_INFO("║ {:<40} │ {:>10} │ {:>10} │ {:>8} ║", "Name", "CV", "Knots", "IsCor");
    MACH3LOG_INFO("╠══════════════════════════════════════════════════════════════════════════════╣");
    for (const auto& param : fileSplinePars) {
        MACH3LOG_INFO("║ {:<40} │ {:>10.3f} │ {:>10} │ {:>8} ║", 
                      param.name, param.cv, param.knots.size(), 
                      param.isCorrection ? "true" : "false");
    }
    MACH3LOG_INFO("╚══════════════════════════════════════════════════════════════════════════════╝");
    MACH3LOG_INFO("Total number of splines in file: {}", fileSplinePars.size());

    // Set up tree reader for event-by-event weight variations
    const std::string treeName = "SystWeights";
    TTreeReader reader(treeName.c_str(), f);

    std::vector<TTreeReaderArray<double>> splineValuesPtrs;
    std::vector<ReusableSpline*> splines;
    
    // Verify required branches exist and bind tree reader arrays for each requested parameter
    for (const SplineParameter& splineParam : splinePars) {
        std::string desiredSplineName = splineParam.name;
        if (!reader.GetTree()->GetBranch(desiredSplineName.c_str())) {
            MACH3LOG_ERROR("Desired spline: {} not found in spline file: {}", desiredSplineName, spline_filename);
            throw MaCh3Exception(__FILE__, __LINE__);
        }

        const auto it = std::find_if(fileSplinePars.begin(), fileSplinePars.end(), [&](const SplineHeader& header) {
                return header.name == desiredSplineName;
            });
        if (it == fileSplinePars.end()) {
            MACH3LOG_ERROR("Desired spline: {} not found in spline headers of file: {}", desiredSplineName, spline_filename);
            throw MaCh3Exception(__FILE__, __LINE__);
        }
        const SplineHeader& header = *it;

        TTreeReaderArray<double> splineValues(reader, desiredSplineName.c_str());
        splineValuesPtrs.emplace_back(std::move(splineValues));
        splines.emplace_back(new ReusableSpline(header.knots));
    }

    uint currentEventIdx = 0;
    size_t nextEventToFindIdx = 0;

    std::vector<std::vector<TSpline3_redDUNE*> > splinesReduced;
    std::vector<RespFuncType> splineTypes;

    MACH3LOG_INFO("Loading splines for {} events...", eventIndices.size());
    const size_t totalEvents = eventIndices.size();
    size_t lastPrintPercent = 0;

    std::map<std::string, size_t> nanCounts;

    // Sequential scan through the TTree matching requested event indices
    while (reader.Next()) {
        if (nextEventToFindIdx >= eventIndices.size()) {
            break;
        }
        if (currentEventIdx == eventIndices[nextEventToFindIdx]) {
            // Process weight variations for all parameters for the current event
            std::vector<TSpline3_redDUNE*> currentEventSplines;
            for (size_t i = 0; i < splines.size(); ++i) {
                const auto& splineValues = splineValuesPtrs[i];
                std::vector<double> values;
                for (auto& val : splineValues) {
                    // Replace invalid values with 1.0 (unit weight) to prevent unphysical calculations
                    if (std::isnan(val)) {
                        const std::string& paramName = splinePars[i].name;
                        if (nanCounts[paramName] == 0) {
                            MACH3LOG_WARN("First NaN detected in spline values for parameter {} at event index {} - replacing with 1.0 (flat spline). Further NaNs for this parameter will be counted quietly.", 
                                          paramName, currentEventIdx);
                        }
                        nanCounts[paramName]++;
                        values.push_back(1.0);
                    } else {
                        values.push_back(val);
                    }
                }
                
                splines[i]->SetVariationWeights(values);
                splines[i]->Refresh(); // Recompute cubic spline polynomial coefficients

                // Optimization: store nullptr for flat splines to bypass evaluation during fits
                TSpline3* splinePtr = splines[i];
                if (splines[i]->IsFlat()) {
                    currentEventSplines.emplace_back(nullptr);
                } else {
                    currentEventSplines.emplace_back(new TSpline3_redDUNE(splinePtr));
                }
            }
            splinesReduced.emplace_back(std::move(currentEventSplines));
            splineTypes.push_back(RespFuncType::kTSpline3_red);
            
            // Handle duplicate event entries by cloning reduced spline objects
            ++nextEventToFindIdx;
            while (nextEventToFindIdx < eventIndices.size() && 
                   eventIndices[nextEventToFindIdx] == currentEventIdx) {
                std::vector<TSpline3_redDUNE*> duplicateEventSplines;
                for (size_t i = 0; i < splines.size(); ++i) {
                    auto* originalSpline = splinesReduced.back()[i];
                    if (originalSpline == nullptr) {
                        duplicateEventSplines.emplace_back(nullptr);
                    } else {
                        duplicateEventSplines.emplace_back(new TSpline3_redDUNE(*originalSpline));
                    }
                }
                splinesReduced.emplace_back(std::move(duplicateEventSplines));
                splineTypes.push_back(RespFuncType::kTSpline3_red);
                ++nextEventToFindIdx;
            }

            // Periodically log loading progress
            size_t currentPercent = (nextEventToFindIdx * 100) / totalEvents;
            if (currentPercent >= lastPrintPercent + 5 || nextEventToFindIdx == totalEvents) {
                M3::Utils::PrintProgressBar(nextEventToFindIdx, totalEvents);
                lastPrintPercent = currentPercent;
            }
        }
        ++currentEventIdx;
    }
    
    std::cout << std::endl;
    MACH3LOG_INFO("Total number of events processed for splines: {}", splinesReduced.size());

    if (!nanCounts.empty()) {
        MACH3LOG_WARN("Summary of NaN spline values replaced with 1.0 (flat spline):");
        for (const auto& [paramName, count] : nanCounts) {
            MACH3LOG_WARN("  - {}: {} NaN values replaced", paramName, count);
        }
    }
    
    // Upcast concrete TSpline3_redDUNE pointers to abstract TResponseFunction_red base pointers
    std::vector<std::vector<TResponseFunction_red*> > splinesReducedGeneric;
    for (const auto& eventSplines : splinesReduced) {
        std::vector<TResponseFunction_red*> genericEventSplines;
        for (const auto& spline : eventSplines) {
            genericEventSplines.push_back(static_cast<TResponseFunction_red*>(spline));
        }
        splinesReducedGeneric.push_back(std::move(genericEventSplines));
    }

    return std::make_pair(splinesReducedGeneric, splineTypes);
}

// Reads the 'systsHeader' tree to extract systematic parameter metadata (names, knots, central values)
std::vector<struct SplineHeader> MonolithSplineHandlerDUNE::GetSplineParametersFromFile(TFile *f, const std::string& treeName) {
    std::vector<struct SplineHeader> splineParams;
    TTreeReader reader(treeName.c_str(), f);
    TTreeReaderArray<char> splineNames(reader, "name");
    TTreeReaderValue<double> splineCV(reader, "cv");
    TTreeReaderValue<bool> splineIsCorrection(reader, "isCorrection");
    TTreeReaderArray<double> splineKnots(reader, "variations");

    while (reader.Next()) {
        SplineHeader param;
        param.name = std::string(&splineNames[0]);
        param.cv = *splineCV;
        param.isCorrection = *splineIsCorrection;
        for (auto& knot : splineKnots) {
            param.knots.push_back(knot);
        }
        splineParams.push_back(param);
    }
    return splineParams;
}

MonolithSplineHandlerDUNE::~MonolithSplineHandlerDUNE() {
}

// Stub implementation for file-based initialization (currently unsupported)
void MonolithSplineHandlerDUNE::InitFromFile(std::string &spline_filename) {
    MACH3LOG_ERROR("MonolithSplineHandlerDUNE::InitFromFile is not implemented. Requested file: {}", spline_filename);
    throw MaCh3Exception(__FILE__, __LINE__);
}
