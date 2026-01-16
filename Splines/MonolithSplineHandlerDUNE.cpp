#include "MonolithSplineHandlerDUNE.h"
#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "Manager/Monitor.h"

MonolithSplineHandlerDUNE::MonolithSplineHandlerDUNE(std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>> initParams) :
    SMonolith(initParams.first, initParams.second)
    { }


MonolithSplineHandlerDUNE::MonolithSplineHandlerDUNE(const std::vector<SplineParameter>& splinePars,
    const std::vector<uint>& eventIndices, const std::string& spline_filename) :
    MonolithSplineHandlerDUNE(GetInitParamsFromConfig(splinePars, eventIndices, spline_filename))
    { }

std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>>
MonolithSplineHandlerDUNE::GetInitParamsFromConfig(
    const std::vector<SplineParameter>& splinePars,
    const std::vector<uint>& eventIndices,
    const std::string& spline_filename) {
    TFile *f = TFile::Open(spline_filename.c_str(),"READ");
    if (!f || f->IsZombie()) {
        MACH3LOG_ERROR("Could not open spline file: {}",spline_filename);
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    std::vector<struct SplineHeader> fileSplinePars = 
        MonolithSplineHandlerDUNE::GetSplineParametersFromFile(f, "systsHeader");

    //Make a nice printout of available splines in the file under the form of an array
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


    const std::string treeName = "SystWeights";
    TTreeReader reader(treeName.c_str(), f);

    //List all available branches (spline names) in the tree

    std::vector<TTreeReaderArray<double>> splineValuesPtrs;
    std::vector<ReusableSpline*> splines;
    
    //Iterate over spline parameters and print their names
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

    while (reader.Next()) {
        if (nextEventToFindIdx >= eventIndices.size()) {
            break;
        }
        if (currentEventIdx == eventIndices[nextEventToFindIdx]) {
            // Process the splines for this event
            std::vector<TSpline3_redDUNE*> currentEventSplines;
            for (size_t i = 0; i < splines.size(); ++i) {
                const auto& splineValues = splineValuesPtrs[i];
                std::vector<double> values;
                for (auto& val : splineValues) {
                    values.push_back(val);
                }
                splines[i]->SetVariationWeights(values);

                splines[i]->Refresh();
                TSpline3* splinePtr = splines[i];
                currentEventSplines.emplace_back(new TSpline3_redDUNE(splinePtr));   
            }
            splinesReduced.emplace_back(std::move(currentEventSplines));
            splineTypes.push_back(RespFuncType::kTSpline3_red);
            
            // Handle duplicates: advance to next different event index or end
            ++nextEventToFindIdx;
            while (nextEventToFindIdx < eventIndices.size() && 
                   eventIndices[nextEventToFindIdx] == currentEventIdx) {
                std::vector<TSpline3_redDUNE*> duplicateEventSplines;
                for (size_t i = 0; i < splines.size(); ++i) {
                    duplicateEventSplines.emplace_back(new TSpline3_redDUNE(*splinesReduced.back()[i]));
                }
                splinesReduced.emplace_back(std::move(duplicateEventSplines));
                splineTypes.push_back(RespFuncType::kTSpline3_red);
                ++nextEventToFindIdx;
            }

            // Print progress bar every 5%
            size_t currentPercent = (nextEventToFindIdx * 100) / totalEvents;
            if (currentPercent >= lastPrintPercent + 5 || nextEventToFindIdx == totalEvents) {
                MaCh3Utils::PrintProgressBar(nextEventToFindIdx, totalEvents);
                lastPrintPercent = currentPercent;
            }
        }
        ++currentEventIdx;
    }
    
    std::cout << std::endl; // New line after progress bar
    MACH3LOG_INFO("Total number of events processed for splines: {}", splinesReduced.size());
    
    std::vector<std::vector<TResponseFunction_red*> > splinesReducedGeneric;
    for (const auto& eventSplines : splinesReduced) {
        std::vector<TResponseFunction_red*> genericEventSplines;
        for (const auto& spline : eventSplines) {
            if(spline->IsFlat()) { //Filtering out flat splines for better performance
                delete spline;
                genericEventSplines.push_back(nullptr);
                continue;
            }
            genericEventSplines.push_back(static_cast<TResponseFunction_red*>(spline));
        }
        splinesReducedGeneric.push_back(std::move(genericEventSplines));
    }

    return std::make_pair(splinesReducedGeneric, splineTypes);

}

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

void MonolithSplineHandlerDUNE::InitFromFile(std::string &spline_filename) {
    // TODO: Implement MonolithSplineHandlerDUNE::InitFromFile to load spline data from file.
    // This method is currently a stub and should not be used in production code.
    MACH3LOG_ERROR("MonolithSplineHandlerDUNE::InitFromFile is not implemented. Requested file: {}", spline_filename);
    throw MaCh3Exception(__FILE__, __LINE__);
}