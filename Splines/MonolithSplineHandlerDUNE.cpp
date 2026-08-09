#include "MonolithSplineHandlerDUNE.h"

#include "TFile.h"
#include "TBranch.h"
#include "TObjArray.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

std::string ResolveSplineBranchName(
    TTree* tree,
    const std::string& requestedName)
{
    if (tree == nullptr) {
        MACH3LOG_ERROR(
            "ResolveSplineBranchName received a null tree"
        );
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    // First try an exact match.
    //
    // Example:
    // CCQETemplateReweight_SBN_v3_HFToCRPA_q0bin0
    if (tree->GetBranch(requestedName.c_str()) != nullptr) {
        return requestedName;
    }

    // Otherwise try suffix matching.
    //
    // Example:
    // requestedName = AhtBY
    //
    // branch =
    // GENIEReWeight_ICARUS_v2_multisigma_AhtBY
    const std::string suffix = "_" + requestedName;

    std::vector<std::string> matches;

    TObjArray* branches = tree->GetListOfBranches();

    if (branches == nullptr) {
        MACH3LOG_ERROR(
            "Could not retrieve the branch list from tree '{}'",
            tree->GetName()
        );
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    for (int i = 0; i < branches->GetEntries(); ++i) {
        const auto* branch =
            dynamic_cast<const TBranch*>(branches->At(i));

        if (branch == nullptr) {
            continue;
        }

        const std::string candidate = branch->GetName();

        if (candidate.size() < suffix.size()) {
            continue;
        }

        const bool matchesSuffix =
            candidate.compare(
                candidate.size() - suffix.size(),
                suffix.size(),
                suffix
            ) == 0;

        if (matchesSuffix) {
            matches.push_back(candidate);
        }
    }

    if (matches.size() == 1) {
        MACH3LOG_INFO(
            "Mapped requested spline '{}' to SystWeights branch '{}'",
            requestedName,
            matches.front()
        );

        return matches.front();
    }

    if (matches.empty()) {
        MACH3LOG_ERROR(
            "No SystWeights branch matches requested spline '{}'",
            requestedName
        );
    } else {
        MACH3LOG_ERROR(
            "Requested spline '{}' ambiguously matches {} branches:",
            requestedName,
            matches.size()
        );

        for (const auto& match : matches) {
            MACH3LOG_ERROR("  {}", match);
        }
    }

    throw MaCh3Exception(__FILE__, __LINE__);
}


int GetBranchIndex(
    TTree* tree,
    const std::string& branchName)
{
    if (tree == nullptr) {
        MACH3LOG_ERROR(
            "GetBranchIndex received a null tree"
        );
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    TObjArray* branches = tree->GetListOfBranches();

    if (branches == nullptr) {
        MACH3LOG_ERROR(
            "Could not retrieve branches from tree '{}'",
            tree->GetName()
        );
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    for (int i = 0; i < branches->GetEntries(); ++i) {
        const auto* branch =
            dynamic_cast<const TBranch*>(branches->At(i));

        if (branch == nullptr) {
            continue;
        }

        if (branchName == branch->GetName()) {
            return i;
        }
    }

    MACH3LOG_ERROR(
        "Branch '{}' was not found in tree '{}'",
        branchName,
        tree->GetName()
    );

    throw MaCh3Exception(__FILE__, __LINE__);
}

} // namespace

MonolithSplineHandlerDUNE::MonolithSplineHandlerDUNE(
    const std::vector<SplineParameter>& splinePars, 
    const std::vector<uint>& eventIndices, 
    const std::string& spline_filename)
: MonolithSplineHandlerDUNE(GetInitParamsFromConfig(splinePars, eventIndices, spline_filename))
{
}

MonolithSplineHandlerDUNE::MonolithSplineHandlerDUNE(
    std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>> initParams)
: UnbinnedSplineHandler(initParams.first, initParams.second)
{
}

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
    
    TTree* systWeightsTree = reader.GetTree();

    if (systWeightsTree == nullptr) {
        MACH3LOG_ERROR(
            "Could not retrieve SystWeights tree from '{}'",
            spline_filename
        );
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    TObjArray* allBranches =
        systWeightsTree->GetListOfBranches();

    if (allBranches == nullptr) {
        MACH3LOG_ERROR(
            "Could not retrieve SystWeights branch list"
        );
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    // The tree contains ordinary event branches before the systematic branches.
    //
    // For your file:
    //
    // total branches - header entries = 14
    const int totalTreeBranches =
        allBranches->GetEntries();

    const int systematicBranchOffset =
        totalTreeBranches -
        static_cast<int>(fileSplinePars.size());

    if (systematicBranchOffset < 0) {
        MACH3LOG_ERROR(
            "Invalid spline file: SystWeights contains {} branches, "
            "but systsHeader contains {} entries",
            totalTreeBranches,
            fileSplinePars.size()
        );
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    MACH3LOG_INFO(
        "Detected {} non-systematic branches before {} systematic branches",
        systematicBranchOffset,
        fileSplinePars.size()
    );

    //Iterate over spline parameters and print their names
    for (const SplineParameter& splineParam : splinePars) {
        // Full unique name supplied by the YAML.
        const std::string requestedSplineName =
            splineParam.name;

        // Find the unique SystWeights branch.
        const std::string branchName =
            ResolveSplineBranchName(
                systWeightsTree,
                requestedSplineName
            );

        // Index among every branch in SystWeights.
        const int treeBranchIndex =
            GetBranchIndex(
                systWeightsTree,
                branchName
            );

        // Convert the full tree-branch index into the corresponding
        // systsHeader entry index.
        const int headerIndex =
            treeBranchIndex - systematicBranchOffset;

        if (headerIndex < 0 ||
            static_cast<std::size_t>(headerIndex) >=
                fileSplinePars.size()) {
            MACH3LOG_ERROR(
                "Invalid spline mapping for '{}': "
                "branch='{}', tree branch index={}, offset={}, "
                "calculated header index={}, header entries={}",
                requestedSplineName,
                branchName,
                treeBranchIndex,
                systematicBranchOffset,
                headerIndex,
                fileSplinePars.size()
            );

            throw MaCh3Exception(__FILE__, __LINE__);
        }

        // The branch ordering and systsHeader ordering correspond after
        // removing the leading non-systematic branches.
        const SplineHeader& header =
            fileSplinePars[
                static_cast<std::size_t>(headerIndex)
            ];

        if (header.knots.empty()) {
            MACH3LOG_ERROR(
                "systsHeader entry {} for branch '{}' contains no knots",
                headerIndex,
                branchName
            );

            throw MaCh3Exception(__FILE__, __LINE__);
        }

        MACH3LOG_INFO(
            "Spline mapping: requested='{}', branch='{}', "
            "treeBranchIndex={}, headerIndex={}, "
            "headerName='{}', CV={}, knots={}",
            requestedSplineName,
            branchName,
            treeBranchIndex,
            headerIndex,
            header.name,
            header.cv,
            header.knots.size()
        );

        TTreeReaderArray<double> splineValues(
            reader,
            branchName.c_str()
        );

        splineValuesPtrs.emplace_back(
            std::move(splineValues)
        );

        splines.emplace_back(
            new ReusableSpline(header.knots)
        );
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
                    if (std::isnan(val)) {
                        MACH3LOG_WARN("NaN detected in spline values for parameter {} at event index {} - replacing with 1.0 (flat spline)", 
                                      splinePars[i].name, currentEventIdx);
                        values.push_back(1.0);
                    } else {
                        values.push_back(val);
                    }
                }
                
                splines[i]->SetVariationWeights(values);

                splines[i]->Refresh();
                TSpline3* splinePtr = splines[i];
                if (splines[i]->IsFlat()) {
                    currentEventSplines.emplace_back(nullptr);
                } else {
                    currentEventSplines.emplace_back(new TSpline3_redDUNE(splinePtr));
                }
            }
            splinesReduced.emplace_back(std::move(currentEventSplines));
            splineTypes.push_back(RespFuncType::kTSpline3_red);
            
            // Handle duplicates: advance to next different event index or end
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

            // Print progress bar every 5%
            size_t currentPercent = (nextEventToFindIdx * 100) / totalEvents;
            if (currentPercent >= lastPrintPercent + 5 || nextEventToFindIdx == totalEvents) {
                M3::Utils::PrintProgressBar(nextEventToFindIdx, totalEvents);
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
            genericEventSplines.push_back(static_cast<TResponseFunction_red*>(spline));
        }
        splinesReducedGeneric.push_back(std::move(genericEventSplines));
    }

    return std::make_pair(splinesReducedGeneric, splineTypes);

}

std::vector<struct SplineHeader>
MonolithSplineHandlerDUNE::GetSplineParametersFromFile(
    TFile* f,
    const std::string& treeName)
{
    std::vector<struct SplineHeader> splineParams;

    if (f == nullptr || f->IsZombie()) {
        MACH3LOG_ERROR(
            "Invalid ROOT file supplied to GetSplineParametersFromFile"
        );
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    TTree* headerTree =
        dynamic_cast<TTree*>(f->Get(treeName.c_str()));

    if (headerTree == nullptr) {
        MACH3LOG_ERROR(
            "Could not find header tree '{}' in file '{}'",
            treeName,
            f->GetName()
        );
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    TTreeReader reader(headerTree);

    // In the new file, name is one std::string per header entry.
    TTreeReaderValue<std::string> splineName(
        reader,
        "name"
    );

    TTreeReaderValue<double> splineCV(
        reader,
        "cv"
    );

    TTreeReaderValue<bool> splineIsCorrection(
        reader,
        "isCorrection"
    );

    TTreeReaderArray<double> splineKnots(
        reader,
        "variations"
    );

    while (reader.Next()) {
        SplineHeader param;

        param.name = *splineName;
        param.cv = *splineCV;
        param.isCorrection = *splineIsCorrection;

        for (const double knot : splineKnots) {
            param.knots.push_back(knot);
        }

        splineParams.push_back(std::move(param));
    }

    MACH3LOG_INFO(
        "Read {} spline entries from header tree '{}'",
        splineParams.size(),
        treeName
    );

    if (splineParams.empty()) {
        MACH3LOG_ERROR(
            "No spline entries were read from header tree '{}'",
            treeName
        );
        throw MaCh3Exception(__FILE__, __LINE__);
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
