#include "MonolithSplineHandlerDUNE.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include <spdlog/fmt/ranges.h>

namespace {
  bool SameHeader(const SplineHeader& a, const SplineHeader& b) {
    return a.name == b.name && a.cv == b.cv &&
      a.isCorrection == b.isCorrection && a.knots == b.knots;
  }
}

MonolithSplineHandlerDUNE::MonolithSplineHandlerDUNE(const std::vector<SplineParameter>& splinePars, const std::vector<uint>& eventIndices, const std::vector<std::string>& spline_filename)
: MonolithSplineHandlerDUNE(GetInitParamsFromConfig(splinePars, eventIndices, spline_filename)) {}

MonolithSplineHandlerDUNE::MonolithSplineHandlerDUNE(std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>> initParams)
: UnbinnedSplineHandler(initParams.first, initParams.second) {}

MonolithSplineHandlerDUNE::MonolithSplineHandlerDUNE(std::vector<std::vector<TResponseFunction_red*> >& splinesPerEvent, const std::vector<RespFuncType>& splineTypes)
: UnbinnedSplineHandler(splinesPerEvent, splineTypes) {}

std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>>
MonolithSplineHandlerDUNE::GetInitParamsFromConfig(const std::vector<SplineParameter>& splinePars, const std::vector<uint>& eventIndices, const std::vector<std::string>& spline_filename) {

  TChain weights("SystWeights");
  for (size_t iFile=0;iFile<spline_filename.size();iFile++) {
    weights.Add(spline_filename[iFile].c_str());
  }
  
  if (!weights.GetListOfFiles() || weights.GetListOfFiles()->GetEntries() == 0) {
    MACH3LOG_ERROR("No SystWeights input files matched: {}", fmt::join(spline_filename, ", "));
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  // Validate every file independently. In particular, do not concatenate
  // repeated header rows and mistake them for distinct parameters.
  std::vector<SplineHeader> fileSplinePars;
  std::map<std::string, SplineHeader> uniqueHeaders;
  auto* files = weights.GetListOfFiles();
  for (int i = 0; i < files->GetEntries(); ++i) {
    auto* element = dynamic_cast<TChainElement*>(files->At(i));
    
    if (!element) {
      MACH3LOG_ERROR("Invalid TChain file element at index {}", i);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    const std::string filename = element->GetTitle();
    std::unique_ptr<TFile> file(TFile::Open(filename.c_str(), "READ"));
    if (!file || file->IsZombie()) {
      MACH3LOG_ERROR("Could not open spline input: {}", filename);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    if (!file->Get("SystWeights") || !file->Get("systsHeader")) {
      MACH3LOG_ERROR("Missing SystWeights or systsHeader in {}", filename);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    const auto headers = GetSplineParametersFromFile(file.get(), "systsHeader");
    if (headers.empty()) {
      MACH3LOG_ERROR("Empty systsHeader in {}", filename);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    std::set<std::string> namesInFile;
    for (const auto& header : headers) {

      if (!namesInFile.insert(header.name).second) {
	MACH3LOG_ERROR("Duplicate spline header '{}' within {}", header.name, filename);
	throw MaCh3Exception(__FILE__, __LINE__);
      }
      
      auto result = uniqueHeaders.emplace(header.name, header);
      if (!result.second && !SameHeader(result.first->second, header)) {
	MACH3LOG_ERROR("Conflicting spline header '{}' in {}", header.name, filename);
	throw MaCh3Exception(__FILE__, __LINE__);
      }
      
      if (i == 0) fileSplinePars.push_back(header);
    }
    
    if (namesInFile.size() != fileSplinePars.size()) {
      MACH3LOG_ERROR("Header parameter count differs in {}", filename);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    for (const auto& expected : fileSplinePars) {
      if (!namesInFile.count(expected.name)) {
	MACH3LOG_ERROR("Missing header '{}' in {}", expected.name, filename);
	throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
    
  }
  
  size_t MaxSizeOfName = 30;
  for (const auto& param : fileSplinePars) {
    MaxSizeOfName = std::max(MaxSizeOfName,param.name.length());
  }
  
  // Total width of the table (including borders and padding)
  const std::size_t TableWidth = MaxSizeOfName + 39;
  std::string Border;
  for (std::size_t i = 0; i < TableWidth; ++i) {
    Border += "═";
  }
  MACH3LOG_INFO("Available splines in input: {} ({} files)", spline_filename, files->GetEntries());
   MACH3LOG_INFO("╔{}╗", Border);
  MACH3LOG_INFO("║ {:^{}} ║", "Spline Parameters from File", TableWidth - 2);
  MACH3LOG_INFO("╠{}╣", Border);
  MACH3LOG_INFO("║ {:<{}} │ {:>10} │ {:>10} │ {:>8} ║", "Name", MaxSizeOfName, "CV", "Knots", "IsCor");
  MACH3LOG_INFO("╠{}╣", Border);  
  for (const auto& param : fileSplinePars) {
    MACH3LOG_INFO("║ {:<{}} │ {:>10.3f} │ {:>10} │ {:>8} ║", param.name, MaxSizeOfName, param.cv, param.knots.size(), param.isCorrection ? "true" : "false");
  }
  MACH3LOG_INFO("╚{}╝", Border);
  
  const Long64_t entries = weights.GetEntries();
  if (entries < 0) {
    MACH3LOG_ERROR("Could not count SystWeights entries in {}", spline_filename);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (!std::is_sorted(eventIndices.begin(), eventIndices.end())) {
    MACH3LOG_ERROR("eventIndices must be sorted for sequential TChain reading");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  if (!eventIndices.empty() && static_cast<Long64_t>(eventIndices.back()) >= entries) {
    MACH3LOG_ERROR("Requested event index {} exceeds SystWeights chain size {}", eventIndices.back(), entries);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  // LoadTree makes the first underlying TTree available for branch checks.
  if (entries > 0 && weights.LoadTree(0) < 0) {
    MACH3LOG_ERROR("Unable to load first SystWeights tree");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  TTreeReader reader(&weights);
  std::vector<std::unique_ptr<TTreeReaderArray<double>>> splineValuesPtrs;
  std::vector<std::unique_ptr<ReusableSpline>> splines;

  for (const SplineParameter& splineParam : splinePars) {
    const std::string& name = splineParam.name;
    if (!weights.GetBranch(name.c_str())) {
      MACH3LOG_ERROR("Desired spline branch '{}' not found in {}", name, spline_filename);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    const auto it = std::find_if(fileSplinePars.begin(), fileSplinePars.end(), [&](const SplineHeader& h) { return h.name == name; });
    
    if (it == fileSplinePars.end()) {
      MACH3LOG_ERROR("Desired spline '{}' not found in headers", name);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    
    splineValuesPtrs.emplace_back(new TTreeReaderArray<double>(reader, name.c_str()));
    splines.emplace_back(new ReusableSpline(it->knots));
  }
  
  Long64_t currentEventIdx = 0;
  size_t nextEventToFindIdx = 0;
  std::vector<std::vector<TSpline3_redDUNE*> > splinesReduced;
  std::vector<RespFuncType> splineTypes;

  MACH3LOG_INFO("Loading splines for {} events from {} chain entries...", eventIndices.size(), entries);
  const size_t totalEvents = eventIndices.size();
  size_t lastPrintPercent = 0;
  std::map<std::string, size_t> nanCounts;
  
  while (nextEventToFindIdx < eventIndices.size() && reader.Next()) {
    if (currentEventIdx == static_cast<Long64_t>(eventIndices[nextEventToFindIdx])) {
      std::vector<TSpline3_redDUNE*> currentEventSplines;

      for (size_t i = 0; i < splines.size(); ++i) {
	const auto& splineValues = *splineValuesPtrs[i];
	std::vector<double> values;
	for (auto val : splineValues) {
	  if (std::isnan(val)) {
	    const std::string& paramName = splinePars[i].name;
	    if (nanCounts[paramName] == 0) {
	      MACH3LOG_WARN("First NaN in spline '{}' at global event index {}; replacing with 1.0", paramName, currentEventIdx);
	    }
	    ++nanCounts[paramName];
	    values.push_back(1.0);
	  } else {
	    values.push_back(val);
	  }
	}
	
	splines[i]->SetVariationWeights(values);
	splines[i]->Refresh();
	TSpline3* splinePtr = splines[i].get();
	currentEventSplines.emplace_back(splines[i]->IsFlat() ? nullptr : new TSpline3_redDUNE(splinePtr));
      }
      
      splinesReduced.emplace_back(std::move(currentEventSplines));
      splineTypes.push_back(RespFuncType::kTSpline3_red);
      ++nextEventToFindIdx;

      while (nextEventToFindIdx < eventIndices.size() && static_cast<Long64_t>(eventIndices[nextEventToFindIdx]) == currentEventIdx) {
	std::vector<TSpline3_redDUNE*> duplicateEventSplines;
	for (auto* originalSpline : splinesReduced.back()) {
	  duplicateEventSplines.emplace_back(originalSpline ? new TSpline3_redDUNE(*originalSpline) : nullptr);
	}

        splinesReduced.emplace_back(std::move(duplicateEventSplines));
	splineTypes.push_back(RespFuncType::kTSpline3_red);
	++nextEventToFindIdx;
      }

      const size_t currentPercent = (nextEventToFindIdx * 100) / totalEvents;
      if (currentPercent >= lastPrintPercent + 5 || nextEventToFindIdx == totalEvents) {
	M3::Utils::PrintProgressBar(nextEventToFindIdx, totalEvents);
	lastPrintPercent = currentPercent;
      }
    }

    ++currentEventIdx;
  }

  if (nextEventToFindIdx != eventIndices.size()) {
    MACH3LOG_ERROR("SystWeights chain ended after resolving {} of {} requested indices", nextEventToFindIdx, eventIndices.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  MACH3LOG_INFO("Total number of events processed for splines: {}", splinesReduced.size());
  
  if (!nanCounts.empty()) {
    MACH3LOG_WARN("Summary of NaN spline values replaced with 1.0:");
    for (const auto& item : nanCounts) {
      MACH3LOG_WARN("  - {}: {} NaN values replaced", item.first, item.second);
    }
  }
  
  std::vector<std::vector<TResponseFunction_red*> > splinesReducedGeneric;
  for (const auto& eventSplines : splinesReduced) {
    std::vector<TResponseFunction_red*> genericEventSplines;
    for (auto* spline : eventSplines) {
      genericEventSplines.push_back(static_cast<TResponseFunction_red*>(spline));
    }
    splinesReducedGeneric.push_back(std::move(genericEventSplines));
  }

  return std::make_pair(std::move(splinesReducedGeneric), std::move(splineTypes));
}

std::vector<struct SplineHeader> MonolithSplineHandlerDUNE::GetSplineParametersFromFile(TFile* file, const std::string& treeName) {
  std::vector<struct SplineHeader> splineParams;
  
  TTreeReader reader(treeName.c_str(), file);
  
  TTreeReaderValue<std::string> splineNames(reader, "name");
  TTreeReaderValue<double> splineCV(reader, "cv");
  TTreeReaderValue<bool> splineIsCorrection(reader, "isCorrection");
  TTreeReaderArray<double> splineKnots(reader, "variations");

  std::vector<std::string> FoundNames;
  
  while (reader.Next()) {
    std::string splineName = *splineNames;
    
    if (std::find(FoundNames.begin(), FoundNames.end(), splineName) == FoundNames.end()) {
      FoundNames.push_back(splineName);
      
      SplineHeader param;
      
      param.name = *splineNames;
      
      if (param.name.empty()) {
	MACH3LOG_ERROR("Empty spline name in {}", file->GetName());
	throw MaCh3Exception(__FILE__, __LINE__);
      }
      
      param.cv = *splineCV;
      param.isCorrection = *splineIsCorrection;
      
      for (auto knot : splineKnots) {
	param.knots.push_back(knot);
      }
      
      splineParams.push_back(std::move(param));
    }
  }
  
  MACH3LOG_INFO("Loaded {} spline headers from {}",splineParams.size(),file->GetName());
  
  return splineParams;
}

MonolithSplineHandlerDUNE::~MonolithSplineHandlerDUNE() {}

void MonolithSplineHandlerDUNE::InitFromFile(std::string &spline_filename) {
  MACH3LOG_ERROR("MonolithSplineHandlerDUNE::InitFromFile is not implemented. Requested file: {}", spline_filename);
  throw MaCh3Exception(__FILE__, __LINE__);
}
