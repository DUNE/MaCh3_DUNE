// Write the binary spline cache used by FastAtmSplineHandler, without building the samples.
//
//   MakeAtmSplineCache <EigenFile> <SplineFile.root> <params.txt|config.yaml> <out.m3spl> [MD5 of SplineFile] [threads]
//
// - eventIDs are read from the last column (Variables::eventID) of the Eigen file (column-major, so one contiguous read);
// - the parameter list is a text file with one ParameterName per line, or a covariance YAML: every uncommented
//   Systematic with Type: Spline is taken (e.g. Configs/CovObjs/xsec_covariance_Atmospherics.yaml for full57);
// - the knot values are read exactly as the original loader does (knot sorting, non-finite values -> 1.0), see FastAtmSplineHandler.h;
// - the MD5 (if given) is checked against the ROOT file (streaming) and stored in the cache; jobs then check
//   InputSplinesMD5Sum against it instead of reading the 23.5 GB file.
// Memory: ~4 bytes per event + the cache size (~3.8 GB for full57); time: a few minutes with 8 threads.
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "Splines/FastAtmSplineHandler.h"
#include "Samples/DUNEUtils.h"
#include "Manager/Manager.h"

int main(int argc, char* argv[]) {
  if (argc < 5) {
    std::cerr << "Usage: " << argv[0] << " <EigenFile> <SplineFile.root> <params.txt|cov.yaml> <out.m3spl> [MD5] [threads]" << std::endl;
    return 1;
  }
  const std::string eigenFile = argv[1], splineFile = argv[2], parFile = argv[3], outFile = argv[4];
  const std::string md5 = argc > 5 ? argv[5] : "";
  const int threads = argc > 6 ? std::atoi(argv[6]) : 0;

  // ---- parameter names
  std::vector<std::string> names;
  if (parFile.size() > 5 && (parFile.substr(parFile.size() - 5) == ".yaml" || parFile.substr(parFile.size() - 4) == ".yml")) {
    YAML::Node cfg = YAML::LoadFile(parFile);
    for (const auto& s : cfg["Systematics"]) {
      const auto& sys = s["Systematic"];
      if (sys["Type"] && sys["Type"].as<std::string>() == "Spline") names.push_back(sys["Names"]["ParameterName"].as<std::string>());
    }
  } else {
    std::ifstream in(parFile);
    std::string n;
    while (in >> n) names.push_back(n);
  }
  if (names.empty()) { MACH3LOG_ERROR("No spline parameters found in {}", parFile); return 1; }
  MACH3LOG_INFO("{} spline parameters from {}", names.size(), parFile);

  // ---- eventIDs from the Eigen file
  std::ifstream ef(eigenFile, std::ios::binary);
  if (!ef) { MACH3LOG_ERROR("Cannot open {}", eigenFile); return 1; }
  int64_t rows = 0, cols = 0;
  ef.read(reinterpret_cast<char*>(&rows), 8);
  ef.read(reinterpret_cast<char*>(&cols), 8);
  // eventID is the last column (SampleHandlerAtm::Variables): 15 in the 16-column file, 16 once inelasticity_reco was added
  if (cols < 16) { MACH3LOG_ERROR("Eigen file has {} columns, expected >= 16", cols); return 1; }
  const int64_t kEventID = cols - 1;
  std::vector<double> col(rows);
  ef.seekg(16 + kEventID*rows*8, std::ios::beg);
  ef.read(reinterpret_cast<char*>(col.data()), rows*8);
  if (!ef) { MACH3LOG_ERROR("Short read from {}", eigenFile); return 1; }
  std::vector<uint32_t> ids(rows);
  for (int64_t i = 0; i < rows; ++i) ids[i] = static_cast<uint32_t>(col[i]);
  std::vector<double>().swap(col);
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
  MACH3LOG_INFO("{}: {} events, {} unique eventIDs", eigenFile, rows, ids.size());

  if (!md5.empty()) {
    if (!CompareMD5Sum(md5, splineFile)) return 1;
  }
  FastAtmSplineHandler::RawSplines raw;
  FastAtmSplineHandler::ReadFromRoot({splineFile}, names, ids, threads, raw);
  raw.SourceMD5 = md5;
  for (size_t p = 0; p < raw.ParNames.size(); ++p) {
    if (raw.NonFinite[p]) MACH3LOG_WARN("  - {}: {} non-finite (NaN/inf) spline values replaced with 1.0", raw.ParNames[p], raw.NonFinite[p]);
  }
  FastAtmSplineHandler::WriteCache(outFile, raw);
  MACH3LOG_INFO("Wrote {} ({} events, {} splines, {} values, {:.2f} GB)", outFile, raw.EventIDs.size(), raw.Par.size(), raw.Delta.size(),
                (raw.Delta.size()*8.0 + raw.Par.size() + raw.Count.size() + raw.EventIDs.size()*4.0)/1e9);
  return 0;
}
