#include "Splines/FastAtmSplineHandler.h"
#include "Splines/MonolithSplineHandlerDUNE.h"   // ReusableSpline: the same TSpline3 builder as the original
#include "Samples/DUNEUtils.h"                  // CompareMD5Sum (streaming)

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <thread>
#include <sys/stat.h>

#ifdef MULTITHREAD
#include <omp.h>
#endif

#include "TChain.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TROOT.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"

namespace {
  constexpr char kCacheMagic[9] = "M3FASPL2";
  constexpr char kCacheEnd[9]   = "M3FAEND1";
  constexpr uint32_t kCacheVersion = 2;

  bool FileExists(const std::string& f) {
    struct stat st;
    return stat(f.c_str(), &st) == 0;
  }

  // ---- small binary I/O helpers ----
  void WriteRaw(FILE* f, const void* p, size_t n, const std::string& name) {
    const char* c = static_cast<const char*>(p);
    const size_t chunk = size_t(1) << 28;
    while (n > 0) {
      const size_t w = std::min(n, chunk);
      if (fwrite(c, 1, w, f) != w) {
        MACH3LOG_ERROR("Failed writing spline cache {}", name);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      c += w; n -= w;
    }
  }
  void ReadRaw(FILE* f, void* p, size_t n, const std::string& name) {
    char* c = static_cast<char*>(p);
    const size_t chunk = size_t(1) << 28;
    while (n > 0) {
      const size_t r = std::min(n, chunk);
      if (fread(c, 1, r, f) != r) {
        MACH3LOG_ERROR("Spline cache {} is truncated or unreadable", name);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      c += r; n -= r;
    }
  }
  template <typename T> void WriteVal(FILE* f, const T& v, const std::string& n) { WriteRaw(f, &v, sizeof(T), n); }
  template <typename T> T ReadVal(FILE* f, const std::string& n) { T v; ReadRaw(f, &v, sizeof(T), n); return v; }
  void WriteStr(FILE* f, const std::string& s, const std::string& n) {
    WriteVal<uint64_t>(f, s.size(), n); WriteRaw(f, s.data(), s.size(), n);
  }
  std::string ReadStr(FILE* f, const std::string& n) {
    const uint64_t l = ReadVal<uint64_t>(f, n);
    if (l > (1u << 20)) { MACH3LOG_ERROR("Corrupt spline cache {}", n); throw MaCh3Exception(__FILE__, __LINE__); }
    std::string s(l, '\0'); ReadRaw(f, &s[0], l, n); return s;
  }
  template <typename T> void WriteVec(FILE* f, const std::vector<T>& v, const std::string& n) {
    WriteVal<uint64_t>(f, v.size(), n); if (!v.empty()) WriteRaw(f, v.data(), v.size()*sizeof(T), n);
  }
  template <typename T> void ReadVec(FILE* f, std::vector<T>& v, const std::string& n) {
    const uint64_t l = ReadVal<uint64_t>(f, n);
    v.clear(); v.shrink_to_fit(); v.resize(l);
    if (l) ReadRaw(f, v.data(), l*sizeof(T), n);
  }

  // Header of one spline parameter in systsHeader (same reading as MonolithSplineHandlerDUNE)
  struct HeaderEntry {
    std::vector<double> knotsFile;  // as stored
    std::vector<size_t> order;      // order[k] = index in file array of sorted knot k
    std::vector<double> knots;      // sorted
  };

  std::map<std::string, HeaderEntry> ReadHeader(TFile* file) {
    std::map<std::string, HeaderEntry> out;
    TTreeReader reader("systsHeader", file);
    TTreeReaderValue<std::string> names(reader, "name");
    TTreeReaderArray<double> knots(reader, "variations");
    while (reader.Next()) {
      const std::string n = *names;
      if (n.empty()) {
        MACH3LOG_ERROR("Empty spline name in {}", file->GetName());
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      if (out.count(n)) continue;  // keep the first occurrence, as the original loader does
      HeaderEntry h;
      for (auto k : knots) h.knotsFile.push_back(k);
      h.order.resize(h.knotsFile.size());
      std::iota(h.order.begin(), h.order.end(), 0);
      std::sort(h.order.begin(), h.order.end(), [&](size_t a, size_t b) { return h.knotsFile[a] < h.knotsFile[b]; });
      for (auto i : h.order) h.knots.push_back(h.knotsFile[i]);
      out[n] = std::move(h);
    }
    return out;
  }
} // namespace

// *****************************************
void FastAtmSplineHandler::Block::Clear() {
  std::vector<uint64_t>().swap(SplOff);
  std::vector<uint64_t>().swap(DatOff);
  std::vector<uint8_t>().swap(Par);
  std::vector<float>().swap(Dat);
  std::vector<double>().swap(DatD);
}
size_t FastAtmSplineHandler::Block::Bytes() const {
  return SplOff.capacity()*sizeof(uint64_t) + DatOff.capacity()*sizeof(uint64_t) + Par.capacity() + Dat.capacity()*sizeof(float) + DatD.capacity()*sizeof(double);
}

// *****************************************
FastAtmSplineHandler::FastAtmSplineHandler(ParameterHandlerGeneric* pars, const std::vector<SplineParameter>& splinePars,
                                           const std::vector<unsigned int>& eventIDs, const std::vector<std::string>& rootFiles,
                                           const FastAtmSplineOptions& opts) : SplineBase(), fPars(pars), fOpts(opts) {
// *****************************************
  TStopwatch clock; clock.Start();
  if (splinePars.size() > 255) {
    MACH3LOG_ERROR("FastAtmSplineHandler supports at most 255 spline parameters per sample, got {}", splinePars.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  nParams = static_cast<short int>(splinePars.size());
  const int P = nParams;
  fNames.resize(P); fGlobalIndex.resize(P); fValPtr.resize(P); fIsFrozen.assign(P, 0);
  for (int p = 0; p < P; ++p) {
    fNames[p] = splinePars[p].name;
    fGlobalIndex[p] = splinePars[p].index;
    fValPtr[p] = fPars->RetPointer(splinePars[p].index);
    if (fOpts.FixedMode > 0 && fPars->IsParameterFixed(splinePars[p].index)) fIsFrozen[p] = 1;
  }

  // ---- unique events
  std::vector<uint32_t> uniq(eventIDs.begin(), eventIDs.end());
  std::sort(uniq.begin(), uniq.end());
  uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
  fUid.resize(eventIDs.size());
  for (size_t i = 0; i < eventIDs.size(); ++i) {
    fUid[i] = static_cast<uint32_t>(std::lower_bound(uniq.begin(), uniq.end(), eventIDs[i]) - uniq.begin());
  }
  MACH3LOG_INFO("FastAtmSplineHandler: {} events, {} unique eventIDs ({:.3f} events per ID), {} spline parameters",
                eventIDs.size(), uniq.size(), uniq.empty() ? 0. : double(eventIDs.size())/double(uniq.size()), P);

  // ---- get the knot values: cache or ROOT
  RawSplines raw;
  bool fromCache = false;
  if (!fOpts.CacheFile.empty() && FileExists(fOpts.CacheFile)) {
    TStopwatch cclock; cclock.Start();
    ReadCache(fOpts.CacheFile, raw, false);
    cclock.Stop();
    bool ok = true;
    if (!fOpts.SourceMD5.empty() && !raw.SourceMD5.empty() && raw.SourceMD5 != fOpts.SourceMD5) {
      MACH3LOG_ERROR("Spline cache {} was made from a spline file with MD5 {}, but the config expects InputSplinesMD5Sum {}",
                     fOpts.CacheFile, raw.SourceMD5, fOpts.SourceMD5);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    for (int p = 0; p < P && ok; ++p) {
      if (std::find(raw.ParNames.begin(), raw.ParNames.end(), fNames[p]) == raw.ParNames.end()) {
        MACH3LOG_WARN("Spline cache {} does not contain parameter {}", fOpts.CacheFile, fNames[p]);
        ok = false;
      }
    }
    if (ok) {
      // every needed eventID must be in the cache (both lists sorted)
      size_t c = 0;
      for (size_t u = 0; u < uniq.size() && ok; ++u) {
        while (c < raw.EventIDs.size() && raw.EventIDs[c] < uniq[u]) ++c;
        if (c == raw.EventIDs.size() || raw.EventIDs[c] != uniq[u]) {
          MACH3LOG_WARN("Spline cache {} does not contain eventID {}", fOpts.CacheFile, uniq[u]);
          ok = false;
        }
      }
    }
    if (ok) {
      fromCache = true;
      MACH3LOG_INFO("Loaded spline cache {} in {:.1f} s ({} events, {} splines, {} parameters, source MD5 {})",
                    fOpts.CacheFile, cclock.RealTime(), raw.EventIDs.size(), raw.Par.size(), raw.ParNames.size(), raw.SourceMD5);
    } else if (rootFiles.empty()) {
      MACH3LOG_ERROR("Spline cache {} does not cover this configuration and no ROOT spline file is available", fOpts.CacheFile);
      throw MaCh3Exception(__FILE__, __LINE__);
    } else {
      MACH3LOG_WARN("Spline cache does not cover this configuration: reading the ROOT spline file instead (cache left unchanged)");
      raw = RawSplines();
    }
  }

  if (!fromCache) {
    if (rootFiles.empty()) {
      MACH3LOG_ERROR("No spline input: neither a usable SplineCacheFile nor InputSplines given");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    if (fOpts.CheckSourceMD5 && !fOpts.SourceMD5.empty() && rootFiles.size() == 1) {
      if (!CompareMD5Sum(fOpts.SourceMD5, rootFiles[0])) {
        MACH3LOG_ERROR("Input Splines MD5Sum does not match the expectation from the config file");
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
    ReadFromRoot(rootFiles, fNames, uniq, fOpts.NThreads, raw);
    raw.SourceMD5 = fOpts.SourceMD5;
    if (!fOpts.CacheFile.empty() && fOpts.WriteCache && !FileExists(fOpts.CacheFile)) {
      TStopwatch wclock; wclock.Start();
      WriteCache(fOpts.CacheFile, raw);
      wclock.Stop();
      MACH3LOG_INFO("Wrote spline cache {} in {:.1f} s", fOpts.CacheFile, wclock.RealTime());
    }
  }

  // ---- knots per parameter (from the file header / cache)
  fNKnots.resize(P); fKnots.resize(P);
  for (int p = 0; p < P; ++p) {
    const size_t q = std::find(raw.ParNames.begin(), raw.ParNames.end(), fNames[p]) - raw.ParNames.begin();
    fKnots[p] = raw.Knots[q];
    fNKnots[p] = static_cast<int>(fKnots[p].size());
    if (fNKnots[p] < 1 || fNKnots[p] > kMaxKnots) {
      MACH3LOG_ERROR("Spline {} has {} knots; FastAtmSplineHandler supports 1-{}", fNames[p], fNKnots[p], kMaxKnots);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    if (raw.NonFinite.size() > q && raw.NonFinite[q] > 0) {
      MACH3LOG_WARN("  - {}: {} non-finite (NaN/inf) spline values replaced with 1.0", fNames[p], raw.NonFinite[q]);
    }
  }

  BuildBasis();

  // ---- fill the blocks for the unique events
  Fill(raw, uniq);
  { RawSplines empty; std::swap(raw, empty); }

  fWeights.assign(uniq.size(), 1.0);
  fLastVal.assign(P, 0.0);
  fFrozenVal.assign(P, 0.0);
  fW.assign(size_t(P)*kMaxKnots, 0.0);
  fSeg.assign(P, 0); fDx.assign(P, 0.f); fXSeg.assign(P, 0.f); fXSegD.assign(P, 0.0);

  if (fHaveFrozen) {
    fFrozenProd.assign(uniq.size(), 1.0);
    bool a = false, f = false;
    UpdateParameterState(true, a, f);
    EvalBlock(fFrozen, nullptr, fFrozenProd.data());
    if (fOpts.FixedMode == 2) {
      for (int p = 0; p < P; ++p) fFrozenVal[p] = fLastVal[p];
      fFrozen.Clear();
      fFrozenDropped = true;
      MACH3LOG_INFO("Splines of fixed parameters evaluated once and dropped (SplineFixedParams: Drop)");
    }
  }

  clock.Stop();
  PrintSummary();
  MACH3LOG_INFO("FastAtmSplineHandler set up in {:.1f} s", clock.RealTime());
}

// *****************************************
FastAtmSplineHandler::~FastAtmSplineHandler() {
// *****************************************
  MACH3LOG_INFO("FastAtmSplineHandler: {} evaluations ({} skipped: no spline parameter changed, {} frozen-block recalculations); event loop {:.3f} s in total, {:.2f} ms per evaluation",
                fNEval, fNSkip, fNFrozenRecalc, fEvalSeconds, (fNEval > fNSkip) ? 1e3*fEvalSeconds/double(fNEval - fNSkip) : 0.);
}

// *****************************************
void FastAtmSplineHandler::ReadFromRoot(const std::vector<std::string>& rootFiles, const std::vector<std::string>& parNames,
                                        const std::vector<uint32_t>& uniqueIDs, int nThreads, RawSplines& out) {
// *****************************************
  TStopwatch clock; clock.Start();
  // ---- headers (every file must agree on the requested parameters)
  std::vector<HeaderEntry> hdr(parNames.size());
  for (size_t iFile = 0; iFile < rootFiles.size(); ++iFile) {
    std::unique_ptr<TFile> f(TFile::Open(rootFiles[iFile].c_str(), "READ"));
    if (!f || f->IsZombie() || !f->Get("SystWeights") || !f->Get("systsHeader")) {
      MACH3LOG_ERROR("Could not open SystWeights/systsHeader in spline input {}", rootFiles[iFile]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    auto h = ReadHeader(f.get());
    for (size_t p = 0; p < parNames.size(); ++p) {
      auto it = h.find(parNames[p]);
      if (it == h.end()) {
        MACH3LOG_ERROR("Desired spline '{}' not found in the systsHeader of {}", parNames[p], rootFiles[iFile]);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      if (iFile == 0) {
        hdr[p] = it->second;
        if (hdr[p].knots != hdr[p].knotsFile) {
          MACH3LOG_WARN("Spline '{}' lists its knots out of order ({}); sorting them", parNames[p], fmt::join(hdr[p].knotsFile, ", "));
        }
      } else if (it->second.knotsFile != hdr[p].knotsFile) {
        MACH3LOG_ERROR("Conflicting spline header '{}' in {}", parNames[p], rootFiles[iFile]);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
  }

  // ---- entries and branch checks
  Long64_t entries = 0;
  {
    TChain ch("SystWeights");
    for (const auto& f : rootFiles) ch.Add(f.c_str());
    entries = ch.GetEntries();
    if (!uniqueIDs.empty() && static_cast<Long64_t>(uniqueIDs.back()) >= entries) {
      MACH3LOG_ERROR("Requested event index {} exceeds SystWeights chain size {}", uniqueIDs.back(), entries);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    ch.LoadTree(0);
    for (size_t p = 0; p < parNames.size(); ++p) {
      TBranch* br = ch.GetBranch(parNames[p].c_str());
      if (!br) {
        MACH3LOG_ERROR("Desired spline branch '{}' not found in {}", parNames[p], fmt::join(rootFiles, ", "));
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      TLeaf* leaf = static_cast<TLeaf*>(br->GetListOfLeaves()->At(0));
      if (!leaf || leaf->GetLenStatic() != static_cast<Int_t>(hdr[p].knots.size()) || leaf->GetLeafCount() != nullptr ||
          std::string(leaf->GetTypeName()) != "Double_t") {
        MACH3LOG_ERROR("Branch '{}' is not a fixed array of {} doubles", parNames[p], hdr[p].knots.size());
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    }
  }

  out.ParNames = parNames;
  out.Knots.resize(parNames.size());
  for (size_t p = 0; p < parNames.size(); ++p) out.Knots[p] = hdr[p].knots;
  out.EventIDs = uniqueIDs;
  out.NonFinite.assign(parNames.size(), 0);

  // ---- threaded read
#ifdef MULTITHREAD
  if (nThreads <= 0) nThreads = omp_get_max_threads();
#endif
  if (nThreads <= 0) nThreads = 1;
  const size_t U = uniqueIDs.size();
  nThreads = static_cast<int>(std::max<size_t>(1, std::min<size_t>(nThreads, U/20000 + 1)));
  ROOT::EnableThreadSafety();
  MACH3LOG_INFO("Reading {} spline branches for {} unique events from {} SystWeights entries with {} threads...",
                parNames.size(), U, entries, nThreads);

  struct Part { std::vector<uint8_t> Count, Par; std::vector<double> Delta; std::vector<uint64_t> NonFinite; };
  std::vector<Part> parts(nThreads);
  std::vector<std::exception_ptr> errors(nThreads);
  std::atomic<size_t> done(0);
  std::mutex printMutex;
  size_t lastPrint = 0;

  auto worker = [&](int t) {
    try {
      const size_t b = U * t / nThreads, e = U * (t + 1) / nThreads;
      Part& part = parts[t];
      part.Count.resize(e - b);
      part.NonFinite.assign(parNames.size(), 0);
      part.Par.reserve((e - b) * 10);
      part.Delta.reserve((e - b) * 64);
      TChain ch("SystWeights");
      for (const auto& f : rootFiles) ch.Add(f.c_str());
      ch.SetBranchStatus("*", 0);
      std::vector<std::vector<double>> buf(parNames.size());
      for (size_t p = 0; p < parNames.size(); ++p) {
        buf[p].assign(hdr[p].knots.size(), 1.0);
        ch.SetBranchStatus(parNames[p].c_str(), 1);
        ch.SetBranchAddress(parNames[p].c_str(), buf[p].data());
      }
      double vals[kMaxKnots];
      for (size_t i = b; i < e; ++i) {
        if (ch.GetEntry(uniqueIDs[i]) <= 0) {
          MACH3LOG_ERROR("Failed to read SystWeights entry {}", uniqueIDs[i]);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
        unsigned int nspl = 0;
        for (size_t p = 0; p < parNames.size(); ++p) {
          const auto& ord = hdr[p].order;
          const size_t n = ord.size();
          bool flat = true;
          for (size_t k = 0; k < n; ++k) {
            double v = buf[p][ord[k]];
            if (!std::isfinite(v)) { v = 1.0; ++part.NonFinite[p]; }  // NaN and inf -> 1
            vals[k] = v;
            if (std::abs(v - 1.0) > 1e-9) flat = false;               // same test as ReusableSpline::IsFlat
          }
          if (flat) continue;
          part.Par.push_back(static_cast<uint8_t>(p));
          for (size_t k = 0; k < n; ++k) part.Delta.push_back(vals[k] - 1.0);
          ++nspl;
        }
        if (nspl > 255) {
          MACH3LOG_ERROR("Event {} has {} non-flat splines (max 255)", uniqueIDs[i], nspl);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
        part.Count[i - b] = static_cast<uint8_t>(nspl);
        const size_t d = ++done;
        if (d % 50000 == 0 || d == U) {
          std::lock_guard<std::mutex> lk(printMutex);
          if (d * 20 / U > lastPrint || d == U) { lastPrint = d * 20 / U; M3::Utils::PrintProgressBar(d, U); }
        }
      }
    } catch (...) {
      errors[t] = std::current_exception();
    }
  };
  std::vector<std::thread> threads;
  for (int t = 0; t < nThreads; ++t) threads.emplace_back(worker, t);
  for (auto& th : threads) th.join();
  for (auto& e : errors) if (e) std::rethrow_exception(e);

  // ---- concatenate
  size_t nSpl = 0, nVal = 0;
  for (auto& p : parts) { nSpl += p.Par.size(); nVal += p.Delta.size(); }
  out.Count.clear(); out.Count.reserve(U);
  out.Par.clear(); out.Par.reserve(nSpl);
  out.Delta.clear(); out.Delta.reserve(nVal);
  for (auto& p : parts) {
    out.Count.insert(out.Count.end(), p.Count.begin(), p.Count.end());
    out.Par.insert(out.Par.end(), p.Par.begin(), p.Par.end());
    out.Delta.insert(out.Delta.end(), p.Delta.begin(), p.Delta.end());
    for (size_t q = 0; q < parNames.size(); ++q) out.NonFinite[q] += p.NonFinite[q];
    Part empty; std::swap(p, empty);
  }
  out.NDelta = out.Delta.size();
  out.DeltaFile = "";
  clock.Stop();
  MACH3LOG_INFO("Read {} non-flat splines ({:.2f} per unique event, {} knot values) in {:.1f} s",
                nSpl, U ? double(nSpl)/double(U) : 0., nVal, clock.RealTime());
}

// *****************************************
void FastAtmSplineHandler::WriteCache(const std::string& fileName, const RawSplines& in) {
// *****************************************
  const std::string tmp = fileName + ".tmp";
  FILE* f = fopen(tmp.c_str(), "wb");
  if (!f) {
    MACH3LOG_ERROR("Cannot write spline cache {}", tmp);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  WriteRaw(f, kCacheMagic, 8, fileName);
  WriteVal<uint32_t>(f, kCacheVersion, fileName);
  WriteStr(f, in.SourceMD5, fileName);
  WriteVal<uint64_t>(f, in.ParNames.size(), fileName);
  for (size_t p = 0; p < in.ParNames.size(); ++p) {
    WriteStr(f, in.ParNames[p], fileName);
    WriteVec(f, in.Knots[p], fileName);
    WriteVal<uint64_t>(f, in.NonFinite.size() > p ? in.NonFinite[p] : 0, fileName);
  }
  WriteVec(f, in.EventIDs, fileName);
  WriteVec(f, in.Count, fileName);
  WriteVec(f, in.Par, fileName);
  WriteVec(f, in.Delta, fileName);
  WriteRaw(f, kCacheEnd, 8, fileName);
  if (fclose(f) != 0 || std::rename(tmp.c_str(), fileName.c_str()) != 0) {
    MACH3LOG_ERROR("Failed to finalise spline cache {}", fileName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

// *****************************************
bool FastAtmSplineHandler::ReadCache(const std::string& fileName, RawSplines& out, bool loadDelta) {
// *****************************************
  FILE* f = fopen(fileName.c_str(), "rb");
  if (!f) return false;
  char magic[8];
  ReadRaw(f, magic, 8, fileName);
  if (std::memcmp(magic, kCacheMagic, 8) != 0 || ReadVal<uint32_t>(f, fileName) != kCacheVersion) {
    fclose(f);
    MACH3LOG_ERROR("{} is not a FastAtmSplineHandler cache (version {})", fileName, kCacheVersion);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  out.SourceMD5 = ReadStr(f, fileName);
  const uint64_t nPar = ReadVal<uint64_t>(f, fileName);
  if (nPar > 255) { fclose(f); MACH3LOG_ERROR("Corrupt spline cache {}", fileName); throw MaCh3Exception(__FILE__, __LINE__); }
  out.ParNames.resize(nPar); out.Knots.resize(nPar); out.NonFinite.resize(nPar);
  for (uint64_t p = 0; p < nPar; ++p) {
    out.ParNames[p] = ReadStr(f, fileName);
    ReadVec(f, out.Knots[p], fileName);
    out.NonFinite[p] = ReadVal<uint64_t>(f, fileName);
  }
  ReadVec(f, out.EventIDs, fileName);
  ReadVec(f, out.Count, fileName);
  ReadVec(f, out.Par, fileName);
  if (loadDelta) {
    ReadVec(f, out.Delta, fileName);
    out.NDelta = out.Delta.size();
    out.DeltaFile = "";
  } else {
    out.NDelta = ReadVal<uint64_t>(f, fileName);
    out.DeltaFile = fileName;
    out.DeltaOffset = static_cast<uint64_t>(ftello(f));
    std::vector<double>().swap(out.Delta);
    if (fseeko(f, static_cast<off_t>(out.NDelta*sizeof(double)), SEEK_CUR) != 0) {
      fclose(f);
      MACH3LOG_ERROR("Spline cache {} is truncated", fileName);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
  ReadRaw(f, magic, 8, fileName);
  fclose(f);
  // consistency
  bool ok = std::memcmp(magic, kCacheEnd, 8) == 0 && out.Count.size() == out.EventIDs.size();
  uint64_t nSpl = 0, nVal = 0;
  for (auto c : out.Count) nSpl += c;
  ok = ok && nSpl == out.Par.size();
  for (size_t s = 0; ok && s < out.Par.size(); ++s) {
    if (out.Par[s] >= nPar) { ok = false; break; }
    nVal += out.Knots[out.Par[s]].size();
  }
  ok = ok && nVal == out.NDelta && std::is_sorted(out.EventIDs.begin(), out.EventIDs.end());
  if (!ok) {
    MACH3LOG_ERROR("Spline cache {} is inconsistent (corrupt or truncated); delete it to rebuild", fileName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  return true;
}

// *****************************************
// For every dial: W_j(v) = delta_sj + dx*(B[s][j] + dx*(C[s][j] + dx*D[s][j])) with B,C,D the TSpline3 coefficients of
// segment s for knot values e_j (unit vectors), built by the same ReusableSpline/TSpline3::BuildCoeff as the original.
void FastAtmSplineHandler::BuildBasis() {
// *****************************************
  const int P = nParams;
  fIsBasis.assign(P, 1);
  fB.assign(P, {}); fC.assign(P, {}); fD.assign(P, {});
  for (int p = 0; p < P; ++p) {
    const int n = fNKnots[p];
    const auto& x = fKnots[p];
    fB[p].assign(size_t(kMaxKnots)*kMaxKnots, 0.0);
    fC[p].assign(size_t(kMaxKnots)*kMaxKnots, 0.0);
    fD[p].assign(size_t(kMaxKnots)*kMaxKnots, 0.0);
    if (n == 1) continue;  // constant spline: W_0 = 1
    bool good = true;
    for (int k = 1; k < n; ++k) if (!(x[k] > x[k-1])) good = false;   // repeated knots
    if (good) {
      ReusableSpline rs(x);
      std::vector<double> y(n, 0.0);
      for (int j = 0; j < n && good; ++j) {
        std::fill(y.begin(), y.end(), 0.0); y[j] = 1.0;
        rs.SetVariationWeights(y);
        rs.Refresh();
        for (int s = 0; s < n - 1; ++s) {
          double xx, yy, b, c, d;
          rs.GetCoeff(s, xx, yy, b, c, d);
          if (!std::isfinite(b) || !std::isfinite(c) || !std::isfinite(d) || yy != (s == j ? 1.0 : 0.0)) { good = false; break; }
          fB[p][s*kMaxKnots + j] = b; fC[p][s*kMaxKnots + j] = c; fD[p][s*kMaxKnots + j] = d;
        }
      }
      // self-test against a direct TSpline3 on generic knot values
      if (good) {
        for (int j = 0; j < n; ++j) y[j] = 1.0 + 0.3*std::sin(1.7*j + 0.4) + 0.05*j;
        rs.SetVariationWeights(y); rs.Refresh();
        double maxdiff = 0;
        for (int s = 0; s < n - 1; ++s) {
          double xx, yy, b, c, d;
          rs.GetCoeff(s, xx, yy, b, c, d);
          double bb = 0, cc = 0, dd = 0;
          for (int j = 0; j < n; ++j) { bb += fB[p][s*kMaxKnots+j]*y[j]; cc += fC[p][s*kMaxKnots+j]*y[j]; dd += fD[p][s*kMaxKnots+j]*y[j]; }
          maxdiff = std::max({maxdiff, std::abs(bb-b)/(1+std::abs(b)), std::abs(cc-c)/(1+std::abs(c)), std::abs(dd-d)/(1+std::abs(d))});
        }
        if (!(maxdiff < 1e-9)) {
          MACH3LOG_WARN("Spline {}: knot basis self-test failed (max rel. diff {}), using per-spline coefficients", fNames[p], maxdiff);
          good = false;
        }
      }
    }
    if (!good) {
      fIsBasis[p] = 0;
      MACH3LOG_WARN("Spline {} has degenerate knots ({}); stored as per-spline (y,b,c,d) like the original monolith",
                    fNames[p], fmt::join(x, ", "));
    }
  }
}


namespace {
  // Sequential access to the knot values, from memory or streamed from the cache file
  class DeltaSource {
   public:
    explicit DeltaSource(const FastAtmSplineHandler::RawSplines& raw) : fRaw(raw) {
      if (!raw.DeltaFile.empty()) {
        fFile = fopen(raw.DeltaFile.c_str(), "rb");
        if (!fFile || fseeko(fFile, static_cast<off_t>(raw.DeltaOffset), SEEK_SET) != 0) {
          MACH3LOG_ERROR("Cannot reopen spline cache {}", raw.DeltaFile);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
        fBuf.resize(size_t(1) << 24);  // 16M values = 128 MB window
      }
    }
    ~DeltaSource() { if (fFile) fclose(fFile); }
    /// pointer to values [pos, pos+n); pos must not decrease between calls
    const double* Get(const uint64_t pos, const size_t n) {
      if (!fFile) return fRaw.Delta.data() + pos;
      if (pos < fStart || pos + n > fEnd) {
        if (pos > fEnd) {  // skip forward
          if (fseeko(fFile, static_cast<off_t>((pos - fEnd)*sizeof(double)), SEEK_CUR) != 0) Fail();
          fEnd = pos;
        } else if (pos < fStart) {
          Fail();
        }
        // keep [pos, fEnd) and read more after it
        const size_t keep = static_cast<size_t>(fEnd - pos);
        if (keep) std::memmove(fBuf.data(), fBuf.data() + (pos - fStart), keep*sizeof(double));
        const size_t want = std::min<uint64_t>(fBuf.size() - keep, fRaw.NDelta - fEnd);
        if (want && fread(fBuf.data() + keep, sizeof(double), want, fFile) != want) Fail();
        fStart = pos; fEnd += want;
        if (pos + n > fEnd) Fail();
      }
      return fBuf.data() + (pos - fStart);
    }
   private:
    [[noreturn]] void Fail() {
      MACH3LOG_ERROR("Error streaming knot values from spline cache {}", fRaw.DeltaFile);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    const FastAtmSplineHandler::RawSplines& fRaw;
    FILE* fFile = nullptr;
    std::vector<double> fBuf;
    uint64_t fStart = 0, fEnd = 0;
  };
}

// *****************************************
// Build the active / frozen blocks for the needed unique events (uniq: sorted, as used for fUid)
void FastAtmSplineHandler::Fill(const RawSplines& raw, const std::vector<uint32_t>& uniq) {
// *****************************************
  const int P = nParams;
  std::vector<int> map(raw.ParNames.size(), -1);
  for (size_t q = 0; q < raw.ParNames.size(); ++q) {
    for (int p = 0; p < P; ++p) if (fNames[p] == raw.ParNames[q]) map[q] = p;
    if (map[q] >= 0 && raw.Knots[q] != fKnots[map[q]]) {
      MACH3LOG_ERROR("Knot mismatch for {}", raw.ParNames[q]);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
  for (int p = 0; p < P; ++p) if (fIsFrozen[p]) fHaveFrozen = true;

  const size_t U = uniq.size();
  Block* blocks[2] = {&fActive, &fFrozen};
  for (auto* b : blocks) { b->Clear(); b->SplOff.assign(U + 1, 0); b->DatOff.assign(U + 1, 0); }
  // pass 1: sizes; pass 2: fill (knot values read sequentially, streamed from the cache file if not in memory)
  std::unique_ptr<DeltaSource> src;
  for (int pass = 0; pass < 2; ++pass) {
    if (pass == 1) src.reset(new DeltaSource(raw));
    size_t c = 0;            // raw event index
    uint64_t rs = 0, rd = 0; // raw spline / value offsets of event c
    uint64_t nSpl[2] = {0, 0}, nDat[2] = {0, 0};
    for (size_t u = 0; u < U; ++u) {
      while (c < raw.EventIDs.size() && raw.EventIDs[c] < uniq[u]) {
        for (unsigned k = 0; k < raw.Count[c]; ++k) { rd += raw.Knots[raw.Par[rs]].size(); ++rs; }
        ++c;
      }
      if (c == raw.EventIDs.size() || raw.EventIDs[c] != uniq[u]) {
        MACH3LOG_ERROR("eventID {} missing from spline data", uniq[u]);
        throw MaCh3Exception(__FILE__, __LINE__);
      }
      if (pass == 0) {
        fActive.SplOff[u] = nSpl[0]; fActive.DatOff[u] = nDat[0];
        fFrozen.SplOff[u] = nSpl[1]; fFrozen.DatOff[u] = nDat[1];
      }
      for (unsigned k = 0; k < raw.Count[c]; ++k) {
        const int q = raw.Par[rs];
        const size_t n = raw.Knots[q].size();
        const int p = map[q];
        if (p >= 0) {
          const int ib = fIsFrozen[p] ? 1 : 0;
          const size_t len = fIsBasis[p] ? n : 4*n;
          if (pass == 1) {
            Block& b = *blocks[ib];
            b.Par[nSpl[ib]] = static_cast<uint8_t>(p);
            const bool dbl = fOpts.DoublePrecision;
            float* dst = dbl ? nullptr : &b.Dat[nDat[ib]];
            double* dstD = dbl ? &b.DatD[nDat[ib]] : nullptr;
            const double* delta = src->Get(rd, n);
            if (fIsBasis[p]) {
              if (dbl) std::memcpy(dstD, delta, n*sizeof(double));
              else for (size_t j = 0; j < n; ++j) dst[j] = static_cast<float>(delta[j]);
            } else {
              // per-spline TSpline3 coefficients, cleaned as TSpline3_redDUNE does
              ReusableSpline rsp(fKnots[p]);
              std::vector<double> y(n);
              for (size_t j = 0; j < n; ++j) y[j] = 1.0 + delta[j];
              rsp.SetVariationWeights(y); rsp.Refresh();
              for (size_t j = 0; j < n; ++j) {
                double xx, yy, bb, cc, dd;
                rsp.GetCoeff(static_cast<int>(j), xx, yy, bb, cc, dd);
                if (!std::isfinite(yy)) yy = 1.0;
                if (!std::isfinite(bb)) bb = 0.0;
                if (!std::isfinite(cc)) cc = 0.0;
                if (!std::isfinite(dd)) dd = 0.0;
                if (dbl) { dstD[4*j] = yy; dstD[4*j+1] = bb; dstD[4*j+2] = cc; dstD[4*j+3] = dd; }
                else { dst[4*j] = float(yy); dst[4*j+1] = float(bb); dst[4*j+2] = float(cc); dst[4*j+3] = float(dd); }
              }
            }
          }
          ++nSpl[ib]; nDat[ib] += len;
        }
        rd += n; ++rs;
      }
      ++c;
    }
    if (pass == 0) {
      fActive.SplOff[U] = nSpl[0]; fActive.DatOff[U] = nDat[0];
      fFrozen.SplOff[U] = nSpl[1]; fFrozen.DatOff[U] = nDat[1];
      fActive.Par.resize(nSpl[0]); fFrozen.Par.resize(nSpl[1]);
      if (fOpts.DoublePrecision) { fActive.DatD.resize(nDat[0]); fFrozen.DatD.resize(nDat[1]); }
      else { fActive.Dat.resize(nDat[0]); fFrozen.Dat.resize(nDat[1]); }
    }
  }
  if (!fHaveFrozen || fFrozen.Par.empty()) {
    if (fHaveFrozen) MACH3LOG_INFO("Fixed spline parameters have no non-flat splines in this sample");
    fHaveFrozen = false;
    fFrozen.Clear();
  }
}

// *****************************************
// Per-parameter state for the current values. Basis dials: segment/dx exactly as core's FindSplineSegment
// (below the first knot: segment 0, above the last: segment n-2, extrapolating the cubic), in double.
void FastAtmSplineHandler::UpdateParameterState(bool force, bool& activeChanged, bool& frozenChanged) {
// *****************************************
  for (int p = 0; p < nParams; ++p) {
    const double v = static_cast<double>(*fValPtr[p]);
    if (!force && v == fLastVal[p]) continue;
    fLastVal[p] = v;
    if (fIsFrozen[p]) frozenChanged = true; else activeChanged = true;
    const int n = fNKnots[p];
    const auto& x = fKnots[p];
    double* W = &fW[size_t(p)*kMaxKnots];
    if (n == 1) { W[0] = 1.0; continue; }
    if (fIsBasis[p]) {
      int s;
      if (v <= x[0]) s = 0;
      else if (v >= x[n-1]) s = n - 2;
      else s = static_cast<int>(std::upper_bound(x.begin(), x.end(), v) - x.begin()) - 1;
      if (s > n - 2) s = n - 2;
      const double dx = v - x[s];
      const double* B = &fB[p][size_t(s)*kMaxKnots];
      const double* C = &fC[p][size_t(s)*kMaxKnots];
      const double* D = &fD[p][size_t(s)*kMaxKnots];
      for (int j = 0; j < n; ++j) W[j] = (j == s ? 1.0 : 0.0) + dx*(B[j] + dx*(C[j] + dx*D[j]));
    } else {
      // fallback dials: same float arithmetic as UnbinnedSplineHandler
      const float xv = static_cast<float>(v);
      int s = 0;
      if (xv <= x[0]) s = 0;
      else if (xv >= x[n-1]) s = n - 1;
      else {
        int kHigh = n - 1;
        while (kHigh - s > 1) { const int kHalf = (s + kHigh)/2; if (xv > x[kHalf]) s = kHalf; else kHigh = kHalf; }
      }
      if (s >= n - 1) s = n - 2;
      fSeg[p] = s;
      fXSeg[p] = static_cast<float>(x[s]);
      fXSegD[p] = x[s];
      fDx[p] = xv - fXSeg[p];
    }
  }
}

// *****************************************
void FastAtmSplineHandler::EvalBlock(const Block& b, const M3::float_t* start, M3::float_t* out) const {
// *****************************************
  if (fOpts.DoublePrecision) EvalBlockT<double>(b, b.DatD.data(), start, out);
  else EvalBlockT<float>(b, b.Dat.data(), start, out);
}

// *****************************************
template <typename T>
void FastAtmSplineHandler::EvalBlockT(const Block& b, const T* dat, const M3::float_t* start, M3::float_t* out) const {
// *****************************************
  const long long U = static_cast<long long>(b.SplOff.size()) - 1;
  if (U <= 0) return;
  const double* W = fW.data();
  const int* nk = fNKnots.data();
  const char* isb = fIsBasis.data();
  const int* seg = fSeg.data();
  const float* dxs = fDx.data();
  const double* lastVal = fLastVal.data();
  const double* xseg = fXSegD.data();
  const uint8_t* par = b.Par.data();
  const uint64_t* so = b.SplOff.data();
  const uint64_t* dof = b.DatOff.data();
  #ifdef MULTITHREAD
  #pragma omp parallel for schedule(static)
  #endif
  for (long long u = 0; u < U; ++u) {
    double w = start ? static_cast<double>(start[u]) : 1.0;
    const T* d = dat + dof[u];
    for (uint64_t s = so[u]; s < so[u+1]; ++s) {
      const int p = par[s];
      const int n = nk[p];
      if (isb[p]) {
        const double* Wp = W + size_t(p)*kMaxKnots;
        double v = 1.0;
        for (int j = 0; j < n; ++j) v += Wp[j]*static_cast<double>(d[j]);
        d += n;
        w *= v;
      } else {
        const T* c = d + 4*seg[p];
        if (sizeof(T) == sizeof(float)) {
          // float store: same float arithmetic as the core monolith
          const float dx = dxs[p];
          w *= static_cast<double>(fmaf(dx, fmaf(dx, fmaf(dx, float(c[3]), float(c[2])), float(c[1])), float(c[0])));
        } else {
          const double dx = lastVal[p] - xseg[p];
          w *= static_cast<double>(c[0]) + dx*(static_cast<double>(c[1]) + dx*(static_cast<double>(c[2]) + dx*static_cast<double>(c[3])));
        }
        d += 4*n;
      }
    }
    out[u] = static_cast<M3::float_t>(w);
  }
}

// *****************************************
void FastAtmSplineHandler::Evaluate() {
// *****************************************
  bool activeChanged = false, frozenChanged = false;
  UpdateParameterState(!fEvaluated, activeChanged, frozenChanged);
  ++fNEval;
  if (frozenChanged && fHaveFrozen) {
    if (fFrozenDropped) {
      bool moved = false;
      for (int p = 0; p < nParams; ++p) {
        if (fIsFrozen[p] && fLastVal[p] != fFrozenVal[p]) {
          MACH3LOG_ERROR("Spline parameter {} was fixed at {} when its splines were dropped (SplineFixedParams: Drop), but now has value {}",
                         fNames[p], fFrozenVal[p], fLastVal[p]);
          moved = true;
        }
      }
      if (moved) {
        MACH3LOG_ERROR("Use SplineFixedParams: Freeze for apps that move fixed parameters");
        throw MaCh3Exception(__FILE__, __LINE__);
      }
    } else {
      EvalBlock(fFrozen, nullptr, fFrozenProd.data());
      ++fNFrozenRecalc;
    }
  }
  if (activeChanged || frozenChanged || !fEvaluated) {
    CalcSplineWeights();
  } else {
    ++fNSkip;
  }
  fEvaluated = true;
}

// *****************************************
void FastAtmSplineHandler::CalcSplineWeights() {
// *****************************************
  const auto t0 = std::chrono::steady_clock::now();
  EvalBlock(fActive, fHaveFrozen ? fFrozenProd.data() : nullptr, fWeights.data());
  fEvalSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// *****************************************
void FastAtmSplineHandler::PrepareSplineFile(std::string FileName) {
// *****************************************
  MACH3LOG_ERROR("FastAtmSplineHandler writes its cache automatically (SplineCacheFile); PrepareSplineFile({}) is not used", FileName);
  throw MaCh3Exception(__FILE__, __LINE__);
}

// *****************************************
void FastAtmSplineHandler::LoadSplineFile(std::string FileName) {
// *****************************************
  MACH3LOG_ERROR("FastAtmSplineHandler reads its cache through SplineCacheFile; LoadSplineFile({}) is not used", FileName);
  throw MaCh3Exception(__FILE__, __LINE__);
}

// *****************************************
size_t FastAtmSplineHandler::GetMemoryBytes() const {
// *****************************************
  return fActive.Bytes() + fFrozen.Bytes() + fUid.capacity()*sizeof(uint32_t) + fWeights.capacity()*sizeof(M3::float_t)
    + fFrozenProd.capacity()*sizeof(M3::float_t);
}

// *****************************************
void FastAtmSplineHandler::PrintSummary() const {
// *****************************************
  const size_t U = fWeights.size();
  MACH3LOG_INFO("--- FastAtmSplineHandler ---");
  MACH3LOG_INFO("{} unique events for {} events; {} active splines ({:.2f}/event, {} values, {} precision), {} frozen splines{}",
                U, fUid.size(), fActive.Par.size(), U ? double(fActive.Par.size())/double(U) : 0., fActive.Dat.size() + fActive.DatD.size(),
                fOpts.DoublePrecision ? "double" : "float",
                fFrozen.Par.size(), fFrozenDropped ? " (dropped)" : "");
  MACH3LOG_INFO("Spline store: {:.1f} MB", double(GetMemoryBytes())/1e6);
  std::vector<uint64_t> cnt(nParams, 0);
  for (auto p : fActive.Par) ++cnt[p];
  for (auto p : fFrozen.Par) ++cnt[p];
  for (int p = 0; p < nParams; ++p) {
    MACH3LOG_INFO("  {:<60} knots {:>2}  splines {:>10}  {}{}", fNames[p], fNKnots[p], cnt[p],
                  fIsBasis[p] ? "basis" : "coeff", fIsFrozen[p] ? (fFrozenDropped ? "  fixed: dropped" : "  fixed: frozen") : "");
  }
}
