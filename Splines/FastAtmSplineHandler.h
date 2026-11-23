#pragma once
// ============================================================================================
// FastAtmSplineHandler: fast evaluation of the atmospheric per-event splines
//
// Drop-in replacement for MonolithSplineHandlerDUNE (event-by-event "PerEvent" splines).
// Same interpolation as the original: for every event and dial, ROOT's TSpline3 (not-a-knot,
// built by the same ReusableSpline code) through the knot weights read from SystWeights, with the
// same knot sorting, the same non-finite -> 1.0 replacement and the same "flat" test
// (all |y-1| <= 1e-9 -> weight 1). What changes is how it is stored and evaluated:
//
//  1. Duplicates removed: every GENIE event appears twice in the atmospheric MC (from the nue and the
//     numu flux) with identical splines. Splines are stored and evaluated once per unique eventID,
//     and both copies point at the same weight (exactly the same number as before).
//  2. Knot-basis storage: a TSpline3 is linear in its knot values, so for a dial with fixed knots
//     x_0..x_{n-1} the spline value at parameter value v is  S(v) = sum_j W_j(v) * y_j , where the
//     W_j(v) depend only on the knots and v (one set per dial per step, from ROOT's own BuildCoeff
//     on unit vectors). Each event spline is stored as its n deviations d_j = y_j - 1 (float), and
//     evaluated as  S = 1 + sum_j W_j d_j  (double). 4x less memory than the (y,b,c,d) per knot of
//     the core monolith, and more precise (rounding is relative to y-1, not to y).
//     Dials whose basis is not finite (repeated knots) fall back to per-spline (y,b,c,d).
//  3. One fused pass per step: the product over an event's splines is formed directly; no per-spline
//     intermediate array.
//  4. Fixed dials: splines of parameters fixed at construction (XsecFix) are kept in a separate
//     block whose product is cached per event and only recomputed if one of them changes value
//     ("Freeze", default). "Drop" evaluates them once and frees them (throws if they ever change).
//  5. Nothing is recomputed if no spline parameter changed since the last call.
//  6. Loading: SystWeights is read with N threads (each with its own TChain), only the needed
//     branches and entries, straight into flat arrays (no per-event heap objects, which is what made
//     the original need ~97 GB). Optionally the result is written to / read from a binary cache
//     file (~2 GB for full57 instead of the 23.5 GB ROOT file).
// ============================================================================================

#include "Splines/SplineBase.h"
#include "Parameters/ParameterHandlerGeneric.h"

#include <cstdint>
#include <string>
#include <vector>

struct FastAtmSplineOptions {
  /// Binary cache file ("" = no cache). Relative paths are resolved against $MACH3.
  std::string CacheFile = "";
  /// Write the cache file when it does not exist yet
  bool WriteCache = true;
  /// MD5 of the ROOT spline file expected by the config (stored in / checked against the cache)
  std::string SourceMD5 = "";
  /// Run the (streaming) MD5 check on the ROOT file before reading it
  bool CheckSourceMD5 = true;
  /// What to do with splines of fixed parameters: 0 = evaluate every step (like the original),
  /// 1 = Freeze (cached product, recomputed only if they change), 2 = Drop (evaluate once, free)
  int FixedMode = 1;
  /// Threads used to read the ROOT file (0 = OpenMP max threads)
  int NThreads = 0;
  /// Store the knot deviations as double (true) or float (false). Double: more precise than the original
  /// everywhere (2x the store, ~3.8 GB for full57); float: half the memory, more precise than the original
  /// within ~1 sigma but up to ~1e-5 relative at +-3 sigma where the cubic extrapolation cancels strongly.
  bool DoublePrecision = true;
};

class FastAtmSplineHandler : public SplineBase {
 public:
  /// @param pars        parameter handler (values are read through RetPointer)
  /// @param splinePars  spline parameters applying to this sample (order defines evaluation order)
  /// @param eventIDs    per-event index into the SystWeights chain (dunemc_atm::eventID), may repeat
  /// @param rootFiles   ROOT files with SystWeights / systsHeader (may be empty if the cache is used)
  FastAtmSplineHandler(ParameterHandlerGeneric* pars, const std::vector<SplineParameter>& splinePars,
                       const std::vector<unsigned int>& eventIDs, const std::vector<std::string>& rootFiles,
                       const FastAtmSplineOptions& opts);
  virtual ~FastAtmSplineHandler();

  void Evaluate() override;
  std::string GetName() const override {return "FastAtmSplineHandler";}
  /// Writes the loaded splines to a cache file
  void PrepareSplineFile(std::string FileName) override;
  void LoadSplineFile(std::string FileName) override;
  void SynchroniseMemTransfer() const override {}

  /// Pointer to the spline weight of event iEvent (shared by events with the same eventID)
  const M3::float_t* RetPointer(const int iEvent) const {return &fWeights[fUid[iEvent]];}
  /// Index of the unique event used by event iEvent
  unsigned int GetUniqueIndex(const int iEvent) const {return fUid[iEvent];}
  size_t GetNUnique() const {return fWeights.size();}
  /// Spline weight of unique event u (after the last Evaluate)
  M3::float_t GetUniqueWeight(const size_t u) const {return fWeights[u];}
  /// Bytes held by the spline store (for the log)
  size_t GetMemoryBytes() const;
  /// Time spent in the event loop so far (seconds) and number of full evaluations
  double GetEvalSeconds() const {return fEvalSeconds;}
  unsigned long long GetNFullEvaluations() const {return fNEval - fNSkip;}

  /// Data as read from the ROOT file or the cache: per unique event, the non-flat splines with
  /// their knot deviations y-1 (float), knots sorted increasingly.
  struct RawSplines {
    std::string SourceMD5;
    std::vector<std::string> ParNames;
    std::vector<std::vector<double>> Knots;
    std::vector<uint64_t> NonFinite;        ///< per parameter: non-finite knot values set to 1
    std::vector<uint32_t> EventIDs;         ///< sorted unique entry indices
    std::vector<uint8_t> Count;             ///< non-flat splines per event
    std::vector<uint8_t> Par;               ///< parameter index (into ParNames) per spline
    std::vector<double> Delta;              ///< y-1 per knot per spline (exact, from the file's doubles)
    std::string DeltaFile;                  ///< if Delta was not loaded: cache file holding it ...
    uint64_t DeltaOffset = 0;               ///< ... at this byte offset ...
    uint64_t NDelta = 0;                    ///< ... with this many values (streamed by Fill)
  };

  /// Read the requested branches for the requested (sorted unique) entries from ROOT files
  static void ReadFromRoot(const std::vector<std::string>& rootFiles, const std::vector<std::string>& parNames,
                           const std::vector<uint32_t>& uniqueIDs, int nThreads, RawSplines& out);
  /// Cache I/O. ReadCache returns false if the file does not exist; throws on a corrupt file
  static void WriteCache(const std::string& fileName, const RawSplines& in);
  /// loadDelta = false leaves the knot values on disk (Fill streams them: no second copy in memory)
  static bool ReadCache(const std::string& fileName, RawSplines& out, bool loadDelta = true);

 protected:
  void CalcSplineWeights() override;

 private:
  /// One storage block (active or frozen splines), CSR over unique events
  struct Block {
    std::vector<uint64_t> SplOff;  ///< U+1 offsets into Par
    std::vector<uint64_t> DatOff;  ///< U+1 offsets into Dat
    std::vector<uint8_t> Par;      ///< local parameter index per spline
    std::vector<float> Dat;        ///< float store: deviations (basis dials) or y,b,c,d per knot (fallback dials)
    std::vector<double> DatD;      ///< double store (DoublePrecision)
    void Clear();
    size_t Bytes() const;
  };

  void BuildBasis();
  void Fill(const RawSplines& raw, const std::vector<uint32_t>& uniq);
  void UpdateParameterState(bool force, bool& activeChanged, bool& frozenChanged);
  void EvalBlock(const Block& b, const M3::float_t* start, M3::float_t* out) const;
  template <typename T> void EvalBlockT(const Block& b, const T* dat, const M3::float_t* start, M3::float_t* out) const;
  void PrintSummary() const;

  static constexpr int kMaxKnots = 16;

  ParameterHandlerGeneric* fPars;
  FastAtmSplineOptions fOpts;

  // Per local parameter
  std::vector<std::string> fNames;
  std::vector<int> fGlobalIndex;
  std::vector<const M3::float_t*> fValPtr;
  std::vector<int> fNKnots;
  std::vector<std::vector<double>> fKnots;
  std::vector<char> fIsBasis;
  std::vector<char> fIsFrozen;
  /// Basis: for segment s, knot j: B,C,D at [s*kMaxKnots + j]
  std::vector<std::vector<double>> fB, fC, fD;
  std::vector<double> fLastVal;
  std::vector<double> fFrozenVal;      ///< values used to evaluate a dropped block
  /// Per-step state used in the event loop
  std::vector<double> fW;              ///< [p*kMaxKnots + j]
  std::vector<int> fSeg;               ///< segment (fallback dials)
  std::vector<float> fDx;              ///< dx (fallback dials, float as in core)
  std::vector<float> fXSeg;            ///< x at segment start (fallback dials, float)
  std::vector<double> fXSegD;          ///< x at segment start (fallback dials, double store)

  std::vector<uint32_t> fUid;          ///< per event -> unique index
  std::vector<M3::float_t> fWeights;   ///< per unique event: spline weight
  std::vector<double> fFrozenProd;     ///< per unique event: product of frozen splines

  Block fActive, fFrozen;
  bool fHaveFrozen = false;
  bool fFrozenDropped = false;
  bool fEvaluated = false;
  unsigned long long fNEval = 0, fNSkip = 0, fNFrozenRecalc = 0;
  double fEvalSeconds = 0;
};
