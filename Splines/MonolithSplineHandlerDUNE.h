#pragma once
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>
#include "TFile.h"  // Used only for per-file systsHeader metadata validation
#include "Splines/UnbinnedSplineHandler.h"

struct SplineHeader {
    std::string name;
    std::vector<double> knots;
    bool isCorrection;
    double cv = 0.0;
};

class MonolithSplineHandlerDUNE : virtual public UnbinnedSplineHandler {
public:
  MonolithSplineHandlerDUNE(const std::vector<SplineParameter>& splinePars, const std::vector<uint>& eventIndices, const std::vector<std::string>& spline_filename);
  /// @brief Construct from per-event splines built by the caller (e.g. straight from the CAF xsec_systs)
  /// @param splinesPerEvent [event][parameter] splines, nullptr for a flat response. Ownership is taken (and the objects deleted) by UnbinnedSplineHandler
  /// @param splineTypes Response function type of each parameter
  MonolithSplineHandlerDUNE(std::vector<std::vector<TResponseFunction_red*> >& splinesPerEvent, const std::vector<RespFuncType>& splineTypes);
  virtual ~MonolithSplineHandlerDUNE();
  
private:
  void InitFromFile(std::string &spline_filename);
  MonolithSplineHandlerDUNE(std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>> initParams);
  static std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>> GetInitParamsFromConfig(const std::vector<SplineParameter>& splinePars,
															 const std::vector<uint>& eventIndices,
															 const std::vector<std::string>& spline_filename);
  static std::vector<struct SplineHeader> GetSplineParametersFromFile(TFile *f, const std::string& treeName);
};

//Forced to write a new version as the M3 core one deleted the splines it is passed................
class TSpline3_redDUNE : public TSpline3_red {
public:
  TSpline3_redDUNE(TSpline3* spline) : TSpline3_red() {
    nPoints = spline->GetNp();
    Par = new M3::float_t*[nPoints];
    XPos = new M3::float_t[nPoints];
    YResp = new M3::float_t[nPoints];
    for (int i = 0; i < nPoints; ++i) {
      // 3 is the size of the TSpline3 coefficients
      Par[i] = new M3::float_t[3];
      double x = -999.99, y = -999.99, b = -999.99, c = -999.99, d = -999.99;
      spline->GetCoeff(i, x, y, b, c, d);
      // Clean up any NaNs or infinities to prevent propagation to event weights
      if (std::isnan(x) || std::isinf(x)) x = 0.0;
      if (std::isnan(y) || std::isinf(y)) y = 1.0;
      if (std::isnan(b) || std::isinf(b)) b = 0.0;
      if (std::isnan(c) || std::isinf(c)) c = 0.0;
      if (std::isnan(d) || std::isinf(d)) d = 0.0;
      XPos[i]   = x;
      YResp[i]  = y;
      Par[i][0] = b;
      Par[i][1] = c;
      Par[i][2] = d;
    }
  }
  
  //Adding a copy constructor because we cannot reuse the same object as M3 core once again loves to delete objects it's given...
  TSpline3_redDUNE(const TSpline3_redDUNE& other) : TSpline3_red() {
    nPoints = other.nPoints;
    Par = new M3::float_t*[nPoints];
    XPos = new M3::float_t[nPoints];
    YResp = new M3::float_t[nPoints];
    for (int i = 0; i < nPoints; ++i) {
      Par[i] = new M3::float_t[3];
      XPos[i] = other.XPos[i];
      YResp[i] = other.YResp[i];
      Par[i][0] = other.Par[i][0];
      Par[i][1] = other.Par[i][1];
      Par[i][2] = other.Par[i][2];
      // Double check values
      if (std::isnan(XPos[i]) || std::isinf(XPos[i])) XPos[i] = 0.0;
      if (std::isnan(YResp[i]) || std::isinf(YResp[i])) YResp[i] = 1.0;
      if (std::isnan(Par[i][0]) || std::isinf(Par[i][0])) Par[i][0] = 0.0;
      if (std::isnan(Par[i][1]) || std::isinf(Par[i][1])) Par[i][1] = 0.0;
      if (std::isnan(Par[i][2]) || std::isinf(Par[i][2])) Par[i][2] = 0.0;
    }
  }
  
  bool IsFlat() const {
    for (int i = 0; i < nPoints; ++i) {
      if (std::abs(YResp[i] - 1.0) > 1e-9) {
	return false;
      }
    }
    return true;
  }
  
};

class ReusableSpline : public TSpline3 {
public:
  // Inherit standard constructors
  ReusableSpline(const std::vector<double>& knots) : TSpline3(){
    fValBeg = 0;
    fValEnd = 0;
    fBegCond = 0;
    fEndCond = 0;
    fNp = static_cast<Int_t>(knots.size());
    fPoly = new TSplinePoly3[fNp];
    for (int i = 0; i < fNp; ++i) {
      fPoly[i].X() = knots[i];
      fPoly[i].Y() = 1.0;  // Initialize Y to 1.0 (neutral weight)
    }
    fXmin = fPoly[0].X();
    fXmax = fPoly[fNp-1].X();
    fName = "ReusableSpline";
  }
  
  bool IsFlat() const {
    for (int i = 0; i < fNp; ++i) {
      if (std::abs(fPoly[i].Y() - 1.0) > 1e-9) {
	return false;
      }
    }
    return true;
  }

  // Expose the protected BuildCoeff method
  void Refresh() {
    this->BuildCoeff();
  }
  
  void SetVariationWeights(const std::vector<double>& weights) {
    if(static_cast<Int_t>(weights.size()) != fNp) {
      MACH3LOG_ERROR("Size of weights does not match number of nodes in spline.");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    for (int i = 0; i < fNp; ++i) {
      fPoly[i].Y() = weights[i];
    }
  }
  
  void PrintPoly() {
    for (int i = 0; i < fNp; ++i) {
      double x = 0, y = 0, b = 0, c = 0, d = 0;
      this->GetCoeff(i, x, y, b, c, d);
      MACH3LOG_INFO("Knot {}: X = {:.2f}, Y = {:.2f}, b = {:.2f}, c = {:.2f}, d = {:.2f}", i, x, y, b, c, d);
    }
  }

};

/// @brief PG Chunked bump allocator for splines built in bulk (hundreds of millions of small objects for a full
///        atmospheric sample): one allocation per chunk instead of several per spline, all released at once.
class SplineArena {
public:
  explicit SplineArena(std::size_t ChunkSize = 64 << 20) : fChunkSize(ChunkSize) {}

  void* Allocate(std::size_t Size, std::size_t Alignment) {
    std::size_t Offset = (fUsed + Alignment - 1) & ~(Alignment - 1);
    if (fChunks.empty() || Offset + Size > fCurrentChunkSize) {
      fCurrentChunkSize = std::max(fChunkSize, Size);
      fChunks.emplace_back(new char[fCurrentChunkSize]);
      fAllocatedBytes += fCurrentChunkSize;
      Offset = 0;
    }
    fUsed = Offset + Size;
    return fChunks.back().get() + Offset;
  }

  /// @brief Release all the memory (every object allocated in the arena must be unused by then)
  void Clear() {
    fChunks.clear();
    fChunks.shrink_to_fit();
    fUsed = fCurrentChunkSize = fAllocatedBytes = 0;
  }

  std::size_t GetAllocatedBytes() const { return fAllocatedBytes; }

private:
  std::size_t fChunkSize;
  std::size_t fCurrentChunkSize = 0;
  std::size_t fUsed = 0;
  std::size_t fAllocatedBytes = 0;
  std::vector<std::unique_ptr<char[]>> fChunks;
};

/// @brief PG TSpline3_red living in a SplineArena, whose knots can be shared by several events (e.g. the nue- and
///        numu-flux copies of a CAF entry). UnbinnedSplineHandler deletes every spline it is given: here that only
///        runs the destructor, the memory goes away with the arena.
class TSpline3_redShared : public TSpline3_red {
public:
  /// @brief Copy the coefficients of a TSpline3 into the arena
  static TSpline3_redShared* Create(SplineArena& Arena, TSpline3* Spline) {
    const int n = Spline->GetNp();
    auto** ParArray = static_cast<M3::float_t**>(Arena.Allocate(n * sizeof(M3::float_t*), alignof(M3::float_t*)));
    //PG x and y of each knot, followed by the b, c, d coefficients of each knot
    auto* Data = static_cast<M3::float_t*>(Arena.Allocate(5 * n * sizeof(M3::float_t), alignof(M3::float_t)));
    for (int i = 0; i < n; ++i) {
      double x = -999.99, y = -999.99, b = -999.99, c = -999.99, d = -999.99;
      Spline->GetCoeff(i, x, y, b, c, d);
      // Clean up any NaNs or infinities to prevent propagation to event weights
      if (!std::isfinite(x)) x = 0.0;
      if (!std::isfinite(y)) y = 1.0;
      if (!std::isfinite(b)) b = 0.0;
      if (!std::isfinite(c)) c = 0.0;
      if (!std::isfinite(d)) d = 0.0;
      Data[i] = x;
      Data[n + i] = y;
      ParArray[i] = Data + 2 * n + 3 * i;
      ParArray[i][0] = b;
      ParArray[i][1] = c;
      ParArray[i][2] = d;
    }
    return new (Arena.Allocate(sizeof(TSpline3_redShared), alignof(TSpline3_redShared))) TSpline3_redShared(n, ParArray, Data, Data + n);
  }

  /// @brief Another spline object using the same knots
  static TSpline3_redShared* Share(SplineArena& Arena, const TSpline3_redShared& Other) {
    return new (Arena.Allocate(sizeof(TSpline3_redShared), alignof(TSpline3_redShared))) TSpline3_redShared(Other.nPoints, Other.Par, Other.XPos, Other.YResp);
  }

  ~TSpline3_redShared() override {
    // The knots belong to the arena: leave nothing for ~TSpline3_red to free
    Par = nullptr;
    XPos = nullptr;
    YResp = nullptr;
    nPoints = 0;
  }

  // The object itself lives in the arena too: it can only be placed there, and deleting it releases nothing
  static void* operator new(std::size_t) = delete;
  static void* operator new(std::size_t, void* Where) noexcept { return Where; }
  static void operator delete(void*) noexcept {}
  static void operator delete(void*, void*) noexcept {}

private:
  TSpline3_redShared(M3::int_t n, M3::float_t** ParArray, M3::float_t* X, M3::float_t* Y) : TSpline3_red() {
    nPoints = n;
    Par = ParArray;
    XPos = X;
    YResp = Y;
  }
};
