// FastCompare: validation and timing app.
// Builds the samples from a normal config, reweights at the nominal point and at NTHROWS reproducible parameter
// throws, and writes the total and per-sample -2lnL, integrals and every bin (%.17g) to a text file, plus timings
// and memory. It uses only public MaCh3/MaCh3_DUNE interfaces, so the same source builds against the original
// code (FastOptions: SplineEngine Legacy, FastReweight false) and the fast engine; diff two output files to compare.
//
//   FastCompare Configs/EventRates_Atmospherics.yaml override.yaml
// Environment: FC_OUT (output text file, default FastCompare.txt), FC_NTHROWS (default 20), FC_SEED (default 12345),
//              FC_SCALE (throw width in prior sigma, default 0.5), FC_REPEAT (timed reweights per point, default 1)
//              FC_POINTS: text file of parameter points to evaluate instead of the throws, one point per
//                line with a value for every parameter in config order ('#' lines are comments). For example a chain's
//                accepted steps, to replay that chain in this build.
//              FC_NOBINS=1: do not write the BINS lines (for many points).
// Each point also gets a PENALTY line (ParameterHandler::GetLikelihood, as the fitter's LogL_systematic_*), and the
// header lists the parameters (PARNAME index name fancyname).
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TRandom3.h"

#include "Samples/MaCh3DUNEFactory.h"
#include "Fitters/MaCh3Factory.h"

namespace {
  double ReadStatusKB(const char* key) {
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line)) {
      if (line.rfind(key, 0) == 0) return std::stod(line.substr(std::string(key).size()));
    }
    return -1;
  }
  int EnvInt(const char* n, int d) { const char* e = std::getenv(n); return (e && *e) ? std::atoi(e) : d; }
  double EnvDouble(const char* n, double d) { const char* e = std::getenv(n); return (e && *e) ? std::atof(e) : d; }
}

int main(int argc, char* argv[]) {
  auto t0 = std::chrono::steady_clock::now();
  auto FitManager = MaCh3ManagerFactory(argc, argv);
  auto [param_handler, samples] = MaCh3DuneFactory(FitManager);
  auto t1 = std::chrono::steady_clock::now();
  const double setupSec = std::chrono::duration<double>(t1 - t0).count();
  const double rssSetup = ReadStatusKB("VmRSS:")/1048576.;
  const double hwmSetup = ReadStatusKB("VmHWM:")/1048576.;

  const std::string outName = std::getenv("FC_OUT") ? std::getenv("FC_OUT") : "FastCompare.txt";
  const int nThrows = EnvInt("FC_NTHROWS", 20);
  const int seed = EnvInt("FC_SEED", 12345);
  const double scale = EnvDouble("FC_SCALE", 0.5);
  const int nRepeat = EnvInt("FC_REPEAT", 1);

  FILE* out = std::fopen(outName.c_str(), "w");
  if (!out) { MACH3LOG_ERROR("Cannot write {}", outName); return 1; }
  std::fprintf(out, "# FastCompare  nThrows %d seed %d scale %g\n", nThrows, seed, scale);
  std::fprintf(out, "SETUP_SECONDS %.3f RSS_GB %.3f HWM_GB %.3f\n", setupSec, rssSetup, hwmSetup);

  const int nPar = param_handler->GetNParameters();
  std::vector<double> nominal(nPar);
  for (int i = 0; i < nPar; ++i) nominal[i] = param_handler->GetParProp(i);
  for (int i = 0; i < nPar; ++i)
    std::fprintf(out, "PARNAME %d %s %s\n", i, param_handler->GetParName(i).c_str(), param_handler->GetParFancyName(i).c_str());

  // FC_POINTS: evaluate the given points (e.g. accepted steps of a chain) instead of random throws
  std::vector<std::vector<double>> points;
  const char* pointsFile = std::getenv("FC_POINTS");
  if (pointsFile && *pointsFile) {
    std::ifstream pf(pointsFile);
    if (!pf) { MACH3LOG_ERROR("FC_POINTS: cannot read {}", pointsFile); return 1; }
    std::string line;
    while (std::getline(pf, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream is(line);
      std::vector<double> p; double v;
      while (is >> v) p.push_back(v);
      if (p.empty()) continue;
      if (static_cast<int>(p.size()) != nPar) {
        MACH3LOG_ERROR("FC_POINTS: line {} has {} values, the config has {} parameters", points.size() + 1, p.size(), nPar);
        return 1;
      }
      points.push_back(p);
    }
    MACH3LOG_INFO("FastCompare: {} points from {}", points.size(), pointsFile);
    std::fprintf(out, "# FC_POINTS %s (%zu points replace the throws)\n", pointsFile, points.size());
  }
  const int nPoints = points.empty() ? nThrows : static_cast<int>(points.size());
  const bool writeBins = !(std::getenv("FC_NOBINS") && std::string(std::getenv("FC_NOBINS")) == "1");

  double sumSec = 0; int nTimed = 0;
  for (int t = 0; t <= nPoints; ++t) {
    // reproducible throw: Gaussian around the starting value, width scale*prior sigma, inside the bounds
    TRandom3 rng(static_cast<UInt_t>(seed + 7919*t));
    for (int i = 0; i < nPar; ++i) {
      double v = nominal[i];
      if (t > 0 && !points.empty()) {
        v = points[t - 1][i];
      } else if (t > 0 && !param_handler->IsParameterFixed(i)) {
        const double lo = param_handler->GetLowerBound(i), hi = param_handler->GetUpperBound(i);
        double err = param_handler->GetDiagonalError(i);
        if (!(err > 0) || err > 0.5*(hi - lo)) err = 0.1*(hi - lo);
        v = nominal[i] + scale*err*rng.Gaus();
        const double eps = 1e-6*(hi - lo);
        if (v < lo + eps) v = lo + eps;
        if (v > hi - eps) v = hi - eps;
      }
      param_handler->SetParProp(i, v);
    }
    std::fprintf(out, "PARS %d", t);
    for (int i = 0; i < nPar; ++i) std::fprintf(out, " %.17g", param_handler->GetParProp(i));
    std::fprintf(out, "\n");

    double llhTot = 0;
    for (auto handler : samples) {
      auto r0 = std::chrono::steady_clock::now();
      for (int r = 0; r < nRepeat; ++r) handler->Reweight();
      auto r1 = std::chrono::steady_clock::now();
      const double sec = std::chrono::duration<double>(r1 - r0).count()/nRepeat;
      if (t > 0) { sumSec += sec; ++nTimed; }
      const double llh = handler->GetLikelihood();
      llhTot += llh;
      std::fprintf(out, "REWEIGHT_SECONDS %d %.6f\n", t, sec);
      for (int s = 0; s < handler->GetNSamples(); ++s) {
        const auto mc = handler->GetMCArray(s);
        double integral = 0;
        for (auto v : mc) integral += v;
        std::fprintf(out, "SAMPLE %d %s LLH %.17g INTEGRAL %.17g\n", t, handler->GetSampleTitle(s).c_str(),
                     handler->GetSampleLikelihood(s), integral);
        if (writeBins) {
          std::fprintf(out, "BINS %d %d", t, s);
          for (auto v : mc) std::fprintf(out, " %.17g", v);
          std::fprintf(out, "\n");
        }
      }
    }
    std::fprintf(out, "PENALTY %d %.17g\n", t, param_handler->GetLikelihood());
    std::fprintf(out, "LLH %d %.17g\n", t, llhTot);
    MACH3LOG_INFO("FastCompare point {}: -2lnL(total) {:.10f}", t, llhTot);
  }
  // same point again: shows the cost when nothing changed
  auto s0 = std::chrono::steady_clock::now();
  for (auto handler : samples) handler->Reweight();
  auto s1 = std::chrono::steady_clock::now();
  std::fprintf(out, "SAMEPOINT_SECONDS %.6f\n", std::chrono::duration<double>(s1 - s0).count());
  std::fprintf(out, "MEAN_REWEIGHT_SECONDS %.6f\n", nTimed ? sumSec/nTimed : 0.);
  std::fprintf(out, "FINAL RSS_GB %.3f HWM_GB %.3f\n", ReadStatusKB("VmRSS:")/1048576., ReadStatusKB("VmHWM:")/1048576.);
  std::fclose(out);
  MACH3LOG_INFO("FastCompare: setup {:.1f} s, mean reweight {:.4f} s over {} throws, RSS {:.2f} GB (peak {:.2f} GB); wrote {}",
                setupSec, nTimed ? sumSec/nTimed : 0., nTimed, ReadStatusKB("VmRSS:")/1048576., ReadStatusKB("VmHWM:")/1048576., outName);
  for (auto s : samples) delete s;
  return 0;
}
