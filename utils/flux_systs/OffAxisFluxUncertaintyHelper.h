#pragma once

// forward declaration
class TFile;
class TH1;
class TAxis;

#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// adapted from CAFAna class from:
//  https://github.com/DUNE/lblpwgtools/commit/195f456add17d2b9af67d6c9fc05cb035b843fe0
class OffAxisFluxUncertaintyHelper {
public:
  ~OffAxisFluxUncertaintyHelper() {}

  // nu_config options
  static int const kND_numu_numode = 0;
  static int const kND_nue_numode = 1;
  static int const kND_numubar_numode = 2;
  static int const kND_nuebar_numode = 3;

  static int const kND_SpecHCRun_numu_numode = 4;
  static int const kND_SpecHCRun_nue_numode = 5;
  static int const kND_SpecHCRun_numubar_numode = 6;
  static int const kND_SpecHCRun_nuebar_numode = 7;

  static int const kND_numu_nubarmode = 8;
  static int const kND_nue_nubarmode = 9;
  static int const kND_numubar_nubarmode = 10;
  static int const kND_nuebar_nubarmode = 11;

  static int const kND_SpecHCRun_numu_nubarmode = 12;
  static int const kND_SpecHCRun_nue_nubarmode = 13;
  static int const kND_SpecHCRun_numubar_nubarmode = 14;
  static int const kND_SpecHCRun_nuebar_nubarmode = 15;

  static int const kFD_numu_numode = 16;
  static int const kFD_nue_numode = 17;
  static int const kFD_numubar_numode = 18;
  static int const kFD_nuebar_numode = 19;

  static int const kFD_numu_nubarmode = 20;
  static int const kFD_nue_nubarmode = 21;
  static int const kFD_numubar_nubarmode = 22;
  static int const kFD_nuebar_nubarmode = 23;

  static int const kUnhandled = 24;

  static int const kInvalidBin = std::numeric_limits<int>::max();

  void Initialize(std::string const &filename, bool verbose = false);

  static OffAxisFluxUncertaintyHelper const &Get();
  size_t GetNFocussingParams() const { return focussing.NDuncerts.size(); }
  std::string GetFocussingParamName(size_t i) const {
    return focussing.UncertLabels.at(i);
  }

  size_t GetNHadProdPCAComponents() const { return hadprod.NDuncerts.size(); }

  int GetNuConfig(int nu_pdg, bool IsND, bool IsNuMode,
                  bool isSpecHCRun = false) const;

  std::vector<std::unique_ptr<TH1>>
  GetNDOffAxisShifts(TFile *f, std::string nd_dir, std::string hname) const;

  int GetFocussingBin(int nu_pdg, double enu_GeV, double off_axis_pos_m,
                      int nu_config) const;

  int GetHadProdBin(int nu_pdg, double enu_GeV, double off_axis_pos_m,
                    int nu_config) const;

  double GetFluxFocussingWeight(size_t param_id, double param_val, int bin,
                                int nu_config) const;

  double GetFluxHadProdWeight(size_t param_id, double param_val, int bin,
                              int nu_config) const;

  struct {
    size_t NParams;

    // param
    std::vector<std::unique_ptr<TAxis>> OffAxisTAxes;
    // param | nucfg | OA bin | E bin
    std::vector<std::vector<std::vector<std::unique_ptr<TH1>>>> NDuncerts;
    // param | nucfg | E bin
    std::vector<std::vector<std::unique_ptr<TH1>>> FDuncerts;

    std::vector<std::string> UncertLabels;

  } focussing;

  struct {
    size_t NPCAComponents;

    // param | nucfg
    std::vector<std::vector<std::unique_ptr<TAxis>>> OffAxisTAxes;
    // param | nucfg | OA bin | E bin
    std::vector<std::vector<std::vector<std::unique_ptr<TH1>>>> NDuncerts;
    // param | nucfg | E bin
    std::vector<std::vector<std::unique_ptr<TH1>>> FDuncerts;

  } hadprod;
};
