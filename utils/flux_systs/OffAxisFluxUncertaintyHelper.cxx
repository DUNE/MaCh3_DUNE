#include "utils/flux_systs/OffAxisFluxUncertaintyHelper.h"

#include "TAxis.h"
#include "TFile.h"
#include "TH1.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

static OffAxisFluxUncertaintyHelper *globalFluxHelper = nullptr;

OffAxisFluxUncertaintyHelper const &OffAxisFluxUncertaintyHelper::Get() {
  if (!globalFluxHelper) {
    globalFluxHelper = new OffAxisFluxUncertaintyHelper();
    globalFluxHelper->Initialize(
        std::string(std::getenv("MACH3")) +
        "/utils/flux_systs/flux_variations_FD_and_PRISM_2023.root");
  }
  return *globalFluxHelper;
}

bool check_axis(TAxis *ref, TAxis *other) {
  if (ref->GetNbins() != other->GetNbins()) {
    std::cout << "-- Mismatched TAxis: Reference NBins = " << ref->GetNbins()
              << ", Other NBins = " << other->GetNbins() << std::endl;
    return false;
  }

  for (int bi = 0; bi < ref->GetNbins(); ++bi) {
    if (std::fabs(ref->GetBinLowEdge(bi + 1) - other->GetBinLowEdge(bi + 1)) >
        1E-8) {
      std::cout << "-- Mismatched TAxis: Low edge of bin " << bi
                << ", Reference low edge: " << ref->GetBinLowEdge(bi + 1)
                << ", other low edge: " << other->GetBinLowEdge(bi + 1)
                << std::endl;
      return false;
    }
    if (std::fabs(ref->GetBinUpEdge(bi + 1) - other->GetBinUpEdge(bi + 1)) >
        1E-8) {
      std::cout << "-- Mismatched TAxis: Up edge of bin " << bi
                << ", Reference up edge: " << ref->GetBinUpEdge(bi + 1)
                << ", other up edge: " << other->GetBinUpEdge(bi + 1)
                << std::endl;
      return false;
    }
  }
  return true;
}

void OffAxisFluxUncertaintyHelper::Initialize(std::string const &filename,
                                              bool verbose) {

  std::string ND_detector_tag = "ND";
  std::string ND_SpecHCRun_detector_tag = "_specrun_";
  std::string FD_detector_tag = "FD";
  std::string nu_mode_beam_tag = "nu";
  std::string nubar_mode_beam_tag = "nubar";
  std::string numu_species_tag = "numu";
  std::string nue_species_tag = "nue";
  std::string numubar_species_tag = "numubar";
  std::string nuebar_species_tag = "nuebar";

  static std::string const location_tags[] = {ND_detector_tag, FD_detector_tag};
  static std::string const beam_mode_tags[] = {nu_mode_beam_tag,
                                               nubar_mode_beam_tag};
  static std::string const species_tags[] = {numu_species_tag, nue_species_tag,
                                             numubar_species_tag,
                                             nuebar_species_tag};
  static std::string const spec_run_tags[] = {"_", ND_SpecHCRun_detector_tag};

  if (verbose) {
    std::cout << "Reading inputs from " << filename << std::endl;
  }

  TFile *inpF = new TFile(filename.c_str());
  if (!inpF || !inpF->IsOpen()) {
    std::cout << "[ERROR]: Couldn't open input file: " << filename << std::endl;
    exit(1);
  }

  {
    focussing.OffAxisTAxes.clear();
    focussing.NDuncerts.clear();
    focussing.FDuncerts.clear();
    focussing.UncertLabels.clear();

    std::string input_dir = "FluxParameters/Focussing";
    TDirectory *d = inpF->GetDirectory(input_dir.c_str());
    if (!d) {
      std::cout << "[ERROR]: Couldn't open directory : " << input_dir
                << " in input file: " << filename << std::endl;
      exit(1);
    }

    int param_id = 0;

    for (TObject *key : *d->GetListOfKeys()) {

      std::string syst(key->GetName());

      if (verbose) {
        std::cout << "  Found parameter directory: " << input_dir << "/" << syst
                  << std::endl;
      }

      int nucfg = kND_numu_numode; // = 0

      std::stringstream input_dir_i("");
      input_dir_i << input_dir << (input_dir.size() ? "/" : "") << syst << "/";
      TDirectory *param_d = inpF->GetDirectory(input_dir_i.str().c_str());

      std::unique_ptr<TAxis> oa_axis =
          std::unique_ptr<TAxis>((TAxis *)param_d->Get("OffAxisTAxis"));
      if (param_id &&
          !check_axis(focussing.OffAxisTAxes.front().get(), oa_axis.get())) {
        exit(1);
      }
      focussing.OffAxisTAxes.emplace_back(std::move(oa_axis));

      focussing.NDuncerts.emplace_back();
      focussing.FDuncerts.emplace_back();

      std::vector<std::vector<std::unique_ptr<TH1>>> param_NDuncerts;
      std::vector<std::unique_ptr<TH1>> param_FDuncerts;

      focussing.UncertLabels.push_back(syst);

      for (size_t lt_it = 0; lt_it < 2; ++lt_it) {
        std::string const &location_tag = location_tags[lt_it];
        for (size_t bm_it = 0; bm_it < 2; ++bm_it) {
          std::string const &beam_mode_tag = beam_mode_tags[bm_it];
          for (size_t sr_it = 0; sr_it < 2; ++sr_it) {
            std::string const &spec_run_tag = spec_run_tags[sr_it];
            for (size_t sp_it = 0; sp_it < 4; ++sp_it) {
              std::string const &species_tag = species_tags[sp_it];

              std::string hname = location_tag + "_" + beam_mode_tag +
                                  spec_run_tag + species_tag;

              if (verbose) {
                std::cout << "    Reading nu configuration: " << hname
                          << std::endl;
              }

              focussing.NDuncerts.back().emplace_back();
              focussing.FDuncerts.back().emplace_back(nullptr);

              if (nucfg < kFD_numu_numode) {
                // Is ND (Any horn current for now)
                // ND dir name is the same as the hist name
                std::string nd_dir = input_dir_i.str() + hname;
                std::vector<std::unique_ptr<TH1>> AllOffAxisShifts =
                    GetNDOffAxisShifts(inpF, nd_dir, hname);

                // check that for a given off axis slice, the energy binning is
                // the same for every parameter
                if (param_id) {
                  if (focussing.NDuncerts.front().at(nucfg).size() !=
                      AllOffAxisShifts.size()) {
                    std::cout << "[ERROR]: Found differing number of off axis "
                                 "slices for two flux focussing parameters: "
                              << focussing.UncertLabels.front() << " and "
                              << syst << std::endl;
                  }
                  for (size_t oab_i = 0; oab_i < AllOffAxisShifts.size();
                       ++oab_i) {
                    if (!check_axis(focussing.NDuncerts.front()
                                        .at(nucfg)
                                        .at(oab_i)
                                        ->GetXaxis(),
                                    AllOffAxisShifts.at(oab_i)->GetXaxis())) {
                      std::cout << "[ERROR]: Found differing energy binning "
                                   "for off axis bin "
                                << oab_i
                                << "slices for two flux focussing parameters: "
                                << focussing.UncertLabels.front() << " and "
                                << syst << std::endl;
                      exit(1);
                    }
                  }
                }

                if (verbose) {
                  std::cout << "      Found " << AllOffAxisShifts.size()
                            << " off axis bins!" << std::endl;
                }

                focussing.NDuncerts.back().at(nucfg) =
                    std::move(AllOffAxisShifts);
                nucfg += 1;
              } else if (spec_run_tag == "_") { // Is FD and not 280kA run
                std::unique_ptr<TH1> h_in =
                    std::unique_ptr<TH1>((TH1 *)param_d->Get(hname.c_str()));
                if (!h_in) {
                  std::cout << "[WARN] Cannot find " << hname << std::endl;
                  continue;
                }

                h_in->SetDirectory(nullptr);

                if (verbose) {
                  std::cout << "      Found an FD variation." << std::endl;
                }

                focussing.FDuncerts.back().at(nucfg) = std::move(h_in);
                nucfg += 1;
              }
            }
          }
        }
      }
      param_id += 1;
    }
    focussing.NParams = param_id;
  }

  { // a horrific stop gap until we have inputs with matching binning

    std::string ND_SpecHCRun_detector_tag = "_280kA_";

    static std::string const spec_run_tags[] = {"_", ND_SpecHCRun_detector_tag};

    hadprod.OffAxisTAxes.clear();
    hadprod.NDuncerts.clear();
    hadprod.FDuncerts.clear();

    std::string input_dir = "FluxParameters/HadronProduction";
    TDirectory *d = inpF->GetDirectory(input_dir.c_str());

    if (!d) {
      std::cout << "[ERROR]: Couldn't open directory : " << input_dir
                << " in input file: " << filename << std::endl;
      exit(1);
    }

    int param_id = 0;

    for (TObject *key : *d->GetListOfKeys()) {

      std::string syst(key->GetName());

      if (verbose) {
        std::cout << "  Found parameter directory: " << input_dir << "/" << syst
                  << std::endl;
      }

      int nucfg = kND_numu_numode; // = 0

      std::stringstream input_dir_i("");
      input_dir_i << input_dir << (input_dir.size() ? "/" : "") << syst << "/";
      TDirectory *param_d = inpF->GetDirectory(input_dir_i.str().c_str());

      hadprod.NDuncerts.emplace_back();
      hadprod.FDuncerts.emplace_back();
      hadprod.OffAxisTAxes.emplace_back();

      std::vector<std::vector<std::unique_ptr<TH1>>> param_NDuncerts;
      std::vector<std::unique_ptr<TH1>> param_FDuncerts;

      for (size_t lt_it = 0; lt_it < 2; ++lt_it) {
        std::string const &location_tag = location_tags[lt_it];
        for (size_t bm_it = 0; bm_it < 2; ++bm_it) {
          std::string const &beam_mode_tag = beam_mode_tags[bm_it];
          for (size_t sr_it = 0; sr_it < 2; ++sr_it) {
            std::string const &spec_run_tag = spec_run_tags[sr_it];
            for (size_t sp_it = 0; sp_it < 4; ++sp_it) {
              std::string const &species_tag = species_tags[sp_it];

              std::string hname = location_tag + "_" + beam_mode_tag +
                                  spec_run_tag + species_tag;

              if (verbose) {
                std::cout << "    Reading nu configuration: " << hname
                          << std::endl;
              }

              hadprod.NDuncerts.back().emplace_back();
              hadprod.FDuncerts.back().emplace_back(nullptr);

              if (nucfg < kFD_numu_numode) {
                // Is ND (Any horn current for now)
                // ND dir name is the same as the hist name
                std::string nd_dir = input_dir_i.str() + hname;
                std::vector<std::unique_ptr<TH1>> AllOffAxisShifts =
                    GetNDOffAxisShifts(inpF, nd_dir, hname);

                std::unique_ptr<TAxis> oa_axis = std::unique_ptr<TAxis>(
                    (TAxis *)inpF->GetDirectory(nd_dir.c_str())
                        ->Get("OffAxisTAxis"));
                if (param_id &&
                    !check_axis(hadprod.OffAxisTAxes.front()[nucfg].get(),
                                oa_axis.get())) {
                  exit(1);
                }
                hadprod.OffAxisTAxes.back().emplace_back(std::move(oa_axis));

                if (param_id &&
                    !check_axis(
                        hadprod.NDuncerts.front().at(nucfg).front()->GetXaxis(),
                        AllOffAxisShifts.front()->GetXaxis())) {
                  exit(1);
                }

                if (verbose) {
                  std::cout << "      Found " << AllOffAxisShifts.size()
                            << " off axis bins!" << std::endl;
                }

                hadprod.NDuncerts.back().at(nucfg) =
                    std::move(AllOffAxisShifts);
                nucfg += 1;
              } else if (spec_run_tag == "_") { // Is FD and not 280kA run
                std::unique_ptr<TH1> h_in =
                    std::unique_ptr<TH1>((TH1 *)param_d->Get(hname.c_str()));
                if (!h_in) {
                  std::cout << "[WARN] Cannot find " << hname << std::endl;
                  continue;
                }

                h_in->SetDirectory(nullptr);

                if (verbose) {
                  std::cout << "      Found an FD variation." << std::endl;
                }

                hadprod.FDuncerts.back().at(nucfg) = std::move(h_in);
                nucfg += 1;
              }
            }
          }
        }
      }
      param_id += 1;
    }
    hadprod.NPCAComponents = param_id;
  }
}

int OffAxisFluxUncertaintyHelper::GetNuConfig(int nu_pdg, bool IsND,
                                              bool IsNuMode,
                                              bool isSpecHCRun) const {

  int nucfg = kUnhandled;

  switch (nu_pdg) {
  case 14: {
    if (IsND) {
      nucfg = IsNuMode
                  ? (isSpecHCRun ? kND_SpecHCRun_numu_numode : kND_numu_numode)
                  : (isSpecHCRun ? kND_SpecHCRun_numu_nubarmode
                                 : kND_numu_nubarmode);
    } else {
      nucfg = IsNuMode ? kFD_numu_numode : kFD_numu_nubarmode;
    }
    break;
  }
  case -14: {
    if (IsND) {
      nucfg = IsNuMode ? (isSpecHCRun ? kND_SpecHCRun_numubar_numode
                                      : kND_numubar_numode)
                       : (isSpecHCRun ? kND_SpecHCRun_numubar_nubarmode
                                      : kND_numubar_nubarmode);
    } else {
      nucfg = IsNuMode ? kFD_numubar_numode : kFD_numubar_nubarmode;
    }
    break;
  }
  case 12: {
    if (IsND) {
      nucfg =
          IsNuMode
              ? (isSpecHCRun ? kND_SpecHCRun_nue_numode : kND_nue_numode)
              : (isSpecHCRun ? kND_SpecHCRun_nue_nubarmode : kND_nue_nubarmode);
    } else {
      nucfg = IsNuMode ? kFD_nue_numode : kFD_nue_nubarmode;
    }
    break;
  }
  case -12: {
    if (IsND) {
      nucfg = IsNuMode ? (isSpecHCRun ? kND_SpecHCRun_nuebar_numode
                                      : kND_nuebar_numode)
                       : (isSpecHCRun ? kND_SpecHCRun_nuebar_nubarmode
                                      : kND_nuebar_nubarmode);
    } else {
      nucfg = IsNuMode ? kFD_nuebar_numode : kFD_nuebar_nubarmode;
    }
    break;
  }
  }

  return nucfg;
}

std::vector<std::unique_ptr<TH1>>
OffAxisFluxUncertaintyHelper::GetNDOffAxisShifts(TFile *f, std::string nd_dir,
                                                 std::string hname) const {
  std::vector<std::unique_ptr<TH1>> OffAxisUncerts;

  TDirectory *d = f->GetDirectory(nd_dir.c_str());
  if (!d) {
    std::cout << "[ERROR]: Failed to open directory: " << nd_dir << std::endl;
    exit(1);
  }
  int n_offaxis = d->GetListOfKeys()->GetSize();

  for (int oa = 0; oa < n_offaxis; oa++) {
    std::string hname_oa = hname + "_" + std::to_string(oa);
    std::unique_ptr<TH1> h_oa =
        std::unique_ptr<TH1>((TH1 *)d->Get(hname_oa.c_str()));
    if (!h_oa) {
      break;
    }
    h_oa->SetDirectory(nullptr);
    OffAxisUncerts.emplace_back(std::move(h_oa));
  }

  return OffAxisUncerts;
}

int OffAxisFluxUncertaintyHelper::GetFocussingBin(int nu_config, double enu_GeV,
                                                  double off_axis_pos_m) const {

  if (nu_config < kFD_numu_numode) {
    // Sign flip in off axis position
    int bin_oa = focussing.OffAxisTAxes.front()->FindFixBin(off_axis_pos_m);
    if ((bin_oa == 0) ||
        (bin_oa == (focussing.OffAxisTAxes.front()->GetNbins() + 1))) {
      // flow bin
      return kInvalidBin;
    }
    bin_oa--;
    auto &ehist = focussing.NDuncerts.front().at(nu_config).at(bin_oa);
    int bin_e = ehist->FindFixBin(enu_GeV);
    if ((bin_e == 0) || (bin_e == (ehist->GetXaxis()->GetNbins() + 1))) {
      // flow bin
      return kInvalidBin;
    }
    return bin_oa * 1000 + bin_e;

  } else {
    if (focussing.FDuncerts.front().size() < 1) {
      std::cout << "[WARN] no param_id" << std::endl;
    }
    if (!focussing.FDuncerts.front().at(nu_config)) {
      std::cout << "[WARN] no nu_config = " << nu_config << std::endl;
    }
    auto &ehist = focussing.FDuncerts.front().at(nu_config);
    int bin_e = ehist->FindFixBin(enu_GeV);
    if ((bin_e == 0) || (bin_e == (ehist->GetXaxis()->GetNbins() + 1))) {
      // flow bin
      return kInvalidBin;
    }
    return bin_e;
  }
}

int OffAxisFluxUncertaintyHelper::GetHadProdBin(int nu_config, double enu_GeV,
                                                double off_axis_pos_m) const {
  if (nu_config < kFD_numu_numode) {
    int bin_oa =
        hadprod.OffAxisTAxes.front()[nu_config]->FindFixBin(off_axis_pos_m);
    if ((bin_oa == 0) ||
        (bin_oa == (hadprod.OffAxisTAxes.front()[nu_config]->GetNbins() + 1))) {
      // flow bin
      return kInvalidBin;
    }
    bin_oa--;
    auto &ehist = hadprod.NDuncerts.front().at(nu_config).at(bin_oa);
    int bin_e = ehist->FindFixBin(enu_GeV);
    if ((bin_e == 0) || (bin_e == (ehist->GetXaxis()->GetNbins() + 1))) {
      // flow bin
      return kInvalidBin;
    }
    return bin_oa * 1000 + bin_e;

  } else {
    if (hadprod.FDuncerts.front().size() < 1) {
      std::cout << "[WARN] no param_id" << std::endl;
    }
    if (!hadprod.FDuncerts.front().at(nu_config)) {
      std::cout << "[WARN] no nu_config = " << nu_config << std::endl;
    }
    auto &ehist = hadprod.FDuncerts.front().at(nu_config);
    int bin_e = ehist->FindFixBin(enu_GeV);
    if ((bin_e == 0) || (bin_e == (ehist->GetXaxis()->GetNbins() + 1))) {
      // flow bin
      return kInvalidBin;
    }
    return bin_e;
  }
}

double OffAxisFluxUncertaintyHelper::GetFluxFocussingWeight(size_t param_id,
                                                            double param_val,
                                                            int nucfg,
                                                            int bin) const {
  if (nucfg == kUnhandled) {
    return 1;
  }

  if (bin == kInvalidBin) {
    return 1;
  }

  if (nucfg < kFD_numu_numode) {
    int bin_oa = bin / 1000;
    int bin_e = bin % 1000;
    return 1 + param_val * focussing.NDuncerts.at(param_id)
                               .at(nucfg)
                               .at(bin_oa)
                               ->GetBinContent(bin_e);
  } else {
    return 1 +
           param_val *
               focussing.FDuncerts.at(param_id).at(nucfg)->GetBinContent(bin);
  }
}

double OffAxisFluxUncertaintyHelper::GetFluxHadProdWeight(size_t param_id,
                                                          double param_val,
                                                          int nucfg,
                                                          int bin) const {
  if (nucfg == kUnhandled) {
    return 1;
  }

  if (bin == kInvalidBin) {
    return 1;
  }

  if (nucfg < kFD_numu_numode) {
    int bin_oa = bin / 1000;
    int bin_e = bin % 1000;
    return 1 + param_val * hadprod.NDuncerts.at(param_id)
                               .at(nucfg)
                               .at(bin_oa)
                               ->GetBinContent(bin_e);
  } else {
    return 1 + param_val *
                   hadprod.FDuncerts.at(param_id).at(nucfg)->GetBinContent(bin);
  }
}

std::vector<double>
OffAxisFluxUncertaintyHelper::GetFluxFocussingOffAxisBinning() {
  std::vector<double> rtn{focussing.OffAxisTAxes[0]->GetBinLowEdge(1)};
  for (int i = 0; i < focussing.OffAxisTAxes[0]->GetNbins(); ++i) {
    rtn.push_back(focussing.OffAxisTAxes[0]->GetBinUpEdge(i + 1));
  }
  return rtn;
}
std::vector<double>
OffAxisFluxUncertaintyHelper::GetFluxHadProdOffAxisBinning(int nu_config) {
  std::vector<double> rtn{
      hadprod.OffAxisTAxes[0][nu_config]->GetBinLowEdge(1)};
  for (int i = 0; i < hadprod.OffAxisTAxes[0][nu_config]->GetNbins(); ++i) {
    rtn.push_back(hadprod.OffAxisTAxes[0][nu_config]->GetBinUpEdge(i + 1));
  }
  return rtn;
}
