// #include "Samples/MaCh3DUNEFactory.h"
// #include "Samples/StructsDUNE.h"
// #include "Samples/SampleHandlerBase.h"
// #include "Fitters/FitterBase.h"
// #include "Manager/Manager.h"
// #include "Parameters/ParameterHandlerBase.h"

// #include <iomanip>
// #include <iostream>
// #include <memory>
// #include <vector>
// #include <string>
// #include <cmath>
// #include <algorithm>
// #include <filesystem>
// #include <map>
// #include <sstream>

// // ---------------------------------------------------------------------------
// // Safely get-or-create a (possibly multi-level) directory inside a TFile.
// // ---------------------------------------------------------------------------
// TDirectory* GetOrMakeDir(TDirectory* top, const std::string& path)
// {
//     TDirectory* dir = top;
//     std::stringstream ss(path);
//     std::string level;
//     while (std::getline(ss, level, '/')) {
//         if (level.empty()) continue;
//         TDirectory* sub = dir->GetDirectory(level.c_str());
//         if (!sub) sub = dir->mkdir(level.c_str());
//         if (!sub) {
//             std::cerr << "[ERROR] Failed to create/get directory '" << level
//                       << "' while building path '" << path << "'" << std::endl;
//             return nullptr;
//         }
//         dir = sub;
//     }
//     return dir;
// }

// // ---------------------------------------------------------------------------
// // Write a 2D histogram (True vs Reco neutrino energy) into fout
// // p < -1 => FakeData Asimov, p < 0 => Asimov, p >= 0 => posterior toy index p
// // ---------------------------------------------------------------------------
// void MakeSpectaVariations2D(SampleHandlerFD*   pdf,
//                              const std::string& varX,
//                              const std::string& varY,
//                              TFile*             fout,
//                              const std::string& ND_or_FD,
//                              const std::string& pdfTitle,
//                              int                p,
//                              int                sample_idx = 0,
//                              std::map<std::string, std::vector<std::vector<std::vector<double>>>>* posteriorStore = nullptr)
// {
//     std::vector<KinematicCut> emptyEventCuts;
//     std::vector<KinematicCut> emptySubEventCuts;

//     TH2* h = nullptr;
//     try {
//         h = pdf->Get2DVarHist(sample_idx, varX, varY,
//                                emptyEventCuts, 0, nullptr, nullptr,
//                                emptySubEventCuts);
//     } catch (const std::exception& e) {
//         std::cerr << "[WARN] Skipping " << varX << " vs " << varY << " for " << pdfTitle
//                    << " (sample_idx=" << sample_idx << "): " << e.what() << std::endl;
//         return;
//     }
//     if (!h) {
//         std::cerr << "[WARN] Null TH2 for " << varX << " vs " << varY
//                    << " from " << pdfTitle << std::endl;
//         return;
//     }

//     TH2D* h2d = dynamic_cast<TH2D*>(h);
//     if (!h2d) {
//         std::cerr << "[ERROR] Could not cast to TH2D for " << varX << " vs " << varY << std::endl;
//         delete h;
//         return;
//     }

//     std::string cleanTitle = pdfTitle;
//     std::replace(cleanTitle.begin(), cleanTitle.end(), ' ', '_');
//     std::replace(cleanTitle.begin(), cleanTitle.end(), '/', '_');

//     std::string baseDir, histName;
//     if (p < -1) {
//         baseDir  = "AsimovFakeData/" + ND_or_FD + "/2D/" + cleanTitle;
//         histName = Form("%s_%s_%s_vs_%s_AsimovFakeData",
//                          ND_or_FD.c_str(), cleanTitle.c_str(), varX.c_str(), varY.c_str());
//     } else if (p < 0) {
//         baseDir  = "Asimov/" + ND_or_FD + "/2D/" + cleanTitle;
//         histName = Form("%s_%s_%s_vs_%s_Asimov",
//                          ND_or_FD.c_str(), cleanTitle.c_str(), varX.c_str(), varY.c_str());
//     } else {
//         baseDir  = ND_or_FD + "/2D/" + cleanTitle;
//         histName = Form("%s_%s_%s_vs_%s_posterior_toy_%03d",
//                          ND_or_FD.c_str(), cleanTitle.c_str(), varX.c_str(), varY.c_str(), p);
//     }

//     TDirectory* dir = GetOrMakeDir(fout, baseDir);
//     if (!dir) {
//         std::cerr << "[WARN] Skipping 2D histogram for " << varX << " vs " << varY
//                    << " from " << pdfTitle << " (could not create directory " << baseDir << ")" << std::endl;
//         delete h2d;
//         return;
//     }
//     dir->cd();

//     TH2D* clone = static_cast<TH2D*>(h2d->Clone(histName.c_str()));
//     clone->SetDirectory(dir);
//     clone->Write();

//     if (p >= 0 && posteriorStore) {
//         std::string key = ND_or_FD + "_" + cleanTitle + "_" + varX + "_vs_" + varY;

//         int nX = clone->GetNbinsX();
//         int nY = clone->GetNbinsY();
//         auto& grid = (*posteriorStore)[key];
//         if (grid.empty()) {
//             grid.resize(nX);
//             for (int ix = 0; ix < nX; ++ix) grid[ix].resize(nY);
//         }
//         for (int ix = 1; ix <= nX; ++ix)
//             for (int iy = 1; iy <= nY; ++iy)
//                 grid[ix-1][iy-1].push_back(clone->GetBinContent(ix, iy));
//     }

//     delete clone;
//     delete h2d;
//     fout->cd();
// }

// // ---------------------------------------------------------------------------
// int main(int argc, char* argv[])
// {
//     if (argc == 1) {
//         std::cout << "Usage: bin/ config.cfg" << std::endl;
//         return 1;
//     }

//     TH1::AddDirectory(kFALSE);

//     // -----------------------------------------------------------------------
//     // Configuration
//     // -----------------------------------------------------------------------
//     auto fitMan = std::make_unique<Manager>(argv[1]);

//     auto PosteriorFile =
//         Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFiles"], __FILE__, __LINE__);
//     auto burn_in =
//         Get<unsigned int>(fitMan->raw()["General"]["MCMC"]["BurnInSteps"], __FILE__, __LINE__);
//     int no_times_sampling_posterior =
//         Get<int>(fitMan->raw()["Predictive"]["SamplePosterior"], __FILE__, __LINE__);
//     auto OscPars =
//         GetFromManager<std::vector<double>>(fitMan->raw()["General"]["OscillationParameters"], {});

//     // -----------------------------------------------------------------------
//     // Output file
//     // -----------------------------------------------------------------------
//     std::filesystem::path outDir = std::filesystem::path(PosteriorFile).parent_path();
//     std::string OutFileName = (outDir / "Posteriorpredictive_out_2D.root").string();
//     std::cout << "[INFO] Output file: " << OutFileName << std::endl;

//     TFile* fOut = new TFile(OutFileName.c_str(), "RECREATE");

//     // -----------------------------------------------------------------------
//     // Build PDFs and parameter handler
//     // -----------------------------------------------------------------------
//     ParameterHandlerGeneric* xsec = nullptr;
//     std::vector<SampleHandlerFD*> DUNEPdfs;
//     MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

//     auto ClassifyNDorFD = [](SampleHandlerFD* pdf) -> std::string {
//         std::string name = pdf->GetName();
//         if (name.find("ND") != std::string::npos) return "ND";
//         if (name.find("FD") != std::string::npos) return "FD";
//         return "FD";
//     };

//     auto ExtractGroup = [&](const std::vector<double>& full, const std::string& group) {
//         std::vector<double> out;
//         out.reserve(xsec->GetNumParFromGroup(group));
//         for (int i = 0; i < (int)full.size(); ++i) {
//             if (xsec->IsParFromGroup(i, group)) out.push_back(full[i]);
//         }
//         return out;
//     };

//     const int nXsecPars = xsec->GetNumParFromGroup("Xsec")
//                          + xsec->GetNumParFromGroup("Flux");

//     // Reset Xsec and Flux to prior, Osc to configured values
//     xsec->SetGroupOnlyParameters("Xsec");
//     xsec->SetGroupOnlyParameters("Flux");
//     xsec->SetGroupOnlyParameters("Osc",  OscPars);

//     int nFakeDataPars = xsec->GetNumParFromGroup("FakeData");
//     std::vector<double> fakeDataZero(nFakeDataPars, 0.0);
//     xsec->SetGroupOnlyParameters("FakeData", fakeDataZero);

//     // -----------------------------------------------------------------------
//     // Asimov 2D reference (nominal parameters)
//     // -----------------------------------------------------------------------
//     std::cout << "[INFO] Generating Asimov 2D spectra..." << std::endl;

//     for (auto& pdf : DUNEPdfs) {
//         std::string ND_or_FD = ClassifyNDorFD(pdf);
//         pdf->Reweight();

//         for (int s = 0; s < pdf->GetNsamples(); ++s) {
//             std::string pdfTitle = pdf->GetSampleTitle(s);
//             if (pdfTitle.empty()) pdfTitle = "FHC_numu_asimov";

//             MakeSpectaVariations2D(pdf, "TrueNeutrinoEnergy", "RecoNeutrinoEnergy",
//                                     fOut, ND_or_FD, pdfTitle, -1, s);
//         }
//     }
//     fOut->cd();

//     xsec->SetGroupOnlyParameters("FakeData", fakeDataZero);
//     std::cout << "[INFO] Reset FakeData parameters to 0 for posterior loop." << std::endl;

//     // -----------------------------------------------------------------------
//     // Open MCMC chain
//     // -----------------------------------------------------------------------
//     TChain* mcmc = new TChain("posteriors");
//     mcmc->Add(PosteriorFile.c_str());
//     Long64_t nEntries = mcmc->GetEntries();
//     if (nEntries == 0) {
//         std::cerr << "[ERROR] MCMC tree has no entries!" << std::endl;
//         return 1;
//     }
//     std::cout << "[INFO] MCMC chain has " << nEntries << " entries." << std::endl;

//     // -----------------------------------------------------------------------
//     // Binary-search for first post-burn-in entry
//     // -----------------------------------------------------------------------
//     UInt_t mcmc_step = 0;
//     mcmc->SetBranchStatus("*",    0);
//     mcmc->SetBranchStatus("step", 1);
//     mcmc->SetBranchAddress("step", &mcmc_step);

//     Long64_t lo = 0, hi = nEntries - 1, firstValid = -1;
//     while (lo <= hi) {
//         Long64_t mid = (lo + hi) / 2;
//         mcmc->GetEntry(mid);
//         if (mcmc_step >= burn_in) { firstValid = mid; hi = mid - 1; }
//         else                      { lo = mid + 1; }
//     }

//     if (firstValid < 0) {
//         std::cerr << "[ERROR] No MCMC entries survive burn-in cut of " << burn_in << std::endl;
//         return 1;
//     }

//     Long64_t nValid = nEntries - firstValid;
//     std::cout << "[INFO] First post-burn-in entry: " << firstValid
//               << "  (" << nValid << " valid entries)" << std::endl;

//     // -----------------------------------------------------------------------
//     // Discover how many param_N branches exist
//     // -----------------------------------------------------------------------
//     mcmc->SetBranchStatus("*", 1);

//     int nBranchesFound = 0;
//     for (int i = 0; i < 2000; ++i) {
//         if (mcmc->GetBranch(Form("param_%d", i))) ++nBranchesFound;
//         else break;
//     }
//     std::cout << "[INFO] Found " << nBranchesFound << " param branches in MCMC." << std::endl;

//     const int nParsToRead = std::min(nBranchesFound, nXsecPars);

//     mcmc->SetBranchStatus("*",    0);
//     mcmc->SetBranchStatus("step", 1);
//     mcmc->SetBranchAddress("step", &mcmc_step);

//     std::vector<Double_t> xsec_tmp(nXsecPars, 0.0);
//     for (int i = 0; i < nParsToRead; ++i) {
//         TString bname = Form("param_%d", i);
//         mcmc->SetBranchStatus(bname, 1);
//         mcmc->SetBranchAddress(bname, &xsec_tmp[i]);
//     }

//     // -----------------------------------------------------------------------
//     // Pre-select random entries (sorted for sequential disk access)
//     // -----------------------------------------------------------------------
//     std::cout << "[INFO] Pre-selecting " << no_times_sampling_posterior
//               << " posterior samples..." << std::endl;

//     auto rnd = std::make_unique<TRandom3>(0);
//     std::vector<Long64_t> chosenEntries(no_times_sampling_posterior);
//     for (int p = 0; p < no_times_sampling_posterior; ++p)
//         chosenEntries[p] = firstValid + (Long64_t)rnd->Integer(nValid);

//     std::sort(chosenEntries.begin(), chosenEntries.end());

//     std::vector<std::vector<Double_t>> sampledParams(no_times_sampling_posterior,
//                                                        std::vector<Double_t>(nXsecPars, 0.0));
//     for (int p = 0; p < no_times_sampling_posterior; ++p) {
//         mcmc->GetEntry(chosenEntries[p]);
//         sampledParams[p].assign(xsec_tmp.begin(), xsec_tmp.end());
//     }

//     std::cout << "[INFO] Finished reading posterior samples. "
//               << "Starting reweighting loop..." << std::endl;

//     // -----------------------------------------------------------------------
//     // Toy loop: combined systematics only, 2D output only
//     // -----------------------------------------------------------------------
//     std::map<std::string, std::vector<std::vector<std::vector<double>>>> posteriorStore2D;

//     for (int p = 0; p < no_times_sampling_posterior; ++p) {
//         if (p % 50 == 0)
//             std::cout << "[INFO] Toy " << p << " / " << no_times_sampling_posterior << std::endl;

//         const auto xsecVals = ExtractGroup(sampledParams[p], "Xsec");
//         const auto fluxVals = ExtractGroup(sampledParams[p], "Flux");

//         xsec->SetGroupOnlyParameters("Xsec", xsecVals);
//         xsec->SetGroupOnlyParameters("Flux", fluxVals);
//         xsec->SetGroupOnlyParameters("Osc",  OscPars);

//         for (auto& pdf : DUNEPdfs) {
//             std::string ND_or_FD = ClassifyNDorFD(pdf);
//             pdf->Reweight();

//             for (int s = 0; s < pdf->GetNsamples(); ++s) {
//                 std::string pdfTitle = pdf->GetSampleTitle(s);
//                 if (pdfTitle.empty()) pdfTitle = "FHC_numu_asimov";

//                 MakeSpectaVariations2D(pdf, "TrueNeutrinoEnergy", "RecoNeutrinoEnergy",
//                                         fOut, ND_or_FD, pdfTitle, p, s, &posteriorStore2D);
//             }
//         }
//         fOut->cd();
//     }

//     // -----------------------------------------------------------------------
//     // Posterior predictive summary (mean + stddev per bin)
//     // -----------------------------------------------------------------------
//     std::cout << "[INFO] Building predictive summary histograms..." << std::endl;

//     for (auto& [key, grid] : posteriorStore2D) {
//         int nX = grid.size();
//         if (nX == 0) continue;
//         int nY = grid[0].size();

//         TH2D hMean(("mean_" + key).c_str(), "", nX, 0, nX, nY, 0, nY);
//         TH2D hStdDev(("stddev_" + key).c_str(), "", nX, 0, nX, nY, 0, nY);

//         for (int ix = 0; ix < nX; ++ix) {
//             for (int iy = 0; iy < nY; ++iy) {
//                 const auto& vals = grid[ix][iy];
//                 double mean = 0.0;
//                 for (double v : vals) mean += v;
//                 mean /= vals.size();

//                 double var = 0.0;
//                 for (double v : vals) var += (v - mean) * (v - mean);
//                 var /= vals.size();

//                 hMean.SetBinContent(ix + 1, iy + 1, mean);
//                 hStdDev.SetBinContent(ix + 1, iy + 1, std::sqrt(var));
//             }
//         }

//         TDirectory* predDir = GetOrMakeDir(fOut, "PosteriorPredictive/2D");
//         if (predDir) {
//             predDir->cd();
//             hMean.Write();
//             hStdDev.Write();
//         }
//         fOut->cd();
//     }

//     // -----------------------------------------------------------------------
//     // Cleanup
//     // -----------------------------------------------------------------------
//     std::cout << "[INFO] Writing output file: " << OutFileName << std::endl;
//     fOut->Write();
//     fOut->Close();

//     delete mcmc;
//     delete xsec;
//     for (auto* sample : DUNEPdfs) delete sample;

//     return 0;
// }


#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Samples/SampleHandlerBase.h"
#include "Fitters/FitterBase.h"
#include "Manager/Manager.h"
#include "Parameters/ParameterHandlerBase.h"

#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <map>
#include <sstream>
#include <regex>

// ---------------------------------------------------------------------------
// Safely get-or-create a (possibly multi-level) directory inside a TFile.
// ---------------------------------------------------------------------------
TDirectory* GetOrMakeDir(TDirectory* top, const std::string& path)
{
    TDirectory* dir = top;
    std::stringstream ss(path);
    std::string level;
    while (std::getline(ss, level, '/')) {
        if (level.empty()) continue;
        TDirectory* sub = dir->GetDirectory(level.c_str());
        if (!sub) sub = dir->mkdir(level.c_str());
        if (!sub) {
            std::cerr << "[ERROR] Failed to create/get directory '" << level
                      << "' while building path '" << path << "'" << std::endl;
            return nullptr;
        }
        dir = sub;
    }
    return dir;
}

// ---------------------------------------------------------------------------
// baseKey is always "<DET>_<sample>" by construction (DET = ND or FD, and
// <sample> may itself contain underscores). Split it so we know which
// detector directory / sample name to use for the posteriorpredictive_plots
// -style summary layout.
// ---------------------------------------------------------------------------
static bool SplitDetAndSample(const std::string& baseKey, std::string& det, std::string& sample)
{
    static const std::regex re("^(ND|FD)_(.+)$");
    std::smatch m;
    if (!std::regex_match(baseKey, m, re)) return false;
    det    = m[1].str();
    sample = m[2].str();
    return true;
}

// ---------------------------------------------------------------------------
// Write a 2D histogram (True vs Reco neutrino energy) into fout
// p < -1 => FakeData Asimov, p == -1 => Asimov, p >= 0 => posterior toy index p
//
// posteriorStore accumulates per-bin values across posterior toys, keyed by
// "<DET>_<sample>", so the mean/stddev can be built at the end without
// re-reading anything from disk.
//
// asimovStore captures a detached clone of the p == -1 Asimov histogram,
// keyed the same way, so it can serve both as the real-binning template
// AND as the reference for the Mean/Asimov ratio plot in the Python
// plotting script -- no separate post-processing step required.
// ---------------------------------------------------------------------------
void MakeSpectaVariations2D(SampleHandlerFD*   pdf,
                             const std::string& varX,
                             const std::string& varY,
                             TFile*             fout,
                             const std::string& ND_or_FD,
                             const std::string& pdfTitle,
                             int                p,
                             int                sample_idx = 0,
                             std::map<std::string, std::vector<std::vector<std::vector<double>>>>* posteriorStore = nullptr,
                             std::map<std::string, TH2D*>* asimovStore = nullptr)
{
    std::vector<KinematicCut> emptyEventCuts;
    std::vector<KinematicCut> emptySubEventCuts;

    TH2* h = nullptr;
    try {
        h = pdf->Get2DVarHist(sample_idx, varX, varY,
                               emptyEventCuts, 0, nullptr, nullptr,
                               emptySubEventCuts);
    } catch (const std::exception& e) {
        std::cerr << "[WARN] Skipping " << varX << " vs " << varY << " for " << pdfTitle
                   << " (sample_idx=" << sample_idx << "): " << e.what() << std::endl;
        return;
    }
    if (!h) {
        std::cerr << "[WARN] Null TH2 for " << varX << " vs " << varY
                   << " from " << pdfTitle << std::endl;
        return;
    }

    TH2D* h2d = dynamic_cast<TH2D*>(h);
    if (!h2d) {
        std::cerr << "[ERROR] Could not cast to TH2D for " << varX << " vs " << varY << std::endl;
        delete h;
        return;
    }

    std::string cleanTitle = pdfTitle;
    std::replace(cleanTitle.begin(), cleanTitle.end(), ' ', '_');
    std::replace(cleanTitle.begin(), cleanTitle.end(), '/', '_');

    // Canonical grouping key: identical for the Asimov (p == -1) histogram
    // and every posterior toy of the same detector/sample/variable pair.
    const std::string key = ND_or_FD + "_" + cleanTitle + "_" + varX + "_vs_" + varY;

    std::string baseDir, histName;
    if (p < -1) {
        baseDir  = "AsimovFakeData/" + ND_or_FD + "/2D/" + cleanTitle;
        histName = Form("%s_%s_%s_vs_%s_AsimovFakeData",
                         ND_or_FD.c_str(), cleanTitle.c_str(), varX.c_str(), varY.c_str());
    } else if (p < 0) {
        baseDir  = "Asimov/" + ND_or_FD + "/2D/" + cleanTitle;
        histName = Form("%s_%s_%s_vs_%s_Asimov",
                         ND_or_FD.c_str(), cleanTitle.c_str(), varX.c_str(), varY.c_str());
    } else {
        baseDir  = ND_or_FD + "/2D/" + cleanTitle;
        histName = Form("%s_%s_%s_vs_%s_posterior_toy_%03d",
                         ND_or_FD.c_str(), cleanTitle.c_str(), varX.c_str(), varY.c_str(), p);
    }

    TDirectory* dir = GetOrMakeDir(fout, baseDir);
    if (!dir) {
        std::cerr << "[WARN] Skipping 2D histogram for " << varX << " vs " << varY
                   << " from " << pdfTitle << " (could not create directory " << baseDir << ")" << std::endl;
        delete h2d;
        return;
    }
    dir->cd();

    TH2D* clone = static_cast<TH2D*>(h2d->Clone(histName.c_str()));
    clone->SetDirectory(dir);
    clone->Write();

    if (p >= 0 && posteriorStore) {
        int nX = clone->GetNbinsX();
        int nY = clone->GetNbinsY();
        auto& grid = (*posteriorStore)[key];
        if (grid.empty()) {
            grid.resize(nX);
            for (int ix = 0; ix < nX; ++ix) grid[ix].resize(nY);
        }
        for (int ix = 1; ix <= nX; ++ix)
            for (int iy = 1; iy <= nY; ++iy)
                grid[ix-1][iy-1].push_back(clone->GetBinContent(ix, iy));
    }

    if (p == -1 && asimovStore) {
        TH2D* asimovClone = (TH2D*)clone->Clone((key + "_asimov_template").c_str());
        asimovClone->SetDirectory(nullptr);
        (*asimovStore)[key] = asimovClone;
    }

    delete clone;
    delete h2d;
    fout->cd();
}

// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    if (argc == 1) {
        std::cout << "Usage: bin/ config.cfg" << std::endl;
        return 1;
    }

    TH1::AddDirectory(kFALSE);

    // -----------------------------------------------------------------------
    // Configuration
    // -----------------------------------------------------------------------
    auto fitMan = std::make_unique<Manager>(argv[1]);

    auto PosteriorFile =
        Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFiles"], __FILE__, __LINE__);
    auto burn_in =
        Get<unsigned int>(fitMan->raw()["General"]["MCMC"]["BurnInSteps"], __FILE__, __LINE__);
    int no_times_sampling_posterior =
        Get<int>(fitMan->raw()["Predictive"]["SamplePosterior"], __FILE__, __LINE__);
    auto OscPars =
        GetFromManager<std::vector<double>>(fitMan->raw()["General"]["OscillationParameters"], {});

    // -----------------------------------------------------------------------
    // Output file
    // -----------------------------------------------------------------------
    std::filesystem::path outDir = std::filesystem::path(PosteriorFile).parent_path();
    std::string OutFileName = (outDir / "Posteriorpredictive_out_2D.root").string();
    std::cout << "[INFO] Output file: " << OutFileName << std::endl;

    TFile* fOut = new TFile(OutFileName.c_str(), "RECREATE");

    // -----------------------------------------------------------------------
    // Build PDFs and parameter handler
    // -----------------------------------------------------------------------
    ParameterHandlerGeneric* xsec = nullptr;
    std::vector<SampleHandlerFD*> DUNEPdfs;
    MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

    auto ClassifyNDorFD = [](SampleHandlerFD* pdf) -> std::string {
        std::string name = pdf->GetName();
        if (name.find("ND") != std::string::npos) return "ND";
        if (name.find("FD") != std::string::npos) return "FD";
        return "FD";
    };

    auto ExtractGroup = [&](const std::vector<double>& full, const std::string& group) {
        std::vector<double> out;
        out.reserve(xsec->GetNumParFromGroup(group));
        for (int i = 0; i < (int)full.size(); ++i) {
            if (xsec->IsParFromGroup(i, group)) out.push_back(full[i]);
        }
        return out;
    };

    const int nXsecPars = xsec->GetNumParFromGroup("Xsec")
                         + xsec->GetNumParFromGroup("Flux");

    // Reset Xsec and Flux to prior, Osc to configured values
    xsec->SetGroupOnlyParameters("Xsec");
    xsec->SetGroupOnlyParameters("Flux");
    xsec->SetGroupOnlyParameters("Osc",  OscPars);

    int nFakeDataPars = xsec->GetNumParFromGroup("FakeData");
    std::vector<double> fakeDataZero(nFakeDataPars, 0.0);
    xsec->SetGroupOnlyParameters("FakeData", fakeDataZero);

    // -----------------------------------------------------------------------
    // Asimov 2D reference (nominal parameters)
    // -----------------------------------------------------------------------
    std::cout << "[INFO] Generating Asimov 2D spectra..." << std::endl;

    // Captures the Asimov histograms in memory, keyed by "<DET>_<sample>",
    // so the summary step at the end can use them for real binning and as
    // the Mean/Asimov reference without touching disk again.
    std::map<std::string, TH2D*> asimovStore2D;

    for (auto& pdf : DUNEPdfs) {
        std::string ND_or_FD = ClassifyNDorFD(pdf);
        pdf->Reweight();

        for (int s = 0; s < pdf->GetNsamples(); ++s) {
            std::string pdfTitle = pdf->GetSampleTitle(s);
            if (pdfTitle.empty()) pdfTitle = "FHC_numu_asimov";

            MakeSpectaVariations2D(pdf, "TrueNeutrinoEnergy", "RecoNeutrinoEnergy",
                                    fOut, ND_or_FD, pdfTitle, -1, s,
                                    /*posteriorStore=*/nullptr, &asimovStore2D);
        }
    }
    fOut->cd();

    xsec->SetGroupOnlyParameters("FakeData", fakeDataZero);
    std::cout << "[INFO] Reset FakeData parameters to 0 for posterior loop." << std::endl;

    // -----------------------------------------------------------------------
    // Open MCMC chain
    // -----------------------------------------------------------------------
    TChain* mcmc = new TChain("posteriors");
    mcmc->Add(PosteriorFile.c_str());
    Long64_t nEntries = mcmc->GetEntries();
    if (nEntries == 0) {
        std::cerr << "[ERROR] MCMC tree has no entries!" << std::endl;
        return 1;
    }
    std::cout << "[INFO] MCMC chain has " << nEntries << " entries." << std::endl;

    // -----------------------------------------------------------------------
    // Binary-search for first post-burn-in entry
    // -----------------------------------------------------------------------
    UInt_t mcmc_step = 0;
    mcmc->SetBranchStatus("*",    0);
    mcmc->SetBranchStatus("step", 1);
    mcmc->SetBranchAddress("step", &mcmc_step);

    Long64_t lo = 0, hi = nEntries - 1, firstValid = -1;
    while (lo <= hi) {
        Long64_t mid = (lo + hi) / 2;
        mcmc->GetEntry(mid);
        if (mcmc_step >= burn_in) { firstValid = mid; hi = mid - 1; }
        else                      { lo = mid + 1; }
    }

    if (firstValid < 0) {
        std::cerr << "[ERROR] No MCMC entries survive burn-in cut of " << burn_in << std::endl;
        return 1;
    }

    Long64_t nValid = nEntries - firstValid;
    std::cout << "[INFO] First post-burn-in entry: " << firstValid
              << "  (" << nValid << " valid entries)" << std::endl;

    // -----------------------------------------------------------------------
    // Discover how many param_N branches exist
    // -----------------------------------------------------------------------
    mcmc->SetBranchStatus("*", 1);

    int nBranchesFound = 0;
    for (int i = 0; i < 2000; ++i) {
        if (mcmc->GetBranch(Form("param_%d", i))) ++nBranchesFound;
        else break;
    }
    std::cout << "[INFO] Found " << nBranchesFound << " param branches in MCMC." << std::endl;

    const int nParsToRead = std::min(nBranchesFound, nXsecPars);

    mcmc->SetBranchStatus("*",    0);
    mcmc->SetBranchStatus("step", 1);
    mcmc->SetBranchAddress("step", &mcmc_step);

    std::vector<Double_t> xsec_tmp(nXsecPars, 0.0);
    for (int i = 0; i < nParsToRead; ++i) {
        TString bname = Form("param_%d", i);
        mcmc->SetBranchStatus(bname, 1);
        mcmc->SetBranchAddress(bname, &xsec_tmp[i]);
    }

    // -----------------------------------------------------------------------
    // Pre-select random entries (sorted for sequential disk access)
    // -----------------------------------------------------------------------
    std::cout << "[INFO] Pre-selecting " << no_times_sampling_posterior
              << " posterior samples..." << std::endl;

    auto rnd = std::make_unique<TRandom3>(0);
    std::vector<Long64_t> chosenEntries(no_times_sampling_posterior);
    for (int p = 0; p < no_times_sampling_posterior; ++p)
        chosenEntries[p] = firstValid + (Long64_t)rnd->Integer(nValid);

    std::sort(chosenEntries.begin(), chosenEntries.end());

    std::vector<std::vector<Double_t>> sampledParams(no_times_sampling_posterior,
                                                       std::vector<Double_t>(nXsecPars, 0.0));
    for (int p = 0; p < no_times_sampling_posterior; ++p) {
        mcmc->GetEntry(chosenEntries[p]);
        sampledParams[p].assign(xsec_tmp.begin(), xsec_tmp.end());
    }

    std::cout << "[INFO] Finished reading posterior samples. "
              << "Starting reweighting loop..." << std::endl;

    // -----------------------------------------------------------------------
    // Toy loop: combined systematics only, 2D output only
    // -----------------------------------------------------------------------
    std::map<std::string, std::vector<std::vector<std::vector<double>>>> posteriorStore2D;

    for (int p = 0; p < no_times_sampling_posterior; ++p) {
        if (p % 50 == 0)
            std::cout << "[INFO] Toy " << p << " / " << no_times_sampling_posterior << std::endl;

        const auto xsecVals = ExtractGroup(sampledParams[p], "Xsec");
        const auto fluxVals = ExtractGroup(sampledParams[p], "Flux");

        xsec->SetGroupOnlyParameters("Xsec", xsecVals);
        xsec->SetGroupOnlyParameters("Flux", fluxVals);
        xsec->SetGroupOnlyParameters("Osc",  OscPars);

        for (auto& pdf : DUNEPdfs) {
            std::string ND_or_FD = ClassifyNDorFD(pdf);
            pdf->Reweight();

            for (int s = 0; s < pdf->GetNsamples(); ++s) {
                std::string pdfTitle = pdf->GetSampleTitle(s);
                if (pdfTitle.empty()) pdfTitle = "FHC_numu_asimov";

                MakeSpectaVariations2D(pdf, "TrueNeutrinoEnergy", "RecoNeutrinoEnergy",
                                        fOut, ND_or_FD, pdfTitle, p, s, &posteriorStore2D);
            }
        }
        fOut->cd();
    }

    // -----------------------------------------------------------------------
    // Posterior predictive summary
    //
    // Writes directly into the layout expected by
    // Utils/plot_posterior_predictive_summary.py and
    // Utils/plot_run_plan_uncertainty_ratio.py -- no separate ROOT macro
    // post-processing step required:
    //
    //   {det}/2D/summary/{sample}_2D_posteriorMean
    //   {det}/2D/summary/{sample}_2D_posteriorErr   (stddev across toys --
    //                                                 nuisance-parameter /
    //                                                 systematic spread)
    //   {det}/2D/summary/{sample}_2D_statErr        (sqrt(posteriorMean) --
    //                                                 genuine Poisson
    //                                                 counting uncertainty)
    //   Asimov/{det}/2D/{sample}_2D_Asimov          (nominal reference, for
    //                                                 the Mean/Asimov ratio)
    // -----------------------------------------------------------------------
    std::cout << "[INFO] Building predictive summary histograms..." << std::endl;

    for (auto& kv : posteriorStore2D) {
        const std::string& key = kv.first;
        auto& grid = kv.second;

        int nX = (int)grid.size();
        if (nX == 0) continue;
        int nY = (int)grid[0].size();

        std::string det, sample;
        if (!SplitDetAndSample(key, det, sample)) {
            std::cerr << "[WARN] Could not classify detector for '" << key
                       << "' (expected an ND_ or FD_ prefix) -- skipping summary." << std::endl;
            continue;
        }

        // Use the Asimov histogram as the binning template so the summary
        // carries real GeV bin edges instead of placeholder bin indices.
        // Falls back to a plain-index histogram (with a warning) if no
        // matching Asimov reference was captured.
        TH2D* templateHist = nullptr;
        auto asimovIt = asimovStore2D.find(key);
        if (asimovIt != asimovStore2D.end()) {
            templateHist = asimovIt->second;
        } else {
            std::cerr << "[WARN] No Asimov reference for '" << key
                       << "' -- summary histograms will use placeholder bin-index axes." << std::endl;
        }

        std::unique_ptr<TH2D> hMean, hPostErr, hStatErr;
        if (templateHist) {
            hMean.reset((TH2D*)templateHist->Clone((sample + "_2D_posteriorMean").c_str()));
            hPostErr.reset((TH2D*)templateHist->Clone((sample + "_2D_posteriorErr").c_str()));
            hStatErr.reset((TH2D*)templateHist->Clone((sample + "_2D_statErr").c_str()));
            hMean->Reset();
            hPostErr->Reset();
            hStatErr->Reset();
        } else {
            hMean.reset(new TH2D((sample + "_2D_posteriorMean").c_str(), "", nX, 0, nX, nY, 0, nY));
            hPostErr.reset(new TH2D((sample + "_2D_posteriorErr").c_str(), "", nX, 0, nX, nY, 0, nY));
            hStatErr.reset(new TH2D((sample + "_2D_statErr").c_str(), "", nX, 0, nX, nY, 0, nY));
        }
        hMean->SetDirectory(nullptr);
        hPostErr->SetDirectory(nullptr);
        hStatErr->SetDirectory(nullptr);

        for (int ix = 0; ix < nX; ++ix) {
            for (int iy = 0; iy < nY; ++iy) {
                const auto& vals = grid[ix][iy];
                if (vals.empty()) continue;

                double mean = 0.0;
                for (double v : vals) mean += v;
                mean /= vals.size();

                double var = 0.0;
                for (double v : vals) var += (v - mean) * (v - mean);
                var /= vals.size();
                double posteriorErr = std::sqrt(var);

                // Genuine Poisson counting uncertainty on the predicted
                // rate -- distinct from the toy-to-toy systematic spread.
                double statErr = (mean > 0.0) ? std::sqrt(mean) : 0.0;

                hMean->SetBinContent(ix + 1, iy + 1, mean);
                hPostErr->SetBinContent(ix + 1, iy + 1, posteriorErr);
                hStatErr->SetBinContent(ix + 1, iy + 1, statErr);
            }
        }

        TDirectory* summaryDir = GetOrMakeDir(fOut, det + "/2D/summary");
        if (summaryDir) {
            summaryDir->cd();
            hMean->Write();
            hPostErr->Write();
            hStatErr->Write();
        }
        fOut->cd();

        if (templateHist) {
            TDirectory* asimovDir = GetOrMakeDir(fOut, "Asimov/" + det + "/2D");
            if (asimovDir) {
                asimovDir->cd();
                std::unique_ptr<TH2D> hAsimov((TH2D*)templateHist->Clone((sample + "_2D_Asimov").c_str()));
                hAsimov->SetDirectory(nullptr);
                hAsimov->Write();
            }
            fOut->cd();
        }
    }

    for (auto& kv : asimovStore2D) delete kv.second;

    // -----------------------------------------------------------------------
    // Cleanup
    // -----------------------------------------------------------------------
    std::cout << "[INFO] Writing output file: " << OutFileName << std::endl;
    fOut->Write();
    fOut->Close();

    delete mcmc;
    delete xsec;
    for (auto* sample : DUNEPdfs) delete sample;

    return 0;
}