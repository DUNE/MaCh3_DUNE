// #include "Samples/MaCh3DUNEFactory.h"
// #include "Samples/StructsDUNE.h"
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

// void MakeSpectaVariations(SampleHandlerFD* pdf, const std::string& var,
//                           TFile* fout, const std::string& ND_or_FD,
//                           const std::string& pdfTitle, int p) {
//     // Refresh internal event weights
//     (void) pdf->GetMCHist(1);  // calls fill1DHist() internally
//     pdf->Reweight();
//     TH1D* h = dynamic_cast<TH1D*>(pdf->Get1DVarHist(var.c_str()));
//     if (!h) {
//         std::cerr << "[WARN] Could not get TH1D for " << var << " from " << pdfTitle << std::endl;
//         return;
//     }

//     std::string cleanTitle = pdfTitle;
//     std::replace(cleanTitle.begin(), cleanTitle.end(), ' ', '_');
//     std::replace(cleanTitle.begin(), cleanTitle.end(), '/', '_');

//     std::string baseDir, histName;
//     if (p < 0) {
//         baseDir = "Asimov/" + ND_or_FD + "/" + var;
//         histName = Form("%s_%s_%s_Asimov",
//                         ND_or_FD.c_str(), cleanTitle.c_str(), var.c_str());
//     } else {
//         baseDir = ND_or_FD + "/" + var;
//         histName = Form("%s_%s_%s_posterior_toy_%03d",
//                         ND_or_FD.c_str(), cleanTitle.c_str(), var.c_str(), p);
//     }

//     TDirectory* dir = fout->GetDirectory(baseDir.c_str());
//     if (!dir) dir = fout->mkdir(baseDir.c_str());
//     dir->cd();

//     TH1D* cloneH = static_cast<TH1D*>(h->Clone(histName.c_str()));
//     cloneH->SetDirectory(dir);
//     cloneH->Write();

//     delete cloneH;
//     //delete h; 
//     fout->cd();
// }


// void MakeSpectaVariations2D(SampleHandlerFD* pdf,
//                             TFile* fout,
//                             const std::string& ND_or_FD,
//                             const std::string& pdfTitle,
//                             int p,
//                         std::map<std::string,
//         std::vector<std::vector<std::vector<double>>>>* posteriorStore = nullptr)
// {
//     // Refresh internal event weights
//     (void) pdf->GetMCHist(2);   // ensure 2D hist is built
//     pdf->Reweight();

//     TH2D* h = dynamic_cast<TH2D*>(pdf->GetMCHist(2));
//     if (!h) {
//         std::cerr << "[WARN] Could not get TH2D from "
//                   << pdfTitle << std::endl;
//         return;
//     }

//     std::string cleanTitle = pdfTitle;
//     std::replace(cleanTitle.begin(), cleanTitle.end(), ' ', '_');
//     std::replace(cleanTitle.begin(), cleanTitle.end(), '/', '_');

//     std::string baseDir, histName;

//     if (p < 0) {
//         baseDir = "Asimov/" + ND_or_FD + "/2D";
//         histName = Form("%s_%s_2D_Asimov",
//                         ND_or_FD.c_str(), cleanTitle.c_str());
//     } else {
//         baseDir = ND_or_FD + "/2D";
//         histName = Form("%s_%s_2D_posterior_toy_%03d",
//                         ND_or_FD.c_str(), cleanTitle.c_str(), p);
//     }

//     TDirectory* dir = fout->GetDirectory(baseDir.c_str());
//     if (!dir) dir = fout->mkdir(baseDir.c_str());
//     dir->cd();

//     TH2D* cloneH = static_cast<TH2D*>(h->Clone(histName.c_str()));
//     if (p >= 0 && posteriorStore) {

//     std::string key = ND_or_FD + "_" + cleanTitle;

//     int nX = cloneH->GetNbinsX();
//     int nY = cloneH->GetNbinsY();

//     if ((*posteriorStore)[key].empty()) {
//         (*posteriorStore)[key].resize(nX);
//         for (int ix = 0; ix < nX; ++ix)
//             (*posteriorStore)[key][ix].resize(nY);
//     }

//     for (int ix = 1; ix <= nX; ++ix) {
//         for (int iy = 1; iy <= nY; ++iy) {

//             (*posteriorStore)[key][ix-1][iy-1]
//                 .push_back(cloneH->GetBinContent(ix, iy));
//         }
//     }
// }
//     cloneH->SetDirectory(dir);
//     cloneH->Write();

//     // Optional: print bin info
//     for (int ix = 1; ix <= cloneH->GetNbinsX(); ++ix) {
//         for (int iy = 1; iy <= cloneH->GetNbinsY(); ++iy) {

//             double rate    = cloneH->GetBinContent(ix, iy);
//             double statErr = cloneH->GetBinError(ix, iy);

//             std::cout
//                 << histName
//                 << " | Bin (" << ix << "," << iy << ")"
//                 << " | Rate = " << rate
//                 << " | StatErr = " << statErr
//                 << std::endl;
//         }
//     }

//     delete cloneH;
//     fout->cd();
// }


// int main(int argc, char* argv[]) {
//     if (argc == 1) {
//         std::cout << "Usage: bin/ config.cfg" << std::endl;
//         return 1;
//     }

//     // --- Manager setup
//     auto fitMan = std::make_unique<manager>(argv[1]);
//     auto PosteriorFile = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFiles"], __FILE__, __LINE__);
//     auto burn_in = Get<unsigned int>(fitMan->raw()["General"]["MCMC"]["BurnInSteps"], __FILE__, __LINE__);
//     int no_times_sampling_posterior = Get<int>(fitMan->raw()["Predictive"]["SamplePosterior"], __FILE__, __LINE__);
//     std::vector<std::string> xsecCovMatrixFile =
//         GetFromManager<std::vector<std::string>>(fitMan->raw()["General"]["Systematics"]["XsecCovFile"], {});
//     auto OscCovFile = GetFromManager<std::vector<std::string>>(fitMan->raw()["General"]["Systematics"]["OscCovFile"], {});
//     auto OscCovName = GetFromManager<std::string>(fitMan->raw()["General"]["Systematics"]["OscCovName"], "osc_cov");
//     auto OscPars = GetFromManager<std::vector<double>>(fitMan->raw()["General"]["OscillationParameters"], {});
//     double shiftAmount = GetFromManager<double>(fitMan->raw()["Predictive"]["ShiftAmount"], 0.2);

//     std::cout << "[INFO] Using shift amount ±" << shiftAmount << std::endl;

//     // --- Setup output file
//     std::filesystem::path dir = std::filesystem::path(PosteriorFile).parent_path();
//     std::string OutFileName = (dir / "Posteriorpredictive_out.root").string();

//     // Add suffix before opening
//     std::string suffix = "_posteriorpredictive_withtoys";
//     size_t dotPos = OutFileName.find_last_of('.');
//     if (dotPos != std::string::npos)
//         OutFileName.insert(dotPos, suffix);
//     else
//         OutFileName += suffix;

//     std::cout << "[INFO] Output file will be: " << OutFileName << std::endl;

//     TFile* fOut = new TFile(OutFileName.c_str(), "RECREATE");
//     MCMCProcessor Processor(PosteriorFile);
//     Processor.Initialise();

//     // --- Covariances
//     ParameterHandlerGeneric* xsec = nullptr;

//     // --- Build PDF instances
//     std::vector<SampleHandlerFD*> DUNEPdfs;
//     MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

//     // key → [ix][iy][values over throws]
//     std::map<std::string,
//     std::vector<std::vector<std::vector<double>>>> posteriorStore;

//     // --- Load posteriors
//     TChain* mcmc = new TChain("posteriors");
//     mcmc->Add(PosteriorFile.c_str());
//     Long64_t nEntries = mcmc->GetEntries();
//     if (nEntries == 0) {
//         std::cerr << "[ERROR] MCMC tree has no entries!" << std::endl;
//         return 1;
//     }

//     // Bind xsec branches
//     std::vector<double> xsec_nominal = xsec->GetPreFitValues();
//     std::vector<Double_t> xsec_tmp(xsec_nominal.size(), 0.0);
//     for (size_t i = 0; i < xsec_nominal.size(); ++i) {
//         TString bname = Form("xsec_%zu", i);
//         if (mcmc->GetBranch(bname))
//             mcmc->SetBranchAddress(bname, &xsec_tmp[i]);
//     }

//     UInt_t mcmc_step = -1;
//     mcmc->SetBranchAddress("step", &mcmc_step);

//     // --- Generate Asimov spectra
//     xsec->SetGroupOnlyParameters("Xsec", xsec_nominal);
//     //xsec->SetGroupOnlyParameters("DetSys", xsec_nominal);
//     xsec->SetGroupOnlyParameters("Osc", OscPars);

//     for (auto& pdf : DUNEPdfs) {
//         std::string pdfTitle = pdf->GetTitle();
//         if (pdfTitle.empty()) pdfTitle = "FHC_numu_asimov";

//         std::string ND_or_FD =
//             (pdfTitle.find("ND") != std::string::npos) ? "ND" :
//             (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";

//         // Create Asimov directory
//         TDirectory* asimovDir = fOut->GetDirectory("Asimov");
//         if (!asimovDir) asimovDir = fOut->mkdir("Asimov");
//         asimovDir->cd();

//         TDirectory* detDir = asimovDir->GetDirectory(ND_or_FD.c_str());
//         if (!detDir) detDir = asimovDir->mkdir(ND_or_FD.c_str());
//         detDir->cd();

//         pdf->Reweight();
//         MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, -1);

//         pdf->Reweight();
//         MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, -1);

//         pdf->Reweight();
//         MakeSpectaVariations(pdf, "Enubias", fOut, ND_or_FD, pdfTitle, -1);
 
//         fOut->cd();
//         pdf->Reweight();
//         MakeSpectaVariations2D(pdf, fOut, ND_or_FD, pdfTitle, -1, nullptr);  // p=-1 → Asimov path

        

//     }


//     // --- ±shift systematic variations
//     std::vector<double> error(xsec_nominal.size(), shiftAmount);
//     std::vector<double> xsec_plus(xsec_nominal.size()), xsec_minus(xsec_nominal.size());
//     for (size_t i = 0; i < xsec_nominal.size(); ++i) {
//         xsec_plus[i] = xsec_nominal[i] + error[i];
//         xsec_minus[i] = xsec_nominal[i] - error[i];
//     }

//     xsec->SetGroupOnlyParameters("Osc", OscPars);
//     for (auto& pdf : DUNEPdfs) {
//         // +shift
//         xsec->SetParameters(xsec_plus);
//         pdf->Reweight();
//         TH1* h_plus = pdf->GetMCHist(1);
//         if (!h_plus) continue;

//         std::string pdfTitle = pdf->GetTitle();
//         std::string ND_or_FD =
//             (pdfTitle.find("ND") != std::string::npos) ? "ND" :
//             (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";

//         TDirectory* shiftDir = fOut->GetDirectory("shift_parameters");
//         if (!shiftDir) shiftDir = fOut->mkdir("shift_parameters");
//         shiftDir->cd();

//         TDirectory* detDir = shiftDir->GetDirectory(ND_or_FD.c_str());
//         if (!detDir) detDir = shiftDir->mkdir(ND_or_FD.c_str());
//         detDir->cd();

//         TH1D* cloneHplus = static_cast<TH1D*>(h_plus->Clone(
//             Form("%s_%s_plus", ND_or_FD.c_str(), pdfTitle.c_str())));
//         cloneHplus->SetDirectory(detDir);
//         cloneHplus->Write();
//         delete cloneHplus;

//         // -shift
//         xsec->SetParameters(xsec_minus);
//         pdf->Reweight();
//         TH1* h_minus = pdf->GetMCHist(1);
//         if (!h_minus) continue;

//         TH1D* cloneHminus = static_cast<TH1D*>(h_minus->Clone(
//             Form("%s_%s_minus", ND_or_FD.c_str(), pdfTitle.c_str())));
//         cloneHminus->SetDirectory(detDir);
//         cloneHminus->Write();
//         delete cloneHminus;

//         fOut->cd();
//     }

//     // --- Posterior predictive draws
//     // --- Posterior predictive draws
//     auto rnd = std::make_unique<TRandom3>(0);
//     const Long64_t maxSampleSteps = 200000;

//     for (int p = 0; p < no_times_sampling_posterior; ++p) {
//         int entry;
//         do {
//             entry = rnd->Integer(std::min<Long64_t>(maxSampleSteps, nEntries));
//             mcmc->GetEntry(entry);
//         } while ((unsigned int)mcmc_step < burn_in);

//         xsec->SetParameters(xsec_tmp);

//         for (auto& pdf : DUNEPdfs) {
//             std::string pdfTitle = pdf->GetTitle();
//             std::string ND_or_FD =
//                 (pdfTitle.find("ND") != std::string::npos) ? "ND" :
//                 (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";


//             // Generate histograms with correct binning for each variable
//             MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
//             MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
//             MakeSpectaVariations(pdf, "Enubias", fOut, ND_or_FD, pdfTitle, p);
//             MakeSpectaVariations2D(pdf,fOut, ND_or_FD, pdfTitle,p, &posteriorStore);
//         }
//         fOut->cd();
//     }
//     std::cout << "\n==============================\n";
// std::cout << "      POSTERIOR ERRORS        \n";
// std::cout << "==============================\n";

// // -----------------------------------------------------------------------
// // Build summary 2D histograms: posterior uncertainty + statistical uncertainty
// // -----------------------------------------------------------------------
// for (auto& [key, grid] : posteriorStore) {

//     std::cout << "\nSample: " << key << "\n";

//     int nX = grid.size();
//     if (nX == 0) continue;
//     int nY = grid[0].size();

//     // ---- Retrieve axis info from the Asimov 2D histogram stored in fOut ----
//     // Key format: "ND_CleanTitle" or "FD_CleanTitle"
//     // Asimov path: "Asimov/<ND_or_FD>/2D/<ND_or_FD>_<CleanTitle>_2D_Asimov"
//     std::string ND_or_FD, cleanTitle;
//     for (const std::string& prefix : {"ND", "FD", "Other"}) {
//         if (key.size() > prefix.size() &&
//             key.substr(0, prefix.size() + 1) == prefix + "_") {
//             ND_or_FD   = prefix;
//             cleanTitle = key.substr(prefix.size() + 1);
//             break;
//         }
//     }
//     if (ND_or_FD.empty()) {
//         std::cerr << "[WARN] Cannot parse key: " << key << " — skipping\n";
//         continue;
//     }

//     // std::string asimovHistPath =
//     //     "Asimov/" + ND_or_FD + "/2D/" +
//     //     ND_or_FD + "_" + cleanTitle + "_2D_Asimov";

//     // TH2D* hRef = dynamic_cast<TH2D*>(fOut->Get(asimovHistPath.c_str()));
//     // if (!hRef) {
//     //     std::cerr << "[WARN] Could not find reference histogram at "
//     //               << asimovHistPath << " — skipping summary hists for "
//     //               << key << "\n";
//     //     continue;
//     // }
//     // REPLACE with explicit directory navigation:
// TH2D* hRef = nullptr;
// {
//     TDirectory* dAsimov = fOut->GetDirectory("Asimov");
//     TDirectory* dDet    = dAsimov ? dAsimov->GetDirectory(ND_or_FD.c_str()) : nullptr;
//     TDirectory* d2D     = dDet    ? dDet->GetDirectory("2D")                : nullptr;
//     if (d2D) {
//         std::string hName = ND_or_FD + "_" + cleanTitle + "_2D_Asimov";
//         hRef = dynamic_cast<TH2D*>(d2D->Get(hName.c_str()));
//         if (!hRef) {
//             std::cerr << "[WARN] Histogram '" << hName
//                       << "' not found in Asimov/" << ND_or_FD << "/2D\n"
//                       << "[DEBUG] Contents of that directory:\n";
//             d2D->ls();
//         }
//     } else {
//         std::cerr << "[WARN] Directory Asimov/" << ND_or_FD << "/2D not found.\n";
//         if (TDirectory* dA = fOut->GetDirectory("Asimov")) dA->ls();
//     }
// }
// if (!hRef) continue;
//     // Write into  <ND_or_FD>/2D/summary/
//     std::string summaryDir = ND_or_FD + "/2D/summary";
//     TDirectory* sDir = fOut->GetDirectory(summaryDir.c_str());
//     if (!sDir) sDir = fOut->mkdir(summaryDir.c_str());
//     sDir->cd();

//     // --- Three summary histograms (clone axes from Asimov ref) ---
//     auto MakeSummaryHist = [&](const std::string& name,
//                                const std::string& zTitle) -> TH2D* {
//         TH2D* h = static_cast<TH2D*>(hRef->Clone(name.c_str()));
//         h->Reset("ICESM");
//         h->SetTitle((name + ";" +
//                      hRef->GetXaxis()->GetTitle() + ";" +
//                      hRef->GetYaxis()->GetTitle() + ";" +
//                      zTitle).c_str());
//         return h;
//     };

//     TH2D* hPostErr = MakeSummaryHist(
//         ND_or_FD + "_" + cleanTitle + "_2D_posteriorErr",
//         "Posterior #sigma");

//     TH2D* hStatErr = MakeSummaryHist(
//         ND_or_FD + "_" + cleanTitle + "_2D_statErr",
//         "Statistical #sigma (#sqrt{mean})");

//     TH2D* hMean = MakeSummaryHist(
//         ND_or_FD + "_" + cleanTitle + "_2D_posteriorMean",
//         "Posterior mean");

//     // --- Fill bins ---
//     for (int ix = 0; ix < nX; ++ix) {
//         for (int iy = 0; iy < nY; ++iy) {

//             auto& values = grid[ix][iy];
//             if (values.empty()) continue;

//             // Posterior mean
//             double mean = 0.0;
//             for (double v : values) mean += v;
//             mean /= static_cast<double>(values.size());

//             // Posterior std-dev
//             double variance = 0.0;
//             for (double v : values)
//                 variance += (v - mean) * (v - mean);
//             variance /= static_cast<double>(values.size());
//             double posteriorErr = std::sqrt(variance);

//             // Statistical uncertainty = sqrt(mean)  [Poisson]
//             double statErr = (mean > 0.0) ? std::sqrt(mean) : 0.0;

//             // ROOT bins are 1-indexed
//             hPostErr->SetBinContent(ix + 1, iy + 1, posteriorErr);
//             hStatErr->SetBinContent(ix + 1, iy + 1, statErr);
//             hMean   ->SetBinContent(ix + 1, iy + 1, mean);

//             std::cout
//                 << "Bin (" << ix+1 << "," << iy+1 << ")"
//                 << " | Mean = "         << mean
//                 << " | PosteriorErr = " << posteriorErr
//                 << " | StatErr = "      << statErr
//                 << "\n";
//         }
//     }

//     hPostErr->SetDirectory(sDir);  hPostErr->Write();
//     hStatErr->SetDirectory(sDir);  hStatErr->Write();
//     hMean   ->SetDirectory(sDir);  hMean   ->Write();

//     delete hPostErr;
//     delete hStatErr;
//     delete hMean;

//     fOut->cd();
// }

//     std::cout << "[INFO] Writing ROOT file: " << OutFileName << std::endl;
//     fOut->Write();
//     fOut->Close();

//     delete xsec;
//     for (auto sample : DUNEPdfs) delete sample;

//     return 0;
// }


#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
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

// ---------------------------------------------------------------------------
// Helper: get or create a single directory level under parent
// ---------------------------------------------------------------------------
static TDirectory* GetOrMkdir(TDirectory* parent, const std::string& name) {
    TDirectory* d = parent->GetDirectory(name.c_str());
    if (!d) d = parent->mkdir(name.c_str());
    return d;
}

void MakeSpectaVariations(SampleHandlerFD* pdf, const std::string& var,
                          TFile* fout, const std::string& ND_or_FD,
                          const std::string& pdfTitle, int p) {
    // Refresh internal event weights
    (void) pdf->GetMCHist(1);  // calls fill1DHist() internally
    pdf->Reweight();
    TH1D* h = dynamic_cast<TH1D*>(pdf->Get1DVarHist(var.c_str()));
    if (!h) {
        std::cerr << "[WARN] Could not get TH1D for " << var << " from " << pdfTitle << std::endl;
        return;
    }

    std::string cleanTitle = pdfTitle;
    std::replace(cleanTitle.begin(), cleanTitle.end(), ' ', '_');
    std::replace(cleanTitle.begin(), cleanTitle.end(), '/', '_');

    std::string histName;
    TDirectory* dir = nullptr;

    if (p < 0) {
        // Asimov / ND_or_FD / var  — build level by level
        TDirectory* d1 = GetOrMkdir(fout,  "Asimov");
        TDirectory* d2 = GetOrMkdir(d1,     ND_or_FD);
        dir             = GetOrMkdir(d2,     var);
        histName = Form("%s_%s_%s_Asimov",
                        ND_or_FD.c_str(), cleanTitle.c_str(), var.c_str());
    } else {
        // ND_or_FD / var
        TDirectory* d1 = GetOrMkdir(fout,  ND_or_FD);
        dir             = GetOrMkdir(d1,    var);
        histName = Form("%s_%s_%s_posterior_toy_%03d",
                        ND_or_FD.c_str(), cleanTitle.c_str(), var.c_str(), p);
    }

    dir->cd();

    TH1D* cloneH = static_cast<TH1D*>(h->Clone(histName.c_str()));
    cloneH->SetDirectory(dir);
    cloneH->Write();

    delete cloneH;
    fout->cd();
}


void MakeSpectaVariations2D(SampleHandlerFD* pdf,
                            TFile* fout,
                            const std::string& ND_or_FD,
                            const std::string& pdfTitle,
                            int p,
                            std::map<std::string,
                                std::vector<std::vector<std::vector<double>>>>* posteriorStore = nullptr)
{
    // Refresh internal event weights
    (void) pdf->GetMCHist(2);   // ensure 2D hist is built
    pdf->Reweight();

    TH2D* h = dynamic_cast<TH2D*>(pdf->GetMCHist(2));
    if (!h) {
        std::cerr << "[WARN] Could not get TH2D from "
                  << pdfTitle << std::endl;
        return;
    }

    std::string cleanTitle = pdfTitle;
    std::replace(cleanTitle.begin(), cleanTitle.end(), ' ', '_');
    std::replace(cleanTitle.begin(), cleanTitle.end(), '/', '_');

    std::string histName;
    TDirectory* dir = nullptr;

    if (p < 0) {
        // Asimov / ND_or_FD / 2D  — build level by level
        TDirectory* d1 = GetOrMkdir(fout,  "Asimov");
        TDirectory* d2 = GetOrMkdir(d1,     ND_or_FD);
        dir             = GetOrMkdir(d2,     "2D");
        histName = Form("%s_%s_2D_Asimov",
                        ND_or_FD.c_str(), cleanTitle.c_str());
    } else {
        // ND_or_FD / 2D
        TDirectory* d1 = GetOrMkdir(fout,  ND_or_FD);
        dir             = GetOrMkdir(d1,    "2D");
        histName = Form("%s_%s_2D_posterior_toy_%03d",
                        ND_or_FD.c_str(), cleanTitle.c_str(), p);
    }

    dir->cd();

    TH2D* cloneH = static_cast<TH2D*>(h->Clone(histName.c_str()));

    if (p >= 0 && posteriorStore) {
        std::string key = ND_or_FD + "_" + cleanTitle;

        int nX = cloneH->GetNbinsX();
        int nY = cloneH->GetNbinsY();

        if ((*posteriorStore)[key].empty()) {
            (*posteriorStore)[key].resize(nX);
            for (int ix = 0; ix < nX; ++ix)
                (*posteriorStore)[key][ix].resize(nY);
        }

        for (int ix = 1; ix <= nX; ++ix)
            for (int iy = 1; iy <= nY; ++iy)
                (*posteriorStore)[key][ix-1][iy-1].push_back(cloneH->GetBinContent(ix, iy));
    }

    cloneH->SetDirectory(dir);
    cloneH->Write();

    // Print bin info
    for (int ix = 1; ix <= cloneH->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= cloneH->GetNbinsY(); ++iy) {
            std::cout
                << histName
                << " | Bin (" << ix << "," << iy << ")"
                << " | Rate = "    << cloneH->GetBinContent(ix, iy)
                << " | StatErr = " << cloneH->GetBinError(ix, iy)
                << std::endl;
        }
    }

    delete cloneH;
    fout->cd();
}


int main(int argc, char* argv[]) {
    if (argc == 1) {
        std::cout << "Usage: bin/ config.cfg" << std::endl;
        return 1;
    }

    // --- Manager setup
    auto fitMan = std::make_unique<manager>(argv[1]);
    auto PosteriorFile = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFiles"], __FILE__, __LINE__);
    auto burn_in = Get<unsigned int>(fitMan->raw()["General"]["MCMC"]["BurnInSteps"], __FILE__, __LINE__);
    int no_times_sampling_posterior = Get<int>(fitMan->raw()["Predictive"]["SamplePosterior"], __FILE__, __LINE__);
    std::vector<std::string> xsecCovMatrixFile =
        GetFromManager<std::vector<std::string>>(fitMan->raw()["General"]["Systematics"]["XsecCovFile"], {});
    auto OscCovFile = GetFromManager<std::vector<std::string>>(fitMan->raw()["General"]["Systematics"]["OscCovFile"], {});
    auto OscCovName = GetFromManager<std::string>(fitMan->raw()["General"]["Systematics"]["OscCovName"], "osc_cov");
    auto OscPars = GetFromManager<std::vector<double>>(fitMan->raw()["General"]["OscillationParameters"], {});
    double shiftAmount = GetFromManager<double>(fitMan->raw()["Predictive"]["ShiftAmount"], 0.2);

    std::cout << "[INFO] Using shift amount ±" << shiftAmount << std::endl;

    // --- Setup output file
    std::filesystem::path dir = std::filesystem::path(PosteriorFile).parent_path();
    std::string OutFileName = (dir / "Posteriorpredictive_out.root").string();

    std::string suffix = "_posteriorpredictive_withtoys";
    size_t dotPos = OutFileName.find_last_of('.');
    if (dotPos != std::string::npos)
        OutFileName.insert(dotPos, suffix);
    else
        OutFileName += suffix;

    std::cout << "[INFO] Output file will be: " << OutFileName << std::endl;

    TFile* fOut = new TFile(OutFileName.c_str(), "RECREATE");
    MCMCProcessor Processor(PosteriorFile);
    Processor.Initialise();

    // --- Covariances
    ParameterHandlerGeneric* xsec = nullptr;

    // --- Build PDF instances
    std::vector<SampleHandlerFD*> DUNEPdfs;
    MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

    // key → [ix][iy][values over throws]
    std::map<std::string,
        std::vector<std::vector<std::vector<double>>>> posteriorStore;

    // --- Load posteriors
    TChain* mcmc = new TChain("posteriors");
    mcmc->Add(PosteriorFile.c_str());
    Long64_t nEntries = mcmc->GetEntries();
    if (nEntries == 0) {
        std::cerr << "[ERROR] MCMC tree has no entries!" << std::endl;
        return 1;
    }

    // Bind xsec branches
    std::vector<double> xsec_nominal = xsec->GetPreFitValues();
    std::vector<Double_t> xsec_tmp(xsec_nominal.size(), 0.0);
    for (size_t i = 0; i < xsec_nominal.size(); ++i) {
        TString bname = Form("xsec_%zu", i);
        if (mcmc->GetBranch(bname))
            mcmc->SetBranchAddress(bname, &xsec_tmp[i]);
    }

    UInt_t mcmc_step = -1;
    mcmc->SetBranchAddress("step", &mcmc_step);

    // --- Generate Asimov spectra
    xsec->SetGroupOnlyParameters("Xsec", xsec_nominal);
    xsec->SetGroupOnlyParameters("Osc", OscPars);

    for (auto& pdf : DUNEPdfs) {
        std::string pdfTitle = pdf->GetTitle();
        if (pdfTitle.empty()) pdfTitle = "FHC_numu_asimov";

        std::string ND_or_FD =
            (pdfTitle.find("ND") != std::string::npos) ? "ND" :
            (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";

        // All directory creation is handled inside MakeSpectra* via GetOrMkdir
        pdf->Reweight();
        MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, -1);

        pdf->Reweight();
        MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, -1);

        pdf->Reweight();
        MakeSpectaVariations(pdf, "Enubias", fOut, ND_or_FD, pdfTitle, -1);

        pdf->Reweight();
        MakeSpectaVariations2D(pdf, fOut, ND_or_FD, pdfTitle, -1, nullptr);

        fOut->cd();
    }

    // --- ±shift systematic variations
    std::vector<double> error(xsec_nominal.size(), shiftAmount);
    std::vector<double> xsec_plus(xsec_nominal.size()), xsec_minus(xsec_nominal.size());
    for (size_t i = 0; i < xsec_nominal.size(); ++i) {
        xsec_plus[i]  = xsec_nominal[i] + error[i];
        xsec_minus[i] = xsec_nominal[i] - error[i];
    }

    xsec->SetGroupOnlyParameters("Osc", OscPars);
    for (auto& pdf : DUNEPdfs) {
        std::string pdfTitle = pdf->GetTitle();
        std::string ND_or_FD =
            (pdfTitle.find("ND") != std::string::npos) ? "ND" :
            (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";

        // +shift
        xsec->SetParameters(xsec_plus);
        pdf->Reweight();
        TH1* h_plus = pdf->GetMCHist(1);
        if (!h_plus) continue;

        TDirectory* shiftDir = GetOrMkdir(fOut, "shift_parameters");
        TDirectory* detDir   = GetOrMkdir(shiftDir, ND_or_FD);
        detDir->cd();

        TH1D* cloneHplus = static_cast<TH1D*>(h_plus->Clone(
            Form("%s_%s_plus", ND_or_FD.c_str(), pdfTitle.c_str())));
        cloneHplus->SetDirectory(detDir);
        cloneHplus->Write();
        delete cloneHplus;

        // -shift
        xsec->SetParameters(xsec_minus);
        pdf->Reweight();
        TH1* h_minus = pdf->GetMCHist(1);
        if (!h_minus) continue;

        TH1D* cloneHminus = static_cast<TH1D*>(h_minus->Clone(
            Form("%s_%s_minus", ND_or_FD.c_str(), pdfTitle.c_str())));
        cloneHminus->SetDirectory(detDir);
        cloneHminus->Write();
        delete cloneHminus;

        fOut->cd();
    }

    // --- Posterior predictive draws
    auto rnd = std::make_unique<TRandom3>(0);
    const Long64_t maxSampleSteps = 200000;

    for (int p = 0; p < no_times_sampling_posterior; ++p) {
        int entry;
        do {
            entry = rnd->Integer(std::min<Long64_t>(maxSampleSteps, nEntries));
            mcmc->GetEntry(entry);
        } while ((unsigned int)mcmc_step < burn_in);

        xsec->SetParameters(xsec_tmp);

        for (auto& pdf : DUNEPdfs) {
            std::string pdfTitle = pdf->GetTitle();
            std::string ND_or_FD =
                (pdfTitle.find("ND") != std::string::npos) ? "ND" :
                (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";

            MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
            MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
            MakeSpectaVariations(pdf, "Enubias",            fOut, ND_or_FD, pdfTitle, p);
            MakeSpectaVariations2D(pdf, fOut, ND_or_FD, pdfTitle, p, &posteriorStore);
        }
        fOut->cd();
    }

    // -----------------------------------------------------------------------
    // Build summary 2D histograms: posterior uncertainty + statistical uncertainty
    // -----------------------------------------------------------------------
    std::cout << "\n==============================\n";
    std::cout << "      POSTERIOR ERRORS        \n";
    std::cout << "==============================\n";

    for (auto& [key, grid] : posteriorStore) {

        std::cout << "\nSample: " << key << "\n";

        int nX = grid.size();
        if (nX == 0) continue;
        int nY = grid[0].size();

        // Parse key: prefix is one of "ND", "FD", "Other"
        std::string ND_or_FD, cleanTitle;
        for (const std::string& prefix : {"ND", "FD", "Other"}) {
            if (key.size() > prefix.size() &&
                key.substr(0, prefix.size() + 1) == prefix + "_") {
                ND_or_FD   = prefix;
                cleanTitle = key.substr(prefix.size() + 1);
                break;
            }
        }
        if (ND_or_FD.empty()) {
            std::cerr << "[WARN] Cannot parse key: " << key << " — skipping\n";
            continue;
        }

        // Navigate to Asimov/ND_or_FD/2D level by level to find the reference histogram
        TH2D* hRef = nullptr;
        {
            TDirectory* dAsimov = fOut->GetDirectory("Asimov");
            TDirectory* dDet    = dAsimov ? dAsimov->GetDirectory(ND_or_FD.c_str()) : nullptr;
            TDirectory* d2D     = dDet    ? dDet->GetDirectory("2D")                : nullptr;
            if (d2D) {
                std::string hName = ND_or_FD + "_" + cleanTitle + "_2D_Asimov";
                hRef = dynamic_cast<TH2D*>(d2D->Get(hName.c_str()));
                if (!hRef) {
                    std::cerr << "[WARN] Histogram '" << hName
                              << "' not found in Asimov/" << ND_or_FD << "/2D\n"
                              << "[DEBUG] Contents:\n";
                    d2D->ls();
                }
            } else {
                std::cerr << "[WARN] Directory Asimov/" << ND_or_FD << "/2D not found.\n";
                if (TDirectory* dA = fOut->GetDirectory("Asimov")) dA->ls();
            }
        }
        if (!hRef) continue;

        // Build summary directory: ND_or_FD / 2D / summary
        TDirectory* sDir = GetOrMkdir(GetOrMkdir(GetOrMkdir(fOut, ND_or_FD), "2D"), "summary");
        sDir->cd();

        // Lambda to clone axis structure from hRef
        auto MakeSummaryHist = [&](const std::string& name,
                                   const std::string& zTitle) -> TH2D* {
            TH2D* h = static_cast<TH2D*>(hRef->Clone(name.c_str()));
            h->Reset("ICESM");
            h->SetTitle((name + ";" +
                         hRef->GetXaxis()->GetTitle() + ";" +
                         hRef->GetYaxis()->GetTitle() + ";" +
                         zTitle).c_str());
            return h;
        };

        TH2D* hPostErr = MakeSummaryHist(ND_or_FD + "_" + cleanTitle + "_2D_posteriorErr",
                                          "Posterior #sigma");
        TH2D* hStatErr = MakeSummaryHist(ND_or_FD + "_" + cleanTitle + "_2D_statErr",
                                          "Statistical #sigma (#sqrt{mean})");
        TH2D* hMean    = MakeSummaryHist(ND_or_FD + "_" + cleanTitle + "_2D_posteriorMean",
                                          "Posterior mean");

        for (int ix = 0; ix < nX; ++ix) {
            for (int iy = 0; iy < nY; ++iy) {

                auto& values = grid[ix][iy];
                if (values.empty()) continue;

                double mean = 0.0;
                for (double v : values) mean += v;
                mean /= static_cast<double>(values.size());

                double variance = 0.0;
                for (double v : values)
                    variance += (v - mean) * (v - mean);
                variance /= static_cast<double>(values.size());
                double posteriorErr = std::sqrt(variance);

                double statErr = (mean > 0.0) ? std::sqrt(mean) : 0.0;

                hPostErr->SetBinContent(ix + 1, iy + 1, posteriorErr);
                hStatErr->SetBinContent(ix + 1, iy + 1, statErr);
                hMean   ->SetBinContent(ix + 1, iy + 1, mean);

                std::cout
                    << "Bin (" << ix+1 << "," << iy+1 << ")"
                    << " | Mean = "         << mean
                    << " | PosteriorErr = " << posteriorErr
                    << " | StatErr = "      << statErr
                    << "\n";
            }
        }

        hPostErr->SetDirectory(sDir);  hPostErr->Write();
        hStatErr->SetDirectory(sDir);  hStatErr->Write();
        hMean   ->SetDirectory(sDir);  hMean   ->Write();

        delete hPostErr;
        delete hStatErr;
        delete hMean;

        fOut->cd();
    }

    std::cout << "[INFO] Writing ROOT file: " << OutFileName << std::endl;
    fOut->Write();
    fOut->Close();

    delete xsec;
    for (auto sample : DUNEPdfs) delete sample;

    return 0;
}