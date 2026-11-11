// #include "Samples/BinningHandler.h"
// #include "Samples/HistogramUtils.h"
// #include "Samples/MaCh3DUNEFactory.h"
// #include "Manager/MaCh3Exception.h"
// #include "Manager/MaCh3Logger.h"
// #include "Splines/SplineMonolith.h"
// #include "Samples/SampleHandlerBase.h"
// #include "Samples/SampleHandlerFD.h"
// #include "Fitters/MaCh3Factory.h"
// #include "Fitters/FitterBase.h"
// #include "Manager/Manager.h"
// #include "Parameters/ParameterHandlerBase.h"
// #include "Parameters/ParameterStructs.h"

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
//     //TH1D* h = dynamic_cast<TH1D*>(pdf->Get1DVarHist(var.c_str()));
//     TH1* h = pdf->Get1DVarHist(
//     0,          // iSample (pick first sample for now)
//     var,        // std::string (NOT var.c_str())
//     {},         // no event cuts
//     0,          // default weight style
//     nullptr,    // default axis
//     {}          // no sub-event cuts
// );
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
//     delete h;  // cleanup to avoid memory leak
//     fout->cd();
// }



// int main(int argc, char* argv[]) {
//     if (argc == 1) {
//         std::cout << "Usage: bin/ config.cfg" << std::endl;
//         return 1;
//     }

//     // --- Manager setup
//     auto fitMan = std::make_unique<Manager>(argv[1]);
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

//     xsec->SetGroupOnlyParameters("Osc", OscPars);

//     // --- Load posteriors
//     TChain* mcmc = new TChain("posteriors");
//     mcmc->Add(PosteriorFile.c_str());
//     Long64_t nEntries = mcmc->GetEntries();
//     if (nEntries == 0) {
//         std::cerr << "[ERROR] MCMC tree has no entries!" << std::endl;
//         return 1;
//     }

//     // Bind xsec branches
//     //std::vector<double> xsec_nominal = xsec->GetPreFitValues();
//     // std::vector<double> xsec_nominal( xsec->GetNumParFromGroup("Xsec"));
//     // std::vector<Double_t> xsec_tmp(xsec_nominal.size(), 0.0);
//     // FIXED CODE - Use exactly 480 parameters
//     // Get the actual number of parameters from the parameter handler
//     //int nXsecPars = xsec->GetNumParFromGroup("Xsec");
//     //std::cout << "[DEBUG] Number of Xsec parameters: " << nXsecPars << std::endl;

//     // Initialize with correct size
//     //std::vector<double> xsec_nominal(nXsecPars);
//     //xsec_nominal = xsec->GetPreFitValues(); // Get actual pre-fit values
//     const int nXsecPars = xsec->GetNumParFromGroup("Xsec");

//     std::vector<double> xsec_nominal(nXsecPars);
//     xsec->SetGroupOnlyParameters("Xsec", xsec_nominal);
//     xsec->SetGroupOnlyParameters("Osc", OscPars);

//     std::vector<Double_t> xsec_tmp(nXsecPars, 0.0);

//     // Debug: Check available branches
// std::cout << "[DEBUG] Available MCMC branches:" << std::endl;
// TObjArray* branches = mcmc->GetListOfBranches();
// for (int i = 0; i < branches->GetEntries(); ++i) {
//     TBranch* branch = static_cast<TBranch*>(branches->At(i));
//     if (std::string(branch->GetName()).find("xsec") != std::string::npos) {
//         std::cout << "  " << branch->GetName() << std::endl;
//     }
// }

// // Check how many xsec branches we actually have
// int nXsecBranches = 0;
// for (size_t i = 0; i < 500; ++i) { // Check up to 500
//     TString bname = Form("xsec_%zu", i);
//     if (mcmc->GetBranch(bname)) {
//         nXsecBranches++;
//     } else {
//         break;
//     }
// }
// std::cout << "[DEBUG] Found " << nXsecBranches << " xsec branches in MCMC" << std::endl;

// // Get pre-fit values but only take first 480
// auto prefitValues = xsec->GetPreFitValues();
// for (int i = 0; i < nXsecBranches && i < prefitValues.size(); ++i) {
//     xsec_nominal[i] = prefitValues[i];
// }

//     for (size_t i = 0; i < xsec_nominal.size(); ++i) {
//         TString bname = Form("xsec_%zu", i);
//         if (mcmc->GetBranch(bname))
//             mcmc->SetBranchAddress(bname, &xsec_tmp[i]);
//     }

//     UInt_t mcmc_step = -1;
//     mcmc->SetBranchAddress("step", &mcmc_step);

//     // --- Generate Asimov spectra
//     xsec->SetGroupOnlyParameters("Xsec", xsec_nominal);
//     xsec->SetGroupOnlyParameters("Osc", OscPars);

//     for (auto& pdf : DUNEPdfs) {
//         std::string pdfTitle = pdf->GetName();
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

//         // --- True Neutrino Energy ---
//         pdf->Reweight();  // refresh weights for Asimov
//         MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, -1);

//         // --- Reco Neutrino Energy ---
//         pdf->Reweight();  // refresh weights for Reco
//         MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, -1);

//         pdf->Reweight();  // refresh weights for Reco
//         MakeSpectaVariations(pdf, "Enubias", fOut, ND_or_FD, pdfTitle, -1);


//         fOut->cd();
//     }


//     // // --- ±shift systematic variations
//     // std::vector<double> error(xsec_nominal.size(), shiftAmount);
//     // std::vector<double> xsec_plus(xsec_nominal.size()), xsec_minus(xsec_nominal.size());
//     // for (size_t i = 0; i < xsec_nominal.size(); ++i) {
//     //     xsec_plus[i] = xsec_nominal[i] + error[i];
//     //     xsec_minus[i] = xsec_nominal[i] - error[i];
//     // }

//     //xsec->SetGroupOnlyParameters("Osc", OscPars);
//     //  xsec->SetGroupOnlyParameters("Xsec", xsec_nominal);
//     // for (auto& pdf : DUNEPdfs) {
//     //     // +shift
//     //     xsec->SetParameters(xsec_plus);
//     //     pdf->Reweight();
//     //     TH1* h_plus = pdf->GetMCHist(1);
//     //     if (!h_plus) continue;

//     //     std::string pdfTitle = pdf->GetName();
//     //     std::string ND_or_FD =
//     //         (pdfTitle.find("ND") != std::string::npos) ? "ND" :
//     //         (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";

//     //     TDirectory* shiftDir = fOut->GetDirectory("shift_parameters");
//     //     if (!shiftDir) shiftDir = fOut->mkdir("shift_parameters");
//     //     shiftDir->cd();

//     //     TDirectory* detDir = shiftDir->GetDirectory(ND_or_FD.c_str());
//     //     if (!detDir) detDir = shiftDir->mkdir(ND_or_FD.c_str());
//     //     detDir->cd();

//     //     TH1D* cloneHplus = static_cast<TH1D*>(h_plus->Clone(
//     //         Form("%s_%s_plus", ND_or_FD.c_str(), pdfTitle.c_str())));
//     //     cloneHplus->SetDirectory(detDir);
//     //     cloneHplus->Write();
//     //     delete cloneHplus;

//     //     // -shift
//     //     xsec->SetParameters(xsec_minus);
//     //     pdf->Reweight();
//     //     TH1* h_minus = pdf->GetMCHist(1);
//     //     if (!h_minus) continue;

//     //     TH1D* cloneHminus = static_cast<TH1D*>(h_minus->Clone(
//     //         Form("%s_%s_minus", ND_or_FD.c_str(), pdfTitle.c_str())));
//     //     cloneHminus->SetDirectory(detDir);
//     //     cloneHminus->Write();
//     //     delete cloneHminus;

//     //     fOut->cd();
//     // }

//     // --- Posterior predictive draws
//     // auto rnd = std::make_unique<TRandom3>(0);
//     // const Long64_t maxSampleSteps = 2000000;

//     // for (int p = 0; p < no_times_sampling_posterior; ++p) {
//     //     int entry;
//     //     do {
//     //         entry = rnd->Integer(std::min<Long64_t>(maxSampleSteps, nEntries));
//     //         mcmc->GetEntry(entry);
//     //     } while ((unsigned int)mcmc_step < burn_in);

//     //     xsec->SetGroupOnlyParameters("Xsec", xsec_tmp);

//     //     for (auto& pdf : DUNEPdfs) {
//     //         pdf->Reweight();
//     //         TH1* h_after = pdf->GetMCHist(1);
//     //         if (!h_after) continue;

//     //         std::string pdfTitle = pdf->GetName();
//     //         std::string ND_or_FD =
//     //             (pdfTitle.find("ND") != std::string::npos) ? "ND" :
//     //             (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";

//     //         TDirectory* detDir = fOut->GetDirectory(ND_or_FD.c_str());
//     //         if (!detDir) detDir = fOut->mkdir(ND_or_FD.c_str());
//     //         detDir->cd();

//     //         TH1D* cloneH = static_cast<TH1D*>(h_after->Clone(
//     //             Form("%s_%s_posterior_toy_%03d", ND_or_FD.c_str(), pdfTitle.c_str(), p)));
//     //         cloneH->SetDirectory(detDir);
//     //         cloneH->Write();
//     //         delete cloneH;

//     //         MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
//     //         MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
//     //     }
//     //     fOut->cd();
//     // }
//     // --- Posterior predictive draws
//     auto rnd = std::make_unique<TRandom3>(0);
//     const Long64_t maxSampleSteps = 2000000;

//     for (int p = 0; p < no_times_sampling_posterior; ++p) {
//         int entry;
//         do {
//             entry = rnd->Integer(std::min<Long64_t>(maxSampleSteps, nEntries));
//             mcmc->GetEntry(entry);
//         } while ((unsigned int)mcmc_step < burn_in);

//        xsec->SetGroupOnlyParameters("Xsec", xsec_tmp);

//         for (auto& pdf : DUNEPdfs) {
//             std::string pdfTitle = pdf->GetName();
//             std::string ND_or_FD =
//                 (pdfTitle.find("ND") != std::string::npos) ? "ND" :
//                 (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";

//             // Generate histograms with correct binning for each variable
//             MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
//             MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
//             MakeSpectaVariations(pdf, "Enubias", fOut, ND_or_FD, pdfTitle, p);
//         }
//         fOut->cd();
//     }

//     std::cout << "[INFO] Writing ROOT file: " << OutFileName << std::endl;
//     fOut->Write();
//     fOut->Close();

//     delete xsec;
//     for (auto sample : DUNEPdfs) delete sample;

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

// ---------------------------------------------------------------------------
// Determine whether a sample title corresponds to ND or FD.
// Matches on current naming conventions:
//   - ND samples contain "OffAxis" (e.g. "OffAxis0m_numuCC_numode")
//   - FD samples contain "FD"      (e.g. "BeamFD_FHC_numu")
// Falls back to FD with a loud warning if neither pattern matches, so a
// misclassification is visible in the logs rather than silent.
// ---------------------------------------------------------------------------
std::string DetermineNDorFD(const std::string& pdfTitle)
{
    if (pdfTitle.find("OffAxis") != std::string::npos) {
        return "ND";
    } else if (pdfTitle.find("FD") != std::string::npos) {
        return "FD";
    } else {
        std::cerr << "[WARN] Could not determine ND/FD for sample title: \""
                   << pdfTitle << "\" — defaulting to FD" << std::endl;
        return "FD";
    }
}

// ---------------------------------------------------------------------------
// Write a 1D histogram for a given kinematic variable into fout
// p < -1 => FakeData Asimov, p < 0 => Asimov, p >= 0 => posterior toy index p
// ---------------------------------------------------------------------------
void MakeSpectaVariations(SampleHandlerFD*   pdf,
                          const std::string& var,
                          TFile*             fout,
                          const std::string& ND_or_FD,
                          const std::string& pdfTitle,
                          int                p,
                          int                sample_idx = 0)
{
    std::vector<KinematicCut> emptyEventCuts;
    std::vector<KinematicCut> emptySubEventCuts;

    TH1* h = pdf->Get1DVarHist(sample_idx, var,
                               emptyEventCuts, 0, nullptr,
                               emptySubEventCuts);
    if (!h) {
        std::cerr << "[WARN] Could not get histogram for " << var
                  << " from " << pdfTitle << std::endl;
        return;
    }

    TH1D* h1d = dynamic_cast<TH1D*>(h);
    if (!h1d) {
        std::cerr << "[ERROR] Could not cast to TH1D for " << var << std::endl;
        delete h;
        return;
    }

    // Build clean directory / histogram names
    std::string cleanTitle = pdfTitle;
    std::replace(cleanTitle.begin(), cleanTitle.end(), ' ', '_');
    std::replace(cleanTitle.begin(), cleanTitle.end(), '/', '_');

    std::string baseDir, histName;
    if (p < -1) {
        baseDir  = "AsimovFakeData/" + ND_or_FD + "/" + var;
        histName = Form("%s_%s_%s_AsimovFakeData",
                    ND_or_FD.c_str(), cleanTitle.c_str(), var.c_str());
    } else if (p < 0) {
        baseDir  = "Asimov/" + ND_or_FD + "/" + var;
        histName = Form("%s_%s_%s_Asimov",
                        ND_or_FD.c_str(), cleanTitle.c_str(), var.c_str());
    } else {
        baseDir  = ND_or_FD + "/" + var;
        histName = Form("%s_%s_%s_posterior_toy_%03d",
                        ND_or_FD.c_str(), cleanTitle.c_str(), var.c_str(), p);
    }

    TDirectory* dir = fout->GetDirectory(baseDir.c_str());
    if (!dir) dir = fout->mkdir(baseDir.c_str());
    dir->cd();

    TH1D* clone = static_cast<TH1D*>(h1d->Clone(histName.c_str()));
    clone->SetDirectory(dir);
    clone->Write();

    delete clone;
    delete h;
    fout->cd();
}


// ---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    if (argc == 1) {
        std::cout << "Usage: bin/ config.cfg" << std::endl;
        return 1;
    }

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
    std::string OutFileName = (outDir / "Posteriorpredictive_out_posteriorpredictive_withtoys.root").string();
    std::cout << "[INFO] Output file: " << OutFileName << std::endl;

    TFile* fOut = new TFile(OutFileName.c_str(), "RECREATE");

    // -----------------------------------------------------------------------
    // Build PDFs and parameter handler
    // -----------------------------------------------------------------------
    ParameterHandlerGeneric* xsec = nullptr;
    std::vector<SampleHandlerFD*> DUNEPdfs;
    MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

    auto ExtractGroup = [&](const std::vector<double>& full, const std::string& group) {
        std::vector<double> out;
        out.reserve(xsec->GetNumParFromGroup(group));
        for (int i = 0; i < (int)full.size(); ++i) {
            if (xsec->IsParFromGroup(i, group)) out.push_back(full[i]);
        }
        return out;
    };

    // -----------------------------------------------------------------------
    // IMPORTANT: use the TOTAL number of parameters in the handler, not just
    // Xsec+Flux. The global parameter indexing used by IsParFromGroup() spans
    // every group (Xsec, Flux, Osc, FakeData, ...), so any buffer we read the
    // MCMC branches into, and any vector we pass to ExtractGroup, must cover
    // the full index range or later groups get silently truncated (this is
    // what caused the "Number of param in group Xsec is 600, while you
    // passed 592" crash).
    // -----------------------------------------------------------------------
    const int nXsecOnly  = xsec->GetNumParFromGroup("Xsec");
    const int nFluxOnly  = xsec->GetNumParFromGroup("Flux");
    const int nOscOnly   = xsec->GetNumParFromGroup("Osc");
    int nFakeDataPars    = xsec->GetNumParFromGroup("FakeData");

    // If your ParameterHandler exposes a direct "total number of parameters"
    // accessor (e.g. GetNumParams()), prefer that instead of summing groups,
    // since summing only known groups can undercount if there are others.
    // const int nTotalPars = xsec->GetNumParams();
    const int nTotalPars = nXsecOnly + nFluxOnly + nOscOnly + nFakeDataPars;

    std::cout << "[INFO] nXsec=" << nXsecOnly
              << " nFlux=" << nFluxOnly
              << " nOsc="  << nOscOnly
              << " nFakeData=" << nFakeDataPars
              << " nTotalPars=" << nTotalPars << std::endl;

    auto prefitValues = xsec->GetPreFitValues();  // still needed later for FakeData extraction

    // Reset Xsec and Flux to their prior values directly — no vector needed
    xsec->SetGroupOnlyParameters("Xsec");
    xsec->SetGroupOnlyParameters("Flux");
    xsec->SetGroupOnlyParameters("Osc",  OscPars);

    std::vector<double> fakeDataZero(nFakeDataPars, 0.0);
    std::vector<double> fakeDataPars(nFakeDataPars);
    xsec->SetGroupOnlyParameters("FakeData", fakeDataZero);

    // -----------------------------------------------------------------------
    // Asimov spectra (nominal parameters)
    // -----------------------------------------------------------------------
    std::cout << "[INFO] Generating Asimov spectra..." << std::endl;

    for (auto& pdf : DUNEPdfs) {
        std::string pdfTitle = pdf->GetSampleTitle(0);
        std::cout << "[DEBUG] Raw sample title: \"" << pdfTitle << "\"" << std::endl;
        if (pdfTitle.empty()) pdfTitle = "FHC_numu_asimov";

        std::string ND_or_FD = DetermineNDorFD(pdfTitle);

        pdf->Reweight();
        MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, -1);
        MakeSpectaVariations(pdf, "RecoNeutrinoEnergy",  fOut, ND_or_FD, pdfTitle, -1);
        MakeSpectaVariations(pdf, "Enubias",             fOut, ND_or_FD, pdfTitle, -1);
        fOut->cd();
    }

    int fd_idx = 0;
    for (int i = 0; i < (int)prefitValues.size(); ++i) {
        if (xsec->IsParFromGroup(i, "FakeData")) {
            fakeDataPars[fd_idx++] = prefitValues[i];
        }
    }

    // --- FakeData Asimov ---
    xsec->SetGroupOnlyParameters("FakeData", fakeDataPars);
    int fpIdx = xsec->GetParIndex("MissingProtonFD"); // or whichever dial name
    MACH3LOG_INFO("Prop vec value for MissingProtonFD just before FD Asimov reweight: {}",
                  xsec->GetParPropVec()[fpIdx]);

    for (auto& pdf : DUNEPdfs) {
        std::string pdfTitle = pdf->GetSampleTitle(0);
        if (pdfTitle.empty()) pdfTitle = "FHC_numu_asimov";

        std::string ND_or_FD = DetermineNDorFD(pdfTitle);

        pdf->Reweight();
        MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, -2);
        MakeSpectaVariations(pdf, "RecoNeutrinoEnergy",  fOut, ND_or_FD, pdfTitle, -2);
        MakeSpectaVariations(pdf, "Enubias",             fOut, ND_or_FD, pdfTitle, -2);
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
    // Binary-search for first post-burn-in entry (read only "step" branch)
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
        std::cerr << "[ERROR] No MCMC entries survive burn-in cut of "
                  << burn_in << std::endl;
        return 1;
    }

    Long64_t nValid = nEntries - firstValid;
    std::cout << "[INFO] First post-burn-in entry: " << firstValid
              << "  (" << nValid << " valid entries)" << std::endl;

    // -----------------------------------------------------------------------
    // Discover how many param_N branches actually exist
    // -----------------------------------------------------------------------
    mcmc->SetBranchStatus("*", 1);   // re-enable everything temporarily

    int nBranchesFound = 0;
    for (int i = 0; i < 2000; ++i) {
        if (mcmc->GetBranch(Form("param_%d", i))) ++nBranchesFound;
        else break;
    }
    std::cout << "[INFO] Found " << nBranchesFound << " param branches in MCMC." << std::endl;
    std::cout << "[INFO] Handler expects " << nTotalPars << " total parameters." << std::endl;
    if (nBranchesFound != nTotalPars) {
        std::cerr << "[WARN] Mismatch between MCMC branches (" << nBranchesFound
                   << ") and handler parameter count (" << nTotalPars
                   << "). Check that PosteriorFile was produced with the same "
                   << "parameter configuration as this handler, and that all "
                   << "relevant groups (Xsec, Flux, Osc, FakeData, ...) are "
                   << "accounted for above." << std::endl;
    }

    // -----------------------------------------------------------------------
    // Enable only the branches we actually need, then bind addresses.
    // Buffer is now sized to nTotalPars so Osc (and any other trailing
    // group) is actually read from the chain instead of being silently cut.
    // -----------------------------------------------------------------------
    const int nParsToRead = std::min(nBranchesFound, nTotalPars);

    mcmc->SetBranchStatus("*",    0);
    mcmc->SetBranchStatus("step", 1);
    mcmc->SetBranchAddress("step", &mcmc_step);

    std::vector<Double_t> xsec_tmp(nTotalPars, 0.0);
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

    std::sort(chosenEntries.begin(), chosenEntries.end());  // sequential I/O

    // Read param values for each chosen entry into memory
    std::vector<std::vector<Double_t>> sampledParams(no_times_sampling_posterior,
                                                     std::vector<Double_t>(nTotalPars, 0.0));
    for (int p = 0; p < no_times_sampling_posterior; ++p) {
        mcmc->GetEntry(chosenEntries[p]);
        sampledParams[p].assign(xsec_tmp.begin(), xsec_tmp.end());
    }

    std::cout << "[INFO] Finished reading posterior samples. "
              << "Starting reweighting loop..." << std::endl;

    // -----------------------------------------------------------------------
    // Posterior predictive loop – pure in-memory reweighting, no more I/O
    // Osc is now sampled per-toy from the chain, matching Xsec/Flux, instead
    // of being pinned to the fixed config-level OscPars for every toy.
    // -----------------------------------------------------------------------
    for (int p = 0; p < no_times_sampling_posterior; ++p) {
        if (p % 50 == 0)
            std::cout << "[INFO] Toy " << p << " / "
                      << no_times_sampling_posterior << std::endl;

        xsec->SetGroupOnlyParameters("Xsec", ExtractGroup(sampledParams[p], "Xsec"));
        xsec->SetGroupOnlyParameters("Flux", ExtractGroup(sampledParams[p], "Flux"));
        xsec->SetGroupOnlyParameters("Osc",  ExtractGroup(sampledParams[p], "Osc"));

        for (auto& pdf : DUNEPdfs) {
            pdf->Reweight();

            std::string pdfTitle = pdf->GetSampleTitle(0);
            if (pdfTitle.empty()) pdfTitle = "FHC_numu_asimov";

            std::string ND_or_FD = DetermineNDorFD(pdfTitle);

            MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
            MakeSpectaVariations(pdf, "RecoNeutrinoEnergy",  fOut, ND_or_FD, pdfTitle, p);
            MakeSpectaVariations(pdf, "Enubias",             fOut, ND_or_FD, pdfTitle, p);
        }
        fOut->cd();
    }

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