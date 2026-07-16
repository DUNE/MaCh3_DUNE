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

// // ---------------------------------------------------------------------------
// // Write a 1D histogram for a given kinematic variable into fout
// // p < 0  => Asimov,  p >= 0 => posterior toy index p
// // ---------------------------------------------------------------------------
// void MakeSpectaVariations(SampleHandlerFD*   pdf,
//                           const std::string& var,
//                           TFile*             fout,
//                           const std::string& ND_or_FD,
//                           const std::string& pdfTitle,
//                           int                p,
//                           int                sample_idx = 0)
// {
//     std::vector<KinematicCut> emptyEventCuts;
//     std::vector<KinematicCut> emptySubEventCuts;

//     TH1* h = pdf->Get1DVarHist(sample_idx, var,
//                                emptyEventCuts, 0, nullptr,
//                                emptySubEventCuts);
//     if (!h) {
//         std::cerr << "[WARN] Could not get histogram for " << var
//                   << " from " << pdfTitle << std::endl;
//         return;
//     }

//     TH1D* h1d = dynamic_cast<TH1D*>(h);
//     if (!h1d) {
//         std::cerr << "[ERROR] Could not cast to TH1D for " << var << std::endl;
//         delete h;
//         return;
//     }

//     // Build clean directory / histogram names
//     std::string cleanTitle = pdfTitle;
//     std::replace(cleanTitle.begin(), cleanTitle.end(), ' ', '_');
//     std::replace(cleanTitle.begin(), cleanTitle.end(), '/', '_');

//     std::string baseDir, histName;
//     if (p < -1) {
//     baseDir  = "AsimovFakeData/" + ND_or_FD + "/" + var;
//     histName = Form("%s_%s_%s_AsimovFakeData",
//                 ND_or_FD.c_str(), cleanTitle.c_str(), var.c_str());
//     } else if (p < 0){
//         baseDir  = "Asimov/" + ND_or_FD + "/" + var;
//         histName = Form("%s_%s_%s_Asimov",
//                         ND_or_FD.c_str(), cleanTitle.c_str(), var.c_str());
//     } else {
//         baseDir  = ND_or_FD + "/" + var;
//         histName = Form("%s_%s_%s_posterior_toy_%03d",
//                         ND_or_FD.c_str(), cleanTitle.c_str(), var.c_str(), p);
//     }

//     TDirectory* dir = fout->GetDirectory(baseDir.c_str());
//     if (!dir) dir = fout->mkdir(baseDir.c_str());
//     dir->cd();

//     TH1D* clone = static_cast<TH1D*>(h1d->Clone(histName.c_str()));
//     clone->SetDirectory(dir);
//     clone->Write();

//     delete clone;
//     delete h;
//     fout->cd();
// }


// // ---------------------------------------------------------------------------
// int main(int argc, char* argv[])
// {
//     if (argc == 1) {
//         std::cout << "Usage: bin/ config.cfg" << std::endl;
//         return 1;
//     }

    

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
//     std::string OutFileName = (outDir / "Posteriorpredictive_out_posteriorpredictive_withtoys.root").string();
//     std::cout << "[INFO] Output file: " << OutFileName << std::endl;

//     TFile* fOut = new TFile(OutFileName.c_str(), "RECREATE");

//     // -----------------------------------------------------------------------
//     // Build PDFs and parameter handler
//     // -----------------------------------------------------------------------
//     ParameterHandlerGeneric* xsec = nullptr;
//     std::vector<SampleHandlerFD*> DUNEPdfs;
//     MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

//     // const int nXsecPars = xsec->GetNumParFromGroup("Xsec")
//     //                     + xsec->GetNumParFromGroup("Flux");
                       
//     // // Nominal = pre-fit values
//     // std::vector<double> xsec_nominal(nXsecPars, 0.0);
//     // auto prefitValues = xsec->GetPreFitValues();
//     // for (int i = 0; i < nXsecPars && i < (int)prefitValues.size(); ++i)
//     //     xsec_nominal[i] = prefitValues[i];

//     // xsec->SetGroupOnlyParameters("Xsec", xsec_nominal);
//     // xsec->SetGroupOnlyParameters("Flux",  xsec_nominal);
//     // xsec->SetGroupOnlyParameters("Osc",  OscPars);

//     auto ExtractGroup = [&](const std::vector<double>& full, const std::string& group) {
//     std::vector<double> out;
//     out.reserve(xsec->GetNumParFromGroup(group));
//     for (int i = 0; i < (int)full.size(); ++i) {
//         if (xsec->IsParFromGroup(i, group)) out.push_back(full[i]);
//     }
//     return out;
//     };
//     const int nXsecPars = xsec->GetNumParFromGroup("Xsec")
//                     + xsec->GetNumParFromGroup("Flux");

//     auto prefitValues = xsec->GetPreFitValues();  // still needed later for FakeData extraction

//     // Reset Xsec and Flux to their prior values directly — no vector needed
//     xsec->SetGroupOnlyParameters("Xsec");
//     xsec->SetGroupOnlyParameters("Flux");
//     xsec->SetGroupOnlyParameters("Osc",  OscPars);

//     int nFakeDataPars = xsec->GetNumParFromGroup("FakeData");
//     std::vector<double> fakeDataZero(nFakeDataPars, 0.0);
//     std::vector<double> fakeDataPars(nFakeDataPars);  
//     xsec->SetGroupOnlyParameters("FakeData", fakeDataZero);
//     // -----------------------------------------------------------------------
//     // Asimov spectra (nominal parameters)
    
//     std::cout << "[INFO] Generating Asimov spectra..." << std::endl;

//     for (auto& pdf : DUNEPdfs) {
//         std::string pdfTitle = pdf->GetSampleTitle(0);
//         std::cout << "[DEBUG] Raw sample title: \"" << pdfTitle << "\"" << std::endl;
//         if (pdfTitle.empty()) pdfTitle = "FHC_numu_asimov";

//         std::string ND_or_FD;
//       if (pdfTitle.find("OffAxis") != std::string::npos) {
//             ND_or_FD = "ND";
//         } else if (pdfTitle.find("FD") != std::string::npos) {
//             ND_or_FD = "FD";
//         } else {
//             std::cerr << "[WARN] Could not determine ND/FD for sample title: \""
//                        << pdfTitle << "\" — defaulting to FD" << std::endl;
//             ND_or_FD = "FD";
//         }

//         pdf->Reweight();
//         MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, -1);
//         MakeSpectaVariations(pdf, "RecoNeutrinoEnergy",  fOut, ND_or_FD, pdfTitle, -1);
//         MakeSpectaVariations(pdf, "Enubias",             fOut, ND_or_FD, pdfTitle, -1);
//         fOut->cd();
//     }

   

//     int fd_idx = 0;
//     for (int i = 0; i < (int)prefitValues.size(); ++i) {
//         if (xsec->IsParFromGroup(i, "FakeData")) {
//             fakeDataPars[fd_idx++] = prefitValues[i];
//         }
//     }


//     // --- FakeData Asimov (new) ---
//     xsec->SetGroupOnlyParameters("FakeData", fakeDataPars);
//     for (auto& pdf : DUNEPdfs) {
        
//         std::string pdfTitle = pdf->GetSampleTitle(0);
//         if (pdfTitle.empty()) pdfTitle = "FHC_numu_asimov";
//         +        std::string ND_or_FD;
// +        if (pdfTitle.find("OffAxis") != std::string::npos) {
// +            ND_or_FD = "ND";
// +        } else if (pdfTitle.find("FD") != std::string::npos) {
// +            ND_or_FD = "FD";
// +        } else {
// +            std::cerr << "[WARN] Could not determine ND/FD for sample title: \""
// +                       << pdfTitle << "\" — defaulting to FD" << std::endl;
// +            ND_or_FD = "FD";
// +        }
//         pdf->Reweight();
//         MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, -2);
//         MakeSpectaVariations(pdf, "RecoNeutrinoEnergy",  fOut, ND_or_FD, pdfTitle, -2);
//         MakeSpectaVariations(pdf, "Enubias",             fOut, ND_or_FD, pdfTitle, -2);
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
//     // Binary-search for first post-burn-in entry (read only "step" branch)
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
//         std::cerr << "[ERROR] No MCMC entries survive burn-in cut of "
//                   << burn_in << std::endl;
//         return 1;
//     }

//     Long64_t nValid = nEntries - firstValid;
//     std::cout << "[INFO] First post-burn-in entry: " << firstValid
//               << "  (" << nValid << " valid entries)" << std::endl;

//     // -----------------------------------------------------------------------
//     // Discover how many param_N branches actually exist
//     // -----------------------------------------------------------------------
//     mcmc->SetBranchStatus("*", 1);   // re-enable everything temporarily

//     int nBranchesFound = 0;
//     for (int i = 0; i < 2000; ++i) {
//         if (mcmc->GetBranch(Form("param_%d", i))) ++nBranchesFound;
//         else break;
//     }
//     std::cout << "[INFO] Found " << nBranchesFound << " param branches in MCMC." << std::endl;

//     const int nParsToRead = std::min(nBranchesFound, nXsecPars);

//     // -----------------------------------------------------------------------
//     // Enable only the branches we actually need, then bind addresses
//     // -----------------------------------------------------------------------
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

//     std::sort(chosenEntries.begin(), chosenEntries.end());  // sequential I/O

//     // Read param values for each chosen entry into memory
//     std::vector<std::vector<Double_t>> sampledParams(no_times_sampling_posterior,
//                                                      std::vector<Double_t>(nXsecPars, 0.0));
//     for (int p = 0; p < no_times_sampling_posterior; ++p) {
//         mcmc->GetEntry(chosenEntries[p]);
//         sampledParams[p].assign(xsec_tmp.begin(), xsec_tmp.end());
//     }

//     std::cout << "[INFO] Finished reading posterior samples. "
//               << "Starting reweighting loop..." << std::endl;

//     // -----------------------------------------------------------------------
//     // Posterior predictive loop – pure in-memory reweighting, no more I/O
//     // -----------------------------------------------------------------------
//     for (int p = 0; p < no_times_sampling_posterior; ++p) {
//         if (p % 50 == 0)
//             std::cout << "[INFO] Toy " << p << " / "
//                       << no_times_sampling_posterior << std::endl;

//         xsec->SetGroupOnlyParameters("Xsec", ExtractGroup(sampledParams[p], "Xsec"));
//         xsec->SetGroupOnlyParameters("Flux", ExtractGroup(sampledParams[p], "Flux"));
//         xsec->SetGroupOnlyParameters("Osc",  OscPars);

//         for (auto& pdf : DUNEPdfs) {
//             pdf->Reweight();

//             std::string pdfTitle = pdf->GetSampleTitle(0);
//             +        std::string ND_or_FD;
// +        if (pdfTitle.find("OffAxis") != std::string::npos) {
// +            ND_or_FD = "ND";
// +        } else if (pdfTitle.find("FD") != std::string::npos) {
// +            ND_or_FD = "FD";
// +        } else {
// +            std::cerr << "[WARN] Could not determine ND/FD for sample title: \""
// +                       << pdfTitle << "\" — defaulting to FD" << std::endl;
// +            ND_or_FD = "FD";
// +        }

//             MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
//             MakeSpectaVariations(pdf, "RecoNeutrinoEnergy",  fOut, ND_or_FD, pdfTitle, p);
//             MakeSpectaVariations(pdf, "Enubias",             fOut, ND_or_FD, pdfTitle, p);
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

// ---------------------------------------------------------------------------
// Determine whether a sample title corresponds to ND or FD.
// Matches on your current naming conventions:
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
// p < 0  => Asimov,  p >= 0 => posterior toy index p
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
    } else if (p < 0){
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
    const int nXsecPars = xsec->GetNumParFromGroup("Xsec")
                    + xsec->GetNumParFromGroup("Flux");

    auto prefitValues = xsec->GetPreFitValues();  // still needed later for FakeData extraction

    // Reset Xsec and Flux to their prior values directly — no vector needed
    xsec->SetGroupOnlyParameters("Xsec");
    xsec->SetGroupOnlyParameters("Flux");
    xsec->SetGroupOnlyParameters("Osc",  OscPars);

    int nFakeDataPars = xsec->GetNumParFromGroup("FakeData");
    std::vector<double> fakeDataZero(nFakeDataPars, 0.0);
    std::vector<double> fakeDataPars(nFakeDataPars);  
    xsec->SetGroupOnlyParameters("FakeData", fakeDataZero);
    // -----------------------------------------------------------------------
    // Asimov spectra (nominal parameters)
    
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


    // --- FakeData Asimov (new) ---
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

    const int nParsToRead = std::min(nBranchesFound, nXsecPars);

    // -----------------------------------------------------------------------
    // Enable only the branches we actually need, then bind addresses
    // -----------------------------------------------------------------------
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

    std::sort(chosenEntries.begin(), chosenEntries.end());  // sequential I/O

    // Read param values for each chosen entry into memory
    std::vector<std::vector<Double_t>> sampledParams(no_times_sampling_posterior,
                                                     std::vector<Double_t>(nXsecPars, 0.0));
    for (int p = 0; p < no_times_sampling_posterior; ++p) {
        mcmc->GetEntry(chosenEntries[p]);
        sampledParams[p].assign(xsec_tmp.begin(), xsec_tmp.end());
    }

    std::cout << "[INFO] Finished reading posterior samples. "
              << "Starting reweighting loop..." << std::endl;

    // -----------------------------------------------------------------------
    // Posterior predictive loop – pure in-memory reweighting, no more I/O
    // -----------------------------------------------------------------------
    for (int p = 0; p < no_times_sampling_posterior; ++p) {
        if (p % 50 == 0)
            std::cout << "[INFO] Toy " << p << " / "
                      << no_times_sampling_posterior << std::endl;

        xsec->SetGroupOnlyParameters("Xsec", ExtractGroup(sampledParams[p], "Xsec"));
        xsec->SetGroupOnlyParameters("Flux", ExtractGroup(sampledParams[p], "Flux"));
        xsec->SetGroupOnlyParameters("Osc",  OscPars);

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

