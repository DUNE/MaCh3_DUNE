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
// Helper: get or create a single directory level under parent
// ---------------------------------------------------------------------------
static TDirectory* GetOrMkdir(TDirectory* parent, const std::string& name) {
    TDirectory* d = parent->GetDirectory(name.c_str());
    if (!d) d = parent->mkdir(name.c_str());
    return d;
}
void MakeSpectaVariations(SampleHandlerFD* pdf,
                          const std::string& var,
                          TFile* fout,
                          const std::string& ND_or_FD,
                          const std::string& pdfTitle,
                          int p,
                          int sample_idx = 0)  {
    
    std::vector<KinematicCut> emptyEventCuts;
    std::vector<KinematicCut> emptySubEventCuts;
    
    TH1* h = pdf->Get1DVarHist(
        sample_idx,           // Sample index
        var,                  // Variable name
        emptyEventCuts,       // Event selection cuts
        0,                    // WeightStyle (0 = use weights)
        nullptr,              // Use default binning
        emptySubEventCuts     // Sub-event selection cuts
    );
    
    TH1D* h1d = dynamic_cast<TH1D*>(h);
    if (!h1d) {
        std::cerr << "ERROR: Could not cast to TH1D" << std::endl;
        return;
    }
    
    //TH1D* h = dynamic_cast<TH1D*>(pdf->Get1DVarHist(var.c_str()));
    if (!h) {
        std::cerr << "[WARN] Could not get TH1D for " << var << " from " << pdfTitle << std::endl;
        return;
    }

    std::string cleanTitle = pdfTitle;
    std::replace(cleanTitle.begin(), cleanTitle.end(), ' ', '_');
    std::replace(cleanTitle.begin(), cleanTitle.end(), '/', '_');

    std::string baseDir, histName;
    if (p < 0) {
        baseDir = "Asimov/" + ND_or_FD + "/" + var;
        histName = Form("%s_%s_%s_Asimov",
                        ND_or_FD.c_str(), cleanTitle.c_str(), var.c_str());
    } else {
        baseDir = ND_or_FD + "/" + var;
        histName = Form("%s_%s_%s_posterior_toy_%03d",
                        ND_or_FD.c_str(), cleanTitle.c_str(), var.c_str(), p);
    }

    TDirectory* dir = fout->GetDirectory(baseDir.c_str());
    if (!dir) dir = fout->mkdir(baseDir.c_str());
    dir->cd();

    TH1D* cloneH = static_cast<TH1D*>(h->Clone(histName.c_str()));
    cloneH->SetDirectory(dir);
    cloneH->Write();

    delete cloneH;
    delete h;  // cleanup to avoid memory leak
    fout->cd();
}


void MakeSpectaVariations2D(SampleHandlerFD* pdf,
                           const std::string& varX,
                           const std::string& varY,
                           TFile* fout,
                           const std::string& ND_or_FD,
                           const std::string& pdfTitle,
                           int p,
                           int sample_idx = 0,
                           std::map<std::string,
                           std::vector<std::vector<std::vector<double>>>>* posteriorStore = nullptr)
{
    std::vector<KinematicCut> emptyEventCuts;
    std::vector<KinematicCut> emptySubEventCuts;
    //pdf->Reweight();

    // Get 2D histogram (same philosophy as 1D)
    TH2* h = pdf->Get2DVarHist(
        sample_idx,
        varX,
        varY,
        emptyEventCuts,
        0,              // weighted
        nullptr,        // default binning
        nullptr,
        emptySubEventCuts
    );

    if (!h) {
    std::cerr << "[WARN] null histogram\n";
    return;
    }

    TH2D* h2d = dynamic_cast<TH2D*>(h);
    if (!h2d) {
        std::cerr << "[ERROR] not TH2D\n";
        delete h;
        return;
    }

    // Clean name
    std::string cleanTitle = pdfTitle;
    std::replace(cleanTitle.begin(), cleanTitle.end(), ' ', '_');
    std::replace(cleanTitle.begin(), cleanTitle.end(), '/', '_');

    std::string baseDir, histName;
    if (p < 0) {
        baseDir = "Asimov/" + ND_or_FD + "/2D";
        histName = Form("%s_%s_%s_vs_%s_Asimov",
                        ND_or_FD.c_str(), cleanTitle.c_str(),
                        varX.c_str(), varY.c_str());
    } else {
        baseDir = ND_or_FD + "/2D";
        histName = Form("%s_%s_%s_vs_%s_posterior_toy_%03d",
                        ND_or_FD.c_str(), cleanTitle.c_str(),
                        varX.c_str(), varY.c_str(), p);
    }

    TDirectory* dir = fout->GetDirectory(baseDir.c_str());
    if (!dir) dir = fout->mkdir(baseDir.c_str());
    dir->cd();

    TH2D* cloneH = static_cast<TH2D*>(h2d->Clone(histName.c_str()));

    // Store posterior values (same as before)
    if (p >= 0 && posteriorStore) {
        std::string key = ND_or_FD + "_" + cleanTitle + "_" + varX + "_vs_" + varY;

        int nX = cloneH->GetNbinsX();
        int nY = cloneH->GetNbinsY();

        if ((*posteriorStore)[key].empty()) {
            (*posteriorStore)[key].resize(nX);
            for (int ix = 0; ix < nX; ++ix)
                (*posteriorStore)[key][ix].resize(nY);
        }

        for (int ix = 1; ix <= nX; ++ix)
            for (int iy = 1; iy <= nY; ++iy)
                (*posteriorStore)[key][ix-1][iy-1]
                    .push_back(cloneH->GetBinContent(ix, iy));
    }

    cloneH->SetDirectory(dir);
    cloneH->Write();

    delete cloneH;
    delete h;   // SAFE here because Get2DVarHist allocates new TH2D

    fout->cd();
}

int main(int argc, char* argv[]) {
    if (argc == 1) {
        std::cout << "Usage: bin/ config.cfg" << std::endl;
        return 1;
    }

    // --- Manager setup
    auto fitMan = std::make_unique<Manager>(argv[1]);
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
        TString bname = Form("param_%zu", i);
        if (mcmc->GetBranch(bname))
            mcmc->SetBranchAddress(bname, &xsec_tmp[i]);
    }

    UInt_t mcmc_step = -1;
    mcmc->SetBranchAddress("step", &mcmc_step);

    // --- Generate Asimov spectra
    xsec->SetGroupOnlyParameters("Xsec", xsec_nominal);
    xsec->SetGroupOnlyParameters("Osc", OscPars);
    
    for (auto& pdf : DUNEPdfs) {
        std::string pdfTitle = pdf->GetName();
        if (pdfTitle.empty()) pdfTitle = "FHC_numu_asimov";

        std::string ND_or_FD =
            (pdfTitle.find("ND") != std::string::npos) ? "ND" :
            (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";

        // All directory creation is handled inside MakeSpectra* via GetOrMkdir
        pdf->Reweight();
        MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, -1);

       
        MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, -1);

       
        MakeSpectaVariations(pdf, "Enubias", fOut, ND_or_FD, pdfTitle, -1);

       MakeSpectaVariations2D(pdf,
                       "TrueNeutrinoEnergy",
                       "Enubias",
                       fOut,
                       ND_or_FD,
                       pdfTitle,
                       -1);
        MakeSpectaVariations2D(pdf,
                       "RecoNeutrinoEnergy",
                       "ELepRec",
                       fOut,
                       ND_or_FD,
                       pdfTitle,
                       -1);
        

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
        std::string pdfTitle = pdf->GetName();
        std::string ND_or_FD =
            (pdfTitle.find("ND") != std::string::npos) ? "ND" :
            (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";

        // // +shift
        // xsec->SetParameters(xsec_plus);
        // pdf->Reweight();
        // TH1* h_plus = pdf->GetMCHist(1);
        // if (!h_plus) continue;

        // TDirectory* shiftDir = GetOrMkdir(fOut, "shift_parameters");
        // TDirectory* detDir   = GetOrMkdir(shiftDir, ND_or_FD);
        // detDir->cd();

        // TH1D* cloneHplus = static_cast<TH1D*>(h_plus->Clone(
        //     Form("%s_%s_plus", ND_or_FD.c_str(), pdfTitle.c_str())));
        // cloneHplus->SetDirectory(detDir);
        // cloneHplus->Write();
        // delete cloneHplus;

        // // -shift
        // xsec->SetParameters(xsec_minus);
        // pdf->Reweight();
        // TH1* h_minus = pdf->GetMCHist(1);
        // if (!h_minus) continue;

        // TH1D* cloneHminus = static_cast<TH1D*>(h_minus->Clone(
        //     Form("%s_%s_minus", ND_or_FD.c_str(), pdfTitle.c_str())));
        // cloneHminus->SetDirectory(detDir);
        // cloneHminus->Write();
        // delete cloneHminus;

        fOut->cd();
    }

   // --- Posterior predictive draws
auto rnd = std::make_unique<TRandom3>(0);
const Long64_t maxSampleSteps = nEntries;

// Build valid entry list once, before the loop
std::vector<Long64_t> validEntries;
Long64_t limit = std::min<Long64_t>(maxSampleSteps, nEntries);

std::cout << "[INFO] Scanning chain for post-burn-in entries...\n";
for (Long64_t i = 0; i < limit; ++i) {
    mcmc->GetEntry(i);
    if (mcmc_step >= burn_in)
        validEntries.push_back(i);
}
std::cout << "[INFO] Found " << validEntries.size()
          << " valid entries after burn-in of " << burn_in << "\n";

if (validEntries.empty()) {
    std::cerr << "[ERROR] No entries pass burn-in cut!\n";
    return 1;
}

for (int p = 0; p < no_times_sampling_posterior; ++p) {
    Long64_t entry = validEntries[rnd->Integer(validEntries.size())];
    mcmc->GetEntry(entry);
    xsec->SetParameters(std::vector<double>(xsec_tmp.begin(), xsec_tmp.end()));

    for (auto& pdf : DUNEPdfs) {
        std::string pdfTitle = pdf->GetName();
        std::string ND_or_FD =
            (pdfTitle.find("ND") != std::string::npos) ? "ND" :
            (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";

        MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
        MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
        MakeSpectaVariations(pdf, "Enubias",            fOut, ND_or_FD, pdfTitle, p);
        MakeSpectaVariations2D(pdf,
                       "TrueNeutrinoEnergy",
                       "Enubias",
                       fOut,
                       ND_or_FD,
                       pdfTitle,
                       p,
                       0,
                       &posteriorStore);
        MakeSpectaVariations2D(pdf,
                       "RecoNeutrinoEnergy",
                       "ELepRec",
                       fOut,
                       ND_or_FD,
                       pdfTitle,
                       p,
                       0,
                       &posteriorStore);
        
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
                std::string hName = key + "_Asimov";
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
                         hRef->GetXaxis()->GetName() + ";" +
                         hRef->GetYaxis()->GetName() + ";" +
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