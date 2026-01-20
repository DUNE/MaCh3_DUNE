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

    // Add suffix before opening
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

        // Create Asimov directory
        TDirectory* asimovDir = fOut->GetDirectory("Asimov");
        if (!asimovDir) asimovDir = fOut->mkdir("Asimov");
        asimovDir->cd();

        TDirectory* detDir = asimovDir->GetDirectory(ND_or_FD.c_str());
        if (!detDir) detDir = asimovDir->mkdir(ND_or_FD.c_str());
        detDir->cd();

        // --- True Neutrino Energy ---
        pdf->Reweight();  // refresh weights for Asimov
        MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, -1);

        // --- Reco Neutrino Energy ---
        pdf->Reweight();  // refresh weights for Reco
        MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, -1);

        fOut->cd();
    }


    // --- ±shift systematic variations
    std::vector<double> error(xsec_nominal.size(), shiftAmount);
    std::vector<double> xsec_plus(xsec_nominal.size()), xsec_minus(xsec_nominal.size());
    for (size_t i = 0; i < xsec_nominal.size(); ++i) {
        xsec_plus[i] = xsec_nominal[i] + error[i];
        xsec_minus[i] = xsec_nominal[i] - error[i];
    }

    xsec->SetGroupOnlyParameters("Osc", OscPars);
    for (auto& pdf : DUNEPdfs) {
        // +shift
        xsec->SetParameters(xsec_plus);
        pdf->Reweight();
        TH1* h_plus = pdf->GetMCHist(1);
        if (!h_plus) continue;

        std::string pdfTitle = pdf->GetTitle();
        std::string ND_or_FD =
            (pdfTitle.find("ND") != std::string::npos) ? "ND" :
            (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";

        TDirectory* shiftDir = fOut->GetDirectory("shift_parameters");
        if (!shiftDir) shiftDir = fOut->mkdir("shift_parameters");
        shiftDir->cd();

        TDirectory* detDir = shiftDir->GetDirectory(ND_or_FD.c_str());
        if (!detDir) detDir = shiftDir->mkdir(ND_or_FD.c_str());
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
    // --- Posterior predictive draws
    auto rnd = std::make_unique<TRandom3>(0);
    const Long64_t maxSampleSteps = 2000000;

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

            // Generate histograms with correct binning for each variable
            MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
            MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
        }
        fOut->cd();
    }
    std::cout << "[INFO] Writing ROOT file: " << OutFileName << std::endl;
    fOut->Write();
    fOut->Close();

    delete xsec;
    for (auto sample : DUNEPdfs) delete sample;

    return 0;
}