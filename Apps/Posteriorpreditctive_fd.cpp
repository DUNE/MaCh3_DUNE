#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"
#include "mcmc/mcmc.h"
#include "manager/manager.h"
#include "covariance/covarianceBase.h"

#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <string>
#include <filesystem> 

int main(int argc, char* argv[]) {
    if (argc == 1) {
        std::cout << "Usage: bin/ config.cfg" << std::endl;
        return 1;
    }

    auto fitMan = std::make_unique<manager>(argv[1]);
    

    auto PosteriorFile = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFiles"], __FILE__, __LINE__);
    int Ntoys = std::stoi(Get<std::string>(fitMan->raw()["Predictive"]["Ntoys"], __FILE__, __LINE__));    
    auto burn_in = Get<unsigned int>(fitMan->raw()["General"]["MCMC"]["BurnInSteps"], __FILE__, __LINE__);
    int no_times_sampling_posterior =Get<int>(fitMan->raw()["Predictive"]["SamplePosterior"], __FILE__, __LINE__);
    std::vector<std::string> xsecCovMatrixFile =GetFromManager<std::vector<std::string>>(fitMan->raw()["General"]["Systematics"]["XsecCovFile"], {});
    auto OscCovFile =GetFromManager<std::vector<std::string>>(fitMan->raw()["General"]["Systematics"]["OscCovFile"], {});
    auto OscCovName =GetFromManager<std::string>(fitMan->raw()["General"]["Systematics"]["OscCovName"], "osc_cov");
    auto OscPars = GetFromManager<std::vector<double>>(fitMan->raw()["General"]["OscillationParameters"], {});



    if(Ntoys ==1){
        std::cout << "NTOYS ==1 which will give useless error bars" << std::endl;
    }
    if(no_times_sampling_posterior < Ntoys ){
        std::cout << "You are sampling the posterior dist fewer times than the number of toys... not sensible" << std::endl;
    }

    // Extract directory path
    std::filesystem::path dir = std::filesystem::path(PosteriorFile).parent_path();
    std::string OutFileName = (dir / "Posteriorpredictive_out.root").string();
    std::cout << "Output file will be: " << OutFileName << std::endl;


    TFile* fOut = new TFile(OutFileName.c_str(), "RECREATE");
    MCMCProcessor Processor(PosteriorFile);
    Processor.Initialise();

    covarianceXsec* xsec = new covarianceXsec(xsecCovMatrixFile, "xsec_cov");
    covarianceOsc* osc = new covarianceOsc(OscCovFile, OscCovName);
    
    osc->setParameters(OscPars);


    //Add extra bit to output file so I cant overwrite my posterior TTree again...
    std::string suffix = "_posteriorpredictive_withtoys";
    size_t dotPos = OutFileName.find_last_of('.');
    if (dotPos != std::string::npos) {
        OutFileName.insert(dotPos, suffix);
    } else {
        OutFileName += suffix;
    }


    std::vector<samplePDFFDBase*> DUNEPdfs;
    MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

    // --- TTree with posteriors in ---
    TChain* mcmc = new TChain("posteriors");
    mcmc->Add(PosteriorFile.c_str());
    Long64_t nEntries = mcmc->GetEntries();
    if (nEntries == 0) {
        std::cerr << "[ERROR] MCMC tree has no entries!" << std::endl;
        return 1;
    }

    // --- Bind xsec branches names to match name to branch correctly
    std::vector<double> xsec_nominal = xsec->getNominalArray();
    std::vector<double> xsec_draw(xsec_nominal.size(), 0.0);
    std::vector<Double_t> xsec_tmp(xsec_nominal.size(), 0.0);

    for (size_t i = 0; i < xsec_nominal.size(); ++i) {
        TString bname = Form("xsec_%zu", i);
        if (mcmc->GetBranch(bname)) {
            mcmc->SetBranchAddress(bname, &xsec_tmp[i]);
        }
    }

    int mcmc_step = -1;
    mcmc->SetBranchAddress("step", &mcmc_step);

    // Drawing from posterior
    auto rnd = std::make_unique<TRandom3>(0);
    std::vector<int> mcmc_draw_entries(no_times_sampling_posterior, -1);

    for (int i = 0; i < no_times_sampling_posterior; ++i) {
        int entry = -1;
        do {
            entry = rnd->Integer(nEntries);
            mcmc->GetEntry(entry);
        } while ((unsigned int)mcmc_step < burn_in);
        mcmc_draw_entries[i] = entry;
    }

    // --- Loop over posterior draws ---
    for (int iSample = 0; iSample < no_times_sampling_posterior; ++iSample) {
        int entry = mcmc_draw_entries[iSample];
        mcmc->GetEntry(entry);

        for (size_t j = 0; j < xsec_draw.size(); ++j)
            xsec_draw[j] = xsec_nominal[j] + (xsec_tmp[j] - xsec_nominal[j]);
        xsec->setParameters(xsec_draw);

        //Loop through each toy to get variations
        for (int iToy = 0; iToy < Ntoys; ++iToy) {
            xsec->throwParameters(); //throw xsec params from TTree
            //osc->throwParameters();  
            
            for (auto pdf : DUNEPdfs) pdf->reweight();

            // Write histograms per toy
            for (size_t j = 0; j < DUNEPdfs.size(); ++j) {
                TH1D* h_after = DUNEPdfs[j]->get1DHist();
                if (!h_after) continue;

                std::string pdfTitle = DUNEPdfs[j]->GetTitle();
                if (pdfTitle.empty()) pdfTitle = Form("FHC_numu_%02zu", j);

                
                auto GetNDorFD = [](const std::string& title) {
                    if (title.find("ND") != std::string::npos) return std::string("ND");
                    if (title.find("FD") != std::string::npos) return std::string("FD");
                    return std::string("Other");
                };

                
                std::string ND_or_FD = GetNDorFD(pdfTitle);
                std::replace(ND_or_FD.begin(), ND_or_FD.end(), ' ', '_');

                TDirectory* detDir = fOut->GetDirectory(ND_or_FD.c_str());
                if (!detDir) detDir = fOut->mkdir(ND_or_FD.c_str());
                detDir->cd();


                TH1D* cloneH = static_cast<TH1D*>(h_after->Clone(
                    Form("%s_posterior_%04d_toy_%03d", pdfTitle.c_str(), iSample, iToy)));
                cloneH->SetDirectory(detDir);
                cloneH->Write();
                delete cloneH;

                fOut->cd();
            }
        }
    }

    // --- Wrap up ---
    std::cout << "[INFO] Writing ROOT file: " << OutFileName << std::endl;
    fOut->Write();
    fOut->Close();

    delete xsec;
    delete osc;
    for (auto sample : DUNEPdfs) delete sample;

    return 0;
}
