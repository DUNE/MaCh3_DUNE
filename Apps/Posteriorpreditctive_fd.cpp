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
#include <vector>
#include <map>

#include "samplePDF/GenericBinningTools.h"

#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "mcmc/mcmc.h"







// void MakeSpectaVariations(samplePDFFDBase* pdf, const std::string& var,
//                           TFile* fout, const std::string& ND_or_FD,
//                           const std::string& pdfTitle, int p) {

//     TH1D* h = dynamic_cast<TH1D*>(pdf->get1DVarHist(var.c_str()));
//     if (!h) {
//         std::cerr << "TH1D for variable '" << var
//                   << "' not found or not TH1D in PDF: " << pdfTitle << std::endl;
//         return;
//     }

    
//     static const std::map<std::string, std::vector<double>> binEdges = {
//         {"TrueNeutrinoEnergy",      {0.0, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0,
//                       3.25, 3.5, 3.75, 4.0, 5.0, 6.0, 10.0}},
//         {"Enubias",  {-2.0, -0.6, -0.581, -0.5595, -0.538, -0.5165, -0.495,
//                       -0.4735, -0.452, -0.4305, -0.409, -0.3875, -0.366,
//                       -0.3445, -0.323, -0.3015, -0.28, -0.2585, -0.237,
//                       -0.2155, -0.194, -0.1725, -0.151, -0.1295, -0.108,
//                       -0.0865, -0.065, -0.0435, -0.022, 0.0, 0.1}},
//         {"RecoNeutrinoEnergy",   {0.0, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0,
//                       3.25, 3.5, 3.75, 4.0, 5.0, 6.0, 10.0}}
//     };

//     auto it = binEdges.find(var);
//     if (it == binEdges.end()) {
//         std::cerr << "You havent added binning for this variable'" << var
//                   << "'Skip it and add later..." << std::endl;
//         return;
//     }

//     // Create custom-binned histogram
//     TH1D* h_clone = ChooseVariableBinning(
//         h,
//         Form("%s_%s_toy_%03d", pdfTitle.c_str(), var.c_str(), p),
//         it->second
//     );

//     // Create output directories
//     TDirectory* varDir = fout->GetDirectory(var.c_str());
//     if (!varDir) varDir = fout->mkdir(var.c_str());
//     TDirectory* detDir = varDir->GetDirectory(ND_or_FD.c_str());
//     if (!detDir) detDir = varDir->mkdir(ND_or_FD.c_str());
//     detDir->cd();

//     h_clone->SetDirectory(detDir);
//     h_clone->Write();
//     delete h_clone;

//     fout->cd();
// }



int main(int argc, char* argv[]) {
    if (argc == 1) {
        std::cout << "Usage: bin/ config.cfg" << std::endl;
        return 1;
    }

    auto fitMan = std::make_unique<manager>(argv[1]);
    

    auto PosteriorFile = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFiles"], __FILE__, __LINE__);
    //int Ntoys = std::stoi(Get<std::string>(fitMan->raw()["Predictive"]["Ntoys"], __FILE__, __LINE__));    
    auto burn_in = Get<unsigned int>(fitMan->raw()["General"]["MCMC"]["BurnInSteps"], __FILE__, __LINE__);
    int no_times_sampling_posterior =Get<int>(fitMan->raw()["Predictive"]["SamplePosterior"], __FILE__, __LINE__);
    std::vector<std::string> xsecCovMatrixFile =GetFromManager<std::vector<std::string>>(fitMan->raw()["General"]["Systematics"]["XsecCovFile"], {});
    auto OscCovFile =GetFromManager<std::vector<std::string>>(fitMan->raw()["General"]["Systematics"]["OscCovFile"], {});
    auto OscCovName =GetFromManager<std::string>(fitMan->raw()["General"]["Systematics"]["OscCovName"], "osc_cov");
    auto OscPars = GetFromManager<std::vector<double>>(fitMan->raw()["General"]["OscillationParameters"], {});
   
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


    for(int p=0 ; p < no_times_sampling_posterior; ++p){

        int entry = -1;
        do {
            entry = rnd->Integer(nEntries);
            mcmc->GetEntry(entry);
        } while ((unsigned int)mcmc_step < burn_in);

        mcmc->GetEntry(entry);

    //for(int p=0 ; p < no_times_sampling_posterior; ++p){

       // mcmc->GetEntry(mcmc_draw_entries[p]); // get the values of the parameters in the chain step p
        xsec->setParameters(xsec_tmp);
        for (size_t i = 0; i < std::min<size_t>(40, xsec_tmp.size()); ++i)
            //std::cout << "xsec_" << i << " = " << xsec_tmp[i] << std::endl;

        for (auto &pdf : DUNEPdfs) pdf->reweight();

        for (size_t j = 0; j < DUNEPdfs.size(); ++j) {
            DUNEPdfs[j]->reweight();
            TH1D* h_after = DUNEPdfs[j]->get1DHist();
            if (!h_after) continue;

            std::string pdfTitle = DUNEPdfs[j]->GetTitle();
            if (pdfTitle.empty()) pdfTitle = Form("FHC_numu_%02zu", j);

            
            std::string ND_or_FD = [](const std::string& title) {
                if (title.find("ND") != std::string::npos) return std::string("ND");
                if (title.find("FD") != std::string::npos) return std::string("FD");
                return std::string("Other");
            }(pdfTitle);
            
            
            TDirectory* detDir = fOut->GetDirectory(ND_or_FD.c_str());
            if (!detDir) detDir = fOut->mkdir(ND_or_FD.c_str());
            detDir->cd();


            TH1D* cloneH = static_cast<TH1D*>(h_after->Clone(
                Form("%s_posterior_%04d_toy_%03d", pdfTitle.c_str(), j, p)));
            cloneH->SetDirectory(detDir);
            cloneH->Write();
            delete cloneH;

            //MakeSpectaVariations(DUNEPdfs[j], "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
            //MakeSpectaVariations(DUNEPdfs[j], "Enubias", fOut, ND_or_FD, pdfTitle, p);
            //MakeSpectaVariations(DUNEPdfs[j], "ELep", fOut, ND_or_FD, pdfTitle, p);
            //MakeSpectaVariations(DUNEPdfs[j], "RecoNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
            //MakeSpectaVariations(DUNEPdfs[j], "q0", fOut, ND_or_FD, pdfTitle, p);
            //MakeSpectaVariations(DUNEPdfs[j], "q3", fOut, ND_or_FD, pdfTitle, p);

            fOut->cd();
        }
    }

    std::cout << "[INFO] Writing ROOT file: " << OutFileName << std::endl;
    fOut->Write();
    fOut->Close();

    delete xsec;
    delete osc;
    for (auto sample : DUNEPdfs) delete sample;

    return 0;
}
