// #include "samplePDFDUNE/MaCh3DUNEFactory.h"
// #include "samplePDFDUNE/StructsDUNE.h"
// #include "mcmc/mcmc.h"
// #include "manager/manager.h"
// #include "covariance/covarianceBase.h"

// #include <iomanip>
// #include <iostream>
// #include <memory>
// #include <vector>
// #include <string>
// #include <cmath>
// #include <algorithm>
// #include <string>
// #include <filesystem> 
// #include <vector>
// #include <map>

// #include "samplePDF/GenericBinningTools.h"

// #include "samplePDFDUNE/MaCh3DUNEFactory.h"
// #include "mcmc/mcmc.h"


// void MakeSpectaVariations(samplePDFFDBase* pdf, const std::string& var,
//                           TFile* fout, const std::string& ND_or_FD,
//                           const std::string& pdfTitle, int p) {

//     TH1D* h = dynamic_cast<TH1D*>(pdf->get1DVarHist(var.c_str()));
//     if (!h) {
//         std::cerr << "TH1D for variable '" << var
//                   << "' not found or not TH1D in PDF: " << pdfTitle << std::endl;
//         return;
//     }

//     // std::string baseDir = "Asimov/" + ND_or_FD;
//     // TDirectory* varDir = fout->GetDirectory(var.c_str());
//     std::string baseDir = "Asimov/" + ND_or_FD + "/" + var;
//     TDirectory* detDir = fout->GetDirectory(baseDir.c_str());
//     if (!detDir) detDir = fout->mkdir(baseDir.c_str());
//     detDir->cd();

//     if (!varDir) varDir = fout->mkdir(var.c_str());
//     TDirectory* detDir = varDir->GetDirectory(baseDir.c_str());
//     if (!detDir) detDir = varDir->mkdir(baseDir.c_str());
//     detDir->cd();
//     h->SetName(Form("%s_%s_%s_%s", ND_or_FD.c_str(), pdfTitle.c_str(), var.c_str(), "Asimov"));


//     h->SetDirectory(detDir);
//     h->Write();
//     delete h;

//     fout->cd();
// }



// int main(int argc, char* argv[]) {
//     if (argc == 1) {
//         std::cout << "Usage: bin/ config.cfg" << std::endl;
//         return 1;
//     }

//     auto fitMan = std::make_unique<manager>(argv[1]);
//     auto PosteriorFile = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFiles"], __FILE__, __LINE__);
//     //int Ntoys = std::stoi(Get<std::string>(fitMan->raw()["Predictive"]["Ntoys"], __FILE__, __LINE__));    
//     auto burn_in = Get<unsigned int>(fitMan->raw()["General"]["MCMC"]["BurnInSteps"], __FILE__, __LINE__);
//     int no_times_sampling_posterior =Get<int>(fitMan->raw()["Predictive"]["SamplePosterior"], __FILE__, __LINE__);
//     std::vector<std::string> xsecCovMatrixFile =GetFromManager<std::vector<std::string>>(fitMan->raw()["General"]["Systematics"]["XsecCovFile"], {});
//     auto OscCovFile =GetFromManager<std::vector<std::string>>(fitMan->raw()["General"]["Systematics"]["OscCovFile"], {});
//     auto OscCovName =GetFromManager<std::string>(fitMan->raw()["General"]["Systematics"]["OscCovName"], "osc_cov");
//     auto OscPars = GetFromManager<std::vector<double>>(fitMan->raw()["General"]["OscillationParameters"], {});

//     double shiftAmount = GetFromManager<double>(fitMan->raw()["Predictive"]["ShiftAmount"], 0.2);
//     std::cout << "[INFO] Using shift amount ±" << shiftAmount << std::endl;

   
//     // Extract directory path
//     std::filesystem::path dir = std::filesystem::path(PosteriorFile).parent_path();
//     std::string OutFileName = (dir / "Posteriorpredictive_out.root").string();
//     std::cout << "Output file will be: " << OutFileName << std::endl;


//     TFile* fOut = new TFile(OutFileName.c_str(), "RECREATE");
//     MCMCProcessor Processor(PosteriorFile);
//     Processor.Initialise();

//     covarianceXsec* xsec = new covarianceXsec(xsecCovMatrixFile, "xsec_cov");
//     covarianceOsc* osc = new covarianceOsc(OscCovFile, OscCovName);
    
//     osc->setParameters(OscPars);


//     //Add extra bit to output file so I cant overwrite my posterior TTree again...
//     std::string suffix = "_posteriorpredictive_withtoys";
//     size_t dotPos = OutFileName.find_last_of('.');
//     if (dotPos != std::string::npos) {
//         OutFileName.insert(dotPos, suffix);
//     } else {
//         OutFileName += suffix;
//     }


//     std::vector<samplePDFFDBase*> DUNEPdfs;
//     MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

//      // --- TTree with posteriors in ---
//     TChain* mcmc = new TChain("posteriors");
//     mcmc->Add(PosteriorFile.c_str());
//     Long64_t nEntries = mcmc->GetEntries();
//     if (nEntries == 0) {
//         std::cerr << "[ERROR] MCMC tree has no entries!" << std::endl;
//         return 1;
//     }

//     // --- Bind xsec branches names to match name to branch correctly
//     std::vector<double> xsec_nominal = xsec->getNominalArray();
//     std::vector<double> xsec_draw(xsec_nominal.size(), 0.0);
//     std::vector<Double_t> xsec_tmp(xsec_nominal.size(), 0.0);

//     for (size_t i = 0; i < xsec_nominal.size(); ++i) {
//         TString bname = Form("xsec_%zu", i);
//         if (mcmc->GetBranch(bname)) {
//             mcmc->SetBranchAddress(bname, &xsec_tmp[i]);
//         }
//     }

//     int mcmc_step = -1;
//     mcmc->SetBranchAddress("step", &mcmc_step);


//     // Reset parameters to nominal values
//     xsec->setParameters(xsec_nominal);
//     osc->setParameters(OscPars);

//     // Reweight all DUNE PDFs for the Asimov prediction
//     for (auto &pdf : DUNEPdfs) {
//         pdf->reweight();

//         TH1D* h_asimov = pdf->get1DHist();
//         if (!h_asimov) continue;

//         std::string pdfTitle = pdf->GetTitle();
//         if (pdfTitle.empty()) pdfTitle = "FHC_numu_asimov";

//         // Identify detector type
//         std::string ND_or_FD = [](const std::string& title) {
//             if (title.find("ND") != std::string::npos) return std::string("ND");
//             if (title.find("FD") != std::string::npos) return std::string("FD");
//             return std::string("Other");
//         }(pdfTitle);

//         // Create Asimov subdirectory (Asimov/ND, Asimov/FD, etc.)
//         TDirectory* asimovDir = fOut->GetDirectory("Asimov");
//         if (!asimovDir) asimovDir = fOut->mkdir("Asimov");
//         asimovDir->cd();

//         TDirectory* detDir = asimovDir->GetDirectory(ND_or_FD.c_str());
//         if (!detDir) detDir = asimovDir->mkdir(ND_or_FD.c_str());
//         detDir->cd();

//         // Clone and write the histogram
//         //TH1D* cloneH = static_cast<TH1D*>(h_asimov->Clone(Form("%s_Asimov", pdfTitle.c_str())));
//         std::string cleanTitle = pdfTitle;
//         std::replace(cleanTitle.begin(), cleanTitle.end(), ' ', '_'); // avoid spaces

//         TH1D* cloneH = static_cast<TH1D*>(h_asimov->Clone(Form("%s_%s_Asimov", ND_or_FD.c_str(), cleanTitle.c_str())));
//         cloneH->SetDirectory(detDir);
//         cloneH->Write();
//         delete cloneH;

//         // Optionally write Asimov spectral variations
//         MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, "Asimov/" + ND_or_FD, pdfTitle, -1);
//         MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, "Asimov/" + ND_or_FD, pdfTitle, -1);

//         fOut->cd();
//     }

//     std::vector<double> error(xsec_nominal.size(), shiftAmount);
//     std::vector<double> xsec_plus(xsec_nominal.size());
//     std::vector<double> xsec_minus(xsec_nominal.size());

//     for (size_t i = 0; i < xsec_nominal.size(); ++i) {
//         xsec_plus[i] = xsec_nominal[i] + error[i];
//         xsec_minus[i] = xsec_nominal[i] - error[i];
//     }
//         // Set oscillation parameters once
//     osc->setParameters(OscPars);

//     for (auto &pdf : DUNEPdfs) {
//         xsec->setParameters(xsec_plus);
//         pdf->reweight();

//         TH1D* h_plus = pdf->get1DHist();
//         if (!h_plus) continue;

//         std::string pdfTitle = pdf->GetTitle();
//         if (pdfTitle.empty()) pdfTitle = "FHC_numu";

//         std::string ND_or_FD = [](const std::string& title) {
//             if (title.find("ND") != std::string::npos) return std::string("ND");
//             if (title.find("FD") != std::string::npos) return std::string("FD");
//             return std::string("Other");
//         }(pdfTitle);

       
//         TDirectory* shiftDir = fOut->GetDirectory("shift_parameters");
//         if (!shiftDir) shiftDir = fOut->mkdir("shift_parameters");
//         shiftDir->cd();

//         TDirectory* detDir = shiftDir->GetDirectory(ND_or_FD.c_str());
//         if (!detDir) detDir = shiftDir->mkdir(ND_or_FD.c_str());
//         detDir->cd();

        
//         TH1D* cloneHplus = static_cast<TH1D*>(h_plus->Clone(Form("%s_%s_plus", ND_or_FD.c_str(), pdfTitle.c_str())));

//         cloneHplus->SetDirectory(detDir);
//         cloneHplus->Write();
//         delete cloneHplus;

//         //MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, "shift_parameters/" + ND_or_FD, pdfTitle + "_plus", -1);
//         //MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, "shift_parameters/" + ND_or_FD, pdfTitle + "_plus", -1);

    
//         xsec->setParameters(xsec_minus);
//         pdf->reweight();

//         TH1D* h_minus = pdf->get1DHist();
//         if (!h_minus) continue;

//         TH1D* cloneHminus = static_cast<TH1D*>(h_minus->Clone(Form("%s_%s_plus", ND_or_FD.c_str(), pdfTitle.c_str())));
//         cloneHminus->SetDirectory(detDir);
//         cloneHminus->Write();
//         delete cloneHminus;

//         //MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, "shift_parameters/" + ND_or_FD, pdfTitle + "_minus", -1);
//         //MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, "shift_parameters/" + ND_or_FD, pdfTitle + "_minus", -1);

//         fOut->cd();
//     }



//       // Drawing from posterior
//     auto rnd = std::make_unique<TRandom3>(0);
//     std::vector<int> mcmc_draw_entries(no_times_sampling_posterior, -1);
//     // Limit sampling to first 2 million steps
//     const Long64_t maxSampleSteps = 2000000;

//     for (int p = 0; p < no_times_sampling_posterior; ++p) {

//         int entry = -1;
//         do {
//             // Draw only from [burn_in, maxSampleSteps)
//             entry = rnd->Integer(std::min<Long64_t>(maxSampleSteps, nEntries));
//             mcmc->GetEntry(entry);
//         } while ((unsigned int)mcmc_step < burn_in);

//         mcmc->GetEntry(entry);


//     // for(int p=0 ; p < no_times_sampling_posterior; ++p){

//     //     int entry = -1;
//     //     do {
//     //         entry = rnd->Integer(nEntries);
//     //         mcmc->GetEntry(entry);
//     //     } while ((unsigned int)mcmc_step < burn_in);

//     //     mcmc->GetEntry(entry);

//     //for(int p=0 ; p < no_times_sampling_posterior; ++p){

//        // mcmc->GetEntry(mcmc_draw_entries[p]); // get the values of the parameters in the chain step p
//         xsec->setParameters(xsec_tmp);
//         for (size_t i = 0; i < std::min<size_t>(40, xsec_tmp.size()); ++i)
//             //std::cout << "xsec_" << i << " = " << xsec_tmp[i] << std::endl;

//         for (auto &pdf : DUNEPdfs) pdf->reweight();

//         for (size_t j = 0; j < DUNEPdfs.size(); ++j) {
//             DUNEPdfs[j]->reweight();
//             TH1D* h_after = DUNEPdfs[j]->get1DHist();
//             if (!h_after) continue;

//             std::string pdfTitle = DUNEPdfs[j]->GetTitle();
//             if (pdfTitle.empty()) pdfTitle = Form("FHC_numu_%02zu", j);

            
//             std::string ND_or_FD = [](const std::string& title) {
//                 if (title.find("ND") != std::string::npos) return std::string("ND");
//                 if (title.find("FD") != std::string::npos) return std::string("FD");
//                 return std::string("Other");
//             }(pdfTitle);
            
            
//             TDirectory* detDir = fOut->GetDirectory(ND_or_FD.c_str());
//             if (!detDir) detDir = fOut->mkdir(ND_or_FD.c_str());
//             detDir->cd();


//             TH1D* cloneH = static_cast<TH1D*>(h_after->Clone(Form("%s_%s_posterior_toy_%03d_hist_%02d", ND_or_FD.c_str(), pdfTitle.c_str(), p, j)));

//             cloneH->SetDirectory(detDir);
//             cloneH->Write();
//             delete cloneH;

//             MakeSpectaVariations(DUNEPdfs[j], "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
//             MakeSpectaVariations(DUNEPdfs[j], "RecoNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
            
//             fOut->cd();
//         }
//     }

//     std::cout << "[INFO] Writing ROOT file: " << OutFileName << std::endl;
//     fOut->Write();
//     fOut->Close();

//     delete xsec;
//     delete osc;
//     for (auto sample : DUNEPdfs) delete sample;

//     return 0;
// }


// #include "samplePDFDUNE/MaCh3DUNEFactory.h"
// #include "samplePDFDUNE/StructsDUNE.h"
// #include "mcmc/mcmc.h"
// #include "manager/manager.h"
// #include "covariance/covarianceBase.h"

// #include <iomanip>
// #include <iostream>
// #include <memory>
// #include <vector>
// #include <string>
// #include <cmath>
// #include <algorithm>
// #include <string>
// #include <filesystem> 
// #include <map>

// #include "samplePDF/GenericBinningTools.h"


// void MakeSpectaVariations(samplePDFFDBase* pdf, const std::string& var,
//                           TFile* fout, const std::string& ND_or_FD,
//                           const std::string& pdfTitle, int p) {

//     TH1D* h = dynamic_cast<TH1D*>(pdf->get1DVarHist(var.c_str()));
//     if (!h) {
//         std::cerr << "TH1D for variable '" << var
//                   << "' not found or not TH1D in PDF: " << pdfTitle << std::endl;
//         return;
//     }

//     std::string baseDir = ND_or_FD;
//     TDirectory* varDir = fout->GetDirectory(var.c_str());
//     if (!varDir) varDir = fout->mkdir(var.c_str());
//     TDirectory* detDir = varDir->GetDirectory(baseDir.c_str());
//     if (!detDir) detDir = varDir->mkdir(baseDir.c_str());
//     detDir->cd();

//     h->SetDirectory(detDir);
//     h->Write();
//     delete h;

//     fout->cd();
// }


// int main(int argc, char* argv[]) {
//     if (argc == 1) {
//         std::cout << "Usage: bin/ config.cfg" << std::endl;
//         return 1;
//     }

//     auto fitMan = std::make_unique<manager>(argv[1]);
    
//     auto PosteriorFile = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFiles"], __FILE__, __LINE__);
//     auto burn_in = Get<unsigned int>(fitMan->raw()["General"]["MCMC"]["BurnInSteps"], __FILE__, __LINE__);
//     int no_times_sampling_posterior = Get<int>(fitMan->raw()["Predictive"]["SamplePosterior"], __FILE__, __LINE__);
//     std::vector<std::string> xsecCovMatrixFile = GetFromManager<std::vector<std::string>>(fitMan->raw()["General"]["Systematics"]["XsecCovFile"], {});
//     auto OscCovFile = GetFromManager<std::vector<std::string>>(fitMan->raw()["General"]["Systematics"]["OscCovFile"], {});
//     auto OscCovName = GetFromManager<std::string>(fitMan->raw()["General"]["Systematics"]["OscCovName"], "osc_cov");
//     auto OscPars = GetFromManager<std::vector<double>>(fitMan->raw()["General"]["OscillationParameters"], {});

//     std::filesystem::path dir = std::filesystem::path(PosteriorFile).parent_path();
//     std::string OutFileName = (dir / "Posteriorpredictive_out.root").string();
//     std::cout << "Output file will be: " << OutFileName << std::endl;

//     TFile* fOut = new TFile(OutFileName.c_str(), "RECREATE");
//     MCMCProcessor Processor(PosteriorFile);
//     Processor.Initialise();

//     covarianceXsec* xsec = new covarianceXsec(xsecCovMatrixFile, "xsec_cov");
//     covarianceOsc* osc = new covarianceOsc(OscCovFile, OscCovName);
//     osc->setParameters(OscPars);

//     std::vector<samplePDFFDBase*> DUNEPdfs;
//     MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

//     std::cout<<"Line 400" <<std::endl;
//     // --- TTree with posteriors in ---
//     TChain* mcmc = new TChain("posteriors");
//     mcmc->Add(PosteriorFile.c_str());
//     Long64_t nEntries = mcmc->GetEntries();
//     if (nEntries == 0) {
//         std::cerr << "[ERROR] MCMC tree has no entries!" << std::endl;
//         return 1;
//     }

//     // --- Bind xsec branches ---
//     std::vector<double> xsec_nominal = xsec->getNominalArray();
//     std::vector<Double_t> xsec_tmp(xsec_nominal.size(), 0.0);

//     for (size_t i = 0; i < xsec_nominal.size(); ++i) {
//         TString bname = Form("xsec_%zu", i);
//         if (mcmc->GetBranch(bname)) {
//             mcmc->SetBranchAddress(bname, &xsec_tmp[i]);
//         }
//     }

//     int mcmc_step = -1;
//     mcmc->SetBranchAddress("step", &mcmc_step);

//     // Reset to nominal
//     xsec->setParameters(xsec_nominal);
//     osc->setParameters(OscPars);

//     // --- Asimov predictions ---
//     std::cout<<"About to do Asimiv" <<std::endl;
//     for (auto &pdf : DUNEPdfs) {
//         pdf->reweight();

//         TH1D* h_asimov = pdf->get1DHist();
//         if (!h_asimov) continue;

//         std::string pdfTitle = pdf->GetTitle();
//         if (pdfTitle.empty()) pdfTitle = "FHC_numu_asimov";

//         std::string ND_or_FD = (pdfTitle.find("ND") != std::string::npos)
//             ? "ND" : (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";

//         TDirectory* asimovDir = fOut->GetDirectory("Asimov");
//         if (!asimovDir) asimovDir = fOut->mkdir("Asimov");
//         TDirectory* detDir = asimovDir->GetDirectory(ND_or_FD.c_str());
//         if (!detDir) detDir = asimovDir->mkdir(ND_or_FD.c_str());
//         detDir->cd();

//         TH1D* cloneH = static_cast<TH1D*>(h_asimov->Clone(Form("%s_Asimov", pdfTitle.c_str())));
//         cloneH->SetDirectory(detDir);
//         cloneH->Write();
//         delete cloneH;

//         MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, "Asimov/" + ND_or_FD, pdfTitle, -1);
//         MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, "Asimov/" + ND_or_FD, pdfTitle, -1);
//     }

//     std::cout<<"About to do toys" <<std::endl;

//     // --- Posterior predictive toys ---
//     auto rnd = std::make_unique<TRandom3>(0);
//     const Long64_t maxSampleSteps = 2000000;

//     for (int p = 0; p < no_times_sampling_posterior; ++p) {
//         int entry = -1;
//         int maxRetries = 1000000;
//         int retries = 0;
//         do {
//             entry = rnd->Integer(std::min<Long64_t>(maxSampleSteps, nEntries));
//             mcmc->GetEntry(entry);
//             retries++;
//             if (retries > maxRetries) {
//                 std::cerr << "[ERROR] Could not find valid posterior entry after burn-in!\n";
//                 return 1;
//             }
//         } while ((unsigned int)mcmc_step < burn_in);

//         xsec->setParameters(xsec_tmp);

//         for (auto &pdf : DUNEPdfs) pdf->reweight();

//         for (size_t j = 0; j < DUNEPdfs.size(); ++j) {
//             TH1D* h_after = DUNEPdfs[j]->get1DHist();
//             if (!h_after) continue;

//             std::string pdfTitle = DUNEPdfs[j]->GetTitle();
//             if (pdfTitle.empty()) pdfTitle = Form("FHC_numu_%02zu", j);

//             std::string ND_or_FD = (pdfTitle.find("ND") != std::string::npos)
//                 ? "ND" : (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";

//             TDirectory* detDir = fOut->GetDirectory(ND_or_FD.c_str());
//             if (!detDir) detDir = fOut->mkdir(ND_or_FD.c_str());
//             detDir->cd();

//             TH1D* cloneH = static_cast<TH1D*>(h_after->Clone(
//                 Form("%s_posterior_toy_%03d", pdfTitle.c_str(), p)));
//             cloneH->SetDirectory(detDir);
//             cloneH->Write();
//             delete cloneH;

//             MakeSpectaVariations(DUNEPdfs[j], "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
//             MakeSpectaVariations(DUNEPdfs[j], "RecoNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
//         }
//     }

//     std::cout << "[INFO] Writing ROOT file: " << OutFileName << std::endl;
//     fOut->Write();
//     fOut->Close();

//     delete xsec;
//     delete osc;
//     for (auto sample : DUNEPdfs) delete sample;

//     return 0;
// }


#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"
#include "mcmc/mcmc.h"
#include "manager/manager.h"
#include "covariance/covarianceBase.h"

#include "samplePDF/GenericBinningTools.h"

#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <map>

// void MakeSpectaVariations(samplePDFFDBase* pdf, const std::string& var,
//                           TFile* fout, const std::string& ND_or_FD,
//                           const std::string& pdfTitle, int p) {
//     TH1D* h = nullptr;

//     // Force regeneration for Reco histograms (since get1DVarHist caches)
//     if (var == "RecoNeutrinoEnergy") {
//         h = pdf->get1DHist();  // Always recomputed via fill1DHist()
//     } else {
//         h = dynamic_cast<TH1D*>(pdf->get1DVarHist(var.c_str()));
//     }

//     if (!h) {
//         std::cerr << "[WARN] TH1D for variable '" << var
//                   << "' not found or not TH1D in PDF: " << pdfTitle << std::endl;
//         return;
//     }

//     // --- Clean title ---
//     std::string cleanTitle = pdfTitle;
//     std::replace(cleanTitle.begin(), cleanTitle.end(), ' ', '_');
//     std::replace(cleanTitle.begin(), cleanTitle.end(), '/', '_');

//     std::string baseDir;
//     std::string histName;

//     if (p < 0) {
//         // Asimov mode
//         baseDir = "Asimov/" + ND_or_FD + "/" + var;
//         histName = Form("%s_%s_%s_Asimov",
//                         ND_or_FD.c_str(), cleanTitle.c_str(), var.c_str());
//     } else {
//         // Posterior toy mode
//         baseDir = ND_or_FD + "/" + var;
//         histName = Form("%s_%s_%s_posterior_toy_%03d",
//                         ND_or_FD.c_str(), cleanTitle.c_str(), var.c_str(), p);
//     }

//     // --- Save histogram ---
//     TDirectory* dir = fout->GetDirectory(baseDir.c_str());
//     if (!dir) dir = fout->mkdir(baseDir.c_str());
//     dir->cd();

//     TH1D* cloneH = static_cast<TH1D*>(h->Clone(histName.c_str()));
//     cloneH->SetDirectory(dir);
//     cloneH->Write();
//     delete cloneH;

//     fout->cd();
// }

void MakeSpectaVariations(samplePDFFDBase* pdf, const std::string& var,
                          TFile* fout, const std::string& ND_or_FD,
                          const std::string& pdfTitle, int p) {
    // Refresh internal event weights
    (void) pdf->get1DHist();  // calls fill1DHist() internally
    pdf->reweight();
    TH1D* h = dynamic_cast<TH1D*>(pdf->get1DVarHist(var.c_str()));
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
    covarianceXsec* xsec = new covarianceXsec(xsecCovMatrixFile, "xsec_cov");
    covarianceOsc* osc = new covarianceOsc(OscCovFile, OscCovName);
    osc->setParameters(OscPars);

    // --- Build PDF instances
    std::vector<samplePDFFDBase*> DUNEPdfs;
    MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

    // --- Load posteriors
    TChain* mcmc = new TChain("posteriors");
    mcmc->Add(PosteriorFile.c_str());
    Long64_t nEntries = mcmc->GetEntries();
    if (nEntries == 0) {
        std::cerr << "[ERROR] MCMC tree has no entries!" << std::endl;
        return 1;
    }

    // Bind xsec branches
    std::vector<double> xsec_nominal = xsec->getNominalArray();
    std::vector<Double_t> xsec_tmp(xsec_nominal.size(), 0.0);
    for (size_t i = 0; i < xsec_nominal.size(); ++i) {
        TString bname = Form("xsec_%zu", i);
        if (mcmc->GetBranch(bname))
            mcmc->SetBranchAddress(bname, &xsec_tmp[i]);
    }

    int mcmc_step = -1;
    mcmc->SetBranchAddress("step", &mcmc_step);

    // --- Generate Asimov spectra
    xsec->setParameters(xsec_nominal);
    osc->setParameters(OscPars);

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
        pdf->reweight();  // refresh weights for Asimov
        MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, -1);

        // --- Reco Neutrino Energy ---
        pdf->reweight();  // refresh weights for Reco
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

    osc->setParameters(OscPars);
    for (auto& pdf : DUNEPdfs) {
        // +shift
        xsec->setParameters(xsec_plus);
        pdf->reweight();
        TH1D* h_plus = pdf->get1DHist();
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
        xsec->setParameters(xsec_minus);
        pdf->reweight();
        TH1D* h_minus = pdf->get1DHist();
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
    const Long64_t maxSampleSteps = 2000000;

    for (int p = 0; p < no_times_sampling_posterior; ++p) {
        int entry;
        do {
            entry = rnd->Integer(std::min<Long64_t>(maxSampleSteps, nEntries));
            mcmc->GetEntry(entry);
        } while ((unsigned int)mcmc_step < burn_in);

        xsec->setParameters(xsec_tmp);

        for (auto& pdf : DUNEPdfs) {
            pdf->reweight();
            TH1D* h_after = pdf->get1DHist();
            if (!h_after) continue;

            std::string pdfTitle = pdf->GetTitle();
            std::string ND_or_FD =
                (pdfTitle.find("ND") != std::string::npos) ? "ND" :
                (pdfTitle.find("FD") != std::string::npos) ? "FD" : "Other";

            TDirectory* detDir = fOut->GetDirectory(ND_or_FD.c_str());
            if (!detDir) detDir = fOut->mkdir(ND_or_FD.c_str());
            detDir->cd();

            TH1D* cloneH = static_cast<TH1D*>(h_after->Clone(
                Form("%s_%s_posterior_toy_%03d", ND_or_FD.c_str(), pdfTitle.c_str(), p)));
            cloneH->SetDirectory(detDir);
            cloneH->Write();
            delete cloneH;

            MakeSpectaVariations(pdf, "TrueNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
            MakeSpectaVariations(pdf, "RecoNeutrinoEnergy", fOut, ND_or_FD, pdfTitle, p);
        }
        fOut->cd();
    }

    std::cout << "[INFO] Writing ROOT file: " << OutFileName << std::endl;
    fOut->Write();
    fOut->Close();

    delete xsec;
    delete osc;
    for (auto sample : DUNEPdfs) delete sample;

    return 0;
}
