#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"
#include "mcmc/mcmc.h"
#include "manager/manager.h"
#include "covariance/covarianceBase.h"

int main(int argc, char* argv[]) {
  /////////// 1. Setup
  if (argc == 1) {
    std::cout << "Usage: bin/ config.cfg" << std::endl;
    return 1;
  }

  auto fitMan = std::unique_ptr<manager>(new manager(argv[1]));
  std::string OutFileName = GetFromManager<std::string>(fitMan->raw()["General"]["OutputFile"], "EventRates.root");

  auto PosteriorFile = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFiles"], __FILE__, __LINE__);
  MCMCProcessor Processor(PosteriorFile);
  Processor.Initialise();

  auto burn_in = Get<unsigned int>(fitMan->raw()["General"]["MCMC"]["BurnInSteps"], __FILE__, __LINE__);
  int no_times_sampling_posterior = Get<int>(fitMan->raw()["Predictive"]["SamplePosterior"], __FILE__, __LINE__);

  TFile* fOut = new TFile(OutFileName.c_str(), "RECREATE");

  /////////// 2. Covariance and PDFs
  std::vector<std::string> xsecCovMatrixFile = Processor.GetXSecCov();
  covarianceXsec* xsec = new covarianceXsec(xsecCovMatrixFile, "xsec_cov");
  covarianceOsc* osc = nullptr;

  std::vector<samplePDFFDBase*> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

  /////////// 3. Set up MCMC tree
  TChain* mcmc = new TChain("posteriors");
  mcmc->Add(PosteriorFile.c_str());

  TObjArray* branches_inmcmc = static_cast<TObjArray*>(mcmc->GetListOfBranches());
  int no_of_branches = branches_inmcmc->GetEntries();
  std::vector<TString> branch_names(no_of_branches);

  for (int i = 0; i < no_of_branches; ++i) {
    TBranch* branch = static_cast<TBranch*>(branches_inmcmc->At(i));
    branch_names[i] = branch->GetName();
   // std::cout << "Branch " << i << ": " << branch_names[i] << std::endl;
  }

  /////////// 4. Bind branches
  std::vector<double> xsec_nominal = xsec->getNominalArray();
  std::vector<double> xsec_draw(xsec_nominal.size(), 0.0);

  int mcmc_step = -1;
  mcmc->SetBranchAddress("step", &mcmc_step);
  
  int no_of_xsec_branches = 0;

  for (const auto& bname : branch_names) {
    if (bname.Contains("LogL") || bname.Contains("_PCA")) continue;

    if (bname.BeginsWith("xsec_")) {
      TString index_str = bname;
      index_str.ReplaceAll("xsec_", "");
      int index = index_str.Atoi();
      if (index >= 0 && index < static_cast<int>(xsec_draw.size())) {
        mcmc->SetBranchAddress(bname, &xsec_draw[index]);
        ++no_of_xsec_branches;
      }
    }
  }

  /////////// 5. Sample steps with burn-in check
  int nEntries = mcmc->GetEntries();
  auto rnd = std::make_unique<TRandom3>(0);
  std::vector<int> mcmc_draw_entries(no_times_sampling_posterior);

  for (int i = 0; i < no_times_sampling_posterior; ++i) {
   int entry = -1;
    do {
      entry = rnd->Integer(nEntries);
      mcmc->GetEntry(entry);
    } while ((unsigned int)mcmc_step < burn_in);

    mcmc_draw_entries[i] = entry;
    //std::cout << "[DEBUG] Entry: " << entry << " -> step: " << mcmc_step << std::endl;
  }

  /////////// 6. Reweight and write predictions
  fOut->cd();

  for (int i = 0; i < no_times_sampling_posterior; ++i) {
    if (i % (no_times_sampling_posterior / 10) == 0) {
      MaCh3Utils::PrintProgressBar(i, no_times_sampling_posterior);
    }

    int entry = mcmc_draw_entries[i];
    mcmc->GetEntry(entry);

    xsec->setParameters(xsec_draw);

    //std::cout << "Sample " << i << " from entry " << entry << " | Step: " << mcmc_step << " | xsec[0]: " << xsec_draw[0] << std::endl;

    for (size_t j = 0; j < DUNEPdfs.size(); ++j) {
      DUNEPdfs[j]->reweight();
      TH1D* h = static_cast<TH1D*>(DUNEPdfs[j]->get1DHist());
      if (!h) continue;

      TH2D* h2 = static_cast<TH2D*>(DUNEPdfs[j]->get2DHist());
      if (!h2) continue;

      // Project onto RecoNeutrinoEnergy (X-axis)
      TH1D* hRecoNuE = h2->ProjectionX();

      // Project onto ELepRec (Y-axis)
      TH1D* hELepRec = h2->ProjectionY();

      /*std::cout << "PDF " << j << " axis: " 
          << h->GetXaxis()->GetTitle() 
          << ", bins: " << h->GetNbinsX() 
          << ", range: [" << h->GetXaxis()->GetXmin() 
          << ", " << h->GetXaxis()->GetXmax() << "]"
          << std::endl;*/

      TString histNameRecoNuE = Form("posterior_prediction_sample%d_projRecoNuE", i);
      hRecoNuE->SetName(histNameRecoNuE);
      hRecoNuE->SetTitle("RecoNeutrinoEnergy projection");
      hRecoNuE->Write();

      TString histNameELepRec = Form("posterior_prediction_sample%d_projELepRec", i);
      hELepRec->SetName(histNameELepRec);
      hELepRec->SetTitle("ELepRec projection");
      hELepRec->Write();

      TString histName = Form("posterior_prediction_sample%d_pdf%d", i, (int)j);
      h->SetName(histName);
      h->SetTitle(histName);
      h->Write();
    }
  }
  /////////// 7. Wrap up
  std::cout << "[INFO] Writing ROOT file: " << OutFileName << "\n";
  std::cout << "[INFO] Number of histograms in memory: " << gDirectory->GetList()->GetSize() << "\n";

  TIter next(gDirectory->GetList());
  TObject* obj;
  while ((obj = next())) {
    std::cout << "[INFO] Writing object: " << obj->GetName() << "\n";
  }

  delete xsec;
  delete osc;
  for (auto sample : DUNEPdfs) delete sample;

  fOut->Write();
  fOut->Close();

  return 0;
}
