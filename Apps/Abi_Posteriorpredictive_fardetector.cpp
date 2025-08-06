//Version for applying oscillations....
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
  //auto var_names = Get<std::string>(fitMan->raw()["Predictive"]["var_names"], __FILE__, __LINE__);
  bool plotVarspectra = true;
  std::vector<std::string> var_names = {"RecoNeutrinoEnergy"};

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

  // Do we want to fix any of the oscillation parameters? Find out from config file!
  auto oscfix = GetFromManager<std::vector<int>>(fitMan->raw()["General"]["Systematics"]["OscFix"], {});
  fOut->cd();
  /////////////////// Variations - which projections to plot

  TH1D* Hist1D = nullptr;
  for (size_t iPDF=0;iPDF<DUNEPdfs.size();iPDF++)
  {
    // Get nominal spectra and event rates
    Hist1D = (TH1D*)DUNEPdfs[iPDF]->get1DHist();

    std::vector<TH1D*> Hist1D_Vars;
    if(plotVarspectra){
      for(size_t var_i = 0 ; var_i < var_names.size() ; var_i++)
      {
        Hist1D_Vars.push_back(static_cast<TH1D*>(DUNEPdfs[iPDF]->get1DVarHist(var_names[var_i])->Clone((DUNEPdfs[iPDF]->GetTitle() + var_names[var_i]).c_str())));
      }
    }

    MACH3LOG_INFO("{:>35} | {}", DUNEPdfs[iPDF]->GetTitle(), Hist1D->Integral());
    std::cout << "[DEBUG] Writing histogram " << Hist1D << " with integral " << Hist1D->Integral() << std::endl;

    Hist1D->Write(DUNEPdfs[iPDF]->GetTitle().c_str());
    if(plotVarspectra)
    {
      for(size_t var_i = 0 ; var_i < var_names.size() ; var_i++)
      {
        Hist1D_Vars[var_i]->Write((DUNEPdfs[iPDF]->GetTitle()+var_names[var_i]).c_str());
      }
    }
  }

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
    //xsec->throwParameters();
    std::cout << "Sample " << i << " from entry " << entry << " | Step: " << mcmc_step << " | xsec[0]: " << xsec_draw[0] << std::endl;

    for (size_t j = 0; j < DUNEPdfs.size(); ++j) {
      DUNEPdfs[j]->reweight();
      TH1D* h = static_cast<TH1D*>(DUNEPdfs[j]->get1DHist());
      if (!h) continue;
      
      TString histName = Form("posterior_prediction_sample%d_pdf%d", i, (int)j);
      h->SetName(histName);
      h->SetTitle(histName);
      h->Write();
    }
    const std::vector<KinematicCut> emptySelectionVec;
    if (plotVarspectra)
    {
        for (size_t iPDF = 0; iPDF < DUNEPdfs.size(); iPDF++){
            //std::cout << "Before reweight: PDF " << iPDF << " integral: " << DUNEPdfs[iPDF]->get1DHist()->Integral() << std::endl;
            DUNEPdfs[iPDF]->reweight();
            //std::cout << "After reweight: PDF " << iPDF << " integral: " << DUNEPdfs[iPDF]->get1DHist()->Integral() << std::endl;
            for (size_t var_i = 0; var_i < var_names.size(); var_i++)
            {

                std::string var = var_names[var_i];

                TAxis tempAxis;
                tempAxis.Set(60, 0.0, 6.0); // Example: 50 bins between 0 and 5
                auto* h = DUNEPdfs[iPDF]->get1DVarHist(var, emptySelectionVec, 0, &tempAxis);

                //std::cout << "Calling get1DVarHist with var: " << var << std::endl;

                if (var.empty()) {
                    std::cerr << "ERROR: var_names[" << var_i << "] is empty!" << std::endl;
                    continue;
                }

                TString Varwritename = DUNEPdfs[iPDF]->GetTitle() + var + var_i;
                Varwritename += var_i;

                //std::cout << "Calling get1DVarHist with var: " << var << std::endl;

                DUNEPdfs[iPDF]->get1DVarHist(var)->Write(Varwritename);
                //std::cout << "Written = " << Varwritename << std::endl;
            }
        }
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
