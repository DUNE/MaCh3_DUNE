#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"
#include "mcmc/mcmc.h"
#include "manager/manager.h"
#include "covariance/covarianceBase.h"



///////////My Attempt at writing my own Posterior predictive function
int main(int argc, char* argv[]){

////////1 DUNE Setup stuff 
if (argc == 1) {
    std::cout << "Usage: bin/ config.cfg" << std::endl;
    return 1;
  }
  auto fitMan = std::unique_ptr<manager>(new manager(argv[1]));
  std::string OutFileName = GetFromManager<std::string>(fitMan->raw()["General"]["OutputFile"], "EventRates.root");

  //// Get the MCMC chain from the previous fit
  auto PosteriorFile = Get<std::string>(fitMan->raw()["Predictive"]["PosteriorFiles"], __FILE__, __LINE__); //This is the ROOT file generated afer the MCMC chain
  MCMCProcessor Processor(PosteriorFile);
  Processor.Initialise();

  //Get the burn in length from the yaml file
  auto burn_in = Get<unsigned int>(fitMan->raw()["General"]["MCMC"]["BurnInSteps"], __FILE__, __LINE__);
  int no_times_sampling_posterior = Get<int>(fitMan->raw()["Predictive"]["SamplePosterior"], __FILE__, __LINE__);
  
  TFile* fOut = new TFile(OutFileName.c_str(), "RECREATE");

  ///Set up the covariance matrices 
  std::vector<std::string> xsecCovMatrixFile = Processor.GetXSecCov();
  covarianceXsec* xsec = new covarianceXsec(xsecCovMatrixFile, "xsec_cov"); 
  covarianceOsc *osc = nullptr;

  std::vector<samplePDFFDBase *> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

  /////////set some other constants
  bool save_reweight = true;
  int mcmc_step = -1;
    
  
  
  //Now setup getting the posteriors from the previous fit
  TChain* mcmc = new TChain("posteriors"); //name of branch containing the posterior dist. of each xsec parameter
  mcmc->Add(PosteriorFile.c_str());
  
  TObjArray* branches_inmcmc = static_cast<TObjArray*>(mcmc->GetListOfBranches());
  int no_of_branches = branches_inmcmc->GetEntries();
  std::vector<TString> branch_names(no_of_branches);
  
  

 /////////////List of branches -useful for debugging
  for (int i = 0; i < no_of_branches; ++i) {
      //TBranch* branch = static_cast<TBranch*>(branches_inmcmc->At(i));
      //branch_names[i] = branches_inmcmc->GetName();
      TBranch* branch = static_cast<TBranch*>(branches_inmcmc->At(i));
      branch_names[i] = branch->GetName();  // CORRECT

      std::cout << "There are  " << no_of_branches << "in the posteriors TTree  " << std::endl;
      std::cout << "Branch " << i << "in the posteriors TTree is called : " << branch_names[i] << std::endl;
  }


  ////////Get the xsec parameters out of the posteriors TTree
std::vector<double> xsec_parameters = xsec->getNominalArray();
int no_of_xsec_branches = 0;
int no_of_steps_in_mcmc_chain = 0;

for (int i = 0; i < no_of_branches; ++i) {
  TString bname = branch_names[i];

  // Skip non-physics branches
  if (bname.Contains("LogL") || bname.Contains("_PCA"))
    continue;

  MACH3LOG_INFO("{} / {}", no_of_xsec_branches, no_of_branches);
  MACH3LOG_INFO("{}", bname);

  if (bname == "step") {
    MACH3LOG_INFO("Found Step branch");
    mcmc->SetBranchAddress("step", &mcmc_step);
    continue;
}
if (bname.Contains("xsec_")) {
    TString index_str = bname;
    index_str.ReplaceAll("xsec_", "");
    int index = index_str.Atoi();
    if (index >= 0 && index < static_cast<int>(xsec_parameters.size())) {
      mcmc->SetBranchAddress(bname, &xsec_parameters[index]);
      ++no_of_xsec_branches;
    }
  }



}
  ///////3 Generate list of random numbers for the steps in the chain

  std::vector<int> mcmc_draw_entries(no_times_sampling_posterior);
  auto rnd = std::make_unique<TRandom3>(0);  

  // Get total number of entries
  int nEntries = mcmc->GetEntries();

// Generate Ntoys random entries with burn-in check
  /*
  for (int i = 0; i < no_times_sampling_posterior; ++i) {
    int step = -1;
    int entry = -1;
    // Sample until Step >= burn_in
    while (static_cast<unsigned int>(step) < burn_in) {
      entry = rnd->Integer(nEntries);
      mcmc->GetEntry(entry);  // Loads Step value
      step = no_of_xsec_branches;
    }
    mcmc_draw_entries[i] = entry;
  }*/
  // Generate Ntoys random entries with burn-in check
  for (int i = 0; i < no_times_sampling_posterior; ++i) {
    int entry = -1;
    mcmc_step = -1;
    /*
    while (static_cast<unsigned int>(mcmc_step) < burn_in) {
      entry = rnd->Integer(nEntries);
      mcmc->GetEntry(entry);  // This fills mcmc_step!
    }*/
   while (static_cast<unsigned int>(mcmc_step) < burn_in) {
      entry = rnd->Integer(nEntries);
      mcmc->GetEntry(entry);
    // step now updated!
    }

    mcmc_draw_entries[i] = entry;
    std::cout << "[DEBUG] Entry: " << entry << " -> step: " << mcmc_step << std::endl;
  }


  //////////4 Sample those steps in the chain
  fOut->cd();
  for (int i = 0; i < no_times_sampling_posterior; ++i) {
  
    if (i % (no_times_sampling_posterior / 10) == 0) {
      MaCh3Utils::PrintProgressBar(i, no_times_sampling_posterior);
    }
    int entry = mcmc_draw_entries[i];
    mcmc->GetEntry(entry);  // Load parameters for this toy
    xsec->setParameters(xsec_parameters);
    std::cout << "Sample " << i << " from entry " << entry << " | Step: " << mcmc_step << " | xsec[0]: " << xsec_parameters[0] << std::endl;
    
    for (size_t j = 0; j < DUNEPdfs.size(); ++j) {
      DUNEPdfs[j]->reweight();
      TH1D* h = static_cast<TH1D*>(DUNEPdfs[j]->get1DHist());
      if (!h) continue;
      TString histName = Form("posterior_prediction_sample%d_pdf%d", i, (int)j);
      h->SetName(histName);
      h->SetTitle(histName);
      h->Write();
    }
  }

  std::cout << "[INFO] Writing ROOT file: " << OutFileName  << "\n";
  std::cout << "[INFO] Number of histograms in memory: " << gDirectory->GetList()->GetSize() << "\n";

  TIter next(gDirectory->GetList());
  TObject* obj;
  while ((obj = next())) {
      std::cout << "[INFO] Writing object: " << obj->GetName() << "\n";
  }
  

  delete xsec;
  delete osc;
 
  for(auto sample: DUNEPdfs) delete sample;
  fOut->Write();

  fOut->Close();
  return 0;


}