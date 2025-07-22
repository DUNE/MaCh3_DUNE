//MaCh3 includes
#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"
#include "mcmc/mcmc.h"
#include "manager/manager.h"
#include "covariance/covarianceBase.h"

int main(int argc, char* argv[])
{
  if (argc < 2) {
    std::cerr << "Usage: Fit config.yaml" << std::endl;
    return 1;
  }

  // Initialise manger responsible for config handling
  //auto FitManager = MaCh3ManagerFactory(argc, argv);

  manager* FitManager = new manager(argv[1]);
  //TFile* sknomfile = new TFile(outputFile.c_str(), "RECREATE");
  //sknomfile->cd();
  

  bool do_by_mode = GetFromManager<bool>(FitManager->raw()["NDOptions"]["PlotByMode"], false);
  bool save_reweight = true;
  
  //KS: We use MCMCProcessor to get names of covariances that were actually used to produce given chain
  auto PosteriorFile = Get<std::string>(FitManager->raw()["Predictive"]["PosteriorFiles"], __FILE__, __LINE__);
  MCMCProcessor Processor(PosteriorFile);
  Processor.Initialise();

  ///Let's ask the manager what are the file with covariance matrix
  std::vector<std::string> xsecCovMatrixFile = Processor.GetXSecCov();

  //Do you want to throw from the SK Detector covariance?
  bool throwSKDetCov = GetFromManager<bool>(FitManager->raw()["Predictive"]["ThrowSKDetCov"], false);

  /// TODO consider passing names from config
  std::vector<std::string> var_names;
  MACH3LOG_INFO("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
  MACH3LOG_INFO("SIZE OF HIST VARS IS {}", var_names.size());

  bool plotVarspectra = false;
  if(var_names.size() > 0){plotVarspectra = true;}

  if(throwSKDetCov){
    MACH3LOG_INFO("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    MACH3LOG_INFO("YOU ASKED TO THROW FROM SK DET COV - SO I WILL!");
    MACH3LOG_INFO("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
  }
  else{
    MACH3LOG_INFO("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
    MACH3LOG_INFO("YOU ASKED NOT TO THROW FROM SK DET COV - SO I WON'T!");
    MACH3LOG_INFO("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
  }
  //Gets xsec covariance directly from the MCMC
  covarianceXsec* xsec = new covarianceXsec(xsecCovMatrixFile, "xsec_cov");
  covarianceOsc* osc = nullptr;
  //covarianceSkDet_joint* skdet = nullptr;

  bool is_skpriorpredictive = Get<bool>(FitManager->raw()["Predictive"]["FixXSec"], __FILE__, __LINE__);
  if(is_skpriorpredictive)
    MACH3LOG_INFO("You've chosen to run SK Prior Predictive Distribution");
  else
    MACH3LOG_INFO("You've chosen to run SK Posterior Predictive Distribution");

  //Load T2K samples, and also the covariance objects
  std::vector<samplePDFFDBase*> T2KPdfs;
  //MakeMaCh3SKInstance(FitManager.get(), T2KPdfs, xsec, skdet, osc);
  MakeMaCh3DuneInstance(FitManager, T2KPdfs, xsec, osc);
  //ApplyFDTuning(FitManager.get(), xsec);
  //ApplyFDTuning(FitManager, xsec);


  //YSP: Will need this anyway I guess. 
  auto oscpars = Get<std::vector<double>>(FitManager->raw()["General"]["OscillationParameters"], __FILE__, __LINE__);

  // Do we want to fix any of the oscillation parameters? Find out from config file!
  auto oscfix = GetFromManager<std::vector<int>>(FitManager->raw()["General"]["Systematics"]["OscFix"], {});
  std::vector<bool> fix_oscpar(oscpars.size(), false);

  for (size_t oscfix_i = 0; oscfix_i < oscfix.size(); oscfix_i++) {
    MACH3LOG_INFO("Fixing oscillation parameter {}", oscfix_i);
    fix_oscpar[oscfix_i] = true;
  }

  osc->acceptStep();
  osc->printPars();

  // One set of osc pars
  for (size_t iPDF = 0; iPDF < T2KPdfs.size(); iPDF++)
  {
    T2KPdfs[iPDF]->reweight();
  }
  xsec->printNominal();

  auto outputFile = Get<std::string>(FitManager->raw()["General"]["OutputFile"], __FILE__, __LINE__);
  TFile* sknomfile = new TFile(outputFile.c_str(), "RECREATE");
  sknomfile->cd();
  TH1D* Hist1D = nullptr;
  for (size_t iPDF=0;iPDF<T2KPdfs.size();iPDF++)
  {
    // Get nominal spectra and event rates
    Hist1D = (TH1D*)T2KPdfs[iPDF]->get1DHist();

    std::vector<TH1D*> Hist1D_Vars;
    if(plotVarspectra){
      for(size_t var_i = 0 ; var_i < var_names.size() ; var_i++)
      {
        Hist1D_Vars.push_back(static_cast<TH1D*>(T2KPdfs[iPDF]->get1DVarHist(var_names[var_i])->Clone((T2KPdfs[iPDF]->GetTitle() + var_names[var_i]).c_str())));
      }
    }

    MACH3LOG_INFO("{:>35} | {}", T2KPdfs[iPDF]->GetTitle(), Hist1D->Integral());
    Hist1D->Write(T2KPdfs[iPDF]->GetTitle().c_str());
    if(plotVarspectra)
    {
      for(size_t var_i = 0 ; var_i < var_names.size() ; var_i++)
      {
        Hist1D_Vars[var_i]->Write((T2KPdfs[iPDF]->GetTitle()+var_names[var_i]).c_str());
      }
    }
  }

  TChain* mcmc = new TChain("posteriors");
  mcmc->Add(PosteriorFile.c_str());
  
  TObjArray* brlis = static_cast<TObjArray*>(mcmc->GetListOfBranches());
  int nbr = brlis->GetEntries();
  std::vector<TString> bnames(nbr);
  int ndraw = 0;
  int Step = 0;
  double reweight = 1; //reweight to reweight prior if "reweight" branch found. Otherwise this is just 1
  std::vector<double> xsecpars = xsec->getNominalArray();
  for(int i = 0; i < nbr; i++)
  {
    TBranch* br = static_cast<TBranch*>(brlis->At(i));
    TString bname = br->GetTitle();
    if(!bname.Contains("LogL") && !bname.Contains("_PCA"))
    {
      MACH3LOG_INFO("{} / {}", ndraw, nbr);
      bnames[ndraw] = bname;
      MACH3LOG_INFO("{}", bnames[ndraw]);

      if(bnames[ndraw].Contains("xsec_"))
      {
        TString index = bnames[ndraw];
        index.ReplaceAll("xsec_","");

        mcmc->SetBranchAddress(bnames[ndraw],&(xsecpars[index.Atoi()]));
      }
      //These branches don't exist for an ND fit
      //They only exist for an SK fit
      //So for an ND chain the SK detector parameters
      //will be fixed unless throwSKDetCov is true.
      //If these branches do exist then set throwSKDetCov
      //to false
      /*if(bnames[ndraw].Contains("skd_joint_"))
      {
        TString index = bnames[ndraw];
        index.ReplaceAll("skd_joint_","");
        mcmc->SetBranchAddress(bnames[ndraw],&(skdpars[index.Atoi()]));
        throwSKDetCov = false;
      }*/
      // KS: Find osc params
      for (int iOsc = 0; iOsc < osc->GetNumParams(); ++iOsc) {
        if (bnames[ndraw] == osc->GetParName(iOsc) && !fix_oscpar[iOsc]) {
          mcmc->SetBranchAddress(bnames[ndraw], &(oscpars[iOsc]));
        }
      }

      if(bnames[ndraw] == "step"){
        MACH3LOG_INFO("Found Step branch");
        mcmc->SetBranchAddress(bnames[ndraw], &Step);
      }

      if(bnames[ndraw]=="RCreweight"){
        MACH3LOG_INFO("Found reweight branch! You've given me a chain with weights to reweight a prior!");
        mcmc->SetBranchAddress(bnames[ndraw], &reweight);
        save_reweight=true;
      }
      ndraw++;
    }
  }


  TTree* reweight_tree = nullptr;
  if(save_reweight){
    reweight_tree = new TTree("reweight", "reweight");
    reweight_tree->Branch("weight", &reweight);
    reweight_tree->Branch("sin2th_13", &(oscpars[2]));
  }
  
  auto rnd = std::make_unique<TRandom3>(0);
  int Ntoys = Get<int>(FitManager->raw()["Predictive"]["Ntoy"], __FILE__, __LINE__);

  //Get the burn-in from the config
  auto burn_in = Get<unsigned int>(FitManager->raw()["General"]["MCMC"]["BurnInSteps"], __FILE__, __LINE__);

  for(int i = 0; i < Ntoys; i++)
  {
      if( i % (Ntoys/10) == 0){MaCh3Utils::PrintProgressBar(i, Ntoys);}
      int entry = 0;
      Step = -999; //YSP: Ensures you get an entry from the mcmc even when burn_in is set to zero (Although not advised :p ).  
      while(static_cast<unsigned int>(Step) < burn_in){//Take 200k burn in steps, WP: Eb C in 1st peaky
        entry = rnd->Integer(mcmc->GetEntries());
        mcmc->GetEntry(entry);
      }
      
      xsec->setParameters(xsecpars);
      //KS: Below line can help you get prior predictive distributions which are helpful for getting pre and post ND fit spectra
      //YSP: If not set in the config, the code runs SK Posterior Predictive distributions by default. If true, then the code runs SK prior predictive.
      if(is_skpriorpredictive) xsec->throwParameters();
      
      //KS: If you want to estiamte error coming only from one systematic (for example flux only) play with below settings
      if (GetFromManager<bool>(FitManager->raw()["Predictive"]["FixXSec"], false, __FILE__, __LINE__)) {
        xsec->SetGroupOnlyParameters("Xsec");
      }
      if (GetFromManager<bool>(FitManager->raw()["Predictive"]["FixFlux"], false, __FILE__, __LINE__)) {
        xsec->SetGroupOnlyParameters("Flux");
      }
      
     
      //osc->setParameters(oscpars);

      for (size_t iPDF=0; iPDF < T2KPdfs.size(); iPDF++)
      {
        T2KPdfs[iPDF]->reweight();
      }
      
      if(do_by_mode)
      {
        for (size_t iPDF = 0; iPDF < T2KPdfs.size(); iPDF++) {
          for(int mode_i= 0; mode_i<=T2KPdfs[iPDF]->GetMaCh3Modes()->GetNModes(); mode_i++) {
            TString mode = T2KPdfs[iPDF]->GetMaCh3Modes()->GetMaCh3ModeName(mode_i);
            TH1D* ModeHist = static_cast<TH1D*>(T2KPdfs[iPDF]->getModeHist1D(1, 1, 2));
            //Reset every time
            ModeHist->Reset();

            for(int samples_i = 0 ; samples_i < T2KPdfs[iPDF]->getNMCSamples() ; samples_i++) {
              ModeHist->Add(T2KPdfs[iPDF]->getModeHist1D(samples_i,mode_i, 2));
            }
            TString writename = T2KPdfs[iPDF]->GetTitle()+"_"+mode;
            writename+=i;
            ModeHist->Write(writename);
          }
        } 
      }
      else{
        for (size_t iPDF = 0;iPDF<T2KPdfs.size(); iPDF++)
        {
          Hist1D = static_cast<TH1D*>(T2KPdfs[iPDF]->get1DHist());
          if (!Hist1D) {
            std::cerr << "[DEBUG] get1DHist returned nullptr for PDF index " << iPDF << std::endl;
            continue; // skip to next iteration
          }
          TString writename = T2KPdfs[iPDF]->GetTitle()+"_";
          writename+=i;
          Hist1D->Write(writename);
        }
      }
      if(plotVarspectra)
      {
        for (size_t iPDF = 0; iPDF < T2KPdfs.size(); iPDF++)
        {
          for(size_t var_i =  0 ; var_i < var_names.size() ; var_i++)
          {
            auto* varHist = T2KPdfs[iPDF]->get1DVarHist(var_names[var_i]);
            if (!varHist) {
    std::cerr << "[DEBUG] get1DVarHist returned nullptr for var: " << var_names[var_i] << std::endl;
    continue;
}
              TString Varwritename = T2KPdfs[iPDF]->GetTitle()+var_names[var_i];
              Varwritename+=i;
            T2KPdfs[iPDF]->get1DVarHist(var_names[i])->Write(Varwritename);
          }
        }
      }
      if(save_reweight)reweight_tree->Fill();
  }//end of toys loop
  if (save_reweight) reweight_tree->Write();
  sknomfile->Write();  // <-- Important: Write all objects
  sknomfile->Close();
  delete sknomfile;

  delete xsec;
  delete osc;
  
  for(auto sample: T2KPdfs) delete sample;
  return 0;
}



