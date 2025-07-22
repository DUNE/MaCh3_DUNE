#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"
#include "mcmc/mcmc.h"
#include "manager/manager.h"
#include "covariance/covarianceBase.h"

int main(int argc, char* argv[])
{
  // Initialise manger responsible for config handling
  manager* FitManager = new manager(argv[1]);

  bool do_by_mode = GetFromManager<bool>(FitManager->raw()["NDOptions"]["PlotByMode"], false);
  bool save_reweight = true;

  // Declare outside loop/scope
  std::string outputFile;
  TFile* sknomfile = nullptr;
  TH1D* Hist1D = nullptr;
    
  // Assign once
  outputFile = Get<std::string>(FitManager->raw()["General"]["OutputFile"], __FILE__, __LINE__);
  sknomfile = new TFile(outputFile.c_str(), "RECREATE");
  if (!sknomfile || sknomfile->IsZombie()) {
    std::cerr << "[ERROR] Failed to open ROOT file for writing: " << outputFile << "\n";
    return 1;
  }

    //KS: We use MCMCProcessor to get names of covariances that were actually used to produce given chain
  auto PosteriorFile = Get<std::string>(FitManager->raw()["Predictive"]["PosteriorFiles"], __FILE__, __LINE__);
  MCMCProcessor Processor(PosteriorFile);
  Processor.Initialise();

  ///Let's ask the manager what are the file with covariance matrix
  std::vector<std::string> xsecCovMatrixFile = Processor.GetXSecCov();

  //Do you want to throw from the SK Detector covariance?
  //bool throwSKDetCov = GetFromManager<bool>(FitManager->raw()["Predictive"]["ThrowSKDetCov"], false);

  /// TODO consider passing names from config
 // std::vector<std::string> var_names;
  std::vector<std::string> var_names = GetFromManager<std::vector<std::string>>(FitManager->raw()["Predictive"]["VarNames"], {});

  MACH3LOG_INFO("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
  MACH3LOG_INFO("SIZE OF HIST VARS IS {}", var_names.size());

  bool plotVarspectra = false;
  if(var_names.size() > 0){plotVarspectra = true;}

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

  //Load DUNE samples, and also the covariance objects
  std::vector<samplePDFFDBase*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec, osc);
  //ApplyFDTuning(FitManager->get(), xsec);

  //YSP: Will need this anyway I guess. 
  //auto oscpars = Get<std::vector<double>>(FitManager->raw()["General"]["OscillationParameters"], __FILE__, __LINE__);

  // Do we want to fix any of the oscillation parameters? Find out from config file!
  auto oscfix = GetFromManager<std::vector<int>>(FitManager->raw()["General"]["Systematics"]["OscFix"], {});
  //std::vector<bool> fix_oscpar(oscpars.size(), false);
 /*
  for (size_t oscfix_i = 0; oscfix_i < oscfix.size(); oscfix_i++) {
    MACH3LOG_INFO("Fixing oscillation parameter {}", oscfix_i);
    fix_oscpar[oscfix_i] = true;
  }*/

  //osc->acceptStep();
  //osc->printPars();

  for (size_t iPDF = 0; iPDF < DUNEPdfs.size(); iPDF++)
{
  //TH1D* Hist1D = (TH1D*)DUNEPdfs[iPDF]->get1DHist();
  TH1D* tempHist = (TH1D*)DUNEPdfs[iPDF]->get1DHist();
  Hist1D = (TH1D*)tempHist->Clone(("hist_" + std::to_string(iPDF)).c_str());
  Hist1D->SetDirectory(sknomfile);  // ensure it's owned by the output file

  if (!Hist1D) {
    MACH3LOG_WARN("Null histogram returned from get1DHist() for index {}", iPDF);
    continue;
  }

  std::string histTitle = DUNEPdfs[iPDF]->GetTitle();
  if (histTitle.empty()) {
    histTitle = "PDF_" + std::to_string(iPDF);
    MACH3LOG_WARN("Empty histogram title at index {} — assigning fallback: {}", iPDF, histTitle);
  }
  sknomfile->cd();
  // Create a unique name, e.g.:
  std::string hist_name = DUNEPdfs[iPDF]->GetTitle();  // base name
  hist_name += "_step" + std::to_string(iPDF);        // make it unique

                                   // write with that name
  Hist1D->SetName(hist_name.c_str());
  Hist1D->Write();

  if (plotVarspectra) {
    for (size_t var_i = 0; var_i < var_names.size(); var_i++) {
      TH1D* varHist = (TH1D*)DUNEPdfs[iPDF]->get1DVarHist(var_names[var_i]);
      if (!varHist) {
        MACH3LOG_WARN("Null var histogram for {} / {}", histTitle, var_names[var_i]);
        continue;
      }
      std::string varName = histTitle + "_" + var_names[var_i];
      TH1D* clonedVar = (TH1D*)varHist->Clone(varName.c_str());
      clonedVar->Write();
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
  TString bname = br->GetName();  // changed from GetTitle() to GetName() for correct branch name

  if(!bname.Contains("LogL") && !bname.Contains("_PCA"))
  {
    MACH3LOG_INFO("{} / {}", ndraw, nbr);
    bnames[ndraw] = bname;
    MACH3LOG_INFO("{}", bnames[ndraw]);

    if(bnames[ndraw].BeginsWith("xsec_"))
    {
      TString indexStr = bnames[ndraw];
      indexStr.ReplaceAll("xsec_", "");
      int index = indexStr.Atoi();

      if (index >= 0 && index < static_cast<int>(xsecpars.size())) {
        mcmc->SetBranchAddress(bnames[ndraw], &xsecpars[index]);
      } else {
        MACH3LOG_WARN("Skipping branch '{}': index {} out of bounds for xsecpars (size {})",
                      bnames[ndraw].Data(), index, xsecpars.size());
      }
    }

    /*
    for (int iOsc = 0; iOsc < osc->GetNumParams(); ++iOsc) {
      if (bnames[ndraw] == osc->GetParName(iOsc) && !fix_oscpar[iOsc]) {
        mcmc->SetBranchAddress(bnames[ndraw], &(oscpars[iOsc]));
      }
    }*/

    if(bnames[ndraw] == "step"){
      MACH3LOG_INFO("Found Step branch");
      mcmc->SetBranchAddress(bnames[ndraw], &Step);
    }

    if(bnames[ndraw] == "RCreweight"){
      MACH3LOG_INFO("Found reweight branch! You've given me a chain with weights to reweight a prior!");
      mcmc->SetBranchAddress(bnames[ndraw], &reweight);
      save_reweight = true;
    }

    ndraw++;
  }
}

  TTree* reweight_tree = nullptr;
  if(save_reweight){
    reweight_tree = new TTree("reweight", "reweight");
    reweight_tree->Branch("weight", &reweight);
    //reweight_tree->Branch("sin2th_13", &(oscpars[2]));
  }
  
  auto rnd = std::make_unique<TRandom3>(0);
  int Ntoys = Get<int>(FitManager->raw()["Predictive"]["Ntoy"], __FILE__, __LINE__);

  //Get the burn-in from the config
  auto burn_in = Get<unsigned int>(FitManager->raw()["General"]["MCMC"]["BurnInSteps"], __FILE__, __LINE__);

  //////////////////////test
  std::vector<int> validEntries;
  for (int e=0; e < mcmc->GetEntries(); e++) {
      mcmc->GetEntry(e);
      if (Step >= burn_in) validEntries.push_back(e);
  }

  std::cout << "about to do toys loop ..." << std::endl;
  for(int i = 0; i < Ntoys; i++)
  {
    std::cout << "in toys loop ..." << std::endl;

      if( i % (Ntoys/10) == 0){MaCh3Utils::PrintProgressBar(i, Ntoys);}
      int entry = 0;
      Step = -999; //YSP: Ensures you get an entry from the mcmc even when burn_in is set to zero (Although not advised :p ).  
      while(static_cast<unsigned int>(Step) < burn_in){
         std::cout << "in while ..." << std::endl;//Take 200k burn in steps, WP: Eb C in 1st peaky
        entry = rnd->Integer(mcmc->GetEntries());

        if (Step < static_cast<int>(burn_in)) {
       std::cerr << "[DEBUG] Step = " << Step << " < burn-in = " << burn_in << "\n";
        }
      int randomIndex = rnd->Integer(validEntries.size());
      mcmc->GetEntry(validEntries[randomIndex]);

        //mcmc->GetEntry(entry);
        std::cout << "Toy " << i << ": Step=" << Step << ", xsecpars[0]=" << xsecpars[0] << std::endl;

      }
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

      for (size_t iPDF=0; iPDF < DUNEPdfs.size(); iPDF++)
      {
        DUNEPdfs[iPDF]->reweight();
      }
      
      if(do_by_mode)
      {
        for (size_t iPDF = 0; iPDF < DUNEPdfs.size(); iPDF++) {
          for(int mode_i= 0; mode_i<=DUNEPdfs[iPDF]->GetMaCh3Modes()->GetNModes(); mode_i++) {
            TString mode = DUNEPdfs[iPDF]->GetMaCh3Modes()->GetMaCh3ModeName(mode_i);
            TH1D* ModeHist = static_cast<TH1D*>(DUNEPdfs[iPDF]->getModeHist1D(1, 1, 2));
            //Reset every time
            ModeHist->Reset();

            for(int samples_i = 0 ; samples_i < DUNEPdfs[iPDF]->getNMCSamples() ; samples_i++) {
              ModeHist->Add(DUNEPdfs[iPDF]->getModeHist1D(samples_i,mode_i, 2));
            }
            TString writename = DUNEPdfs[iPDF]->GetTitle()+"_"+mode;
            writename+=i;
            ModeHist->Write(writename);
          }
        } 
      }
      else{
        for (size_t iPDF = 0;iPDF<DUNEPdfs.size(); iPDF++)
        {
          Hist1D = static_cast<TH1D*>(DUNEPdfs[iPDF]->get1DHist());
          TString writename = DUNEPdfs[iPDF]->GetTitle()+"_";
          writename+=i;
          Hist1D->Write(writename);
        }
      }
      if(plotVarspectra)
      {
        for (size_t iPDF = 0; iPDF < DUNEPdfs.size(); iPDF++)
        {
          for(size_t var_i =  0 ; var_i < var_names.size() ; var_i++)
          {
            TString Varwritename = DUNEPdfs[iPDF]->GetTitle()+var_names[var_i];
            Varwritename+=i;
            //DUNEPdfs[iPDF]->get1DVarHist(var_names[i])->Write(Varwritename);
            if (var_i < var_names.size()) {
            auto* varHist = DUNEPdfs[iPDF]->get1DVarHist(var_names[var_i]);
            if (varHist) {
              varHist->Write(Varwritename);
            } else {
              std::cerr << "[WARNING] get1DVarHist returned nullptr for var '" 
                        << var_names[var_i] << "' in PDF: " << DUNEPdfs[iPDF]->GetTitle() << "\n";
            }
      }

          }
        }
      }
      if(save_reweight)reweight_tree->Fill();
  }//end of toys loop
  if (!sknomfile || sknomfile->IsZombie()) {
    std::cerr << "[ERROR] Failed to open ROOT file for writing: " << outputFile << "\n";
    return 1;
}

  if(save_reweight) reweight_tree->Write();

  std::cout << "[INFO] Writing ROOT file: " << outputFile << "\n";
std::cout << "[INFO] Number of histograms in memory: " << gDirectory->GetList()->GetSize() << "\n";

TIter next(gDirectory->GetList());
TObject* obj;
while ((obj = next())) {
    std::cout << "[INFO] Writing object: " << obj->GetName() << "\n";
}
  sknomfile->Write(); 
  sknomfile->Close();
  delete sknomfile;

  delete xsec;
  delete osc;
 
  for(auto sample: DUNEPdfs) delete sample;
  return 0;
}