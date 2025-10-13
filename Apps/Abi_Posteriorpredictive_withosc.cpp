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
  std::vector<std::string> xsecCovMatrixFile =  GetFromManager<std::vector<std::string>>(fitMan->raw()["General"]["Systematics"]["XsecCovFile"], {}); // Processor.GetXSecCov();
  //std::vector<std::string> xsecCovMatrixFile = {"/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/CovObjs/Erecyrec_xsec_closuretest.yaml"};//{"/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/configs/CovObjs/TrueNeutrinoEnergy_Enubias_relative.yaml"};//Processor.GetXSecCov();
  covarianceXsec* xsec = new covarianceXsec(xsecCovMatrixFile, "xsec_cov");
  auto OscCovFile = GetFromManager<std::vector<std::string>>(fitMan->raw()["General"]["Systematics"]["OscCovFile"], {});
  auto OscCovName = GetFromManager<std::string>(fitMan->raw()["General"]["Systematics"]["OscCovName"], "osc_cov");

  covarianceOsc* osc = new covarianceOsc(OscCovFile, OscCovName);

  auto OscPars = GetFromManager<std::vector<double>>(fitMan->raw()["General"]["OscillationParameters"], {});
  osc->setParameters(OscPars);

  std::vector<samplePDFFDBase*> DUNEPdfs;
  //MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);
  MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

  
  for (auto pdf : DUNEPdfs) {
    //std::cout << pdf->GetTitle() << " integral: " << pdf->get1DHist()->Integral() << std::endl;
  }
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
    //std::cout << "[DEBUG] Writing histogram " << Hist1D << " with integral " << Hist1D->Integral() << std::endl;

    Hist1D->Write(DUNEPdfs[iPDF]->GetTitle().c_str());
    if(plotVarspectra)
    {
      for(size_t var_i = 0 ; var_i < var_names.size() ; var_i++)
      {
        Hist1D_Vars[var_i]->Write((DUNEPdfs[iPDF]->GetTitle()+var_names[var_i]).c_str());
      }
    }
  }
  for (auto pdf : DUNEPdfs) {
    //std::cout << pdf->GetTitle() << " integral: " << pdf->get1DHist()->Integral() << std::endl;
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

  // --- print all branch names (for inspection)
std::cout << "[DEBUG] All branches in tree (" << no_of_branches << "):\n";
for (int i = 0; i < no_of_branches; ++i) {
  std::cout << "  " << i << ": " << branch_names[i] << "\n";
}

// ------------------------------------------------------------
// Diagnostic + safety guard for posterior variation
// ------------------------------------------------------------
{
    std::cout << "\n[DEBUG] Checking variation in xsec_* branches...\n";

    int nBranches = mcmc->GetListOfBranches()->GetEntries();
    Long64_t nEntries = mcmc->GetEntries();
    Long64_t nCheckEntries = std::min<Long64_t>(10000, nEntries);

    double val = 0.0;
    int countVarying = 0;
    int countZero = 0;

    for (int i = 0; i < nBranches; ++i) {
        TBranch* br = (TBranch*)mcmc->GetListOfBranches()->At(i);
        TString name = br->GetName();

        if (!name.BeginsWith("xsec_")) continue;

        mcmc->SetBranchAddress(name, &val);

        double minVal = 1e30;
        double maxVal = -1e30;
        double sumVal = 0.0;

        for (Long64_t j = 0; j < nCheckEntries; ++j) {
            mcmc->GetEntry(j);
            if (val < minVal) minVal = val;
            if (val > maxVal) maxVal = val;
            sumVal += val;
        }

        double meanVal = sumVal / nCheckEntries;

        std::cout << std::setw(10) << name
                  << " | min: " << std::setw(12) << minVal
                  << " | max: " << std::setw(12) << maxVal
                  << " | mean: " << std::setw(12) << meanVal << "\n";

        if (fabs(maxVal - minVal) > 1e-12) ++countVarying;
        if (fabs(maxVal) < 1e-12 && fabs(minVal) < 1e-12) ++countZero;
    }

    std::cout << "[DEBUG] Summary: " << countVarying
              << " varying xsec_* branches, "
              << countZero << " constant zero branches.\n";

    if (countVarying == 0) {
        std::cerr << "\n[WARNING] No variation detected in posterior samples!\n"
                     "           All xsec_* branches appear constant.\n"
                     "           Skipping posterior predictive generation to avoid\n"
                     "           producing 1000 identical histograms.\n\n";
        return 0; // <<< Safely exit before the sampling loop >>>
    }

    std::cout << "[DEBUG] Posterior variation detected. Proceeding with sampling.\n\n";
}


// --- 1. Enable all branches ---
mcmc->SetBranchStatus("*", 1);

// --- 2. Step variable ---
int mcmc_step = -1;
mcmc->SetBranchAddress("step", &mcmc_step);

// --- 3. Prepare working xsec arrays ---
std::vector<double> xsec_nominal = xsec->getNominalArray();
std::vector<double> xsec_draw(xsec_nominal.size(), 0.0);  // for PDFs

// Temporary storage for ROOT binding
std::vector<Double_t> xsec_tmp(xsec_nominal.size(), 0.0);

// --- 4. Bind each xsec_* branch individually ---
int no_of_xsec_branches = 0;
for (size_t i = 0; i < xsec_nominal.size(); ++i) {
    TString bname = Form("xsec_%zu", i);

    // Check branch exists
    if (mcmc->GetBranch(bname)) {
        mcmc->SetBranchAddress(bname, &xsec_tmp[i]);
        ++no_of_xsec_branches;
        std::cout << "[DEBUG] Bound branch " << bname << " -> xsec_tmp[" << i << "]\n";
    } else {
        std::cerr << "[WARN] Branch " << bname << " does not exist in tree\n";
    }
}

std::cout << "[DEBUG] Number of xsec branches bound: " << no_of_xsec_branches << "\n";

// --- 5. Test one entry ---
mcmc->GetEntry(0);

// Copy to working array for PDFs
for (size_t i = 0; i < xsec_draw.size(); ++i) xsec_draw[i] = xsec_tmp[i];

/*
std::cout << "[DEBUG] After GetEntry(0): mcmc_step=" << mcmc_step
          << " | xsec_draw[0..4]="
          << xsec_draw[0] << " " << xsec_draw[1] << " "
          << xsec_draw[2] << " " << xsec_draw[3] << " " << xsec_draw[4]
          << std::endl;

// --- Test one entry ---
mcmc->GetEntry(0);
std::cout << "[DEBUG] After GetEntry(0): mcmc_step=" << mcmc_step
          << " | xsec_draw[0..4]="
          << xsec_draw[0] << " " << xsec_draw[1] << " "
          << xsec_draw[2] << " " << xsec_draw[3] << " " << xsec_draw[4]
          << std::endl;*/

  // ===============================================================
// Detect which xsec parameters actually vary in the posterior
// ===============================================================
std::cout << "[DEBUG] Scanning xsec parameter ranges to identify active parameters...\n";

std::vector<double> xsec_min(xsec_draw.size(),  1e9);
std::vector<double> xsec_max(xsec_draw.size(), -1e9);

// Quick scan over the MCMC chain to find min/max per parameter
Long64_t n_scan = std::min<Long64_t>(5000, mcmc->GetEntries());
for (Long64_t i = 0; i < n_scan; ++i) {
    mcmc->GetEntry(i);
    for (size_t j = 0; j < xsec_tmp.size(); ++j) {
        xsec_min[j] = std::min(xsec_min[j], xsec_tmp[j]);
        xsec_max[j] = std::max(xsec_max[j], xsec_tmp[j]);
    }
}

// Detect active (non-flat) parameters
constexpr double VARIATION_THRESHOLD = 1e-6;
std::vector<int> active_xsecs;
for (size_t j = 0; j < xsec_draw.size(); ++j) {
    double range = std::fabs(xsec_max[j] - xsec_min[j]);
    if (range > VARIATION_THRESHOLD && !(xsec_min[j] == 0.0 && xsec_max[j] == 0.0)) {
        active_xsecs.push_back(static_cast<int>(j));
    }
}

std::cout << "[INFO] Active xsec parameters detected: " << active_xsecs.size()
          << " out of " << xsec_draw.size() << std::endl;

if (active_xsecs.empty()) {
    std::cerr << "[WARNING] No active xsec parameters found! "
              << "All reweighting will use nominal values.\n";
}

// ===============================================================
// End of detection block
// ===============================================================



  /////////// 5. Sample steps with burn-in check
  int nEntries = mcmc->GetEntries();
  auto rnd = std::make_unique<TRandom3>(0);
  std::vector<int> mcmc_draw_entries(no_times_sampling_posterior);

  // Quick check: read first entry (or a few entries) and print xsec_draw contents
if (nEntries > 0) {
  mcmc->GetEntry(0);
 // std::cout << "[DEBUG] After GetEntry(0): mcmc_step=" << mcmc_step << " ; xsec_draw[0..4] = ";
  for (int k = 0; k < std::min((size_t)5, xsec_draw.size()); ++k) std::cout << xsec_draw[k] << " ";
  std::cout << "\n";
} else {
  std::cerr << "[ERROR] MCMC tree has zero entries!\n";
}

// If you have root interactive access, run:
// mcmc->Scan("xsec_0:xsec_1", "", "colsize=15")


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
    // ----------- 6. Reweight and write predictions -------------
fOut->cd();

for (int i = 0; i < no_times_sampling_posterior; ++i) {
    if (i % (no_times_sampling_posterior / 10) == 0) {
        MaCh3Utils::PrintProgressBar(i, no_times_sampling_posterior);
    }

    // --- Draw a random MCMC entry
    int entry = -1;
    do {
        entry = rnd->Integer(nEntries);
        mcmc->GetEntry(entry);   // fills xsec_tmp and mcmc_step
    } while ((unsigned int)mcmc_step < burn_in);

    /*
    // --- Update working xsec_draw array from ROOT branch values
    for (size_t j = 0; j < xsec_draw.size(); ++j) {
        xsec_draw[j] = xsec_tmp[j];
    }*/

    //test to see if making the parameter variations bigger gives a noticable sized error band
    double inflateFactor = 100.0; 
    for (size_t j = 0; j < xsec_draw.size(); ++j) {
        double nominal = xsec_nominal[j];
        double delta = xsec_tmp[j] - nominal;
        xsec_draw[j] = nominal + inflateFactor * delta;
    }
    xsec->setParameters(xsec_draw);


    
    /*std::cout << "[DEBUG] Sample " << i << " from entry " << entry
              << " | Step: " << mcmc_step
              << " | xsec[0..4]: "
              << xsec_draw[0] << " " << xsec_draw[1] << " "
              << xsec_draw[2] << " " << xsec_draw[3] << " "
              << xsec_draw[4] << "\n";*/

    // --- Apply sample parameters to xsec object
    xsec->setParameters(xsec_draw);
    osc->setParameters(OscPars);

    std::cout << "[DEBUG] First few active xsecs before reweight: ";
    for (int k = 0; k < 5; ++k)
        std::cout << xsec_tmp[active_xsecs[k]] << " ";
    std::cout << std::endl;


    // --- Reweight PDFs and write histograms
    for (size_t j = 0; j < DUNEPdfs.size(); ++j) {
        if (!DUNEPdfs[j]) continue;

        // Reweight using current sample
        DUNEPdfs[j]->reweight();

        // Get 1D projection
        TH1D* h = DUNEPdfs[j]->get1DHist();
        if (!h) continue;

        TString histName = Form("%s_posterior_sample_%04d", 
                                DUNEPdfs[j]->GetTitle(), i);
        TH1D* writeHist = static_cast<TH1D*>(h->Clone(histName));
        writeHist->SetDirectory(fOut);
        writeHist->Write();
        delete writeHist;

        
        std::cout << "[DEBUG] PDF " << j 
                  << " after reweight, integral: " << h->Integral() << "\n";
    }
}


    // inside the sample loop after xsec->setParameters(xsec_draw) and DUNEPdfs[j]->reweight()
    const std::vector<KinematicCut> emptySelectionVec;
    std::string var = "RecoNeutrinoEnergy";

    for (size_t j = 0; j < DUNEPdfs.size(); ++j) {
        if (!DUNEPdfs[j]) continue;

    // reweight once
    DUNEPdfs[j]->reweight();

    // 1) get the 1D projected RecoNeutrinoEnergy histogram for this pdf
    // Try the overload that projects variable: check both signatures to be safe
    TH1D* varHist = nullptr;
    // If you have an overload that takes selection & axis:
    TAxis tempAxis; // optional if you need to override binning
    // tempAxis.Set(nbins, xmin, xmax); // only if needed
    // varHist = static_cast<TH1D*>(DUNEPdfs[j]->get1DVarHist(var, emptySelectionVec, 0, &tempAxis));

    // Fallback to the simpler overload if available
    if (!varHist) varHist = static_cast<TH1D*>(DUNEPdfs[j]->get1DVarHist(var));

    if (!varHist) {
      std::cerr << "[WARN] get1DVarHist returned nullptr for var '" << var
                << "' pdf " << j << " (title=" << DUNEPdfs[j]->GetTitle() << ")\n";
      continue;
    }

    // Clone so we have a uniquely-named histogram to write
    TString sampleName = DUNEPdfs[j]->GetTitle();
    TString outName = Form("%s_%s_posterior_sample_%04d_pdf%zu",
                           sampleName.Data(), var.c_str(), i, j);
    TH1D* writeHist = static_cast<TH1D*>(varHist->Clone(outName));
    writeHist->SetDirectory(fOut); // keep in output file
    writeHist->Write();
    delete writeHist; // safe because written to file
}

    //const std::vector<KinematicCut> emptySelectionVec;
    //std::string var = "RecoNeutrinoEnergy";
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

                TString Varwritename = Form("%s_%s_%d", DUNEPdfs[iPDF]->GetTitle().c_str(), var.c_str(), (int)var_i);


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
