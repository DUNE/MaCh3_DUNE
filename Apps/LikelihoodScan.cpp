#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>

#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Fitters/MaCh3Factory.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include "TIterator.h"
#include <fstream>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include <fstream>
#include <iomanip>

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TKey.h"
#include <fstream>
#include <iomanip>

void DumpDirRecursive(TDirectory* dir, const std::string& path, std::ofstream& out) {
  TIter next(dir->GetListOfKeys());
  TKey* key;
  while ((key = (TKey*)next())) {
    TObject* obj = key->ReadObj();
    std::string fullName = path + "/" + obj->GetName();

    if (obj->InheritsFrom(TDirectory::Class())) {
      DumpDirRecursive((TDirectory*)obj, fullName, out);
    }
    else if (obj->InheritsFrom(TH2D::Class())) {
      TH2D* h2 = (TH2D*)obj;
      out << "# " << fullName << "\n# param_x\tparam_y\tlogL\n";
      for (int i = 1; i <= h2->GetNbinsX(); ++i)
        for (int j = 1; j <= h2->GetNbinsY(); ++j)
          out << h2->GetXaxis()->GetBinCenter(i) << "\t"
              << h2->GetYaxis()->GetBinCenter(j) << "\t"
              << h2->GetBinContent(i, j) << "\n";
      out << "\n";
    }
    else if (obj->InheritsFrom(TH1D::Class())) {
      TH1D* h = (TH1D*)obj;
      out << "# " << fullName << "\n# param_value\tlogL\n";
      for (int i = 1; i <= h->GetNbinsX(); ++i)
        out << h->GetBinCenter(i) << "\t" << h->GetBinContent(i) << "\n";
      out << "\n";
    }
  }
}

void DumpLLHScanToTxt(const std::string& rootFileName,
                       const std::string& txtFileName) {
  TFile* f = TFile::Open(rootFileName.c_str(), "READ");
  if (!f || f->IsZombie()) {
    MACH3LOG_ERROR("Could not open {}", rootFileName);
    return;
  }

  std::ofstream out(txtFileName);
  out << std::scientific << std::setprecision(8);

  DumpDirRecursive(f, "", out);

  out.close();
  f->Close();
  MACH3LOG_INFO("Wrote LLH scan values to {}", txtFileName);
}
int main(int argc, char * argv[]) {
  auto FitManager = MaCh3ManagerFactory(argc, argv);
  {
    TFile precreate("llhfsdcov.root", "RECREATE");
    precreate.Close();
  }
  // 1D scan on by default, and 2D off
  const bool do_1d_llhscan = GetFromManager(FitManager->raw()["General"]["1DLLHScan"], true);
  const bool do_2d_llhscan = GetFromManager(FitManager->raw()["General"]["2DLLHScan"], false);

  if (!do_1d_llhscan && !do_2d_llhscan) {
    MACH3LOG_ERROR("Neither 1D or 2D llhscan enabled");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  //###############################################################################################################################
  //Create samplePDFFD objects
  
  ParameterHandlerGeneric* xsec = nullptr;
  
  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  //###############################################################################################################################
  //Perform reweight, print total integral, and set data
  int idx = xsec->GetParIndex("MissingProtonFD");

  xsec->SetSingleParameter(idx, 0.2);
  xsec->SetFixParameter(idx);
  std::vector<TH1*> DUNEHists;
  for(auto handler : DUNEPdfs){
    for (unsigned iSample = 0; iSample < handler->GetNsamples(); ++iSample) {
      handler->Reweight();
      DUNEHists.push_back(handler->GetMCHist(iSample));
      MACH3LOG_INFO("Event rate for {} : {:<5.2f}", handler->GetSampleTitle(iSample), handler->GetMCHist(iSample)->Integral());

      if (handler->GetNDim(iSample) == 1) {
        handler->AddData(iSample, (TH1D*)DUNEHists.back());
      } else if (handler->GetNDim(iSample) == 2) {
        handler->AddData(iSample, (TH2D*)DUNEHists.back());
      }
    }
  }
  auto MaCh3Fitter = MaCh3FitterFactory(FitManager.get());
  xsec->SetSingleParameter(idx, 0.0);
  //###############################################################################################################################
  //Lets benefit from the core code utilities 
  
  //Add samples to FitterBase
  for(auto Sample : DUNEPdfs){
    MaCh3Fitter->AddSampleHandler(Sample);
  }

  MaCh3Fitter->AddSystObj(xsec);
  
 if (do_1d_llhscan) {
  MaCh3Fitter->RunLLHScan();
  MaCh3Fitter->GetStepScaleBasedOnLLHScan();
}
if (do_2d_llhscan) {
  MaCh3Fitter->Run2DLLHScan();
}

MaCh3Fitter.reset();  // flush & close twstcov.root

DumpLLHScanToTxt("twstcov.root", "llh_scan_output.txt");
}
