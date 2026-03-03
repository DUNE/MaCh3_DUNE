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
#include <TH3.h>
#include "TError.h"

#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Fitters/MaCh3Factory.h"


void PlotSigmaVariations(const std::string& filename,
                         const std::vector<double>& sigmaVariations)
{
  TFile* file = TFile::Open(filename.c_str(), "READ");
  if (!file || file->IsZombie()) {
    MACH3LOG_ERROR("Could not open {}", filename);
    return;
  }

  auto canvas = std::make_unique<TCanvas>("c", "c", 1000, 800);
  canvas->SetGrid();

  std::string pdfName = filename.substr(0, filename.find(".root")) + "_SigmaVar.pdf";
  canvas->Print((pdfName + "[").c_str());

  TIter parIter(file->GetListOfKeys());
  TKey* parKey = nullptr;

  while ((parKey = static_cast<TKey*>(parIter()))) {
    if (std::string(parKey->GetClassName()) != "TDirectoryFile") continue;

    TDirectory* parDir = file->GetDirectory(parKey->GetName());
    std::string parName = parKey->GetName();

    TIter sampIter(parDir->GetListOfKeys());
    TKey* sampKey = nullptr;

    while ((sampKey = static_cast<TKey*>(sampIter()))) {
      if (std::string(sampKey->GetClassName()) != "TDirectoryFile") continue;

      TDirectory* sampDir =
          parDir->GetDirectory(sampKey->GetName());
      std::string sampName = sampKey->GetName();

      std::vector<TH1*> hists;
      double maxY = 0;

      for (size_t i = 0; i < sigmaVariations.size(); ++i) {
        auto* h = sampDir->Get<TH1>(Form("Variation_%zu", i));
        if (!h) continue;

        h = static_cast<TH1*>(h->Clone());
        h->SetDirectory(nullptr);

        h->SetLineWidth(2);
        h->SetLineColor(kBlack + i);
        h->SetLineStyle(i == 3 ? kSolid : kDashed);

        maxY = std::max(maxY, h->GetMaximum());
        hists.push_back(h);
      }

      if (hists.empty()) continue;

      canvas->Clear();
      hists[0]->SetTitle(Form("%s – %s", parName.c_str(), sampName.c_str()));
      hists[0]->SetMaximum(maxY * 1.2);
      hists[0]->Draw("HIST");

      for (size_t i = 1; i < hists.size(); ++i)
        hists[i]->Draw("HIST SAME");

      auto leg = std::make_unique<TLegend>(0.6, 0.6, 0.88, 0.88);
      for (size_t i = 0; i < hists.size(); ++i)
        leg->AddEntry(hists[i],
                      Form("%+.0f#sigma", sigmaVariations[i]),
                      "l");

      leg->Draw();
      canvas->Print(pdfName.c_str());

      for (auto* h : hists) delete h;
    }
  }

  canvas->Print((pdfName + "]").c_str());
  file->Close();
}


int main(int argc, char * argv[]) {
  gErrorIgnoreLevel = kFatal;
  auto FitManager = MaCh3ManagerFactory(argc, argv);

  //###############################################################################################################################

  //DB Sigma variations in units of each parameters Sigma
  std::vector<double> sigmaVariations = {-3,-2,-1,0, 1,2,3};

  //###############################################################################################################################
  //Create samplePDFFD objects
  
  ParameterHandlerGeneric* xsec = nullptr;

  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  //###############################################################################################################################
  //Perform reweight and print total integral

  MACH3LOG_INFO("=======================================================");
  for(SampleHandlerFD* Sample: DUNEPdfs){
    Sample->Reweight();
    MACH3LOG_INFO("Event rate for {} : {:<5.2f}", Sample->GetTitle(), Sample->GetMCHist(Sample->GetNDim())->Integral());
  }
  
  //###############################################################################################################################
  //DB Can't use the core sigma variations as it's entirely set up around the concept of multiple selections per samplePDF object
  //   Thats not the case in the FD code, which has one selection per samplePDF object
  //   Consequently have to write out own code
  
  std::vector<ParameterHandlerBase*> CovObjs;
  CovObjs.emplace_back(xsec);

  MACH3LOG_INFO("=======================================================");

  std::string OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();
  TFile* File = TFile::Open(OutputFileName.c_str(),"RECREATE");

  for (ParameterHandlerBase* CovObj : CovObjs) {
    MACH3LOG_INFO("Starting Variations for covarianceBase Object: {}", CovObj->GetName());

    int nPars = CovObj->GetNumParams();
    for (int iPar = 0; iPar < nPars; iPar++) {
        std::string ParName = CovObj->GetParFancyName(iPar);
        double VarInit = CovObj->GetParInit(iPar);
        double VarSigma = CovObj->GetDiagonalError(iPar);

        MACH3LOG_INFO("\tParameter : {:<30} - Variations around value : {:<10.7f} , in units of 1 Sigma : {:<10.7f}",
                      ParName, VarInit, VarSigma);

        File->cd();
        File->mkdir(ParName.c_str());
        File->cd(ParName.c_str());

        for (size_t iSigVar = 0; iSigVar < sigmaVariations.size(); iSigVar++) {
            double VarVal = VarInit + sigmaVariations[iSigVar] * VarSigma;
            if (VarVal < CovObj->GetLowerBound(iPar)) VarVal = CovObj->GetLowerBound(iPar);
            if (VarVal > CovObj->GetUpperBound(iPar)) VarVal = CovObj->GetUpperBound(iPar);

            MACH3LOG_INFO("\t\tVariation {:<5.3f} - Parameter Value : {:<10.7f}", sigmaVariations[iSigVar], VarVal);
            CovObj->SetParProp(iPar, VarVal);

            for (size_t iSample = 0; iSample < DUNEPdfs.size(); iSample++) {
                std::string SampleName = DUNEPdfs[iSample]->GetTitle();

                File->cd(ParName.c_str());
                if (iSigVar == 0) {
                    File->mkdir((ParName + "/" + SampleName).c_str());
                }
                File->cd((ParName + "/" + SampleName).c_str());

                DUNEPdfs[iSample]->Reweight();
                TH1* HistOrig = DUNEPdfs[iSample]->GetMCHist(DUNEPdfs[iSample]->GetNDim());
                TH1* HistToWrite = nullptr;

                // Project to X if 2D or 3D
                if (HistOrig->InheritsFrom("TH3")) {
                    HistToWrite = static_cast<TH3*>(HistOrig)->ProjectionX(Form("%s_projX", HistOrig->GetName()));
                } else if (HistOrig->InheritsFrom("TH2")) {
                    HistToWrite = static_cast<TH2*>(HistOrig)->ProjectionX(Form("%s_projX", HistOrig->GetName()));
                } else {
                    HistToWrite = static_cast<TH1*>(HistOrig->Clone());
                    HistToWrite->SetDirectory(nullptr);
                }

                MACH3LOG_INFO("\t\t\tSample : {:<30} - Integral : {:<10}", SampleName, HistToWrite->Integral());

                HistToWrite->Write(Form("Variation_%i", (int)iSigVar));
                delete HistToWrite;  // Prevent memory leak
            }
        }

        // Restore original parameter value
        CovObj->SetParProp(iPar, VarInit);
    }

    MACH3LOG_INFO("=======================================================");
}
  
  // for (ParameterHandlerBase* CovObj: CovObjs) {
  //   MACH3LOG_INFO("Starting Variations for covarianceBase Object: {}",CovObj->GetName());
    
  //   int nPars = CovObj->GetNumParams();
  //   for (int iPar=0;iPar<nPars;iPar++) {
  //     std::string ParName = CovObj->GetParFancyName(iPar);
  //     double VarInit = CovObj->GetParInit(iPar);
  //     double VarSigma = CovObj->GetDiagonalError(iPar);
      
  //     MACH3LOG_INFO("\tParameter : {:<30} - Variations around value : {:<10.7f} , in units of 1 Sigma : {:<10.7f}",ParName,VarInit,VarSigma);

  //     File->cd();
  //     File->mkdir(ParName.c_str());
  //     File->cd(ParName.c_str());
      
  //     for (size_t iSigVar=0;iSigVar<sigmaVariations.size();iSigVar++) {
	// double VarVal = VarInit + sigmaVariations[iSigVar]*VarSigma;
	// if (VarVal < CovObj->GetLowerBound(iPar)) VarVal = CovObj->GetLowerBound(iPar);
	// if (VarVal > CovObj->GetUpperBound(iPar)) VarVal = CovObj->GetUpperBound(iPar);
	
	// MACH3LOG_INFO("\t\tVariation {:<5.3f} - Parameter Value : {:<10.7f}",sigmaVariations[iSigVar],VarVal);
	// CovObj->SetParProp(iPar,VarVal);

	// for (size_t iSample=0;iSample<DUNEPdfs.size();iSample++) {
	//   std::string SampleName = DUNEPdfs[iSample]->GetTitle();
	  
	//   File->cd(ParName.c_str());
	//   if (iSigVar == 0) {
	//     File->mkdir((ParName+"/"+SampleName).c_str());
	//   }
	//   File->cd((ParName+"/"+SampleName).c_str());
	  
	//   DUNEPdfs[iSample]->Reweight();
	//   TH1* Hist = DUNEPdfs[iSample]->GetMCHist(DUNEPdfs[iSample]->GetNDim());
	//   MACH3LOG_INFO("\t\t\tSample : {:<30} - Integral : {:<10}",SampleName,Hist->Integral());
	  
	//   Hist->Write(Form("Variation_%i",(int)iSigVar));
	// }
  //     }

  //     CovObj->SetParProp(iPar,VarInit);
  //   }

  //   MACH3LOG_INFO("=======================================================");
  // }

  File->Close();
  

  PlotSigmaVariations(OutputFileName, sigmaVariations);

  
}
