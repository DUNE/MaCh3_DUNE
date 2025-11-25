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

#include "mcmc/mcmc.h"
#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"

int main(int argc, char * argv[]) {
  if(argc == 1){
    MACH3LOG_ERROR("Usage: bin/EventRatesDUNEBeam config.cfg");
    return 1;
  }
  manager* FitManager = new manager(argv[1]);

  //###############################################################################################################################

  //DB Sigma variations in units of each parameters Sigma
  std::vector<double> sigmaVariations = {-3,-2, -1, 0, 1,2, 3};

  //###############################################################################################################################
  //Create samplePDFFD objects
  
  covarianceXsec* xsec = nullptr;
  covarianceOsc* osc = nullptr;

  std::vector<samplePDFFDBase*> DUNEPdfs;
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec, osc);

  //###############################################################################################################################
  //Perform reweight and print total integral

  MACH3LOG_INFO("=======================================================");
  for(samplePDFFDBase* Sample: DUNEPdfs){
    Sample->reweight();
    MACH3LOG_INFO("Event rate for {} : {:<5.2f}", Sample->GetTitle(), Sample->get1DHist()->Integral());
  }
  
  //###############################################################################################################################
  //DB Can't use the core sigma variations as it's entirely set up around the concept of multiple selections per samplePDF object
  //   Thats not the case in the FD code, which has one selection per samplePDF object
  //   Consequently have to write out own code
  
  std::vector<covarianceBase*> CovObjs;
  CovObjs.emplace_back(xsec);
  //CovObjs.emplace_back(osc);

  MACH3LOG_INFO("=======================================================");

  std::string OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();
  TFile* File = TFile::Open(OutputFileName.c_str(),"RECREATE");
  
  /*
  for (covarianceBase* CovObj: CovObjs) {
    MACH3LOG_INFO("Starting Variations for covarianceBase Object: {}",CovObj->getName());
    
    int nPars = CovObj->getNpars();
    for (int iPar=0;iPar<nPars;iPar++) {
      std::string ParName = CovObj->GetParFancyName(iPar);
      double VarInit = CovObj->getParInit(iPar);
      double VarSigma = CovObj->getDiagonalError(iPar);
      
      MACH3LOG_INFO("\tParameter : {:<30} - Variations around value : {:<10.7f} , in units of 1 Sigma : {:<10.7f}",ParName,VarInit,VarSigma);

      File->cd();
      File->mkdir(ParName.c_str());
      File->cd(ParName.c_str());
      
      for (size_t iSigVar=0;iSigVar<sigmaVariations.size();iSigVar++) {
	double VarVal = VarInit + sigmaVariations[iSigVar]*VarSigma;
	if (VarVal < CovObj->GetLowerBound(iPar)) VarVal = CovObj->GetLowerBound(iPar);
	if (VarVal > CovObj->GetUpperBound(iPar)) VarVal = CovObj->GetUpperBound(iPar);
	
	MACH3LOG_INFO("\t\tVariation {:<5.3f} - Parameter Value : {:<10.7f}",sigmaVariations[iSigVar],VarVal);
	CovObj->setParProp(iPar,VarVal);

	for (size_t iSample=0;iSample<DUNEPdfs.size();iSample++) {
	  std::string SampleName = DUNEPdfs[iSample]->GetTitle();
	  
	  File->cd(ParName.c_str());
	  if (iSigVar == 0) {
	    File->mkdir((ParName+"/"+SampleName).c_str());
	  }
	  File->cd((ParName+"/"+SampleName).c_str());
	  
	  DUNEPdfs[iSample]->reweight();
	  TH1* Hist;
    if (DUNEPdfs[iSample]->GetNDim() == 1)
      Hist = DUNEPdfs[iSample]->get1DHist();
    else if (DUNEPdfs[iSample]->GetNDim() == 2)
      Hist = DUNEPdfs[iSample]->get2DHist();
	  MACH3LOG_INFO("\t\t\tSample : {:<30} - Integral : {:<10}",SampleName,Hist->Integral());
	  
	  Hist->Write(Form("Variation_%i",(int)iSigVar));
	}
      }

      CovObj->setParProp(iPar,VarInit);
    }

    MACH3LOG_INFO("=======================================================");
  }*/
  std::string pdfFileName = "SigmaVariations.pdf";
bool firstCanvas = true;

for (covarianceBase* CovObj: CovObjs) {
    int nPars = CovObj->getNpars();
    for (int iPar = 0; iPar < nPars; iPar++) {
        std::string ParName = CovObj->GetParFancyName(iPar);
        double VarInit = CovObj->getParInit(iPar);
        double VarSigma = CovObj->getDiagonalError(iPar);

        File->cd();
        File->mkdir(ParName.c_str());
        File->cd(ParName.c_str());

        // Canvas per parameter
        TCanvas* cPar = new TCanvas(Form("c_%s", ParName.c_str()), Form("Systematic :  %s", ParName.c_str()), 1000, 700);
        cPar->cd();

        TLegend* legend = new TLegend(0.7,0.7,0.9,0.9);
        THStack* hsAll = new THStack("hsAll", Form("Systematic :  %s", ParName.c_str()));

       
      std::vector<int> lineColors = {
          kBlue+4, // -3, lightest
          kBlue+2, 
          kBlue,   
          kBlack,  // 0σ (nominal)
          kBlue,   
          kBlue+2, 
          kBlue+4  // +3, lightest
      };

      int colorIndex = 0;


        for (size_t iSample = 0; iSample < DUNEPdfs.size(); iSample++) {
            std::string SampleName = DUNEPdfs[iSample]->GetTitle();
            File->cd((ParName + "/" + SampleName).c_str());

            for (size_t iSigVar = 0; iSigVar < sigmaVariations.size(); iSigVar++) {
                double VarVal = VarInit + sigmaVariations[iSigVar] * VarSigma;
                if (VarVal < CovObj->GetLowerBound(iPar)) VarVal = CovObj->GetLowerBound(iPar);
                if (VarVal > CovObj->GetUpperBound(iPar)) VarVal = CovObj->GetUpperBound(iPar);
                CovObj->setParProp(iPar, VarVal);

                DUNEPdfs[iSample]->reweight();
                if (DUNEPdfs[iSample]->GetNDim() != 1) continue;

                //double w = DUNEPdfs[iSample]->reweight();
                //std::cout << "Weight for " << sigmaVariations[iSigVar] << " sigma = " << w << std::endl;

                TH1* Hist = (TH1*)DUNEPdfs[iSample]->get1DHist()->Clone(Form("Hist_%s_%i", SampleName.c_str(), (int)iSigVar));
                int colorIndex = static_cast<int>(iSigVar); // iSigVar = 0..6 if you redefine sigmaVariations as {-3,-2,-1,0,1,2,3}

                // Set line color
                Hist->SetLineColor(lineColors[colorIndex]);
                Hist->SetLineWidth(2);
                  
            
                std::cout << "sigmaVariations[iSigVar] = " << sigmaVariations[iSigVar] << "  VarVal =  " <<  VarVal << std::endl;               // Label axes
                Hist->GetXaxis()->SetTitle("E^{Reco.}_{#nu} [GeV]");
                Hist->GetYaxis()->SetTitle("Event Rate");
                // Compute color index
           
                hsAll->Add(Hist);
                legend->AddEntry(Hist, Form("%s %+g #sigma", SampleName.c_str(), sigmaVariations[iSigVar]), "l");

                // --- Print bin contents for debugging ---
                std::cout << "=== " << ParName 
                          << " | Sample: " << SampleName 
                          << " | Variation: " << sigmaVariations[iSigVar] 
                          << " sigma (VarVal=" << VarVal << ") ===" 
                          << std::endl;

                for (int b = 1; b <= Hist->GetNbinsX(); b++) {
                    std::cout << "  Bin " << b 
                              << "  x=" << Hist->GetBinCenter(b)
                              << "  content=" << Hist->GetBinContent(b)
                              << std::endl;
                }
                std::cout << std::endl;


                Hist->Write(Form("Variation_%i", (int)iSigVar));
            }
        }

        hsAll->Draw("NOSTACK HIST");
        legend->Draw();

        // Save to single PDF
        if (firstCanvas)
            cPar->SaveAs(Form("%s(", pdfFileName.c_str())); // Open PDF
        else
            cPar->SaveAs(pdfFileName.c_str());              // Append page
        firstCanvas = false;

        cPar->Write();
        delete cPar;

        CovObj->setParProp(iPar, VarInit);
        std::cout << "DEBUG param " << ParName 
          << " internal = " << CovObj->getParProp(iPar) 
          << std::endl;

    }
}

// Close the PDF
TCanvas* cEnd = new TCanvas();
cEnd->SaveAs(Form("%s)", pdfFileName.c_str()));
delete cEnd;




  File->Close();
  
}
