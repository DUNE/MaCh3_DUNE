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
#include <TError.h>

#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Fitters/MaCh3Factory.h"

#include <TFile.h>
#include <TMatrixDSym.h>
#include <TMinuit.h>


struct MyChi2Blob {

    std::vector<SampleHandlerBase*> samples;
    ParameterHandlerGeneric* systematics;

    // *******************
    double CalcChi2(const double* x) {
    // *******************

        int ParCounter = 0;
        double llh = 0;

        std::vector<double> pars;
        const int NumPar = systematics->GetNumParams();
        //KS: Avoid push back as they are slow
        pars.resize(NumPar);
        for(int i = 0; i < NumPar; ++i, ++ParCounter)
        {
            double ParVal = x[ParCounter];

            pars[i] = ParVal;
            //std::cout << "Param[" << i << "] = " << pars[i] << std::endl;
        }


        static int call = 0;
        if (call < 5) {
          std::cout << "Minuit call " << call << "\n";
          for (int i = 0; i < 5; ++i) {
            std::cout << "x[" << i << "] = " << x[i] << "\n";
          }
        }
        call++;

        
        systematics->SetParameters(pars);
        // systematics->AcceptStep();
        //systematics->SetParCurrProp(pars);

        
        llh += systematics->CalcLikelihood();
        
        //std::cout<< "llh += systematics->CalcLikelihood() = " << llh << std::endl;
       
        //systematics->PrintNominalCurrProp();
        // Could multi-thread this
        // But since sample reweight is multi-threaded it's probably better to do that
    
        for (size_t i = 0; i < samples.size(); i++)
        {
            samples[i]->Reweight();
            
        }

        //DB for atmospheric event by event sample migration, need to fully reweight all samples to allow event passing prior to likelihood evaluation
        for (size_t i = 0; i < samples.size(); i++) {
            // Get the sample likelihoods and add them
            llh += samples[i]->GetLikelihood();
            //std::cout<< "samples[i]->GetLikelihood() = " << llh << std::endl;

        }
        //std::cout<< "llh = " << llh << std::endl;
        llh = 2.0*llh;
        //std::cout<< "llh = " << llh << std::endl;
        return llh;
    }

};

int main(int argc, char * argv[]) {
  gErrorIgnoreLevel = kFatal;
  auto FitManager = MaCh3ManagerFactory(argc, argv);

  //###############################################################################################################################
  //Create samplePDFFD objects

  ParameterHandlerGeneric* xsec = nullptr;
  std::vector<SampleHandlerFD*> DUNEPdfs;
  std::vector<TH1*> DUNEHists;//akp
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  MyChi2Blob myblob;
   myblob.systematics = xsec;// dont need mayber?
    for(auto Sample : DUNEPdfs){
      double weight = 1.0;
      //xsec->SetParCurrProp(0, weight); nuwro
        
      //Reweight
      Sample->Reweight();
      //Fix it
      DUNEHists.push_back(Sample->GetMCHist(Sample->GetNDim()));
      //xsec->SetParCurrProp(0, 0.0); nuwro
      //xsec->ToggleFixParameter(0); nuwro
      if (Sample->GetNDim() == 1) {
          Sample->AddData((TH1D*)DUNEHists.back());
      } else if (Sample->GetNDim() == 2) {
          Sample->AddData((TH2D*)DUNEHists.back());
      }
      myblob.samples.push_back(Sample);
      
      
    }

   
  std::cout << "Total events in histogram = " << DUNEHists.back()->Integral() << std::endl;

    
  
  std::cout << "Number of PDFs: " << DUNEPdfs.size() << "\n";
  if(!DUNEPdfs.empty()) {
      std::cout << "First PDF pointer: " << DUNEPdfs[0] << "\n";
  }

    {
      std::cout << "=== PRIOR CURVATURE TEST ===\n";

      std::vector<double> pars(xsec->GetNumParams());
      for (int i = 0; i < xsec->GetNumParams(); ++i)
          pars[i] = xsec->GetParInit(i);

      xsec->SetParameters(pars);
      double llh0 = 2.0 * xsec->CalcLikelihood();

      pars[0] += 0.5;  // large shift
      xsec->SetParameters(pars);
      double llh1 = 2.0 * xsec->CalcLikelihood();

      std::cout << "Prior chi2 nominal = " << llh0 << "\n";
      std::cout << "Prior chi2 shifted = " << llh1 << "\n";
  }

  {
      std::cout << "=== SAMPLE SENSITIVITY TEST ===\n";

      std::vector<double> pars(xsec->GetNumParams());
      for (int i = 0; i < xsec->GetNumParams(); ++i)
          pars[i] = xsec->GetParInit(i);

      xsec->SetParameters(pars);
      for (auto* s : myblob.samples)
          s->Reweight();

      double llh0 = myblob.samples[0]->GetLikelihood();

      pars[0] += 0.5;
      xsec->SetParameters(pars);
      for (auto* s : myblob.samples)
          s->Reweight();

      double llh1 = myblob.samples[0]->GetLikelihood();

      std::cout << "Sample LLH nominal = " << llh0 << "\n";
      std::cout << "Sample LLH shifted = " << llh1 << "\n";
  }



  auto minuit = std::unique_ptr<ROOT::Math::Minimizer>(
  ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));

    //KS: For none PCA this will be equal to normal parameters
  const int NparsMinuitFull = xsec->GetNumParams(); //NPars;

  ROOT::Math::Functor fChi2(&myblob, &MyChi2Blob::CalcChi2, xsec->GetNumParams());
  minuit->SetFunction(fChi2);
  minuit->SetStrategy(2);

  //KS: add config or something
  minuit->SetPrintLevel(4);
  
  minuit->SetMaxFunctionCalls(FitManager->raw()["General"]["Minuit2"]["NSteps"].as<unsigned>());
  minuit->SetMaxIterations(100000);
  double tolerance = FitManager->raw()["General"]["Minuit2"]["ToleranceLevel"].as<double>();
  minuit->SetTolerance(tolerance);
  MACH3LOG_INFO("Preparing Minuit");
  int ParCounter = 0;

  for (int i = 0; i < xsec->GetNumParams(); ++i) {
    std::cout << i << " "
              << xsec->GetParName(i)
              << " error = " << xsec->GetDiagonalError(i)
              << " fixed = " << xsec->IsParameterFixed(i)
              << "\n";
}




      for(int i = 0; i < xsec->GetNumParams(); ++i, ++ParCounter)
      {
        //KS: Index, name, prior, step scale [different to MCMC],
        minuit->SetVariable(ParCounter, (xsec->GetParName(i)), xsec->GetParInit(i), xsec->GetDiagonalError(i)/100);
        minuit->SetVariableValue(ParCounter, xsec->GetParInit(i));
        //KS: lower bound, upper bound, if Mirroring enabled then ignore
        minuit->SetVariableLimits(ParCounter, xsec->GetLowerBound(i), xsec->GetUpperBound(i));
        if(xsec->IsParameterFixed(i))
        { std::cout<< "Found fixed parameter..." << std::endl;
          minuit->FixVariable(ParCounter);
        }
      }
    
    MACH3LOG_INFO("Starting MIGRAD");
    minuit->Minimize(); //nuwro

    MACH3LOG_INFO("Starting HESSE");
    minuit->Hesse();

    std::cout << "Status: " << minuit->Status() << "\n";
    std::cout << "CovStatus: " << minuit->CovMatrixStatus() << "\n";

    
   std::cout << "Minuit status: " << minuit->Status() << std::endl;
  std::cout << "Covariance matrix status: " 
          << minuit->CovMatrixStatus() << std::endl;

    
    TMatrixDSym* Postmatrix = new TMatrixDSym(NparsMinuitFull);

    for(int i = 0; i < NparsMinuitFull; ++i)
    {
      for(int j = 0; j < NparsMinuitFull; ++j)
      {
        (*Postmatrix)(i,j) = 0;
        (*Postmatrix)(i,j) = minuit->CovMatrix(i,j);
      }
    }
    // Open a ROOT file
    TFile* outFile = new TFile("CovMatrix.root", "RECREATE");
    //Postmatrix->GetName("CovarianceMatrix");
    Postmatrix->Write();  // Write the matrix to the file
    outFile->Close();

    delete Postmatrix;
    delete outFile;
}

