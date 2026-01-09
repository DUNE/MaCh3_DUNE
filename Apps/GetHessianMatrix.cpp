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
        systematics->SetParameters(pars);

        systematics->AcceptStep();
        

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
        }
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
  
  std::cout << "Number of PDFs: " << DUNEPdfs.size() << "\n";
  if(!DUNEPdfs.empty()) {
      std::cout << "First PDF pointer: " << DUNEPdfs[0] << "\n";
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
  minuit->SetMaxIterations(10000);
  double tolerance = FitManager->raw()["General"]["Minuit2"]["ToleranceLevel"].as<double>();
  minuit->SetTolerance(tolerance);
  MACH3LOG_INFO("Preparing Minuit");
  int ParCounter = 0;


      for(int i = 0; i < xsec->GetNumParams(); ++i, ++ParCounter)
      {
        //KS: Index, name, prior, step scale [different to MCMC],
        minuit->SetVariable(ParCounter, (xsec->GetParName(i)), xsec->GetParInit(i), xsec->GetDiagonalError(i)/10);
        minuit->SetVariableValue(ParCounter, xsec->GetParInit(i));
        //KS: lower bound, upper bound, if Mirroring enabled then ignore
        minuit->SetVariableLimits(ParCounter, xsec->GetLowerBound(i), xsec->GetUpperBound(i));
        if(xsec->IsParameterFixed(i))
        {
          minuit->FixVariable(ParCounter);
        }
      }
    
    MACH3LOG_INFO("Starting MIGRAD");
    minuit->Minimize(); //nuwro

    MACH3LOG_INFO("Starting HESSE");
    minuit->Hesse();
    

    
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

