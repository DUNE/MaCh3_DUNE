#include <iostream>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>
#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDF/GenericBinningTools.h"
#include "mcmc/mcmc.h"

std::string getNameNoExt(std::string name, std::string ext)  
{                                                            
  std::size_t pos ;                                          
  pos = name.find(ext);                                      
  name = name.substr(0,pos);                                 
  return name ;                                              
}                                                             

void saveCanvas(TCanvas* canvas, std::string name, std::string legend)                                                                  
{                                                            
  name = getNameNoExt(name, ".root") ;                       
  name = name + legend + ".root" ;                     
  canvas -> SaveAs(name.c_str()) ;                           

  name = getNameNoExt(name, ".root") ;                       
  name = name + ".png" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           

  name = getNameNoExt(name, ".png") ;                        
  name = name + ".pdf" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           

  name = getNameNoExt(name, ".pdf") ;                        
  name = name + ".eps" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           
} 

int main(int argc, char * argv[]) {

  if(argc == 1){
    std::cout << "Usage: bin/mini_MCMC config.cfg" << std::endl;
    return 1;
  }

  
    bool do_1d_llhscan = true;
    bool do_2d_llhscan = false;

    auto fitMan = std::unique_ptr<manager>(new manager(argv[1]));

    covarianceXsec *xsec = nullptr;
    covarianceOsc *osc = nullptr;

    std::vector<samplePDFFDBase *> DUNEPdfs;
    MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

    gStyle->SetOptStat(false);

    std::vector<double> xsecpar = xsec->getNominalArray();  
    xsec->setParameters();

    std::vector<std::string> sample_names;
    for (auto Sample : DUNEPdfs) {
      Sample->reweight();
      std::string name = Sample->GetSampleName();
      sample_names.push_back(name);

      if (Sample -> GetNDim() == 1) {
        TH1D *Asimov_1D = (TH1D*)Sample->get1DHist()->Clone((name+"_asimov").c_str());
        Sample -> addData(Asimov_1D); 
      }
        
      if (Sample -> GetNDim() == 2) {
        TH2D *Asimov_2D = (TH2D*)Sample->get2DHist()->Clone((name+"_asimov").c_str());
        Sample -> addData(Asimov_2D); 
      }
    }

    std::string OutFileName = fitMan->raw()["General"]["OutputFile"].as<std::string>();
    auto OutFile = std::unique_ptr<TFile>(new TFile(OutFileName.c_str(), "RECREATE"));

  // // ---------------- 2D LLH Scan Setup ---------------- //
  // unsigned int parX = 26; // index of first parameter to scan
  // unsigned int parY = 29; // index of second parameter to scan
  // unsigned int nPoints = fitMan->raw()["LLHScans"]["ParPoints"].as<int>();
  // double nSigma = 3.0;

  // double xMin = xsec->getNominal(parX) - nSigma * xsec->GetError(parX);
  // double xMax = xsec->getNominal(parX) + nSigma * xsec->GetError(parX);
  // double yMin = xsec->getNominal(parY) - nSigma * xsec->GetError(parY);
  // double yMax = xsec->getNominal(parY) + nSigma * xsec->GetError(parY);

  // TH2D* h2DScan = new TH2D("h2DScan", "2D LLH Scan",
  //                          nPoints, xMin, xMax,
  //                          nPoints, yMin, yMax);

  // double dx = (xMax - xMin)/(nPoints-1);
  // double dy = (yMax - yMin)/(nPoints-1);

  // std::vector<double> xsecValues = xsec->getNominalArray();

  // for(unsigned int ix = 0; ix < nPoints; ++ix){
  //   xsecValues[parX] = xMin + ix*dx;

  //   for(unsigned int iy = 0; iy < nPoints; ++iy){
  //     xsecValues[parY] = yMin + iy*dy;

  //     // Set parameters
  //     xsec->setParameters(xsecValues);

  //     // Reweight and calculate total LLH
  //     double totalLLH = 0;
  //     for(auto Sample : DUNEPdfs){
  //         Sample->reweight();
  //         totalLLH += 2 * Sample->GetLikelihood();
  //     }

  //     totalLLH += 2 * xsec->GetLikelihood(); // Systematic penalty

  //     h2DScan->SetBinContent(ix+1, iy+1, totalLLH);
  //   }
  // }

  // OutFile->cd();
  // h2DScan->Write();

  // std::cout << "Finished 2D LLH Scan!" << std::endl;


  unsigned int nPoints = fitMan->raw()["LLHScans"]["ParPoints"].as<int>();
  double nSigma = 3.0;

  std::vector<double> nominalVals = xsec->getNominalArray();


    

  // // ---------------- 2D LLH Scan Setup ---------------- //
  // unsigned int parX = 26; // index of first parameter to scan
  // unsigned int parY = 29; // index of second parameter to scan
  // unsigned int nPoints = fitMan->raw()["LLHScans"]["ParPoints"].as<int>();
  // double nSigma = 3.0;

  // double xMin = xsec->getNominal(parX) - nSigma * xsec->GetError(parX);
  // double xMax = xsec->getNominal(parX) + nSigma * xsec->GetError(parX);
  // double yMin = xsec->getNominal(parY) - nSigma * xsec->GetError(parY);
  // double yMax = xsec->getNominal(parY) + nSigma * xsec->GetError(parY);

  // TH2D* h2DScan = new TH2D("h2DScan", "2D LLH Scan",
  //                          nPoints, xMin, xMax,
  //                          nPoints, yMin, yMax);

  // double dx = (xMax - xMin)/(nPoints-1);
  // double dy = (yMax - yMin)/(nPoints-1);

  // std::vector<double> xsecValues = xsec->getNominalArray();

  // for(unsigned int ix = 0; ix < nPoints; ++ix){
  //   xsecValues[parX] = xMin + ix*dx;

  //   for(unsigned int iy = 0; iy < nPoints; ++iy){
  //     xsecValues[parY] = yMin + iy*dy;

  //     // Set parameters
  //     xsec->setParameters(xsecValues);

  //     // Reweight and calculate total LLH
  //     double totalLLH = 0;
  //     for(auto Sample : DUNEPdfs){
  //         Sample->reweight();
  //         totalLLH += 2 * Sample->GetLikelihood();
  //     }

  //     totalLLH += 2 * xsec->GetLikelihood(); // Systematic penalty

  //     h2DScan->SetBinContent(ix+1, iy+1, totalLLH);
  //   }
  // }

  // OutFile->cd();
  // h2DScan->Write();

  // std::cout << "Finished 2D LLH Scan!" << std::endl;



  // Define which parameters to scan
  std::vector<unsigned int> paramList = {12, 18, 26, 29, 441, 200, 443,444}; // For 1D scans
  std::vector<std::pair<unsigned int, unsigned int>> paramPairs = {
      {26, 29},
      {12, 15},
      {18, 20},
      {441,442},
      {443,444},
      {441,443}
  }; // For 2D scans


  // ---------------- 1D LLH Scans ---------------- //
  if (do_1d_llhscan) {
      for (auto par : paramList) {

          double parNom = xsec->getNominal(par);
          double parErr = xsec->GetError(par);
          double pMin = parNom - nSigma * parErr;
          double pMax = parNom + nSigma * parErr;
          double dp = (pMax - pMin) / (nPoints - 1);

          std::string histName = Form("h1DScan_par%d", par);
          std::string histTitle = Form("1D LLH Scan for Parameter %d", par);
          TH1D* h1DScan = new TH1D(histName.c_str(), histTitle.c_str(),
                                  nPoints, pMin, pMax);

          std::vector<double> xsecValues = nominalVals;

          for (unsigned int i = 0; i < nPoints; ++i) {
              xsecValues[par] = pMin + i * dp;

              xsec->setParameters(xsecValues);

              double totalLLH = 0;
              for (auto Sample : DUNEPdfs) {
                  Sample->reweight();
                  totalLLH += 2 * Sample->GetLikelihood();
              }
              totalLLH += 2 * xsec->GetLikelihood();

              h1DScan->SetBinContent(i + 1, totalLLH);
          }

          OutFile->cd();
          h1DScan->Write();
          std::cout << "Finished 1D LLH scan for parameter " << par << std::endl;

          delete h1DScan;
      }
  }


  // ---------------- 2D LLH Scans ---------------- //
  if (do_2d_llhscan) {
      for (auto [parX, parY] : paramPairs) {

          double xMin =0; //xsec->getNominal(parX) - nSigma * xsec->GetError(parX);
          double xMax = 3;//xsec->getNominal(parX) + nSigma * xsec->GetError(parX);
          double yMin =0.; //xsec->getNominal(parY) - nSigma * xsec->GetError(parY);
          double yMax = 3;//xsec->getNominal(parY) + nSigma * xsec->GetError(parY);

          double dx = (xMax - xMin) / (nPoints - 1);
          double dy = (yMax - yMin) / (nPoints - 1);

          std::string histName = Form("h2DScan_par%d_par%d", parX, parY);
          std::string histTitle = Form("2D LLH Scan for Par %d vs Par %d", parX, parY);

          TH2D* h2DScan = new TH2D(histName.c_str(), histTitle.c_str(),
                                    nPoints, xMin, xMax,
                                    nPoints, yMin, yMax);

          std::vector<double> xsecValues = nominalVals;

          for (unsigned int ix = 0; ix < nPoints; ++ix) {
              xsecValues[parX] = xMin + ix * dx;

              for (unsigned int iy = 0; iy < nPoints; ++iy) {
                  xsecValues[parY] = yMin + iy * dy;

                  xsec->setParameters(xsecValues);

                  double totalLLH = 0;
                  for (auto Sample : DUNEPdfs) {
                      Sample->reweight();
                      totalLLH += 2 * Sample->GetLikelihood();
                  }
                  totalLLH += 2 * xsec->GetLikelihood();

                  h2DScan->SetBinContent(ix + 1, iy + 1, totalLLH);
              }
          }

          OutFile->cd();
          h2DScan->Write();
          std::cout << "Finished 2D LLH scan for parameters "
                    << parX << " and " << parY << std::endl;

          delete h2DScan;
      }
  }

  std::cout << "All LLH scans completed successfully!" << std::endl;

  return 0;


//   // Define which parameters to scan
//   std::vector<unsigned int> paramList = {12, 18, 26, 29, 441, 200, 443,444}; // For 1D scans
//   std::vector<std::pair<unsigned int, unsigned int>> paramPairs = {
//       {26, 29},
//       {12, 15},
//       {18, 20},
//       {441,442},
//       {443,444},
//       {441,443}
//   }; // For 2D scans


//   // ---------------- 1D LLH Scans ---------------- //
//   if (do_1d_llhscan) {
//       for (auto par : paramList) {

//           double parNom = xsec->getNominal(par);
//           double parErr = xsec->GetError(par);
//           double pMin = parNom - nSigma * parErr;
//           double pMax = parNom + nSigma * parErr;
//           double dp = (pMax - pMin) / (nPoints - 1);

//           std::string histName = Form("h1DScan_par%d", par);
//           std::string histTitle = Form("1D LLH Scan for Parameter %d", par);
//           TH1D* h1DScan = new TH1D(histName.c_str(), histTitle.c_str(),
//                                   nPoints, pMin, pMax);

//           std::vector<double> xsecValues = nominalVals;

//           for (unsigned int i = 0; i < nPoints; ++i) {
//               xsecValues[par] = pMin + i * dp;

//               xsec->setParameters(xsecValues);

//               double totalLLH = 0;
//               for (auto Sample : DUNEPdfs) {
//                   Sample->reweight();
//                   totalLLH += 2 * Sample->GetLikelihood();
//               }
//               totalLLH += 2 * xsec->GetLikelihood();

//               h1DScan->SetBinContent(i + 1, totalLLH);
//           }

//           OutFile->cd();
//           h1DScan->Write();
//           std::cout << "Finished 1D LLH scan for parameter " << par << std::endl;

//           delete h1DScan;
//       }
//   }


//   // ---------------- 2D LLH Scans ---------------- //
//   if (do_2d_llhscan) {
//       for (auto [parX, parY] : paramPairs) {

//           double xMin =0; //xsec->getNominal(parX) - nSigma * xsec->GetError(parX);
//           double xMax = 3;//xsec->getNominal(parX) + nSigma * xsec->GetError(parX);
//           double yMin =0.; //xsec->getNominal(parY) - nSigma * xsec->GetError(parY);
//           double yMax = 3;//xsec->getNominal(parY) + nSigma * xsec->GetError(parY);

//           double dx = (xMax - xMin) / (nPoints - 1);
//           double dy = (yMax - yMin) / (nPoints - 1);

//           std::string histName = Form("h2DScan_par%d_par%d", parX, parY);
//           std::string histTitle = Form("2D LLH Scan for Par %d vs Par %d", parX, parY);

//           TH2D* h2DScan = new TH2D(histName.c_str(), histTitle.c_str(),
//                                     nPoints, xMin, xMax,
//                                     nPoints, yMin, yMax);

//           std::vector<double> xsecValues = nominalVals;

//           for (unsigned int ix = 0; ix < nPoints; ++ix) {
//               xsecValues[parX] = xMin + ix * dx;

//               for (unsigned int iy = 0; iy < nPoints; ++iy) {
//                   xsecValues[parY] = yMin + iy * dy;

//                   xsec->setParameters(xsecValues);

//                   double totalLLH = 0;
//                   for (auto Sample : DUNEPdfs) {
//                       Sample->reweight();
//                       totalLLH += 2 * Sample->GetLikelihood();
//                   }
//                   totalLLH += 2 * xsec->GetLikelihood();

//                   h2DScan->SetBinContent(ix + 1, iy + 1, totalLLH);
//               }
//           }

//           OutFile->cd();
//           h2DScan->Write();
//           std::cout << "Finished 2D LLH scan for parameters "
//                     << parX << " and " << parY << std::endl;

//           delete h2DScan;
//       }
//   }

// std::cout << "All LLH scans completed successfully!" << std::endl;

//   return 0;
}
