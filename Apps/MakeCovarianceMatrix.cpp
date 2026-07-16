// #include "Fitters/MaCh3Factory.h"
// #include "Samples/MaCh3DUNEFactory.h"

// #include "TFile.h"

// #include "yaml-cpp/yaml.h"

// #include "Eigen/Dense"

// #include <cmath>
// #include <iostream>
// #include <vector>

// // #define DEBUG_MAKE_COV_MX

// int main(int argc, char *argv[]) {

//   if (argc <= 5) {
//     std::cout << "[ERROR]: Runlike " << argv[0]
//               << " <Mach3_Config.yml> <Parameter.yml> <NThrows> <output.root>"
//               << std::endl;
//   }

//   auto conf_yml = M3OpenConfig(argv[1]);
//   conf_yml["General"]["Systematics"] =
//       YAML::Load(std::string("XsecCovFile:\n\t[") + argv[2] +
//                  "]\nXsecCovName: xsec_cov\nXsecStepScale: 1.0");
//   auto fitMan = std::make_unique<Manager>(conf_yml);

//   auto nThrows = std::stol(argv[3]);
//   auto OutputFileName = std::string(argv[4]);

//   ParameterHandlerGeneric *xsec = nullptr;
//   std::vector<SampleHandlerFD *> DUNEPdfs;
//   MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

  

//   if (DUNEPdfs.size() != 1) {
//     throw std::runtime_error("only want a single samplehandler");
//   }

//   auto &pdf = *DUNEPdfs.front();
  

//   auto get_pred = [&]() -> Eigen::VectorXd {
//     pdf.Reweight();
//     auto mcvals = pdf.GetMCArray();

// #ifdef DEBUG_MAKE_COV_MX
//     std::cout << "[ ";
//     for (auto &v : mcvals) {
//       std::cout << v << ", ";
//     }
//     std::cout << "]" << std::endl;
// #endif

//     return Eigen::Map<Eigen::VectorXd>(mcvals.data(), mcvals.size());
//   };
//   Eigen::VectorXd cvpred = get_pred();

//     std::vector<TH1D> binHists;
// binHists.reserve(cvpred.size());
// for (int i = 0; i < cvpred.size(); ++i) {
//     binHists.emplace_back(
//     Form("bin_%d", i),
//     Form("Bin %d throw distribution", i),
//     100, 0.0, 1e6
// );
// }

//   // Set all parameters to prior values before getting nominal
//   for (int i = 0; i < xsec->GetNumParams(); ++i) {
//     xsec->SetPar(i, xsec->GetParInit(i));
//   }
  
//   Eigen::VectorXd mean = Eigen::VectorXd::Zero(cvpred.size());
//   MACH3LOG_INFO("cvpred.size() = {}", cvpred.size());

//   Eigen::MatrixXd CovMatrix = Eigen::MatrixXd::Zero(cvpred.size(), cvpred.size());

//   double rms_accum =0;
//   for (int iThrow = 0; iThrow < nThrows; ++iThrow) {

//     xsec->ThrowParameters();
//     if (iThrow == 0) {
//     std::cout << "\nThrown parameters:\n";

//     for (int i = 0; i < xsec->GetNumParams(); ++i) {
//       std::cout<< i<< " "<< xsec->GetParName(i)<< " = "<< xsec->GetCorrThrows(i)<< std::endl;}
//     }
//     Eigen::VectorXd throw_pred = get_pred();
//     mean += throw_pred;
//     Eigen::VectorXd diff = (throw_pred - cvpred);
//     CovMatrix += (diff * diff.transpose());

//     double rmsShift = diff.norm();
//     rms_accum += rmsShift;
//     if ((iThrow & 1023) == 0) {
//         MACH3LOG_INFO("Throw {}/{}", iThrow, nThrows);
//     }
//     for (int i = 0; i < cvpred.size(); ++i) {
//       binHists[i].Fill(throw_pred(i));
//     }   
//   }
//   std::cout<< "Average throw norm = "<< rms_accum / nThrows<< std::endl;

//   mean      *= 1.0 / double(nThrows);
//   CovMatrix *= 1.0 / double(nThrows - 1);

//   Eigen::VectorXd meanDiff = mean - cvpred;

//   std::cout<< "||mean-cv|| = "<< meanDiff.norm()<< std::endl;

//   std::cout << "\nBin  CV  Sigma  FracErr\n";

//   for (int i = 0; i < cvpred.size(); ++i) {
//       double sigma = std::sqrt(CovMatrix(i,i));

//       std::cout
//         << i
//         << "  " << cvpred(i)
//         << "  " << sigma
//         << "  " << sigma/cvpred(i)
//         << "\n";
//   }
//   // Reset to nominal before getting W2
//   for (int i = 0; i < xsec->GetNumParams(); ++i) {
//       xsec->SetPar(i, xsec->GetParInit(i));
//   }
//   //pdf.Reweight();
//   //auto w2vals = pdf.GetW2Array();  // std::vector<double>, sum(w_i^2) per bin
//   //auto events_permcbin = pdf.GetMCArray();
//   // Add  stability MC stat uncertainty diagonal
//   for (int i = 0; i < cvpred.size(); ++i) {
//        //std::cout<< "w2vals[i] = " << w2vals[i] << std::endl;
//        //std::cout<< "sqrt w2vals[i] = " << sqrt(w2vals[i]) << std::endl;
//        //std::cout<< "events_permcbin[i] = " << events_permcbin[i] << std::endl;
//        //std::cout << "w2vals[i]/events_permcbin[i] = "  << w2vals[i]/events_permcbin[i] << std::endl;
//        CovMatrix(i, i) +=  1e-15; //+ w2vals[i]/events_permcbin[i];
//        //std::cout << "xsec throw variance / w2 = " 
//          // << CovMatrix(i,i) / w2vals[i] << std::endl;
//    } 

//    double avgFrac = 0.0;

//   for (int i = 0; i < cvpred.size(); ++i) {
//       avgFrac += std::sqrt(CovMatrix(i,i)) / cvpred(i);
//   }

//   avgFrac /= cvpred.size();

//   std::cout
//     << "Average fractional uncertainty = "
//     << avgFrac
//     << std::endl;

//   auto OutputFile =
//       std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
//   OutputFile->cd();

//   TMatrixD cvpred_root(cvpred.size(), 1, cvpred.data());
//   TMatrixD mean_root(mean.size(), 1, mean.data());
//   TMatrixD CovMatrix_root(CovMatrix.rows(), CovMatrix.cols(), CovMatrix.data());
 
//   Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(CovMatrix);

// std::cout
//   << "Min eigenvalue = "
//   << eig.eigenvalues().minCoeff()
//   << std::endl;

// std::cout
//   << "Max eigenvalue = "
//   << eig.eigenvalues().maxCoeff()
//   << std::endl;

//   std::cout
//   << "Diag min = "
//   << CovMatrix.diagonal().minCoeff()
//   << std::endl;

// std::cout
//   << "Diag max = "
//   << CovMatrix.diagonal().maxCoeff()
//   << std::endl;

// cvpred_root.Write("NDCV");
// mean_root.Write("NDThrowsMean");
// CovMatrix_root.Write("NDCovMatrix");

// for (auto &h : binHists) {
//     h.Write();
// }

// OutputFile->Write();
// OutputFile->Close();
// }


// //   cvpred_root.Write("NDCV");
// //   mean_root.Write("NDThrowsMean");
// //   CovMatrix_root.Write("NDCovMatrix");

// //   // write PDF of covmatrix (correlation matrix)

// //   OutputFile->Write();
// //   OutputFile->Close();
// // mean      *= 1.0 / double(nThrows);
// // CovMatrix *= 1.0 / double(nThrows - 1);

// // // Reset to nominal before getting W2
// // for (int i = 0; i < xsec->GetNumParams(); ++i) {
// //     xsec->SetPar(i, xsec->GetParInit(i));
// // }
// // pdf.Reweight();
// // auto w2vals = pdf.GetW2Array();

// // // --- Identify non-empty bins ---
// // std::vector<int> activeBins;
// // for (int i = 0; i < cvpred.size(); ++i) {
// //     if (cvpred(i) > 0 || w2vals[i] > 0) {
// //         activeBins.push_back(i);
// //     } else {
// //         MACH3LOG_INFO("Skipping empty bin {}: MC={}, W2={}", i, cvpred(i), w2vals[i]);
// //     }
// // }
// // MACH3LOG_INFO("Active bins: {}/{}", activeBins.size(), cvpred.size());

// // // --- Build reduced covariance matrix over active bins only ---
// // int nActive = activeBins.size();
// // Eigen::MatrixXd CovReduced = Eigen::MatrixXd::Zero(nActive, nActive);
// // for (int i = 0; i < nActive; ++i) {
// //     for (int j = 0; j < nActive; ++j) {
// //         CovReduced(i, j) = CovMatrix(activeBins[i], activeBins[j]);
// //     }
// // }

// // // --- Add W2 diagonal to reduced matrix ---
// // for (int i = 0; i < nActive; ++i) {
// //     MACH3LOG_DEBUG("w2vals[{}] = {}", activeBins[i], w2vals[activeBins[i]]);
// //     CovReduced(i, i) += w2vals[activeBins[i]];
// // }

// // // --- Sanity checks ---
// // Eigen::FullPivLU<Eigen::MatrixXd> lu(CovReduced);
// // MACH3LOG_INFO("Reduced matrix rank: {}/{}", lu.rank(), nActive);

// // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(CovReduced);
// // MACH3LOG_INFO("Min eigenvalue: {}", solver.eigenvalues().minCoeff());
// // MACH3LOG_INFO("Max eigenvalue: {}", solver.eigenvalues().maxCoeff());
// // if (solver.eigenvalues().minCoeff() <= 0) {
// //     MACH3LOG_WARN("Matrix is not positive definite!");
// // }

// // // --- Also store reduced cvpred and mean ---
// // Eigen::VectorXd cvpredReduced(nActive);
// // Eigen::VectorXd meanReduced(nActive);
// // for (int i = 0; i < nActive; ++i) {
// //     cvpredReduced(i) = cvpred(activeBins[i]);
// //     meanReduced(i)   = mean(activeBins[i]);
// // }

// // // --- Write to ROOT file ---
// // auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
// // OutputFile->cd();

// // // Store active bin indices for later reference
// // TVectorD activeBins_root(nActive);
// // for (int i = 0; i < nActive; ++i) activeBins_root(i) = activeBins[i];
// // activeBins_root.Write("ActiveBins");

// // TMatrixD cvpred_root(nActive, 1, cvpredReduced.data());
// // TMatrixD mean_root(nActive, 1, meanReduced.data());

// // // Fix Eigen column-major -> ROOT row-major
// // Eigen::MatrixXd CovReduced_rm = CovReduced;
// // TMatrixD CovMatrix_root(nActive, nActive, CovReduced_rm.data());


// // // Check correlation matrix
// // for (int i = 0; i < nActive; ++i) {
// //     for (int j = 0; j < nActive; ++j) {
// //         double rho = CovReduced(i,j) / std::sqrt(CovReduced(i,i) * CovReduced(j,j));
// //         if (std::abs(rho) > 1.0 + 1e-6) {
// //             MACH3LOG_WARN("Unphysical correlation bin ({},{}): rho={}", i, j, rho);
// //         }
// //     }
// // }



// /////////////////////
// //adding stats uncert
//   // mean     *= 1.0 / double(nThrows);
//   // CovMatrix *= 1.0 / double(nThrows - 1);

//   // // Add Poisson statistical term 
//   // for (int i = 0; i < cvpred.size(); ++i) {
//   //     CovMatrix(i, i) += cvpred(i); 
//   // }

// //   for (int iThrow = 0; iThrow < nThrows; ++iThrow) {

// //     xsec->ThrowParameters();

// //     Eigen::VectorXd throw_pred = get_pred();

// //     mean += throw_pred;

// //     Eigen::VectorXd diff = (throw_pred - cvpred);

// // #ifdef DEBUG_MAKE_COV_MX
// //     std::cout << "diff: [ ";
// //     for (int i = 0; i < diff.rows(); ++i) {
// //       std::cout << diff(i, 0) << ", ";
// //     }
// //     std::cout << "]" << std::endl;
// // #endif

// //     CovMatrix += (diff * diff.transpose());
// //     for (int i = 0; i < cvpred.size(); ++i) {
// //     CovMatrix(i, i) += cvpred(i);
// // }

// //     if (iThrow % 100 == 0) {
// //       MACH3LOG_INFO("Throw {}/{}", iThrow, nThrows);
// //     }
// //   }

// //   mean *= 1.0 / double(nThrows);
// //   CovMatrix *= 1.0 / double(nThrows - 1);


#include "Fitters/MaCh3Factory.h"
#include "Samples/MaCh3DUNEFactory.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1D.h"

#include "yaml-cpp/yaml.h"

#include "Eigen/Dense"

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

// #define DEBUG_MAKE_COV_MX

int main(int argc, char *argv[]) {

  if (argc <= 5) {
    std::cout << "[ERROR]: Runlike " << argv[0]
              << " <Mach3_Config.yml> <Parameter.yml> <NThrows> <output.root>"
              << std::endl;
  }

  auto conf_yml = M3OpenConfig(argv[1]);
  conf_yml["General"]["Systematics"] =
      YAML::Load(std::string("XsecCovFile:\n\t[") + argv[2] +
                 "]\nXsecCovName: xsec_cov\nXsecStepScale: 1.0");
  auto fitMan = std::make_unique<Manager>(conf_yml);

  auto nThrows = std::stol(argv[3]);
  auto OutputFileName = std::string(argv[4]);

  ParameterHandlerGeneric *xsec = nullptr;
  std::vector<SampleHandlerFD *> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

  if (DUNEPdfs.size() != 1) {
    throw std::runtime_error("only want a single samplehandler");
  }

  auto &pdf = *DUNEPdfs.front();

  auto get_pred = [&]() -> Eigen::VectorXd {
    pdf.Reweight();
    auto mcvals = pdf.GetMCArray();

#ifdef DEBUG_MAKE_COV_MX
    std::cout << "[ ";
    for (auto &v : mcvals) {
      std::cout << v << ", ";
    }
    std::cout << "]" << std::endl;
#endif

    return Eigen::Map<Eigen::VectorXd>(mcvals.data(), mcvals.size());
  };
  Eigen::VectorXd cvpred = get_pred();

  std::vector<TH1D> binHists;
  binHists.reserve(cvpred.size());
  for (int i = 0; i < cvpred.size(); ++i) {
    binHists.emplace_back(
        Form("bin_%d", i),
        Form("Bin %d throw distribution", i),
        100, 0.0, 1e6);
  }

  // Set all parameters to prior values before getting nominal
  for (int i = 0; i < xsec->GetNumParams(); ++i) {
    xsec->SetPar(i, xsec->GetParInit(i));
  }

  Eigen::VectorXd mean = Eigen::VectorXd::Zero(cvpred.size());
  MACH3LOG_INFO("cvpred.size() = {}", cvpred.size());

  Eigen::MatrixXd CovMatrix = Eigen::MatrixXd::Zero(cvpred.size(), cvpred.size());

  // --- before the throw loop ---
  int nExampleThrows = 5; // how many spectra you want to save
  std::vector<int> exampleThrowIdx;
  for (int k = 0; k < nExampleThrows; ++k) {
    exampleThrowIdx.push_back((nThrows / nExampleThrows) * k); // evenly spaced
  }

  // Real reco-energy-axis histograms, one slot per saved example throw.
  // Filled via Get1DVarHist() below, so no fixed binning is assumed here.
  std::vector<std::unique_ptr<TH1>> exampleSpectra(nExampleThrows);

  // CV spectrum on the real reco-energy axis (grabbed once, at nominal params)
  std::unique_ptr<TH1> cvHist(
      (TH1 *)pdf.Get1DVarHist(0, pdf.GetXBinVarName(0), {})->Clone("cv_spectrum"));
  cvHist->SetDirectory(nullptr);

  double rms_accum = 0;
  for (int iThrow = 0; iThrow < nThrows; ++iThrow) {

    xsec->ThrowParameters();
    if (iThrow == 0) {
      std::cout << "\nThrown parameters:\n";

      for (int i = 0; i < xsec->GetNumParams(); ++i) {
        std::cout << i << " " << xsec->GetParName(i) << " = "
                  << xsec->GetCorrThrows(i) << std::endl;
      }
    }
    Eigen::VectorXd throw_pred = get_pred();
    mean += throw_pred;
    Eigen::VectorXd diff = (throw_pred - cvpred);
    CovMatrix += (diff * diff.transpose());

    // save this throw's spectrum if selected
    auto it = std::find(exampleThrowIdx.begin(), exampleThrowIdx.end(), iThrow);
    if (it != exampleThrowIdx.end()) {
      int k = std::distance(exampleThrowIdx.begin(), it);
      TH1 *h = pdf.Get1DVarHist(0, pdf.GetXBinVarName(0), {});
      exampleSpectra[k].reset((TH1 *)h->Clone(Form("throw_spectrum_%d", iThrow)));
      exampleSpectra[k]->SetDirectory(nullptr);
    }

    double rmsShift = diff.norm();
    rms_accum += rmsShift;
    if ((iThrow & 1023) == 0) {
      MACH3LOG_INFO("Throw {}/{}", iThrow, nThrows);
    }
    for (int i = 0; i < cvpred.size(); ++i) {
      binHists[i].Fill(throw_pred(i));
    }
  }
  std::cout << "Average throw norm = " << rms_accum / nThrows << std::endl;

  mean *= 1.0 / double(nThrows);
  CovMatrix *= 1.0 / double(nThrows - 1);

  Eigen::VectorXd meanDiff = mean - cvpred;

  std::cout << "||mean-cv|| = " << meanDiff.norm() << std::endl;

  std::cout << "\nBin  CV  Sigma  FracErr\n";

  for (int i = 0; i < cvpred.size(); ++i) {
    double sigma = std::sqrt(CovMatrix(i, i));

    std::cout
        << i
        << "  " << cvpred(i)
        << "  " << sigma
        << "  " << sigma / cvpred(i)
        << "\n";
  }
  // Reset to nominal before getting W2
  for (int i = 0; i < xsec->GetNumParams(); ++i) {
    xsec->SetPar(i, xsec->GetParInit(i));
  }
  // pdf.Reweight();
  // auto w2vals = pdf.GetW2Array();  // std::vector<double>, sum(w_i^2) per bin
  // auto events_permcbin = pdf.GetMCArray();
  // Add stability MC stat uncertainty diagonal
  for (int i = 0; i < cvpred.size(); ++i) {
    CovMatrix(i, i) += 1e-15; // + w2vals[i]/events_permcbin[i];
  }

  double avgFrac = 0.0;

  for (int i = 0; i < cvpred.size(); ++i) {
    avgFrac += std::sqrt(CovMatrix(i, i)) / cvpred(i);
  }

  avgFrac /= cvpred.size();

  std::cout
      << "Average fractional uncertainty = "
      << avgFrac
      << std::endl;

  auto OutputFile =
      std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();

  TMatrixD cvpred_root(cvpred.size(), 1, cvpred.data());
  TMatrixD mean_root(mean.size(), 1, mean.data());
  TMatrixD CovMatrix_root(CovMatrix.rows(), CovMatrix.cols(), CovMatrix.data());

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(CovMatrix);

  std::cout
      << "Min eigenvalue = "
      << eig.eigenvalues().minCoeff()
      << std::endl;

  std::cout
      << "Max eigenvalue = "
      << eig.eigenvalues().maxCoeff()
      << std::endl;

  std::cout
      << "Diag min = "
      << CovMatrix.diagonal().minCoeff()
      << std::endl;

  std::cout
      << "Diag max = "
      << CovMatrix.diagonal().maxCoeff()
      << std::endl;

  TCanvas c("c_spectra", "spectra", 900, 700);

  cvHist->SetLineColor(kBlack);
  cvHist->SetLineWidth(3);
  cvHist->Draw("HIST");

  int colors[] = {kRed, kBlue, kGreen + 2, kMagenta, kOrange + 1};
  TLegend leg(0.7, 0.7, 0.9, 0.9);
  leg.AddEntry(cvHist.get(), "CV", "l");

  for (int k = 0; k < nExampleThrows; ++k) {
    if (!exampleSpectra[k]) continue; // guard in case a slot never got filled
    exampleSpectra[k]->SetLineColor(colors[k % 5]);
    exampleSpectra[k]->SetLineWidth(2);
    exampleSpectra[k]->Draw("HIST SAME");
    leg.AddEntry(exampleSpectra[k].get(), Form("Throw %d", exampleThrowIdx[k]), "l");
  }

  leg.Draw();
  c.Print("varied_spectra.pdf");

  cvpred_root.Write("NDCV");
  mean_root.Write("NDThrowsMean");
  CovMatrix_root.Write("NDCovMatrix");

  cvHist->Write();
  for (int k = 0; k < nExampleThrows; ++k) {
    if (exampleSpectra[k]) exampleSpectra[k]->Write();
  }

  for (auto &h : binHists) {
    h.Write();
  }

  OutputFile->Write();
  OutputFile->Close();
}