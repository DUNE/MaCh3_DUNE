#include <TFile.h>
#include <TKey.h>
#include <TMatrixTSym.h>
#include <TVectorT.h>
#include <iostream>
#include <fstream>

void dump_cov() {

  const char* filename =
    "/scratch/abipeake/DecemberChains/"
    "Saturday100126_templateparams_adaptiveshorter_2D/"
    "Saturday100126_templateparams_adaptiveshorter_2D_adaptive_chain_0_job_0.root";

  TFile f(filename, "READ");
  if (f.IsZombie()) {
    std::cerr << "Error opening file\n";
    return;
  }

  // Get covariance matrix (first key)
  TKey* k0 = (TKey*) f.GetListOfKeys()->At(0);
  auto* cov = (TMatrixTSym<double>*) k0->ReadObj();

  // Dump to text
  std::ofstream out("adaptive_cov_dump.txt");

  int n = cov->GetNrows();
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      out << (*cov)(i,j) << " ";
    }
    out << "\n";
  }

  out.close();

  std::cout << "Wrote adaptive covariance to adaptive_cov_dump.txt\n";
}
