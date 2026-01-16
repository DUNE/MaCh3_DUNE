// void ImproveProposalMatrix() {
//   TFile f("/scratch/abipeake/MaCh3_DUNE_Nov2025/MaCh3_DUNE/CovMatrix_templateparameters.root", "READ");
//   auto* C = (TMatrixTSym<double>*)f.Get("TMatrixTSym<double>");
//   if (!C) {
//     std::cerr << "Matrix not found!" << std::endl;
//     return;
//   }

//   const int N = C->GetNrows();
//   TMatrixDSym Csym(N);
//   Csym = *C;

//   // Eigen decomposition
//   TMatrixDSymEigen eig(Csym);
//   TVectorD evals = eig.GetEigenValues();
//   TMatrixD evecs = eig.GetEigenVectors();

//   // Floor eigenvalues
//   const double lambda_floor = 1.0;
//   for (int i = 0; i < N; ++i) {
//     if (evals[i] < lambda_floor)
//       evals[i] = lambda_floor;
//   }

//   // Reconstruct matrix: V * D * V^T
//   TMatrixD D(N, N);
//   for (int i = 0; i < N; ++i)
//     D(i,i) = evals[i];

//   TMatrixD Cnew = evecs * D * TMatrixD(TMatrixD::kTransposed, evecs);

//   // Save
//   TFile fout("/scratch/abipeake/MaCh3_DUNE_Nov2025/MaCh3_DUNE/CovMatrix_templateparameters_floored.root", "RECREATE");
//   TMatrixTSym<double> Csave(N);
//   for (int i = 0; i < N; ++i)
//     for (int j = 0; j <= i; ++j)
//       Csave(i,j) = Cnew(i,j);

//   Csave.Write("TMatrixTSym");
//   fout.Close();

//   std::cout << "Saved floored proposal matrix." << std::endl;
// }

void ImproveProposalMatrix() {

  const char* inFile  = "/scratch/abipeake/MaCh3_DUNE_Nov2025/MaCh3_DUNE/CovMatrix_templateparameters.root";
  const char* outFile = "/scratch/abipeake/MaCh3_DUNE_Nov2025/MaCh3_DUNE/CovMatrix_templateparameters_floored.root";
  const char* matName = "TMatrixTSym<double>";

  const char* meanName = "xsec_cov_means";

  TFile f(inFile, "READ");
  auto* C = dynamic_cast<TMatrixTSym<double>*>(f.Get(matName));
  if (!C) {
    std::cerr << "Input matrix not found!" << std::endl;
    return;
  }

  const int N = C->GetNrows();

  // Convert to DSym
  TMatrixDSym Csym(N);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j <= i; ++j)
      Csym(i,j) = (*C)(i,j);

  // Eigen decomposition
  TMatrixDSymEigen eig(Csym);
  TVectorD evals = eig.GetEigenValues();
  TMatrixD evecs = eig.GetEigenVectors();

  // Relative eigenvalue floor
  double maxEval = evals.Max();
  double lambda_floor = 1e-6 * maxEval;

  for (int i = 0; i < N; ++i)
    if (evals[i] < lambda_floor)
      evals[i] = lambda_floor;

  // Reconstruct covariance
  TMatrixD D(N,N);
  for (int i = 0; i < N; ++i)
    D(i,i) = evals[i];

  TMatrixD Cnew = evecs * D * TMatrixD(TMatrixD::kTransposed, evecs);

  // Save outputs
  TFile fout(outFile, "RECREATE");

  TMatrixTSym<double> Csave(N);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j <= i; ++j)
      Csave(i,j) = Cnew(i,j);

  Csave.Write(matName);

  // Means vector (CRITICAL)
  TVectorD means(N);
  for (int i = 0; i < N; ++i)
    means(i) = 1.0;

  means.Write(meanName);

  fout.Close();

  std::cout << "Saved improved proposal + means vector" << std::endl;
}
