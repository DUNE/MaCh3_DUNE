void scale_covariance() {
    // Open original file
    TFile *oldFile = TFile::Open("/scratch/abipeake/October_chains/Tuesday_041125_enubuas_postadaptivetestnew/combined_drawCorr.root");
    
    // Get the covariance matrix
    auto *Covariance = (TMatrixTSym<double>*)oldFile->Get("Covariance");
    int dim = Covariance->GetNrows();  // matrix dimension

    // Compute scaling factor
    double scale = pow(2.38, 2) / dim;
    std::cout << "Scaling factor = " << scale << std::endl;

    // Create a scaled copy
    TMatrixTSym<double> Cov_scaled(*Covariance);
    Cov_scaled *= scale;

    // Write to new file
    TFile *newFile = new TFile("/scratch/abipeake/October_chains/Tuesday_041125_enubuas_postadaptivetestnew/scaled_covariance.root", "RECREATE");
    Cov_scaled.Write("Covariance_scaled");

    newFile->Close();
    oldFile->Close();
}
