void drawCorrelationMatrix_hessian()
{
    // --- Open file and get covariance matrix ---
    TFile* f = TFile::Open("/scratch/abipeake/MaCh3_DUNE_Nov2025/MaCh3_DUNE/CovMatrix_templateparamssubset.root", "READ");
    if (!f || f->IsZombie()) {
        Error("DrawCovAndCorr", "Cannot open CovMatrix.root");
        return;
    }

    TMatrixTSym<double>* cov = nullptr;
    f->GetObject("TMatrixTSym<double>", cov);
    if (!cov) {
        Error("DrawCovAndCorr", "Covariance matrix not found");
        return;
    }

    int N = cov->GetNrows();

    // --- Create histograms ---
    TH2D* hcov = new TH2D(
        "hcov",
        "Covariance Matrix;Parameter ; Parameter ",
        N, 0, N,
        N, 0, N
    );

    TH2D* hcorr = new TH2D(
        "hcorr",
        "Correlation Matrix;Parameter ; Parameter ",
        N, 0, N,
        N, 0, N
    );

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            double cij = (*cov)(i,j);
            hcov->SetBinContent(i+1, j+1, cij);

            double sii = (*cov)(i,i);
            double sjj = (*cov)(j,j);
            hcorr->SetBinContent(i+1, j+1, cij / std::sqrt(sii * sjj));
        }
    }

    // --- Red–white–blue palette ---
    const Int_t NRGBs = 3;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = {0.00, 0.50, 1.00};
    Double_t red[NRGBs]   = {0.00, 1.00, 1.00};
    Double_t green[NRGBs] = {0.00, 1.00, 0.00};
    Double_t blue[NRGBs]  = {1.00, 1.00, 0.00};

    TColor::CreateGradientColorTable(
        NRGBs, stops, red, green, blue, NCont
    );
    gStyle->SetNumberContours(NCont);
    gStyle->SetOptStat(0);

    // --- Canvas ---
    TCanvas* c = new TCanvas("c", "Covariance and Correlation", 900, 800);

    // --- Multi-page PDF ---
    TString pdfName = "CovarianceAndCorrelation_templatesubset.pdf";
    c->Print(pdfName + "[");

    // Covariance
    double max = hcov->GetMaximum();
    hcov->SetMinimum(-max);
    hcov->Draw("COLZ");
    c->Print(pdfName);

    // Correlation
    hcorr->SetMinimum(-1.0);
    hcorr->SetMaximum( 1.0);
    hcorr->Draw("COLZ");
    c->Print(pdfName);

    c->Print(pdfName + "]");

    // --- Cleanup ---
    delete c;
    f->Close();
}
