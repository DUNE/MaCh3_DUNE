#include <TH2.h>

// This macro should be fed 2D binned true x reco for some variable binned very finely.
// This macro will try to combine bins in such a way as so "threshold" proportion of
// truth will be in each reco bin.
// Adjust "threshold" as needed

const double threshold = 0.5;

double efficiency(const TH2D *h, const int startBin, const int endBin) {
    double inBin = 0;
    double outBin = 0;
    for (int i = startBin; i <= endBin; i++) {
        for (int j = startBin; j <= endBin; j++) {
            inBin += h->GetBinContent(i, j);
        }
    }
    for (int i = startBin; i <= endBin; i++) {
        for (int j = 1; j < startBin; j++) {
            outBin += h->GetBinContent(i, j);
        }
    }
    for (int i = startBin; i <= endBin; i++) {
        for (int j = endBin+1; j <= h->GetNbinsY(); j++) {
            outBin += h->GetBinContent(i, j);
        }
    }
    return inBin / (inBin + outBin);
}

TH2D* flipHist(TH2D *h) {
    TH2D *h2 = new TH2D("{h.GetName()}_flip", "{h.GetName()}_flip", h->GetNbinsY(),
                        h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax(), h->GetNbinsX(),
                        h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());

    for(int i = 1; i <= h->GetNbinsX(); i++) {
        for(int j = 1; j <= h->GetNbinsY(); j++) {
            h2->SetBinContent(j, i, h->GetBinContent(i, j));
        }
    }
    return h2;
}

void chooseBinningAtm(std::string inFile, std::string histName, bool flip=false) {
    TFile *file = TFile::Open(inFile.c_str());
    TH2D *h = (TH2D *) file->Get(histName.c_str());

    // Need this because issue https://github.com/mach3-software/MaCh3/issues/245
    // Need to switch ordering of x-y to avoid segfaults for some variables
    if (flip) {
        h = flipHist(h);
    }

    TAxis *axis = h->GetXaxis();
    
    for (int i = 1; i <= h->GetNbinsX(); i++) {
        for (int j = i; j <= h->GetNbinsX(); j++) {
            double e = efficiency(h, i, j);
            if (j == h->GetNbinsX() || e >= threshold) {
                std::cout << axis->GetBinLowEdge(i) << " -> " << axis->GetBinLowEdge(j) + axis->GetBinWidth(j)
                        << " : " << 100 * e << "%" << std::endl;
                i = j + 1;
            }
        }
    }
}
