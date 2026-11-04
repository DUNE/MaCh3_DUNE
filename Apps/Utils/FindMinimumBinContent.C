
#include <iostream>
#include <algorithm>
#include "TFile.h"
#include "TH1D.h"

void FindMinimumBinContent()
{
    // Open ROOT file
    TFile *file = TFile::Open("DuneEventRates_Y.root", "READ");

    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open DuneEventRates.root" << std::endl;
        return;
    }

    // List of histogram names
    const char* histNames[] = {
        "hFC_numuselec",
        "hFC_nueselec",
        "hPC_numuselec",
        "hPC_nueselec"
    };

    // Loop over histograms
    for (const char* name : histNames) {

        TH1D *hist = (TH1D*)file->Get(name);

        if (!hist) {
            std::cerr << "Histogram " << name
                      << " not found!" << std::endl;
            continue;
        }

        // Find minimum bin content (excluding underflow and overflow)
        double minContent = hist->GetBinContent(1);
        int minBin = 1;

        for (int bin = 2; bin <= hist->GetNbinsX(); bin++) {

            double content = hist->GetBinContent(bin);

            if (content < minContent) {
                minContent = content;
                minBin = bin;
            }
        }

        std::cout << "Histogram: " << name << std::endl;
        std::cout << "  Minimum bin content = " << minContent << std::endl;
        std::cout << "  Bin number          = " << minBin << std::endl;
        std::cout << "  Bin centre          = "
                  << hist->GetBinCenter(minBin) << std::endl;
    }

    file->Close();
    delete file;
}
