// AdaptiveTailMerge2D_macro.C
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
//"/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/Chain_testing/CERN/Eventrate_q0q310e4.root


TH2D* AdaptiveTailMerge2D_macro(const TH2D* h2, double minCounts) {
    int nx = h2->GetNbinsX();
    int ny = h2->GetNbinsY();

    std::vector<double> xEdges;
    std::vector<double> yEdges;

    xEdges.push_back(h2->GetXaxis()->GetBinLowEdge(1));
    yEdges.push_back(h2->GetYaxis()->GetBinLowEdge(1));

    // Adaptive binning in X
    int ix = 1;
    while (ix <= nx) {
        int ix_end = ix;
        double sum = 0;
        do {
            for (int iy = 1; iy <= ny; iy++) {
                sum += h2->GetBinContent(ix_end, iy);
            }
            ix_end++;
        } while (ix_end <= nx && sum < minCounts);
        xEdges.push_back(h2->GetXaxis()->GetBinUpEdge(ix_end - 1));
        ix = ix_end;
    }

    // Adaptive binning in Y
    int iy = 1;
    while (iy <= ny) {
        int iy_end = iy;
        double sum = 0;
        do {
            for (int ix2 = 1; ix2 <= nx; ix2++) {
                sum += h2->GetBinContent(ix2, iy_end);
            }
            iy_end++;
        } while (iy_end <= ny && sum < minCounts);
        yEdges.push_back(h2->GetYaxis()->GetBinUpEdge(iy_end - 1));
        iy = iy_end;
    }

    // Create new histogram
    TH2D* h2_new = new TH2D("h2_adaptive", h2->GetTitle(),
                            xEdges.size() - 1, &xEdges[0],
                            yEdges.size() - 1, &yEdges[0]);

    // Fill new histogram
    for (int ix_new = 1; ix_new <= h2_new->GetNbinsX(); ix_new++) {
        for (int iy_new = 1; iy_new <= h2_new->GetNbinsY(); iy_new++) {
            double xlow = h2_new->GetXaxis()->GetBinLowEdge(ix_new);
            double xup  = h2_new->GetXaxis()->GetBinUpEdge(ix_new);
            double ylow = h2_new->GetYaxis()->GetBinLowEdge(iy_new);
            double yup  = h2_new->GetYaxis()->GetBinUpEdge(iy_new);

            double sum = 0;
            for (int ix_old = 1; ix_old <= nx; ix_old++) {
                double xcen = h2->GetXaxis()->GetBinCenter(ix_old);
                if (xcen < xlow || xcen >= xup) continue;
                for (int iy_old = 1; iy_old <= ny; iy_old++) {
                    double ycen = h2->GetYaxis()->GetBinCenter(iy_old);
                    if (ycen < ylow || ycen >= yup) continue;
                    sum += h2->GetBinContent(ix_old, iy_old);
                }
            }
            h2_new->SetBinContent(ix_new, iy_new, sum);
        }
    }

    return h2_new;
}

// Macro entry point
void AdaptiveTailMerge2D_macro(const char* infile = "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/Chain_testing/CERN/Eventrate_q0q310e4.root", const char* histname = "myhist2", double minCounts = 1000, const char* outfile = "output.root") {
    TFile* fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: cannot open input file." << std::endl;
        return;
    }
    TH2D* h2 = dynamic_cast<TH2D*>(fin->Get(histname));
    if (!h2) {
        std::cerr << "Error: histogram not found." << std::endl;
        return;
    }

    TH2D* h2_new = AdaptiveTailMerge2D_macro(h2, minCounts);

    // Print new X bin edges
    std::cout << "New X bin edges:" << std::endl;
    for (int i = 1; i <= h2_new->GetNbinsX() + 1; ++i) {
        std::cout << h2_new->GetXaxis()->GetBinLowEdge(i) << " ";
    }
    std::cout << std::endl;

    // Print new Y bin edges
    std::cout << "New Y bin edges:" << std::endl;
    for (int i = 1; i <= h2_new->GetNbinsY() + 1; ++i) {
        std::cout << h2_new->GetYaxis()->GetBinLowEdge(i) << " ";
    }
    std::cout << std::endl;


    TFile* fout = TFile::Open(outfile, "RECREATE");
    h2_new->Write();
    fout->Close();
    fin->Close();

    std::cout << "Adaptive histogram saved to " << outfile << std::endl;
}
