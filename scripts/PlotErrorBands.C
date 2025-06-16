#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

void PlotErrorBands(const std::string& filename, const std::string& histBaseName = "posterior_prediction_sample", int ntoys = 500, int pdfIndex = 0) {
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "[ERROR] Failed to open file: " << filename << std::endl;
        return;
    }

    std::vector<TH1D*> toys;
    int nFound = 0;

    for (int i = 0; i < ntoys; ++i) {
        std::string histName = histBaseName + std::to_string(i) + "_" + std::to_string(pdfIndex); //posterior_prediction_sample487_pdf0

        TH1D* h = dynamic_cast<TH1D*>(file->Get(histName.c_str()));
        
        if (!h) {
            std::cerr << "[WARNING] Histogram not found: " << histName << std::endl;
            continue;
        }

        toys.push_back((TH1D*)h->Clone());
        ++nFound;
        if (nFound <= 3) {
            std::cout << "[DEBUG] Loaded " << histName << " with integral: " << h->Integral() << std::endl;
        }
    }

    if (toys.empty()) {
        std::cerr << "[ERROR] No valid toy histograms loaded." << std::endl;
        return;
    }

    std::cout << "[INFO] Successfully loaded " << nFound << " toy histograms." << std::endl;

    int nbins = toys[0]->GetNbinsX();
    std::vector<double> mean(nbins, 0.0);
    std::vector<double> stddev(nbins, 0.0);

    // Calculate per-bin mean and variance
    for (int bin = 1; bin <= nbins; ++bin) {
        ////////Check the toys
        std::cout << "Bin " << bin << " values:\n";
        for (int i = 0; i < std::min((int)toys.size(), 5); ++i) {
            std::cout << "  Toy " << i << ": " << toys[i]->GetBinContent(bin) << std::endl;
        }
                double sum = 0.0;
        for (auto& h : toys) sum += h->GetBinContent(bin);
        double avg = sum / toys.size();
        mean[bin - 1] = avg;

        double var = 0.0;
        for (auto& h : toys) {
            double diff = h->GetBinContent(bin) - avg;
            var += diff * diff;
        }
        var /= (toys.size() - 1);  // unbiased variance
        stddev[bin - 1] = std::sqrt(var);

        //std::cout << "Bin = " << bin << ": mean = " << avg << ", stddev = " << std::sqrt(var) << std::endl;
    }

    // Create error band graph
    TGraphAsymmErrors* band = new TGraphAsymmErrors(nbins);
    for (int bin = 1; bin <= nbins; ++bin) {
        double x = toys[0]->GetBinCenter(bin);
        double y = mean[bin - 1];
        double err = stddev[bin - 1];

        band->SetPoint(bin - 1, x, y);
        band->SetPointError(bin - 1, 0, 0, err, err);
    }

    band->SetFillColorAlpha(kBlue, 0.3);
    band->SetLineColor(kBlue + 2);
    band->SetMarkerStyle(0);
    band->SetTitle("±1σ Error Band");

    // Plotting
    TCanvas* c = new TCanvas("c", "Error Bands", 800, 600);
    toys[0]->SetLineColor(kRed + 1);
    toys[0]->SetLineWidth(2);
    toys[0]->SetTitle("Prediction with Error Bands");
    toys[0]->GetXaxis()->SetTitle("Reconstructed Neutrino Energy");
    toys[0]->GetXaxis()->SetRangeUser(0, 30);
    toys[0]->Draw("HIST");
    band->Draw("E3 SAME");

    TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(toys[0], "Nominal Prediction", "l");
    leg->AddEntry(band, "#pm1#sigma Error Band", "f");
    leg->Draw();

    c->SaveAs("error_bands.pdf");
    delete c;
    file->Close();
}
