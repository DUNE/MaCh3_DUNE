#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <iostream>

void make1DNormalisedHistRatio() {
    // Open the ROOT file
    TFile *file = TFile::Open("./outputs/NDGAr_EventRates_CC_BField0_5_FidRad_160_FidLen_209_InstrRad_249_45_InstrLen_259_WithECal_Enu_1to5GeV.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open the file!" << std::endl;
        return;
    }

    // Get the histograms
    TH1 *hist1 = dynamic_cast<TH1*>(file->Get("FHC_numu_Interaction Mode (bad caf events)"));
    TH1 *hist2 = dynamic_cast<TH1*>(file->Get("FHC_numu_Interaction Mode"));

    if (!hist1 || !hist2) {
        std::cerr << "Error: Could not find 'hist1' or 'hist2' in the file." << std::endl;
        return;
    }

    // Normalize the histograms
    if (hist1->Integral() != 0)
        hist1->Scale(1.0 / hist1->Integral());
    else
        std::cerr << "Warning: hist1 has zero integral and cannot be normalized." << std::endl;

    if (hist2->Integral() != 0)
        hist2->Scale(1.0 / hist2->Integral());
    else
        std::cerr << "Warning: hist2 has zero integral and cannot be normalized." << std::endl;

    // Create the ratio histogram
    TH1 *ratio = dynamic_cast<TH1*>(hist1->Clone("ratio"));
    ratio->Divide(hist2);

    // Draw the ratio
    TCanvas *canvas = new TCanvas("canvas", "Bad CAF events / All events (normalised)", 800, 600);
    ratio->SetTitle("Bad CAF events / All events (normalised);Mode;Ratio");
    ratio->Draw("E1");

    // Save to PDF
    canvas->SaveAs("Interaction_mode_ratio.pdf");

    // Clean up
    file->Close();
    delete file;
    delete canvas;

    return;
}
