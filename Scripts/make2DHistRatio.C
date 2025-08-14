#include <TFile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <iostream>
#include "DUNEStyle.h"

void make2DHistRatio() {
  // Open the ROOT files
  TFile *file1 = TFile::Open("./Outputs/projections_outputs/Projections_CC_BField0_5_FidRad_160_FidLen_209_InstRad_249_45_InstLen_259_WithECal_seccurve_2layers_2MeV.root");
  TFile *file2 = TFile::Open("./Outputs/projections_outputs/Projections_CC_BField0_5_FidRad_160_FidLen_209_InstRad_249_45_InstLen_259_WithECal_2layers_2MeV.root");
  if (!file1 || file1->IsZombie() || !file2 || file2->IsZombie()) {
    std::cerr << "Error: Cannot open the file!" << std::endl;
    return;
  }

  // Get the histograms
  // TH2D *hist1 = dynamic_cast<TH2D*>(file1->Get("FHC_numu_NDGAr_Accepted_TrueQ3_vs_TrueQ0"));
  // TH2D *hist2 = dynamic_cast<TH2D*>(file2->Get("FHC_numu_NDGAr_Accepted_TrueQ3_vs_TrueQ0"));
  TH2D *hist1 = dynamic_cast<TH2D*>(file1->Get("FHC_numu_NDGAr_Accepted_Particle_BAngle_Momentum_Pi0"));
  TH2D *hist2 = dynamic_cast<TH2D*>(file2->Get("FHC_numu_NDGAr_Accepted_Particle_BAngle_Momentum_Pi0"));

  if (!hist1 || !hist2) {
    std::cerr << "Error: Could not find 'hist1' or 'hist2' in the file." << std::endl;
    return;
  }

  // Create the ratio histogram
  TH1 *ratio = dynamic_cast<TH2D*>(hist1->Clone("ratio"));
  ratio->Divide(hist2);

  // Draw the ratio
  TCanvas *canvas = new TCanvas("canvas", "Sec Curvature Acceptance Improvement pi+", 800, 600);
  ratio->SetTitle("Acceptance With Sec Curvature / Without Sec Curvature;Angle to B-Field [#circ];Momentum [GeV/c]");
  ratio->Draw("COLZ");

  // Save to PDF
  canvas->SaveAs("SecCurveImprovement_pi0bangle.pdf");

  // Clean up
  file1->Close();
  delete file1;
  file2->Close();
  delete file2;

  delete canvas;

  return;
}
