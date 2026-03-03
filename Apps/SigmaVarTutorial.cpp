// #include <iostream>
// #include <chrono>
// #include <iomanip>
// #include <vector>

// #include <TH1D.h>
// #include <THStack.h>
// #include <TStyle.h>
// #include <TCanvas.h>
// #include <TRint.h>
// #include <TLegend.h>
// #include <TColor.h>
// #include <TMath.h>
// #include "TError.h"

// #include "Samples/MaCh3DUNEFactory.h"
// #include "Samples/StructsDUNE.h"
// #include "Fitters/MaCh3Factory.h"


// void PlotSigmaVariations(const std::string& filename,
//                          const std::vector<double>& sigmaVariations,
//                          const std::string& projAxis = "X") // "X" or "Y"
// {
//     TFile* file = TFile::Open(filename.c_str(), "READ");
//     if (!file || file->IsZombie()) {
//         MACH3LOG_ERROR("Could not open {}", filename);
//         return;
//     }

//     auto canvas = std::make_unique<TCanvas>("c", "c", 1000, 800);
//     canvas->SetGrid();

//     std::string pdfName = filename.substr(0, filename.find(".root")) + "_SigmaVar.pdf";
//     canvas->Print((pdfName + "[").c_str());

//     TIter parIter(file->GetListOfKeys());
//     TKey* parKey = nullptr;

//     while ((parKey = static_cast<TKey*>(parIter()))) {
//         if (std::string(parKey->GetClassName()) != "TDirectoryFile") continue;

//         TDirectory* parDir = file->GetDirectory(parKey->GetName());
//         std::string parName = parKey->GetName();

//         TIter sampIter(parDir->GetListOfKeys());
//         TKey* sampKey = nullptr;

//         while ((sampKey = static_cast<TKey*>(sampIter()))) {
//             if (std::string(sampKey->GetClassName()) != "TDirectoryFile") continue;

//             TDirectory* sampDir = parDir->GetDirectory(sampKey->GetName());
//             std::string sampName = sampKey->GetName();

//             std::vector<TH1*> hists;
//             double maxY = 0;

//             for (size_t i = 0; i < sigmaVariations.size(); ++i) {
//                 std::string histName = Form("Variation_%zu_%s", i, projAxis.c_str());
//                 TH1* hOrig = sampDir->Get<TH1>(histName.c_str());
//                 if (!hOrig) continue;

//                 TH1* h = nullptr;
//                 if (dynamic_cast<TH2*>(hOrig)) {
//                     if (projAxis == "X") h = static_cast<TH2*>(hOrig)->ProjectionX();
//                     else if (projAxis == "Y") h = static_cast<TH2*>(hOrig)->ProjectionY();
//                     h->SetName(Form("%s_proj%s", hOrig->GetName(), projAxis.c_str()));
//                 } else {
//                     h = static_cast<TH1*>(hOrig->Clone());
//                 }

//                 h->SetDirectory(nullptr);
//                 h->SetLineWidth(2);
//                 h->SetLineColor(kBlack + i);
//                 h->SetLineStyle(i == 3 ? kSolid : kDashed);

//                 maxY = std::max(maxY, h->GetMaximum());
//                 hists.push_back(h);
//             }

//             if (hists.empty()) continue;

//             canvas->Clear();
//             hists[0]->SetTitle(Form("%s – %s (Projection %s)", parName.c_str(), sampName.c_str(), projAxis.c_str()));
//             hists[0]->SetMaximum(maxY * 1.2);
//             hists[0]->Draw("HIST");

//             for (size_t i = 1; i < hists.size(); ++i)
//                 hists[i]->Draw("HIST SAME");

//             auto leg = std::make_unique<TLegend>(0.6, 0.6, 0.88, 0.88);
//             for (size_t i = 0; i < hists.size(); ++i)
//                 leg->AddEntry(hists[i],
//                               Form("%+.0f#sigma", sigmaVariations[i]),
//                               "l");

//             leg->Draw();
//             canvas->Print(pdfName.c_str());

//             for (auto* h : hists) delete h;
//         }
//     }

//     canvas->Print((pdfName + "]").c_str());
//     file->Close();
// }

// int main(int argc, char * argv[]) {
//     gErrorIgnoreLevel = kFatal;
//     auto FitManager = MaCh3ManagerFactory(argc, argv);

//     // Sigma variations in units of each parameter sigma
//     std::vector<double> sigmaVariations = {-3, -2, -1, 0, 1, 2, 3};

//     // Create SamplePDF objects
//     ParameterHandlerGeneric* xsec = nullptr;
//     std::vector<SampleHandlerFD*> DUNEPdfs;
//     MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

//     // Reweight & print integrals
//     MACH3LOG_INFO("=======================================================");
//     for (SampleHandlerFD* Sample : DUNEPdfs) {
//         Sample->Reweight();
//         MACH3LOG_INFO("Event rate for {} : {:<5.2f}", Sample->GetTitle(), Sample->GetMCHist(Sample->GetNDim())->Integral());
//     }

//     // Sigma variation loop
//     std::vector<ParameterHandlerBase*> CovObjs = { xsec };
//     MACH3LOG_INFO("=======================================================");
//     std::string OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();
//     TFile* File = TFile::Open(OutputFileName.c_str(), "RECREATE");

//     for (ParameterHandlerBase* CovObj : CovObjs) {
//         MACH3LOG_INFO("Starting Variations for covarianceBase Object: {}", CovObj->GetName());

//         int nPars = CovObj->GetNumParams();
//         for (int iPar = 0; iPar < nPars; iPar++) {
//             std::string ParName = CovObj->GetParFancyName(iPar);
//             double VarInit = CovObj->GetParInit(iPar);
//             double VarSigma = CovObj->GetDiagonalError(iPar);

//             MACH3LOG_INFO("\tParameter : {:<30} - Variations around value : {:<10.7f} , in units of 1 Sigma : {:<10.7f}",
//                           ParName, VarInit, VarSigma);

//             File->cd();
//             File->mkdir(ParName.c_str());
//             File->cd(ParName.c_str());

//             for (size_t iSigVar = 0; iSigVar < sigmaVariations.size(); iSigVar++) {
//                 double VarVal = VarInit + sigmaVariations[iSigVar] * VarSigma;
//                 VarVal = std::clamp(VarVal, CovObj->GetLowerBound(iPar), CovObj->GetUpperBound(iPar));
//                 MACH3LOG_INFO("\t\tVariation {:<5.3f} - Parameter Value : {:<10.7f}", sigmaVariations[iSigVar], VarVal);
//                 CovObj->SetParProp(iPar, VarVal);

//                 for (size_t iSample = 0; iSample < DUNEPdfs.size(); iSample++) {
//                     std::string SampleName = DUNEPdfs[iSample]->GetTitle();
//                     File->cd(ParName.c_str());
//                     if (iSigVar == 0) File->mkdir((ParName + "/" + SampleName).c_str());
//                     File->cd((ParName + "/" + SampleName).c_str());

//                     DUNEPdfs[iSample]->Reweight();
//                     TH2* Hist2D = dynamic_cast<TH2*>(DUNEPdfs[iSample]->GetMCHist(DUNEPdfs[iSample]->GetNDim()));
// if (!Hist2D) {
//     MACH3LOG_ERROR("Histogram is not 2D for sample {}", DUNEPdfs[iSample]->GetTitle());
//     continue;
// }

// // Save X projection
// TH1* histX = Hist2D->ProjectionX();
// histX->SetName(Form("Variation_%zu_X", iSigVar));
// histX->Write();
// delete histX;

// // Save Y projection
// TH1* histY = Hist2D->ProjectionY();
// histY->SetName(Form("Variation_%zu_Y", iSigVar));
// histY->Write();
// delete histY;

// MACH3LOG_INFO("\t\t\tSample : {:<30} - Integral X/Y : {:<10.2f} / {:<10.2f}",
//               SampleName, Hist2D->ProjectionX()->Integral(), Hist2D->ProjectionY()->Integral());
//                     MACH3LOG_INFO("\t\t\tSample : {:<30} - Integral X/Y : {:<10.2f} / {:<10.2f}",
//                                   SampleName, Hist2D->ProjectionX()->Integral(), Hist2D->ProjectionY()->Integral());
//                 }
//             }
//             CovObj->SetParProp(iPar, VarInit);
//         }
//         MACH3LOG_INFO("=======================================================");
//     }

//     File->Close();

//     // Plot X and Y projections
//     PlotSigmaVariations(OutputFileName, sigmaVariations, "X");
//     PlotSigmaVariations(OutputFileName, sigmaVariations, "Y");
// }

#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TSpline.h>
#include <TCanvas.h>
#include <TString.h>
#include <vector>
#include <iostream>

void MakeMAQESigmaRatios() {

    const char* splineFileName = "OffAxisND_splines.root";
    const char* nominalFileName = "offaxis_fakedataeventrate.root";
    const char* nominalHistName = "h_nominal"; // replace with your actual histogram name
    std::vector<double> sigmas = {1.0, 1.5, 2.0};

    TFile* fSpline = TFile::Open(splineFileName);
    if (!fSpline) { std::cout << "Cannot open spline file" << std::endl; return; }

    TFile* fNom = TFile::Open(nominalFileName);
    TH2D* hNom = (TH2D*)fNom->Get(nominalHistName);
    if (!hNom) { std::cout << "Cannot load nominal histogram" << std::endl; return; }

    int nx = hNom->GetNbinsX();
    int ny = hNom->GetNbinsY();

    TCanvas* c = new TCanvas("c","c",900,700);
    TString pdfName = "MAQE_SigmaRatios.pdf";
    c->Print(pdfName + "[");

    for (double sigmaVal : sigmas) {

        TH2D* hRatio = (TH2D*)hNom->Clone(Form("hRatio_%g", sigmaVal));
        hRatio->Reset();

        for (int x = 1; x <= nx; x++) {
            for (int y = 1; y <= ny; y++) {
                double wNom = hNom->GetBinContent(x,y);
                if (wNom == 0) { hRatio->SetBinContent(x,y,1.0); continue; }

                TString splineName = Form("dev_MAQE_sp_0_%d_%d", x-1, y-1); // adjust to your naming
                TSpline3* sp = (TSpline3*)fSpline->Get(splineName);
                double wVar = sp ? sp->Eval(sigmaVal) : 1.0;

                hRatio->SetBinContent(x,y, wVar ); // ratio = spline(σ)/1, multiplicative
            }
        }

        // X projection
        TH1D* projX = hRatio->ProjectionX(Form("projX_%g", sigmaVal));
        projX->SetLineColor(kRed);
        projX->SetLineWidth(2);
        projX->SetTitle(Form("MAQE Sigma Ratio = %.1f (X proj)", sigmaVal));
        c->Clear();
        projX->Draw("HIST");
        c->Print(pdfName);

        // Y projection
        TH1D* projY = hRatio->ProjectionY(Form("projY_%g", sigmaVal));
        projY->SetLineColor(kBlue);
        projY->SetLineWidth(2);
        projY->SetTitle(Form("MAQE Sigma Ratio = %.1f (Y proj)", sigmaVal));
        c->Clear();
        projY->Draw("HIST");
        c->Print(pdfName);

        delete projX;
        delete projY;
        delete hRatio;
    }

    c->Print(pdfName + "]");
    fSpline->Close();
    fNom->Close();

    std::cout << "Done! Sigma ratios plotted in " << pdfName << std::endl;
}