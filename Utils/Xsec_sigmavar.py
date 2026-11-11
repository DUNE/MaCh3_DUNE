// Plot_SigmaVariations2D_thesis.C
//
// Produces individual, thesis-quality PDF + PNG plots (one pair of files
// per parameter/sample) showing:
//   - the +1 sigma varied 2D histogram
//   - the ratio (+1 sigma / nominal)
//
// Designed for direct \includegraphics use in a LaTeX thesis:
//   \includegraphics[width=0.48\textwidth]{plots/<par>_<sample>_var.pdf}
//   \includegraphics[width=0.48\textwidth]{plots/<par>_<sample>_ratio.pdf}
//
// Usage (in ROOT):
//   root -l -b -q 'Plot_SigmaVariations2D_thesis.C("/path/to/sigmavar_recovtrue.root")'

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TLine.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <utility>

//------------------------------------------------------------------------
// Set up a clean, consistent "thesis" style. Called once at the top.
//------------------------------------------------------------------------
void SetThesisStyle()
{
    gStyle->Reset();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);

    // Font 42 = Helvetica, precision 2 (scales cleanly at any canvas size)
    const int kFont = 42;
    gStyle->SetTextFont(kFont);
    gStyle->SetLabelFont(kFont, "XYZ");
    gStyle->SetTitleFont(kFont, "XYZ");
    gStyle->SetTitleFont(kFont, "");
    gStyle->SetLegendFont(kFont);

    gStyle->SetLabelSize(0.045, "XYZ");
    gStyle->SetTitleSize(0.050, "XYZ");
    gStyle->SetTitleOffset(1.15, "X");
    gStyle->SetTitleOffset(1.35, "Y");
    gStyle->SetTitleOffset(1.35, "Z");

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadBorderMode(0);

    // Colourblind-safe, perceptually uniform palette
    gStyle->SetPalette(kViridis);
    gStyle->SetNumberContours(99);

    gStyle->SetTitleAlign(23);
    gStyle->SetTitleX(0.5);
    gStyle->SetTitleBorderSize(0);
    gStyle->SetTitleFontSize(0.055);
}

//------------------------------------------------------------------------
// Explicit physics-variable axis titles, overriding whatever is (or
// isn't) baked into the histograms themselves. Matches the MaCh3
// binning convention: VarStr = ["TrueNeutrinoEnergy", "RecoNeutrinoEnergy"]
// i.e. X = true, Y = reco, both in GeV.
//------------------------------------------------------------------------
static const char* kXAxisTitle = "True Neutrino Energy (GeV)";
static const char* kYAxisTitle = "Reconstructed Neutrino Energy (GeV)";

void SetPhysicsAxisTitles(TH2 *h)
{
    if (!h) return;
    h->GetXaxis()->SetTitle(kXAxisTitle);
    h->GetYaxis()->SetTitle(kYAxisTitle);
}

//------------------------------------------------------------------------
// Turn an arbitrary ROOT object/directory name into a safe filename.
//------------------------------------------------------------------------
TString SanitizeName(const TString &in)
{
    TString out = in;
    out.ReplaceAll(" ", "_");
    out.ReplaceAll("/", "-");
    out.ReplaceAll("\\", "-");
    out.ReplaceAll("#", "");
    out.ReplaceAll("(", "");
    out.ReplaceAll(")", "");
    out.ReplaceAll(":", "-");
    return out;
}

//------------------------------------------------------------------------
// Turn a raw sample/parameter name into a human-readable label,
// e.g. "OnAxis_numuCC_numode_true1" -> "OnAxis numuCC numode true1"
//------------------------------------------------------------------------
TString PrettifyLabel(const TString &in)
{
    TString out = in;
    out.ReplaceAll("_", " ");
    return out;
}

//------------------------------------------------------------------------
// Escape LaTeX special characters so raw names can go safely into
// \caption{} / section text.
//------------------------------------------------------------------------
TString EscapeLatex(const TString &in)
{
    TString out = in;
    out.ReplaceAll("_", "\\_");
    out.ReplaceAll("&", "\\&");
    out.ReplaceAll("%", "\\%");
    out.ReplaceAll("#", "\\#");
    return out;
}

//------------------------------------------------------------------------
// Write one (or more, if there are many samples) LaTeX figure block(s)
// for a single systematic parameter, with one subfigure per sample.
// samples: vector of (fileBaseName, prettySampleLabel) pairs, in the
// order they were found in the ROOT file.
//------------------------------------------------------------------------
void WriteLatexFigureBlock(std::ofstream &tex,
                            const TString &parName,
                            const std::vector<std::pair<TString,TString>> &samples,
                            const TString &plotDirForLatex,
                            int perFigure = 6,
                            double subfigWidth = 0.31)
{
    if (samples.empty()) return;

    const TString parLabelEsc = EscapeLatex(parName);
    const int nChunks = (samples.size() + perFigure - 1) / perFigure;

    for (int chunk = 0; chunk < nChunks; ++chunk) {
        const int start = chunk * perFigure;
        const int end   = std::min((int)samples.size(), start + perFigure);

        tex << "\\begin{figure}[htbp]\n";
        tex << "    \\centering\n";

        for (int i = start; i < end; ++i) {
            const TString &base   = samples[i].first;
            const TString &sample = samples[i].second;

            tex << "    \\begin{subfigure}[b]{" << subfigWidth << "\\textwidth}\n";
            tex << "        \\centering\n";
            tex << "        \\includegraphics[width=\\textwidth]{"
                << plotDirForLatex << "/" << base << "_ratio.pdf}\n";
            tex << "        \\caption{" << EscapeLatex(sample) << "}\n";
            tex << "        \\label{fig:" << base << "_ratio}\n";
            tex << "    \\end{subfigure}\n";

            // Line break every 3 subfigures for a tidy grid, unless
            // it's the last one in this chunk.
            if ((i - start + 1) % 3 == 0 && i != end - 1) {
                tex << "    \\\\[1em]\n";
            } else if (i != end - 1) {
                tex << "    \\hfill\n";
            }
        }

        tex << "    \\caption{Effect of a $+1\\sigma$ shift in the " << parLabelEsc
            << " cross-section systematic on the reconstructed-vs-true neutrino "
               "energy distribution, shown as the ratio to nominal for each sample"
            << (nChunks > 1 ? Form(" (%d/%d)", chunk + 1, nChunks) : "") << ".}\n";
        tex << "    \\label{fig:" << SanitizeName(parName) << "_ratio_all"
            << (nChunks > 1 ? Form("_%d", chunk + 1) : "") << "}\n";
        tex << "\\end{figure}\n\n";
    }
}

void Plot_SigmaVariations2D_thesis(
    const char* filename    = "/scratch/abipeake/MaCh3_DUNE_Nov2025/MaCh3_DUNE/sigmavar_recovtrue.root",
    const char* outdir      = "plots_thesis",
    int subfiguresPerFigure = 6,   // how many sample subfigures per LaTeX figure before it wraps to a new one
    const char* texRelativePlotDir = "plots_thesis") // path to use inside \includegraphics{} - adjust to match your thesis directory structure
{
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cout << "File not found: " << filename << "\n";
        return;
    }

    gSystem->mkdir(outdir, /*recursive=*/true);

    SetThesisStyle();

    TString texPath = Form("%s/systematics_subfigures.tex", outdir);
    std::ofstream tex(texPath.Data());
    tex << "% Auto-generated by Plot_SigmaVariations2D_thesis.C\n";
    tex << "% One figure per cross-section systematic, one subfigure per sample,\n";
    tex << "% showing the +1 sigma / nominal ratio.\n";
    tex << "%\n";
    tex << "% Requires in your thesis preamble:\n";
    tex << "%   \\usepackage{graphicx}\n";
    tex << "%   \\usepackage{subcaption}\n";
    tex << "% Recommended if figures get separated from their section by other floats:\n";
    tex << "%   \\usepackage{placeins}  % then add \\FloatBarrier between sections\n";
    tex << "%\n";
    tex << "% Include in your thesis with:\n";
    tex << "%   \\input{" << texPath << "}\n\n";

    // Canvas: tall, split into a main panel (variation) and a smaller
    // ratio panel underneath, in the style of ATLAS/CMS/DUNE plots.
    TCanvas *c = new TCanvas("c", "c", 800, 900);
    c->cd();

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 1.0);
    pad1->SetTopMargin(0.10);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.14);
    pad1->SetRightMargin(0.17);
    pad1->Draw();

    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.35);
    pad2->SetTopMargin(0.04);
    pad2->SetBottomMargin(0.30);
    pad2->SetLeftMargin(0.14);
    pad2->SetRightMargin(0.17);
    pad2->Draw();

    TIter parDirIter(f->GetListOfKeys());
    TKey *parKey;

    int nPlots = 0;

    while ((parKey = (TKey*)parDirIter())) {

        TString parName = parKey->GetName();
        TDirectory *parDir = (TDirectory*)f->Get(parName);
        if (!parDir) continue;

        // Collected as we go, then written out as one LaTeX figure
        // (or several, if there are many samples) once this parameter
        // is done.
        std::vector<std::pair<TString,TString>> samplesForLatex;

        TIter sampleIter(parDir->GetListOfKeys());
        TKey *sampleKey;

        while ((sampleKey = (TKey*)sampleIter())) {

            TString sampleName = sampleKey->GetName();
            TDirectory *sampleDir = (TDirectory*)parDir->Get(sampleName);
            if (!sampleDir) continue;

            // Nominal (sigma = 0) and +1 sigma
            TH2 *hNom  = (TH2*)sampleDir->Get("Variation_2");
            TH2 *hPlus = (TH2*)sampleDir->Get("Variation_3");
            if (!hNom || !hPlus) continue;

            TH2 *hVar   = (TH2*)hPlus->Clone("hVar");
            TH2 *hRatio = (TH2*)hPlus->Clone("hRatio");
            hVar->SetDirectory(0);
            hRatio->SetDirectory(0);
            hRatio->Divide(hNom);

            TString base = SanitizeName(parName) + "_" + SanitizeName(sampleName);

            // ---------------- Variation panel ----------------
            pad1->cd();
            pad1->Clear();
            pad1->SetLogz(0);

            SetPhysicsAxisTitles(hVar);
            hVar->SetTitle(Form("%s: %s (+1#sigma)", parName.Data(), sampleName.Data()));
            hVar->GetXaxis()->SetLabelSize(0); // hide x labels, ratio panel shows them
            hVar->GetXaxis()->SetTitleSize(0);
            hVar->GetZaxis()->SetTitle("Events");
            hVar->GetZaxis()->SetTitleOffset(1.25);
            hVar->Draw("COLZ");
            pad1->Update();

            // ---------------- Ratio panel ----------------
            pad2->cd();
            pad2->Clear();

            SetPhysicsAxisTitles(hRatio);
            hRatio->SetTitle("");
            hRatio->GetXaxis()->SetLabelSize(0.09);
            hRatio->GetXaxis()->SetTitleSize(0.10);
            hRatio->GetYaxis()->SetLabelSize(0.09);
            hRatio->GetYaxis()->SetTitleSize(0.10);
            hRatio->GetYaxis()->SetTitleOffset(0.65);
            hRatio->GetZaxis()->SetTitle("+1#sigma / Nominal");
            hRatio->GetZaxis()->SetLabelSize(0.08);
            hRatio->GetZaxis()->SetTitleSize(0.09);
            hRatio->GetZaxis()->SetTitleOffset(0.65);

            // Symmetric range centred on 1, sized to the actual spread
            // of the ratio (falls back to a sensible default if the
            // ratio is essentially flat).
            double zMin = hRatio->GetMinimum();
            double zMax = hRatio->GetMaximum();
            double halfRange = std::max({std::abs(1.0 - zMin), std::abs(zMax - 1.0), 0.02});
            hRatio->SetMinimum(1.0 - halfRange);
            hRatio->SetMaximum(1.0 + halfRange);

            hRatio->Draw("COLZ");
            pad2->Update();

            // Save individual files
            TString pdfVar   = Form("%s/%s_var.pdf", outdir, base.Data());
            TString pngVar   = Form("%s/%s_var.png", outdir, base.Data());
            TString pdfRatio = Form("%s/%s_ratio.pdf", outdir, base.Data());
            TString pngRatio = Form("%s/%s_ratio.png", outdir, base.Data());

            // Save the variation panel alone
            {
                TCanvas cVar("cVar", "cVar", 800, 650);
                cVar.SetTopMargin(0.10);
                cVar.SetBottomMargin(0.13);
                cVar.SetLeftMargin(0.14);
                cVar.SetRightMargin(0.17);
                TH2 *hVarStandalone = (TH2*)hVar->Clone("hVarStandalone");
                hVarStandalone->SetDirectory(0);
                SetPhysicsAxisTitles(hVarStandalone);
                hVarStandalone->GetXaxis()->SetLabelSize(0.045);
                hVarStandalone->GetXaxis()->SetTitleSize(0.050);
                hVarStandalone->Draw("COLZ");
                cVar.Print(pdfVar);
                cVar.Print(pngVar);
                delete hVarStandalone;
            }

            // Save the ratio panel alone
            {
                TCanvas cRatio("cRatio", "cRatio", 800, 650);
                cRatio.SetTopMargin(0.10);
                cRatio.SetBottomMargin(0.13);
                cRatio.SetLeftMargin(0.14);
                cRatio.SetRightMargin(0.17);
                TH2 *hRatioStandalone = (TH2*)hRatio->Clone("hRatioStandalone");
                hRatioStandalone->SetDirectory(0);
                SetPhysicsAxisTitles(hRatioStandalone);
                hRatioStandalone->SetTitle(Form("%s: %s (+1#sigma / Nominal)", parName.Data(), sampleName.Data()));
                hRatioStandalone->GetXaxis()->SetLabelSize(0.045);
                hRatioStandalone->GetXaxis()->SetTitleSize(0.050);
                hRatioStandalone->GetYaxis()->SetLabelSize(0.045);
                hRatioStandalone->GetYaxis()->SetTitleSize(0.050);
                hRatioStandalone->GetYaxis()->SetTitleOffset(1.35);
                hRatioStandalone->Draw("COLZ");
                cRatio.Print(pdfRatio);
                cRatio.Print(pngRatio);
                delete hRatioStandalone;
            }

            // Also save the combined (stacked) view, matching the original layout
            TString pdfCombined = Form("%s/%s_combined.pdf", outdir, base.Data());
            c->Print(pdfCombined);

            samplesForLatex.push_back({base, PrettifyLabel(sampleName)});

            delete hVar;
            delete hRatio;
            ++nPlots;
        }

        // Write this parameter's LaTeX figure block(s) now that all its
        // sample ratio plots exist.
        WriteLatexFigureBlock(tex, parName, samplesForLatex,
                               texRelativePlotDir, subfiguresPerFigure);
    }

    tex.close();

    std::cout << "Wrote " << nPlots << " plot set(s) to \"" << outdir << "/\"\n";
    std::cout << "Wrote LaTeX subfigures to \"" << texPath << "\"\n";
    f->Close();
}