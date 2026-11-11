#include "Samples/BinningHandler.h"
#include "Samples/SampleHandlerBase.h"
#include "Samples/SampleHandlerFD.h"
#include "Samples/MaCh3DUNEFactory.h"
#include "Manager/MaCh3Logger.h"
#include "Manager/Manager.h"
#include "Parameters/ParameterHandlerBase.h"

#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <filesystem>
#include <set>
#include <sstream>
#include <ctime>
#include <random>
#include <fstream>

#include "yaml-cpp/yaml.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLatex.h"

struct BinDef {
    int index;
    double q0_min, q0_max;
    double q3_min, q3_max;
};

struct BinningResult {
    std::vector<BinDef> binDefs;
    std::vector<double> q0_edges;
    std::vector<double> q3_edges;
};

std::string AddTimestampToROOTFilename(const std::string& baseName) {
    time_t now = time(0);
    tm* ltm = localtime(&now);
    std::ostringstream oss;
    oss << baseName.substr(0, baseName.find_last_of(".")) << "_"
        << 1900 + ltm->tm_year
        << 1 + ltm->tm_mon
        << ltm->tm_mday << "_"
        << ltm->tm_hour
        << ltm->tm_min
        << ltm->tm_sec << ".root";
    return oss.str();
}

// -------------------------
// Palettes
// -------------------------
void setViridisPalette() {
    gStyle->SetPalette(kViridis);
    gStyle->SetNumberContours(255);
}

void setBlueWhiteRedPalette(double whiteStop) {
    const Int_t nRGBs = 3;
    Double_t stops[nRGBs] = {0.00, whiteStop, 1.00};
    Double_t red[nRGBs]   = {0.00, 1.00, 1.00};
    Double_t green[nRGBs] = {0.00, 1.00, 0.00};
    Double_t blue[nRGBs]  = {1.00, 1.00, 0.00};
    const Int_t nColors = 255;
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, nColors);
    gStyle->SetNumberContours(nColors);
}

// -------------------------
// Thesis-quality style
// -------------------------
void setThesisStyle() {
    gStyle->Reset();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);           // titles drawn manually via TLatex for full control

    // Fonts: 132 = Times, precision 2 (scales cleanly at any zoom).
    // Use 42 instead of 132 everywhere below if your thesis uses a sans-serif font.
    const int font = 132;
    gStyle->SetTextFont(font);
    gStyle->SetLabelFont(font, "XYZ");
    gStyle->SetTitleFont(font, "XYZ");
    gStyle->SetLegendFont(font);

    gStyle->SetLabelSize(0.035, "XYZ");
    gStyle->SetTitleSize(0.045, "XYZ");
    gStyle->SetTitleOffset(1.1, "X");
    gStyle->SetTitleOffset(1.3, "Y");
    gStyle->SetTitleOffset(1.4, "Z");

    // Margins: leave room for the z-axis colour bar title
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.18);

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTickLength(0.02, "XYZ");

    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
}

// Draw a manual title above the frame (replaces SetOptTitle(1) titles)
void drawTitle(const std::string& text) {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(132);
    latex.SetTextSize(0.045);
    latex.SetTextAlign(21); // centered
    latex.DrawLatex(0.5, 0.94, text.c_str());
}

// -------------------------
// Style helper for 1D error-bar plots (thesis-quality)
// -------------------------
void styleThesis1D(TH1D* h) {
    h->SetLineColor(kBlack);
    h->SetLineWidth(2);
    h->SetMarkerStyle(20);      // filled circle
    h->SetMarkerSize(0.9);
    h->SetMarkerColor(kBlack);
    h->SetStats(0);
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetNdivisions(505);
    // Give a bit of vertical headroom above the tallest bin (incl. error bar)
    double ymax = h->GetMaximum() + h->GetBinError(h->GetMaximumBin());
    h->GetYaxis()->SetRangeUser(0.0, 1.15 * ymax);
}

// Divide every bin (content AND error) by its bin width, turning an event
// count histogram into a differential ("per unit x") distribution. This is
// the standard treatment for variable-width bins so that bin height is
// directly comparable across bins of different widths.
void normalizeByBinWidth(TH1D* h) {
    h->Scale(1.0, "width");
}

// -------------------------
// LaTeX-style axis labels for known MaCh3-DUNE kinematic variables
// -------------------------
std::string toLatexLabel(const std::string& var) {
    if (var == "TrueNeutrinoEnergy") return "E_{#nu}^{True} [GeV]";
    if (var == "Enubias")            return "E_{#nu}^{Bias} [GeV]";
    // add further mappings here as new variables come up
    return var; // fallback: raw variable name
}

// -------------------------
// Helper: get latest non-empty posteriors tree
// -------------------------
TTree* getLatestNonEmpty(TFile* f) {
    int cycle = 9999;
    while (cycle > 0) {
        TString name = Form("posteriors;%d", cycle);
        TTree* t = (TTree*) f->Get(name);
        if (t && t->GetEntries() > 0) return t;
        cycle--;
    }
    return nullptr;
}

// ---------------------------------------------------------------------------
// Quantile-based rebinning support
// ---------------------------------------------------------------------------
// Given a filled histogram, return nBins+1 edges such that each bin holds
// ~1/nBins of the total (weighted) content, using ROOT's TH1::GetQuantiles.
// ---------------------------------------------------------------------------
std::vector<double> ComputeQuantileEdges(TH1D* h, int nBins) {
    int nProbs = nBins + 1;
    std::vector<double> probs(nProbs);
    for (int i = 0; i < nProbs; ++i) probs[i] = static_cast<double>(i) / nBins;

    std::vector<double> edges(nProbs);
    h->GetQuantiles(nProbs, edges.data(), probs.data());

    // Snap the outer edges to the histogram's actual axis range, since
    // GetQuantiles' interpolation can slightly under/overshoot depending
    // on the granularity of the input histogram.
    edges.front() = h->GetXaxis()->GetXmin();
    edges.back()  = h->GetXaxis()->GetXmax();
    return edges;
}

// ---------------------------------------------------------------------------
// Build a finely-binned 1D histogram of `var` (summed over all DUNEPdfs and
// all samples), used purely as input to ComputeQuantileEdges. `nFineBins`
// should be much larger than the number of final quantile bins so the
// quantile estimate isn't limited by histogram granularity.
// ---------------------------------------------------------------------------
TH1D* BuildFine1DHist(std::vector<SampleHandlerFD*>& DUNEPdfs,
                       const std::string& var,
                       double lo, double hi,
                       int nFineBins = 2000) {
    // IMPORTANT: build the axis from an explicit edges array, matching the
    // pattern used for axisX/axisY elsewhere in this file (via
    // TAxis(edges.size()-1, edges.data())), NOT the uniform
    // TAxis(nBins, lo, hi) constructor. Get1DVarHist/Get2DVarHist appear to
    // rely on an explicit bin-edges array being present (e.g. checking
    // axis->GetXbins()->GetSize()); a uniformly-constructed TAxis has no
    // such array and was observed to collapse to effectively one bin.
    std::vector<double> fineEdges(nFineBins + 1);
    double width = (hi - lo) / nFineBins;
    for (int i = 0; i <= nFineBins; ++i) fineEdges[i] = lo + i * width;

    TAxis fineAxis(static_cast<int>(fineEdges.size()) - 1, fineEdges.data());

    auto* hFine = new TH1D("hFineQuantileInput", "", nFineBins, lo, hi);
    hFine->SetDirectory(nullptr);
    hFine->Sumw2();

    std::vector<KinematicCut> SelectionVector;
    for (auto& pdf : DUNEPdfs) {
        for (int s = 0; s < pdf->GetNsamples(); ++s) {
            TH1* h = pdf->Get1DVarHist(s, var, SelectionVector, 0, &fineAxis);
            if (!h) {
                std::cerr << "[Warning] Get1DVarHist returned null for sample "
                          << s << ", variable " << var << "\n";
                continue;
            }
            TH1D* h1d = dynamic_cast<TH1D*>(h);
            if (!h1d) { delete h; continue; }

            // Sanity check: bail loudly if the returned histogram doesn't
            // actually have the fine binning we asked for, rather than
            // silently producing garbage quantiles again.
            if (h1d->GetNbinsX() != nFineBins) {
                std::cerr << "[Warning] Get1DVarHist returned " << h1d->GetNbinsX()
                          << " bins for sample " << s << ", expected " << nFineBins
                          << " -- binning request may not have been honoured.\n";
            }

            hFine->Add(h1d);
            delete h1d;
        }
    }

    // Extra diagnostic: warn if the content is suspiciously concentrated,
    // so a future mis-binning is caught immediately rather than silently
    // producing bad quantile edges again.
    int maxBin = hFine->GetMaximumBin();
    double maxFrac = (hFine->Integral() > 0)
        ? hFine->GetBinContent(maxBin) / hFine->Integral() : 0.0;
    if (maxFrac > 0.5) {
        std::cerr << "[Warning] " << maxFrac * 100.0 << "% of " << var
                  << " content is in a single fine bin ([" 
                  << hFine->GetBinLowEdge(maxBin) << ", "
                  << hFine->GetBinLowEdge(maxBin) + hFine->GetBinWidth(maxBin)
                  << "]) -- this usually indicates a binning bug rather than\n"
                  << "                  genuine physics. Check Get1DVarHist's handling of the axis argument.\n";
    }

    return hFine;
}

// ---------------------------------------------------------------------------
// Write new quantile-based bin edges to a standalone text file, formatted
// as copy-pasteable YAML KinematicCuts blocks (Enubias/var2-style), for
// dropping into a regenerated templateparams.yaml.
// ---------------------------------------------------------------------------
void WriteQuantileEdgesToFile(const std::vector<double>& edges,
                               const std::string& varName,
                               const std::string& outPath) {
    std::ofstream out(outPath);
    if (!out) {
        std::cerr << "[Warning] Could not open " << outPath << " for writing.\n";
        return;
    }
    out << "# Quantile-based bin edges for " << varName
        << " (" << edges.size() - 1 << " bins), generated by mcmc_diagnostics.\n";
    out << "# Plain edge list:\n";
    out << "# [";
    for (size_t i = 0; i < edges.size(); ++i) {
        out << std::setprecision(6) << edges[i];
        if (i + 1 < edges.size()) out << ", ";
    }
    out << "]\n\n";

    out << "# YAML KinematicCuts blocks (one per bin):\n";
    for (size_t i = 0; i + 1 < edges.size(); ++i) {
        out << "    - " << varName << ":\n"
            << "      - " << std::setprecision(6) << edges[i]   << "\n"
            << "      - " << std::setprecision(6) << edges[i+1] << "\n";
    }
    out.close();
    std::cout << "[Info] Wrote " << (edges.size() - 1) << " quantile bin edges for "
              << varName << " to " << outPath << "\n";
}

BinningResult extract_2D_bins_from_yaml(const std::string& yaml_file,
                                        const std::string& xsecvar1,
                                        const std::string& xsecvar2) {
    BinningResult result;
    std::set<double> q0_set;
    std::set<double> q3_set;

    try {
        YAML::Node config = YAML::LoadFile(yaml_file);
        const auto& systs = config["Systematics"];
        int index = 0;

        for (const auto& systematic : systs) {
            const auto& sys = systematic["Systematic"];
            const auto& cuts = sys["KinematicCuts"];

            double q0_min = 0, q0_max = 0;
            double q3_min = 0, q3_max = 0;
            bool has_q0 = false, has_q3 = false;

            for (const auto& cut : cuts) {
                if (cut[xsecvar1] && cut[xsecvar1].IsSequence() && cut[xsecvar1].size() == 2) {
                    double q0_lo = cut[xsecvar1][0].as<double>();
                    double q0_hi = cut[xsecvar1][1].as<double>();
                    q0_set.insert(q0_lo);
                    q0_set.insert(q0_hi);
                    q0_min = q0_lo;
                    q0_max = q0_hi;
                    has_q0 = true;
                }
                if (cut[xsecvar2] && cut[xsecvar2].IsSequence() && cut[xsecvar2].size() == 2) {
                    double q3_lo = cut[xsecvar2][0].as<double>();
                    double q3_hi = cut[xsecvar2][1].as<double>();
                    q3_set.insert(q3_lo);
                    q3_set.insert(q3_hi);
                    q3_min = q3_lo;
                    q3_max = q3_hi;
                    has_q3 = true;
                }
            }

            if (has_q0 || has_q3) {
                result.binDefs.push_back({index++, q0_min, q0_max, q3_min, q3_max});
            }
        }

        result.q0_edges.assign(q0_set.begin(), q0_set.end());
        result.q3_edges.assign(q3_set.begin(), q3_set.end());

        std::cout << "Parsed " << result.binDefs.size() << " 2D bins\n";
        std::cout << "   q0 edges: " << (result.q0_edges.empty() ? 0 : result.q0_edges.size() - 1) << "\n";
        std::cout << "   q3 edges: " << (result.q3_edges.empty() ? 0 : result.q3_edges.size() - 1) << "\n";

    } catch (const std::exception& e) {
        std::cerr << "YAML parsing failed: " << e.what() << "\n";
    }

    return result;
}


int main(int argc, char* argv[]) {

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <mcmc_output.root> <config.yaml> [--quantile-bins N] [--quantile-var var1|var2]"
                  << std::endl;
        std::cerr << "  --quantile-bins N     Also compute an N-bin quantile-based binning for\n"
                  << "                        the chosen variable from the raw MC distribution,\n"
                  << "                        and produce a population-balance comparison against\n"
                  << "                        the current (YAML) binning. Does NOT require a\n"
                  << "                        posterior refit -- it is an MC-only diagnostic.\n";
        std::cerr << "  --quantile-var        Which variable to rebin: 'var1' (xsec_var1) or\n"
                  << "                        'var2' (xsec_var2, e.g. Enubias). Default: var2.\n";
        return 1;
    }

    setThesisStyle();

    std::string mcmc_file   = argv[1];
    std::string config_file = argv[2];

    int quantileBins = 0;              // 0 => quantile rebinning disabled
    std::string quantileVarChoice = "var2";
    for (int i = 3; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--quantile-bins" && i + 1 < argc) {
            quantileBins = std::atoi(argv[++i]);
        } else if (arg == "--quantile-var" && i + 1 < argc) {
            quantileVarChoice = argv[++i];
        }
    }

    std::filesystem::path mcmcPath(mcmc_file);
std::string outDir = mcmcPath.parent_path().string();
if (!outDir.empty()) outDir += "/";

std::string rootOut = outDir + AddTimestampToROOTFilename("mcmc_diagnostics.root");
    std::string pdfOut  = mcmc_file.substr(0, mcmc_file.find_last_of('.')) + "_diagnostics.pdf";

    // --- Open MCMC file ---
    TFile* f_post = TFile::Open(mcmc_file.c_str(), "READ");
    if (!f_post || f_post->IsZombie()) {
        std::cerr << "Cannot open MCMC file: " << mcmc_file << std::endl;
        return 1;
    }

    TTree* post = getLatestNonEmpty(f_post);
    if (!post) {
        std::cerr << "No non-empty 'posteriors' tree found in file." << std::endl;
        return 1;
    }
    std::cout << "Using " << post->GetName()
              << " with " << post->GetEntries() << " entries." << std::endl;

    // --- Load YAML config ---
    auto FitManager = std::make_unique<Manager>(config_file);
    auto xsec_var1 = FitManager->raw()["General"]["Systematics"]["xsec_var1"].as<std::string>();
    auto xsec_var2 = FitManager->raw()["General"]["Systematics"]["xsec_var2"].as<std::string>();
    auto xsec_yaml = FitManager->raw()["General"]["Systematics"]["XsecCovFile"][0].as<std::string>();

    if (xsec_yaml.empty()) {
        std::cerr << "Error: xsec_yaml path is empty or missing in config.\n";
        return 1;
    }

    // --- Parse binning ---
    BinningResult binning = extract_2D_bins_from_yaml(xsec_yaml, xsec_var1, xsec_var2);
    std::vector<BinDef>& binDefs  = binning.binDefs;
    std::vector<double>& q0_edges = binning.q0_edges;
    std::vector<double>& q3_edges = binning.q3_edges;

    if (q0_edges.empty() || q3_edges.empty()) {
        std::cerr << "Bin edges could not be parsed from YAML file." << std::endl;
        return 1;
    }

    TAxis* axisX = new TAxis(q0_edges.size() - 1, q0_edges.data());
    TAxis* axisY = new TAxis(q3_edges.size() - 1, q3_edges.data());

    // --- Setup MaCh3 PDFs ---
    ParameterHandlerGeneric* xsec = nullptr;
    std::vector<SampleHandlerFD*> DUNEPdfs;
    MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

    // --- Build 2D event-rate histogram ---
    TH2D* h_xsecvar_eventratehist = nullptr;
    for (auto& pdf : DUNEPdfs) {
        std::vector<KinematicCut> SelectionVector;
        TH2* h = pdf->Get2DVarHist(0, xsec_var1, xsec_var2, SelectionVector, 0, axisX, axisY);
        if (!h) {
            std::cerr << "Warning: Get2DVarHist returned null.\n";
            continue;
        }
        if (!h_xsecvar_eventratehist) {
            h_xsecvar_eventratehist = (TH2D*) h->Clone("h_xsecvar_eventratehist");
            h_xsecvar_eventratehist->SetDirectory(nullptr);
        } else {
            h_xsecvar_eventratehist->Add(h);
        }
    }

    if (!h_xsecvar_eventratehist) {
        std::cerr << "Error: No 2D histograms were produced.\n";
        return 1;
    }

    // Make sure per-bin errors are tracked so that ProjectionX/Y propagate
    // statistical uncertainty correctly. If Sumw2 was never called upstream
    // (e.g. the PDF fills unweighted), this call sets each bin's error to
    // sqrt(content), i.e. the standard Poisson counting uncertainty.
    if (h_xsecvar_eventratehist->GetSumw2N() == 0) {
        h_xsecvar_eventratehist->Sumw2();
    }

    h_xsecvar_eventratehist->GetXaxis()->SetTitle(toLatexLabel(xsec_var1).c_str());
    h_xsecvar_eventratehist->GetYaxis()->SetTitle(toLatexLabel(xsec_var2).c_str());
    h_xsecvar_eventratehist->GetZaxis()->SetTitle("Events");

    // --- 1D projections of the event-rate histogram, with stat. uncertainty ---
    // ProjectionX sums over q3 (var2) for each q0 (var1) bin; ProjectionY does the reverse.
    // Errors are propagated automatically from the 2D Sumw2 array above.
    TH1D* h_proj_q0 = h_xsecvar_eventratehist->ProjectionX("h_proj_q0");
    TH1D* h_proj_q3 = h_xsecvar_eventratehist->ProjectionY("h_proj_q3");

    // Divide by bin width (Scale propagates the errors correctly too),
    // since the q0/q3 binning is non-uniform.
    normalizeByBinWidth(h_proj_q0);
    normalizeByBinWidth(h_proj_q3);

    h_proj_q0->GetXaxis()->SetTitle(toLatexLabel(xsec_var1).c_str());
    h_proj_q0->GetYaxis()->SetTitle("Events / bin width");
    styleThesis1D(h_proj_q0);

    h_proj_q3->GetXaxis()->SetTitle(toLatexLabel(xsec_var2).c_str());
    h_proj_q3->GetYaxis()->SetTitle("Events / bin width");
    styleThesis1D(h_proj_q3);

    // -------------------------------------------------------------------
    // Optional: quantile-based rebinning diagnostic (MC population only).
    //
    // IMPORTANT: this computes a candidate new binning and shows how much
    // more evenly the event rate would be distributed across bins, but it
    // does NOT rerun the fit. The posterior mean/stddev maps further below
    // are always evaluated on the ORIGINAL YAML binning, since the MCMC
    // chain's param_N branches only exist for the parameters that were
    // actually floated in that binning. To get posterior mean/stddev on
    // the new binning, regenerate templateparams.yaml using the edges
    // written out below, then rerun the fit and this script on the new
    // posterior file.
    // -------------------------------------------------------------------
    TH2D* h_xsecvar_eventratehist_quantile = nullptr;
    TH1D* h_proj_quantilevar = nullptr;
    std::vector<double> quantileEdges;

    if (quantileBins > 0) {
        bool rebinVar1 = (quantileVarChoice == "var1");
        const std::string& quantVarName = rebinVar1 ? xsec_var1 : xsec_var2;
        std::vector<double>& fixedEdges  = rebinVar1 ? q3_edges  : q0_edges;   // the OTHER axis, kept as-is
        double loEdge = rebinVar1 ? q0_edges.front() : q3_edges.front();
        double hiEdge = rebinVar1 ? q0_edges.back()  : q3_edges.back();

        std::cout << "\n[Info] Computing " << quantileBins << "-bin quantile binning for "
                  << quantVarName << " (range " << loEdge << " to " << hiEdge << ")...\n";

        TH1D* hFine = BuildFine1DHist(DUNEPdfs, quantVarName, loEdge, hiEdge, 2000);
        quantileEdges = ComputeQuantileEdges(hFine, quantileBins);
        delete hFine;

        WriteQuantileEdgesToFile(
            quantileEdges, quantVarName,
            outDir + quantVarName + "_quantile_bins.yaml.txt");

        // Rebuild the 2D event-rate histogram with the new quantile axis on
        // the chosen variable, keeping the other axis fixed at its current
        // (YAML) edges, so the comparison isolates the effect of rebinning
        // just this one variable.
        TAxis* axisQuant = new TAxis(quantileEdges.size() - 1, quantileEdges.data());
        TAxis* axisFixed = new TAxis(fixedEdges.size() - 1, fixedEdges.data());
        TAxis* axisXQ = rebinVar1 ? axisQuant : axisFixed;
        TAxis* axisYQ = rebinVar1 ? axisFixed : axisQuant;

        for (auto& pdf : DUNEPdfs) {
            std::vector<KinematicCut> SelectionVector;
            TH2* h = pdf->Get2DVarHist(0, xsec_var1, xsec_var2, SelectionVector, 0, axisXQ, axisYQ);
            if (!h) continue;
            if (!h_xsecvar_eventratehist_quantile) {
                h_xsecvar_eventratehist_quantile = (TH2D*) h->Clone("h_xsecvar_eventratehist_quantile");
                h_xsecvar_eventratehist_quantile->SetDirectory(nullptr);
            } else {
                h_xsecvar_eventratehist_quantile->Add((TH2D*)h);
            }
        }

        if (h_xsecvar_eventratehist_quantile) {
            if (h_xsecvar_eventratehist_quantile->GetSumw2N() == 0)
                h_xsecvar_eventratehist_quantile->Sumw2();

            h_xsecvar_eventratehist_quantile->GetXaxis()->SetTitle(toLatexLabel(xsec_var1).c_str());
            h_xsecvar_eventratehist_quantile->GetYaxis()->SetTitle(toLatexLabel(xsec_var2).c_str());
            h_xsecvar_eventratehist_quantile->GetZaxis()->SetTitle("Events");

            h_proj_quantilevar = rebinVar1
                ? h_xsecvar_eventratehist_quantile->ProjectionX("h_proj_quantilevar")
                : h_xsecvar_eventratehist_quantile->ProjectionY("h_proj_quantilevar");
            normalizeByBinWidth(h_proj_quantilevar);
            h_proj_quantilevar->GetXaxis()->SetTitle(toLatexLabel(quantVarName).c_str());
            h_proj_quantilevar->GetYaxis()->SetTitle("Events / bin width");
            styleThesis1D(h_proj_quantilevar);

            std::cout << "[Info] Quantile-binned event-rate histogram built ("
                      << quantileEdges.size() - 1 << " bins on " << quantVarName << ").\n";
        } else {
            std::cerr << "[Warning] Quantile rebinning requested but no histogram was produced.\n";
        }

        delete axisQuant;
        delete axisFixed;
    }

    // --- Setup posterior parameter reading ---
    // Only read the param_N branches that correspond to parameters defined
    // in templateparams.yaml (i.e. those present in binDefs).
    int nParams = 0;
    for (const auto& bin : binDefs) {
        nParams = std::max(nParams, bin.index + 1);
    }
    std::cout << "[Info] Restricting to " << binDefs.size()
              << " parameters from templateparams.yaml (param_0 .. param_"
              << nParams - 1 << ")\n";

    std::vector<double> xsec_vals(nParams, 0.0);
    post->SetBranchStatus("*", 0);
    for (const auto& bin : binDefs) {
        std::string bname = "param_" + std::to_string(bin.index);
        if (!post->GetBranch(bname.c_str())) {
            std::cerr << "[Warning] Branch " << bname
                      << " not found in tree, skipping.\n";
            continue;
        }
        post->SetBranchStatus(bname.c_str(), 1);
        post->SetBranchAddress(bname.c_str(), &xsec_vals[bin.index]);
    }

    // --- Compute posterior mean and stddev per bin ---
    TH2D* h_mean = new TH2D("h_param_mean", "",
                            q0_edges.size()-1, q0_edges.data(),
                            q3_edges.size()-1, q3_edges.data());
    TH2D* h_stddev = new TH2D("h_param_stddev", "",
                              q0_edges.size()-1, q0_edges.data(),
                              q3_edges.size()-1, q3_edges.data());

    h_mean->GetXaxis()->SetTitle(toLatexLabel(xsec_var1).c_str());
    h_mean->GetYaxis()->SetTitle(toLatexLabel(xsec_var2).c_str());
    h_mean->GetZaxis()->SetTitle("Posterior mean (relative to prior)");

    h_stddev->GetXaxis()->SetTitle(toLatexLabel(xsec_var1).c_str());
    h_stddev->GetYaxis()->SetTitle(toLatexLabel(xsec_var2).c_str());
    h_stddev->GetZaxis()->SetTitle("Posterior std. dev.");

    h_mean->GetXaxis()->SetNdivisions(505);
    h_mean->GetYaxis()->SetNdivisions(505);
    h_stddev->GetXaxis()->SetNdivisions(505);
    h_stddev->GetYaxis()->SetNdivisions(505);

    std::vector<double> sum(nParams, 0.0);
    std::vector<double> sq_sum(nParams, 0.0);

    Long64_t nEntries = post->GetEntries();

    // --- Burn-in cut ---
    Long64_t burnIn = 0;
    if (FitManager->raw()["General"]["MCMC"]["BurnInSteps"]) {
        burnIn = FitManager->raw()["General"]["MCMC"]["BurnInSteps"].as<Long64_t>();
    }
    if (burnIn <= 0 || burnIn >= nEntries) {
        std::cout << "[Info] BurnInSteps not set, zero, or exceeds chain length; "
                     "defaulting to 20% of the chain.\n";
        burnIn = static_cast<Long64_t>(0.2 * nEntries);
    }
    std::cout << "[Info] Applying burn-in cut: discarding first " << burnIn
              << " of " << nEntries << " entries ("
              << 100.0 * burnIn / nEntries << "%).\n";

    Long64_t nUsed = nEntries - burnIn;
    for (Long64_t i = burnIn; i < nEntries; ++i) {
        post->GetEntry(i);
        for (const auto& bin : binDefs) {
            if (bin.index >= nParams) continue;
            double val = xsec_vals[bin.index];
            sum[bin.index]    += val;
            sq_sum[bin.index] += val * val;
        }
    }

    for (const auto& bin : binDefs) {
        if (bin.index >= nParams) continue;
        double mean   = sum[bin.index] / nUsed;
        double var    = std::max(0.0, (sq_sum[bin.index] / nUsed) - mean * mean);
        double stddev = std::sqrt(var);

        double q0_center = 0.5 * (bin.q0_min + bin.q0_max);
        double q3_center = 0.5 * (bin.q3_min + bin.q3_max);

        int bin_q0 = h_mean->GetXaxis()->FindBin(q0_center);
        int bin_q3 = h_mean->GetYaxis()->FindBin(q3_center);

        h_mean->SetBinContent(bin_q0, bin_q3, mean);
        h_stddev->SetBinContent(bin_q0, bin_q3, stddev);
    }

    // --- Write output ---
    auto OutputFile = std::unique_ptr<TFile>(TFile::Open(rootOut.c_str(), "RECREATE"));
    if (!OutputFile || OutputFile->IsZombie()) {
        std::cerr << "Error creating output file: " << rootOut << std::endl;
        return 1;
    }
    OutputFile->cd();

    h_xsecvar_eventratehist->Write("xsec_eventrate_histo");
    h_proj_q0->Write("xsec_eventrate_proj_q0");
    h_proj_q3->Write("xsec_eventrate_proj_q3");
    h_mean->Write("xsec_param_mean");
    h_stddev->Write("xsec_param_stddev");

    if (h_xsecvar_eventratehist_quantile) {
        h_xsecvar_eventratehist_quantile->Write("xsec_eventrate_histo_quantile");
    }
    if (h_proj_quantilevar) {
        h_proj_quantilevar->Write("xsec_eventrate_proj_quantilevar");
    }

    // --- Save to PDF ---
    TCanvas* c = new TCanvas("c", "c", 1000, 750);
    c->Print((pdfOut + "[").c_str());  // open PDF

    setViridisPalette();
    h_xsecvar_eventratehist->Draw("COLZ");
    drawTitle("Event Rate");
    c->Print(pdfOut.c_str());
    c->Clear();

    // 1D projections: histogram bars with statistical error bars on top
    h_proj_q0->Draw("HIST");
    h_proj_q0->Draw("E1 SAME");
    drawTitle((toLatexLabel(xsec_var1)).c_str());
    c->Print(pdfOut.c_str());
    c->Clear();

    h_proj_q3->Draw("HIST");
    h_proj_q3->Draw("E1 SAME");
    drawTitle(("Event rate projection: " + toLatexLabel(xsec_var2)).c_str());
    c->Print(pdfOut.c_str());
    c->Clear();

    // Quantile-binning comparison pages, if computed above.
    if (h_xsecvar_eventratehist_quantile) {
        setViridisPalette();
        h_xsecvar_eventratehist_quantile->Draw("COLZ");
        drawTitle("Event Rate (quantile-rebinned)");
        c->Print(pdfOut.c_str());
        c->Clear();
    }
    if (h_proj_quantilevar) {
        // Overlay old (fixed-width, current-YAML) vs new (quantile) 1D
        // projections of the rebinned variable for a direct visual check
        // of how much more balanced the event population is per bin.
        bool rebinVar1 = (quantileVarChoice == "var1");
        TH1D* hOldProj = rebinVar1 ? h_proj_q0 : h_proj_q3;

        h_proj_quantilevar->SetLineColor(kAzure + 1);
        h_proj_quantilevar->SetMarkerColor(kAzure + 1);
        h_proj_quantilevar->Draw("HIST");
        h_proj_quantilevar->Draw("E1 SAME");

        TH1D* hOldClone = (TH1D*) hOldProj->Clone("hOldProjClone");
        hOldClone->SetLineColor(kOrange + 7);
        hOldClone->SetMarkerColor(kOrange + 7);
        hOldClone->SetLineStyle(2);
        hOldClone->Draw("HIST SAME");
        hOldClone->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.55, 0.72, 0.90, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(h_proj_quantilevar, "Quantile binning", "l");
        leg->AddEntry(hOldClone, "Current (YAML) binning", "l");
        leg->Draw();

        drawTitle(("Event rate: current vs quantile binning ("
                   + toLatexLabel(quantileVarChoice == "var1" ? xsec_var1 : xsec_var2) + ")").c_str());
        c->Print(pdfOut.c_str());
        c->Clear();
        delete hOldClone;
        delete leg;
    }

    double meanZmin = 0.0, meanZmax = 2.0, whiteCenter = 1.0;
    double whiteStop = (whiteCenter - meanZmin) / (meanZmax - meanZmin);
    setBlueWhiteRedPalette(whiteStop);
    h_mean->GetZaxis()->SetRangeUser(meanZmin, meanZmax);
    h_mean->Draw("COLZ");
    drawTitle("Posterior post-fit mean, #LT#theta#GT");
    c->Print(pdfOut.c_str());
    c->Clear();

    setViridisPalette();
    h_stddev->GetZaxis()->SetRangeUser(0.0, 0.5);
    h_stddev->Draw("COLZ");
    drawTitle("Posterior post-fit std. dev., #sigma_{#theta}");
    c->Print(pdfOut.c_str());
    c->Clear();

    c->Print((pdfOut + "]").c_str());  // close PDF

    // --- Cleanup ---
    delete h_xsecvar_eventratehist;
    delete h_proj_q0;
    delete h_proj_q3;
    delete h_mean;
    delete h_stddev;
    delete h_xsecvar_eventratehist_quantile;
    delete h_proj_quantilevar;
    delete axisX;
    delete axisY;
    delete c;
    f_post->Close();
    delete f_post;

    OutputFile->Write();
    OutputFile->Close();

    std::cout << "Saved diagnostics to: " << rootOut << std::endl;
    std::cout << "PDF written to: " << pdfOut << std::endl;
    return 0;
}