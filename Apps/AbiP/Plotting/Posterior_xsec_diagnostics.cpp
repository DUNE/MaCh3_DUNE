#include "Samples/BinningHandler.h"
#include "Samples/HistogramUtils.h"
#include "Fitters/MaCh3Factory.h"
#include "Samples/SampleHandlerBase.h"
#include "Samples/SampleHandlerFD.h"
#include "Samples/MaCh3DUNEFactory.h"
#include "Manager/MaCh3Exception.h"
#include "Manager/MaCh3Logger.h"
#include "Splines/SplineMonolith.h"
#include "Fitters/FitterBase.h"
#include "Manager/Manager.h"
#include "Parameters/ParameterHandlerBase.h"
#include "Parameters/ParameterStructs.h"

#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <map>
#include <ctime>
#include <random>
#include <fstream>
#include <sstream>
#include <set>
#include <chrono>

#include "yaml-cpp/yaml.h"

#include "TLatex.h"
#include "TH1D.h"
#include "THStack.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TRint.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TMath.h"

// ============================================================
// Data structures
// ============================================================
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

struct HistAxis {
    int bins;
    double min;
    double max;
    HistAxis(int b, double mn, double mx) : bins(b), min(mn), max(mx) {}
};

struct PosteriorSample {
    std::vector<double> xsec_draw;
};

// ============================================================
// Filenames
// ============================================================
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

// ============================================================
// Thesis-quality style (matches PosteriorMean.cpp)
// ============================================================
void setThesisStyle() {
    gStyle->Reset();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);           // titles drawn manually via TLatex

    const int font = 132;             // Times, precision 2. Use 42 for Helvetica.
    gStyle->SetTextFont(font);
    gStyle->SetLabelFont(font, "XYZ");
    gStyle->SetTitleFont(font, "XYZ");
    gStyle->SetLegendFont(font);

    gStyle->SetLabelSize(0.035, "XYZ");
    gStyle->SetTitleSize(0.045, "XYZ");
    gStyle->SetTitleOffset(1.1, "X");
    gStyle->SetTitleOffset(1.3, "Y");
    gStyle->SetTitleOffset(1.4, "Z");

    gStyle->SetPadTopMargin(0.10);
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
    gStyle->SetPaintTextFormat("4.2f");
}

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

void drawTitle(const std::string& text) {
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(132);
    latex.SetTextSize(0.045);
    latex.SetTextAlign(21); // centered
    latex.DrawLatex(0.5, 0.95, text.c_str());
}

// ============================================================
// LaTeX-style axis labels for known MaCh3-DUNE kinematic variables
// ============================================================
std::string formatAxisTitle(const std::string& var) {
    if (var == "q0")                       return "q_{0}";
    if (var == "q3")                       return "q_{3}";
    if (var == "TrueNeutrinoEnergy")       return "E_{#nu}^{True} [GeV]";
    if (var == "ENuProxy_minus_Enutrue")   return "E_{#nu}^{True} - E_{#nu}^{Proxy} [GeV]";
    if (var == "yRec")                     return "y_{Rec}";
    if (var == "RecoNeutrinoEnergy")       return "E_{#nu}^{Reco} [GeV]";
    if (var == "isRelativeEnubias")        return "E_{#nu}^{Bias} / E_{#nu} [GeV]";
    if (var == "Enubias")                  return "E_{#nu}^{Bias} [GeV]";
    // add more mappings as needed
    return var; // fallback: raw variable name
}

// ============================================================
// Correlation matrix helpers
// ============================================================
TMatrixD CovToCorr(const TMatrixD& cov) {
    int n = cov.GetNrows();
    TMatrixD corr(n, n);
    for (int i = 0; i < n; ++i) {
        double sigma_i = std::sqrt(cov(i, i));
        for (int j = 0; j < n; ++j) {
            double sigma_j = std::sqrt(cov(j, j));
            corr(i, j) = (sigma_i > 0 && sigma_j > 0) ? cov(i, j) / (sigma_i * sigma_j) : 0.0;
        }
    }
    return corr;
}

std::unique_ptr<TMatrixD> LoadCorrMatrix(TFile* f) {
    if (auto corr = dynamic_cast<TMatrixDSym*>(f->Get("Correlation"))) {
        std::cout << "[Info] Using correlation matrix, size: "
                  << corr->GetNrows() << "x" << corr->GetNcols() << "\n";
        return std::make_unique<TMatrixD>(*corr);
    } else if (auto cov = dynamic_cast<TMatrixDSym*>(f->Get("covariance_matrix"))) {
        std::cout << "[Info] Using covariance matrix, size: "
                  << cov->GetNrows() << "x" << cov->GetNcols()
                  << " -> converting to correlation.\n";
        return std::make_unique<TMatrixD>(CovToCorr(*cov));
    }
    std::cerr << "[Error] No correlation or covariance matrix found in file.\n";
    return nullptr;
}

// ============================================================
// Binning / tree helpers
// ============================================================
BinningResult extract_2D_bins_from_yaml(const std::string& yaml_file,
                                        const std::string& xsecvar1,
                                        const std::string& xsecvar2) {
    BinningResult result;
    std::set<double> q0_set, q3_set;

    try {
        YAML::Node config = YAML::LoadFile(yaml_file);
        const auto& systs = config["Systematics"];
        int index = 0;

        for (const auto& systematic : systs) {
            const auto& sys  = systematic["Systematic"];
            const auto& cuts = sys["KinematicCuts"];

            double q0_min = 0, q0_max = 0;
            double q3_min = 0, q3_max = 0;
            bool has_q0 = false, has_q3 = false;

            for (const auto& cut : cuts) {
                if (cut[xsecvar1] && cut[xsecvar1].IsSequence() && cut[xsecvar1].size() == 2) {
                    q0_min = cut[xsecvar1][0].as<double>();
                    q0_max = cut[xsecvar1][1].as<double>();
                    q0_set.insert(q0_min);
                    q0_set.insert(q0_max);
                    has_q0 = true;
                }
                if (cut[xsecvar2] && cut[xsecvar2].IsSequence() && cut[xsecvar2].size() == 2) {
                    q3_min = cut[xsecvar2][0].as<double>();
                    q3_max = cut[xsecvar2][1].as<double>();
                    q3_set.insert(q3_min);
                    q3_set.insert(q3_max);
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

// ============================================================
// Posterior predictive
// ============================================================
TH1D* MakePosteriorPredictiveHist(SampleHandlerFD* pdf,
                                  const std::vector<double>& xsec_draw,
                                  ParameterHandlerGeneric* xsec,
                                  const HistAxis& axis) {
    if (!pdf) return nullptr;

    if (xsec) xsec->SetParameters(xsec_draw);
    pdf->Reweight();

    TH1* h_pdf = pdf->GetMCHist(0);
    if (!h_pdf) {
        std::cerr << "PDF has no 1D histogram!" << std::endl;
        return nullptr;
    }

    TH1D* h_out = new TH1D(Form("posterior_predictive_%s", pdf->GetName()),
                           "", axis.bins, axis.min, axis.max);

    for (int bin = 1; bin <= axis.bins; ++bin) {
        double x = h_out->GetBinCenter(bin);
        int pdfBin = h_pdf->GetXaxis()->FindBin(x);
        h_out->SetBinContent(bin, h_pdf->GetBinContent(pdfBin));
    }

    return h_out;
}

// ============================================================
// Slicing / correlation helpers
// ============================================================
std::vector<int> GetSliceIndices(const std::vector<BinDef>& binDefs,
                                 int sliceBin,
                                 bool sliceAlong_xvar_2,
                                 TH2D* h_mean,
                                 int nCorrRows) {
    std::vector<int> sliceIndices;
    for (const auto& bin : binDefs) {
        int bin_index = sliceAlong_xvar_2
                        ? h_mean->GetYaxis()->FindBin(0.5*(bin.q3_min + bin.q3_max)) - 1
                        : h_mean->GetXaxis()->FindBin(0.5*(bin.q0_min + bin.q0_max)) - 1;

        if (bin_index == sliceBin - 1) {
            if (bin.index >= 0 && bin.index < nCorrRows) {
                sliceIndices.push_back(bin.index);
            } else {
                std::cerr << "[Warning] bin.index " << bin.index
                          << " exceeds corr matrix size, skipping.\n";
            }
        }
    }
    return sliceIndices;
}

TMatrixD GetSliceCorrelation(const TMatrixD& corr, const std::vector<int>& sliceIndices) {
    int m = sliceIndices.size();
    TMatrixD sliceCorr(m, m);
    int nrows = corr.GetNrows();
    int ncols = corr.GetNcols();

    for (int a = 0; a < m; ++a) {
        for (int b = 0; b < m; ++b) {
            if (sliceIndices[a] < nrows && sliceIndices[b] < ncols) {
                sliceCorr(a, b) = corr(sliceIndices[a], sliceIndices[b]);
            } else {
                std::cerr << "[Warning] slice index out of bounds ("
                          << sliceIndices[a] << ", " << sliceIndices[b] << ")\n";
                sliceCorr(a, b) = 0.0;
            }
        }
    }
    return sliceCorr;
}

TH2D* PlotSliceCorrelation(const TMatrixD& sliceCorr,
                           const std::vector<std::string>& paramNames,
                           const std::vector<int>& sliceIndices,
                           int sliceNumber) {
    int m = sliceIndices.size();
    TH2D* h_corr = new TH2D(Form("h_corr_slice_%d", sliceNumber), "",
                            m, 0, m, m, 0, m);

    for (int a = 0; a < m; ++a)
        for (int b = 0; b < m; ++b)
            h_corr->SetBinContent(a + 1, b + 1, sliceCorr(a, b));

    for (int a = 0; a < m; ++a) {
        std::string label = (sliceIndices[a] < (int)paramNames.size())
                                ? paramNames[sliceIndices[a]]
                                : Form("param_%d", sliceIndices[a]);
        h_corr->GetXaxis()->SetBinLabel(a + 1, label.c_str());
        h_corr->GetYaxis()->SetBinLabel(a + 1, label.c_str());
    }

    h_corr->GetZaxis()->SetRangeUser(-1, 1);
    h_corr->SetStats(0);
    return h_corr;
}

// Combined page: event rate slice + posterior slice + correlation matrix for that slice
void Save1DSlicesWithEventRateAndCorr(
    TH2D* h_event, TH2D* h_event_err,
    TH2D* h_post,  TH2D* h_post_err,
    const TMatrixD& corr,
    const std::vector<BinDef>& binDefs,
    const std::vector<std::string>& param_names,
    TCanvas* c,
    const std::string& pdfOut,
    bool sliceAlong_xvar_2,
    const std::string& xvar1,
    const std::string& xvar2,
    TFile* fOut) {

    if (!h_post || !h_event || !h_post_err || !h_event_err) return;
    fOut->cd();

    int nSlices = sliceAlong_xvar_2 ? h_post->GetNbinsY() : h_post->GetNbinsX();
    int nBinsX  = h_post->GetNbinsX();
    int nBinsY  = h_post->GetNbinsY();

    setBlueWhiteRedPalette(0.5);

    const double fillThreshold = 1e-12;
    const double corrMaskThreshold = 0.0;
    const double sentinel = -999.0; // marker for masked values

    for (int i = 1; i <= nSlices; ++i) {
        c->Clear();
        std::vector<int> sliceInds = GetSliceIndices(binDefs, i, sliceAlong_xvar_2, h_post, corr.GetNrows());
        if (sliceInds.empty()) continue;

        const BinDef& firstBin = binDefs[sliceInds.front()];
        const BinDef& lastBin  = binDefs[sliceInds.back()];

        auto formatEdge = [](double val) {
            if (fabs(val - round(val)) < 1e-6)        return Form("%.0f", val);
            else if (fabs(val*10 - round(val*10)) < 1e-6) return Form("%.1f", val);
            else                                        return Form("%.3f", val);
        };

        std::string sliceAxisLabel, orthogonalAxisLabel;
        double minEdge, maxEdge;
        if (sliceAlong_xvar_2) {
            sliceAxisLabel = xvar1;
            orthogonalAxisLabel = xvar2;
            minEdge = firstBin.q3_min;
            maxEdge = lastBin.q3_max;
        } else {
            sliceAxisLabel = xvar2;
            orthogonalAxisLabel = xvar1;
            minEdge = firstBin.q0_min;
            maxEdge = lastBin.q0_max;
        }

        std::string pageTitle = Form("Slice %d: %s #leq %s < %s (%s vs %s)",
                             i, formatEdge(minEdge), formatAxisTitle(sliceAxisLabel).c_str(),
                             formatEdge(maxEdge), formatAxisTitle(orthogonalAxisLabel).c_str(),
                             formatAxisTitle(sliceAxisLabel).c_str());

        c->Divide(2, 1);
        TPad* padLeft  = (TPad*) c->cd(1);
        TPad* padRight = (TPad*) c->cd(2);

        TH1D* h_event_slice = sliceAlong_xvar_2
            ? h_event->ProjectionX(Form("h_event_slice_%d", i), i, i)
            : h_event->ProjectionY(Form("h_event_slice_%d", i), i, i);
        TH1D* h_post_slice = sliceAlong_xvar_2
            ? new TH1D(Form("h_post_slice_%d", i), "", nBinsX, h_post->GetXaxis()->GetXbins()->GetArray())
            : new TH1D(Form("h_post_slice_%d", i), "", nBinsY, h_post->GetYaxis()->GetXbins()->GetArray());

        for (int b = 1; b <= (sliceAlong_xvar_2 ? nBinsX : nBinsY); ++b) {
            double content = sliceAlong_xvar_2 ? h_post->GetBinContent(b, i) : h_post->GetBinContent(i, b);
            double err     = sliceAlong_xvar_2 ? h_post_err->GetBinContent(b, i) : h_post_err->GetBinContent(i, b);
            h_post_slice->SetBinContent(b, content);
            h_post_slice->SetBinError(b, err);
        }

        padLeft->Divide(1, 2, 0, 0);
        int nbins = h_post_slice->GetNbinsX();

        int firstFilled = 0, lastFilled = 0;
        for (int b = 1; b <= nbins; ++b) {
            double vpost = h_post_slice->GetBinContent(b);
            double vev = (b <= h_event_slice->GetNbinsX()) ? h_event_slice->GetBinContent(b) : 0.0;
            if ((std::isfinite(vpost) && vpost > fillThreshold) ||
                (std::isfinite(vev) && vev > fillThreshold)) { firstFilled = b; break; }
        }
        for (int b = nbins; b >= 1; --b) {
            double vpost = h_post_slice->GetBinContent(b);
            double vev = (b <= h_event_slice->GetNbinsX()) ? h_event_slice->GetBinContent(b) : 0.0;
            if ((std::isfinite(vpost) && vpost > fillThreshold) ||
                (std::isfinite(vev) && vev > fillThreshold)) { lastFilled = b; break; }
        }
        if (firstFilled > 0 && lastFilled > 0) {
            double xlow  = h_post_slice->GetXaxis()->GetBinLowEdge(firstFilled);
            double xhigh = h_post_slice->GetXaxis()->GetBinUpEdge(lastFilled);
            h_event_slice->GetXaxis()->SetRangeUser(xlow, xhigh);
            h_post_slice->GetXaxis()->SetRangeUser(xlow, xhigh);
        }
        h_post_slice->GetYaxis()->SetRangeUser(0.2, 1.8);

        // xvar1 is always the axis remaining after projection (ProjectionX when
        // sliceAlong_xvar_2, ProjectionY otherwise) - see call sites in main().
        std::string projectedAxisLabel = formatAxisTitle(xvar1);
        padLeft->cd(1); gPad->SetBottomMargin(0.12);
        h_event_slice->GetXaxis()->SetTitle(projectedAxisLabel.c_str());
        h_event_slice->GetYaxis()->SetTitle("Events");
        h_event_slice->Draw("E1 HIST");

        padLeft->cd(2); gPad->SetBottomMargin(0.12);
        h_post_slice->GetXaxis()->SetTitle(projectedAxisLabel.c_str());
        h_post_slice->GetYaxis()->SetTitle("Posterior value");
        h_post_slice->Draw("E1 HIST");

        h_event_slice->Write();
        h_post_slice->Write();

        c->cd();
        drawTitle(pageTitle);
        c->Update();

        // --- Correlation slice ---
        TMatrixD sliceCorr = GetSliceCorrelation(corr, sliceInds);
        for (int r = 0; r < sliceCorr.GetNrows(); ++r)
            for (int ccol = 0; ccol < sliceCorr.GetNcols(); ++ccol) {
                double v = sliceCorr(r, ccol);
                if (!std::isfinite(v) || fabs(v) < corrMaskThreshold) sliceCorr(r, ccol) = sentinel;
            }

        if (sliceCorr.GetNrows() > 0) {
            TH2D* h_corr_slice = PlotSliceCorrelation(sliceCorr, param_names, sliceInds, i);
            if (h_corr_slice) {
                std::string var1 = xvar1.empty() ? "q0" : formatAxisTitle(xvar1);
                std::string var2 = xvar2.empty() ? "q3" : formatAxisTitle(xvar2);

                auto fmt = [](double vmin, double vmax) {
                    TString s;
                    if (fabs(vmin - round(vmin)) < 1e-6 && fabs(vmax - round(vmax)) < 1e-6)
                        s.Form("%.0f-%.0f", vmin, vmax);
                    else
                        s.Form("%.2f-%.2f", vmin, vmax);
                    return s;
                };

                for (int idx = 0; idx < (int)sliceInds.size(); ++idx) {
                    const BinDef& bin = binDefs[sliceInds[idx]];
                    TString label = Form("(%s, %s) = (%s, %s)",
                                        var1.c_str(), var2.c_str(),
                                        fmt(bin.q0_min, bin.q0_max).Data(),
                                        fmt(bin.q3_min, bin.q3_max).Data());
                    h_corr_slice->GetXaxis()->SetBinLabel(idx + 1, label);
                    h_corr_slice->GetYaxis()->SetBinLabel(idx + 1, label);
                }

                int m = h_corr_slice->GetNbinsX();
                double labsz = std::max(0.008, 0.025 - 0.0006 * m);
                h_corr_slice->SetLabelSize(labsz, "X");
                h_corr_slice->SetLabelSize(labsz, "Y");
                h_corr_slice->GetXaxis()->LabelsOption("v");
                h_corr_slice->GetYaxis()->LabelsOption("h");

                padRight->cd();
                padRight->SetRightMargin(0.20);
                padRight->SetLeftMargin(0.25);
                padRight->SetBottomMargin(0.30);
                h_corr_slice->SetMinimum(-1.0);
                h_corr_slice->SetMaximum(1.0);

                for (int ix = 1; ix <= h_corr_slice->GetNbinsX(); ++ix)
                    for (int iy = 1; iy <= h_corr_slice->GetNbinsY(); ++iy)
                        if (h_corr_slice->GetBinContent(ix, iy) <= sentinel/2)
                            h_corr_slice->SetBinContent(ix, iy, -1.1);

                h_corr_slice->Draw("COLZ");
                gPad->Modified();
                gPad->Update();
                h_corr_slice->Write();
                c->Modified();
                c->Update();
                c->Print(pdfOut.c_str());
                delete h_corr_slice;
            }
        }

        delete h_post_slice;
        delete h_event_slice;
    }
}

// ============================================================
// main
// ============================================================
int main(int argc, char* argv[]) {

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <mcmc_output.root> <config.yaml> <corr_matrix.root>" << std::endl;
        return 1;
    }

    setThesisStyle();

    std::string mcmc_file   = argv[1];
    std::string config_file = argv[2];
    std::string corr_file   = argv[3];

    std::string rootOut = AddTimestampToROOTFilename("mcmc_diagnostics.root");
    std::string pdfOut  = mcmc_file.substr(0, mcmc_file.find_last_of('.')) + "_diagnostics.pdf";

    TCanvas* c_master = new TCanvas("c", "c", 1200, 800);

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

    // --- Load correlation/covariance matrix and parameter names ---
    TFile* fCorr = TFile::Open(corr_file.c_str(), "READ");
    if (!fCorr || fCorr->IsZombie()) {
        std::cerr << "Error: cannot open " << corr_file << "\n";
        return 1;
    }
    auto corrMatrix = LoadCorrMatrix(fCorr);
    if (!corrMatrix) return 1;

    std::vector<std::string> paramNames;
    TObject* namesObj = fCorr->Get("branch_names");
    if (auto arr = dynamic_cast<TObjArray*>(namesObj); arr && arr->GetEntries() > 0) {
        for (int i = 0; i < arr->GetEntries(); ++i) {
            if (auto s = dynamic_cast<TObjString*>(arr->At(i))) {
                paramNames.push_back(s->GetString().Data());
            }
        }
        std::cout << "[Info] Loaded " << paramNames.size()
                  << " paramNames from " << corr_file << std::endl;
    }
    if (paramNames.empty()) {
        std::cerr << "[Warning] paramNames empty. Creating placeholder names.\n";
        paramNames.resize(corrMatrix->GetNrows());
        for (int i = 0; i < (int)paramNames.size(); ++i) paramNames[i] = Form("param_%d", i);
    }

    // --- Binning from YAML ---
    auto xsec_var1 = FitManager->raw()["General"]["Systematics"]["xsec_var1"].as<std::string>();
    auto xsec_var2 = FitManager->raw()["General"]["Systematics"]["xsec_var2"].as<std::string>();
    auto xsec_yaml = FitManager->raw()["General"]["Systematics"]["XsecCovFile"][0].as<std::string>();

    if (xsec_yaml.empty()) {
        std::cerr << "Error: xsec_yaml path is empty or missing in config.\n";
        return 1;
    }

    BinningResult binning = extract_2D_bins_from_yaml(xsec_yaml, xsec_var1, xsec_var2);
    std::vector<BinDef>&  binDefs  = binning.binDefs;
    std::vector<double>&  q0_edges = binning.q0_edges;
    std::vector<double>&  q3_edges = binning.q3_edges;

    if (q0_edges.empty() || q3_edges.empty()) {
        std::cerr << "Bin edges could not be parsed from YAML file." << std::endl;
        return 1;
    }

    TAxis* axisX = new TAxis(q0_edges.size() - 1, q0_edges.data());
    TAxis* axisY = new TAxis(q3_edges.size() - 1, q3_edges.data());

    TH2D* h_mean = new TH2D("h_param_mean", "",
                            q0_edges.size()-1, q0_edges.data(),
                            q3_edges.size()-1, q3_edges.data());
    TH2D* h_stddev = new TH2D("h_param_stddev", "",
                              q0_edges.size()-1, q0_edges.data(),
                              q3_edges.size()-1, q3_edges.data());
    h_mean->GetXaxis()->SetTitle(formatAxisTitle(xsec_var1).c_str());
    h_mean->GetYaxis()->SetTitle(formatAxisTitle(xsec_var2).c_str());
    h_mean->GetZaxis()->SetTitle("Posterior mean (relative to prior)");
    h_stddev->GetXaxis()->SetTitle(formatAxisTitle(xsec_var1).c_str());
    h_stddev->GetYaxis()->SetTitle(formatAxisTitle(xsec_var2).c_str());
    h_stddev->GetZaxis()->SetTitle("Posterior std. dev.");

    // --- Output ROOT file ---
    auto OutputFile = std::unique_ptr<TFile>(TFile::Open(rootOut.c_str(), "RECREATE"));
    if (!OutputFile || OutputFile->IsZombie()) {
        std::cerr << "Error creating output file " << rootOut << std::endl;
        return 1;
    }
    OutputFile->cd();

    // --- MaCh3 setup ---
    ParameterHandlerGeneric* xsec = nullptr;
    std::vector<SampleHandlerFD*> DUNEPdfs;
    MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);
    // Note: MaCh3FitterFactory is not currently used beyond construction here.
    // Keeping the call in case downstream code relies on side effects of building it;
    // remove if not needed.
    auto MaCh3Fitter = MaCh3FitterFactory(FitManager.get());

    // --- Build combined 2D event-rate histogram ---
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
    h_xsecvar_eventratehist->GetXaxis()->SetTitle(formatAxisTitle(xsec_var1).c_str());
    h_xsecvar_eventratehist->GetYaxis()->SetTitle(formatAxisTitle(xsec_var2).c_str());
    h_xsecvar_eventratehist->GetZaxis()->SetTitle("Events");
    h_xsecvar_eventratehist->Write("xsec_eventrate_histo");

    // --- Labeled clone: bin index overlaid on each cell ---
    TH2D* h_eventrate_labels = (TH2D*) h_xsecvar_eventratehist->Clone("h_eventrate_labels");
    h_eventrate_labels->Write();

    int nBinsX = h_eventrate_labels->GetNbinsX();
    int nBinsY = h_eventrate_labels->GetNbinsY();
    double textSize = std::min(0.05, 0.5 * std::min(1.0/nBinsX, 1.0/nBinsY));

    TLatex binLabels;
    binLabels.SetTextAlign(22);
    binLabels.SetTextSize(textSize);

    double maxContent = h_eventrate_labels->GetMaximum();
    for (const auto& bin : binDefs) {
        double q0_center = 0.5 * (bin.q0_min + bin.q0_max);
        double q3_center = 0.5 * (bin.q3_min + bin.q3_max);
        int binX = h_eventrate_labels->GetXaxis()->FindBin(q0_center);
        int binY = h_eventrate_labels->GetYaxis()->FindBin(q3_center);
        double binContent = h_eventrate_labels->GetBinContent(binX, binY);

        int color = (maxContent > 0 && binContent / maxContent > 0.5) ? kWhite : kBlack;
        binLabels.SetTextColor(color);
        binLabels.DrawLatex(q0_center, q3_center, std::to_string(bin.index).c_str());
    }

    // --- Determine which param_N branches to read ---
    int nParams = 0;
    for (const auto& bin : binDefs) nParams = std::max(nParams, bin.index + 1);
    std::cout << "[Info] Reading " << binDefs.size()
              << " parameters from templateparams.yaml (param_0 .. param_"
              << nParams - 1 << ")\n";

    std::vector<double> xsec_vals(nParams, 0.0);
    post->SetBranchStatus("*", 0);
    for (int i = 0; i < nParams; ++i) {
        std::string bname = "param_" + std::to_string(i);
        if (!post->GetBranch(bname.c_str())) {
            std::cerr << "[Warning] Branch " << bname << " not found in tree, skipping.\n";
            continue;
        }
        post->SetBranchStatus(bname.c_str(), 1);
        post->SetBranchAddress(bname.c_str(), &xsec_vals[i]);
    }

    // --- Burn-in cut ---
    // Applied consistently to both the direct posterior mean/stddev below
    // and the posterior-predictive event-rate sampling further down.
    Long64_t nEntries = post->GetEntries();
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

    // --- Posterior mean & stddev per bin (post burn-in) ---
    std::vector<double> sum(nParams, 0.0);
    std::vector<double> sq_sum(nParams, 0.0);

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
    h_mean->Write("xsec_param_mean");
    h_stddev->Write("xsec_param_stddev");

    // --- Posterior sampling for posterior-predictive event rate ---
    int nSamples = 1000;
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<Long64_t> dis(burnIn, nEntries - 1);

    std::vector<PosteriorSample> posteriorSamples;
    posteriorSamples.reserve(nSamples);
    for (int i = 0; i < nSamples; ++i) {
        Long64_t idx = dis(gen);
        post->GetEntry(idx);
        PosteriorSample sample;
        sample.xsec_draw.assign(xsec_vals.begin(), xsec_vals.end());
        posteriorSamples.push_back(std::move(sample));
    }
    std::cout << "Loaded " << posteriorSamples.size()
              << " random posterior samples from " << nEntries << " entries." << std::endl;

    HistAxis fineAxis(1000, 0, 10); // GeV

    TH2D* h_eventrate_sum    = (TH2D*) h_xsecvar_eventratehist->Clone("h_eventrate_sum");
    TH2D* h_eventrate_sum2   = (TH2D*) h_xsecvar_eventratehist->Clone("h_eventrate_sum2");
    TH2D* h_eventrate_mean   = (TH2D*) h_xsecvar_eventratehist->Clone("h_eventrate_mean");
    TH2D* h_eventrate_stddev = (TH2D*) h_xsecvar_eventratehist->Clone("h_eventrate_stddev");
    h_eventrate_sum->Reset();
    h_eventrate_sum2->Reset();
    h_eventrate_mean->Reset();
    h_eventrate_stddev->Reset();

    for (const auto& sample : posteriorSamples) {
        for (auto* pdf : DUNEPdfs) {
            if (!pdf) continue;
            TH1D* h_pred = MakePosteriorPredictiveHist(pdf, sample.xsec_draw, xsec, fineAxis);
            if (!h_pred) continue;

            for (int bx = 1; bx <= h_eventrate_sum->GetNbinsX(); ++bx) {
                for (int by = 1; by <= h_eventrate_sum->GetNbinsY(); ++by) {
                    // Simple approximation: map 2D bin to nearest predictive bin
                    int predBin = std::min(std::max(1, bx), h_pred->GetNbinsX());
                    double val = h_pred->GetBinContent(predBin);

                    double evSum  = h_eventrate_sum->GetBinContent(bx, by);
                    double evSum2 = h_eventrate_sum2->GetBinContent(bx, by);
                    h_eventrate_sum->SetBinContent(bx, by, evSum + val);
                    h_eventrate_sum2->SetBinContent(bx, by, evSum2 + val * val);
                }
            }
            delete h_pred;
        }
    }

    for (int bx = 1; bx <= h_eventrate_sum->GetNbinsX(); ++bx) {
        for (int by = 1; by <= h_eventrate_sum->GetNbinsY(); ++by) {
            double evSum  = h_eventrate_sum->GetBinContent(bx, by);
            double evSum2 = h_eventrate_sum2->GetBinContent(bx, by);
            double mean   = evSum / nSamples;
            double var    = evSum2 / nSamples - mean * mean;
            double stddev = (var > 0) ? std::sqrt(var) : 0.0;

            h_eventrate_mean->SetBinContent(bx, by, mean);
            h_eventrate_stddev->SetBinContent(bx, by, stddev);
        }
    }

    // --- h_stddev with the bottom row removed (often a low-stats overflow-like bin) ---
    const TArrayD* xEdges = h_stddev->GetXaxis()->GetXbins();
    const TArrayD* yEdges = h_stddev->GetYaxis()->GetXbins();
    if (xEdges->GetSize() == 0 || yEdges->GetSize() == 0) {
        std::cerr << "[Error] Histogram does not use variable binning.\n";
    }

    std::vector<double> newYEdges;
    for (int i = 1; i < yEdges->GetSize(); ++i) newYEdges.push_back(yEdges->At(i));

    TH2D* h_stddev_noBottom = new TH2D("h_param_stddev_noBottom", "",
        h_stddev->GetNbinsX(), xEdges->GetArray(),
        h_stddev->GetNbinsY() - 1, newYEdges.data());
    h_stddev_noBottom->SetDirectory(nullptr);
    h_stddev_noBottom->GetXaxis()->SetTitle(formatAxisTitle(xsec_var1).c_str());
    h_stddev_noBottom->GetYaxis()->SetTitle(formatAxisTitle(xsec_var2).c_str());
    h_stddev_noBottom->GetZaxis()->SetTitle("Posterior std. dev.");

    for (int bx = 1; bx <= h_stddev->GetNbinsX(); ++bx) {
        for (int by = 2; by <= h_stddev->GetNbinsY(); ++by) { // skip bottom row
            h_stddev_noBottom->SetBinContent(bx, by - 1, h_stddev->GetBinContent(bx, by));
            h_stddev_noBottom->SetBinError(bx, by - 1, h_stddev->GetBinError(bx, by));
        }
    }

    // ========================================================
    // PDF pages
    // ========================================================
    c_master->Print((pdfOut + "[").c_str()); // open PDF

    // Page 1: event rate with bin-index labels
    c_master->Clear();
    setViridisPalette();
    h_eventrate_labels->GetXaxis()->SetTitle(formatAxisTitle(xsec_var1).c_str());
    h_eventrate_labels->GetYaxis()->SetTitle(formatAxisTitle(xsec_var2).c_str());
    h_eventrate_labels->GetZaxis()->SetTitle("Events");
    h_eventrate_labels->Draw("COLZ TEXT");
    drawTitle("DUNE-PRISM event rate (bin index labelled)");
    c_master->Print(pdfOut.c_str());

    // Page 2: event rate
    c_master->Clear();
    h_xsecvar_eventratehist->Draw("COLZ");
    drawTitle("DUNE-PRISM event rate");
    c_master->Update();
    c_master->Print(pdfOut.c_str());

    // Page 3: posterior mean, white centred at 1.0
    c_master->Clear();
    {
        double meanZmin = 0.5, meanZmax = 1.5, whiteCenter = 1.0;
        double whiteStop = (whiteCenter - meanZmin) / (meanZmax - meanZmin);
        setBlueWhiteRedPalette(whiteStop);
        h_mean->GetZaxis()->SetRangeUser(meanZmin, meanZmax);
    }
    h_mean->Draw("COLZ");
    drawTitle("Posterior post-fit mean");
    c_master->Update();
    c_master->Print(pdfOut.c_str());

    // Page 4: posterior stddev
    c_master->Clear();
    setViridisPalette();
    h_stddev->GetZaxis()->SetRangeUser(0.0, 0.5);
    h_stddev->Draw("COLZ");
    drawTitle("Posterior post-fit std. dev.");
    c_master->Update();
    c_master->Print(pdfOut.c_str());

    // Page 5: posterior stddev with bottom row removed
    c_master->Clear();
    h_stddev_noBottom->GetZaxis()->SetRangeUser(0.0, 0.5);
    h_stddev_noBottom->Draw("COLZ");
    drawTitle("Posterior post-fit std. dev. (bottom row removed)");
    h_stddev_noBottom->Write("xsec_param_stddev_noBottom");
    c_master->Update();
    c_master->Print(pdfOut.c_str());

    std::cout << "xsec_var1 = " << xsec_var1 << std::endl;
    std::cout << "xsec_var2 = " << xsec_var2 << std::endl;

    // Combined slice pages, both slicing directions
    Save1DSlicesWithEventRateAndCorr(
        h_eventrate_labels, h_eventrate_stddev,
        h_mean, h_stddev,
        *corrMatrix, binDefs, paramNames,
        c_master, pdfOut,
        /*sliceAlong_xvar_2=*/true,
        xsec_var1, xsec_var2,
        OutputFile.get());

    Save1DSlicesWithEventRateAndCorr(
        h_eventrate_labels, h_eventrate_stddev,
        h_mean, h_stddev,
        *corrMatrix, binDefs, paramNames,
        c_master, pdfOut,
        /*sliceAlong_xvar_2=*/false,
        xsec_var2, xsec_var1,
        OutputFile.get());

    c_master->Print((pdfOut + "]").c_str()); // close PDF

    // --- Cleanup ---
    delete h_xsecvar_eventratehist;
    delete h_eventrate_labels;
    delete h_eventrate_sum;
    delete h_eventrate_sum2;
    delete h_eventrate_mean;
    delete h_eventrate_stddev;
    delete h_mean;
    delete h_stddev;
    delete h_stddev_noBottom;
    delete axisX;
    delete axisY;
    delete c_master;
    f_post->Close();
    delete f_post;
    fCorr->Close();
    delete fCorr;

    OutputFile->Write();
    OutputFile->Close();

    std::cout << "Saved diagnostics to: " << rootOut << std::endl;
    std::cout << "Finished diagnostics. Output written to: " << pdfOut << std::endl;
    return 0;
}
