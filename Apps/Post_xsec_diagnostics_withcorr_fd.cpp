#include <iostream>
#include <chrono>
#include "TLatex.h"
#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include "TPaveText.h"
#include <TColor.h>
#include <TMath.h>
#include <cmath> 
#include <fstream>
#include <sstream>
#include <set>
#include <string>
#include <yaml-cpp/yaml.h>
#include <vector>
#include <ctime>
#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"

#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"
#include "mcmc/mcmc.h"
#include "manager/manager.h"
#include "covariance/covarianceBase.h"
#include <vector>
#include <string>
#include <TH1D.h>
#include <TString.h>
#include <iostream>
#include <random>

#include <iostream>
#include <chrono>
#include "TLatex.h"
#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include "TPaveText.h"
#include <TColor.h>
#include <TMath.h>
#include <cmath> 
#include <fstream>
#include <sstream>
#include <set>
#include <string>
#include <yaml-cpp/yaml.h>
#include <vector>
#include <ctime>
#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"

#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"
#include "mcmc/mcmc.h"
#include "manager/manager.h"
#include "covariance/covarianceBase.h"
#include <vector>
#include <string>
#include <TH1D.h>
#include <TString.h>
#include <iostream>
#include <random>



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

std::string title;

// Forward declaration of your PDF type
class samplePDFFDBase;
struct PosteriorSample {
    std::vector<double> xsec_draw; // vector of xsec parameters for this posterior sample
};

TMatrixD CovToCorr(const TMatrixD& cov) {
    int n = cov.GetNrows();
    TMatrixD corr(n, n);

    for (int i = 0; i < n; ++i) {
        double sigma_i = std::sqrt(cov(i, i));
        for (int j = 0; j < n; ++j) {
            double sigma_j = std::sqrt(cov(j, j));
            if (sigma_i > 0 && sigma_j > 0) {
                corr(i, j) = cov(i, j) / (sigma_i * sigma_j);
            } else {
                corr(i, j) = 0.0;
            }
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
                  << " → converting to correlation.\n";
        TMatrixD corrMat = CovToCorr(*cov);
        return std::make_unique<TMatrixD>(corrMat);
    } else {
        std::cerr << "[Error] No correlation or covariance matrix found in file.\n";
        return nullptr;
    }
}


// Example histogram axis struct
struct HistAxis {
    int bins;
    double min;
    double max;
    HistAxis(int b, double mn, double mx) : bins(b), min(mn), max(mx) {}
};


std::string AddTimestampToFilename(const std::string& baseName) {
    time_t now = time(0);
    tm* ltm = localtime(&now);
    std::ostringstream oss;
    oss << baseName.substr(0, baseName.find_last_of(".")) << "_"
        << 1900 + ltm->tm_year
        << 1 + ltm->tm_mon
        << ltm->tm_mday << "_"
        << ltm->tm_hour
        << ltm->tm_min
        << ltm->tm_sec << ".pdf";
    return oss.str();
}

BinningResult extract_2D_bins_from_yaml(const std::string& yaml_file, const std::string& xsecvar1, const std::string& xsecvar2) {
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
                if (cut[xsecvar1]) {
                    const auto& q0 = cut[xsecvar1];
                    if (q0.IsSequence() && q0.size() == 2) {
                        q0_min = q0[0].as<double>();
                        q0_max = q0[1].as<double>();
                        has_q0 = true;
                    }
                }
                if (cut[xsecvar2]) {
                    const auto& q3 = cut[xsecvar2];
                    if (q3.IsSequence() && q3.size() == 2) {
                        q3_min = q3[0].as<double>();
                        q3_max = q3[1].as<double>();
                        has_q3 = true;
                    }
                }
            }

            if (has_q0 && has_q3) {
                q0_set.insert(q0_min);
                q0_set.insert(q0_max);
                q3_set.insert(q3_min);
                q3_set.insert(q3_max);

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

  std::string formatAxisTitle(const std::string& var) {
    if (var == "q0") return "q_{0}";
    if (var == "q3") return "q_{3}";
    if (var == "TrueNeutrinoEnergy") return "E^_{#nu} [GeV]";
    if (var == "ENuProxy_minus_Enutrue") return "E^{True}_{#nu} -E^{Proxy}_{#nu} [GeV]";
    if (var == "yRec") return "y_{Rec}";
    if (var == "RecoNeutrinoEnergy") return "E^{Reco.}_{#nu} [GeV]";
    if (var == "isRelativeEnubias") return "E^{Bias}_{#nu}/E^_{#nu}  [GeV]";
    // Add more mappings as needed

    // Fallback: just return var unchanged if no match found
    return var;
}




void SavePosterior1D(const std::vector<PosteriorSample>& posteriorSamples,
                     int paramIndex,
                     const std::string& name,
                     TCanvas* c,
                     const std::string& pdfOut) 
{
    TH1D* h_post = new TH1D(Form("posterior_%s_%d", name.c_str(), paramIndex),
                            Form("Posterior %s (param %d)", name.c_str(), paramIndex),
                            50, 0, 10); // adjust bins/range

    for (const auto& sample : posteriorSamples)
        h_post->Fill(sample.xsec_draw[paramIndex]);

    h_post->SetLineColor(kBlack);
    h_post->SetLineWidth(2);
    h_post->Draw("HIST");
    h_post->Write();
    
    c->Print(pdfOut.c_str());
    delete h_post;
}


std::pair<TH1D*, TH1D*> GetSliceWithErrors(TH2D* h_mean,
                                           TH2D* h_stddev,
                                           int slice,
                                           bool sliceAlongY)
{
    TH1D* h_slice_mean = sliceAlongY
        ? h_mean->ProjectionX(Form("h_mean_slice_%d", slice), slice, slice)
        : h_mean->ProjectionY(Form("h_mean_slice_%d", slice), slice, slice);

    TH1D* h_slice_std = sliceAlongY
        ? h_stddev->ProjectionX(Form("h_std_slice_%d", slice), slice, slice)
        : h_stddev->ProjectionY(Form("h_std_slice_%d", slice), slice, slice);

    // Attach bin errors
    for (int b = 1; b <= h_slice_mean->GetNbinsX(); ++b)
        h_slice_mean->SetBinError(b, h_slice_std->GetBinContent(b));

    return {h_slice_mean, h_slice_std};
}


void Save1DSlicesWithErrors(TH2D* h_mean, TH2D* h_error,
                            const std::string& sliceLabel,
                            TCanvas* c,
                            const std::string& pdfOut,
                            bool sliceAlongQ3,
                            const std::string& xLabel)
{
    int nSlices = sliceAlongQ3 ? h_mean->GetNbinsY() : h_mean->GetNbinsX();

    for (int i = 1; i <= nSlices; ++i) {
        c->Clear();

        // Project mean and error
        TH1D* h_slice_mean  = sliceAlongQ3 ?
            h_mean->ProjectionX(Form("%s_mean_slice_%d", h_mean->GetName(), i), i, i) :
            h_mean->ProjectionY(Form("%s_mean_slice_%d", h_mean->GetName(), i), i, i);

        TH1D* h_slice_error = sliceAlongQ3 ?
            h_error->ProjectionX(Form("%s_error_slice_%d", h_error->GetName(), i), i, i) :
            h_error->ProjectionY(Form("%s_error_slice_%d", h_error->GetName(), i), i, i);

        // Attach error to mean slice
        for (int b = 1; b <= h_slice_mean->GetNbinsX(); ++b)
            h_slice_mean->SetBinError(b, h_slice_error->GetBinContent(b));

        // Titles
        double fixedLow  = sliceAlongQ3 ? h_mean->GetYaxis()->GetBinLowEdge(i) : h_mean->GetXaxis()->GetBinLowEdge(i);
        double fixedHigh = sliceAlongQ3 ? h_mean->GetYaxis()->GetBinUpEdge(i)  : h_mean->GetXaxis()->GetBinUpEdge(i);
        h_slice_mean->SetTitle(Form("%s: %.2f < %s < %.2f", sliceLabel.c_str(), fixedLow, xLabel.c_str(), fixedHigh));

        // Axis labels
        h_slice_mean->GetXaxis()->SetTitle(xLabel.c_str());
        h_slice_mean->GetYaxis()->SetTitle("Value");

        // Style
        h_slice_mean->SetLineColor(kBlack);
        h_slice_mean->SetLineWidth(2);
        h_slice_mean->SetMarkerSize(0);
        h_slice_mean->SetFillColorAlpha(kAzure+1, 0.35); // optional: shaded band

        // Draw shaded error and mean
        h_slice_mean->Draw("E1");        // non shaded ±error
        h_slice_mean->Draw("HIST SAME"); // mean line on top
        h_slice_mean->Write();
        h_slice_error->Write();
        c->Update();

        // Save one page
        c->Print(pdfOut.c_str());

        delete h_slice_mean;
        delete h_slice_error;
    }
}

    // -------------------------
    // Helper: get latest non-empty posteriors tree
    // -------------------------
    TTree* getLatestNonEmpty(TFile* f) {
        int cycle = 9999;  // start from a high number
        while (cycle > 0) {
            TString name = Form("posteriors;%d", cycle);
            TTree* t = (TTree*) f->Get(name);
            if (t && t->GetEntries() > 0) return t;
            cycle--;
        }
        return nullptr; // not found
    }


void Save1DSlicePair(TH2D* h_event_mean, TH2D* h_event_stddev,
                     TH2D* h_param_mean, TH2D* h_param_stddev,
                     const std::string& sliceLabel,
                     TCanvas* c,
                     const std::string& pdfOut,
                     bool sliceAlongQ3,
                     const std::string& xLabel)
{
    int nSlices = sliceAlongQ3 ? h_event_mean->GetNbinsY() : h_event_mean->GetNbinsX();

    for (int i = 1; i <= nSlices; ++i) {
        c->Clear();
        c->Divide(1,2); // stacked pads: top = event rate, bottom = posterior

        // --- Slice 1: Event rate ---
        TH1D* h_ev_mean = sliceAlongQ3 ?
            h_event_mean->ProjectionX(Form("ev_mean_slice_%d",i), i,i) :
            h_event_mean->ProjectionY(Form("ev_mean_slice_%d",i), i,i);

        TH1D* h_ev_std = sliceAlongQ3 ?
            h_event_stddev->ProjectionX(Form("ev_std_slice_%d",i), i,i) :
            h_event_stddev->ProjectionY(Form("ev_std_slice_%d",i), i,i);

        for (int b = 1; b <= h_ev_mean->GetNbinsX(); ++b)
            h_ev_mean->SetBinError(b, h_ev_std->GetBinContent(b));

        h_ev_mean->SetLineColor(kBlue+2);
        h_ev_mean->SetLineWidth(2);
        h_ev_mean->SetMarkerStyle(20);
        h_ev_mean->SetMarkerColor(kBlue+2);
        h_ev_mean->GetXaxis()->SetTitle(xLabel.c_str());
        h_ev_mean->GetYaxis()->SetTitle("Event rate");
        h_ev_mean->SetTitle(Form("%s Slice %d: Event rate", sliceLabel.c_str(), i));

        // auto y-axis range
        h_ev_mean->SetMinimum(0);
        h_ev_mean->SetMaximum(1.2*h_ev_mean->GetBinContent(h_ev_mean->GetMaximumBin()));

        c->cd(1);
        h_ev_mean->Draw("E1");
        h_ev_mean->Write();

        // --- Slice 2: Posterior values ---
        TH1D* h_pa_mean = sliceAlongQ3 ?
            h_param_mean->ProjectionX(Form("pa_mean_slice_%d",i), i,i) :
            h_param_mean->ProjectionY(Form("pa_mean_slice_%d",i), i,i);

        TH1D* h_pa_std = sliceAlongQ3 ?
            h_param_stddev->ProjectionX(Form("pa_std_slice_%d",i), i,i) :
            h_param_stddev->ProjectionY(Form("pa_std_slice_%d",i), i,i);

        for (int b = 1; b <= h_pa_mean->GetNbinsX(); ++b)
            h_pa_mean->SetBinError(b, h_pa_std->GetBinContent(b));

        h_pa_mean->SetLineColor(kRed+1);
        h_pa_mean->SetLineWidth(2);
        h_pa_mean->SetMarkerStyle(21);
        h_pa_mean->SetMarkerColor(kRed+1);
        h_pa_mean->GetXaxis()->SetTitle(xLabel.c_str());
        h_pa_mean->GetYaxis()->SetTitle("Posterior value");
        h_pa_mean->SetTitle(Form("%s Slice %d: Posterior", sliceLabel.c_str(), i));

        // auto y-axis range
        h_pa_mean->SetMinimum(0);
        h_pa_mean->SetMaximum(1.2*h_pa_mean->GetBinContent(h_pa_mean->GetMaximumBin()));

        c->cd(2);
        h_pa_mean->Draw("E1");
        h_pa_mean->Write();

        // export page
        c->Print(pdfOut.c_str());

        // cleanup
        delete h_ev_mean; delete h_ev_std;
        delete h_pa_mean; delete h_pa_std;
    }
}


TH1D* MakePosteriorPredictiveHist(
    samplePDFFDBase* pdf,                     // your PDF object
    const std::vector<double>& xsec_draw,  
    covarianceXsec* xsec, // pointer to your cross-section object         // cross-section parameters for this posterior sample
    const std::vector<std::string>& xsec_vars,// names of xsec parameters
    const HistAxis& axis                       // histogram axis info
) {
    if (!pdf) return nullptr;

    // Apply xsec parameters to the reweighter
    if (xsec) xsec->setParameters(xsec_draw);
    // Update PDF with new weights
    pdf->reweight();
    // Create a new histogram for this posterior predictive
    TString histName = Form("posterior_predictive_%s", pdf->GetTitle());
    // Get the PDF histogram
    TH1D* h_pdf = pdf->get1DHist();
    if (!h_pdf) {
        std::cerr << "PDF has no 1D histogram!" << std::endl;
        return nullptr;
    }

    TH1D* h_out = new TH1D(Form("posterior_predictive_%s", pdf->GetTitle()),
                           Form("posterior_predictive_%s", pdf->GetTitle()),
                           axis.bins, axis.min, axis.max);


    // Loop over bins and fill the posterior predictive histogram
    for (int bin = 1; bin <= axis.bins; ++bin) {
        double x = h_out->GetBinCenter(bin);
        int pdfBin = h_pdf->GetXaxis()->FindBin(x);
        double val = h_pdf->GetBinContent(pdfBin);
        h_out->SetBinContent(bin, val);
    }

    return h_out;
}


std::vector<int> GetSliceIndices(const std::vector<BinDef>& binDefs,
                                 int sliceBin,
                                 bool sliceAlongQ3,
                                 TH2D* h_mean,
                                 int nCorrRows) {  // pass corr.GetNrows()
    std::vector<int> sliceIndices;
    for (const auto& bin : binDefs) {
        int bin_index = sliceAlongQ3
                        ? h_mean->GetYaxis()->FindBin(0.5*(bin.q3_min + bin.q3_max)) - 1
                        : h_mean->GetXaxis()->FindBin(0.5*(bin.q0_min + bin.q0_max)) - 1;

        if (bin_index == sliceBin - 1) {
            std::cout << "[Debug] Slice " << sliceBin 
                      << " matched bin.index " << bin.index
                      << " (bin_index=" << bin_index << ")\n";

            if (bin.index >= 0 && bin.index < nCorrRows) {
                sliceIndices.push_back(bin.index);
            } else {
                std::cerr << "[Warning] bin.index " << bin.index
                          << " exceeds corr matrix size, skipping.\n";
            }
        }
    }

    std::cout << "[Debug] GetSliceIndices found " << sliceIndices.size()
          << " indices for slice " << sliceBin << std::endl;

    return sliceIndices;

}


TH2D* PlotSliceCorrelation(const TMatrixD& sliceCorr,
                           const std::vector<std::string>& paramNames,
                           const std::vector<int>& sliceIndices,
                           int sliceNumber)
{
    int m = sliceIndices.size();
    TH2D* h_corr = new TH2D(Form("h_corr_slice_%d", sliceNumber),
                            Form(" %d", sliceNumber),
                            m, 0, m, m, 0, m);

    for (int a = 0; a < m; ++a) {
        for (int b = 0; b < m; ++b) {
            h_corr->SetBinContent(a + 1, b + 1, sliceCorr(a, b));
        }
    }

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
                std::cerr << "Warning: slice index out of bounds ("
                          << sliceIndices[a] << ", " << sliceIndices[b] << ")\n";
                sliceCorr(a, b) = 0.0;
            }
        }
    }
    return sliceCorr;
}




TH2D* PlotSliceCorrelation_threshold(const TMatrixD& sliceCorr,
                                    const std::vector<std::string>& paramNames,
                                    const std::vector<int>& sliceIndices,
                                    int sliceNumber,
                                    double threshold_correlation) {
    int m = static_cast<int>(sliceIndices.size());

    // Safety checks
    if (m <= 0) {
        std::cerr << "[PlotSliceCorrelation_threshold] Warning: empty sliceIndices (m=0). Skipping.\n";
        return nullptr;
    }
    if (sliceCorr.GetNrows() != m || sliceCorr.GetNcols() != m) {
        std::cerr << "[PlotSliceCorrelation_threshold] Error: sliceCorr dim mismatch: "
                  << sliceCorr.GetNrows() << "x" << sliceCorr.GetNcols()
                  << " vs m=" << m << "\n";
        return nullptr;
    }

    // verify indices are within bounds of paramNames and sliceCorr source indices are valid
    int maxParamIndex = static_cast<int>(paramNames.size()) - 1;
    int maxSliceIndex = -1;
    for (int idx : sliceIndices) maxSliceIndex = std::max(maxSliceIndex, idx);
    if (maxParamIndex < 0) {
        std::cerr << "[PlotSliceCorrelation_threshold] Error: paramNames empty.\n";
        return nullptr;
    }
    if (maxSliceIndex < 0) {
        std::cerr << "[PlotSliceCorrelation_threshold] Error: sliceIndices contains negative index.\n";
        return nullptr;
    }
    if (maxSliceIndex >= (int)paramNames.size()) {
        std::cerr << "[PlotSliceCorrelation_threshold] Error: sliceIndices refer to paramNames index "
                  << maxSliceIndex << " but paramNames.size()=" << paramNames.size() << "\n";
        return nullptr;
    }

    // Create histogram
    TH2D* h_corr = new TH2D(Form("h_corr_slice_%d", sliceNumber),
                            Form("", sliceNumber),
                            m, 0, m, m, 0, m);

    // Fill contents with threshold applied
    for (int a = 0; a < m; ++a) {
        for (int b = 0; b < m; ++b) {
            if (a==0 && b==0) {
                //std::cerr << "[Debug] First entry in slice " << sliceNumber << " = " << value << std::endl;
            }

            double value = sliceCorr(a, b);
            if (!std::isfinite(value)) value = 0.0;
                h_corr->SetBinContent(a+1, b+1, value);

            // Optional: blank label for small correlations
            if (std::fabs(value) < threshold_correlation) {
                h_corr->SetBinContent(a+1, b+1, value); // still draw the color
                // but don’t add text if you later use "TEXT" option
            }
        }
    }

    // Set labels safely: use temporary std::string copies to ensure lifetime
    for (int a = 0; a < m; ++a) {
        int pi = sliceIndices[a];
        if (pi >= 0 && pi < (int)paramNames.size()) {
            // SetBinLabel stores label internally, but to be safe convert to TString
            TString label(paramNames[pi]);
            h_corr->GetXaxis()->SetBinLabel(a+1, label.Data());
            h_corr->GetYaxis()->SetBinLabel(a+1, label.Data());
        } else {
            h_corr->GetXaxis()->SetBinLabel(a+1, Form("%d", pi));
            h_corr->GetYaxis()->SetBinLabel(a+1, Form("%d", pi));
            std::cerr << "[PlotSliceCorrelation_threshold] Warning: label index " << pi
                      << " out-of-range, using numeric label.\n";
        }
    }

    h_corr->GetZaxis()->SetRangeUser(-1, 1);
    h_corr->SetStats(0);

    return h_corr;
}


void Save1DSlicesWithEventRateAndCorr(
    TH2D* h_event, TH2D* h_event_err,
    TH2D* h_post,  TH2D* h_post_err,
    const TMatrixD& corr,
    const std::vector<BinDef>& binDefs,
    const std::vector<std::string>& param_names,
    TCanvas* c,
    const std::string& pdfOut,
    bool sliceAlongQ3,
    const std::string& xvar1,
    const std::string& xvar2,
    TFile* fOut)
{
    if (!h_post || !h_event || !h_post_err || !h_event_err) return;

    fOut->cd();

    int nSlices = sliceAlongQ3 ? h_post->GetNbinsY() : h_post->GetNbinsX();
    int nBinsX  = h_post->GetNbinsX();
    int nBinsY  = h_post->GetNbinsY();

    // --- Red-white-blue palette ---
    const Int_t nRGBs = 3;
    Double_t stops[nRGBs] = {0.0, 0.5, 1.0};
    Double_t red[nRGBs]   = {0.0, 1.0, 1.0};
    Double_t green[nRGBs] = {0.0, 1.0, 0.0};
    Double_t blue[nRGBs]  = {1.0, 1.0, 0.0};
    TColor::CreateGradientColorTable(nRGBs, stops, red, green, blue, 255);
    gStyle->SetNumberContours(255);
    gStyle->SetPaintTextFormat("4.2f");

    const int contextBins = 3;
    const double fillThreshold = 1e-12;
    const double corrMaskThreshold = 0.1;
    const double sentinel = -999.0; // marker for masked values

    for (int i = 1; i <= nSlices; ++i) {
        c->Clear();
        std::vector<int> sliceInds = GetSliceIndices(binDefs, i, sliceAlongQ3, h_post, corr.GetNrows());
        if (sliceInds.empty()) continue;

        const BinDef& firstBin = binDefs[sliceInds.front()];
        const BinDef& lastBin  = binDefs[sliceInds.back()];

        auto formatEdge = [](double val) {
            if (fabs(val - round(val)) < 1e-6) return Form("%.0f", val);
            else if (fabs(val*10 - round(val*10)) < 1e-6) return Form("%.1f", val);
            else return Form("%.3f", val);
        };

        std::string pageTitle = sliceAlongQ3
            ? Form("Slice %d: %s #leq q3 < %s GeV", i, formatEdge(firstBin.q3_min), formatEdge(lastBin.q3_max))
            : Form("Slice %d: %s #leq q0 < %s GeV", i, formatEdge(firstBin.q0_min), formatEdge(lastBin.q0_max));

        // Left/right pads
        c->Divide(2,1);
        TPad* padLeft  = (TPad*) c->cd(1);
        TPad* padRight = (TPad*) c->cd(2);

        // --- Event slices ---
        TH1D* h_event_slice = sliceAlongQ3
            ? h_event->ProjectionX(Form("h_event_slice_%d", i), i, i)
            : h_event->ProjectionY(Form("h_event_slice_%d", i), i, i);
        TH1D* h_post_slice = sliceAlongQ3
            ? new TH1D(Form("h_post_slice_%d", i), "", nBinsX, h_post->GetXaxis()->GetXbins()->GetArray())
            : new TH1D(Form("h_post_slice_%d", i), "", nBinsY, h_post->GetYaxis()->GetXbins()->GetArray());

        // Fill posterior slice & errors
        for (int b = 1; b <= (sliceAlongQ3 ? nBinsX : nBinsY); ++b) {
            double content = sliceAlongQ3 ? h_post->GetBinContent(b,i) : h_post->GetBinContent(i,b);
            double err     = sliceAlongQ3 ? h_post_err->GetBinContent(b,i) : h_post_err->GetBinContent(i,b);
            h_post_slice->SetBinContent(b, content);
            h_post_slice->SetBinError(b, err);
        }
        h_event_slice->SetTitle(pageTitle.c_str());
        // Left pads draw
        padLeft->Divide(1,2,0,0);
        //int nbins = h_post_slice->GetNbinsX();
        int nbins = h_post_slice->GetNbinsX();
        h_post_slice->SetTitle(pageTitle.c_str());

        int firstFilled = 0;
        int lastFilled = 0;

        // find first non-empty bin
        for (int b = 1; b <= nbins; ++b) {
            double vpost = h_post_slice->GetBinContent(b);
            double vev = (b <= h_event_slice->GetNbinsX()) ? h_event_slice->GetBinContent(b) : 0.0;
            if ((std::isfinite(vpost) && vpost > fillThreshold) || 
                (std::isfinite(vev) && vev > fillThreshold)) {
                firstFilled = b;
                break;
            }
        }

        // find last non-empty bin
        for (int b = nbins; b >= 1; --b) {
            double vpost = h_post_slice->GetBinContent(b);
            double vev = (b <= h_event_slice->GetNbinsX()) ? h_event_slice->GetBinContent(b) : 0.0;
            if ((std::isfinite(vpost) && vpost > fillThreshold) || 
                (std::isfinite(vev) && vev > fillThreshold)) {
                lastFilled = b;
                break;
            }
        }

        if (firstFilled > 0 && lastFilled > 0) {
            double xlow  = h_post_slice->GetXaxis()->GetBinLowEdge(firstFilled);
            double xhigh = h_post_slice->GetXaxis()->GetBinUpEdge(lastFilled);

            h_event_slice->GetXaxis()->SetRangeUser(xlow, xhigh);
            h_post_slice->GetXaxis()->SetRangeUser(xlow, xhigh);
        }

        
        h_post_slice->GetYaxis()->SetRangeUser(0.9, 1.1);

        padLeft->cd(1); gPad->SetBottomMargin(0.12); h_event_slice->GetXaxis()->SetTitle(xvar1.c_str()); h_event_slice->Draw("E1 HIST");
        padLeft->cd(2); gPad->SetBottomMargin(0.12); h_post_slice->GetXaxis()->SetTitle(xvar1.c_str()); h_post_slice->SetYTitle("Posterior value"); h_post_slice->Draw("E1 HIST");

        h_event_slice->Write();
        h_post_slice->Write();

// --- draw the page title BEFORE printing the canvas ---
        c->cd(); // go to canvas (not a child pad)
        // ensure there is room at the top for the title
        c->SetTopMargin(0.08); // enlarge top margin if necessary
        TLatex latex; 
        latex.SetNDC(kTRUE); 
        latex.SetTextAlign(22); 
        latex.SetTextSize(0.045); 
        latex.DrawLatex(0.5, 0.96, pageTitle.c_str()); // centered title near top
        c->Update(); 

        // --- Correlation slice ---
        TMatrixD sliceCorr = GetSliceCorrelation(corr, sliceInds);
        int nonMaskedCount = 0;
        for (int r = 0; r < sliceCorr.GetNrows(); ++r) {
            for (int ccol = 0; ccol < sliceCorr.GetNcols(); ++ccol) {
                double v = sliceCorr(r, ccol);
                if (!std::isfinite(v) || fabs(v) < corrMaskThreshold) sliceCorr(r, ccol) = sentinel;
                else ++nonMaskedCount;
            }
        }

        if (sliceCorr.GetNrows() > 0) {
            TH2D* h_corr_slice = PlotSliceCorrelation(sliceCorr, param_names, sliceInds, i);
            if (h_corr_slice) {
                // Compact axis labels
                for (int idx = 0; idx < (int)sliceInds.size(); ++idx) {
                    const BinDef& bin = binDefs[sliceInds[idx]];
                    TString label = Form("[%.2f-%.2f, %.2f-%.2f]",
                     bin.q0_min, bin.q0_max,
                     bin.q3_min, bin.q3_max);
                    h_corr_slice->GetXaxis()->SetBinLabel(idx+1, label);
                    h_corr_slice->GetYaxis()->SetBinLabel(idx+1, label);
                }

                int m = h_corr_slice->GetNbinsX();
                double labsz = std::max(0.006, 0.02 - 0.0005 * m);
                h_corr_slice->SetLabelSize(labsz, "X");
                h_corr_slice->SetLabelSize(labsz, "Y");
                h_corr_slice->GetXaxis()->LabelsOption("v");

                padRight->cd();
                padRight->SetRightMargin(0.20); padRight->SetLeftMargin(0.25); padRight->SetBottomMargin(0.30);
                h_corr_slice->SetMinimum(-1.0); h_corr_slice->SetMaximum(1.0);

                // NaNs/masked bins white
                for (int ix = 1; ix <= h_corr_slice->GetNbinsX(); ++ix)
                    for (int iy = 1; iy <= h_corr_slice->GetNbinsY(); ++iy)
                        if (h_corr_slice->GetBinContent(ix, iy) <= sentinel/2) h_corr_slice->SetBinContent(ix, iy, -1.1);

                h_corr_slice->Draw("COLZ");
                gPad->Modified(); gPad->Update();
                h_corr_slice->Write();
                c->Modified(); c->Update(); c->Print(pdfOut.c_str());
                delete h_corr_slice;
            }
        }

        
        delete h_post_slice;
        delete h_event_slice;
    }
}



int main(int argc, char* argv[]) {

    // --- Setup ROOT canvas ---
    std::string pdfOut = AddTimestampToFilename("mcmc_summary_plots.pdf");
    std::string rootOut = AddTimestampToROOTFilename("mcmc_diagnostics.root");
    
    // --- Setup ROOT canvas ---
    TCanvas* c_master = new TCanvas("c", "c", 1200, 800);

    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <mcmc_output.root> <config.yaml> <corr_matrix.root>" << std::endl;
        return 1;
    }

    std::string mcmc_file  = argv[1];
    std::string config_file = argv[2];
    std::string corr_file  = argv[3];

    // -------------------------
    // Open original MCMC file
    // -------------------------
    TFile* f_post = TFile::Open(mcmc_file.c_str(), "READ");
    if (!f_post || f_post->IsZombie()) {
        std::cerr << "Error: cannot open MCMC file " << mcmc_file << std::endl;
        return 1;
    }

    TTree* post = getLatestNonEmpty(f_post); // user-defined
    if (!post) {
        std::cerr << "Error: No non-empty 'posteriors' tree found in file." << std::endl;
        return 1;
    }

    // -------------------------
    // Load YAML config
    // -------------------------
    YAML::Node config = YAML::LoadFile(config_file);
    manager* FitManager = new manager(config_file);
    if (!FitManager) {
        std::cerr << "Error: Failed to create FitManager from YAML config." << std::endl;
        return 1;
    }

    std::cout << "Using " << post->GetName() 
              << " with " << post->GetEntries() << " entries." << std::endl;

    // -------------------------
    // Open correlation file
    // -------------------------
    TFile* fCorr = TFile::Open(corr_file.c_str(), "READ");
    if (!fCorr || fCorr->IsZombie()) {
        std::cerr << "Error: cannot open " << corr_file << "\n";
        return 1;
    }

    auto corrMatrix = LoadCorrMatrix(fCorr); // user-defined
    if (!corrMatrix) return 1;

    // --- Load parameter names ---
    std::vector<std::string> paramNames;
    TDirectory* namesDir = (TDirectory*) fCorr->Get("Names");
    if (namesDir) {
        TIter nextkey(namesDir->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*) nextkey())) {
            TObject* obj = key->ReadObj();
            if (obj->InheritsFrom("TObjString")) {
                paramNames.push_back(((TObjString*)obj)->GetString().Data());
            }
        }
    }

    TObject* obj = fCorr->Get("branch_names");
    TObjArray* arr = dynamic_cast<TObjArray*>(obj);
    if (arr && arr->GetEntries() > 0) {
        paramNames.clear();
        for (int i = 0; i < arr->GetEntries(); ++i) {
            TObjString* s = dynamic_cast<TObjString*>(arr->At(i));
            if (s) paramNames.push_back(std::string(s->GetString()));
        }
        std::cerr << "[Info] Loaded " << paramNames.size() << " paramNames from " << corr_file << std::endl;
    }

    if (paramNames.empty()) {
        std::cerr << "[Warning] paramNames empty. Creating placeholder names.\n";
        paramNames.resize(corrMatrix->GetNrows());
        for (int i = 0; i < (int)paramNames.size(); ++i) paramNames[i] = Form("param_%d", i);
    }

    // -------------------------
    // Setup binning from YAML
    // -------------------------
    auto xsec_var1 = FitManager->raw()["General"]["Systematics"]["xsec_var1"].as<std::string>();
    auto xsec_var2 = FitManager->raw()["General"]["Systematics"]["xsec_var2"].as<std::string>();
    auto xsec_yaml = FitManager->raw()["General"]["Systematics"]["XsecCovFile"][0].as<std::string>();

    if (xsec_yaml.empty()) {
        std::cerr << "Error: xsec_yaml path is empty or missing in config.\n";
        return 1;
    }

    BinningResult binning = extract_2D_bins_from_yaml(xsec_yaml, xsec_var1, xsec_var2);
    std::vector<BinDef> binDefs   = binning.binDefs;
    std::vector<double> q0_edges  = binning.q0_edges;
    std::vector<double> q3_edges  = binning.q3_edges;

    if (q0_edges.empty() || q3_edges.empty()) {
        std::cerr << "Error: Bin edges could not be parsed from YAML file." << std::endl;
        return 1;
    }

    TAxis* axisX = new TAxis(q0_edges.size() - 1, q0_edges.data());
    TAxis* axisY = new TAxis(q3_edges.size() - 1, q3_edges.data());

    // -------------------------
    // Posterior mean/stddev histograms
    // -------------------------
    TH2D* h_mean = new TH2D("h_param_mean", "Posterior Mean;q0 [GeV];q3 [GeV]",
                            q0_edges.size()-1, &q0_edges[0],
                            q3_edges.size()-1, &q3_edges[0]);
    TH2D* h_stddev = new TH2D("h_param_stddev", "Posterior StdDev;",
                              q0_edges.size()-1, &q0_edges[0],
                              q3_edges.size()-1, &q3_edges[0]);

    // -------------------------
    // Slice example
    // -------------------------
    int sliceBin = 1;          // pick first slice bin
    bool sliceAlongQ3 = true;  // slice along q3 direction

    auto sliceIndices = GetSliceIndices(binDefs, sliceBin, sliceAlongQ3,
                                            h_mean, corrMatrix->GetNrows());
    if (!sliceIndices.empty()) {
        auto sliceCorr = GetSliceCorrelation(*corrMatrix, sliceIndices);
        auto hCorrPlot = PlotSliceCorrelation(sliceCorr, paramNames, sliceIndices, sliceBin);
        hCorrPlot->Draw("COLZ ");
        if (!hCorrPlot) {
        std::cerr << "[Error] PlotSliceCorrelation returned nullptr!" << std::endl;
        }
        else {
         std::cerr << "[Debug] hCorrPlot has " << hCorrPlot->GetNbinsX() << "x" << hCorrPlot->GetNbinsY() << " bins" << std::endl;
        }
    }

    


    // -------------------------
    // Derive output filenames
    // -------------------------
    std::string outputRootFile = mcmc_file.substr(0, mcmc_file.find_last_of('.')) + "_diagnostics.root";
    std::string outputPdfFile  = mcmc_file.substr(0, mcmc_file.find_last_of('.')) + "_diagnostics.pdf";

    auto OutputFile = std::unique_ptr<TFile>(TFile::Open(outputRootFile.c_str(), "RECREATE"));
    if (!OutputFile || OutputFile->IsZombie()) {
        std::cerr << "Error creating output file " << outputRootFile << std::endl;
        return 1;
    }
    OutputFile->cd();

 
    // -------------------------
    // MaCh3 setup
    // -------------------------
    covarianceXsec* xsec = nullptr;
    //covarianceOsc* osc = nullptr;

    auto OscCovFile =GetFromManager<std::vector<std::string>>(FitManager->raw()["General"]["Systematics"]["OscCovFile"], {});
    auto OscCovName =GetFromManager<std::string>(FitManager->raw()["General"]["Systematics"]["OscCovName"], "osc_cov");
    auto OscPars = GetFromManager<std::vector<double>>(FitManager->raw()["General"]["OscillationParameters"], {});

    covarianceOsc* osc = new covarianceOsc(OscCovFile, OscCovName);
    osc->setParameters(OscPars);

    std::vector<samplePDFFDBase*> DUNEPdfs;
    MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec, osc);
    std::unique_ptr<mcmc> MaCh3Fitter = std::make_unique<mcmc>(FitManager);

    // -------------------------
    // Build combined 2D histogram
    // -------------------------
    TH2D* h_q0q3_combined = nullptr;
    for (auto& pdf : DUNEPdfs) {
        auto* dunePdf = dynamic_cast<samplePDFDUNEBeamFD*>(pdf);
        if (!dunePdf) continue;

        std::vector<KinematicCut> SelectionVector;
        TH2* h = dunePdf->get2DVarHist(xsec_var1, xsec_var2, SelectionVector, 0, axisX, axisY);
        if (!h) {
            std::cerr << "Warning: get2DVarHist returned null.\n";
            continue;
        }
        if (!h_q0q3_combined) {
            h_q0q3_combined = (TH2D*)h->Clone("h_q0q3_combined");
            h_q0q3_combined->SetDirectory(nullptr);
        } else {
            h_q0q3_combined->Add(h);
        }
    }

    if (!h_q0q3_combined) {
        std::cerr << "Error: No 2D histograms were produced.\n";
        return 1;
    }

    // Write combined histogram to output diagnostics file
    h_q0q3_combined->GetXaxis()->SetTitle(xsec_var1.c_str());
    h_q0q3_combined->GetYaxis()->SetTitle(xsec_var2.c_str());
    h_q0q3_combined->Write("xsec_eventrate_histo");

    // --- Event rate histogram with adaptive label colors ---
    TH2D* h_eventrate_labels = (TH2D*)h_q0q3_combined->Clone("h_eventrate_labels");
    h_eventrate_labels->SetTitle("");

    //TCanvas* c_eventrate_labels = new TCanvas("c_eventrate_labels", "", 900, 700);
    //c_eventrate_labels->SetRightMargin(0.15);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kPastel);
    h_eventrate_labels->GetXaxis()->SetTitle(xsec_var1.c_str());
    h_eventrate_labels->GetYaxis()->SetTitle(xsec_var2.c_str());
    
    h_eventrate_labels->Draw("COLZ");
    h_eventrate_labels->Write();

    // Dynamic text size
    int nBinsX = h_eventrate_labels->GetNbinsX();
    int nBinsY = h_eventrate_labels->GetNbinsY();
    double textSize = 0.5 * std::min(1.0/nBinsX, 1.0/nBinsY);
    if (textSize > 0.05) textSize = 0.05;

    TLatex latex;
    latex.SetTextAlign(22); // center
    latex.SetTextSize(textSize);

    for (const auto& bin : binDefs) {
        double q0_center = 0.5 * (bin.q0_min + bin.q0_max);
        double q3_center = 0.5 * (bin.q3_min + bin.q3_max);

        int binX = h_eventrate_labels->GetXaxis()->FindBin(q0_center);
        int binY = h_eventrate_labels->GetYaxis()->FindBin(q3_center);
        double binContent = h_eventrate_labels->GetBinContent(binX, binY);

        // Adapt label color depending on bin content
        int color = kBlack;
        double maxContent = h_eventrate_labels->GetMaximum();
        if (binContent / maxContent > 0.5) color = kWhite;

        latex.SetTextColor(color);
        latex.DrawLatex(q0_center, q3_center, std::to_string(bin.index).c_str());
    }



    // -------------------------
    // Setup xsec branch reading
    // -------------------------
    // --- Setup xsec branch reading (revised) ---
    int nParams = 0;
    TObjArray* branches = post->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); ++i) {
        std::string name = branches->At(i)->GetName();
        if (name.rfind("xsec_", 0) == 0) {
            int idx = std::stoi(name.substr(5));
            nParams = std::max(nParams, idx + 1);
        }
    }

    std::vector<double> xsec_vals(nParams, 0.0);
    std::vector<double*> xsec_ptrs(nParams);
    for (int i = 0; i < nParams; ++i) xsec_ptrs[i] = &xsec_vals[i];

    post->SetBranchStatus("*", 0);
    // Enable ALL xsec_* for predictive sampling correctness
    for (int i = 0; i < nParams; ++i) {
        std::string bname = "xsec_" + std::to_string(i);
        post->SetBranchStatus(bname.c_str(), 1);
        post->SetBranchAddress(bname.c_str(), xsec_ptrs[i]);
    }



    std::vector<double> sum(nParams, 0.0);
    std::vector<double> sq_sum(nParams, 0.0);

    Long64_t nEntries = post->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        post->GetEntry(i);
        for (const auto& bin : binDefs) {
            if (bin.index >= nParams) continue;
            double val = xsec_vals[bin.index];
            sum[bin.index] += val;
            sq_sum[bin.index] += val * val;
        }
    }

    for (const auto& bin : binDefs) {
        if (bin.index >= nParams) continue;
        double mean = sum[bin.index] / nEntries;
        double var = std::max(0.0, (sq_sum[bin.index] / nEntries) - mean * mean);
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
    h_mean->GetZaxis()->SetRangeUser(0.95,1.05);
    h_stddev->GetZaxis()->SetRangeUser(0, 0.1);
    

    // -------------------------
    // Random subset of posterior samples
    // -------------------------
    int nSamples = 1000; // how many random samples
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<Long64_t> dis(0, nEntries-1);
    int burnIn = 0.1 * nEntries;

    std::vector<PosteriorSample> posteriorSamples;
    posteriorSamples.reserve(nSamples);

    for (int i = burnIn; i < nSamples; ++i) {
        Long64_t idx = dis(gen);
        post->GetEntry(idx);

        PosteriorSample sample;
        sample.xsec_draw.resize(nParams);
        for (int p = 0; p < nParams; ++p)
            sample.xsec_draw[p] = xsec_vals[p];

        posteriorSamples.push_back(sample);
    }

    std::cout << "Loaded " << posteriorSamples.size() << " random posterior samples from "  << nEntries << " entries." << std::endl;

    // Define your histogram axis
    HistAxis fineAxis(1000, 0, 10); // 0–5 GeV, 500 bins

    // Example: vector of cross-section variables
    std::vector<std::string> xsec_vars = {xsec_var1, xsec_var2};
    // --- Initialize histograms for event rates ---
    TH2D* h_eventrate_sum    = (TH2D*) h_q0q3_combined->Clone("h_eventrate_sum");
    TH2D* h_eventrate_sum2   = (TH2D*) h_q0q3_combined->Clone("h_eventrate_sum2");
    TH2D* h_eventrate_mean   = (TH2D*) h_q0q3_combined->Clone("h_eventrate_mean");
    TH2D* h_eventrate_stddev = (TH2D*) h_q0q3_combined->Clone("h_eventrate_stddev");

    h_eventrate_sum->Reset();
    h_eventrate_sum2->Reset();
    h_eventrate_mean->Reset();
    h_eventrate_stddev->Reset();

    // --- Loop over posterior samples ---
    for (size_t i = 0; i < posteriorSamples.size(); ++i) {
        for (size_t j = 0; j < DUNEPdfs.size(); ++j) {
            auto* dunePdf = dynamic_cast<samplePDFFDBase*>(DUNEPdfs[j]);
            if (!dunePdf) continue;

            // Posterior predictive histogram (1D)
            TH1D* h_pred = MakePosteriorPredictiveHist(dunePdf, 
                                                    posteriorSamples[i].xsec_draw,
                                                    xsec, 
                                                    xsec_vars, 
                                                    fineAxis);
            if (!h_pred) continue;

            // Accumulate into 2D event-rate histograms
            for (int bx = 1; bx <= h_eventrate_sum->GetNbinsX(); ++bx) {
                for (int by = 1; by <= h_eventrate_sum->GetNbinsY(); ++by) {
                    // Map 2D bin to 1D prediction bin (simple approximation: take central bin)
                    int predBin = std::min(std::max(1, bx), h_pred->GetNbinsX());
                    double val = h_pred->GetBinContent(predBin);

                    double sum   = h_eventrate_sum->GetBinContent(bx, by);
                    double sum2  = h_eventrate_sum2->GetBinContent(bx, by);
                    h_eventrate_sum->SetBinContent(bx, by, sum + val);
                    h_eventrate_sum2->SetBinContent(bx, by, sum2 + val*val);
                }
            }
            delete h_pred; // clean up
        }
    }

    // --- Compute mean & stddev ---
    //int nSamples = posteriorSamples.size();
    for (int bx = 1; bx <= h_eventrate_sum->GetNbinsX(); ++bx) {
        for (int by = 1; by <= h_eventrate_sum->GetNbinsY(); ++by) {
            double sum   = h_eventrate_sum->GetBinContent(bx, by);
            double sum2  = h_eventrate_sum2->GetBinContent(bx, by);

            double mean  = sum / nSamples;
            double var   = sum2/nSamples - mean*mean;
            double stddev = (var > 0) ? sqrt(var) : 0.0;

            h_eventrate_mean->SetBinContent(bx, by, mean);
            h_eventrate_stddev->SetBinContent(bx, by, stddev);
        }
    }

        // --- Open multipage PDF ---
    //c_master->Print((outputFileName + "[").c_str()); 
    c_master->Print((outputPdfFile + "[").c_str());  // open PDF


    // 1️⃣ Event rate with labels
    c_master->Clear();
    h_eventrate_labels->Draw("COLZ TEXT");
    h_eventrate_labels->Write();
    c_master->Print(outputPdfFile.c_str());          // add a page

    h_eventrate_labels->Write(); // save to ROOT

    // 2️⃣ Event rate heatmap without labels
    c_master->Clear();
    h_q0q3_combined->GetXaxis()->SetTitle(xsec_var1.c_str());
    h_q0q3_combined->GetYaxis()->SetTitle(xsec_var2.c_str());

    h_q0q3_combined->Draw("COLZ");
    h_q0q3_combined->Write();
    c_master->Update();
    //c_master->Print(outputFileName.c_str());
    c_master->Print(outputPdfFile.c_str());          // add a page


    // 3️⃣ Posterior mean 2D
    c_master->Clear();
    
    h_mean->GetXaxis()->SetTitle(xsec_var1.c_str());
    h_mean->GetYaxis()->SetTitle(xsec_var2.c_str());
    h_mean->SetTitle("Posterior post-fit mean");
    h_mean->Draw("COLZ");
    h_mean->Write();
    c_master->Update();
    //c_master->Print(outputFileName.c_str());
    c_master->Print(outputPdfFile.c_str());          // add a page


    // 4️⃣ Posterior stddev 2D
    c_master->Clear();
    h_stddev->GetXaxis()->SetTitle(xsec_var1.c_str());
    h_stddev->GetYaxis()->SetTitle(xsec_var2.c_str());

    h_stddev->SetTitle("Posterior post-fit. std dev.");

    h_stddev->Draw("COLZ");
    h_stddev->Write();
    c_master->Update();
    //c_master->Print(outputFileName.c_str());
    c_master->Print(outputPdfFile.c_str());          // add a page

    std::cout << "xsec_var1  = " << xsec_var1 << std::endl;
    std::cout << "xsec_var2  = " << xsec_var2 << std::endl;

    // Combined slices: event rate, posterior, and correlation on each page
    Save1DSlicesWithEventRateAndCorr(
        h_eventrate_labels, h_eventrate_stddev,  // event rate and its error
        h_mean, h_stddev,                        // posterior mean and stddev
        *corrMatrix,                                   // correlation matrix
        binDefs,                                 // bin definitions
        paramNames,                             // parameter names
        c_master,                                // canvas
        outputPdfFile,                                  // output PDF
        true,                                    // slice along q3
        xsec_var2,                               // x-axis variable for 1D slices
        xsec_var1 ,
        OutputFile.get()                               // y-axis variable for correlation/2D
    );


    // Combined slices: event rate, posterior, and correlation on each page
    Save1DSlicesWithEventRateAndCorr(
        h_eventrate_labels, h_eventrate_stddev,  // event rate and its error
        h_mean, h_stddev,                        // posterior mean and stddev
        *corrMatrix,                                   // correlation matrix
        binDefs,                                 // bin definitions
        paramNames,                             // parameter names
        c_master,                                // canvas
        outputPdfFile,                                  // output PDF
        false,                                    // slice along q3
        xsec_var1,                               // x-axis variable for 1D slices
        xsec_var2  ,
        OutputFile.get()                              // y-axis variable for correlation/2D
    );

    c_master->Print((outputPdfFile + "]").c_str());  // close PDF

    // -------------------------
    // Clean up
    // -------------------------
    delete h_q0q3_combined;
    delete h_mean;
    delete h_stddev;
    delete axisX;
    delete axisY;
    delete FitManager;
    delete osc;
    delete f_post;
    fCorr->Close();
    delete fCorr;

    OutputFile->Write();
    OutputFile->Close();
    std::cout << "Saved diagnostics to: " << rootOut << std::endl;
    std::cout << "Finished diagnostics. Output written to: " << outputPdfFile << std::endl;
    return 0;
}