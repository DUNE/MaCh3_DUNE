#include <iostream>
#include <chrono>
#include "TLatex.h"
#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
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

std::string title;

// Forward declaration of your PDF type
class samplePDFFDBase;
struct PosteriorSample {
    std::vector<double> xsec_draw; // vector of xsec parameters for this posterior sample
};

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
    if (var == "TrueNeutrinoEnergy") return "E^{True}_{#nu} [GeV]";
    if (var == "ENuProxy_minus_Enutrue") return "E^{True}_{#nu} -E^{Proxy}_{#nu} [GeV]";
    if (var == "yRec") return "y_{Rec}";
    if (var == "RecoNeutrinoEnergy") return "E^{Reco.}_{#nu} [GeV]";
    // Add more mappings as needed

    // Fallback: just return var unchanged if no match found
    return var;
}


////////////////////////////////////////////////////////////////////
/*
void SavePosteriorWithError(samplePDFFDBase* pdf,
                            const std::vector<int>& mcmc_draw_entries,
                            TTree* posterior_tree,
                            covarianceXsec* xsec,
                            const std::vector<double>& xsec_draw,
                            const std::string& varName,
                            const std::string& plotTitle,
                            const std::string& xLabel,
                            const std::string& pdfOutFile)
{
    // Number of posterior samples
    int nSamples = mcmc_draw_entries.size();
    const std::vector<KinematicCut> emptySelectionVec;

    // Get template histogram to define binning
    TH1* hTemplateTH1 = pdf->get1DVarHist(varName, emptySelectionVec);
    TH1D* hTemplate = dynamic_cast<TH1D*>(hTemplateTH1);
    if (!hTemplate) {
        std::cerr << "[ERROR] get1DVarHist did not return TH1D for variable " << varName << "\n";
        return;
    }

    int nBins = hTemplate->GetNbinsX();
    double xMin = hTemplate->GetXaxis()->GetXmin();
    double xMax = hTemplate->GetXaxis()->GetXmax();

    // Fill TH2: X = variable bins, Y = posterior sample index
    TH2D* h2 = new TH2D("h2Var", "Posterior samples", nBins, xMin, xMax, nSamples, 0.5, nSamples + 0.5);

    for (int i = 0; i < nSamples; ++i) {
        int entry = mcmc_draw_entries[i];
        posterior_tree->GetEntry(entry);

        // Apply current xsec parameters
        xsec->setParameters(xsec_draw);
        pdf->reweight();

        TH1* hSampleTH1 = pdf->get1DVarHist(varName, emptySelectionVec);
        TH1D* hSample = dynamic_cast<TH1D*>(hSampleTH1);
        if (!hSample) continue;

        for (int b = 1; b <= nBins; ++b)
            h2->SetBinContent(b, i + 1, hSample->GetBinContent(b));
    }

    // Compute mean ± stddev per bin
    TH1D* hMean = (TH1D*)hTemplate->Clone((varName + "_mean").c_str());
    TH1D* hStdDev = (TH1D*)hTemplate->Clone((varName + "_stddev").c_str());
    hMean->Reset();
    hStdDev->Reset();

    for (int b = 1; b <= nBins; ++b) {
        double sum = 0, sum2 = 0;
        for (int y = 1; y <= nSamples; ++y) {
            double val = h2->GetBinContent(b, y);
            sum += val;
            sum2 += val * val;
        }
        double mean = sum / nSamples;
        double stddev = sqrt(sum2 / nSamples - mean * mean);
        hMean->SetBinContent(b, mean);
        hMean->SetBinError(b, stddev); // for E2 drawing
        hStdDev->SetBinContent(b, stddev);
    }

    // Draw plot
    TCanvas* c = new TCanvas(("c_" + varName).c_str(), plotTitle.c_str(), 800, 600);
    hMean->SetLineColor(kBlack);
    hMean->SetLineWidth(2);
    hMean->SetMarkerSize(0);

    hMean->SetTitle(plotTitle.c_str());
    hMean->GetXaxis()->SetTitle(xLabel.c_str());
    hMean->GetYaxis()->SetTitle("Event Rate");

    hMean->Draw("E2");        // shaded ±1σ
    hMean->Draw("HIST SAME"); // mean line

    c->Print(pdfOutFile.c_str());

    delete h2;
    delete hMean;
    delete hStdDev;
    delete c;
}
*/

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
    c->Print(pdfOut.c_str());
    delete h_post;
}

/*
void Save1DSlices(TH2D* h_mean, TH2D* h_error,
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
        h_slice_mean->Draw("E2");        // shaded ±error
        h_slice_mean->Draw("HIST SAME"); // mean line on top
        c->Update();

        // Save one page
        c->Print(pdfOut.c_str());

        delete h_slice_mean;
        delete h_slice_error;
    }
}
*/
/*
void Save1DSlicesWithErrors(TH2D* h_mean, TH2D* h_stddev,
                            bool sliceAlongQ3,
                            const std::string& sliceTitle,
                            TCanvas* c,
                            const std::string& pdfOut) {
    int nSlices = sliceAlongQ3 ? h_mean->GetNbinsY() : h_mean->GetNbinsX();
    for (int i = 1; i <= nSlices; ++i) {
        c->Clear();

        TH1D* h_slice = sliceAlongQ3 ? h_mean->ProjectionX("slice", i, i)
                                     : h_mean->ProjectionY("slice", i, i);

        // Assign errors from h_stddev
        for (int bin = 1; bin <= h_slice->GetNbinsX(); ++bin) {
            int binX = sliceAlongQ3 ? bin : i;
            int binY = sliceAlongQ3 ? i : bin;
            double err = h_stddev->GetBinContent(binX, binY);
            h_slice->SetBinError(bin, err);
        }

        h_slice->SetTitle(Form("%s slice %d", sliceTitle.c_str(), i));
        h_slice->GetYaxis()->SetTitle("Event Rate");
        h_slice->Draw("E2");
        h_slice->Draw("HIST SAME");

        c->Print(pdfOut.c_str());
        delete h_slice;
    }
}
*/
/////////////////////////////////////////////////////////////////////////

/*
void Save1DSlices(TH2D* h2d, bool slice_along_q3, const std::string& label, TCanvas* c, const std::string& pdfOut,const std::string& tag,const std::string& xsec_var1,
                  const std::string& xsec_var2) {

  
    int nSlices = slice_along_q3 ? h2d->GetNbinsY() : h2d->GetNbinsX();
    std::string axisTitle = slice_along_q3 ? formatAxisTitle(xsec_var1): formatAxisTitle(xsec_var2);
    
    
    // Modify bin errors here
    for (int bin = 1; bin <= h_slice->GetNbinsX(); ++bin) {
        double newErr = your calculation;
        h_slice->SetBinError(bin, newErr);
    }


    for (int i = 1; i <= nSlices; ++i) {
        std::string slice_name = Form("%s_slice_%d", h2d->GetName(), i);
        TH1D* h_slice = slice_along_q3 ? h2d->ProjectionX(slice_name.c_str(), i, i)
                                       : h2d->ProjectionY(slice_name.c_str(), i, i);

        // Set axis labels
        h_slice->SetTitle(Form("%s Slice %d", label.c_str(), i));
        h_slice->GetXaxis()->SetTitle(axisTitle.c_str());
        h_slice->GetYaxis()->SetTitle("Posterior Value");

        h_slice->Draw("HIST E");
        c->Print(pdfOut.c_str());
    }
}*/

/*
void SaveSlices(TH2D* h_mean, TH2D* h_stddev,
                bool slice_along_q3,
                const std::string& label,
                TCanvas* c, const std::string& pdfOut,
                const std::string& tag,
                const std::string& xsec_var1,
                const std::string& xsec_var2,
                bool predictiveBand = false) 
{
    if (!h_mean) {
        std::cerr << "Error: SaveSlices called with null histogram!" << std::endl;
        return;
    }

    int nSlices = slice_along_q3 ? h_mean->GetNbinsY() : h_mean->GetNbinsX();
    std::string axisTitle = slice_along_q3 ? formatAxisTitle(xsec_var1) : formatAxisTitle(xsec_var2);

    // Color palette
    std::vector<int> colors = {kBlue, kRed, kGreen+2, kMagenta, kCyan+2, kOrange, kViolet, kSpring, kTeal};

    for (int i = 1; i <= nSlices; ++i) {
        // --- Build slice ---
        std::string slice_name = Form("%s_slice_%d", h_mean->GetName(), i);
        TH1D* h_slice = slice_along_q3 ? h_mean->ProjectionX(slice_name.c_str(), i, i)
                                       : h_mean->ProjectionY(slice_name.c_str(), i, i);

        // --- Assign errors from stddev ---
        if (h_stddev) {
            for (int bin = 1; bin <= h_slice->GetNbinsX(); ++bin) {
                int binX = slice_along_q3 ? bin : i;
                int binY = slice_along_q3 ? i : bin;
                double err = h_stddev->GetBinContent(binX, binY);
                h_slice->SetBinError(bin, err);
            }
        }

        // --- Titles ---
        double fixedLow  = slice_along_q3 ? h_mean->GetYaxis()->GetBinLowEdge(i) : h_mean->GetXaxis()->GetBinLowEdge(i);
        double fixedHigh = slice_along_q3 ? h_mean->GetYaxis()->GetBinUpEdge(i)  : h_mean->GetXaxis()->GetBinUpEdge(i);

        if (slice_along_q3) {
            h_slice->SetTitle(Form("%s: %.2f < q3 < %.2f", label.c_str(), fixedLow, fixedHigh));
        } else {
            h_slice->SetTitle(Form("%s: %.2f < q0 < %.2f", label.c_str(), fixedLow, fixedHigh));
        }

        h_slice->GetXaxis()->SetTitle(axisTitle.c_str());
        h_slice->GetYaxis()->SetTitle("Posterior Value");

        h_mean_slice->SetLineColor(kBlack);
        h_mean_slice->SetLineWidth(2);
        // --- Style ---
        int color = colors[(i-1) % colors.size()];
        h_slice->SetLineColor(color);
        h_slice->SetFillColorAlpha(color, 0.35);

        // --- Draw fresh on empty canvas ---
        c->Clear();
        if (predictiveBand) {
            h_slice->Draw("E2");        // shaded posterior predictive error band
            h_slice->Draw("HIST SAME"); // overlay mean line
        } else {
            h_slice->Draw("E2");
            h_slice->Draw("HIST SAME");
        }
        c->Update();

        // --- Save one page ---
        c->Print(pdfOut.c_str());

        delete h_slice;
    }
}
*/
/*
void SaveSlices(TH2* h_mean, TH2* h_stddev, bool sliceAlongQ3,
                const std::string& title, TCanvas* c, const std::string& pdfOut,
                const std::string& tag, const std::string& xLabel, const std::string& yLabel) 
{
    int nSlices = sliceAlongQ3 ? h_mean->GetNbinsX() : h_mean->GetNbinsY();

    for (int i = 1; i <= nSlices; i++) {
        c->Clear();

        // --- mean line ---
        TH1D* h_mean_slice = sliceAlongQ3 ? h_mean->ProjectionY("mean_slice", i, i)
                                          : h_mean->ProjectionX("mean_slice", i, i);

        // --- error band ---
        TH1D* h_err_slice  = sliceAlongQ3 ? h_stddev->ProjectionY("err_slice", i, i)
                                          : h_stddev->ProjectionX("err_slice", i, i);

        // Apply styles
        h_mean_slice->SetLineColor(kBlack);
        h_mean_slice->SetLineWidth(2);

        h_err_slice->SetFillColorAlpha(kAzure+1, 0.35); // transparent blue band
        h_err_slice->SetLineColor(kAzure+1);
        h_err_slice->SetMarkerSize(0);

        // Attach errors from stddev to mean slice
        for (int bin = 1; bin <= h_mean_slice->GetNbinsX(); ++bin) {
            double err = h_err_slice->GetBinContent(bin);
            h_mean_slice->SetBinError(bin, err);
        }

        // Titles
        h_mean_slice->SetTitle((title + " slice " + std::to_string(i)).c_str());
        h_mean_slice->GetXaxis()->SetTitle(sliceAlongQ3 ? yLabel.c_str() : xLabel.c_str());
        h_mean_slice->GetYaxis()->SetTitle("Value");

        // Draw
        h_mean_slice->Draw("E2");        // draw mean with error band (shaded)
        h_mean_slice->Draw("HIST SAME"); // black line on top

        c->Print((pdfOut + "(").c_str());
    }
}*/


/*
TH1D* MakePosteriorPredictiveHist(
    const std::string& name,
    int nbins, double xmin, double xmax,
    const std::vector<std::vector<double>>& posteriorSamples, // [sample][parameter]
    const std::vector<std::string>& xsec_vars,                // knob names
    samplePDFFDBase* pdf,                                           // your MaCh3 pdf
    const std::vector<int>& selectionVec = {}
) {
    TH1D* h_mean = new TH1D((name + "_mean").c_str(), name.c_str(), nbins, xmin, xmax);
    TH1D* h_var  = new TH1D((name + "_var").c_str(),  name.c_str(), nbins, xmin, xmax);

    for (size_t s = 0; s < posteriorSamples.size(); ++s) {
        // 1. set parameters for this sample
        for (size_t j = 0; j < xsec_vars.size(); ++j) {
            pdf->SetParameter(xsec_vars[j], posteriorSamples[s][j]);
            //pdf->reweight();
        }

        // 2. reweight with these parameters
        pdf->reweight();

        // 3. get a histogram prediction
        TH1D* h_tmp = pdf->MakePrediction(name + Form("_tmp%d", (int)s), nbins, xmin, xmax, selectionVec);

        // 4. accumulate mean & variance online
        for (int b = 1; b <= nbins; ++b) {
            double val = h_tmp->GetBinContent(b);
            double prev_mean = h_mean->GetBinContent(b);
            double prev_var  = h_var->GetBinContent(b);

            double new_mean = prev_mean + (val - prev_mean) / (s + 1);
            double new_var  = prev_var + (val - prev_mean) * (val - new_mean);

            h_mean->SetBinContent(b, new_mean);
            h_var->SetBinContent(b, new_var);
        }

        delete h_tmp;
    }

    // finalize variance → stddev
    for (int b = 1; b <= nbins; ++b) {
        if (posteriorSamples.size() > 1) {
            double var = h_var->GetBinContent(b) / (posteriorSamples.size() - 1);
            h_var->SetBinContent(b, sqrt(var));
        }
    }

    return h_mean; // keep h_var alongside it
}
*/


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
        c->Update();

        // Save one page
        c->Print(pdfOut.c_str());

        delete h_slice_mean;
        delete h_slice_error;
    }
}



// Skeleton for MakePosteriorPredictiveHist
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



int main(int argc, char* argv[]) {

   // std::string pdfOut = "EventRates.pdf";

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <mcmc_output.root> <bin_config.yaml>" << std::endl;
        return 1;
    }

    // Match usage string: first arg is MCMC ROOT, second is YAML config
    std::string mcmc_file = argv[1]; // ROOT file
    std::string yaml_file = argv[2]; // YAML config

    // === MCMC setup ===
    std::vector<BinDef> binDefs;
    std::vector<double> q0_edges;
    std::vector<double> q3_edges;

    manager* FitManager = new manager(yaml_file);
    if (!FitManager) {
        std::cerr << "Error: Failed to create FitManager from YAML config: " << yaml_file << std::endl;
        return 1;
    }
     std::string pdfOut = AddTimestampToFilename("mcmc_summary_plots.pdf");
    auto OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();
    auto OutputFile = std::unique_ptr<TFile>(TFile::Open(OutputFileName.c_str(), "RECREATE"));
    if (!OutputFile || OutputFile->IsZombie()) {
        std::cerr << "Failed to open output file: " << OutputFileName << std::endl;
        return 1;
    }

    auto xsec_var1 = FitManager->raw()["General"]["Systematics"]["xsec_var1"].as<std::string>();
    auto xsec_var2 = FitManager->raw()["General"]["Systematics"]["xsec_var2"].as<std::string>();
    auto xsec_yaml = FitManager->raw()["General"]["Systematics"]["XsecCovFile"][0].as<std::string>();
    if (xsec_yaml.empty()) {
        std::cerr << "Error: xsec_yaml path is empty or missing in config.\n";
        return 1;
    }

    BinningResult binning = extract_2D_bins_from_yaml(xsec_yaml, xsec_var1, xsec_var2);
    binDefs = binning.binDefs;
    q0_edges = binning.q0_edges;
    q3_edges = binning.q3_edges;
    if (q0_edges.empty() || q3_edges.empty()) {
        std::cerr << "Error: Bin edges could not be parsed from YAML file." << std::endl;
        return 1;
    }

    TAxis* axisX = new TAxis(q0_edges.size() - 1, q0_edges.data());
    TAxis* axisY = new TAxis(q3_edges.size() - 1, q3_edges.data());

    // MaCh3 setup
    covarianceXsec* xsec = nullptr;
    covarianceOsc* osc = nullptr;
    std::vector<samplePDFFDBase*> DUNEPdfs;
    MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec, osc);
    std::unique_ptr<mcmc> MaCh3Fitter = std::make_unique<mcmc>(FitManager);

    // Build combined histogram
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

    // Write combined histogram once
    if (h_q0q3_combined) {
        TFile outFile("my_eventratehistogram_templateparams_q0q3.root", "RECREATE");
        h_q0q3_combined->GetXaxis()->SetTitle(xsec_var1.c_str());
        h_q0q3_combined->GetYaxis()->SetTitle(xsec_var2.c_str());
        h_q0q3_combined->Write("xsec_eventrate_histo");
        outFile.Close();
    } else {
        std::cerr << "Error: No 2D histograms were produced.\n";
        return 1;
    }

    // --- Event rate histogram with adaptive label colors ---
    TH2D* h_eventrate_labels = (TH2D*)h_q0q3_combined->Clone("h_eventrate_labels");
    h_eventrate_labels->SetTitle("Event Rate with Bin Labels;q_{0} [GeV];q_{3} [GeV]");

    TCanvas* c_eventrate_labels = new TCanvas("c_eventrate_labels", "", 900, 700);
    c_eventrate_labels->SetRightMargin(0.15);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    h_eventrate_labels->Draw("COLZ");

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

    c_eventrate_labels->Print(pdfOut.c_str());
    h_eventrate_labels->Write();


    // Load posteriors
    TFile* f_post = new TFile(mcmc_file.c_str());
    TTree* post = (TTree*)f_post->Get("posteriors");
    if (!post) {
        std::cerr << "Error: Could not find 'posteriors' TTree in file " << mcmc_file << std::endl;
        return 1;
    }

    // Count xsec_* params
    int nParams = 0;
    TObjArray* branches = post->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); ++i) {
        std::string name = branches->At(i)->GetName();
        if (name.find("xsec_") == 0) {
            int idx = std::stoi(name.substr(5));
            if (idx >= nParams) nParams = idx + 1;
        }
    }
    std::cout << "Detected nParams = " << nParams << std::endl;

    // Prepare arrays for reading posteriors
    std::vector<double> xsec_vals(nParams, 0.0);
    std::vector<double*> xsec_ptrs(nParams, nullptr);
    for (int i = 0; i < nParams; ++i) xsec_ptrs[i] = &xsec_vals[i];

    post->SetBranchStatus("*", 0);
    for (const auto& bin : binDefs) {
        if (bin.index >= nParams) {
            std::cerr << "Skipping bin.index = " << bin.index << " (>= nParams)\n";
            continue;
        }
        std::string name = "xsec_" + std::to_string(bin.index);
        post->SetBranchStatus(name.c_str(), 1);
        post->SetBranchAddress(name.c_str(), xsec_ptrs[bin.index]);
    }

    // Histograms for mean/stddev
    TH2D* h_mean = new TH2D("h_param_mean", "Posterior Mean;q_{0} [GeV];q_{3} [GeV]",
                            q0_edges.size()-1, &q0_edges[0],
                            q3_edges.size()-1, &q3_edges[0]);
    TH2D* h_stddev = new TH2D("h_param_stddev", "Posterior StdDev;",
                              q0_edges.size()-1, &q0_edges[0],
                              q3_edges.size()-1, &q3_edges[0]);

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
    std::cout << "h_mean bins: X=" << h_mean->GetNbinsX() << " Y=" << h_mean->GetNbinsY() << std::endl;

    // Save histograms
    h_mean->Write("xsec_param_mean");
    h_stddev->Write("xsec_param_stddev");


    //////////////////Now for the Posteriors in 1D!


    ////////////////////////////////
    // Build posterior samples (with thinning)
    int thinFactor = 1000; // <-- adjust as needed
    std::vector<PosteriorSample> posteriorSamples;

    for (Long64_t i = 0; i < nEntries; i += thinFactor) {
        post->GetEntry(i);

        PosteriorSample sample;
        sample.xsec_draw.resize(nParams);
        for (int p = 0; p < nParams; ++p) {
            sample.xsec_draw[p] = xsec_vals[p];
        }
        posteriorSamples.push_back(sample);
    }

    std::cout << "Loaded " << posteriorSamples.size()
              << " posterior samples (thinned by " << thinFactor << ")." << std::endl;

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
    int nSamples = posteriorSamples.size();
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
    // Plotting setup
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);

    TCanvas* c = new TCanvas("c", "", 900, 700);
    c->SetRightMargin(0.15);

    // --- Open multipage PDF with first canvas (event rate with labels) ---
    c_eventrate_labels->Print((pdfOut + "[").c_str());
    c_eventrate_labels->Print(pdfOut.c_str());

    // --- Event rate (2D heatmap without labels) ---
    c->Clear();
    h_q0q3_combined->Draw("COLZ");
    c->Update();
    c->Print(pdfOut.c_str());

    // --- Posterior mean (2D) ---
    c->Clear();
    h_mean->Draw("COLZ");
    c->Update();
    c->Print(pdfOut.c_str());

    // --- Posterior stddev (2D) ---
    c->Clear();
    h_stddev->Draw("COLZ");
    c->Update();
    // --- Open multipage PDF ---
c->Print((pdfOut + "[").c_str());

// 1️⃣ Event rate 2D with labels
c_eventrate_labels->Print(pdfOut.c_str());

// 2️⃣ Event rate 2D without labels
c->Clear();
h_q0q3_combined->Draw("COLZ");
c->Update();
c->Print(pdfOut.c_str());

// 3️⃣ Posterior mean 2D
c->Clear();
h_mean->Draw("COLZ");
c->Update();
c->Print(pdfOut.c_str());

// 4️⃣ Posterior stddev 2D
c->Clear();
h_stddev->Draw("COLZ");
c->Update();
c->Print(pdfOut.c_str());

// 5️⃣ 1D slices of posterior mean & stddev
//Save1DSlices(h_mean, h_stddev, "Posterior Mean (q3 slices)", c, pdfOut, true, xsec_var1);
//Save1DSlices(h_mean, h_stddev, "Posterior Mean (q0 slices)", c, pdfOut, false, xsec_var2);

// 6️⃣ 1D slices of posterior predictive event rates
//Save1DSlices(h_eventrate_mean, h_eventrate_stddev, "Event Rate (q3 slices)", c, pdfOut, true, xsec_var1);
//Save1DSlices(h_eventrate_mean, h_eventrate_stddev, "Event Rate (q0 slices)", c, pdfOut, false, xsec_var2);

// Posterior mean slices
//Save1DSlices(h_mean, h_stddev, "Posterior Mean (q3 slices)", c, pdfOut, true, xsec_var1);
//Save1DSlices(h_mean, h_stddev, "Posterior Mean (q0 slices)", c, pdfOut, false, xsec_var2);

// Event rate slices with posterior predictive error
//Save1DSlices(h_eventrate_mean, h_eventrate_stddev, "Event Rate (q3 slices)", c, pdfOut, true, xsec_var1);
//Save1DSlices(h_eventrate_mean, h_eventrate_stddev, "Event Rate (q0 slices)", c, pdfOut, false, xsec_var2);

//Save1DSlicesWithErrors(h_eventrate_mean, h_eventrate_stddev, true, "Event rate along q3", c, pdfOut);
//Save1DSlicesWithErrors(h_eventrate_mean, h_eventrate_stddev, false, "Event rate along q0", c, pdfOut);

// Posterior mean slices
Save1DSlicesWithErrors(h_mean, h_stddev, "Posterior Mean (q3 slices)", c, pdfOut, true, xsec_var1);
Save1DSlicesWithErrors(h_mean, h_stddev, "Posterior Mean (q0 slices)", c, pdfOut, false, xsec_var2);

// Event rate slices with posterior predictive error
Save1DSlicesWithErrors(h_eventrate_labels, h_eventrate_stddev, "Event Rate (Relative Enu bias slices)", c, pdfOut, true, xsec_var1);
Save1DSlicesWithErrors(h_eventrate_labels, h_eventrate_stddev, "Event Rate (True Enu slices)", c, pdfOut, false, xsec_var2);

// 7️⃣ 1D posterior distributions
for (int p = 0; p < nParams; ++p) {
    SavePosterior1D(posteriorSamples, p, "xsec", c, pdfOut);
}

// --- Close multipage PDF ---
c->Print((pdfOut + "]").c_str());



    OutputFile->Close();
    std::cout << "Histograms saved to: " << pdfOut << std::endl;
      
    // Cleanup
    delete h_q0q3_combined;
    delete h_mean;
    delete h_stddev;
    delete axisX;
    delete axisY;
    delete c;
    delete f_post;
    delete FitManager;

    return 0;
}
