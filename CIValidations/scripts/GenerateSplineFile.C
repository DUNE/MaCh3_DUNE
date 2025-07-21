#include "TFile.h"
#include "TH1D.h"
#include "TString.h"

#include <iostream>
#include <vector>
#include <string>

TString GenerateSplineName(TString systematic, TString mode, int etrue_bin, int erec_bin)
{
    TString spline_name;
    spline_name.Form("dev_%s_%s_sp_%d_%d", systematic.Data(), mode.Data(), etrue_bin, erec_bin);
    return spline_name;
}

void GenerateSplines(TFile *infile, TString systematic, std::vector<TString> modes, int etrue_bins, int erec_bins, std::vector<double> y_knots)
{
    infile->cd();
    int n_knots = static_cast<int>(y_knots.size());
    double min_x_knot = (static_cast<double>(n_knots) - 1.0) / 2.0;
    std::vector<double> x_knots(n_knots);
    for (int i = 0; i < n_knots; ++i)
    {
        x_knots[i] = min_x_knot + i;
    }

    for (const auto &mode : modes)
    {
        for (int etrue_bin = 0; etrue_bin < etrue_bins; ++etrue_bin)
        {
            for (int erec_bin = 0; erec_bin < erec_bins; ++erec_bin)
            {
                TString spline_name = GenerateSplineName(systematic, mode, etrue_bin, erec_bin);

                TSpline3* spline = new TSpline3(spline_name.Data(), x_knots.data(), y_knots.data(), n_knots);
                spline->SetName(spline_name.Data());
                spline->Write();
                delete spline;
            }
        }
    }
}

void GenerateSplineFile(TString output_file, TString input_spline_file)
{
    std::vector<TString> modes = {"ccqe", "ccmec", "ccdis", "ccres", "cccoh", "ccdiff", "ccnueel", "unknown", "ccamnugamma", "unknown",
                                  "cccohel", "ccibd", "ccglasres", "ccimdannihilation", "ncqe", "ncdis", "ncres", "nccoh", "ncdiff",
                                  "ncnueel", "ncamnugamma", "ncmec", "nccohel", "ncibd", "ncglasres", "ncimdannihilation"};

    std::vector<TString> cc_modes = {"ccqe"};

    std::vector<TString> nc_modes = {"ncqe"};

    // Literally just want to copy dev tmp
    TFile *input_file = TFile::Open(input_spline_file, "READ");
    if (!input_file || !input_file->IsOpen())
    {
        std::cerr << "Error opening file: " << input_spline_file << std::endl;
        return;
    }
    TH2D *h = (TH2D *)input_file->Get("dev_tmp_0_0");
    TH2D *dev_tmp_0_0 = static_cast<TH2D *>(h->Clone());
    dev_tmp_0_0->SetDirectory(nullptr); // Remove directory to avoid conflicts
    input_file->Close();

    TFile *infile = TFile::Open(output_file, "RECREATE");
    if (!infile || !infile->IsOpen())
    {
        std::cerr << "Error opening file: " << output_file << std::endl;
        return;
    }
    infile->cd();

    dev_tmp_0_0->Write();

    int x_bins = dev_tmp_0_0->GetNbinsX();
    int y_bins = dev_tmp_0_0->GetNbinsY();

    // CC splines
    std::vector<double> cc_y_knots = {1.0, 0.5, 2.0, 0.5, 1.0};                             // Example values
    std::vector<double> nc_y_knots = {2.0, 0.5, 2.0, 0.5, 2.0};                             // Example values
    std::vector<double> all_par_knots = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6}; // Example values

    // now we make splines
    std::cout << "Generating splines for systematic: CC" << std::endl;
    GenerateSplines(infile, "ccdummy", cc_modes, x_bins, y_bins, cc_y_knots);

    std::cout << "Generating splines for systematic: NC" << std::endl;
    GenerateSplines(infile, "ncdummy", nc_modes, x_bins, y_bins, nc_y_knots);
    std::cout << "Generating splines for systematic: ALL" << std::endl;

    // GenerateSplines(infile, "allmodesdummy", modes, x_bins, y_bins, all_par_knots);
    // std::cout << "Splines generated and written to file: " << output_file << std::endl;
    // infile->Close();
}
