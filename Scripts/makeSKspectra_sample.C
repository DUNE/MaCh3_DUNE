#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TRandom3.h"
#include <vector>
#include <iostream>
#include <cmath>

// Make sure these are declared somewhere:
// extern std::vector<TString> samplenameuser;
// extern std::vector<int> maxysample;
// extern TString M3Mode[11];
// extern TString SplineMode[11];
// extern const int nSamples;
//YSP: The sample names are slightly different between OA2023 and OA2024
std::vector<TString> samplenameuser = {"ND_FHC_CCnumu"};
std::vector<int> maxysample = {2, 2, 3};
const int nSamples = 3;
TString M3Mode[11] = {"ccqe0", "cc1pi0", "ccnp0", "ccdis0", "nc0", "ccqe1", "cc1pi1", "ccnp1", "ccdis1", "nc1", "other"};
TString SplineMode[11] = {"ccqe", "cc1pi", "ccnp", "ccdis", "nc", "ccqe", "cc1pi", "ccnp", "ccdis", "nc", "other"};
//YSP: This sets the max y axis range for the posterior distributions. WIP to remove this altogether
//int maxysample[nSamples] = {150,150,150,25,25,5};

//Added this array to reduce all the if statements and string comparisons look cleaner now

void makeSKspectra_sample_by_mode(const char *infile = "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/Chain_testing/newconfig_eventrate.root", TString name = "hist", TString outname = "output.root", int maxy = 1, int nthrows = 1){
  TFile* file = new TFile(infile, "OPEN");
  std::cout << "infile is " << infile << std::endl;

  int nbins = maxy * 1000;
  //TString tempname = name + "_ccqe0";
  TString tempname = name + "_0";  // results in ND_FHC_CCnumu_0
 

  std::cout << "Got temp name as " << tempname << std::endl;
  std::cout << "Got name as " << name << std::endl;

  TH1D* temp = (TH1D*)file->Get(tempname);
  int nbinsx = temp->GetNbinsX();
  double binedges[10000];

  std::cout << "BLAH" << std::endl;
  printf("%i \n", nbinsx);

  for (int i = 0; i < nbinsx; i++) {
    binedges[i] = temp->GetXaxis()->GetBinLowEdge(i + 1);  // ROOT bins are 1-indexed
}
binedges[nbinsx] = temp->GetXaxis()->GetBinUpEdge(nbinsx);  // Add the final upper edge


  std::cout << "Made Th2D spectra" << std::endl;
  std::vector<TH2D*> spectra;
  std::vector<TH1D*> hAvg;
  TString mode = "";

  for (int mode_i = 0; mode_i < 11; mode_i++) {
    mode = M3Mode[mode_i];
    spectra.push_back(new TH2D("spectra_" + mode, "", nbinsx, binedges, nbins + 1, 0, maxy));
    hAvg.push_back(new TH1D("hAvg" + mode, "", 1000, 0, 500));
  }

  double average = 0;
  int nsteps = nthrows;
  std::cout << "Number of throws is " << nthrows << std::endl;

  TRandom3* rnd = new TRandom3(0);
  TFile* outfile = new TFile(outname, "RECREATE");

  for (int mode_i = 0; mode_i < spectra.size(); mode_i++) {
    mode = M3Mode[mode_i];

    for (int i = 0; i < nsteps; i++) {
      TString name1 = name + "_" + mode;
      name1 += i;

      temp = (TH1D*)file->Get(name1);
      if (!temp) {
        std::cout << "ERROR: Histogram '" << name1 << "' not found in file!" << std::endl;
        file->ls();
        continue;  // or return if this is critical
        }
        average += temp->Integral();
        hAvg[mode_i]->Fill(temp->Integral());

      for (int j = 1; j <= temp->GetNbinsX(); j++)
        spectra[mode_i]->Fill(temp->GetXaxis()->GetBinCenter(j), temp->GetBinContent(j));
    }

    std::cout << "Filled spectra " << std::endl;
  }

  ///////////////////////////////
  // Now make projections
  ///////////////////////////////

  std::cout << "Now looping through to make projections" << std::endl;
  double x[100], y[100], ex[100], ey[100];
  TH1D* proj = nullptr;
  std::vector<TGraphErrors*> errorbars;
  TString dummyname;
  TString histname;
  std::vector<TH1D*> hist1d;

  for (int mode_i = 0; mode_i < spectra.size(); mode_i++) {
    mode = SplineMode[mode_i];
    if (name == "sk_nue") dummyname = "nue_" + mode;
    if (name == "sk_nuebar") dummyname = "nuebar_" + mode;
    if (name == "sk_nue1pi") dummyname = "nue1pi_" + mode;
    if (name == "sk_numu") dummyname = "numu_" + mode;
    if (name == "sk_numubar") dummyname = "numubar_" + mode;

    std::cout << "Starting mode " << mode << " " << mode_i << " / " << spectra.size() << std::endl;

    histname = "hist_1d_" + dummyname;
    hist1d.push_back((TH1D*)temp->Clone(histname));
    hist1d[mode_i]->Reset();
    hist1d[mode_i]->GetYaxis()->SetRange(0, 20);

    if (name == "sk_numu" || name == "sk_numubar") {
      if (mode_i != 0 && mode_i != 4) {
        spectra[mode_i]->RebinX(4);
        hist1d[mode_i]->RebinX(4);
      } else if (mode_i == 4) {
        spectra[mode_i]->RebinX(8);
        hist1d[mode_i]->RebinX(8);
      }
    } else {
      if (mode_i > 3) {
        spectra[mode_i]->RebinX(2);
        hist1d[mode_i]->RebinX(2);
      }
    }

    for (int i = 1; i <= spectra[mode_i]->GetNbinsX(); i++) {
      TString namepy = name + "_" + mode + "_py";
      namepy += i;
      proj = spectra[mode_i]->ProjectionY(namepy, i, i);

      int first_bin = proj->GetNbinsX() + 1;
      int last_bin = 0;

      for (int bin_i = 1; bin_i <= proj->GetNbinsX(); ++bin_i) {
        double val = proj->GetBinContent(bin_i);
        if (val > 0 && bin_i < first_bin) first_bin = bin_i;
        if (val > 0 && bin_i > last_bin) last_bin = bin_i;
      }

      double diff = proj->GetBinCenter(last_bin) - proj->GetBinCenter(first_bin);
      double opt_bin_width = diff / 25;
      double current_bin_width = proj->GetBinWidth(last_bin);
      int rebin = 1;
      if (diff > 0) rebin = int(opt_bin_width / current_bin_width);
      if (rebin == 0) rebin = 1;
      while ((proj->GetNbinsX() % rebin) != 0 && rebin > 1) rebin--;

      if (rebin > 1) proj->Rebin(rebin);

      x[i - 1] = spectra[mode_i]->GetXaxis()->GetBinCenter(i);
      ex[i - 1] = spectra[mode_i]->GetXaxis()->GetBinWidth(i) / 2.0;
      y[i - 1] = proj->GetMean();
      ey[i - 1] = proj->GetRMS();

      if (y[i - 1] > 1E-4) {
        if (mode_i == 0) std::cout << "Filling bin " << i << " with " << y[i - 1] << std::endl;
        hist1d[mode_i]->SetBinContent(i, y[i - 1]);
        hist1d[mode_i]->SetBinError(i, ey[i - 1]);
      } else {
        hist1d[mode_i]->SetBinContent(i, 0);
        hist1d[mode_i]->SetBinError(i, 0);
      }

      proj->Write(namepy);
    }

    errorbars.push_back(new TGraphErrors(nbinsx, x, y, ex, ey));
    errorbars[mode_i]->Write("graph_" + dummyname);
    spectra[mode_i]->Write("spectra2d_" + dummyname);
    hist1d[mode_i]->Write("hist_1d_" + dummyname);
    hAvg[mode_i]->Write("hist_1d_" + dummyname + "_events");
  }

  file->Close();
  outfile->Close();
}


void makeSKspectra(const char *infile = "/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/Chain_testing/newconfig_eventrate.root", TString outname = "output.root", int nthrows = 1, int do_var = -1, bool do_by_mode = false, int which_oa = 2023) {
  if (do_by_mode) {
    for (int i = 0; i < nSamples; ++i) {
      std::cout << std::endl << " ---------------- " << samplenameuser[i] << " ---------------- " << std::endl << std::endl;
      makeSKspectra_sample_by_mode(infile, samplenameuser[i], outname, maxysample[i], nthrows);
    }
  }
}


// ROOT entry point
void makeSKspectra_sample() {
  makeSKspectra("/scratch/abipeake/MaCh3_DUNE_merged/MaCh3_DUNE/Chain_testing/newconfig_eventrate.root", "output.root", 1, -1, true, 2023);
}
