#include <vector>
#include <utility>
#include <string>
#include <map>
#include <limits>
#include <iostream>
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"

// ############## CHANGE AS NEEDED ############## //
const int NParams = 1920; // Number of energy normalisation parameters
const int NEBins = 20; // Number of probability bins in the colour plots
std::vector<std::string> Sample = {"FHC_numu", "FHC_nue", "RHC_numu", "RHC_nue"}; // Creating vector of sample titles
//std::vector<std::string> Sample = {"FHC_numu"};

// ############## FIXED PARAMS/VALUES ############## //
const int EBinsChan = 40; // Number of energy bins per channel
std::vector<double> BinStart = {0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0,4.25,4.5,4.75,5.0,5.25,5.5,5.75,6.0,6.25,6.5,6.75,7.0,7.25,7.5,7.75,8.0,8.25,8.5,8.75,9.0,9.25,9.5,9.75,10.0};
std::vector<double> BinMid = {0.125,0.375,0.625,0.875,1.125,1.375,1.625,1.875,2.125,2.375,2.625,2.875,3.125,3.375,3.625,3.875,4.125,4.375,4.625,4.875,5.125,5.375,5.625,5.875,6.125,6.375,6.625,6.875,7.125,7.375,7.625,7.875,8.125,8.375,8.625,8.875,9.125,9.375,9.625,9.875};
std::vector<double> NormParams(NParams, 0.0); // Creating empty vector to store all energy normalisation parameters
std::map<int, std::string> NuType = {{12, "nue"}, {14, "numu"}, {16, "nutau"}, {-12, "nuebar"}, {-14, "numubar"}, {-16, "nutaubar"}}; // Creating map between neutrino "value" and "name"
std::vector<std::pair<int,int>> OscChan = {{14, 14}, {12, 12}, {-14, -14}, {-12, -12}, {14, 12}, {-14, -12}, {12, 14}, {-12, -14}, {14, 16}, {12, 16}, {-14, -16}, {-12, -16}}; // Creating all oscillation channels as pairs of neutrino "value"
std::vector<TH2D*> Hists; // Creating empty vector to store histograms

void NormaliseAndPlot(TH2D* h, TString HistName) // Create a script to column normalise and then plot a given histogram with a given name
{
    for (int x = 1; x <= h->GetNbinsX(); x++) { // For all energy bins, avoiding underflow bin
        bool FixedCol = true; // Assume a column is for a fixed parameter to start
        for(int y = 2; y <= h->GetNbinsY(); y++){ // For every bin except the first (if fixed, all events will be in bin 1)
            if(h->GetBinContent(x,y) != 0.0){ // If bin content is not zero
                FixedCol = false; // Not a fixed parameter - only need one non-zero bin at not zero to prove this
                break;
            }
        }
        if (FixedCol) { // If not fixed then skip this, if fixed then;
        for (int y = 1; y <= h->GetNbinsY(); y++) {
            h->SetBinContent(x, y, std::numeric_limits<double>::quiet_NaN()); // Set all bin contents to be NaN, will all show up dark blue
        }
        continue; 
        }
        double ColSum = 0;
        for (int y = 1; y <= h->GetNbinsY(); y++) { // For all oscillation probability bins, avoiding underflow bin
            ColSum += h->GetBinContent(x,y); // Running total of the number of events in the column
        }
        if (ColSum == 0) continue;
        for (int y = 1; y <= h->GetNbinsY(); y++) { // For all oscillation probability bins, avoiding underflow bin
            double Content = h->GetBinContent(x,y); // Get the content of each bin
            h->SetBinContent(x, y, Content/ColSum); // Normalise each bin by the total number of events in the column
        }
    }

    auto c = new TCanvas(); // Create a canvas to plot
    c->Draw(); // Draw the canvas
    c->cd();
    //c->SetLogy(); // Set log scale
    h->Draw("COLZ"); // Draw colour scale histogram
    h->GetXaxis()->SetTitle("Energy (GeV)"); // Setting axes properties
    h->GetYaxis()->SetTitle("Oscillation Probability");
    h->GetZaxis()->SetTitle("Number of Events");
    h->GetZaxis()->SetTitleOffset(1.5);
    h->SetStats(0); // Cosmetics of plot
    c->SetRightMargin(0.15);
    c->Update();
    c->Print(HistName + "_OscProbCol.pdf"); // Print colour histogram to pdf
}

void ChanColPlot(TString file_name) // Create a script to produce colour plots for each oscillation channel from a given root file
{
    auto f = TFile::Open(file_name); // Open the specified file

    if(f->IsZombie() || !f){ // Check if file exists
        std::cerr<<"Error cannot find "<<file_name<<std::endl;
        throw;
    }

    auto t1 = f->Get<TTree>("posteriors"); // Read in all posteriors from file

    for(int i = 0; i < NParams; i++){ // For all energy normalisation parameters
        std::string name = "xsec_" + std::to_string(i); // Create relevent xsec name
        t1->SetBranchAddress(name.c_str(), &NormParams[i]); // Reassign address to correct name and store values in given vector
    }

    for (const auto& sample : Sample) { // For each sample in our vector of sample names
        for (const auto& chan : OscChan) { // For each channel in our pairing of oscillation channels
            std::string unosc_name = NuType[chan.first]; // Unoscillated neutrino type is first value in pairing
            std::string osc_name   = NuType[chan.second]; // Oscillated neutrino type is second value in pairing
            std::string HistName = sample + "_" + unosc_name + "_" + osc_name; // Create name for each histogram using sample name and oscillation channel
            Hists.push_back(new TH2D(HistName.c_str(), ("Oscillation Probability for " + HistName).c_str(), EBinsChan, BinStart.data(), NEBins, 0.0, 0.1)); // Create histogram for each channel in each sample
        }
    }

    Long64_t nentries = t1->GetEntries(); // Get number of entries (steps) in given file
    for (Long64_t i=4000000; i<nentries; i++) { // Start at some burn-in cut, run for rest of entries
        t1->GetEntry(i); // Get the relevant information
        for(int p = 0; p < NParams; p++){ // For each parameter
            int k = p % EBinsChan; // Find the value of param index modulo 40 (each value corresponds to a different energy bin (40))
            int g = p / EBinsChan; // Find the integer division of param index / 48 (each value corresponds to a different histogram (48))
            Hists[g]->Fill(BinMid[k], NormParams[p]); // Fill each histogram with the relevant parameter from the desired energy bin
        }
    }

    // for (auto* h : Hists){ // For each histogram in our list
    //     NormaliseAndPlot(h, h->GetName()); // Normalise by column and plot, getting names as required
    // }

    //NormaliseAndPlot(Hists[0], Hists[0]->GetName()); // FHC_numu: numu -> numu
    NormaliseAndPlot(Hists[16], Hists[16]->GetName()); // FHC_nue: numu -> nue 
    //NormaliseAndPlot(Hists[26], Hists[26]->GetName()); // RHC_numu: numubar -> numubar
    NormaliseAndPlot(Hists[41], Hists[41]->GetName()); // RHC_nue: numubar -> nuebar
}