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

// ############## INITIAL VALUES ############## //

const int NParams = 480; // Number of energy normalisation parameters

const int NEBins = 40; // Number of energy bins per channel
const double EMin = 0.0;
const double EMax = 10.0;
const double BinWidth = (EMax - EMin) / NEBins;

const int NProbBins = 100; // Number of probability bins in the colour plots
const double ProbMin = 0.0;
const double ProbMax = 1.0;
const double ProbBinWidth = (ProbMax - ProbMin)/ NProbBins;

// ############## INITIAL VECTORS/MAPPINGS ############## //

std::vector<double> NormParams(NParams, 0.0); // Creating empty vector to store all energy normalisation parameters

std::map<int, std::string> NuTypes = {{12, "nue"}, {14, "numu"}, {16, "nutau"}, {-12, "nuebar"}, {-14, "numubar"}, {-16, "nutaubar"}}; // Creating map between neutrino "value" and "name"
std::vector<std::pair<int,int>> OscChans = {{14, 14}, {12, 12}, {-14, -14}, {-12, -12}, {14, 12}, {-14, -12}, {12, 14}, {-12, -14}, {14, 16}, {12, 16}, {-14, -16}, {-12, -16}}; // Creating all oscillation channels as pairs of neutrino "value"

std::vector<std::string> SampleTitles = {"FHC_numu"};
//std::vector<std::string> Sample = {"FHC_numu", "FHC_nue", "RHC_numu", "RHC_nue"}; // Creating vector of sample titles

std::vector<TH2D*> EventHistos; // Creating empty vector to store histograms
std::vector<double> BinLowEdges;
std::vector<double> ProbBinLowEdges;
std::vector<double> BinCentres;

// ############## CREATION OF ENERGY AND PROBABILITY BINNING SCHEMES ############## //

void SetupBinning(){
    for (int a = 0; a <= NEBins; a++){
        double BinLow = EMin + a * BinWidth;
        BinLowEdges.push_back(BinLow);
    }
    for (int b = 0; b < NEBins; b++){
        double BinMid = (BinLowEdges[b] + BinLowEdges[b+1]) / 2;
        BinCentres.push_back(BinMid);
    }
    for (int c = 0; c <= NProbBins; c++){
        double ProbBinLow = ProbMin + c * ProbBinWidth;
        ProbBinLowEdges.push_back(ProbBinLow);
    }
}

// ############## READING DATA AND CREATING PLOTS TO PRINT ############## //

void NormaliseAndPlot(TH2D* Histo, TString HistoName) // Create a script to column normalise and then plot a given histogram with a given name
{
    for (int x = 1; x <= Histo->GetNbinsX(); x++) { // For all energy bins, avoiding underflow bin
        bool FixedCol = true; // Assume a column is for a fixed parameter to start
        for (int y = 2; y <= Histo->GetNbinsY(); y++) { // For every bin except the first (if fixed, all events will be in bin 1)
            if (Histo->GetBinContent(x,y) != 0.0) { // If bin content is not zero
                FixedCol = false; // Not a fixed parameter - only need one non-zero bin at not zero to prove this
                break;
            }
        }
        if (FixedCol) { // If not fixed then skip this, if fixed then;
            for (int y = 1; y <= Histo->GetNbinsY(); y++) {
                Histo->SetBinContent(x, y, std::numeric_limits<double>::quiet_NaN()); // Set all bin contents to be NaN, will all show up dark blue
            }
            continue; 
        }
        double ColSum = 0;
        for (int y = 1; y <= Histo->GetNbinsY(); y++) { // For all oscillation probability bins, avoiding underflow bin
            ColSum += Histo->GetBinContent(x,y); // Running total of the number of events in the column
        }
        if (ColSum == 0) continue;
        for (int y = 1; y <= Histo->GetNbinsY(); y++) { // For all oscillation probability bins, avoiding underflow bin
            double Content = Histo->GetBinContent(x,y); // Get the content of each bin
            Histo->SetBinContent(x, y, Content/ColSum); // Normalise each bin by the total number of events in the column
        }
    }

    // TFile* OscData = TFile::Open("CCIndChanOsc.root"); // Loading in histograms
    // TFile* UnoscData = TFile::Open("CCIndChanUnosc.root");

    // auto HistoOsc = OscData->Get<TH1D>("FHC_numu_FHC_numu_x_numu");
    // auto HistoUnosc = UnoscData->Get<TH1D>("FHC_numu_FHC_numu_x_numu");

    // TH1D* NominalLine = new TH1D("Nominals", "Ideal Values for Params", NEBins, BinLowEdges.data());

    // for (int j = 1; j <= HistoOsc->GetNbinsX(); j++){
    //     NominalLine->SetBinContent(j, HistoOsc->GetBinContent(j) / HistoUnosc->GetBinContent(j));
    // }

    auto Canvas = new TCanvas(); // Create a canvas to plot
    Canvas->Draw(); // Draw the canvas
    Canvas->cd();
    Canvas->SetLogy(); // Set log scale
    Histo->Draw("COLZ"); // Draw colour scale histogram
    Histo->GetXaxis()->SetTitle("Energy (GeV)"); // Setting axes properties
    Histo->GetYaxis()->SetTitle("Oscillation Probability");
    Histo->GetZaxis()->SetTitle("Number of Events");
    Histo->GetZaxis()->SetTitleOffset(1.5);
    Histo->SetStats(0); // Cosmetics of plot
    // NominalLine->Draw("SAME");
    // NominalLine->SetStats(0);
    // NominalLine->SetLineColor(kRed); 
    // NominalLine->SetLineWidth(2);
    Canvas->SetRightMargin(0.15);
    Canvas->Update();
    Canvas->Print(HistoName + "_OscProbCol.pdf"); // Print colour histogram to pdf
}

// ############## LOADING BIN SCHEME, DATA, AND PLOTTING ############## //

void ChanColPlot(TString File_To_Plot) // Create a script to produce colour plots for each oscillation channel from a given root file
{
    SetupBinning();

    auto RootFile = TFile::Open(File_To_Plot); // Open the specified file

    if (RootFile->IsZombie() || !RootFile) { // Check if file exists
        throw std::runtime_error("File does not exist");
    }

    auto Posteriors = RootFile->Get<TTree>("posteriors"); // Read in all posteriors from file

    for (int i = 0; i < NParams; i++) { // For all energy normalisation parameters
        std::string ChainVarName = "param_" + std::to_string(i); // Create relevent variable name
        Posteriors->SetBranchAddress(ChainVarName.c_str(), &NormParams[i]); // Reassign address to correct name and store values in given vector
    }

    for (const auto& sample : SampleTitles) { // For each sample in our vector of sample names
        for (const auto& channel : OscChans) { // For each channel in our pairing of oscillation channels
            std::string UnoscName = NuTypes[channel.first]; // Unoscillated neutrino type is first value in pairing
            std::string OscName   = NuTypes[channel.second]; // Oscillated neutrino type is second value in pairing
            std::string HistName = sample + "_" + UnoscName + "_" + OscName; // Create name for each histogram using sample name and oscillation channel
            EventHistos.push_back(new TH2D(HistName.c_str(), ("Oscillation Probability for " + HistName).c_str(), NEBins, BinLowEdges.data(), NProbBins, ProbBinLowEdges.data())); // Create histogram for each channel in each sample
        }
    }

    Long64_t MCEntries = Posteriors->GetEntries(); // Get number of entries (steps) in given file
    const int BurnInCut = 0.25 * MCEntries; // Number of steps to cut out for burn in
    for (Long64_t i = BurnInCut; i < MCEntries; i++) { // Start at some burn-in cut, run for rest of entries
        Posteriors->GetEntry(i); // Get the relevant information
        for (int p = 0; p < NParams; p++) { // For each parameter
            int k = p % NEBins; // Find the value of param index modulo 40 (each value corresponds to a different energy bin (40))
            int g = p / NEBins; // Find the integer division of param index / 48 (each value corresponds to a different histogram (48))
            EventHistos[g]->Fill(BinCentres[k], NormParams[p]); // Fill each histogram with the relevant parameter from the desired energy bin
        }
    }

    // for (auto* IndChan : EventHistos){ // For each histogram in our list
    //     NormaliseAndPlot(IndChan, IndChan->GetName()); // Normalise by column and plot, getting names as required
    // }

    NormaliseAndPlot(EventHistos[0], EventHistos[0]->GetName()); // FHC_numu: numu -> numu
    //NormaliseAndPlot(EventHistos[4], EventHistos[4]->GetName());
    //NormaliseAndPlot(Hists[16], Hists[16]->GetName()); // FHC_nue: numu -> nue 
    //NormaliseAndPlot(Hists[26], Hists[26]->GetName()); // RHC_numu: numubar -> numubar
    //NormaliseAndPlot(Hists[41], Hists[41]->GetName()); // RHC_nue: numubar -> nuebar
}

// ###################################################################### //