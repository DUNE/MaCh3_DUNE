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

std::vector<double> Params;
std::vector<double> AvgParams(480, 0.0);
std::vector<int> Count(480, 0);

void ChanOnlyParams()
{
    TFile* Osc = TFile::Open("TrueIndChanOsc.root"); // Loading in histograms
    TFile* Unosc = TFile::Open("TrueIndChanUnosc.root"); // Note that any "Reweight" file has the energy normalisation parameters already applied

    std::ofstream out("ParamsOutput.txt"); // Creating output text file

    TIter next(Osc->GetListOfKeys()); // Producing list of keys (histograms) from file
    TKey* key; 
    while ((key = (TKey*)next())) { 
        auto HistoOsc = Osc->Get<TH1D>(key->GetName()); // Get titles of individual histograms
        auto HistoUnosc = Unosc->Get<TH1D>(key->GetName());

        int NumBins = HistoOsc->GetNbinsX(); // Find number of bins (energy normalisation parameters) for each histogram
        for(int j = 1; j <= NumBins; j++) {
            double BinSizeOsc = HistoOsc->GetBinContent(j); // Find number of events for each energy parameter
            double BinSizeUnosc = HistoUnosc->GetBinContent(j);
            double Param;
            if(BinSizeUnosc == 0) { // If no unoscillated data, set value to 0, print output statement
                Param = 0;
                Params.push_back(Param);
            }
            else { // If unoscillated data, find ratio, print output statement
                Param = BinSizeOsc / BinSizeUnosc;
                Params.push_back(Param);
            }
        }
    }
    for(int p = 0; p < 1920; p++){
        int sample = p / 480;
        int channel = (p % 480) / 40;
        int bin = p % 40;
        int Index = channel * 40 + bin;
        AvgParams[Index] += Params[p];
        Count[Index] += 1;
    }
    for(int i = 0; i < 480; i++){
        AvgParams[i] /= Count[i];
    }
    for(const auto& val : AvgParams){
        out << val << "\n";
    }
}

void BinVals() // Create a script to print ratio of bins (params) for specific criteria
{
    TFile* Osc = TFile::Open("TrueIndChanOsc.root"); // Loading in histograms
    TFile* Unosc = TFile::Open("TrueIndChanUnosc.root"); // Note that any "Reweight" file has the energy normalisation parameters already applied

    std::ofstream out("ParamsOutput.txt"); // Creating output text file

    TIter next(Osc->GetListOfKeys()); // Producing list of keys (histograms) from file
    TKey* key; 
    while ((key = (TKey*)next())) { 
        auto HistoOsc = Osc->Get<TH1D>(key->GetName()); // Get titles of individual histograms
        auto HistoUnosc = Unosc->Get<TH1D>(key->GetName());

        out << "Printing All Params for Sample " << key->GetName() << "\n";

        int NumBins = HistoOsc->GetNbinsX(); // Find number of bins (energy normalisation parameters) for each histogram
        for(int j = 1; j <= NumBins; j++) {
            double BinSizeOsc = HistoOsc->GetBinContent(j); // Find number of events for each energy parameter
            double BinSizeUnosc = HistoUnosc->GetBinContent(j);
            double Param;
            if(BinSizeUnosc == 0) { // If no unoscillated data, set value to 0, print output statement
                Param = 0;
                Params.push_back(Param);
                //out << "Param in range " << HistoOsc->GetXaxis()->GetBinLowEdge(j) << " to " << HistoOsc->GetXaxis()->GetBinUpEdge(j) << " has a value of " << Param << " where the ratio is given by " << BinSizeOsc << " divided by " << BinSizeUnosc << "\n";
            }
            else { // If unoscillated data, find ratio, print output statement
                Param = BinSizeOsc / BinSizeUnosc;
                Params.push_back(Param);
                out << "Param in range " << HistoOsc->GetXaxis()->GetBinLowEdge(j) << " to " << HistoOsc->GetXaxis()->GetBinUpEdge(j) << " has a value of " << Param << " where the ratio is given by " << BinSizeOsc << " divided by " << BinSizeUnosc << "\n";
                //out << Param << "\n";
                // if((Param > 1.1) || (Param < 0.9)) // If param is greater than 10% away from 1 (ideal)
                // out << "Param for sample " << key->GetName() << " in range " << HistoOsc->GetXaxis()->GetBinLowEdge(j) << " to " << HistoOsc->GetXaxis()->GetBinUpEdge(j) << " has a value of " << Param << "\n";
                }
            }; 
        out << "\n";
    }
    for(const auto& x : Params){
        out << x << "\n";
    }
}

void IntegralComp() // Create a script to print ratio of integrals for two histograms
{
    TFile* Osc = TFile::Open("CCIndChanOsc.root"); // Loading in histograms
    TFile* Unosc = TFile::Open("TrueCCIndChanReweight.root"); // Note that any "Reweight" file has the energy normalisation parameters already applied

    std::ofstream out("IntComp.txt"); // Creating output text file

    TIter next(Osc->GetListOfKeys()); // Producing list of keys (histograms) from file
    TKey* key; 
    while ((key = (TKey*)next())) { 
        auto HistoOsc = Osc->Get<TH1D>(key->GetName()); // Get titles of individual histograms
        auto HistoUnosc = Unosc->Get<TH1D>(key->GetName());

        out << "Comparing Integrals for Sample " << key->GetName() << "\n";

        double IntOsc = HistoOsc->Integral(); // Calculate integrals of each histogram
        double IntUnosc = HistoUnosc->Integral();
        double IntRatio;
        if(IntUnosc == 0) { // If no integral possible for unoscillated data, set to 0
            IntRatio = 0;
            out << IntRatio << "\n";
        }
        else { // If integral possible, calculate ratio of these
            IntRatio = IntOsc / IntUnosc; 
            out << IntRatio << "\n";
            }
    }
}