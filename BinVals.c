#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TRandom.h"

void BinVals()
{
    TFile* Osc = TFile::Open("TrueIndChanOsc.root");
    TFile* Unosc = TFile::Open("TrueIndChanReweight.root");

    std::ofstream out("ParamsOutput.txt");

    TIter next(Osc->GetListOfKeys()); 
    TKey* key; 
    while ((key = (TKey*)next())) { 
        auto HistoOsc = Osc->Get<TH1D>(key->GetName());
        auto HistoUnosc = Unosc->Get<TH1D>(key->GetName());

        out << "Printing All Params for Sample " << key->GetName() << "\n";

        int NumBins = HistoOsc->GetNbinsX();

        for(int j = 1; j <= NumBins; j++) {
            double BinSizeOsc = HistoOsc->GetBinContent(j);
            double BinSizeUnosc = HistoUnosc->GetBinContent(j);

            double Param;
            if(BinSizeUnosc == 0) {
                Param = 0;
                //out << "Param in range " << HistoOsc->GetXaxis()->GetBinLowEdge(j) << " to " << HistoOsc->GetXaxis()->GetBinUpEdge(j) << " has a value of " << Param << " where the ratio is given by " << BinSizeOsc << " divided by " << BinSizeUnosc << "\n";
            }
            else {
                Param = BinSizeOsc / BinSizeUnosc; 
                out << "Param in range " << HistoOsc->GetXaxis()->GetBinLowEdge(j) << " to " << HistoOsc->GetXaxis()->GetBinUpEdge(j) << " has a value of " << Param << " where the ratio is given by " << BinSizeOsc << " divided by " << BinSizeUnosc << "\n";
                //out << Param << "\n";
                // if((Param > 1.1) || (Param < 0.9))
                // out << "Param for sample " << key->GetName() << " in range " << HistoOsc->GetXaxis()->GetBinLowEdge(j) << " to " << HistoOsc->GetXaxis()->GetBinUpEdge(j) << " has a value of " << Param << "\n";
                }
            }; 

        out << "\n";
    }
}

void IntegralComp()
{
    TFile* Osc = TFile::Open("CCIndChanOsc.root");
    TFile* Unosc = TFile::Open("TrueCCIndChanReweight.root");

    std::ofstream out("IntComp.txt");

    TIter next(Osc->GetListOfKeys()); 
    TKey* key; 

    while ((key = (TKey*)next())) { 
        auto HistoOsc = Osc->Get<TH1D>(key->GetName());
        auto HistoUnosc = Unosc->Get<TH1D>(key->GetName());

        out << "Comparing Integrals for Sample " << key->GetName() << "\n";

        int NumBins = HistoOsc->GetNbinsX();

        double IntOsc = HistoOsc->Integral();
        double IntUnosc = HistoUnosc->Integral();

        double IntRatio;
        if(IntUnosc == 0) {
            IntRatio = 0;
            out << IntRatio << "\n";
        }
        else {
            IntRatio = IntOsc / IntUnosc; 
            out << IntRatio << "\n";
            }
        //out << IntOsc << "\n";
    }
}