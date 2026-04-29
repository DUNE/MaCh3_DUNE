#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <memory>

#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>
#include <TFile.h>

#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Fitters/MaCh3Factory.h"


// ============================================================================
// Write histograms to ROOT file
// ============================================================================

void Write1DHistogramsToFile(
    const std::string& OutFileName,
    const std::vector<TH1*>& Histograms)
{
    auto OutputFile = std::unique_ptr<TFile>(
        TFile::Open(OutFileName.c_str(), "RECREATE")
    );

    OutputFile->cd();

    for (auto Hist : Histograms) {
        Hist->Write();
    }

    OutputFile->Close();
}


// ============================================================================
// Write histograms to PDF
// ============================================================================

void Write1DHistogramsToPdf(
    std::string OutFileName,
    const std::vector<TH1*>& Histograms)
{
    // Replace .root with .pdf
    OutFileName.erase(OutFileName.find('.'));
    OutFileName += ".pdf";

    auto c1 = std::unique_ptr<TCanvas>(
        new TCanvas("c1", "c1", 800, 600)
    );

    c1->cd();

    c1->Print((OutFileName + "[").c_str());

    for (auto Hist : Histograms) {
        Hist->Draw("HIST");
        c1->Print(OutFileName.c_str());
    }

    c1->Print((OutFileName + "]").c_str());
}



void Print2DHistogramTable(TH2* hist)
{
    double totalEvents = hist->Integral();

    std::string xVarName = hist->GetXaxis()->GetTitle();
    std::string yVarName = hist->GetYaxis()->GetTitle();

    if (xVarName.empty()) xVarName = "X Variable";
    if (yVarName.empty()) yVarName = "Y Variable";

    
    std::cout << "Histogram: " << hist->GetName() << "\n";
    std::cout << "Total events: "
              << std::fixed << std::setprecision(2)
              << totalEvents << "\n";
    std::cout << "Only showing template bins with > 0.01% contribution\n";
    

    std::cout
        << std::left
        << std::setw(25) << (xVarName + " Bin Range")
        << std::setw(25) << (yVarName + " Bin Range")
        << std::setw(15) << "Events"
        << std::setw(15) << "% of Total"
        << "\n";

    std::cout << std::string(80, '-') << "\n";

    for (int ix = 1; ix <= hist->GetNbinsX(); ++ix)
    {
        double xLow  = hist->GetXaxis()->GetBinLowEdge(ix);
        double xHigh = hist->GetXaxis()->GetBinUpEdge(ix);

        for (int iy = 1; iy <= hist->GetNbinsY(); ++iy)
        {
            double yLow  = hist->GetYaxis()->GetBinLowEdge(iy);
            double yHigh = hist->GetYaxis()->GetBinUpEdge(iy);

            double content = hist->GetBinContent(ix, iy);

            double percent = 0.0;
            if (totalEvents > 0.0)
            {
                percent = 100.0 * content / totalEvents;
            }

            // Only print bins contributing more than 1%
            if (percent < 0.01) continue;

            std::stringstream xRange;
            xRange << "["
                   << std::fixed << std::setprecision(2)
                   << xLow << ", " << xHigh << "]";

            std::stringstream yRange;
            yRange << "["
                   << std::fixed << std::setprecision(2)
                   << yLow << ", " << yHigh << "]";

            std::cout
                << std::left
                << std::setw(25) << xRange.str()
                << std::setw(25) << yRange.str()
                << std::setw(15) << std::fixed << std::setprecision(2)
                << content
                << std::setw(15) << std::fixed << std::setprecision(2)
                << percent
                << "\n";
        }
    }

    std::cout << "========================================================================\n";
}


int main(int argc, char* argv[])
{
    MaCh3Utils::MaCh3Usage(argc, argv);
    auto fitMan = MaCh3ManagerFactory(argc, argv);

    
    ParameterHandlerGeneric* xsec = nullptr;

    std::vector<SampleHandlerFD*> DUNEPdfs;
    MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);


    
    std::vector<TH1*> DUNEHists;

    for (auto handler : DUNEPdfs) {

        handler->Reweight();

        for (int iSample = 0; iSample < handler->GetNsamples(); ++iSample) {

            TH1* hist = handler->GetMCHist(iSample);

            
            double totalEvents = hist->Integral(); // True total event count

            
            TH1* plotHist = (TH1*)hist->Clone(
                Form("%s_binWidthNorm", hist->GetName())
            );

            plotHist->Scale(1.0, "width");
            plotHist->GetYaxis()->SetTitle("Events / GeV");

           
            hist->SetTitle(
                Form("%s (Raw Counts)", hist->GetTitle())
            );

            plotHist->SetTitle(
                Form("%s (Bin Width Normalized)", hist->GetTitle())
            );

            
            //DUNEHists.push_back(hist);
            DUNEHists.push_back(plotHist);

            std::string EventRateString =
                fmt::format("{:.2f}", totalEvents);

            MACH3LOG_INFO(
                "Event rate for {} : {:<5}",
                handler->GetSampleTitle(iSample),
                EventRateString
            );

            handler->PrintIntegral(iSample);

            
            std::vector<KinematicCut> SelectionVec;

            

            TH2* Hist2D = handler->Get2DVarHist(
                iSample,
                "TrueNeutrinoEnergy",   // X variable
                "Enubias",   // Y variable
                SelectionVec
            );

            Print2DHistogramTable(Hist2D);

            delete Hist2D;
        }
    }

    std::string OutFileName =
        GetFromManager<std::string>(
            fitMan->raw()["General"]["OutputFile"],
            "EventRatesOutput.root"
        );

    Write1DHistogramsToFile(OutFileName, DUNEHists);
    Write1DHistogramsToPdf(OutFileName, DUNEHists);


    // =========================================================================
    // Oscillation Mode Breakdown
    // =========================================================================

    MACH3LOG_INFO("========================================================================");
    MACH3LOG_INFO("Oscillation Mode Breakdown:");

    for (auto handler : DUNEPdfs) {

        for (int iSample = 0; iSample < handler->GetNsamples(); ++iSample) {

            MACH3LOG_INFO("======================");

            int nOscChannels = handler->GetNOscChannels(iSample);

            for (int iOscChan = 0; iOscChan < nOscChannels; ++iOscChan) {

                std::vector<KinematicCut> SelectionVec;

                KinematicCut SelecChannel;
                SelecChannel.ParamToCutOnIt =
                    handler->ReturnKinematicParameterFromString(
                        "OscillationChannel"
                    );

                SelecChannel.LowerBound = iOscChan;
                SelecChannel.UpperBound = iOscChan + 1;

                SelectionVec.push_back(SelecChannel);

                TH1* Hist = handler->Get1DVarHist(
                    iSample,
                    handler->GetXBinVarName(iSample),
                    SelectionVec
                );

                MACH3LOG_INFO(
                    "{:<20} : {:<20} : {:<20.2f}",
                    handler->GetSampleTitle(iSample),
                    handler->GetFlavourName(iSample, iOscChan),
                    Hist->Integral()
                );
            }
        }
    }


    // =========================================================================
    // Interaction Mode Breakdown
    // =========================================================================

    MACH3LOG_INFO("========================================================================");
    MACH3LOG_INFO("Interaction Mode Breakdown:");

    for (auto handler : DUNEPdfs) {

        for (int iSample = 0; iSample < handler->GetNsamples(); ++iSample) {

            MACH3LOG_INFO("======================");

            MaCh3Modes* Modes = handler->GetMaCh3Modes();
            int nModeChannels = Modes->GetNModes();

            for (int iModeChan = 0; iModeChan < nModeChannels; ++iModeChan) {

                std::vector<KinematicCut> SelectionVec;

                KinematicCut SelecChannel;
                SelecChannel.ParamToCutOnIt =
                    handler->ReturnKinematicParameterFromString("Mode");

                SelecChannel.LowerBound = iModeChan;
                SelecChannel.UpperBound = iModeChan + 1;

                SelectionVec.push_back(SelecChannel);

                TH1* Hist = handler->Get1DVarHist(
                    iSample,
                    handler->GetXBinVarName(iSample),
                    SelectionVec
                );

                MACH3LOG_INFO(
                    "{:<20} : {:<20} : {:<20.2f}",
                    handler->GetSampleTitle(iSample),
                    Modes->GetMaCh3ModeName(iModeChan),
                    Hist->Integral()
                );
            }
        }
    }

    return 0;
}