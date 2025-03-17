#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <string>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>

#include "samplePDFDUNE/MaCh3DUNEFactory.h"
#include "samplePDFDUNE/StructsDUNE.h"
#include "tests/Comparison.h"

// Function to write to both file and terminal
void writeToBoth(std::ofstream& outFile, const std::string& text) {
    outFile << text << std::endl; // Write to file with std::endl
    MACH3LOG_INFO(text); // Write to terminal with std::endl
}

int main(int argc, char * argv[]) {
    MaCh3Utils::MaCh3Usage(argc, argv);

    auto fitMan = std::unique_ptr<manager>(new manager(argv[1]));
    
    bool skip_checks = false;

    if (CheckNodeExists(fitMan->raw(), "General", "Tests", "SkipChecks")) {
        skip_checks = fitMan->raw()["General"]["Tests"]["SkipChecks"].as<bool>();
    }

    //###############################################################################################################################
    //Create samplePDFFD objects
    
    covarianceXsec* xsec = nullptr;
    covarianceOsc* osc = nullptr;
    
    std::vector<samplePDFFDBase*> DUNEPdfs;
    MakeMaCh3DuneInstance(fitMan.get(), DUNEPdfs, xsec, osc);

    //###############################################################################################################################
    //Perform reweight and print total integral
    //###############################################################################################################################
    //Make oscillation channel breakdown
    // Initialise output file
    std::ofstream outFile("TestNewSampleOut.txt");

    writeToBoth(outFile, "======");
    writeToBoth(outFile, "Total Event Rates");
    writeToBoth(outFile, "======");
    for(auto Sample : DUNEPdfs){
        Sample->reweight();
        
        std::string EventRateString = fmt::format("{:.2f}", Sample->get1DHist()->Integral());
        writeToBoth(outFile, "Event rate for " + Sample->GetName() + ": " + EventRateString);
    }

    writeToBoth(outFile, "======");
    writeToBoth(outFile, "Event Rates by osc channel");
    writeToBoth(outFile, "======");

    for(auto Sample : DUNEPdfs) {
        int nOscChannels = Sample->getNMCSamples();
        for (int iOscChan=0;iOscChan<nOscChannels;iOscChan++) {
            std::vector< std::vector<double> > SelectionVec;

            std::vector<double> SelecChannel(3);
            SelecChannel[0] = Sample->ReturnKinematicParameterFromString("OscChannel");
            SelecChannel[1] = iOscChan;
            SelecChannel[2] = iOscChan+1;
            SelectionVec.push_back(SelecChannel);
            
            TH1* Hist = Sample->get1DVarHist("TrueNeutrinoEnergy",SelectionVec);
            writeToBoth(outFile, Sample->GetName() + " " + Sample->getFlavourName(iOscChan) + ": " + std::to_string(Hist->Integral()));
        }

        TH1* Hist = Sample->get1DVarHist("TrueNeutrinoEnergy");
        writeToBoth(outFile, Sample->GetName() + " " + std::to_string(Hist->Integral()));
    }

    //###############################################################################################################################
    //Make interaction channel breakdown
    writeToBoth(outFile, "======");
    writeToBoth(outFile, "Event Rates by interaction mode");
    writeToBoth(outFile, "======");

    for(auto Sample : DUNEPdfs) {
      int nModeChannels = Sample->GetMaCh3Modes()->GetNModes();
      for (int iModeChan=0;iModeChan<nModeChannels;iModeChan++) {
	std::vector< std::vector<double> > SelectionVec;
	
	std::vector<double> SelecChannel(3);
	SelecChannel[0] = Sample->ReturnKinematicParameterFromString("Mode");
	SelecChannel[1] = iModeChan;
	SelecChannel[2] = iModeChan+1;
	SelectionVec.push_back(SelecChannel);
	
	TH1* Hist = Sample->get1DVarHist("TrueNeutrinoEnergy",SelectionVec);
	writeToBoth(outFile, Sample->GetName() + " " + Sample->GetMaCh3Modes()->GetMaCh3ModeName(iModeChan) + ": " + std::to_string(Hist->Integral()));
      }
      
      TH1* Hist = Sample->get1DVarHist("TrueNeutrinoEnergy");
      writeToBoth(outFile, Sample->GetName() + ": " + std::to_string(Hist->Integral()));
    }

    // Do you want to gener
    if(skip_checks){
        return 0;
    }
    std::string results_file = "tests/test_results/EventRatesV2TDR.txt";
    if(CheckNodeExists(fitMan->raw()["Tests"]["TestResultsFile"])){
      results_file = fitMan->raw()["Tests"]["TestResultsFile"].as<std::string>();
    }

    // Okay now we've written we need to compare 
    bool TheSame = CompareTwoFiles(results_file, "TestNewSampleOut.txt");

    // Are we right?
    if(!TheSame) {
        MACH3LOG_CRITICAL("Different event rates detected");
        throw MaCh3Exception(__FILE__ , __LINE__ );
    } else {
        MACH3LOG_INFO("Event rates match");
    }

    return 0;
}
