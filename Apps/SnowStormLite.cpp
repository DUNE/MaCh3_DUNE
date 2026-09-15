#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <filesystem>

#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"


// ============================================================
// Event information
//
// Add the quantities from your CAF that are needed by your
// functional parameters.
// ============================================================

struct EventData
{
    // ============================================================
    // Reconstructed energies
    //
    // These are read directly from the input TTree. Depending on
    // whether this is a reco-numu or reco-nue file, different CAF
    // branches will be connected to these variables.
    // ============================================================

    double recoEnergy;
    double recoHadEnergy;
    double recoLepEnergy;
    double recoProtonEnergy;
    double recoNeutronEnergy;
    double recoPipEnergy;
    double recoPimEnergy;
    double recoPi0Energy;
    double recoHadEnergySum;


    // ============================================================
    // Derived quantities
    //
    // These are NOT read from branches. They are calculated after
    // reading each event.
    // ============================================================

    double recoHadEnergySqrt;
    double recoLepEnergySqrt;
    double recoNeutronEnergySqrt;
    double recoPi0EnergySqrt;
    double recoHadEnergySumSqrt;
    


    // ============================================================
    // True particle energies
    //
    // These are read directly from the input TTree.
    // ============================================================

    double trueProtonEnergy;
    double truePipEnergy;
    double truePimEnergy;

    // Generic true charged-lepton energy. This corresponds to
    // rw_LepE in the functions you provided.
    double trueLepEnergy;

    double trueNeutronEnergy;
    double truePi0Energy;


    // ============================================================
    // Interaction information
    //
    // Used to identify true CC numu/nue events.
    // ============================================================

    int isCC;
    int nuPDG;
};


// ============================================================
// Functional parameter
//
// Each parameter has:
//   - a name
//   - a Gaussian distribution
//   - a function giving its contribution to reco energy
// ============================================================

struct FunctionalParameter
{
    std::string name;

    double mean;
    double sigma;

    // Returns the change in reconstructed energy
    std::function<double(const EventData&, double)> response;
};

enum class RecoChannel
{
    Numu,
    Nue
};

bool IsTrueCCNumu(const EventData& e)
{
    return e.isCC == 1 && e.nuPDG == 14;
}

bool IsTrueCCNue(const EventData& e)
{
    return e.isCC == 1 && e.nuPDG == 12;
}


// ============================================================
// Main
// Usage: SnowStormLite input.root
// ============================================================

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr
            << "Usage: "
            << argv[0]
            << " input.root"
            << std::endl;

        return 1;
    }


    // ========================================================
    // Input/output filenames
    // ========================================================

    const std::string inputFilename = argv[1];
  
    RecoChannel recoChannel;
    
    if (inputFilename.find("numuselec") != std::string::npos)
    {
        recoChannel = RecoChannel::Numu;
    }
    else if (inputFilename.find("nueselec") != std::string::npos)
    {
        recoChannel = RecoChannel::Nue;
    }
    else
    {
        std::cerr
            << "ERROR: Could not determine reconstruction channel "
            << "from filename.\n"
            << "Expected filename to contain either 'numu' or 'nue'."
            << std::endl;
    
        return 1;
    }

    std::filesystem::path inputPath(inputFilename);

    const std::filesystem::path outputDirectory = "/work4/ppd/scarf1407/DUNE_2026_FD_TDR_CAFs_SnowStormLite_TotalEScaleSqrt";

    const std::string outputFilename =
        (outputDirectory /
         (inputPath.stem().string() +
          "_SnowStormLite_TotalEScaleSqrt.root")).string();


    // ========================================================
    // Open input CAF
    // ========================================================

    TFile inputFile(
        inputFilename.c_str(),
        "READ"
    );

    if (inputFile.IsZombie())
    {
        std::cerr
            << "ERROR: Could not open "
            << inputFilename
            << std::endl;

        return 1;
    }


    // CHANGE THIS TO THE ACTUAL TREE NAME
    TTree* inputTree =
        dynamic_cast<TTree*>(
            inputFile.Get("caf")
        );

    if (!inputTree)
    {
        std::cerr
            << "ERROR: Could not find input TTree"
            << std::endl;

        return 1;
    }


    // ========================================================
    // Set up input branches
    // ========================================================

    EventData event{};

    if (recoChannel == RecoChannel::Numu)
    {

      inputTree->SetBranchAddress(
          "Ev_reco_numu",
          &event.recoEnergy
      );

      inputTree->SetBranchAddress(
          "RecoHadEnNumu",
          &event.recoHadEnergy
      );

      inputTree->SetBranchAddress(
          "RecoLepEnNumu",
          &event.recoLepEnergy
      );
    }
    else {

      inputTree->SetBranchAddress(
          "Ev_reco_nue",
          &event.recoEnergy
      );

      inputTree->SetBranchAddress(
          "RecoHadEnNue",
          &event.recoHadEnergy
      );

      inputTree->SetBranchAddress(
          "RecoLepEnNue",
          &event.recoLepEnergy
      );
    }

    inputTree->SetBranchAddress(
        "eRecoP",
        &event.recoProtonEnergy
    );

    inputTree->SetBranchAddress(
        "eRecoN",
        &event.recoNeutronEnergy
    );

    inputTree->SetBranchAddress(
        "eRecoPip",
        &event.recoPipEnergy
    );

    inputTree->SetBranchAddress(
        "eRecoPim",
        &event.recoPimEnergy
    );

    inputTree->SetBranchAddress(
        "eRecoPi0",
        &event.recoPi0Energy
    );

    inputTree->SetBranchAddress(
        "eP",
        &event.trueProtonEnergy
    );

    inputTree->SetBranchAddress(
        "ePip",
        &event.truePipEnergy
    );

    inputTree->SetBranchAddress(
        "ePim",
        &event.truePimEnergy
    );

    inputTree->SetBranchAddress(
        "LepE",
        &event.trueLepEnergy
    );

    inputTree->SetBranchAddress(
        "eN",
        &event.trueNeutronEnergy
    );

    inputTree->SetBranchAddress(
        "ePi0",
        &event.truePi0Energy
   );

    inputTree->SetBranchAddress(
        "isCC",
        &event.isCC
   );

    inputTree->SetBranchAddress(
        "nuPDG",
        &event.isCC
    );


    // ========================================================
    // Define functional parameters
    //
    // Add all of your systematic functions here.
    // ========================================================

    std::vector<FunctionalParameter> parameters;
    std::vector<FunctionalParameter> parametersDummy;

    // ============================================================
    // TOTAL ENERGY SCALE
    //
    // TotalEScale + TotalEScaleNotCCNumu
    // ============================================================
    
    parameters.push_back({
        "TotalEScale",
        0.0, // mean
        0.02, // sigma
    
        [recoChannel](const EventData& e, double par)
        {
            double deltaE = par * e.recoHadEnergy;
    
            // TotalEScaleNotCCNumu component
            if (!IsTrueCCNumu(e) && recoChannel == RecoChannel::Nue)
                deltaE += par * e.recoLepEnergy;
    
            return deltaE;
        }
    });
    
    
    // ============================================================
    // TOTAL ENERGY SCALE SQRT
    //
    // TotalEScaleSqrt + TotalEScaleSqrtNotCCNumu
    // ============================================================
    
    parametersDummy.push_back({
        "TotalEScaleSqrt",
        0.0,
        0.01,
    
        [recoChannel](const EventData& e, double par)
        {
            double deltaE =
                par * e.recoHadEnergy * e.recoHadEnergySqrt;
    
            if (!IsTrueCCNumu(e) && recoChannel == RecoChannel::Nue)
            {
                deltaE +=
                    par * e.recoLepEnergy * e.recoLepEnergySqrt;
            }
    
            return deltaE;
        }
    });
    
    
    // ============================================================
    // TOTAL ENERGY SCALE INVERSE SQRT
    //
    // TotalEScaleInvSqrt + TotalEScaleInvSqrtNotCCNumu
    // ============================================================
    
    parameters.push_back({
        "TotalEScaleInvSqrt",
        0.0,
        0.02,
    
        [recoChannel](const EventData& e, double par)
        {
            double deltaE = par * e.recoHadEnergySqrt;
    
            if (!IsTrueCCNumu(e) && recoChannel == RecoChannel::Nue)
                deltaE += par * e.recoLepEnergySqrt;
    
            return deltaE;
        }
    });
    
    
    // ============================================================
    // HADRONIC ENERGY SCALE
    // ============================================================
    
    parameters.push_back({
        "HadEScale",
        0.0,
        0.05,
    
        [recoChannel](const EventData& e, double par)
        {
            return par * e.recoHadEnergySum;
        }
    });
    
    
    // ============================================================
    // HADRONIC ENERGY SCALE SQRT
    // ============================================================
    
    parameters.push_back({
        "HadEScaleSqrt",
        0.0,
        0.05,
    
        [recoChannel](const EventData& e, double par)
        {
            return par *
                   e.recoHadEnergySum *
                   e.recoHadEnergySumSqrt;
        }
    });
    
    
    // ============================================================
    // HADRONIC ENERGY SCALE INVERSE SQRT
    // ============================================================
    
    parameters.push_back({
        "HadEScaleInvSqrt",
        0.0,
        0.05,
    
        [recoChannel](const EventData& e, double par)
        {
            return par * e.recoHadEnergySumSqrt;
        }
    });
    
    
    // ============================================================
    // MUON ENERGY SCALE
    //
    // Original function calls TotalEScaleNotCCNumu, but this is an
    // INDEPENDENT parameter, so the implementation is written directly.
    // ============================================================
    
    parameters.push_back({
        "MuEScale",
        0.0,
        0.02,
    
        [recoChannel](const EventData& e, double par)
        {
            if (IsTrueCCNumu(e))
                return par * e.recoLepEnergy;
           
            return 0.0;
        }
    });
    
    
    // ============================================================
    // MUON ENERGY SCALE SQRT
    // ============================================================
    
    parameters.push_back({
        "MuEScaleSqrt",
        0.0,
        0.005,
    
        [recoChannel](const EventData& e, double par)
        {
            if (IsTrueCCNumu(e))
                return par *
                   e.recoLepEnergy *
                   e.recoLepEnergySqrt;
           
            return 0.0;
        }
    });
    
    
    // ============================================================
    // MUON ENERGY SCALE INVERSE SQRT
    // ============================================================
    
    parameters.push_back({
        "MuEScaleInvSqrt",
        0.0,
        0.02,
    
        [recoChannel](const EventData& e, double par)
        {
            if (IsTrueCCNumu(e))
                return par *
                   e.recoLepEnergySqrt;
           
            return 0.0;
        }
    });
    
    
    // ============================================================
    // NEUTRON ENERGY SCALE
    // ============================================================
    
    parameters.push_back({
        "NEScale",
        0.0,
        0.2,
    
        [recoChannel](const EventData& e, double par)
        {
            return par * e.recoNeutronEnergy;
        }
    });
    
    
    // ============================================================
    // NEUTRON ENERGY SCALE SQRT
    // ============================================================
    
    parameters.push_back({
        "NEScaleSqrt",
        0.0,
        0.3,
    
        [recoChannel](const EventData& e, double par)
        {
            return par *
                   e.recoNeutronEnergy *
                   e.recoNeutronEnergySqrt;
        }
    });
    
    
    // ============================================================
    // NEUTRON ENERGY SCALE INVERSE SQRT
    // ============================================================
    
    parameters.push_back({
        "NEScaleInvSqrt",
        0.0,
        0.3,
    
        [recoChannel](const EventData& e, double par)
        {
            return par * e.recoNeutronEnergySqrt;
        }
    });
    
    
    // ============================================================
    // ELECTROMAGNETIC ENERGY SCALE
    //
    // EMEScale + EMEScaleCCNue
    // ============================================================
    
    parameters.push_back({
        "EMEScale",
        0.0,
        0.025,
    
        [recoChannel](const EventData& e, double par)
        {
            double deltaE = par * e.recoPi0Energy;
    
            // EMEScaleCCNue component
            if (IsTrueCCNue(e) && recoChannel == RecoChannel::Nue)
                deltaE += par * e.recoLepEnergy;
    
            return deltaE;
        }
    });
    
    
    // ============================================================
    // ELECTROMAGNETIC ENERGY SCALE SQRT
    //
    // EMEScaleSqrt + EMEScaleSqrtCCNue
    // ============================================================
    
    parameters.push_back({
        "EMEScaleSqrt",
        0.0,
        0.025,
    
        [recoChannel](const EventData& e, double par)
        {
            double deltaE =
                par *
                e.recoPi0Energy *
                e.recoPi0EnergySqrt;
    
            if (IsTrueCCNue(e) && recoChannel == RecoChannel::Nue)
            {
                deltaE +=
                    par *
                    e.recoLepEnergy *
                    e.recoLepEnergySqrt;
            }
    
            return deltaE;
        }
    });
    
    
    // ============================================================
    // ELECTROMAGNETIC ENERGY SCALE INVERSE SQRT
    //
    // EMEScaleInvSqrt + EMEScaleInvSqrtCCNue
    // ============================================================
    
    parameters.push_back({
        "EMEScaleInvSqrt",
        0.0,
        0.025,
    
        [recoChannel](const EventData& e, double par)
        {
            double deltaE = par * e.recoPi0EnergySqrt;
    
            if (IsTrueCCNue(e) && recoChannel == RecoChannel::Nue)
                deltaE += par * e.recoLepEnergySqrt;
    
            return deltaE;
        }
    });
    
    
    // ============================================================
    // HADRONIC RESOLUTION
    // ============================================================
    
    parameters.push_back({
        "HadRes",
        0.0,
        0.02,
    
        [recoChannel](const EventData& e, double par)
        {
            const double trueHadEnergy =
                e.trueProtonEnergy +
                e.truePipEnergy +
                e.truePimEnergy;
    
            return par *
                   (trueHadEnergy - e.recoHadEnergySum);
        }
    });
    
    
    // ============================================================
    // MUON RESOLUTION
    // ============================================================
    
    parameters.push_back({
        "MuRes",
        0.0,
        0.02,
    
        [recoChannel](const EventData& e, double par)
        {
            if (IsTrueCCNumu(e) && recoChannel == RecoChannel::Numu)
                 return par *
                        (e.trueLepEnergy - e.recoLepEnergy);

	    return 0.0;
        }
    });
    
    
    // ============================================================
    // NEUTRON RESOLUTION
    // ============================================================
    
    parameters.push_back({
        "NRes",
        0.0,
        0.1,
    
        [recoChannel](const EventData& e, double par)
        {
            return par *
                   (e.trueNeutronEnergy -
                    e.recoNeutronEnergy);
        }
    });
    
    
    // ============================================================
    // ELECTROMAGNETIC RESOLUTION
    //
    // EMRes + EMResCCNue
    // ============================================================
    
    parameters.push_back({
        "EMRes",
        0.0,
        0.02,
    
        [recoChannel](const EventData& e, double par)
        {
            double deltaE =
                par *
                (e.truePi0Energy - e.recoPi0Energy);
    
            // EMResCCNue calls MuRes in the original implementation.
            // Written out explicitly here, using the SAME EMRes parameter.
            if (IsTrueCCNue(e) && recoChannel == RecoChannel::Nue)
            {
                deltaE +=
                    par *
                    (e.trueLepEnergy - e.recoLepEnergy);
            }
    
            return deltaE;
        }
    });

    // ========================================================
    // Random number generator
    // ========================================================

    // Use a fixed seed for reproducibility while developing.
    // Change to 0 if you want ROOT to generate a new seed each run.
    TRandom3 random(12345);


    // ========================================================
    // Create output file
    // ========================================================

    TFile outputFile(
        outputFilename.c_str(),
        "RECREATE"
    );


    // Clone the structure of the original tree.
    //
    // All original CAF information will therefore be retained.
    TTree* outputTree =
        inputTree->CloneTree(0);

    outputTree->SetName(
        inputTree->GetName()
    );


    // ========================================================
    // New output branches
    // ========================================================

    // The final varied reconstructed energy
    double variedRecoEnergy = 0.0;

    outputTree->Branch(
        "rw_erec_snowstorm",
        &variedRecoEnergy,
        "rw_erec_snowstorm/D"
    );


    // --------------------------------------------------------
    // Create ONE branch per functional parameter
    // --------------------------------------------------------

    // This vector stores the current event's parameter values.
    //
    // IMPORTANT:
    // We resize it BEFORE creating branches so that the memory
    // addresses used by ROOT remain stable.
    std::vector<double> parameterThrows(
        parametersDummy.size(),
        0.0
    );


    for (std::size_t i = 0;
         i < parametersDummy.size();
         ++i)
    {
        const std::string branchName =
            "snowstorm_" + parametersDummy[i].name;

        const std::string leafName =
            branchName + "/D";

        outputTree->Branch(
            branchName.c_str(),
            &parameterThrows[i],
            leafName.c_str()
        );
    }


    // ========================================================
    // Event loop
    // ========================================================

    const Long64_t nEvents =
        inputTree->GetEntries();

    std::cout
        << "Processing "
        << nEvents
        << " events..."
        << std::endl;


    for (Long64_t iEvent = 0;
         iEvent < nEvents;
         ++iEvent)
    {
        // Read original CAF event
        inputTree->GetEntry(iEvent);

        event.recoHadEnergySqrt = std::sqrt(event.recoHadEnergy);

        event.recoLepEnergySqrt = std::sqrt(event.recoLepEnergy);

        event.recoNeutronEnergySqrt = std::sqrt(event.recoNeutronEnergy);

        event.recoPi0EnergySqrt = std::sqrt(event.recoPi0Energy);
        
        event.recoHadEnergySum = event.recoProtonEnergy + event.recoPipEnergy + event.recoPimEnergy;

        event.recoHadEnergySumSqrt = std::sqrt(event.recoHadEnergySum);


        // Start with nominal reconstructed energy
        variedRecoEnergy =
            event.recoEnergy;


        // ----------------------------------------------------
        // Throw every systematic parameter for this event
	// Edit this loop to select a subset of the systematic parameters
        // ----------------------------------------------------
        for (std::size_t iPar = 0;
             iPar < parametersDummy.size();
             ++iPar)
        {
            const FunctionalParameter& parameter =
                parametersDummy[iPar];


            // Draw the SnowStorm value
            parameterThrows[iPar] =
                random.Gaus(
                    parameter.mean,
                    parameter.sigma
                );


            // Apply the functional response if appropriate
            variedRecoEnergy +=
                parameter.response(
                    event,
                    parameterThrows[iPar]
                    );
        }


        // ----------------------------------------------------
        // Write:
        //
        //   - all original CAF information
        //   - all parameter throws
        //   - varied reconstructed energy
        // ----------------------------------------------------

        outputTree->Fill();


        // Progress
        if (iEvent % 100000 == 0)
        {
            std::cout
                << "Processed "
                << iEvent
                << " / "
                << nEvents
                << std::endl;
        }
    }


    // ========================================================
    // Write output
    // ========================================================

    outputFile.cd();

    outputTree->Write();

    outputFile.Close();
    inputFile.Close();


    std::cout
        << "\nFinished successfully.\n"
        << "Output: "
        << outputFilename
        << std::endl;


    return 0;
}
