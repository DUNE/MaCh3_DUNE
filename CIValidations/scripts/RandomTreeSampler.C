#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TRandom3.h>
#include <iostream>
#include <set>

/*
HW: Script that grabs all CAFs in a folder, randomly samples entries from each of them and makes a "Franken"-caf file.

Run:
root -l -b -q 'RandomTreeSampler.C("path/to/input/dir", "output_file.root")'
*/

void RandomTreeSampler(const char* inputDir, const char* outputFile = "sampled_output.root") {
    // Create a random number generator
    TRandom3 randGen(0); // Seed with 0 (uses UUID for unique seed)

    // Prepare the output file
    TFile* outFile = new TFile(outputFile, "RECREATE");
    TTree* outTree = nullptr;
    TH1D* normHist = nullptr;
    bool normCopied = false;

    // Get list of files in the input directory
    void* dirp = gSystem->OpenDirectory(inputDir);
    if (!dirp) {
        std::cerr << "Error opening directory: " << inputDir << std::endl;
        return;
    }

    const char* entry;
    while ((entry = gSystem->GetDirEntry(dirp))) {
        TString filename = entry;
        if (!filename.EndsWith(".root")) continue;

        TString fullPath = TString(inputDir) + "/" + filename;

        // Open the input file
        TFile* inFile = TFile::Open(fullPath);
        if (!inFile || inFile->IsZombie()) {
            std::cerr << "Error opening file: " << fullPath << std::endl;
            continue;
        }

        // Get the input tree
        TTree* inTree = (TTree*)inFile->Get("caf");
        if (!inTree) {
            std::cerr << "No 'caf' tree in file: " << fullPath << std::endl;
            inFile->Close();
            continue;
        }

        // Copy the "norm" histogram from the first valid file
        if (!normCopied) {
            normHist = (TH1D*)inFile->Get("norm");
            if (normHist) {
                outFile->cd();
                TH1D* histCopy = (TH1D*)normHist->Clone("norm");
                histCopy->Write();
                normCopied = true;
                std::cout << "Copied norm histogram from " << filename << std::endl;
            } else {
                std::cerr << "No 'norm' histogram in file: " << fullPath << std::endl;
            }
        }

        Long64_t max_entries = 500;
        Long64_t nEntries = inTree->GetEntries();
        Long64_t total_entries = std::min(max_entries, nEntries);

        // Create the output tree on first successful file
        if (!outTree) {
            outTree = inTree->CloneTree(0); // Clone structure but don't copy any entries
            outTree->SetName("caf");
            std::cout << "Created output tree from first file: " << fullPath << std::endl;
        }

        // IMPORTANT: Connect the input tree to the output tree's branch addresses
        outTree->CopyAddresses(inTree);

        // Select random entries
        std::set<Long64_t> selectedEntries;
        while (selectedEntries.size() < total_entries) {
            Long64_t entryNum = randGen.Integer(nEntries);
            selectedEntries.insert(entryNum);
        }

        // Copy the selected entries
        for (Long64_t entryNum : selectedEntries) {
            inTree->GetEntry(entryNum);
            outTree->Fill();
        }

        std::cout << "Copied " << total_entries << " entries from " << filename << std::endl;
        inFile->Close();
    }

    // Write and close the output file
    if (outTree) {
        outFile->cd();
        outTree->Write();
        std::cout << "Successfully created output file: " << outputFile
                  << " with " << outTree->GetEntries() << " entries" << std::endl;
        
        // If we never found a norm histogram, create an empty one
        if (!normCopied) {
            std::cerr << "Warning: No 'norm' histogram found in any input file. Creating empty one." << std::endl;
            TH1D* emptyNorm = new TH1D("norm", "Normalization histogram", 1, 0, 1);
            emptyNorm->Write();
            delete emptyNorm;
        }
    } else {
        std::cerr << "No valid trees were processed!" << std::endl;
    }

    outFile->Close();
}