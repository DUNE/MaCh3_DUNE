#include <TSystem.h>
#include <TFile.h>
#include <TError.h>
#include <TString.h>
#include <iostream>
#include <vector>

/*
    HW: ROOT macro script, literally just loops over files in a directory and opens then, if thye're borken raises error. 
    Goes through full directory to ensure nothing is broken. Developed since CVMFS can have some issues with the CI
*/


void check_corrupt(const char* filename) {
    TFile* file = TFile::Open(filename);
    if (file && !file->IsZombie()) {
        std::cout << "\033[1;34m" << filename << "\033[0m \033[1;32mFINE\033[0m" << std::endl;
        file->Close();
        delete file;
    } else {
        throw std::runtime_error("Could not open file properly");
    }
}

void check_root_files(const char* directory) {
    int fail = 0;
    
    void* dirp = gSystem->OpenDirectory(directory);
    if (!dirp) {
        std::cerr << "Error opening directory: " << directory << std::endl;
        return;
    }
    
    const char* entry;
    std::vector<TString> root_files;
    
    // Collect all .root files
    while ((entry = gSystem->GetDirEntry(dirp))) {
        TString filename = entry;
        if (filename.EndsWith(".root")) {
            TString fullpath = TString::Format("%s/%s", directory, filename.Data());
            root_files.push_back(fullpath);
        }
    }
    gSystem->FreeDirectory(dirp);
    
    // Check each file
    for (const auto& filepath : root_files) {
        try {
            check_corrupt(filepath.Data());
        } catch (...) {
            fail++;
            std::cout << "\033[1;33m" << filepath << "\033[0m \033[1;31mCORRUPT!\033[0m" << std::endl;
        }
    }
    
    if (fail > 0) {
        Throw("\033[1;33m Warning, failed check on \033[0m \033[1;31m" << fail << " files\033[0m")
        // To match Python's behavior of raising an exception
    }
}

void root_file_corrupt(const char* directory = nullptr) {
    // Set error handling
    gErrorIgnoreLevel = kWarning;
    if (!directory) {
        std::cerr << "Usage: .root -l -b -q root_file_corrupt.C+\"<directory>\"" << std::endl;
        throw std::runtime_error("No directory provided");
    }
    
    check_root_files(directory);
}