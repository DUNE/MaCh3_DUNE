#include "MonolithSplineHandlerDUNE.h"
#include "TFile.h"

MonolithSplineHandlerDUNE::MonolithSplineHandlerDUNE(std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>> initParams) :
    SMonolith(initParams.first, initParams.second)
    { }


MonolithSplineHandlerDUNE::MonolithSplineHandlerDUNE(ParameterHandlerGeneric* xsec) :
    MonolithSplineHandlerDUNE(GetInitParamsFromConfig(xsec))
    { }

std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>> MonolithSplineHandlerDUNE:: GetInitParamsFromConfig(ParameterHandlerGeneric* xsec) {
    YAML::Node ConfigCurrent = xsec->GetConfig();
//     //Dump the whole config
    std::cout << ConfigCurrent << std::endl;

    // std::string spline_filename = xsec->Get<std::string>("SplineFile", __FILE__, __LINE__);
    // // For now, we will just initialize from file
    // MonolithSplineHandlerDUNE tempHandler(xsec);
    // tempHandler.InitFromFile(spline_filename);
    // // Here we would extract the init parameters from the loaded splines
    // // This is a placeholder; actual implementation would depend on how splines are stored
    // std::pair<std::vector<std::vector<TResponseFunction_red*> >, std::vector<RespFuncType>> initParams;
    // return initParams;
}

MonolithSplineHandlerDUNE::~MonolithSplineHandlerDUNE() {
}

void MonolithSplineHandlerDUNE::InitFromFile(std::string &spline_filename) {
    // Implementation for initializing from file goes here
    TFile* file = TFile::Open(spline_filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        MACH3LOG_ERROR("Failed to open spline file: {}", spline_filename);
        throw MaCh3Exception(__FILE__, __LINE__);
    }
    const std::string treeName = "SystWeights";

    TTree *tree;
    file->GetObject(treeName.c_str(), tree);
    if (!tree) {
        MACH3LOG_ERROR("Spline tree not found in file: {}", spline_filename);
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    TObjArray* list = tree->GetListOfBranches();
    for (int i = 0; i < list->GetEntries(); ++i) {
        TBranch* branch = (TBranch*)list->At(i);
        std::string splineName = branch->GetName();
        MACH3LOG_INFO("Found spline: {}", splineName);
        // Further processing to load the spline data can be added here
    }
  
}