#include "Samples/MaCh3DUNEFactory.h"
#include "Fitters/MaCh3Factory.h"
#include "Fitters/PredictiveThrower.h"

void RunPredictive(std::vector<SampleHandlerFD*> DUNEPdfs, ParameterHandlerGeneric* xsec, manager* PredictiveManager) {
  std::unique_ptr<PredictiveThrower> MaCh3Fitter = std::make_unique<PredictiveThrower>(PredictiveManager);
  MaCh3Fitter->AddSystObj(xsec);
  
  for (size_t i = 0; i < DUNEPdfs.size(); ++i) {
    MaCh3Fitter->AddSampleHandler(DUNEPdfs[i]);
  }
  
  MaCh3Fitter->ProduceToys();
}

int main(int argc, char * argv[]) {
  //###############################################################################################################################
  MaCh3Utils::MaCh3Usage(argc, argv);
  auto FitManager = MaCh3ManagerFactory(argc, argv);

  //###############################################################################################################################
  //Create SampleHandlerFD objects
   
  ParameterHandlerGeneric* xsec = nullptr;
  std::vector<SampleHandlerFD*> DUNEPdfs;
  
  MakeMaCh3DuneInstance(FitManager, DUNEPdfs, xsec);

  //###############################################################################################################################
  //Perform the predictive study

  if (FitManager->raw()["Predictive"]["LoopOverParameterSets"]) {
    std::vector<std::vector<std::string>> ParameterSets = Get<std::vector<std::vector<std::string>>>(FitManager->raw()["Predictive"]["LoopOverParameterSets"],__FILE__,__LINE__);
    
    for (size_t iParameterSet=0;iParameterSet<ParameterSets.size();++iParameterSet) {
      
      YAML::Node Config = YAML::Node(FitManager->raw());
      Config["General"]["OutputFile"] = Form("PredictiveOutput_ParameterSet_%i.root",iParameterSet);
      Config["Predictive"]["ThrowSinlgeParams"] = ParameterSets[iParameterSet];
      manager* PredictiveManager = new manager(Config);

      RunPredictive(DUNEPdfs, xsec, PredictiveManager);
    }
  } else {
    RunPredictive(DUNEPdfs, xsec, FitManager.get());
  }
  
  //###############################################################################################################################
}
