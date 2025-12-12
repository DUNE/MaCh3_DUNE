#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/SampleHandlerAtm.h"
#include "Fitters/MaCh3Factory.h"

int main(int argc, char * argv[]) {

  MaCh3Utils::MaCh3Usage(argc, argv);
  auto fitMan = MaCh3ManagerFactory(argc, argv);

  //###############################################################################################################################
  //Create SampleHandlerFD objects
  
  ParameterHandlerGeneric* xsec = nullptr;
  std::vector<SampleHandlerFD*> DUNEPdfs;
  MakeMaCh3DuneInstance(fitMan, DUNEPdfs, xsec);

  for (size_t iSample=0;iSample<DUNEPdfs.size();iSample++) {
    SampleHandlerAtm* Sample = dynamic_cast<SampleHandlerAtm*>(DUNEPdfs[iSample]);
    if (!Sample) {
      MACH3LOG_ERROR("Can only convert SampleHandlerAtm CAFs to Eigen matrix binary input file");
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    std::string FileName = "AtmSample_"+Sample->GetTitle()+".eig";
    Sample->TransferToEigen(FileName);
  }
}
