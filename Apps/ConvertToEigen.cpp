#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/SampleHandlerAtm.h"
#include "Fitters/MaCh3Factory.h"

int main(int argc, char * argv[]) {

  M3::Utils::MaCh3Usage(argc, argv);
  auto FitManager = MaCh3ManagerFactory(argc, argv);

  //###############################################################################################################################
  //Create SampleHandlerFD objects
  
  ParameterHandlerGeneric* xsec = nullptr;

  //Create SampleHandlerBase objects
  auto [param_handler, DUNEPdfs] = MaCh3DuneFactory(FitManager);

  for (size_t iSample=0;iSample<DUNEPdfs.size();iSample++) {
    SampleHandlerAtm* Sample = dynamic_cast<SampleHandlerAtm*>(DUNEPdfs[iSample]);
    if (!Sample) {
      MACH3LOG_ERROR("Can only convert SampleHandlerAtm CAFs to Eigen matrix binary input file");
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    std::string FileName = "AtmSample.eig";
    Sample->TransferToEigen(FileName);
  }
}
