#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/SampleHandlerAtm.h"
#include "Samples/DUNEUtils.h"
#include "Fitters/MaCh3Factory.h"

int main(int argc, char * argv[]) {
  M3::Utils::MaCh3Usage(argc, argv);
  auto FitManager = MaCh3ManagerFactory(argc, argv);

  std::string HistogramFileName = Get<std::string>(FitManager->raw()["PlotSliceHistograms"]["InputHistogramFileName"],__FILE__,__LINE__);
  int SliceAxis = Get<int>(FitManager->raw()["PlotSliceHistograms"]["SliceAxis"],__FILE__,__LINE__);
  bool LogX = Get<bool>(FitManager->raw()["PlotSliceHistograms"]["LogX"],__FILE__,__LINE__);
  bool LogY = Get<bool>(FitManager->raw()["PlotSliceHistograms"]["LogY"],__FILE__,__LINE__);
  int nColumns = Get<int>(FitManager->raw()["PlotSliceHistograms"]["NumColumns"],__FILE__,__LINE__);

  std::unordered_map<std::string, std::string> FancyTitles;
  auto FancyTitlesNode = FitManager->raw()["PlotSliceHistograms"]["FancyTitles"];
  if (FancyTitlesNode) {
    for (const auto& entry : FancyTitlesNode) {
      FancyTitles[entry.first.as<std::string>()] =
	entry.second.as<std::string>();
    }
  }
  
  //###############################################################################################################################
  //Create SampleHandlerFD objects
  
  ParameterHandlerGeneric* xsec = nullptr;

  //Create SampleHandlerBase objects
  auto [param_handler, DUNEPdfs] = MaCh3DuneFactory(FitManager);

  TFile* HistogramFile = new TFile(HistogramFileName.c_str());
  if (!HistogramFile || HistogramFile->IsZombie()) {
    MACH3LOG_ERROR("Can't find Histogram file: {}",HistogramFileName);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  for (size_t iSampleHandler=0;iSampleHandler<DUNEPdfs.size();iSampleHandler++) {
    SampleHandlerAtm* Sample = dynamic_cast<SampleHandlerAtm*>(DUNEPdfs[iSampleHandler]);
    if (!Sample) {
      MACH3LOG_ERROR("Can only convert SampleHandlerAtm CAFs to Eigen matrix binary input file");
      throw MaCh3Exception(__FILE__, __LINE__);
    }

    auto BinningHandler = Sample->GetBinningHandler();

    for (int iSample=0;iSample<Sample->GetNSamples();iSample++) {
      std::string HistogramName = "h"+Sample->GetSampleTitle(iSample);
      TH1* Hist = HistogramFile->Get<TH1>(HistogramName.c_str());

      auto BinningVariables = Sample->GetSampleBinningVariables(iSample);
      for (size_t iVar=0;iVar<BinningVariables.size();iVar++) {
	if (auto it = FancyTitles.find(BinningVariables[iVar]); it != FancyTitles.end()) {
	  BinningVariables[iVar] = it->second;
	}
      }
      
      auto slices = Make2DSlices(Hist, *BinningHandler, iSample, BinningVariables, SliceAxis);
      auto canvas = Plot2DSlices(slices, HistogramName+"_Canvas", LogX, LogY, nColumns);
      canvas->Print((HistogramName+".pdf").c_str());
    }
  }
}
