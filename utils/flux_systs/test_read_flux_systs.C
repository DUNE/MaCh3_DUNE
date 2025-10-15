{
  gROOT->ProcessLine(".include ../../");
  gROOT->ProcessLine(".L OffAxisFluxUncertaintyHelper.cxx+");

  gROOT->ProcessLine(R"(
    OffAxisFluxUncertaintyHelper flux_helper;
    flux_helper.Initialize("flux_variations_FD_and_PRISM_2023.root", true);


    for(size_t i = 0; i < flux_helper.GetNFocussingParams(); ++i){
      std::cout << i << ": " << flux_helper.GetFocussingParamName(i) 
                << std::endl;
    }

    std::cout << "Focussing Off Axis Binning: [ ";
    for(auto be : flux_helper.GetFluxFocussingOffAxisBinning()){
      std::cout << be << ", ";
    }
    std::cout << "]" << std::endl;
    std::cout << "Have " << flux_helper.GetNHadProdPCAComponents() 
              << " hadron production PCA components." << std::endl;

    for(int nucfg = 0; nucfg < OffAxisFluxUncertaintyHelper::kFD_numu_numode; ++nucfg){
      std::cout << "HadProd Off Axis Binning for nucfg: " << nucfg  << ", [ ";
      for(auto be : flux_helper.GetFluxFocussingOffAxisBinning()){
        std::cout << be << ", ";
      }
      std::cout << "]" << std::endl;
    }
  )");
}