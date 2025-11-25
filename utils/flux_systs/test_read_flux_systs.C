{
  gROOT->ProcessLine("#define VERBOSE_MODE");
  gROOT->ProcessLine(".L OffAxisFluxUncertaintyHelper.cxx+");

  gROOT->ProcessLine(R"(
    OffAxisFluxUncertaintyHelper flux_helper;
    flux_helper.Initialize("flux_variations_FD_and_PRISM_2023.root", true);


    for(size_t i = 0; i < flux_helper.GetNFocussingParams(); ++i){
      std::cout << i << ": " << flux_helper.GetFocussingParamName(i) << std::endl;
    }
  )");

  
}