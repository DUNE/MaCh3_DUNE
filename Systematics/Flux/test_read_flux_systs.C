{
  gROOT->ProcessLine(".include ../../");
  gROOT->ProcessLine(".L OffAxisFluxUncertaintyHelper.cxx+");

  gROOT->ProcessLine(R"(
    OffAxisFluxUncertaintyHelper flux_helper;
    flux_helper.Initialize("flux_variations_FD_and_PRISM_2023.root", true);


    for(size_t i = 0; i < flux_helper.GetNFocussingParams(); ++i){
      std::cout << i << ": " << flux_helper.GetFocussingParamName(i) 
                << std::endl;

      std::cout << "    Binning: [ ";
      for(auto be : flux_helper.GetFluxFocussingOffAxisBinning()){
        std::cout << be << ", ";
      }
      std::cout << "]" << std::endl;
    }

    std::cout << "Have " << flux_helper.GetNHadProdPCAComponents() 
              << " hadron production PCA components." << std::endl;

    for(int nucfg = 0; nucfg < OffAxisFluxUncertaintyHelper::kFD_numu_numode; ++nucfg){
      std::cout << "HadProd Off Axis Binning for nucfg: " << nucfg  << ", [ ";
      for(auto be : flux_helper.GetFluxHadProdOffAxisBinning(nucfg)){
        std::cout << be << ", ";
      }
      std::cout << "]" << std::endl;
    }

    TFile fout("flux_1sigma_weights.root","RECREATE");

    for(size_t i = 0; i < flux_helper.GetNFocussingParams(); ++i){
      auto dir = fout.mkdir(flux_helper.GetFocussingParamName(i).c_str());

      for(double oa_m : std::vector<double>{0.0, 12.0}){
        TH1D hist(("ND_Numu_" + std::to_string(int(oa_m)) + "m").c_str(),
          (flux_helper.GetFocussingParamName(i) +";Enu;Weight").c_str(), 
          1000,0,5);

        int nucfg = 0;
        for(int bi = 0; bi < 1000; ++bi){
          int syst_bi = flux_helper.GetFocussingBin(nucfg,
                          hist.GetXaxis()->GetBinCenter(bi+1), oa_m);
          hist.SetBinContent(bi+1, flux_helper.GetFluxFocussingWeight(i, 1, 
                                      nucfg, syst_bi));
        }
        dir->WriteTObject(&hist, hist.GetName());
      }
    }

    for(size_t i = 0; i < flux_helper.GetNHadProdPCAComponents(); ++i){
      auto dir = fout.mkdir(("HadronProduction_pca_" + std::to_string(i)).c_str());

      for(double oa_m : std::vector<double>{0.0, 12.0}){
        TH1D hist(("ND_Numu_" + std::to_string(int(oa_m)) + "m").c_str(),
          ("HadronProduction_pca_" + std::to_string(i)).c_str(), 
          1000,0,5);

        int nucfg = 0;
        for(int bi = 0; bi < 1000; ++bi){
          int syst_bi = flux_helper.GetHadProdBin(nucfg,
                          hist.GetXaxis()->GetBinCenter(bi+1), oa_m);
          hist.SetBinContent(bi+1, flux_helper.GetFluxHadProdWeight(i, 1, 
                                      nucfg, syst_bi));
        }
        dir->WriteTObject(&hist, hist.GetName());
      }
    }

    fout.Write();
    fout.Close();

  )");
}