void MakeOscBinning() 
{

  TFile *BinFile = new TFile("OscBinning.root", "NEW");
  std::vector<double> Binning;
  for (double i = 0.0; i <= 10.0; i+=0.01) 
  {
    Binning.push_back(i);
  }

  Binning.push_back(15.0);  
  Binning.push_back(20.0);  
  Binning.push_back(40.0);  
  Binning.push_back(100.0);  
  Binning.push_back(1e8);


  TH1D* BinningHist = new TH1D("OscBinning", "OscBinning", Binning.size()-1, Binning.data());
  BinningHist->Write();
  BinFile->Close();
}
	 
