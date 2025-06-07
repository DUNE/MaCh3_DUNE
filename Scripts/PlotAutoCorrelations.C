#include <TString.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>

void PlotAutoCorrelations(std::string FileName) {

  TFile* File = new TFile(FileName.c_str(),"READ");

  gStyle->SetOptStat(0);
  TDirectory* AutoCorrDir2 = (TDirectory*)File->Get("Trace");
  TIter next2(AutoCorrDir2->GetListOfKeys());
  std::vector<TString> Titles;
  TString HistName2;
  TObject* obj2;
  while (obj2 = next2()) {
    HistName2 = obj2->GetName();
	Titles.push_back(HistName2);
  }

  TDirectory* AutoCorrDir = (TDirectory*)File->Get("Auto_corr");
  TIter next(AutoCorrDir->GetListOfKeys());

  TString HistName;
  TH1D* Hist;

  TCanvas* Canvas = new TCanvas("Canv","");
  Canvas->Print("AutoCorrNDlong.pdf[");

  TObject* obj;
  int i = 0;
  while (obj = next()) {
    HistName = obj->GetName();
    std::cout << Titles[i] << std::endl;
    Hist = (TH1D*)AutoCorrDir->Get(HistName);
    //Hist->GetXaxis()->SetRangeUser(0, 2.7e6;
    Hist->SetTitle(Titles[i]);
    Hist->Draw();
    Canvas->Print("AutoCorrNDlong.pdf");
	i++;
  }
  Canvas->Print("AutoCorrNDlong.pdf]");
  
  /*
  for (int i=0;i<1;i++) {
    std::cout << KeyList[i] << std::endl;
  }
  */
}

