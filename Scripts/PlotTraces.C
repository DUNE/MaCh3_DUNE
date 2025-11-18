#include <TString.h>
#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>

void PlotTraces(std::string FileName) {

  TFile* File = new TFile(FileName.c_str(),"READ");
  TDirectory* TraceDir = (TDirectory*)File->Get("Trace");
  TIter next(TraceDir->GetListOfKeys());
  gStyle->SetOptStat(0);

  TString HistName;
  TH1D* Hist;

  TCanvas* Canvas = new TCanvas("Canv","");
  Canvas->Print("TracesNDlong.pdf[");

  TObject* obj;
  while (obj = next()) {
    HistName = obj->GetName();
    std::cout << HistName << std::endl;
	Hist = (TH1D*)TraceDir->Get(HistName);
	Hist->SetBinContent(Hist->GetXaxis()->GetNbins(), Hist->GetBinContent(Hist->GetXaxis()->GetNbins()-1));
	std::cout << "Min = " << Hist->GetMinimum() << std::endl;
	std::cout << "Max = " << Hist->GetMaximum() << std::endl;
	//Hist->GetYaxis()->SetRangeUser(0.95*Hist->GetMinimum(), 1.05*Hist->GetMaximum());

    Hist->Draw();
    //Hist->GetXaxis()->SetRangeUser(0, 2.7e6);
    Canvas->Print("TracesNDlong.pdf");
  }
  Canvas->Print("TracesNDlong.pdf]");
  
  /*
  for (int i=0;i<1;i++) {
    std::cout << KeyList[i] << std::endl;
  }
  */
}

