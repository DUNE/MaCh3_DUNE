#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraphErrors.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMath.h"


#include <iostream>

void makePredictivePlots()
{

//   TFile* datafile = new TFile("Run1_4DataHist.root");
   TFile* datafile = new TFile("nd280fit2013output.root");
   TH2D* cc0pidata = (TH2D*)datafile->Get("cc0pidata");
   TH2D* cc1pidata = (TH2D*)datafile->Get("cc1pidata");
   TH2D* ccmpidata = (TH2D*)datafile->Get("ccmpidata");   

   
   TFile* file = new TFile("posteriorPredictive_v3_nomcerr.root");
   TH2D* cc0pi = (TH2D*)file->Get("cc0pi_0");
   TH2D* cc1pi = (TH2D*)file->Get("cc1pi_0");
   TH2D* ccmpi = (TH2D*)file->Get("ccmpi_0");
   
   TH2D* cc0piprobhi = (TH2D*)cc0pi->Clone("cc0piprobhi");
   cc0piprobhi->Reset();
   TH2D* cc1piprobhi = (TH2D*)cc1pi->Clone("cc1piprobhi");
   cc1piprobhi->Reset();
   TH2D* ccmpiprobhi = (TH2D*)ccmpi->Clone("ccmpiprobhi");
   ccmpiprobhi->Reset();

   int nzbins=5000;
   double zbins[5001];
   double zlo=0, zhi=1000;
   for(int i=0; i<=nzbins; i++)
      zbins[i]=zlo+(zhi-zlo)/double(nzbins)*i;


   TH3D* cc0pipp = new TH3D("cc0pipp","cc0pipp",cc0pi->GetXaxis()->GetNbins(),cc0pi->GetXaxis()->GetXbins()->GetArray(),
			    cc0pi->GetYaxis()->GetNbins(),cc0pi->GetYaxis()->GetXbins()->GetArray(),nzbins,zbins);
   TH3D* cc1pipp = new TH3D("cc1pipp","cc1pipp",cc1pi->GetXaxis()->GetNbins(),cc1pi->GetXaxis()->GetXbins()->GetArray(),
			    cc1pi->GetYaxis()->GetNbins(),cc1pi->GetYaxis()->GetXbins()->GetArray(),nzbins,zbins);
   TH3D* ccmpipp = new TH3D("ccmpipp","ccmpipp",ccmpi->GetXaxis()->GetNbins(),ccmpi->GetXaxis()->GetXbins()->GetArray(),
			    ccmpi->GetYaxis()->GetNbins(),ccmpi->GetYaxis()->GetXbins()->GetArray(),nzbins,zbins);

   
   TList* list = file->GetListOfKeys();
   
   TH1D* cc0pisum = new TH1D("cc0pisum","cc0pisum",100,16500,17500);
   TH1D* cc1pisum = new TH1D("cc1pisum","cc1pisum",100,3800,4200);
   TH1D* ccmpisum = new TH1D("ccmpisum","ccmpisum",100,3800,4200);


   for(int i=0; i<list->GetEntries(); i++)
   {
      TString keyname = list->At(i)->GetName();
      if(keyname.Contains("cc0pi"))
      {
	 cc0pi = (TH2D*)file->Get(keyname);
	 cc0pisum->Fill(cc0pi->Integral());
	 for(int j=1; j<=cc0pi->GetNbinsX(); j++)
	 {
	    for(int k=1; k<=cc0pi->GetNbinsY(); k++)
	    {
	       cc0pipp->Fill(cc0pi->GetXaxis()->GetBinCenter(j),cc0pi->GetYaxis()->GetBinCenter(k),cc0pi->GetBinContent(j,k));
	    }
	 }
      }
      if(keyname.Contains("cc1pi"))
      {
	 cc1pi = (TH2D*)file->Get(keyname);
	 cc1pisum->Fill(cc1pi->Integral());
	 for(int j=1; j<=cc1pi->GetNbinsX(); j++)
	 {
	    for(int k=1; k<=cc1pi->GetNbinsY(); k++)
	    {
	       cc1pipp->Fill(cc1pi->GetXaxis()->GetBinCenter(j),cc1pi->GetYaxis()->GetBinCenter(k),cc1pi->GetBinContent(j,k));
	    }
	 }
      }
      if(keyname.Contains("ccmpi"))
      {
	 ccmpi = (TH2D*)file->Get(keyname);
	 ccmpisum->Fill(ccmpi->Integral());
	 for(int j=1; j<=ccmpi->GetNbinsX(); j++)
	 {
	    for(int k=1; k<=ccmpi->GetNbinsY(); k++)
	    {
	       ccmpipp->Fill(ccmpi->GetXaxis()->GetBinCenter(j),ccmpi->GetYaxis()->GetBinCenter(k),ccmpi->GetBinContent(j,k));
	    }
	 }
      }
   }

   TCanvas* c0 = new TCanvas("c0","c0",0,0,1800,900);
   c0->Divide(3,1);
   c0->cd(1);
   cc0pisum->Draw();
   cout << cc0pidata->Integral() << " " << cc0pisum->GetMinimum() << " " << cc0pisum->GetMaximum() << endl;
   double x[2]={cc0pidata->Integral(),cc0pidata->Integral()};
   double y[2]={cc0pisum->GetMinimum(),cc0pisum->GetMaximum()};
   TGraph* g0 = new TGraph(2,x,y);
   g0->SetLineColor(kRed);
   g0->SetLineWidth(3);
   g0->Draw("L SAME");
   c0->cd(2);
   cc1pisum->Draw();
   x[0]=cc1pidata->Integral(); y[0]=cc1pisum->GetMinimum();
   x[1]=cc1pidata->Integral(); y[1]=cc1pisum->GetMaximum();
   TGraph* g1 = new TGraph(2,x,y);
   g1->SetLineColor(kRed);
   g1->SetLineWidth(3);
   g1->Draw("L SAME");
   c0->cd(3);
   ccmpisum->Draw();
   x[0]=ccmpidata->Integral(); y[0]=ccmpisum->GetMinimum();
   x[1]=ccmpidata->Integral(); y[1]=ccmpisum->GetMaximum();
   TGraph* gm = new TGraph(2,x,y);
   gm->SetLineColor(kRed);
   gm->SetLineWidth(3);
   gm->Draw("L SAME");


   TFile* outfile = new TFile("posteriorPredictive_v3_nomcerr_output.root","RECREATE");
   cc0pisum->Write("cc0pisum");
   cc1pisum->Write("cc1pisum");
   ccmpisum->Write("ccmpisum");
//   outfile->Close();

   c0->SaveAs("posteriorpredictivebins_v3_nomcerr.pdf(");


   TCanvas* c1 = new TCanvas("c1","c1",0,0,800,800);
   c1->cd();

   TF1* poisson = new TF1("poisson","TMath::PoissonI(x,[0])*[1]",0,1000);
   poisson->SetLineColor(kBlue);
   poisson->SetNpx(3000);

   double negLogL = 0;

   vector<double> means;
   vector<double> rms;
   vector<double> index;
   vector<double> ierr;
   vector<double> dataval;

   int inde=0;
   
   for(int j=1; j<=cc0pi->GetNbinsX(); j++)
   {
      for(int k=1; k<=cc0pi->GetNbinsY(); k++)
      {
	 TString name="pz0pi";
	 name+=j;
	 name+=k;
	 TH1D* projz0 = cc0pipp->ProjectionZ(name,j,j,k,k);
	 projz0->SetAxisRange(projz0->GetMean()-16*projz0->GetRMS(),projz0->GetMean()+16*projz0->GetRMS());
	 means.push_back(projz0->GetMean()); rms.push_back(projz0->GetRMS()); index.push_back(inde); ierr.push_back(0); inde++;
	 double ndata=cc0pidata->GetBinContent(j,k);
	 dataval.push_back(ndata);
	 double xx[2]={ndata,ndata};
	 double yy[2]={projz0->GetMinimum(),projz0->GetMaximum()};
	 TGraph* g = new TGraph(2,xx,yy);
	 g->SetLineColor(kRed);
	 g->SetLineWidth(3);
	 poisson->SetParameters(ndata+1e-9,1);
	 if(ndata==0)
	    poisson->SetRange(0,5);
	 else if(ndata-5*sqrt(ndata)<0)
	    poisson->SetRange(0,ndata+5*sqrt(ndata));
	 else
	    poisson->SetRange(ndata-5*sqrt(ndata),ndata+5*sqrt(ndata));

	 char title[70];
	 sprintf(title,"CC0#pi p(%5.0f-%5.0f) cos#theta(%1.2f-%1.2f)",cc0pi->GetXaxis()->GetBinLowEdge(j),cc0pi->GetXaxis()->GetBinLowEdge(j+1),
		 cc0pi->GetYaxis()->GetBinLowEdge(k),cc0pi->GetYaxis()->GetBinLowEdge(k+1));
	 poisson->SetTitle(title);
	 projz0->SetTitle(title);
	 projz0->Scale(1/projz0->Integral("width"));
	 projz0->Draw("");
	 poisson->Draw("SAME");
	 g->Draw("L SAME");
	 c1->Update();
	 c1->SaveAs("posteriorpredictivebins_v3_nomcerr.pdf");

	 double mean = projz0->GetMean();
	 negLogL+=(mean - xx[0] + xx[0] * TMath::Log(xx[0]/mean));

	 double prophi = projz0->Integral(projz0->FindBin(ndata),projz0->GetNbinsX())/projz0->Integral();
	 if(prophi==0)
	    prophi=1E-7;
	 cc0piprobhi->SetBinContent(j,k,prophi);
      }
   }

   cout << negLogL << endl;
   for(int j=1; j<=cc1pi->GetNbinsX(); j++)
   {
      for(int k=1; k<=cc1pi->GetNbinsY(); k++)
      {
	 TString name="pz1pi";
	 name+=j;
	 name+=k;
	 TH1D* projz1 = cc1pipp->ProjectionZ(name,j,j,k,k);
	 projz1->SetAxisRange(projz1->GetMean()-16*projz1->GetRMS(),projz1->GetMean()+16*projz1->GetRMS());
	 means.push_back(projz1->GetMean()); rms.push_back(projz1->GetRMS()); index.push_back(inde); ierr.push_back(0); inde++;
	 double ndata=cc1pidata->GetBinContent(j,k);
	 dataval.push_back(ndata);
	 double xx[2]={ndata,ndata};
	 double yy[2]={projz1->GetMinimum(),projz1->GetMaximum()};
	 TGraph* g = new TGraph(2,xx,yy);
	 g->SetLineColor(kRed);
	 g->SetLineWidth(3);
	 poisson->SetParameters(ndata+1e-9,1);
	 if(ndata==0)
	    poisson->SetRange(0,5);
	 else  if(ndata-5*sqrt(ndata)<0)
	    poisson->SetRange(0,ndata+5*sqrt(ndata));
	 else
	    poisson->SetRange(ndata-5*sqrt(ndata),ndata+5*sqrt(ndata));

	 char title[70];
	 sprintf(title,"CC1#pi p(%5.0f-%5.0f) cos#theta(%1.2f-%1.2f)",cc1pi->GetXaxis()->GetBinLowEdge(j),cc1pi->GetXaxis()->GetBinLowEdge(j+1),
		 cc1pi->GetYaxis()->GetBinLowEdge(k),cc1pi->GetYaxis()->GetBinLowEdge(k+1));
	 poisson->SetTitle(title);
	 poisson->Draw();
	 projz1->DrawNormalized("SAME");
	 g->Draw("L SAME");
	 c1->Update();
	 c1->SaveAs("posteriorpredictivebins_v3_nomcerr.pdf");

//	 getchar();
	 double mean = projz1->GetMean();
	 if(xx[0]>0)
	    negLogL+=(mean - xx[0] + xx[0] * TMath::Log(xx[0]/mean));
	 else 
	    negLogL += mean;
	    
	 cout << xx[0] << " " << mean << ": " << negLogL << endl;
	 
	 double prophi = projz1->Integral(projz1->FindBin(ndata),projz1->GetNbinsX())/projz1->Integral();
	 if(prophi==0)
	    prophi=1E-7;
	 cc1piprobhi->SetBinContent(j,k,prophi);
      }
   }
   cout << negLogL << endl;
   for(int j=1; j<=ccmpi->GetNbinsX(); j++)
   {
      for(int k=1; k<=ccmpi->GetNbinsY(); k++)
      {
	 TString name="pzmpi";
	 name+=j;
	 name+=k;
	 TH1D* projzm = ccmpipp->ProjectionZ(name,j,j,k,k);
	 projzm->SetAxisRange(projzm->GetMean()-16*projzm->GetRMS(),projzm->GetMean()+16*projzm->GetRMS());
	 means.push_back(projzm->GetMean()); rms.push_back(projzm->GetRMS()); index.push_back(inde); ierr.push_back(0); inde++;
	 double ndata=ccmpidata->GetBinContent(j,k);
	 dataval.push_back(ndata);
	 double xx[2]={ndata,ndata};
	 double yy[2]={projzm->GetMinimum(),projzm->GetMaximum()};
	 TGraph* g = new TGraph(2,xx,yy);
	 g->SetLineColor(kRed);
	 g->SetLineWidth(3);
	 poisson->SetParameters(ndata+1e-9,1);
	 if(ndata==0)
	    poisson->SetRange(0,5);
	 else	 if(ndata-5*sqrt(ndata)<0)
	    poisson->SetRange(0,ndata+5*sqrt(ndata));
	 else
	    poisson->SetRange(ndata-5*sqrt(ndata),ndata+5*sqrt(ndata));

	 char title[70];
	 sprintf(title,"CC oth p(%5.0f-%5.0f) cos#theta(%1.2f-%1.2f)",ccmpi->GetXaxis()->GetBinLowEdge(j),ccmpi->GetXaxis()->GetBinLowEdge(j+1),
		 ccmpi->GetYaxis()->GetBinLowEdge(k),ccmpi->GetYaxis()->GetBinLowEdge(k+1));
	 poisson->SetTitle(title);
	 poisson->Draw();
	 projzm->DrawNormalized("SAME");
	 g->Draw("L SAME");
	 c1->Update();
	 c1->SaveAs("posteriorpredictivebins_v3_nomcerr.pdf");

//	 getchar();

	 double mean = projzm->GetMean();
	 if(xx[0]>0)
	    negLogL+=(mean - xx[0] + xx[0] * TMath::Log(xx[0]/mean));
	 else 
	    negLogL += mean;

	 double prophi = projzm->Integral(projzm->FindBin(ndata),projzm->GetNbinsX())/projzm->Integral();
	 if(prophi==0)
	    prophi=1E-7;
	 ccmpiprobhi->SetBinContent(j,k,prophi);
      }
   }
   cout << negLogL << endl;

   c0->cd(1);
   cc0piprobhi->Draw("colz");
   c0->cd(2);
   cc1piprobhi->Draw("colz");
   c0->cd(3);
   ccmpiprobhi->Draw("colz");

   c0->SaveAs("posteriorpredictivebins_v3_nomcerr.pdf");

   TGraphErrors* gMeans = new TGraphErrors(means.size(),&index[0],&means[0],&ierr[0],&rms[0]);
   TGraph* gData = new TGraph(dataval.size(),&index[0],&dataval[0]);
   gData->SetMarkerStyle(25);
   gData->SetMarkerColor(kBlue);
   c1->cd();
   gMeans->Draw("AP");
   gData->Draw("P SAME");
   c1->SaveAs("posteriorpredictivebins_v3_nomcerr.pdf)");
   
   outfile->cd();
   gMeans->Write("Means");
   gData->Write("Data");

//    cc0pipp->Draw("box");
   
//    cc0pipp->ProjectionZ("pz",5,5,5,5,"d");

}
