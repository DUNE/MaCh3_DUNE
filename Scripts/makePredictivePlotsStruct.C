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
using namespace std;

struct sample
{
   TH2D* data;
   TH2D* data_norm;
   TH2D* nominal; 
   TH2D* nom_norm;
   TH3D* pp;
   TH1D* sum;
   TH2D* mean;
   TH2D* mean_norm;
   TH2D* lnL;
   TCanvas* canvas;
   int pad;
   TString key;
};

int nzbins=5000;
double zbins[5001];
double zlo=0, zhi=1000;

TF1* poisson = new TF1("poisson","TMath::PoissonI(x,[0])*[1]",0,1000);

void getDataNominal(sample& s, TFile* file, TString dataname, TString nomname)
{
   s.data = (TH2D*)file->Get(dataname);
   s.data->GetXaxis()->SetTitle("p_{#mu} (MeV)");
   s.data->GetYaxis()->SetTitle("cos(#theta_{#mu})");
   s.data->GetZaxis()->SetTitle("Counts");
   s.nominal = (TH2D*)file->Get(nomname);
   s.nominal->GetXaxis()->SetTitle("p_{#mu} (MeV)");
   s.nominal->GetYaxis()->SetTitle("cos(#theta_{#mu})");
   s.nominal->GetZaxis()->SetTitle("Counts");
   TString ppname = s.key;
   ppname+="pp";
   s.pp = new TH3D(ppname,ppname,
		    s.data->GetXaxis()->GetNbins(),s.data->GetXaxis()->GetXbins()->GetArray(),
		    s.data->GetYaxis()->GetNbins(),s.data->GetYaxis()->GetXbins()->GetArray(),
		    nzbins,zbins);
   TString meanname = s.key;
   meanname+="_mean";
   s.mean = (TH2D*)s.nominal->Clone(meanname);
   TString sumname = s.key;
   sumname+="_sum";
   s.sum = new TH1D(sumname,sumname,200,s.data->Integral()*0.8,s.data->Integral()*1.2);
   s.sum->GetXaxis()->SetTitle("N Events");
   s.sum->GetYaxis()->SetTitle("Counts");
   TString lnLname = s.key;
   lnLname+="_lnL";
   s.lnL = (TH2D*)s.nominal->Clone(lnLname);
}

void addThrow(sample& s, TFile* file, TString keyname)
{
   TH2D* hist  = (TH2D*)file->Get(keyname);
   s.sum->Fill(hist->Integral());
   for(int j=1; j<=hist->GetNbinsX(); j++)
   {
      for(int k=1; k<=hist->GetNbinsY(); k++)
      {
	 s.pp->Fill(hist->GetXaxis()->GetBinCenter(j),
		     hist->GetYaxis()->GetBinCenter(k),
		     hist->GetBinContent(j,k));
      }
   }
   
   hist->Delete();

}

void drawSum(sample& s)
{
   s.canvas->cd(s.pad);
   s.sum->Draw();
   cout << s.sum->GetMean() << " " << s.sum->GetMinimum() << " " << s.sum->GetMaximum() << endl;
   double x[2]={s.data->Integral(),s.data->Integral()};
   double y[2]={s.sum->GetMinimum(),s.sum->GetMaximum()};
   TGraph* g0 = new TGraph(2,x,y);
   g0->SetLineColor(kRed);
   g0->SetLineWidth(3);
   g0->Draw("L SAME");
   s.canvas->Update();
//   getchar();
}

void makeMean(sample& s, double& logL, TString pdfname)
{
   double negLogL=0;
   for(int j=1; j<=s.data->GetNbinsX(); j++)
   {
      for(int k=1; k<=s.data->GetNbinsY(); k++)
      {
	 TString name=s.key;
	 name+="_px_";
	 name+=j;
	 name+=k;
	 TH1D* projz = s.pp->ProjectionZ(name,j,j,k,k);
	 projz->SetAxisRange(projz->GetMean()-16*projz->GetRMS(),
			     projz->GetMean()+16*projz->GetRMS());
	 // means.push_back(projz->GetMean()); 
	 // rms.push_back(projz->GetRMS()); 
	 // index.push_back(inde); 
	 // ierr.push_back(0); inde++;
	 
	 double ndata=s.data->GetBinContent(j,k);
	 // dataval.push_back(ndata);

	 double xx[2]={ndata,ndata};
	 double yy[2]={projz->GetMinimum(),projz->GetMaximum()};

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
	 sprintf(title,"%s p(%5.0f-%5.0f) cos#theta(%1.2f-%1.2f)",s.key.Data(),
		 s.data->GetXaxis()->GetBinLowEdge(j),s.data->GetXaxis()->GetBinLowEdge(j+1),
		 s.data->GetYaxis()->GetBinLowEdge(k),s.data->GetYaxis()->GetBinLowEdge(k+1));
	 poisson->SetTitle(title);
	 projz->SetTitle(title);
	 projz->Scale(1/projz->Integral("width"));

	 // s.canvas->cd(s.pad);
	 // projz->Draw("");
	 // poisson->Draw("SAME");
	 // g->Draw("L SAME");
	 // s.canvas->Update();
	 // s.canvas->SaveAs(pdfname);

	 double mean = projz->GetMean();
	 double ll=0;
	 if(xx[0]>0)
	    ll=(mean - ndata + ndata * TMath::Log(ndata/mean));
	 else 
	    ll=mean;

	 negLogL+=ll;
	 s.mean->SetBinContent(j,k,mean);
	 if(mean-ndata>=0)
	    ll*=2;
	 else
	    ll*=-2;
	 s.lnL->SetBinContent(j,k,ll);
      }
   }

   logL+=negLogL;
   cout << "This sample -lnL: " << negLogL << "\tTotal -lnL: " << logL << endl;

}

void drawMean(sample& s)
{
   s.canvas->cd(s.pad);
   s.mean->Draw("colz");
}

void drawData(sample& s)
{
   s.canvas->cd(s.pad);
   s.data->Draw("colz");
}

void drawNominal(sample& s)
{
   s.canvas->cd(s.pad);
   s.nominal->Draw("colz");
}

void drawlnL(sample& s)
{
   s.canvas->cd(s.pad);
   s.lnL->GetZaxis()->SetTitle("-2lnL#times sign(MC-data)");
   s.lnL->Draw("colz");
}

void saveHistograms(sample& s, TFile* file)
{
   file->cd();
   s.data->Write();
   s.nominal->Write();
   s.mean->Write();
   s.lnL->Write();
   s.sum->Write();
}

void makeNormalized(sample& s)
{
   TString name = s.key;
   name+="data_norm";
   s.data_norm=(TH2D*)s.data->Clone(name);
   name = s.key;
   name+="nom_norm";
   s.nom_norm=(TH2D*)s.nominal->Clone(name);
   name = s.key;
   name+="mean_norm";
   s.mean_norm=(TH2D*)s.mean->Clone(name);


   for(int i=1; i<=s.data_norm->GetNbinsX(); i++)
   {
      for(int j=1; j<=s.data_norm->GetNbinsY(); j++)
      {
	 double dataval=s.data_norm->GetBinContent(i,j);
	 double binwidthx=s.data_norm->GetXaxis()->GetBinWidth(i);
 	 double binwidthy=s.data_norm->GetYaxis()->GetBinWidth(j);
	 s.data_norm->SetBinContent(i,j,dataval*100.0/binwidthx*0.01/binwidthy);
	 dataval = s.nom_norm->GetBinContent(i,j);
	 s.nom_norm->SetBinContent(i,j,dataval*100.0/binwidthx*0.01/binwidthy);
	 dataval = s.mean_norm->GetBinContent(i,j);
	 s.mean_norm->SetBinContent(i,j,dataval*100.0/binwidthx*0.01/binwidthy);

	 double dataerr=s.data_norm->GetBinError(i,j);
	 s.data_norm->SetBinError(i,j,dataerr*100.0/binwidthx*0.01/binwidthy);
	    
      }
   }
   s.data_norm->GetZaxis()->SetTitle("Events/100 MeV/0.01");
   s.nom_norm->GetZaxis()->SetTitle("Events/100 MeV/0.01");
   s.mean_norm->GetZaxis()->SetTitle("Events/100 MeV/0.01");

}

void drawNominalNorm(sample& s)
{
   s.canvas->cd(s.pad);
   s.nom_norm->GetXaxis()->SetRange(1,s.data->GetNbinsX()-1);
   s.nom_norm->GetYaxis()->SetRange(2,s.data->GetNbinsY());
   s.nom_norm->Draw("colz");
}
void drawDataNorm(sample& s)
{
   s.canvas->cd(s.pad);
   s.data_norm->GetXaxis()->SetRange(1,s.data->GetNbinsX()-1);
   s.data_norm->GetYaxis()->SetRange(2,s.data->GetNbinsY());
   s.data_norm->Draw("colz");
}
void drawMeanNorm(sample& s)
{
   s.canvas->cd(s.pad);
   s.mean_norm->GetXaxis()->SetRange(1,s.data->GetNbinsX()-1);
   s.mean_norm->GetYaxis()->SetRange(2,s.data->GetNbinsY());
   s.mean_norm->Draw("colz");
}


void projectXAndDrawNormalized(sample s, int lobin, int hibin,  TString pdfname="", bool drawleg=false)
{
   TString dataname = s.data->GetName();
   dataname+=hibin;
   TString premcname = s.nominal->GetName();
   premcname+=hibin;
   TString postmcname = s.mean->GetName();
   postmcname+=hibin;
   TH1D* pxdata = s.data->ProjectionX(dataname,lobin,hibin);
   TH1D* pxpremc = s.nominal->ProjectionX(premcname,lobin,hibin);
   TH1D* pxpostmc = s.mean->ProjectionX(postmcname,lobin,hibin);

   char title[1000];
   double lo = s.data->GetYaxis()->GetBinLowEdge(lobin);
   double hi = s.data->GetYaxis()->GetBinUpEdge(hibin);

   if(hibin!=-1)
   {
      sprintf(title,"%s, %2.2f < cos#theta < %2.2f",pxdata->GetTitle(),lo,hi);
      pxdata->SetTitle(title);
   }
   
   double lnLpre=0;
   double lnLpost=0;

   for(int i=1; i<=pxdata->GetNbinsX(); i++)
   {
      double dataval=pxdata->GetBinContent(i);
      double dataerr=pxdata->GetBinError(i);
      double premcval=pxpremc->GetBinContent(i);
      double postmcval=pxpostmc->GetBinContent(i);
      double binwidthx=pxdata->GetXaxis()->GetBinWidth(i);

      if(dataval>0)
	 lnLpre+=premcval-dataval+dataval*TMath::Log(dataval/premcval);
      else
	 lnLpre+=premcval;
      if(dataval>0)
	 lnLpost+=postmcval-dataval+dataval*TMath::Log(dataval/postmcval);
      else
	 lnLpost+=postmcval;
      
      pxdata->SetBinContent(i,dataval*100.0/binwidthx);
      pxdata->SetBinError(i,dataerr*100.0/binwidthx);
      pxpremc->SetBinContent(i,premcval*100.0/binwidthx);
      pxpostmc->SetBinContent(i,postmcval*100.0/binwidthx);
   } 
   std::cout << title << " prefit -2lnL: " << 2*lnLpre << " postfit -2lnL: " << 2*lnLpost << std::endl;
//   std::cout << "v3-v2: " << 2*(lnLpost-lnLbanff) << std::endl;
   
   double maxval=0.0;
   if(pxdata->GetMaximum()>maxval)
      maxval=pxdata->GetMaximum();
   if(pxpremc->GetMaximum()>maxval)
	 maxval=pxpremc->GetMaximum(); 
   if(pxpostmc->GetMaximum()>maxval)
	 maxval=pxpostmc->GetMaximum(); 
   pxdata->SetMinimum(0.0);
   pxdata->SetMaximum(1.1*maxval);
   
   
   s.canvas->cd(s.pad);
//   s.canvas->GetPad(s.pad)->SetLogy();
   pxdata->GetXaxis()->SetRange(1,pxdata->GetXaxis()->GetNbins()-1);
   pxdata->GetYaxis()->SetTitle("Events/100 MeV");
   pxdata->Draw("E");
   pxpremc->SetLineColor(kGray+2);
   pxpremc->Draw("SAME");
   pxpostmc->SetLineColor(kBlue);
   pxpostmc->Draw("SAME");
   
   if(drawleg)
   {
      TLegend* leg = new TLegend(0.5,0.5,0.8,0.8);
      leg->SetFillColor(0);
      leg->AddEntry(pxdata,"Data","PL");
      leg->AddEntry(pxpremc,"Prefit MC","L");
      leg->AddEntry(pxpostmc,"Postfit MC","L");
	 leg->Draw("SAME");
   }
      
   if(pdfname!="")
      s.canvas->SaveAs(pdfname);
   
}

void projectYAndDrawNormalized(sample s, int lobin, int hibin,  TString pdfname="", bool drawleg=false)
{
   TString dataname = s.data->GetName();
   dataname+=hibin;
   TString premcname = s.nominal->GetName();
   premcname+=hibin;
   TString postmcname = s.mean->GetName();
   postmcname+=hibin;
   TH1D* pydata = s.data->ProjectionY(dataname,lobin,hibin);
   TH1D* pypremc = s.nominal->ProjectionY(premcname,lobin,hibin);
   TH1D* pypostmc=s.mean->ProjectionY(postmcname,lobin,hibin);
      
   char title[1000];
   double lo = pydata->GetXaxis()->GetBinLowEdge(lobin);
   double hi = pydata->GetXaxis()->GetBinUpEdge(hibin);
//   std::cout << lo << " " << hi << std::endl;
   if(hibin!=-1)
   {
      sprintf(title,"%s, %2.2f < p < %2.2f",pydata->GetTitle(),lo,hi);
      pydata->SetTitle(title);
   }
   double lnLpre=0;
   double lnLpost=0;
   for(int i=1; i<=pydata->GetNbinsX(); i++)
   {
      double dataval=pydata->GetBinContent(i);
      double dataerr=pydata->GetBinError(i);
      double premcval=pypremc->GetBinContent(i);
      double binwidthx=pydata->GetXaxis()->GetBinWidth(i);

      if(dataval>0)
	 lnLpre+=premcval-dataval+dataval*TMath::Log(dataval/premcval);
      else
	 lnLpre+=premcval;

      pydata->SetBinContent(i,dataval*0.01/binwidthx);
      pydata->SetBinError(i,dataerr*0.01/binwidthx);
      pypremc->SetBinContent(i,premcval*0.01/binwidthx);
      
      double postmcval=pypostmc->GetBinContent(i);
      pypostmc->SetBinContent(i,postmcval*0.01/binwidthx);
      
      if(dataval>0)
	 lnLpost+=postmcval-dataval+dataval*TMath::Log(dataval/postmcval);
      else
	 lnLpost+=postmcval;
      
   }

   std::cout << title << " prefit -2lnL: " << 2*lnLpre << " postfit -2lnL: " << 2*lnLpost << std::endl;

   double maxval=0.0;
   if(pydata->GetMaximum()>maxval)
      maxval=pydata->GetMaximum();
   if(pypremc->GetMaximum()>maxval)
      maxval=pypremc->GetMaximum(); 
   if(pypostmc->GetMaximum()>maxval)
      maxval=pypostmc->GetMaximum(); 
   pydata->SetMinimum(0.0);
   pydata->SetMaximum(1.1*maxval);
   
//   pydata->GetYaxis()->SetTitle("Events");
   pydata->GetYaxis()->SetTitle("Events/0.01");
   pydata->GetXaxis()->SetRange(2,pydata->GetXaxis()->GetNbins());
   pydata->Draw("E");
   pypremc->SetLineColor(kGray+2);
   pypremc->Draw("SAME");
   pypostmc->SetLineColor(kBlue);
   pypostmc->Draw("SAME");
   
   if(drawleg)
   {
      TLegend* leg = new TLegend(0.2,0.5,0.5,0.8);
      leg->SetFillColor(0);
      leg->AddEntry(pydata,"Data","PL");
      leg->AddEntry(pypremc,"Prefit MC","L");
      leg->AddEntry(pypostmc,"Postfit MC","L");
//       leg->AddEntry(pybanff,"No CCQE Priors Postfit MC","L");
      leg->Draw("SAME");

   }

   if(pdfname!="")
      s.canvas->SaveAs(pdfname);

}

void drawLogLDist(sample& s)
{
   TString name = s.key;
   name+="lnl_1D";
   TH1D* tmp = new TH1D(name,name,75,-15,15);
   for(int i=1; i<=s.lnL->GetNbinsX(); i++)
      for(int j=1; j<=s.lnL->GetNbinsY(); j++)
	 tmp->Fill(s.lnL->GetBinContent(i,j));

   s.canvas->cd(s.pad);
   tmp->GetXaxis()->SetTitle("-2lnL#times sign(MC-data)");
   tmp->Draw();
   std::cout << s.key << ": " << tmp->GetMean() << " " << tmp->GetRMS() <<  " " << tmp->GetRMS()/sqrt(tmp->Integral()) << std::endl;

}


void makePredictivePlots(TString infile, TString outprefix)
{
   poisson->SetLineColor(kBlue);
   poisson->SetNpx(3000);
   for(int i=0; i<=nzbins; i++)
      zbins[i]=zlo+(zhi-zlo)/double(nzbins)*i;

   const int nsamples=14;
   TString samplekeys[nsamples] = {"cc0pi","cc1pi","ccmpi",
   				   "anu_ccqe","anu_ccnqe","wsbkg_ccqe","wsbkg_ccnqe",
   				   "cc0pi_2","cc1pi_2","ccmpi_2",
   				   "anu_ccqe_2","anu_ccnqe_2","wsbkg_ccqe_2","wsbkg_ccnqe_2"};
   TString datakeys[nsamples] = {"data_cc0pi_1","data_cc1pi_1","data_ccmpi_1",
   				 "data_antinu_ccqe_1","data_antinu_ccnqe_1",
   				 "data_wsbkg_ccqe_1","data_wsbkg_ccnqe_1",
   				 "data_cc0pi_2","data_cc1pi_2","data_ccmpi_2",
   				 "data_antinu_ccqe_2","data_antinu_ccnqe_2",
   				 "data_wsbkg_ccqe_2","data_wsbkg_ccnqe_2"};
   TString nomkeys[nsamples] = {"nu_cc0pi_1","nu_cc1pi_1","nu_ccmpi_1",
   				"antinu_ccqe_1","antinu_ccnqe_1",
   				"wsbkg_ccqe_1","wsbkg_ccnqe_1",
   				"nu_cc0pi_2","nu_cc1pi_2","nu_ccmpi_2",
   				"antinu_ccqe_2","antinu_ccnqe_2",
   				"wsbkg_ccqe_2","wsbkg_ccnqe_2"};
   sample samples[nsamples];
   int pad[nsamples]={1,2,3,5,6,7,8,1,2,3,5,6,7,8};


   // const int nsamples=7;
   // TString samplekeys[nsamples] = {"cc0pi","cc1pi","ccmpi",
   // 				   "anu_ccqe","anu_ccnqe","wsbkg_ccqe","wsbkg_ccnqe"};
   // TString datakeys[nsamples] = {"data_cc0pi_1","data_cc1pi_1","data_ccmpi_1",
   // 				 "data_antinu_ccqe_1","data_antinu_ccnqe_1",
   // 				 "data_wsbkg_ccqe_1","data_wsbkg_ccnqe_1"};
   // TString nomkeys[nsamples] = {"nu_cc0pi_1","nu_cc1pi_1","nu_ccmpi_1",
   // 				"antinu_ccqe_1","antinu_ccnqe_1",
   // 				"wsbkg_ccqe_1","wsbkg_ccnqe_1"};
   // sample samples[nsamples];
   // int pad[nsamples]={1,2,3,5,6,7,8};


   // const int nsamples=7;
   // TString samplekeys[nsamples] = {"cc0pi_2","cc1pi_2","ccmpi_2",
   // 				   "anu_ccqe_2","anu_ccnqe_2","wsbkg_ccqe_2","wsbkg_ccnqe_2"};
   // TString datakeys[nsamples] = {"data_cc0pi_2","data_cc1pi_2","data_ccmpi_2",
   // 				 "data_antinu_ccqe_2","data_antinu_ccnqe_2",
   // 				 "data_wsbkg_ccqe_2","data_wsbkg_ccnqe_2"};
   // TString nomkeys[nsamples] = {"nu_cc0pi_2","nu_cc1pi_2","nu_ccmpi_2",
   // 				"antinu_ccqe_2","antinu_ccnqe_2",
   // 				"wsbkg_ccqe_2","wsbkg_ccnqe_2"};
   // sample samples[nsamples];
   // int pad[nsamples]={1,2,3,5,6,7,8};


   TFile* datafile = new TFile("antinunominal.root");

   for(int i=0; i<nsamples; i++)
   {
      samples[i].key=samplekeys[i];
      samples[i].pad=pad[i];
      getDataNominal(samples[i],datafile,datakeys[i],nomkeys[i]);
   }
   
   TFile* file = new TFile(infile);
   TList* list = file->GetListOfKeys();

   TCanvas* c1 = new TCanvas("c1","FGD1",0,0,1800,900);
   c1->Divide(4,2);
   TCanvas* c2 = new TCanvas("c2","FGD2",0,0,1800,900);
   c2->Divide(4,2);

   for(int i=0; i<list->GetEntries(); i++)
   {
      TString keyname = list->At(i)->GetName();
      for(int j=0; j<nsamples; j++)
      {
   	 if(j<7 && keyname.Contains(samples[j].key) && !keyname.Contains("_2_"))
	 {
   	    addThrow(samples[j],file,keyname);
	    samples[j].canvas=c1;
	 }
   	 if(j>=7 && keyname.Contains(samples[j].key) && keyname.Contains("_2_"))
//   	 if(j>=0 && keyname.Contains(samples[j].key) && keyname.Contains("_2_"))
	 {
   	    addThrow(samples[j],file,keyname);
	    samples[j].canvas=c2;
	 }

      }

   }
   
   

   for(int i=0; i<nsamples; i++)
      drawSum(samples[i]);

   c1->Update();
   c2->Update();
   
   TString pdfname = outprefix;
   pdfname+=".pdf(";
   c1->SaveAs(pdfname);
   pdfname.ReplaceAll("(","");
   c2->SaveAs(pdfname);
   
   double negLogL=0;
   for(int i=0; i<nsamples; i++)
   {
      c1->Clear("D");
      c2->Clear("D");
      makeMean(samples[i],negLogL,pdfname);
   }

   for(int i=0; i<nsamples; i++)
      drawMean(samples[i]);
  
   c1->SaveAs(pdfname);
   c2->SaveAs(pdfname);

   for(int i=0; i<nsamples; i++)
      drawlnL(samples[i]);
   
   c1->SaveAs(pdfname);
   c2->SaveAs(pdfname);

   for(int i=0; i<nsamples; i++)
   {
      makeNormalized(samples[i]);
      drawDataNorm(samples[i]);
   }
   c1->SaveAs(pdfname);
   c2->SaveAs(pdfname);

   for(int i=0; i<nsamples; i++)
      drawNominalNorm(samples[i]);
   
   c1->SaveAs(pdfname);
   c2->SaveAs(pdfname);

   for(int i=0; i<nsamples; i++)
      drawMeanNorm(samples[i]);
   
   c1->SaveAs(pdfname);
   c2->SaveAs(pdfname);


   for(int i=0; i<nsamples; i++)
   {
     c1->Clear("D");
     c2->Clear("D");
     for(int j=1; j<=samples[i].data->GetNbinsY(); j++)
	projectXAndDrawNormalized(samples[i],j,j,pdfname);
   }

   for(int i=0; i<nsamples; i++)
   {
      bool leg=false;
      if(i==0 || i==7)
	 leg=true;
      projectXAndDrawNormalized(samples[i],0,-1,"",leg);
   }

   c1->SaveAs(pdfname);
   c2->SaveAs(pdfname);
   
   for(int i=0; i<nsamples; i++)
   {
      bool leg=false;
      if(i==0 || i==7)
	 leg=true;
      projectYAndDrawNormalized(samples[i],0,-1,"",leg);
   }

   c1->SaveAs(pdfname);
   c2->SaveAs(pdfname);

   for(int i=0; i<nsamples; i++)
      drawLogLDist(samples[i]);

   c1->SaveAs(pdfname);
   c2->SaveAs(pdfname);


   TString outfilename = outprefix;
   outfilename+=".root";
   TFile* outfile = new TFile(outfilename,"RECREATE");

   for(int i=0; i<nsamples; i++)
      saveHistograms(samples[i],outfile);
   
   pdfname+="]";
   c1->SaveAs(pdfname);

}
