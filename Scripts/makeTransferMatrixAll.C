#include "TObjArray.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TString.h"
#include "TStyle.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TF2.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include <iostream>

using namespace std;

void makeTransferMatrixAll(TString addstring, TString outprefix)
{
   TChain* c = new TChain("posteriors","");
   c->Add(addstring);


   char stepcut[300];
//   sprintf(stepcut,"(LocalEntry$>=20000&&LogL<5000)*exp(-0.5*pow((xsec_4-1)/0.3333,2))*exp(-0.5*pow((xsec_5-0.27)/0.35228,2))");
   sprintf(stepcut,"(LocalEntry$>=20000&&LogL<5000)");
  
   //c->Draw("Accb_0",stepcut);
   //c->Draw("Accb_1",stepcut);
   std::cout << "Get List of Branches " << std::endl;
   TObjArray* brlis = (TObjArray*)c->GetListOfBranches();
   std::cout << "Get Number Branches " << std::endl;
   int nbr = brlis->GetEntries();
   std::cout << " nbr " << nbr << std::endl;
   TString bnames[850];
   int ndraw=0;
   for(int i=0; i<nbr; i++)
   {
      TBranch* br = (TBranch*)brlis->At(i);
      TString bname = br->GetName();
      if(bname.Contains("b_"))
      {
	 TString index = bname;
	 index.ReplaceAll("b_","");
	 if(!bname.Contains("LogL"))
	 {
	    bnames[ndraw]=bname;
	    cout << bnames[ndraw] << endl;
	    ndraw++;
	 }
      }
      if(bname.Contains("fsi_"))
      {
	 TString index = bname;
	 index.ReplaceAll("fsi_","");
	 if(!bname.Contains("LogL"))
	 {
	    bnames[ndraw]=bname;
	    cout << bnames[ndraw] << endl;
	    ndraw++;
	 }
      }

      if(bname.Contains("xsec_"))
      {
	 TString index = bname;
	 index.ReplaceAll("xsec_","");
	 if(!bname.Contains("LogL"))
	 {
	    bnames[ndraw]=bname;
	    cout << bname << endl;
	    ndraw++;
	 }
      }
      if(bname.Contains("ndd_"))
      {
	 TString index = bname;
	 index.ReplaceAll("ndd_","");
	 if(!bname.Contains("LogL"))
	 {
	    bnames[ndraw]=bname;
	    cout << bname << endl;
	    ndraw++;
	 }
      }

//       if(bname.Contains("b_") && !bname.Contains("Prob") && !bname.Contains("LogL"))
//       {
// 	 bnames[ndraw]=bname;
// 	 ndraw++;
//      } 
   }

   gStyle->SetOptFit(111);
   TCanvas* c0 = new TCanvas("c0","c0",0,0,800,800);
   std::cout << "Print pdf file  " << std::endl;

   TString canvasname=outprefix;
   canvasname+=".pdf[";
   c0->Print(canvasname);				
   canvasname.ReplaceAll("[","");
   TH1F *hpost = new TH1F("hpost","hpost",100,0,2);    
   TH2F *hpost2 = new TH2F("hpost2","hpost2",200,-2,2,200,-2,2);    
   TF1 *gauss = new TF1("gauss","[0]/sqrt(2.0*3.14159)/[2]*TMath::Exp(-0.5*pow(x-[1],2)/[2]/[2])",-5,5);
   TF1* gaus = new TF1("gaus","gaus",-5,5);
   TVectorD* mean_vec = new TVectorD(ndraw);
   TVectorD* err_vec = new TVectorD(ndraw); 
   TVectorD* gaus_mean_vec = new TVectorD(ndraw);
   TVectorD* gaus_err_vec = new TVectorD(ndraw); 
   TMatrixT<double>* covariance = new TMatrixT<double>(ndraw,ndraw);

   for(int i=0; i<ndraw; i++)
   {
      for(int j=0; j<=i; j++)
      {
//	 cout << bnames[i] << "," << bnames[j] << endl;
	 if(i==j)
	 {
	    cout << bnames[i] << std::endl;
	    hpost->SetBins(100,0.0,2.0);
            // If beam systematic
	    if (bnames[i].Contains("b_")) {
	       hpost->SetBins(50,0.5,1.5);
            // If MACCQE (== xsec_0)
            } else if (bnames[i]=="xsec_0") {
	       hpost->SetBins(100,0.7,1.3);
            } else if (bnames[i]=="xsec_13" ||  bnames[i]=="xsec_24" ||  bnames[i]=="xsec_25" ||  bnames[i]=="xsec_20" ||  bnames[i]=="xsec_21" ||  bnames[i]=="xsec_22" ||  bnames[i]=="xsec_23" || bnames[i].Contains("fsi")) {
               hpost->SetBins(100,-1,1);
            }
	    hpost->SetTitle(bnames[i]);
	    if( bnames[i]=="xsec_11" ||  bnames[i]=="xsec_22" ||  bnames[i]=="xsec_23" ||  bnames[i]=="xsec_18" ||  bnames[i]=="xsec_19" 
		||  bnames[i]=="xsec_20" ||  bnames[i]=="xsec_21" || bnames[i].Contains("fsi"))
	       hpost->SetBins(100,-1,1);

	    c->Project("hpost", bnames[i], stepcut);

	    double mean = hpost->GetMean(); 
	    double rms = hpost->GetRMS();
	    double peakval = hpost->GetBinCenter(hpost->GetMaximumBin() ); 
	    gauss->SetRange(mean-1.5*rms,mean+1.5*rms);
	    gauss->SetParameters(hpost->GetMaximum()*rms*sqrt(2*3.14),peakval,rms);
	    hpost->Fit("gauss","Rq");

	    (*mean_vec)(i)=mean;
	    (*err_vec)(i)=rms;
	    (*gaus_mean_vec)(i)=gauss->GetParameter(1);
	    (*gaus_err_vec)(i)=gauss->GetParameter(2);
	    (*covariance)(i,i)=rms*rms;

	    hpost->Draw();
	    c0->Print(canvasname);

	 } else {/*
	    TString drawcmd = bnames[j];
	    drawcmd+=":";
	    drawcmd+=bnames[i];
	    c->Project("hpost2",drawcmd,stepcut);
	    hpost2->Draw("colz");
	    cout << "covariance: " <<   hpost2->GetCovariance() << endl;
	    (*covariance)(i,j)=hpost2->GetCovariance();
	    (*covariance)(j,i)=hpost2->GetCovariance();
	    
	    int xb,yb,zb;
	    hpost2->GetBinXYZ(hpost2->GetMaximumBin(),xb,yb,zb);
	    cout << hpost2->GetXaxis()->GetBinCenter(xb) << " " << hpost2->GetYaxis()->GetBinCenter(yb) << endl;
	    c0->Update();
	  */
	 }
      }

   }


   canvasname+=")";
   c0->Print(canvasname);

   vector<double> val;
   vector<double> vale;
   vector<double> ind;
   vector<double> inde;

   for(int i=0; i<mean_vec->GetNrows(); i++) {
      val.push_back((*mean_vec)(i));
      vale.push_back((*err_vec)(i));
      ind.push_back(i);
      inde.push_back(0.0);
   }
   
   TGraphErrors* gParams = new TGraphErrors(val.size(),&ind[0],&val[0],&inde[0],&vale[0]);

   TH2D* hCov = new TH2D("hCov","hCov",covariance->GetNrows(),0,covariance->GetNrows(),
			 covariance->GetNcols(),0,covariance->GetNcols());
   TH2D* hCovSq = new TH2D("hCovSq","hCovSq",covariance->GetNrows(),0,covariance->GetNrows(),
			   covariance->GetNcols(),0,covariance->GetNcols());

   for(int i=0; i<covariance->GetNrows(); i++) {
      for(int j=0; j<covariance->GetNcols(); j++) {
	 double va = (*covariance)(i,j);
	 hCov->SetBinContent(i+1,j+1,va);
	 hCovSq->SetBinContent(i+1,j+1,((va>0)-(va<0))*sqrt(fabs(va)));
      }
   }


   TString rootfilename=outprefix;
   rootfilename+=".root";
   TFile* file = new TFile(rootfilename,"RECREATE");
   mean_vec->Write("FitParameters");
   err_vec->Write("FitErrors");
   gaus_mean_vec->Write("FitParameters_gauss");
   gaus_err_vec->Write("FitErrors_gauss");
   covariance->Write("postfit_cov");
   gParams->Write("mach3params");
   hCov->Write();
   hCovSq->Write();
   file->Close();

}
