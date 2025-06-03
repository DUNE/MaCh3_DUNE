// Quick script to add 'norm' hists to SK minituples (as required in samplePDFNue.cpp and samplePDFNumu.cpp. Should be run in root from the MaCh3 directory:
// " root
//   .L utils/MakeNormHists.cxx++
//   MakeNormHists() "

// Will check for 'norm' hist in the file. If it exists, will update it with the numbers given below. If it doesn't exists, will create it.


#include <stdlib.h>
#include "iostream"

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"


void MakeNormHists()
{
  // Numbers for normalisation
  double numerators[16];
  numerators[0]  = 1e21; // NHC numu
  numerators[1]  = 1e21; // NHC nue
  numerators[2]  = 1e21; // NHC numubar
  numerators[3]  = 1e21; // NHC nuebar
  numerators[4]  = 1e21; // NHC signue
  numerators[5]  = 1e21; // NHC signuebar
  numerators[6]  = 1e21; // NHC signnumu
  numerators[7]  = 1e21; // NHC signnumubar
  numerators[8]  = 1e21; //RHC numu
  numerators[9]  = 1e21; //RHC nue
  numerators[10]  = 1e21; //RHC numubar
  numerators[11]  = 1e21; //RHC nuebar
  numerators[12]  = 1e21; //RHC signue
  numerators[13]  = 1e21; //RHC signuebar
  numerators[14]  = 1e21; //RHC signumu
  numerators[15]  = 1e21; //RHC signumubar
  // numerators[6]  = 205.213957312; // RHC numu
  // numerators[7]  = 8.79315608601; // RHC nue
  // numerators[8]  = 338.496583643; // RHC numubar
  // numerators[9]  = 6.52180351576; // RHC nuebar
  // numerators[10] = 211.81697101 ; // RHC signue
  // numerators[11] = 357.710318015; // RHC signuebar

  double denominators[16];
  denominators[0] = 1.0; // NHC numu
  denominators[1] = 1.0; // NHC nue
  denominators[2] = 1.0; // NHC numubar
  denominators[3] = 1.0; // NHC nuebar
  denominators[4] = 1.0; // NHC signue
  denominators[5] = 1.0; // NHC signuebar
  denominators[6] = 1.0; // NHC signumu
  denominators[7] = 1.0; // NHC signumubar
  denominators[8] = 1.0; // RHC numu
  denominators[9] = 1.0; // RHC nue
  denominators[10] = 1.0; // RHC numubar
  denominators[11] = 1.0; // RHC nuebar
  denominators[12] = 1.0; // RHC signue
  denominators[13] = 1.0; // RHC signuebar
  denominators[14] = 1.0; // RHC signumu
  denominators[15] = 1.0; // RHC signumubar
  // denominators[6] = 634250; // RHC numu
  // denominators[7] = 126829; // RHC nue
  // denominators[8] = 634817; // RHC numubar
  // denominators[9] = 126982; // RHC nuebar
  // denominators[10] = 126923; // RHC signue
  // denominators[11] = 126928; // RHC signuebar

  // POT used to create MC
  //double pot = 1E21;
  //double pot_anu = 1E21;

  //folder containing mtuples
  char* mtuple_folder = "/vols/t2k/users/ljw20/data/DUNE_2021/DUNE_2021_CAFs";
  char* unsplit_mtuple_folder = "/vols/t2k/users/ljw20/data/DUNE_2021/DUNE_2021_HaddedCAFs";
  // Names of SK files 
  // Note: it's important that the first index follows the same order as the indices in 'numerators' and 'denominators' (ie. numu / nue / numubar / nuebar / signue / signuebar)
  char *sk_filenames[6][8];
  char *unsplit_mtuples[6];

  //These files get their POT norms from FHC non swap
  sk_filenames[0][0]  = (char *)"FD_FHC_ger_numu_x_numu_numuselec.root"; // NHC numu
  sk_filenames[0][1]  = (char *)"FD_FHC_ger_numu_x_numu_nueselec.root";
  sk_filenames[0][2]  = (char *)"FD_FHC_ger_nue_x_nue_numuselec.root"; // NHC nue
  sk_filenames[0][3]  = (char *)"FD_FHC_ger_nue_x_nue_nueselec.root"; //
  sk_filenames[0][4]  = (char *)"FD_FHC_ger_numubar_x_numubar_numuselec.root"; // NHC numubar
  sk_filenames[0][5]  = (char *)"FD_FHC_ger_numubar_x_numubar_nueselec.root";
  sk_filenames[0][6]  = (char *)"FD_FHC_ger_nuebar_x_nuebar_numuselec.root";// NHC nuebar
  sk_filenames[0][7]  = (char *)"FD_FHC_ger_nuebar_x_nuebar_nueselec.root";
  //These files get their POT norms from FHC nue swap
  sk_filenames[1][0]  = (char *)"FD_FHC_ger_numu_x_nue_numuselec.root"; // NHC signnue
  sk_filenames[1][1]  = (char *)"FD_FHC_ger_numu_x_nue_nueselec.root";
  sk_filenames[1][2]  = (char *)"FD_FHC_ger_nue_x_nutau_numuselec.root"; // NHC signutau
  sk_filenames[1][3]  = (char *)"FD_FHC_ger_nue_x_nutau_nueselec.root";
  sk_filenames[1][4]  = (char *)"FD_FHC_ger_numubar_x_nuebar_numuselec.root"; // NHC signuebar
  sk_filenames[1][5]  = (char *)"FD_FHC_ger_numubar_x_nuebar_nueselec.root"; // NHC signuebar
  sk_filenames[1][6]  = (char *)"FD_FHC_ger_nuebar_x_nutaubar_numuselec.root"; // NHC signutaubar
  sk_filenames[1][7]  = (char *)"FD_FHC_ger_nuebar_x_nutaubar_nueselec.root"; // NHC signutaubar
  //These files get their POT norms from FHC tauswap
  sk_filenames[2][0]  = (char *)"FD_FHC_ger_nue_x_numu_numuselec.root"; // NHC signumu
  sk_filenames[2][1]  = (char *)"FD_FHC_ger_nue_x_numu_nueselec.root"; // NHC signumu
  sk_filenames[2][2]  = (char *)"FD_FHC_ger_numu_x_nutau_numuselec.root"; // NHC signutau
  sk_filenames[2][3]  = (char *)"FD_FHC_ger_numu_x_nutau_nueselec.root"; // NHC signutau
  sk_filenames[2][4]  = (char *)"FD_FHC_ger_nuebar_x_numubar_numuselec.root"; // NHC signumubar
  sk_filenames[2][5]  = (char *)"FD_FHC_ger_nuebar_x_numubar_nueselec.root"; // NHC signumubar
  sk_filenames[2][6]  = (char *)"FD_FHC_ger_numubar_x_nutaubar_numuselec.root"; // NHC signutaubar
  sk_filenames[2][7]  = (char *)"FD_FHC_ger_numubar_x_nutaubar_nueselec.root"; // NHC signutaubar
  //These files get their POT norms from RHC non swap
  sk_filenames[3][0]  = (char *)"FD_RHC_ger_numu_x_numu_numuselec.root"; // NHC numu
  sk_filenames[3][1]  = (char *)"FD_RHC_ger_numu_x_numu_nueselec.root";
  sk_filenames[3][2]  = (char *)"FD_RHC_ger_nue_x_nue_numuselec.root"; // NHC nue
  sk_filenames[3][3]  = (char *)"FD_RHC_ger_nue_x_nue_nueselec.root"; //
  sk_filenames[3][4]  = (char *)"FD_RHC_ger_numubar_x_numubar_numuselec.root"; // NHC numubar
  sk_filenames[3][5]  = (char *)"FD_RHC_ger_numubar_x_numubar_nueselec.root";
  sk_filenames[3][6]  = (char *)"FD_RHC_ger_nuebar_x_nuebar_numuselec.root";// NHC nuebar
  sk_filenames[3][7]  = (char *)"FD_RHC_ger_nuebar_x_nuebar_nueselec.root";
  //These files get their POT norms from RHC nue swap
  sk_filenames[4][0]  = (char *)"FD_RHC_ger_numu_x_nue_numuselec.root"; // NHC signnue
  sk_filenames[4][1]  = (char *)"FD_RHC_ger_numu_x_nue_nueselec.root";
  sk_filenames[4][2]  = (char *)"FD_RHC_ger_nue_x_nutau_numuselec.root"; // NHC signutau
  sk_filenames[4][3]  = (char *)"FD_RHC_ger_nue_x_nutau_nueselec.root";
  sk_filenames[4][4]  = (char *)"FD_RHC_ger_numubar_x_nuebar_numuselec.root"; // NHC signuebar
  sk_filenames[4][5]  = (char *)"FD_RHC_ger_numubar_x_nuebar_nueselec.root"; // NHC signuebar
  sk_filenames[4][6]  = (char *)"FD_RHC_ger_nuebar_x_nutaubar_numuselec.root"; // NHC signutaubar
  sk_filenames[4][7]  = (char *)"FD_RHC_ger_nuebar_x_nutaubar_nueselec.root"; // NHC signutaubar
  //These files get their POT norms from RHC tauswap
  sk_filenames[5][0]  = (char *)"FD_RHC_ger_nue_x_numu_numuselec.root"; // NHC signumu
  sk_filenames[5][1]  = (char *)"FD_RHC_ger_nue_x_numu_nueselec.root"; // NHC signumu
  sk_filenames[5][2]  = (char *)"FD_RHC_ger_numu_x_nutau_numuselec.root"; // NHC signutau
  sk_filenames[5][3]  = (char *)"FD_RHC_ger_numu_x_nutau_nueselec.root"; // NHC signutau
  sk_filenames[5][4]  = (char *)"FD_RHC_ger_nuebar_x_numubar_numuselec.root"; // NHC signumubar
  sk_filenames[5][5]  = (char *)"FD_RHC_ger_nuebar_x_numubar_nueselec.root"; // NHC signumubar
  sk_filenames[5][6]  = (char *)"FD_RHC_ger_numubar_x_nutaubar_numuselec.root"; // NHC signutaubar
  sk_filenames[5][7]  = (char *)"FD_RHC_ger_numubar_x_nutaubar_nueselec.root"; // NHC signutaubar


  unsplit_mtuples[0] = (char*)"FD_FHC_ger_nonswap.root";
  unsplit_mtuples[1] = (char*)"FD_FHC_ger_nueswap.root";
  unsplit_mtuples[2] = (char*)"FD_FHC_ger_tauswap.root";
  unsplit_mtuples[3] = (char*)"FD_RHC_ger_nonswap.root";
  unsplit_mtuples[4] = (char*)"FD_RHC_ger_nueswap.root";
  unsplit_mtuples[5] = (char*)"FD_RHC_ger_tauswap.root";

  double pot;
  double fPOT;
  TTree* trPot;
  // For each file, check if TH1D* norm exists, and replace (or create) it

  //Number of original files i.e. nonswap, nueswap etc.
  for (int swaptype=0; swaptype<6; swaptype++)
  {
	// there are 8 "selections" because there are 4 beam compenents and 2 selections, i.e. blah_x_blah_numuselec and blah_x_blah_nueslsec
	for (int selec=0; selec<8; selec++)
	{
	  //if(selec>2 && swaptype>5){continue;}
	  std::cout << " on swaptype " << swaptype << std::endl;
	  char FileName[200];
	  char FileName2[200];
	  sprintf(FileName, "%s/%s", mtuple_folder, sk_filenames[swaptype][selec]);
	  sprintf(FileName2, "%s/%s", unsplit_mtuple_folder, unsplit_mtuples[swaptype]);
	  std::cout << "on file " << FileName << std::endl;

	  TFile *sk_file = new TFile(FileName,"UPDATE"); 

	  TH1D *norm_old = (TH1D*)sk_file->Get("norm");
          if(norm_old) // delete histogram from file
	  {
		sk_file->Delete("norm;*");
	  }

	  TH1D *norm = new TH1D("norm","norm",10,0,10);

	  // Set norm bin 1 content (neutrino type normalisation)
	  double val = numerators[swaptype]/denominators[swaptype];
	  norm->SetBinContent(2,val);

	  // Set norm bin 2 content (POT used to generate MC)
	  if (swaptype<8) // Neutrino mode (normal horn current)
	  {
		//selec == 0 is here just so that we don't have to resum all the entries in the same swap file
		if (selec == 0)
		{
		  TFile *unsplit_file = new TFile(FileName2,"UPDATE"); 
		  //if (unsplit_file->GetListOfKeys()->Contains("cafmaker/meta"))
		  //{
			trPot = (TTree*)unsplit_file->Get("cafmaker/meta");
			double pot;
			fPOT = 0;
			trPot->SetBranchAddress("pot", &pot);
			for(int n = 0; n < trPot->GetEntries(); ++n){
			  trPot->GetEntry(n);
			  fPOT += pot;
			}
			unsplit_file->Close();
		  //}
		  //else{std::cout << "[ERROR:] Couldn't find meta tree in file " << FileName2 << std::endl;}
		}
		else  {
		  std::cout << "No need to recalculate fPOT as this file corresponds to the same unsplit file " << std::endl; 
		}
		sk_file->cd();
		norm->SetBinContent(1,1e21/fPOT);
		norm->Write();
	  }
	  else // Antuneutrino mode (reverse horn current)
	  {
		//norm->SetBinContent(2,pot_anu);
	  }

	  // Save norm hist to file
	  sk_file->Close();
	}
  }
}
