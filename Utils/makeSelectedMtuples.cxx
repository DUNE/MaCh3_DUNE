#include <TFile.h>
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "EVGCore/EventRecord.h"
#include "nusystematics/artless/response_helper.hh"
#include <stdio.h>


void Split(std::string infile, int initial_nuPDG, int final_nuPDG, bool numuselec) 
{  

   const char* path  = "/vols/t2k/users/ljw20/software/m3_dune/Inputs/mtuples/GenieCAF/HaddedFiles";
   const char* outpath  = "/vols/t2k/users/ljw20/software/m3_dune/Inputs/mtuples/GenieSampleCuts2";

   int nutype;
   int nutypeunosc;
   double cvn;
   float cvn_threshold;
   
   std::string initial_nu_name; 
   std::string final_nu_name; 
   std::string selec_name;
   std::string outfile_name; 
   std::string cvn_name;

   if (initial_nuPDG == 14) {initial_nu_name = "numu_x_";}    
   if (initial_nuPDG == -14) {initial_nu_name = "numubar_x_";}    
   if (initial_nuPDG == 12) {initial_nu_name = "nue_x_";}    
   if (initial_nuPDG == -12) {initial_nu_name = "nuebar_x_";}    
   if (initial_nuPDG == 16) {initial_nu_name = "nutau_x_";}    
   if (initial_nuPDG == -16) {initial_nu_name = "nutaubar_x_";}    

   if (final_nuPDG == 14) {final_nu_name = "numu_";}    
   if (final_nuPDG == -14) {final_nu_name = "numubar_";}    
   if (final_nuPDG == 12) {final_nu_name = "nue_";}    
   if (final_nuPDG == -12) {final_nu_name = "nuebar_";}    
   if (final_nuPDG == 16) {final_nu_name = "nutau_";}    
   if (final_nuPDG == -16) {final_nu_name = "nutaubar_";}    

   if(numuselec) {
     cvn_name = "cvnnumu";
     cvn_threshold = 0.5; 
     selec_name = "numuselec";}
   else {
     cvn_name = "cvnnue";
     cvn_threshold = 0.85; 
     selec_name = "nueselec";}

   outfile_name = infile.substr(0, 11) + initial_nu_name + final_nu_name + selec_name + ".root";
   
   char filepath[200];
   char outfilepath[200];
   sprintf(filepath, "%s/%s", path, infile.c_str());
   sprintf(outfilepath, "%s/%s", outpath, outfile_name.c_str());


   TFile * tf = new TFile(filepath, "READ");
   TFile * outf = new TFile(outfilepath, "RECREATE");     

   TTree * tree = (TTree*) tf->Get( "cafmaker/caf" );
   TTree * gtree = (TTree*) tf->Get( "cafmaker/genieEvt" );
   TTree * newtree = tree->CloneTree(0);
   TTree * gnewtree = gtree->CloneTree(0);

   tree->SetBranchAddress("nuPDGunosc",&nutypeunosc);
   tree->SetBranchAddress("nuPDG",&nutype);
   tree->SetBranchAddress(cvn_name.c_str(),&cvn);

   int N = tree->GetEntries();
   std::cout << "N = " << N << std::endl; 
   for (int i = 0; i < N; ++i) {
     tree->GetEntry(i);
     // if( i % 100 == 0 ) printf( "Event %d of %d...\n", i, N );
     
     if(nutypeunosc == initial_nuPDG and nutype == final_nuPDG and cvn > cvn_threshold) {
       newtree->Fill();
       gtree->GetEntry(i);
       gnewtree->Fill(); }  
     
     }

   int n_selection =  newtree->GetEntries();
   int gn_selection =  gnewtree->GetEntries();
   std::cout << "Check: number of events in " << outfile_name << " is: " << n_selection << std::endl;
   std::cout << "Check: number of events in GENIE TREE is: " << gn_selection << std::endl;
   
   if (n_selection != 0) {
     newtree->Write("caf");
     gnewtree->Write("genieEvt");
     std::cout << "Writing output ROOT File" << std::endl;
     outf->Write();
     outf->Close(); }
   
   else {
     std::cout << "No events inside " << outfile_name << " , so no output ROOT file " << std::endl;
     remove(outfilepath);
   }


}

int main( int argc, char const *argv[] )
{

  if( (argc == 2) && ((std::string("--help") == argv[1]) || (std::string("-h") == argv[1])) ) {
    std::cout << "Help yourself by looking at the source code to see what the options are." << std::endl;
    return 0;
  }

  std::string infile;
  infile = argv[1];

  
  // NUMU_X_NUMU
  Split(infile, 14, 14, true); //numu_x_numu_numuselec
  Split(infile, 14, 14, false); //numu_x_numu_nueselec
  Split(infile, -14, -14, true); //numubar_x_numubar_numuselec
  Split(infile, -14, -14, false); //numubar_x_numubar_nueselec
  
  // NUMU_X_NUE
  Split(infile, 14, 12, true); //numu_x_nue_numuselec
  Split(infile, 14, 12, false); //numu_x_nue_nueselec
  Split(infile, -14, -12, true); //numubar_x_nuebar_numuselec
  Split(infile, -14, -12, false); //numubar_x_nuebar_nueselec
  
  // NUMU_X_NUTAU
  Split(infile, 14, 16, true); //numu_x_nutau_numuselec
  Split(infile, 14, 16, false); //numu_x_nutau_nueselec
  Split(infile, -14, -16, true); //numubar_x_nutaubar_numuselec
  Split(infile, -14, -16, false); //numubar_x_nutaubar_nueselec
  
  // NUE_X_NUE
  Split(infile, 12, 12, true); //nue_x_nue_numuselec 
  Split(infile, 12, 12, false); //nue_x_nue_nueselec
  Split(infile, -12, -12, true); //nuebar_x_nuebar_numuselec
  Split(infile, -12, -12, false); //nuebar_x_nuebar_nueselec
  
  // NUE_X_NUMU
  Split(infile, 12, 14, true); //nue_x_numu_numuselec
  Split(infile, 12, 14, false); //nue_x_numu_nueselec
  Split(infile, -12, -14, true); //nuebar_x_numubar_numuselec
  Split(infile, -12, -14, false); //nuebar_x_numubar_nueselec
  
  //NUE_X_NUTAU
  Split(infile, 12, 16, true); //nue_x_nutau_numuselec
  Split(infile, 12, 16, false); //nue_x_nutau_nueselec
  Split(infile, -12, -16, true); //nuebar_x_nutaubar_numuselec
  Split(infile, -12, -16, false); //nuebar_x_nutaubar_nueselec

}




  
