#include <iostream>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TNamed.h"

int main(int argc, char *argv[])
{
  if (argc != 3 && argc != 4)
    {
      std::cout << "usage: reduceDataSet \"/some/files/blah*.root\" outputfile.root \"some information string\"" << std::endl;
      std::cout << "note: the last argument is optional, enter a sentence to describe whats in the file." << std::endl;
      return 1;
    }

  double t_sin2th_23;
  double t_sin2th_13;
  double t_delm2_23;
  double t_delta_cp;
  double t_LogL;
  int t_step;
  TNamed info;

  if (argc == 4)
    info = TNamed("info", argv[3]);

  TChain* myChain = new TChain("posteriors","");
  myChain->Add(argv[1]);

  myChain->SetBranchStatus("*",0);
  myChain->SetBranchStatus("sin2th_23",1);
  myChain->SetBranchStatus("sin2th_13",1);
  myChain->SetBranchStatus("delm2_23",1);
  myChain->SetBranchStatus("delta_cp",1);
  myChain->SetBranchStatus("LogL",1);
  myChain->SetBranchStatus("step",1);
  
  myChain->SetBranchAddress("sin2th_23", &t_sin2th_23);
  myChain->SetBranchAddress("sin2th_13", &t_sin2th_13);
  myChain->SetBranchAddress("delm2_23", &t_delm2_23);
  myChain->SetBranchAddress("delta_cp", &t_delta_cp);
  myChain->SetBranchAddress("LogL", &t_LogL);
  myChain->SetBranchAddress("step", &t_step);
  
  TFile *newfile = TFile::Open(argv[2],"recreate");
  TTree *newtree = new TTree("osc_posteriors","oscillation parameter posteriors, marginalized");

  newtree->Branch("theta23", &t_sin2th_23, "theta23/D");
  newtree->Branch("theta13", &t_sin2th_13, "theta13/D");
  newtree->Branch("dm23", &t_delm2_23, "dm23/D");
  newtree->Branch("dcp", &t_delta_cp, "dcp/D");
  newtree->Branch("step", &t_step, "step/I");
  newtree->Branch("LogL", &t_LogL, "LogL/D");

  for (int i = 0; i < myChain->GetEntries(); ++i)
    {
      if (i % 10000 == 0)
	std::cout << i << std::endl;
      myChain->GetEntry(i);
      newtree->Fill();
    }

  if (argc == 4)
    info.Write();

  newtree->Print();
  newtree->AutoSave();
  
  delete newfile;

 return 0;
}
