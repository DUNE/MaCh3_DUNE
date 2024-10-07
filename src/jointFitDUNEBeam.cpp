#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>

#include <TH1D.h>
#include <THStack.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TRint.h>
#include <TLegend.h>
#include <TColor.h>
#include <TMath.h>

#include "mcmc/mcmc.h"
#include "samplePDFDUNE/MaCh3DUNEFactory.h"


std::string getNameNoExt(std::string name, std::string ext)  
{                                                            
  std::size_t pos ;                                          
  pos = name.find(ext);                                      
  name = name.substr(0,pos);                                 
  return name ;                                              
}                                                            
                                                             
void saveCanvas(TCanvas* canvas, std::string name, std::string legend)                                                                  
{                                                            
  name = getNameNoExt(name, ".root") ;                       
  name = name + legend + ".root" ;                     
  canvas -> SaveAs(name.c_str()) ;                           
                                                             
  name = getNameNoExt(name, ".root") ;                       
  name = name + ".png" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           
                                                             
  name = getNameNoExt(name, ".png") ;                        
  name = name + ".pdf" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           
                                                             
  name = getNameNoExt(name, ".pdf") ;                        
  name = name + ".eps" ;                                     
  canvas -> SaveAs(name.c_str()) ;                           
                                                             
} 

void DrawModeBreakDown(samplePDFDUNEBeamFDBase*& Sample, std::unique_ptr<TFile> &File){

  File->cd();
  TH1D* ModeHistogram = nullptr;
  THStack* StackedModes = new THStack(Sample->GetName().c_str(), "stack");
  TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);

  std::ofstream outfile;
  std::string FileName=Sample->GetName()+"_EventRates.txt";
  outfile.open(FileName);

  outfile << "\\begin{table}[ht]" << std::endl;
  outfile << "\\begin{center}" << std::endl;
  outfile << "\\caption{Integral breakdown for sample: " << Sample->GetName() << "}" << std::endl;
  outfile << "\\label{" << Sample->GetName() << "-EventRate}" << std::endl;


  int space = 14;
  TString nColumns;
  nColumns += "|c|";
  outfile << "\\begin{tabular}{|l" << nColumns.Data() << "}" << std::endl;
  outfile << "\\hline" << std::endl;
  outfile << "&" << std::setw(space) << Sample->GetName() << " " << std::endl;

  double total = 0;
  double mode_rate = 0;
  for (unsigned int iMode=0;iMode<kMaCh3_nModes ;iMode++) {

	////
	std::string name = Sample->GetName();
	name += MaCh3mode_ToDUNEString(static_cast<MaCh3_Mode>(iMode));
	ModeHistogram = (TH1D*)Sample->get1DVarHist(kRecoNeutrinoEnergy, iMode);
	ModeHistogram->SetName(name.c_str());
	mode_rate = ModeHistogram->Integral();
	total += mode_rate; 

	outfile << MaCh3mode_ToDUNEString(static_cast<MaCh3_Mode>(iMode)).c_str() << std::setw(space) << " & " << mode_rate << " \\\\" << std::endl;

	//Add this to the stacked histogram
	if(ModeHistogram->Integral() > 0){
	  ModeHistogram->SetFillColor(kRed+iMode);
	  ModeHistogram->SetLineColor(kRed+iMode);
	  StackedModes->Add(ModeHistogram);	
	  leg->AddEntry(ModeHistogram, MaCh3mode_ToDUNEString(static_cast<MaCh3_Mode>(iMode)).c_str(), "f");
	  ModeHistogram->Write();
	}
  }
  outfile << "\\hline " << std::endl;
  outfile << "&" << std::setw(space) << "Total:" << total << "\\\\ \\hline" << std::endl;

  std::string PdfName = Sample->GetName();
  PdfName += "_ModeStack.pdf";
  StackedModes->Write();
  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();
  StackedModes->Draw("HIST");
  leg->Draw("SAMES");
  c1->Print(PdfName.c_str());
  delete c1;
  

  return;
}

int main(int argc, char * argv[]) {

  // ----------------------- OPTIONS ---------------------------------------- //
  if(argc == 1){
    std::cout << "Usage: bin/jointFitDUNEBeam configs/config.yaml" << std::endl;
    return 1;
  }

  manager *FitManager = new manager(argv[1]);
  auto OutputFileName = FitManager->raw()["General"]["OutputFile"].as<std::string>();

  covarianceXsec* xsec = nullptr;
  covarianceOsc* osc = nullptr;

  //####################################################################################
  //Create samplePDFSKBase Objs
  std::cout << "Loading T2K samples.." << "\n" << std::endl;
  std::vector<samplePDFDUNEBeamFDBase*> DUNEPdfs;
  MakeMaCh3DuneBeamInstance(FitManager, DUNEPdfs, xsec, osc); 
  //Setup the cross-section parameters
  //This should get the prior values.
  std::vector<double> XsecParVals = xsec->getNominalArray();

  xsec->setParameters(XsecParVals);
  xsec->setStepScale(FitManager->raw()["General"]["Systematics"]["XsecStepScale"].as<double>());

  //Some place to store the histograms
  std::vector<TH1D*> PredictionHistograms;
  std::vector<std::string> sample_names;

  auto OutputFile = std::unique_ptr<TFile>(new TFile(OutputFileName.c_str(), "RECREATE"));
  OutputFile->cd();

  for (unsigned sample_i = 0 ; sample_i < DUNEPdfs.size() ; ++sample_i) {
    	
	std::string name = DUNEPdfs[sample_i]->GetName();
	sample_names.push_back(name);
	TString NameTString = TString(name.c_str());

	DUNEPdfs[sample_i] -> SetupOscCalc(osc->GetPathLength(), osc->GetDensity());
	osc->setParameters();
	DUNEPdfs[sample_i] -> reweight();
	TH1D *SampleHistogram = (TH1D*)DUNEPdfs[sample_i] -> get1DHist() -> Clone(NameTString+"_unosc");
	PredictionHistograms.push_back(SampleHistogram);

	DUNEPdfs[sample_i]->addData(PredictionHistograms[sample_i]);
  }
 
  //Now print out some event rates, we'll make a nice latex table at some point 
  for (unsigned iPDF = 0; iPDF < DUNEPdfs.size() ; ++iPDF) {
	std::cout << "Integrals of nominal hists: " << std::endl;
	std::cout << sample_names[iPDF].c_str() << ": " << PredictionHistograms[iPDF]-> Integral() << std::endl;
	std::cout << "~~~~~~~~~~~~~~~~" << std::endl;
  }

  //###########################################################################################################
  // Set covariance objects equal to output of previous chain
  
  std::vector<double> oscparstarts;
  std::map<TString,std::vector<double> > parstarts;
  int lastStep = 0;

  bool StartFromPreviousChain = GetFromManager(FitManager->raw()["General"]["StartFromPos"], false);

  if(StartFromPreviousChain) {//Start from values at the end of an already run chain
    //Read in paramter names and values from file
    std::cout << "MCMC getting starting position from " << FitManager->raw()["General"]["PosFileName"].as<std::string>() << std::endl;
    TFile *infile = new TFile(FitManager->raw()["General"]["PosFileName"].as<std::string>().c_str(), "READ");
    TTree *posts = (TTree*)infile->Get("posteriors");
    TObjArray* brlis = (TObjArray*)posts->GetListOfBranches();
    int nbr = brlis->GetEntries();
    TString branch_names[nbr];
    double branch_vals[nbr];
    int step_val;
    for (int i = 0; i < nbr; ++i) {
      TBranch *br = (TBranch*)brlis->At(i);
      TString bname = br->GetName();
      branch_names[i] = bname;
      if(branch_names[i] == "step") {
        posts->SetBranchAddress("step",&step_val);
        continue;
      }
      std::cout << " * Loading " << bname << std::endl;
      posts->SetBranchAddress(branch_names[i], &branch_vals[i]);
    }
    posts->GetEntry(posts->GetEntries()-1);
    std::map<TString,double> init_pars;
    for (int i = 0; i < nbr; ++i) {
      init_pars.insert( std::pair<TString, double>(branch_names[i], branch_vals[i]));
    }
    infile->Close();
    delete infile;
    
    //Make vectors of parameter value for each covariance
    std::vector<TString> covtypes;
    covtypes.push_back("xsec");
 
    for(unsigned icov=0;icov<covtypes.size();icov++){
      std::vector<double> covparstarts;
      std::map<TString, double>::const_iterator it;
      int iPar=0;
      while(it != init_pars.end()){
  	it = init_pars.find(covtypes[icov]+"_"+TString::Format("%d",iPar));
  	if (it != init_pars.end()) {
  	  covparstarts.push_back(it->second);
  	}
  	iPar++;
      }
      if(covparstarts.size()!=0) parstarts.insert(std::pair<TString,std::vector<double> >(covtypes[icov],covparstarts));
      else std::cout<<"Did not find any parameters in posterior tree for: "<<covtypes[icov]<<std::endl<<"assuming previous chain didn't use them"<<std::endl;
    }

    std::map<TString, double>::const_iterator itt;

    // set the oscillation parameters
    itt = init_pars.find("sin2th_12");
    oscparstarts.push_back(itt->second);
    itt = init_pars.find("sin2th_23");
    oscparstarts.push_back(itt->second);
    itt = init_pars.find("sin2th_13");
    oscparstarts.push_back(itt->second);
    itt = init_pars.find("delm2_12");
    oscparstarts.push_back(itt->second);
    itt = init_pars.find("delm2_23");
    oscparstarts.push_back(itt->second);
    itt = init_pars.find("delta_cp");
    oscparstarts.push_back(itt->second);

    lastStep = step_val;

  }

  //###########################################################################################################
  // Back to actual nominal for fit
  // If starting from end values of a previous chain set everything to the values from there
  // and do acceptStep() to update fParCurr with these values
  //
  if(GetFromManager(FitManager->raw()["General"]["StartFromPos"], false)) {
    if(parstarts.find("xsec")!=parstarts.end()) {
      xsec->setParameters(parstarts["xsec"]);
      xsec->acceptStep();
    }
    else {xsec->setParameters();}
  }
  else {
    xsec->setParameters();
  }

  //###########################################################################################################
  // MCMC

  //mcmc *MaCh3Fitter = new mcmc(FitManager);

  std::unique_ptr<FitterBase> MaCh3Fitter = std::make_unique<mcmc>(FitManager);

  //if(lastStep > 0) MaCh3Fitter->setInitialStepNumber(lastStep+1);

  // add samples
  for(auto Sample : DUNEPdfs){
	MaCh3Fitter->addSamplePDF(Sample);
  }

  //start chain from random position
  xsec->throwParameters();
  osc->throwParameters();

  // add systematic objects
  if (GetFromManager(FitManager->raw()["General"]["StatOnly"], false)){
    MaCh3Fitter->addSystObj(osc);
	MACH3LOG_INFO("Running a stat-only fit so no systematics will be applied");
  }
  else {
    MaCh3Fitter->addSystObj(osc);
    MaCh3Fitter->addSystObj(xsec);
  }
  
  // run!
  MaCh3Fitter->runMCMC();


  return 0;
 }