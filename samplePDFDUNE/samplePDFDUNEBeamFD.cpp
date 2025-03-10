#include "samplePDFDUNEBeamFD.h"

#include "manager/manager.h"

#include "splines/splinesDUNE.h"

#include "duneanaobj/StandardRecord/StandardRecord.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"


samplePDFDUNEBeamFD::samplePDFDUNEBeamFD(std::string mc_version_,
                                         covarianceXsec *xsec_cov_)
    : samplePDFFDBase(mc_version_, xsec_cov_) {
  // Call insitialise in samplePDFFD
  Initialise();
}

samplePDFDUNEBeamFD::~samplePDFDUNEBeamFD() {}

void samplePDFDUNEBeamFD::Init() {
  dunemcSamples.resize(nSamples, dunemc_base()); 
  pot = SampleManager->raw()["POT"].as<double>();
  gen_pot = SampleManager->raw()["GEN_POT"].as<double>();

  numu_cut = SampleManager->raw()["numu_cut"].as<double>();
  nue_cut = SampleManager->raw()["nue_cut"].as<double>();
  std::cout<< "POT IS !!! ================ " <<pot <<std::endl;
  std::cout<< "numu_cut IS !!! ================ " << numu_cut <<std::endl;
  std::cout<< "nue_cut IS !!! ================ " << nue_cut <<std::endl;
 
    

  MACH3LOG_INFO(
      "-------------------------------------------------------------------");
}


double samplePDFDUNEBeamFD::CalculatePOT() {
  TChain calc_pot_chain("cafmaker/meta");  // Use correct tree name

  std::string pot_branch = "pot";  // Use the correct branch name

  for (int i = 0; i < (int)dunemcSamples.size(); ++i) {
      calc_pot_chain.AddFile((mc_files[i]).c_str());
  }

  // Check if the branch exists before proceeding
  if (!calc_pot_chain.GetBranch(pot_branch.c_str())) {
      std::cerr << "Error: Branch " << pot_branch << " not found in the tree!" << std::endl;
      return 0.0;
  }

  double pot_value = 0.0;
  calc_pot_chain.SetBranchAddress(pot_branch.c_str(), &pot_value);

  double sum_pot = 0.0;
  Long64_t nEntries = calc_pot_chain.GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
      calc_pot_chain.GetEntry(i);
      sum_pot += pot_value;
  }

  std::cout << "Summed POT: " << sum_pot << std::endl;
  return sum_pot;
}

/*
void samplePDFDUNEBeamFD::SetupSplines() {
  ///@todo move all of the spline setup into core
  if (XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline) > 0) {
    MACH3LOG_INFO(
        "Found {} splines for this sample so I will create a spline object",
        XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    splinesDUNE *DUNESplines = new splinesDUNE(XsecCov);
    splineFile = (splineFDBase *)DUNESplines;
    InitialiseSplineObject();
  } else {
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or "
                  "evaluate splines",
                  XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    splineFile = nullptr;
  }
}*/


void samplePDFDUNEBeamFD::SetupSplines() {
  ///@todo move all of the spline setup into core
  if (XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline) > 0) {
    MACH3LOG_INFO(
        "Found {} splines for this sample so I will create a spline object",
        XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    splinesDUNE *DUNESplines = new splinesDUNE(XsecCov);
    splineFile = (splineFDBase *)DUNESplines;
    InitialiseSplineObject();
  } else {
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or "
                  "evaluate splines",
                  XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    splineFile = nullptr;
  }
}


void samplePDFDUNEBeamFD::SetupWeightPointers() {
  for (int i = 0; i < (int)dunemcSamples.size(); ++i) {
    for (int ev_i = 0; ev_i < dunemcSamples[i].nEvents; ++ev_i) {
      MCSamples[i].ntotal_weight_pointers[ev_i] = 6;
      MCSamples[i].total_weight_pointers[ev_i] =
          new const double *[MCSamples[i].ntotal_weight_pointers[ev_i]];
      MCSamples[i].total_weight_pointers[ev_i][0] = &(dunemcSamples[i].pot_s);
      MCSamples[i].total_weight_pointers[ev_i][1] = &(dunemcSamples[i].norm_s);
      MCSamples[i].total_weight_pointers[ev_i][2] =
          MCSamples[i].osc_w_pointer[ev_i];
      MCSamples[i].total_weight_pointers[ev_i][3] =
          &(dunemcSamples[i].rw_berpaacvwgt[ev_i]);
      MCSamples[i].total_weight_pointers[ev_i][4] =
          &(dunemcSamples[i].flux_w[ev_i]);
      MCSamples[i].total_weight_pointers[ev_i][5] =
          &(MCSamples[i].xsec_w[ev_i]);
    }
  }
}

double summed_pot = 0.0; //set the sum of the pot from each file to be 0,befor any are read in



int samplePDFDUNEBeamFD::setupExperimentMC(int iSample) {
  std::cout << "Summed POT, after sample "<< iSample << "is loaded = " << summed_pot << std::endl;
   
  double production_pot = 1.3628319e+23;

  dunemc_base *duneobj = &(dunemcSamples[iSample]);

  MACH3LOG_INFO(
      "-------------------------------------------------------------------");
  MACH3LOG_INFO("input file: {}", mc_files[iSample]);

  std::unique_ptr<TFile> in_file{
      TFile::Open(mc_files[iSample].c_str(), "READ")};
  if (!in_file) {
    MACH3LOG_ERROR("Could not open file: {}", mc_files[iSample]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  auto in_tree = in_file->Get<TTree>("cafmaker/cafTree");
  if (!in_tree) {
    MACH3LOG_ERROR("Could not find \"cafmaker/cafTree\" TTree in {}",
                   mc_files[iSample]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
 auto meta_tree = in_file->Get<TTree>("cafmaker/meta");
  if (!meta_tree) {
    MACH3LOG_ERROR("Could not find \"cafmaker/meta\" TTree in {}",
                   mc_files[iSample]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  TTreeReader tree_rdr(in_tree);
  TTreeReader tree_rdrmeta(meta_tree);
  // The branch "px" contains floats; access them as myPx.
  TTreeReaderValue<caf::StandardRecord> sr(tree_rdr, "rec");
  TTreeReaderValue<caf::StandardRecord> meta(tree_rdrmeta, "meta");

  //double analysis_POT = fitMan->raw()["General"]["POT"].as<double>(); 
  duneobj->nEvents = tree_rdr.GetEntries();

  MACH3LOG_INFO("Found \"cafmaker/cafTree\" TTree in {}", mc_files[iSample]);
  MACH3LOG_INFO("With number of entries: {}", duneobj->nEvents);

 
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_erec_shifted = new double[duneobj->nEvents];
  duneobj->rw_vtx_x = new double[duneobj->nEvents];
  duneobj->rw_vtx_y = new double[duneobj->nEvents];
  duneobj->rw_vtx_z = new double[duneobj->nEvents];
  duneobj->rw_cvnnumu_shifted = new double[duneobj->nEvents];
  duneobj->rw_cvnnue_shifted = new double[duneobj->nEvents];
  duneobj->mode = new double[duneobj->nEvents];
  duneobj->Target = new int[duneobj->nEvents];
  duneobj->rw_isCCevent = new double[duneobj->nEvents];

  duneobj->rw_isnumu = new int[duneobj->nEvents];
  duneobj->rw_isnue = new int[duneobj->nEvents];
  duneobj->mc_isnumu = new double[duneobj->nEvents];
  duneobj->mc_isnue = new double[duneobj->nEvents];
  duneobj->mc_isnutau = new double[duneobj->nEvents];

  duneobj->numuefficency =  new double[duneobj->nEvents];
  duneobj->nueefficency =  new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->selected_numuCCevent_energy  =  new double[duneobj->nEvents];
  duneobj->selected_numuCCevent_vertexpos_x  =  new double[duneobj->nEvents];
  duneobj->selected_numuCCevent_vertexpos_y=  new double[duneobj->nEvents];
  duneobj->selected_numuCCevent_vertexpos_z  =  new double[duneobj->nEvents];

  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents];

  duneobj->nu_pdg = new double[duneobj->nEvents];

  duneobj->selected_nueCCevent_energy  =  new double[duneobj->nEvents];
  duneobj->selected_nueCCevent_vertexpos_x  =  new double[duneobj->nEvents];
  duneobj->selected_nueCCevent_vertexpos_y=  new double[duneobj->nEvents];
  duneobj->selected_nueCCevent_vertexpos_z  =  new double[duneobj->nEvents];
  duneobj->nuflavour = new double[duneobj->nEvents];
  
  TTreeReader metatree_rdr(in_file->Get<TTree>("cafmaker/meta"));
  TTreeReaderValue<double> pot(metatree_rdr, "pot");
  double total_pot = 0;
  for (auto entryi : metatree_rdr) {
     total_pot += *pot;
  }
  std::cout << "POT in file =  " << total_pot << std::endl;

  summed_pot = summed_pot + total_pot;

  std::cout << "Updated summed POT =  " << summed_pot << std::endl;
 
  duneobj->nutype = sample_nutype[iSample];
  duneobj->oscnutype = sample_oscnutype[iSample];
  duneobj->signal = sample_signal[iSample];
  
 
  // after allocation, set all event entries to a bad bad value
  std::fill_n(duneobj->rw_etru, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_erec_shifted, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_vtx_x, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_vtx_y, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_vtx_z, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_cvnnumu_shifted, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_cvnnue_shifted, duneobj->nEvents, -99999);
  std::fill_n(duneobj->mode, duneobj->nEvents, -99999);
  std::fill_n(duneobj->Target, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_isCCevent, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_isnumu, duneobj->nEvents, -99999);
  std::fill_n(duneobj->rw_isnue, duneobj->nEvents, -99999);
  std::fill_n(duneobj->mc_isnumu, duneobj->nEvents, -99999);
  std::fill_n(duneobj->mc_isnue, duneobj->nEvents, -99999);
  std::fill_n(duneobj->mc_isnutau, duneobj->nEvents, -99999);
  std::fill_n(duneobj->numuefficency, duneobj->nEvents, -99999);
  std::fill_n(duneobj->nueefficency, duneobj->nEvents, -99999);
  std::fill_n(duneobj->selected_numuCCevent_energy , duneobj->nEvents, -99999);
  std::fill_n(duneobj->selected_numuCCevent_vertexpos_x , duneobj->nEvents, -99999);
  std::fill_n(duneobj->selected_numuCCevent_vertexpos_y , duneobj->nEvents, -99999);
  std::fill_n(duneobj->selected_numuCCevent_vertexpos_z , duneobj->nEvents, -99999);
  std::fill_n(duneobj->nuflavour , duneobj->nEvents, -99999);
  std::fill_n(duneobj->nu_pdg, duneobj->nEvents, -99999);
  
  std::fill_n(duneobj->selected_nueCCevent_energy , duneobj->nEvents, -99999);
  std::fill_n(duneobj->selected_nueCCevent_vertexpos_x , duneobj->nEvents, -99999);
  std::fill_n(duneobj->selected_nueCCevent_vertexpos_y , duneobj->nEvents, -99999);
  std::fill_n(duneobj->selected_nueCCevent_vertexpos_z , duneobj->nEvents, -99999);
  

  //TH1D numerator("numerator_isnumuCC_isselnumuCC",";E_{#nu};Rate (IsSelNumuCC && IsTrueNumuCC)",100,0,10);
  //TH1D denominator("denominator_isnumuCC",";E_{#nu};Rate IsSelNumuCC",100,0,10);

   
    double newpot = CalculatePOT();

  for (auto entryi : tree_rdr) {

    if (sr->mc.nu.size() < 1) { // no true neutrinos here
      continue;
    }

    //bool isnumuCC = duneobj->mc_isnumu[entryi];
    //if(!isnumuCC) { continue; } // don't need to worry about this event

    //bool issel_numuCC = (sr->common.ixn.pandora[0].nuhyp.cvn.numu > numu_cut);

    /*if(issel_numuCC){
      numerator.Fill(duneobj->enutrue[entryi]);
    }*/
    if (sr->mc.nu[0].pdg == 14){
    duneobj->mc_isnumu[entryi]==sr->mc.nu[0].pdg;
    };

    if (sr->mc.nu[0].pdg == 12){
    duneobj->mc_isnue[entryi]=sr->mc.nu[0].pdg;
    };

    if (sr->mc.nu[0].pdg == 16){
    duneobj->mc_isnutau[entryi]=sr->mc.nu[0].pdg;
    };

    duneobj->nuflavour[entryi]=sr->mc.nu[0].pdg;
    duneobj->rw_etru[entryi] = sr->mc.nu[0].E;
    duneobj->mode[entryi] = sr->mc.nu[0].mode;
    duneobj->Target[entryi] = (sr->mc.nu[0].targetPDG / 1000) % 1000;
    duneobj->rw_isCCevent[entryi] = sr->mc.nu[0].iscc;

    duneobj->flux_w[entryi] = 1.0;

    

    //std::cout<< "rw_etru ==  " <<  sr->mc.nu[0].E <<  std::endl;
    //std::cout<< "rw_isCCevent ==  " <<   sr->mc.nu[0].iscc <<  std::endl;

    /* 
    if (sr->mc.nu[0].pdg == 14) {
        duneobj->mc_isnumu[entryi] = sr->mc.nu[0].pdgorig;
    }
    if (sr->mc.nu[0].pdg == 12) {
        duneobj->mc_isnue[entryi] = sr->mc.nu[0].pdgorig;
    }*/
    
    
    if (sr->common.ixn.pandora.size() < 1) { // no reconstructed objects
      continue;
    }
    duneobj->nu_pdg[entryi] = sr->mc.nu[0].pdg;////////////////////////////pdg code of reco. particle
    duneobj->rw_erec_shifted[entryi] = sr->common.ixn.pandora[0].Enu.lep_calo;
    duneobj->rw_vtx_x[entryi] = sr->common.ixn.pandora[0].vtx.x;
    duneobj->rw_vtx_y[entryi] = sr->common.ixn.pandora[0].vtx.y;
    duneobj->rw_vtx_z[entryi] = sr->common.ixn.pandora[0].vtx.z;
    duneobj->rw_cvnnumu_shifted[entryi] =
        sr->common.ixn.pandora[0].nuhyp.cvn.numu;
    duneobj->rw_cvnnue_shifted[entryi] =
        sr->common.ixn.pandora[0].nuhyp.cvn.nue;                
    //std::cout << "sr->common.ixn.pandora[0].nuhyp.cvn.numu; = " << sr->common.ixn.pandora[0].nuhyp.cvn.numu <<std::endl;
    if (sr->mc.nu[0].pdg == 14  && sr->mc.nu[0].iscc ==1 ) { 
        duneobj->selected_numuCCevent_energy[entryi] = sr->mc.nu[0].E;
        duneobj->selected_numuCCevent_vertexpos_x[entryi] = sr->mc.nu[0].vtx.x;
        duneobj->selected_numuCCevent_vertexpos_y[entryi] = sr->mc.nu[0].vtx.y;
        duneobj->selected_numuCCevent_vertexpos_z[entryi] = sr->mc.nu[0].vtx.z;
    }
    if (sr->mc.nu[0].pdg == 12  && sr->mc.nu[0].iscc ==1 ) { 
        duneobj->selected_nueCCevent_energy[entryi] = sr->mc.nu[0].E;
        duneobj->selected_nueCCevent_vertexpos_x[entryi] = sr->mc.nu[0].vtx.x;
        duneobj->selected_nueCCevent_vertexpos_y[entryi] = sr->mc.nu[0].vtx.y;
        duneobj->selected_nueCCevent_vertexpos_z[entryi] = sr->mc.nu[0].vtx.z;
    }
     duneobj->rw_berpaacvwgt[entryi] = 1.0;

    //std::cout << "DEBUG: numu_cut = " << numu_cut << std::endl;
    
    /*
     if (sr->mc.nu[0].pdg == 12 && sr->mc.nu[0].iscc==1 && sr->common.ixn.pandora[0].nuhyp.cvn.nue  > nue_cut) { //////for efficiency/purity
        eventsthatpasscut += 1;
        //std::cout << "eventsthatpasscut = " << eventsthatpasscut << std::endl;
    }*/
    
    /*
    if (sr->mc.nu[0].pdg == 12 && sr->mc.nu[0].iscc==1 ) { //////for efficiency/purity
        total_true_ccnumu  += 1;
        //std::cout << "total_true_ccnumu  = " << total_true_ccnumu  << std::endl;
    }

    
    
    if (sr->common.ixn.pandora[0].nuhyp.cvn.nue  > nue_cut) { //////for efficiency/purity
        total_events_incut  += 1;
        //std::cout << "total_events_incut  = " << total_events_incut  << std::endl;
    }

    if (sr->mc.nu[0].pdg == 12 && sr->mc.nu[0].iscc == 1 && sr->common.ixn.pandora[0].nuhyp.cvn.nue > nue_cut) {
      eventsthatpasscut += 1;
      std::cout << "Passed event: PDG=" << sr->mc.nu[0].pdg 
                << ", isCC=" << sr->mc.nu[0].iscc 
                << ", nue_score=" << sr->common.ixn.pandora[0].nuhyp.cvn.nue 
                << " > " << nue_cut << std::endl;
  }*/

  } //end of for loop
  duneobj->norm_s = 1.0/newpot;
  duneobj->pot_s = *pot; //new double[duneobj->nEvents];
  return duneobj->nEvents;
}

double samplePDFDUNEBeamFD::getEfficiency_nue() {
  /*if (mc_true_total == 0) {
      std::cerr << "Error: Input cannot be zero." << std::endl;
      return -1; // Return -1 to indicate an error
  }*/
  int no_samples = (int)dunemcSamples.size();

  double eventsthatpasscut = 0;
  double total_true_ccnumu = 0;
  double total_events_incut = 0; // Declare the missing variable

  for (int sample = 0; sample < no_samples; ++sample) {
      for (int event = 0; event < dunemcSamples[sample].nEvents; ++event) {
          if (dunemcSamples[sample].nu_pdg[event] == 12 && dunemcSamples[sample].rw_isCCevent[event] == 1) {
              total_true_ccnumu += *(MCSamples[sample].osc_w_pointer[event]); //+= 1;
              //std::cout<<"pdg of event in nu e efficiency cut = " << dunemcSamples[sample].nu_pdg[event] << std::endl;
          }
          if (dunemcSamples[sample].rw_cvnnue_shifted[event] > nue_cut) {
              total_events_incut += *(MCSamples[sample].osc_w_pointer[event]); //+= 1;
          }
          if (dunemcSamples[sample].nu_pdg[event] == 12 && dunemcSamples[sample].rw_isCCevent[event] == 1 &&
              dunemcSamples[sample].rw_cvnnue_shifted[event] > nue_cut) {
               // std::cout<<"pdg of event in nu e efficiency cut = " << dunemcSamples[sample].nu_pdg[event] << std::endl;
              eventsthatpasscut += *(MCSamples[sample].osc_w_pointer[event]);// += 1;
          }
      }
  }
  std::cout << "total_true_ccnumu" << total_true_ccnumu << std::endl;
  std::cout << "total_events_incut" << total_events_incut << std::endl;
  std::cout << "eventsthatpasscut" << eventsthatpasscut << std::endl;

  if (total_true_ccnumu == 0) { // Prevent division by zero
      std::cerr << "Error: No CC νμ events found." << std::endl;
      return -1;
  }

  double efficiency = eventsthatpasscut / total_true_ccnumu;
  std::cout<< "efficiency = " << efficiency << std::endl;
  return efficiency;
}

double samplePDFDUNEBeamFD::getCut(){
  double cut = nue_cut;
  return cut;
}


double samplePDFDUNEBeamFD::getPurity_nue() {
  /*if (events_incut == 0) {  // Ensure 'events_incut' is defined elsewhere
      std::cerr << "Error: Input energy cannot be zero." << std::endl;
      return -1; // Return -1 to indicate an error
  }*/

  double total_events_incut_purity = 0;
  double eventsthatpasscut_purity = 0;

  int no_samples = (int)dunemcSamples.size();

  for (int sample = 0; sample < no_samples; ++sample) {
      for (int event = 0; event < dunemcSamples[sample].nEvents; ++event) {
          if (dunemcSamples[sample].rw_cvnnue_shifted[event] > nue_cut) { 
              total_events_incut_purity += *(MCSamples[sample].osc_w_pointer[event]);//+= 1;
              
          }
          if (dunemcSamples[sample].nu_pdg[event] == 12 && dunemcSamples[sample].rw_isCCevent[event] == 1 &&
            dunemcSamples[sample].rw_cvnnue_shifted[event] > nue_cut) {
            eventsthatpasscut_purity += *(MCSamples[sample].osc_w_pointer[event]);//+= 1;
            //std::cout<<"pdg of event in nue purity cut = " << dunemcSamples[sample].nu_pdg[event] << std::endl;
        }
      }
  }

  if (total_events_incut_purity == 0) {  // Prevent division by zero
      std::cerr << "Error: No events in cut found." << std::endl;
      return -1;
  }
  std::cout << "eventsthatpasscut_purity" << eventsthatpasscut_purity << std::endl;
  std::cout << "total_events_incut_purity" << total_events_incut_purity << std::endl;

  double purity = eventsthatpasscut_purity / total_events_incut_purity;
  std::cout<< "purity = " << purity << std::endl;
  return purity;
}


double samplePDFDUNEBeamFD::getEfficiency_numu() {
 
  int no_samples = (int)dunemcSamples.size();

  double eventsthatpasscut = 0;
  double total_true_ccnumu = 0;
  double total_events_incut = 0; // Declare the missing variable

  for (int sample = 0; sample < no_samples; ++sample) {
      for (int event = 0; event < dunemcSamples[sample].nEvents; ++event) {
          if (dunemcSamples[sample].nu_pdg[event] == 14 && dunemcSamples[sample].rw_isCCevent[event] == 1) {
            //std::cout<<"pdg of event in numu efficiency cut = " << dunemcSamples[sample].nu_pdg[event] << std::endl;
              total_true_ccnumu += *(MCSamples[sample].osc_w_pointer[event]);//+= 1;
          }
          if (dunemcSamples[sample].rw_cvnnumu_shifted[event] > numu_cut) {
              total_events_incut += *(MCSamples[sample].osc_w_pointer[event]);// += 1;
          }
          if (dunemcSamples[sample].nu_pdg[event] == 14 && dunemcSamples[sample].rw_isCCevent[event] == 1 &&
              dunemcSamples[sample].rw_cvnnumu_shifted[event] > numu_cut) {
                //std::cout<<"pdg of event in numu efficiency cut = " << dunemcSamples[sample].nu_pdg[event] << std::endl;
              eventsthatpasscut += *(MCSamples[sample].osc_w_pointer[event]); //+= 1;
          }
      }
  }
  std::cout << "total_true_ccnumu" << total_true_ccnumu << std::endl;
  std::cout << "total_events_incut" << total_events_incut << std::endl;
  std::cout << "eventsthatpasscut" << eventsthatpasscut << std::endl;

  if (total_true_ccnumu == 0) { // Prevent division by zero
      std::cerr << "Error: No CC νμ events found." << std::endl;
      return -1;
  }

  double efficiency = eventsthatpasscut / total_true_ccnumu;
  std::cout<< "efficiency = " << efficiency << std::endl;
  return efficiency;
}


double samplePDFDUNEBeamFD::getPurity_numu() {


  double total_events_incut_purity = 0;
  double eventsthatpasscut_purity = 0;

  int no_samples = (int)dunemcSamples.size();

  for (int sample = 0; sample < no_samples; ++sample) {
      for (int event = 0; event < dunemcSamples[sample].nEvents; ++event) {
          if (dunemcSamples[sample].rw_cvnnumu_shifted[event] > numu_cut) { 
              total_events_incut_purity += *(MCSamples[sample].osc_w_pointer[event]); //+= 1;
          }
          if (dunemcSamples[sample].nu_pdg[event] == 14 && dunemcSamples[sample].rw_isCCevent[event] == 1 &&
            dunemcSamples[sample].rw_cvnnumu_shifted[event] > numu_cut) {
             // std::cout<<"pdg of event in numu epurity cut = " << dunemcSamples[sample].nu_pdg[event] << std::endl;
            eventsthatpasscut_purity += *(MCSamples[sample].osc_w_pointer[event]);//+= 1;
        }
      }
  }

  if (total_events_incut_purity == 0) {  // Prevent division by zero
      std::cerr << "Error: No events in cut found." << std::endl;
      return -1;
  }
  std::cout << "eventsthatpasscut_purity" << eventsthatpasscut_purity << std::endl;
  std::cout << "total_events_incut_purity" << total_events_incut_purity << std::endl;

  double purity = eventsthatpasscut_purity / total_events_incut_purity;
  std::cout<< "purity = " << purity << std::endl;
  return purity;
}


TH1D *samplePDFDUNEBeamFD::get1DVarHist(KinematicTypes Var1, int kModeToFill,
                                        int kChannelToFill, int WeightStyle,
                                        TAxis *Axis) {
  bool fChannel;
  bool fMode;

  if (kChannelToFill != -1) {
    if (kChannelToFill > dunemcSamples.size()) {
      MACH3LOG_ERROR("Required channel is not available. kChannelToFill should "
                     "be between 0 and {}",
                     dunemcSamples.size());
      MACH3LOG_ERROR("kChannelToFill given: {}", kChannelToFill);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    fChannel = true;
  } else {
    fChannel = false;
  }

  if (kModeToFill != -1) {
    if (kModeToFill > kMaCh3_nModes) {
      MACH3LOG_ERROR("Required mode is not available. kModeToFill should be "
                     "between 0 and {}",
                     kMaCh3_nModes);
      MACH3LOG_ERROR("kModeToFill given: {}", kModeToFill);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    fMode = true;
  } else {
    fMode = false;
  }

  std::vector<std::vector<double>> SelectionVec;

  if (fMode) {
    std::vector<double> SelecMode(3);
    SelecMode[0] = kM3Mode;
    SelecMode[1] = kModeToFill;
    SelecMode[2] = kModeToFill + 1;
    SelectionVec.push_back(SelecMode);
  }

  if (fChannel) {
    std::vector<double> SelecChannel(3);
    SelecChannel[0] = kOscChannel;
    SelecChannel[1] = kChannelToFill;
    SelecChannel[2] = kChannelToFill + 1;
    SelectionVec.push_back(SelecChannel);
  }

  return get1DVarHist(Var1, SelectionVec, WeightStyle, Axis);
}

/*! DB New version of get1DVarHist which only fills histogram with events
 * passing IsEventSelected This works by having the Selection vector, where each
 * component of Selection is a 2 or 3 length vector If Selection[i].size()==3,
 * Selection[i][0] is the ND280KinematicType which is being cut, and only eventsnuhyp
 * with ND280KinematicType values between Selection[i][1] and Selection[i][2]
 * are accepted
 */
TH1D *
samplePDFDUNEBeamFD::get1DVarHist(KinematicTypes Var1,
                                  std::vector<std::vector<double>> SelectionVec,
                                  int WeightStyle, TAxis *Axis) {

  Selection = SelectionVec;

  for (unsigned int iStoredSelection = 0;
       iStoredSelection < StoredSelection.size(); iStoredSelection++) {
    Selection.push_back(StoredSelection[iStoredSelection]);
  }

  for (unsigned int iSelection = 0; iSelection < Selection.size();
       iSelection++) {
    if (Selection[iSelection].size() != 3) {
      MACH3LOG_ERROR("Selection Vector[{}] is not formed correctly. Expect "
                     "size == 3, given: {}",
                     iSelection, Selection[iSelection].size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }

  // DB Cut on OscChannel in this function due to speed increase from
  // considering duneSamples structure (ie. Array of length NChannels)
  bool fChannel = false;
  int kChannelToFill = -1;
  for (unsigned int iSelection = 0; iSelection < Selection.size();
       iSelection++) {
    if (Selection[iSelection][0] == kOscChannel) {
      fChannel = true;
      kChannelToFill = Selection[iSelection][1];
    }
  }

  if (fChannel && kChannelToFill > dunemcSamples.size()) {
    MACH3LOG_ERROR("Required channel is not available. kChannelToFill should "
                   "be between 0 and {}",
                   dunemcSamples.size());
    MACH3LOG_ERROR("kChannelToFill given: {}", kChannelToFill);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  TH1D *_h1DVar;
  std::vector<double> xBinEdges =
      ReturnKinematicParameterBinning(ReturnStringFromKinematicParameter(Var1));
  _h1DVar = new TH1D("", "", xBinEdges.size() - 1, xBinEdges.data());

  // This should be the same as FillArray in core basically, except that
  // events will end up in different bins
  for (int i = 0; i < dunemcSamples.size(); i++) {
    if (fChannel && (i != kChannelToFill)) {
      continue;
    }
    for (int ev_i = 0; ev_i < dunemcSamples[i].nEvents; ev_i++) {

      // DB Determine which events pass selection
      if (!IsEventSelected(i, ev_i)) {
        continue;
      }

      double Weight = GetEventWeight(i, ev_i);
      if (WeightStyle == 1) {
        Weight = *(MCSamples[i].osc_w_pointer[ev_i]) * dunemcSamples[i].pot_s *
                 dunemcSamples[i].norm_s * dunemcSamples[i].flux_w[ev_i];
      }

      // ETA - not sure about this
      if (MCSamples[i].xsec_w[ev_i] == 0.)
        continue;

      double Var1_Val;

      Var1_Val = ReturnKinematicParameter(Var1, i, ev_i);
      if (Var1_Val != _DEFAULT_RETURN_VAL_) {
        _h1DVar->Fill(Var1_Val, Weight);
      }
    }
  }

  /* DB: This is commented out be default
  // This code shifts the histogram meaning to Events/Bin Width but this affects
  the overall integral of the histogram so it should not be used anywhere we
  care about event rates
  // We could use Hist->Integral("width") but it would require a lot of
  modification throughout the code

  if (Var1!=kPDFBinning) {
    //_h1DVar->SetBinContent(1,_h1DVar->GetBinContent(0)+_h1DVar->GetBinContent(1));
    //_h1DVar->SetBinContent(_h1DVar->GetNbinsX(),_h1DVar->GetBinContent(_h1DVar->GetNbinsX())+_h1DVar->GetBinContent(_h1DVar->GetNbinsX()+1));

    for (int x=1;x<=_h1DVar->GetNbinsX();x++) {
      _h1DVar->SetBinContent(x,_h1DVar->GetBinContent(x)/_h1DVar->GetXaxis()->GetBinWidth(x));
    }

    _h1DVar->GetYaxis()->SetTitle("Events/Bin Width");
  }
  */

  return _h1DVar;
}

double samplePDFDUNEBeamFD::ReturnKinematicParameter(double KinematicVariable,
                                                     int iSample, int iEvent) {
  return *GetPointerToKinematicParameter(KinematicVariable, iSample, iEvent);
}

double
samplePDFDUNEBeamFD::ReturnKinematicParameter(std::string KinematicParameter,
                                              int iSample, int iEvent) {
  return ReturnKinematicParameter(
      ReturnKinematicParameterFromString(KinematicParameter), iSample, iEvent);
}

const double *samplePDFDUNEBeamFD::GetPointerToKinematicParameter(
    std::string KinematicParameter, int iSample, int iEvent) {
  return GetPointerToKinematicParameter(
      ReturnKinematicParameterFromString(KinematicParameter), iSample, iEvent);
}

int samplePDFDUNEBeamFD::ReturnKinematicParameterFromString(
    std::string KinematicParameterStr) {
  if (KinematicParameterStr.find("TrueNeutrinoEnergy") != std::string::npos) {
    return kTrueNeutrinoEnergy;
  }
  if (KinematicParameterStr.find("RecoNeutrinoEnergy") != std::string::npos) {
    return kRecoNeutrinoEnergy;
  }
  if (KinematicParameterStr.find("TrueXPos") != std::string::npos) {
    return kTrueXPos;
  }
  if (KinematicParameterStr.find("TrueYPos") != std::string::npos) {
    return kTrueYPos;
  }
  if (KinematicParameterStr.find("TrueZPos") != std::string::npos) {
    return kTrueZPos;
  }
  if (KinematicParameterStr.find("CVNNumu") != std::string::npos) {
    return kCVNNumu;
  }
  if (KinematicParameterStr.find("CVNNue") != std::string::npos) {
    return kCVNNue;
  }
  if (KinematicParameterStr.find("M3Mode") != std::string::npos) {
    return kM3Mode;
  }
  if (KinematicParameterStr.find("Numu_efficiency") != std::string::npos) {
    return kNumu_efficiency;
  }
  if (KinematicParameterStr.find("Nue_efficiency") != std::string::npos) {
    return kNue_efficiency;
  }
  if (KinematicParameterStr.find("selected_numuCCevent_energy") != std::string::npos) {
    return kselected_numuCCevent_energy;
  }
  if (KinematicParameterStr.find("selected_nueCCevent_energy") != std::string::npos) {
    return kselected_nueCCevent_energy;
  }
  if (KinematicParameterStr.find("selected_numuCCevent_vertexpos_x") != std::string::npos) {
    return kselected_numuCCevent_vertexpos_x;
  }
  if (KinematicParameterStr.find("selected_nueCCevent_vertexpos_x") != std::string::npos) {
    return kselected_nueCCevent_vertexpos_x;
  }
  if (KinematicParameterStr.find("selected_numuCCevent_vertexpos_y") != std::string::npos) {
    return kselected_numuCCevent_vertexpos_y;
  }
  if (KinematicParameterStr.find("selected_nueCCevent_vertexpos_y") != std::string::npos) {
    return kselected_nueCCevent_vertexpos_y;
  }
  if (KinematicParameterStr.find("selected_numuCCevent_vertexposz") != std::string::npos) {
    return kselected_numuCCevent_vertexpos_z;
  }
  if (KinematicParameterStr.find("selected_nueCCevent_vertexpos_z") != std::string::npos) {
    return kselected_nueCCevent_vertexpos_z;
  }
  if (KinematicParameterStr.find("ismc_numu") != std::string::npos) {
    return kismc_numu;
  }
  if (KinematicParameterStr.find("ismc_nue") != std::string::npos) {
    return kismc_nue;
  }
  if (KinematicParameterStr.find("ismc_nutau") != std::string::npos) {
    return kismc_nutau;
  }
  if (KinematicParameterStr.find("isCC") != std::string::npos) {
    return kisCC;
  }
  if (KinematicParameterStr.find("nuflavour") != std::string::npos) {
    return knuflavour;
  }
  if (KinematicParameterStr.find("reco_pdg") != std::string::npos) {
    return krecopdg;
  }
}

/*
kselected_numuCCevent_energy,
    kselected_numuCCevent_vertexpos_x,
    kselected_numuCCevent_vertexpos_y,
    kselected_numuCCevent_vertexpos_z,
    kselected_nueCCevent_energy,
    kselected_nueCCevent_vertexpos_x,
    kselected_nueCCevent_vertexpos_y,
    kselected_nueCCevent_vertexpos_z*/

const double *
samplePDFDUNEBeamFD::GetPointerToKinematicParameter(double KinematicVariable,
                                                    int iSample, int iEvent) {
  KinematicTypes KinPar = (KinematicTypes)std::round(KinematicVariable);

  switch (KinPar) {
  case kTrueNeutrinoEnergy:
    return &dunemcSamples[iSample].rw_etru[iEvent];
  case kRecoNeutrinoEnergy:
    return &dunemcSamples[iSample].rw_erec_shifted[iEvent];
  case kTrueXPos:
    return &dunemcSamples[iSample].rw_vtx_x[iEvent];
  case kTrueYPos:
    return &dunemcSamples[iSample].rw_vtx_y[iEvent];
  case kTrueZPos:
    return &dunemcSamples[iSample].rw_vtx_z[iEvent];
  case kCVNNumu:
    return &dunemcSamples[iSample].rw_cvnnumu_shifted[iEvent];
  case kCVNNue:
    return &dunemcSamples[iSample].rw_cvnnue_shifted[iEvent];
  case kNumu_efficiency:
    return &dunemcSamples[iSample].numuefficency[iEvent];
  case kNue_efficiency:
    return &dunemcSamples[iSample].nueefficency[iEvent];
  case kselected_numuCCevent_energy:
      //std::cout << "in selected_numuCCevent_energy case, the energy is = " << std::endl;
    return &dunemcSamples[iSample].selected_numuCCevent_energy[iEvent];
  case kselected_nueCCevent_energy:
    return &dunemcSamples[iSample].selected_nueCCevent_energy[iEvent];
  case kselected_numuCCevent_vertexpos_x:
    return &dunemcSamples[iSample].selected_numuCCevent_vertexpos_x[iEvent];
  case kselected_nueCCevent_vertexpos_x:
    return &dunemcSamples[iSample].selected_nueCCevent_vertexpos_x[iEvent];
  case kselected_numuCCevent_vertexpos_y:
    return &dunemcSamples[iSample].selected_numuCCevent_vertexpos_y[iEvent];
  case kselected_nueCCevent_vertexpos_y:
    return &dunemcSamples[iSample].selected_nueCCevent_vertexpos_y[iEvent];
  case kselected_numuCCevent_vertexpos_z:
    return &dunemcSamples[iSample].selected_numuCCevent_vertexpos_z[iEvent];
  case kselected_nueCCevent_vertexpos_z:
    return &dunemcSamples[iSample].selected_nueCCevent_vertexpos_z[iEvent];
  //default:
  case kismc_numu:
    return &dunemcSamples[iSample].mc_isnumu[iEvent];
  case kismc_nue:
    return &dunemcSamples[iSample].mc_isnue[iEvent];
  case kismc_nutau:
    return &dunemcSamples[iSample].mc_isnutau[iEvent];
  case kisCC:
    return &dunemcSamples[iSample].rw_isCCevent[iEvent];
  case knuflavour:
    return &dunemcSamples[iSample].nuflavour[iEvent];
  case krecopdg:
    return &dunemcSamples[iSample].nu_pdg[iEvent];
  default:
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type...");
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}

inline std::string samplePDFDUNEBeamFD::ReturnStringFromKinematicParameter(
    int KinematicParameter) {
  std::string KinematicString = "";

  switch (KinematicParameter) {
  case kRecoNeutrinoEnergy:
    KinematicString = "RecoNeutrinoEnergy";
    break;
  case kTrueNeutrinoEnergy:
    KinematicString = "TrueNeutrinoEnergy";
    break;
  case kTrueXPos:
    KinematicString = "TrueXPos";
    break;
  case kTrueYPos:
    KinematicString = "TrueYPos";
    break;
  case kTrueZPos:
    KinematicString = "TrueZPos";
    break;
  case kCVNNumu:
    KinematicString = "CVNNumu";
    break;
  case kCVNNue:
    KinematicString = "CVNNue";
    break;
  case kM3Mode:
    KinematicString = "M3Mode";
    break;
  case kselected_numuCCevent_energy:
    KinematicString = "selected_numuCCevent_energy";
    break;
  case kselected_nueCCevent_energy:
    KinematicString = "selected_nueCCevent_energy";
    break;
  case kselected_numuCCevent_vertexpos_x:
    KinematicString ="selected_numuCCevent_vertexpos_x";
    break;
  case kselected_nueCCevent_vertexpos_x:
    KinematicString ="selected_nueCCevent_vertexpos_x";
    break;
  case kselected_numuCCevent_vertexpos_y:
    KinematicString ="selected_numuCCevent_vertexpos_y";
    break;
  case kselected_nueCCevent_vertexpos_y:
    KinematicString ="selected_nueCCevent_vertexpos_y";
    break;
  case kselected_numuCCevent_vertexpos_z:
    KinematicString ="selected_numuCCevent_vertexpos_z";
    break;
  case kselected_nueCCevent_vertexpos_z:
    KinematicString ="selected_nueCCevent_vertexpos_z";
    break;
  case kismc_numu:
    KinematicString ="ismc_numu";
    break;
  case kismc_nue:
    KinematicString ="ismc_nue";
    break;
  case kismc_nutau:
    KinematicString ="ismc_nutau";
    break;
  case kisCC:
    KinematicString ="isCC";
    break;
  case knuflavour:
    KinematicString ="nuflavour";
    break;
  case krecopdg:
    KinematicString ="reco_pdg";
    break;
  default:
    break;
  }

  return KinematicString;
}

void samplePDFDUNEBeamFD::setupFDMC(int iSample) {
  dunemc_base *duneobj = &(dunemcSamples[iSample]);
  fdmc_base *fdobj = &(MCSamples[iSample]);

  fdobj->nutype = duneobj->nutype;
  fdobj->oscnutype = duneobj->oscnutype;
  fdobj->signal = duneobj->signal;
  fdobj->SampleDetID = SampleDetID;

  for (int iEvent = 0; iEvent < fdobj->nEvents; ++iEvent) {
    fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
    fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
    fdobj->Target[iEvent] = &(duneobj->Target[iEvent]);
    fdobj->isNC[iEvent] = !(duneobj->rw_isCCevent[iEvent]);
  }
}

void samplePDFDUNEBeamFD::applyShifts(int iSample, int iEvent) {}

std::vector<double> samplePDFDUNEBeamFD::ReturnKinematicParameterBinning(
    std::string KinematicParameterStr) {
  std::vector<double> binningVector;
  KinematicTypes KinematicParameter = static_cast<KinematicTypes>(
      ReturnKinematicParameterFromString(KinematicParameterStr));

  int nBins = 0;
  double bin_width = 0;
  switch (KinematicParameter) {
  case (kRecoNeutrinoEnergy):
    nBins = 20;
    bin_width = 0.5; // GeV
    break;
  case (kTrueNeutrinoEnergy):
    nBins = 20;
    bin_width = 0.5; // GeV
    break;
  default:
    nBins = 10;
    bin_width = 1.0;
    break;
  }

  for (int bin_i = 0; bin_i < nBins; bin_i++) {
    binningVector.push_back(bin_i * bin_width);
  }

  return binningVector;
}
