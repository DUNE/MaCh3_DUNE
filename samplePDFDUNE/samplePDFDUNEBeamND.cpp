#include <TROOT.h>

#include "samplePDFDUNEBeamND.h"
#include "StructsDUNE.h"
#include "TString.h"
#include <assert.h>
#include <manager/MaCh3Logger.h>
#include <stdexcept>
#include <string>
#include "TMath.h"
#include "manager/manager.h"

samplePDFDUNEBeamND::samplePDFDUNEBeamND(std::string mc_version_, covarianceXsec* xsec_cov_, covarianceOsc* osc_cov_) : samplePDFFDBase(mc_version_, xsec_cov_, osc_cov_) {  
  Initialise();
}

samplePDFDUNEBeamND::~samplePDFDUNEBeamND() {
}

void samplePDFDUNEBeamND::Init() {
  dunendmcSamples.resize(nSamples,dunemc_base());
  funcParsGrid.resize(nSamples);

  if (CheckNodeExists(SampleManager->raw(), "POT")) {
    pot = SampleManager->raw()["POT"].as<double>();
  } else{
    MACH3LOG_ERROR("POT not defined in {}, please add this!", SampleManager->GetFileName());
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  IsRHC = SampleManager->raw()["SampleBools"]["isrhc"].as<bool>();
  SampleDetID = SampleManager->raw()["DetID"].as<int>();
  iselike = SampleManager->raw()["SampleBools"]["iselike"].as<bool>();

  // SetupFunctionalParameters();

  /*
  tot_escale_nd_pos = -999;
  tot_escale_sqrt_nd_pos = -999;
  tot_escale_invsqrt_nd_pos = -999;
  had_escale_nd_pos = -999;
  had_escale_sqrt_nd_pos = -999;
  had_escale_invsqrt_nd_pos = -999;
  mu_escale_nd_pos = -999;
  mu_escale_sqrt_nd_pos = -999;
  mu_escale_invsqrt_nd_pos = -999;
  n_escale_nd_pos = -999;
  n_escale_sqrt_nd_pos = -999;
  n_escale_invsqrt_nd_pos = -999;
  em_escale_nd_pos = -999;
  em_escale_sqrt_nd_pos = -999;
  em_escale_invsqrt_nd_pos = -999;
  had_res_nd_pos = -999;
  mu_res_nd_pos = -999;
  n_res_nd_pos = -999;
  em_res_nd_pos = -999;

  std::vector<std::string> funcParsNames = XsecCov->GetParsNamesFromDetID(SampleDetID, SystType::kFunc);
  std::vector<int> funcParsIndex = XsecCov->GetParsIndexFromDetID(SampleDetID, SystType::kFunc);

  nNDDetectorSystPointers = funcParsIndex.size();
  NDDetectorSystPointers = std::vector<const double*>(nNDDetectorSystPointers);

  int func_it = 0;
  for (std::vector<int>::iterator it = funcParsIndex.begin(); it != funcParsIndex.end(); ++it, ++func_it) {
    std::string name = funcParsNames.at(func_it);
    
    if (name == "TotalEScaleND") {
      tot_escale_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(tot_escale_nd_pos);
    }
    else if (name == "TotalEScaleSqrtND") {
      tot_escale_sqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(tot_escale_sqrt_nd_pos);
    }
    else if (name == "TotalEScaleInvSqrtND") {
      tot_escale_invsqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(tot_escale_invsqrt_nd_pos);
    }
    else if (name == "HadEScaleND") {
      had_escale_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(had_escale_nd_pos);
    }
    else if (name == "HadEScaleSqrtND") {
      had_escale_sqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(had_escale_sqrt_nd_pos);
    }
    else if (name == "HadEScaleInvSqrtND") {
      had_escale_invsqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(had_escale_invsqrt_nd_pos);
    }
    else if (name == "MuEScaleND") {
      mu_escale_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_escale_nd_pos);
    }
    else if (name == "MuEScaleSqrtND") {
      mu_escale_sqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_escale_sqrt_nd_pos);
    }
    else if (name == "MuEScaleInvSqrtND") {
      mu_escale_invsqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_escale_invsqrt_nd_pos);
    }
    else if (name == "NEScaleND") {
      n_escale_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(n_escale_nd_pos);
    }
    else if (name == "NEScaleSqrtND") {
      n_escale_sqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(n_escale_sqrt_nd_pos);
    }
    else if (name == "NEScaleInvSqrtND") {
      n_escale_invsqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(n_escale_invsqrt_nd_pos);
    }
    else if (name == "EMEScaleND") {
      em_escale_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(em_escale_nd_pos);
    }
    else if (name == "EMEScaleSqrtND") {
      em_escale_sqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(em_escale_sqrt_nd_pos);
    }
    else if (name == "EMEScaleInvSqrtND") {
      em_escale_invsqrt_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(em_escale_invsqrt_nd_pos);
    }
    else if (name == "HadResND") {
      had_res_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(had_res_nd_pos);
    }
    else if (name == "MuResND") {
      mu_res_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(mu_res_nd_pos);
    }
    else if (name == "NResND") {
      n_res_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(n_res_nd_pos);
    }
    else if (name == "EMResND") {
      em_res_nd_pos = *it;
      NDDetectorSystPointers[func_it] = XsecCov->retPointer(em_res_nd_pos);
    }
    else { 
      MACH3LOG_ERROR("Found a functional parameter which wasn't specified in the xml | samplePDFDUNEBeamND: {}",name);
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
  */
  
  std::cout << "-------------------------------------------------------------------" <<std::endl;
}

void samplePDFDUNEBeamND::SetupSplines() {
  ///@todo move all of the spline setup into core
  if(XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline) > 0){
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    SplineHandler = std::unique_ptr<splineFDBase>(new splinesDUNE(XsecCov));
    InitialiseSplineObject();
  }
  else{
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or evaluate splines", XsecCov->GetNumParamsFromDetID(SampleDetID, kSpline));
    SplineHandler = nullptr;
  }

  return;
}

// === HH: Functional parameters ===
void samplePDFDUNEBeamND::TotalEScale(const double * par, std::size_t iSample, std::size_t iEvent) {
  // Total energy scale uncertainties for anything but CC Numu, see: 
  // https://github.com/DUNE/lblpwgtools/blob/3d475f50a998fbfa6266df9a0c4eb3056c0cdfe5/CAFAna/Systs/EnergySysts.h#L39
  
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_erec_had[iEvent];
}

void samplePDFDUNEBeamND::TotalEScaleNotCCNumu(const double * par, std::size_t iSample, std::size_t iEvent) {
  // A special case for Not (CC Numu), where we also scale Erec by lepton energy
  // Since we reconstruct muon energy in a different way, see:
  // https://github.com/DUNE/lblpwgtools/blob/3d475f50a998fbfa6266df9a0c4eb3056c0cdfe5/CAFAna/Systs/EnergySysts.h#L39
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_erec_lep[iEvent];
}

void samplePDFDUNEBeamND::TotalEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_erec_had[iEvent] * dunendmcSamples[iSample].rw_erec_had_sqrt[iEvent];
}

void samplePDFDUNEBeamND::TotalEScaleSqrtNotCCNumu(const double * par, std::size_t iSample, std::size_t iEvent) {
  // See comments in TotalEScaleNotCCNumu
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_erec_lep[iEvent] * dunendmcSamples[iSample].rw_erec_lep_sqrt[iEvent];
}

void samplePDFDUNEBeamND::TotalEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  // Erec/sqrt(Erec) = sqrt(Erec)
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_erec_had_sqrt[iEvent];
}

void samplePDFDUNEBeamND::TotalEScaleInvSqrtNotCCNumu(const double * par, std::size_t iSample, std::size_t iEvent) {
  // See comments in TotalEScaleNotCCNumu
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_erec_lep_sqrt[iEvent];
}

void samplePDFDUNEBeamND::HadEScale(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_sum_ehad[iEvent];
}

void samplePDFDUNEBeamND::HadEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_sum_ehad[iEvent] * dunendmcSamples[iSample].rw_sum_ehad_sqrt[iEvent];
}

void samplePDFDUNEBeamND::HadEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_sum_ehad_sqrt[iEvent];
}

void samplePDFDUNEBeamND::MuEScale(const double * par, std::size_t iSample, std::size_t iEvent) {
  // HH TODO: Functionally this is the same as TotalEScaleNotCCNumu, not sure if this function is even needed
  TotalEScaleNotCCNumu(par, iSample, iEvent);
}

void samplePDFDUNEBeamND::MuEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  // See comments in MuEScale
  TotalEScaleSqrtNotCCNumu(par, iSample, iEvent);
}

void samplePDFDUNEBeamND::MuEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  // See comments in MuEScale
  TotalEScaleInvSqrtNotCCNumu(par, iSample, iEvent);
}

void samplePDFDUNEBeamND::NEScale(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_eRecoN[iEvent];
}

void samplePDFDUNEBeamND::NEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_eRecoN[iEvent] * dunendmcSamples[iSample].rw_eRecoN_sqrt[iEvent];
}

void samplePDFDUNEBeamND::NEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_eRecoN_sqrt[iEvent];
}

void samplePDFDUNEBeamND::EMEScale(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_eRecoPi0[iEvent];
}

void samplePDFDUNEBeamND::EMEScaleCCNue(const double * par, std::size_t iSample, std::size_t iEvent) {
  // Again this is the same as TotalEScaleNotCCNumu, not sure if this function is needed
  TotalEScaleNotCCNumu(par, iSample, iEvent);
}

void samplePDFDUNEBeamND::EMEScaleSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_eRecoPi0[iEvent] * dunendmcSamples[iSample].rw_eRecoPi0_sqrt[iEvent];
}

void samplePDFDUNEBeamND::EMEScaleSqrtCCNue(const double * par, std::size_t iSample, std::size_t iEvent) {
  // See comments in EMEScaleCCNue
  TotalEScaleSqrtNotCCNumu(par, iSample, iEvent);
}

void samplePDFDUNEBeamND::EMEScaleInvSqrt(const double * par, std::size_t iSample, std::size_t iEvent) {
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * dunendmcSamples[iSample].rw_eRecoPi0_sqrt[iEvent];
}

void samplePDFDUNEBeamND::EMEScaleInvSqrtCCNue(const double * par, std::size_t iSample, std::size_t iEvent) {
  // See comments in EMEScaleCCNue
  TotalEScaleInvSqrtNotCCNumu(par, iSample, iEvent);
}

void samplePDFDUNEBeamND::HadRes(const double * par, std::size_t iSample, std::size_t iEvent) {
  // True sum - reco sum
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * (dunendmcSamples[iSample].rw_eP[iEvent] 
    + dunendmcSamples[iSample].rw_ePip[iEvent] 
    + dunendmcSamples[iSample].rw_ePim[iEvent] 
    - dunendmcSamples[iSample].rw_sum_ehad[iEvent]);
}

void samplePDFDUNEBeamND::MuRes(const double * par, std::size_t iSample, std::size_t iEvent) {
  // True muon energy - reco muon energy
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * (dunendmcSamples[iSample].rw_LepE[iEvent] - dunendmcSamples[iSample].rw_erec_lep[iEvent]);
}

void samplePDFDUNEBeamND::NRes(const double * par, std::size_t iSample, std::size_t iEvent) {
  // True neutron energy - reco neutron energy
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * (dunendmcSamples[iSample].rw_eN[iEvent] - dunendmcSamples[iSample].rw_eRecoN[iEvent]);
}

void samplePDFDUNEBeamND::EMRes(const double * par, std::size_t iSample, std::size_t iEvent) {
  // True pi0 energy - reco pi0 energy
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] += (*par) * (dunendmcSamples[iSample].rw_ePi0[iEvent] - dunendmcSamples[iSample].rw_eRecoPi0[iEvent]);
}

void samplePDFDUNEBeamND::EMResCCNue(const double * par, std::size_t iSample, std::size_t iEvent) {
  // This is the same as MuRes, again not sure if this function is needed
  MuRes(par, iSample, iEvent);
}

void samplePDFDUNEBeamND::DebugShift(const double * par, std::size_t iSample, std::size_t iEvent) {
  if (dunendmcSamples[iSample].rw_erec[iEvent] < 2) {
    dunendmcSamples[iSample].rw_erec_shifted[iEvent] = 4.0;
  }
}

void samplePDFDUNEBeamND::RegisterFunctionalParameters() {
  MACH3LOG_INFO("Registering functional parameters");
  // This function manually populates the map of functional parameters
  // Maps the name of the functional parameter to the pointer of the function
  
  // This is the part where we manually enter things
  // A lambda function has to be used so we can refer to a non-static member function
  RegisterIndividualFuncPar("DebugNothing", 
                            kDebugNothing, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) {});

  RegisterIndividualFuncPar("DebugShift",
                            kDebugShift, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->DebugShift(par, iSample, iEvent); });

  RegisterIndividualFuncPar("TotalEScaleND",
                            kTotalEScale, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->TotalEScale(par, iSample, iEvent); });

  RegisterIndividualFuncPar("TotalEScaleNotCCNumuND",
                            kTotalEScaleNotCCNumu, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->TotalEScaleNotCCNumu(par, iSample, iEvent); });

  RegisterIndividualFuncPar("TotalEScaleSqrtND",
                            kTotalEScaleSqrt, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->TotalEScaleSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("TotalEScaleSqrtNotCCNumuND",
                            kTotalEScaleSqrtNotCCNumu, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->TotalEScaleSqrtNotCCNumu(par, iSample, iEvent); });

  RegisterIndividualFuncPar("TotalEScaleInvSqrtND",
                            kTotalEScaleInvSqrt, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->TotalEScaleInvSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("TotalEScaleInvSqrtNotCCNumuND",
                            kTotalEScaleInvSqrtNotCCNumu, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->TotalEScaleInvSqrtNotCCNumu(par, iSample, iEvent); });

  RegisterIndividualFuncPar("HadEScaleND",
                            kHadEScale, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->HadEScale(par, iSample, iEvent); });

  RegisterIndividualFuncPar("HadEScaleSqrtND",
                            kHadEScaleSqrt, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->HadEScaleSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("HadEScaleInvSqrtND",
                            kHadEScaleInvSqrt, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->HadEScaleInvSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("MuEScaleND",
                            kMuEScale, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->MuEScale(par, iSample, iEvent); });

  RegisterIndividualFuncPar("MuEScaleSqrtND",
                            kMuEScaleSqrt, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->MuEScaleSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("MuEScaleInvSqrtND",
                            kMuEScaleInvSqrt, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->MuEScaleInvSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("NEScaleND",
                            kNEScale, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->NEScale(par, iSample, iEvent); });

  RegisterIndividualFuncPar("NEScaleSqrtND",
                            kNEScaleSqrt, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->NEScaleSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("NEScaleInvSqrtND",
                            kNEScaleInvSqrt, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->NEScaleInvSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMEScaleND",
                            kEMEScale, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMEScale(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMEScaleCCNueND",
                            kEMEScaleCCNue, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMEScaleCCNue(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMEScaleSqrtND",
                            kEMEScaleSqrt, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMEScaleSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMEScaleSqrtCCNueND",
                            kEMEScaleSqrtCCNue, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMEScaleSqrtCCNue(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMEScaleInvSqrtND",
                            kEMEScaleInvSqrt, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMEScaleInvSqrt(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMEScaleInvSqrtCCNueND",
                            kEMEScaleInvSqrtCCNue, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMEScaleInvSqrtCCNue(par, iSample, iEvent); });

  RegisterIndividualFuncPar("HadResND",
                            kHadRes, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->HadRes(par, iSample, iEvent); });

  RegisterIndividualFuncPar("MuResND",
                            kMuRes, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->MuRes(par, iSample, iEvent); });

  RegisterIndividualFuncPar("NResND",
                            kNRes, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->NRes(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMResND",
                            kEMRes, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMRes(par, iSample, iEvent); });

  RegisterIndividualFuncPar("EMResCCNueND",
                            kEMResCCNue, 
                            [this](const double * par, std::size_t iSample, std::size_t iEvent) { this->EMResCCNue(par, iSample, iEvent); });

  MACH3LOG_INFO("Finished registering functional parameters");
}

// HH: Reset the shifted values to the original values
void samplePDFDUNEBeamND::resetShifts(int iSample, int iEvent) {
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] = dunendmcSamples[iSample].rw_erec[iEvent];
}
// =================================

void samplePDFDUNEBeamND::SetupWeightPointers() {
  for (int i = 0; i < (int)dunendmcSamples.size(); ++i) {
    for (int j = 0; j < dunendmcSamples[i].nEvents; ++j) {
      MCSamples[i].ntotal_weight_pointers[j] = 6;
      MCSamples[i].total_weight_pointers[j].resize(MCSamples[i].ntotal_weight_pointers[j]);
      MCSamples[i].total_weight_pointers[j][0] = &(dunendmcSamples[i].pot_s);
      MCSamples[i].total_weight_pointers[j][1] = &(dunendmcSamples[i].norm_s);
      MCSamples[i].total_weight_pointers[j][2] = MCSamples[i].osc_w_pointer[j];
      MCSamples[i].total_weight_pointers[j][3] = &(dunendmcSamples[i].rw_berpaacvwgt[j]);
      MCSamples[i].total_weight_pointers[j][4] = &(dunendmcSamples[i].flux_w[j]);
      MCSamples[i].total_weight_pointers[j][5] = &(MCSamples[i].xsec_w[j]);
    }
  }
}

int samplePDFDUNEBeamND::setupExperimentMC(int iSample) {

  dunemc_base *duneobj = &(dunendmcSamples[iSample]);
  int nupdgUnosc = sample_nupdgunosc[iSample];
  int nupdg = sample_nupdg[iSample];
  
  MACH3LOG_INFO("-------------------------------------------------------------------");
  MACH3LOG_INFO("input file: {}", mc_files.at(iSample));
  
  _sampleFile = TFile::Open(mc_files.at(iSample).c_str(), "READ");
  _data = (TTree*)_sampleFile->Get("caf");

  if(_data){
    MACH3LOG_INFO("Found \"caf\" tree in {}", mc_files[iSample]);
    MACH3LOG_INFO("With number of entries: {}", _data->GetEntries());
  }
  else{
    MACH3LOG_ERROR("Could not find \"caf\" tree in {}", mc_files[iSample]);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  _data->SetBranchStatus("*", 0);
  _data->SetBranchStatus("Ev", 1);
  _data->SetBranchAddress("Ev", &_ev);
  _data->SetBranchStatus("Ev_reco", 1);
  _data->SetBranchAddress("Ev_reco", &_erec);
  _data->SetBranchStatus("Elep_reco", 1);
  _data->SetBranchAddress("Elep_reco", &_erec_lep);
  _data->SetBranchStatus("mode",1);
  _data->SetBranchAddress("mode",&_mode);
  _data->SetBranchStatus("isCC", 1);
  _data->SetBranchAddress("isCC", &_isCC);
  _data->SetBranchStatus("isFHC", 1);
  _data->SetBranchAddress("isFHC", &_isFHC);
  _data->SetBranchStatus("reco_q", 1);
  _data->SetBranchAddress("reco_q", &_reco_q);
  _data->SetBranchStatus("nuPDGunosc", 1);
  _data->SetBranchAddress("nuPDGunosc", &_nuPDGunosc);
  _data->SetBranchStatus("nuPDG", 1);
  _data->SetBranchAddress("nuPDG", &_nuPDG);
  _data->SetBranchStatus("BeRPA_A_cvwgt", 1);
  _data->SetBranchAddress("BeRPA_A_cvwgt", &_BeRPA_cvwgt);

  _data->SetBranchStatus("eRecoP", 1);
  _data->SetBranchAddress("eRecoP", &_eRecoP);
  _data->SetBranchStatus("eRecoPip", 1);
  _data->SetBranchAddress("eRecoPip", &_eRecoPip);
  _data->SetBranchStatus("eRecoPim", 1);
  _data->SetBranchAddress("eRecoPim", &_eRecoPim);
  _data->SetBranchStatus("eRecoPi0", 1);
  _data->SetBranchAddress("eRecoPi0", &_eRecoPi0);
  _data->SetBranchStatus("eRecoN", 1);
  _data->SetBranchAddress("eRecoN", &_eRecoN);

  _data->SetBranchStatus("LepE", 1);
  _data->SetBranchAddress("LepE", &_LepE);
  _data->SetBranchStatus("eP", 1);
  _data->SetBranchAddress("eP", &_eP);
  _data->SetBranchStatus("ePip", 1);
  _data->SetBranchAddress("ePip", &_ePip);
  _data->SetBranchStatus("ePim", 1);
  _data->SetBranchAddress("ePim", &_ePim);
  _data->SetBranchStatus("ePi0", 1);
  _data->SetBranchAddress("ePi0", &_ePi0);
  _data->SetBranchStatus("eN", 1);
  _data->SetBranchAddress("eN", &_eN);

  _data->SetBranchStatus("reco_nue", 1);
  _data->SetBranchAddress("reco_nue", &_reco_nue);
  _data->SetBranchStatus("reco_numu", 1);
  _data->SetBranchAddress("reco_numu", &_reco_numu);

  if (!IsRHC) { 
    duneobj->norm_s = (1e21/1.5e21);
  } else {
    duneobj->norm_s = (1e21/1.905e21);
  }
  duneobj->pot_s = (pot)/1e21;
  
  std::cout << "pot: " << pot << std::endl;
  std::cout << "pot_s: " << duneobj->pot_s << std::endl;
  std::cout << "norm_s: " << duneobj->norm_s << std::endl;

  duneobj->nEvents = _data->GetEntries();

  // HH: Downsampling by choosing the first X% of the events
  // HH TODO: Instead of choosing the first X% of the events, we should randomly sample X% of the events 
  if (CheckNodeExists(SampleManager->raw(), "Downsample")) {
    double downsample = SampleManager->raw()["Downsample"].as<double>();
    MACH3LOG_INFO("Downsample found in {}, will sample only {} of each file and scale POT correspondingly!", SampleManager->GetFileName(), downsample);
    duneobj->nEvents = (int)(duneobj->nEvents*downsample);
    duneobj->pot_s = duneobj->pot_s/downsample;
    MACH3LOG_INFO("New number of events: {}", duneobj->nEvents);
    MACH3LOG_INFO("New POT: {}", pot);
    MACH3LOG_INFO("New pot_s: {}", duneobj->pot_s);
    MACH3LOG_INFO("New norm_s: {}", duneobj->norm_s);
  } else{
    MACH3LOG_INFO("Downsample not defined in {}, continuing without downsampling!", SampleManager->GetFileName());
  }



  duneobj->rw_yrec = new double[duneobj->nEvents];
  duneobj->rw_erec_lep = new double[duneobj->nEvents];
  duneobj->rw_erec_had = new double[duneobj->nEvents];
  duneobj->rw_etru = new double[duneobj->nEvents];
  duneobj->rw_erec = new double[duneobj->nEvents];
  duneobj->rw_erec_shifted = new double[duneobj->nEvents];
  duneobj->rw_theta = new double[duneobj->nEvents];
  duneobj->flux_w = new double[duneobj->nEvents];
  duneobj->rw_isCC = new int[duneobj->nEvents];
  duneobj->rw_isFHC = new double[duneobj->nEvents];
  duneobj->rw_reco_q = new double[duneobj->nEvents];
  duneobj->rw_nuPDGunosc = new int[duneobj->nEvents];
  duneobj->rw_nuPDG = new int[duneobj->nEvents];
  duneobj->rw_berpaacvwgt = new double[duneobj->nEvents]; 

  duneobj->rw_eRecoP = new double[duneobj->nEvents];
  duneobj->rw_eRecoPip = new double[duneobj->nEvents];
  duneobj->rw_eRecoPim = new double[duneobj->nEvents];
  duneobj->rw_eRecoPi0 = new double[duneobj->nEvents];
  duneobj->rw_eRecoN = new double[duneobj->nEvents];

  duneobj->rw_LepE = new double[duneobj->nEvents];
  duneobj->rw_eP = new double[duneobj->nEvents];
  duneobj->rw_ePip = new double[duneobj->nEvents];
  duneobj->rw_ePim = new double[duneobj->nEvents];
  duneobj->rw_ePi0 = new double[duneobj->nEvents];
  duneobj->rw_eN = new double[duneobj->nEvents];

  duneobj->nupdg = new int[duneobj->nEvents];
  duneobj->nupdgUnosc = new int[duneobj->nEvents];
  duneobj->mode = new double[duneobj->nEvents];
  duneobj->Target = new int[duneobj->nEvents];

  duneobj->rw_reco_numu = new int[duneobj->nEvents];
  duneobj->rw_reco_nue = new int[duneobj->nEvents];

  duneobj->rw_erec_had_sqrt = new double[duneobj->nEvents];
  duneobj->rw_erec_lep_sqrt = new double[duneobj->nEvents];
  duneobj->rw_eRecoN_sqrt = new double[duneobj->nEvents];
  duneobj->rw_eRecoPi0_sqrt = new double[duneobj->nEvents];

  duneobj->rw_sum_ehad = new double[duneobj->nEvents];
  duneobj->rw_sum_ehad_sqrt = new double[duneobj->nEvents];

  _data->GetEntry(0);

  //FILL DUNE STRUCT
  for (int i = 0; i < duneobj->nEvents; ++i) { // Loop through tree
    _data->GetEntry(i);

    duneobj->nupdg[i] = sample_nupdg[iSample];
    duneobj->nupdgUnosc[i] = sample_nupdgunosc[iSample];
    
    duneobj->rw_erec[i] = (double)_erec;
    duneobj->rw_erec_shifted[i] = (double)_erec;
    duneobj->rw_erec_lep[i] = (double)_erec_lep;
    duneobj->rw_erec_had[i] = (double)(_erec - _erec_lep);
    duneobj->rw_yrec[i] = (double)((_erec - _erec_lep)/_erec);
    duneobj->rw_etru[i] = (double)_ev; // in GeV
    duneobj->rw_theta[i] = (double)_LepNuAngle;
    duneobj->rw_isCC[i] = _isCC;
    duneobj->rw_isFHC[i] = (double)_isFHC;
    duneobj->rw_reco_q[i] = _reco_q;
    duneobj->rw_nuPDGunosc[i] = _nuPDGunosc;
    duneobj->rw_nuPDG[i] = _nuPDG;
    duneobj->rw_berpaacvwgt[i] = (double)_BeRPA_cvwgt;
    
    duneobj->rw_eRecoP[i] = (double)_eRecoP; 
    duneobj->rw_eRecoPip[i] = (double)_eRecoPip; 
    duneobj->rw_eRecoPim[i] = (double)_eRecoPim; 
    duneobj->rw_eRecoPi0[i] = (double)_eRecoPi0; 
    duneobj->rw_eRecoN[i] = (double)_eRecoN; 
    
    duneobj->rw_LepE[i] = (double)_LepE; 
    duneobj->rw_eP[i] = (double)_eP; 
    duneobj->rw_ePip[i] = (double)_ePip; 
    duneobj->rw_ePim[i] = (double)_ePim; 
    duneobj->rw_ePi0[i] = (double)_ePi0; 
    duneobj->rw_eN[i] = (double)_eN; 

    duneobj->rw_reco_numu[i] = static_cast<int>(_reco_numu);
    duneobj->rw_reco_nue[i] = static_cast<int>(_reco_nue);

    duneobj->rw_erec_had_sqrt[i] = sqrt(duneobj->rw_erec_had[i]);
    duneobj->rw_erec_lep_sqrt[i] = sqrt(duneobj->rw_erec_lep[i]);
    duneobj->rw_eRecoN_sqrt[i] = sqrt(duneobj->rw_eRecoN[i]);
    duneobj->rw_eRecoPi0_sqrt[i] = sqrt(duneobj->rw_eRecoPi0[i]);

    duneobj->rw_sum_ehad[i] = duneobj->rw_eRecoP[i] + duneobj->rw_eRecoPip[i] + duneobj->rw_eRecoPim[i];
    duneobj->rw_sum_ehad_sqrt[i] = sqrt(duneobj->rw_sum_ehad[i]);

    //Assume everything is on Argon for now....
    duneobj->Target[i] = 40;
    
    int mode= TMath::Abs(_mode);       
    duneobj->mode[i]=(double)GENIEMode_ToMaCh3Mode(mode, _isCC);
    
    duneobj->flux_w[i] = 1.0;
  }
  
  _sampleFile->Close();
  return duneobj->nEvents;
}


const double* samplePDFDUNEBeamND::GetPointerToKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent) {
  double* KinematicValue;
  
  switch(KinPar){
  case kTrueNeutrinoEnergy:
    KinematicValue = &dunendmcSamples[iSample].rw_etru[iEvent]; 
    break;
  case kRecoQ:
    KinematicValue = &dunendmcSamples[iSample].rw_reco_q[iEvent];
    break;
  case kRecoNeutrinoEnergy:
    // HH: Changed this to match BeamFD
    KinematicValue = &dunendmcSamples[iSample].rw_erec_shifted[iEvent];
    break;
  case kIsFHC:
    KinematicValue = &dunendmcSamples[iSample].rw_isFHC[iEvent];
    break;
  case kRecoY:
    KinematicValue = &dunendmcSamples[iSample].rw_yrec[iEvent];
    break;
  default:
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type...");
    std::cout << KinPar << ReturnStringFromKinematicParameter(KinPar) << std::endl;
    throw MaCh3Exception(__FILE__, __LINE__);
  }
  
  return KinematicValue;
}

const double* samplePDFDUNEBeamND::GetPointerToKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = (KinematicTypes) std::round(KinematicVariable);
  return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

const double* samplePDFDUNEBeamND::GetPointerToKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return GetPointerToKinematicParameter(KinPar,iSample,iEvent);
}

double samplePDFDUNEBeamND::ReturnKinematicParameter(double KinematicVariable, int iSample, int iEvent) {
  KinematicTypes KinPar = (KinematicTypes) std::round(KinematicVariable);
  return ReturnKinematicParameter(KinPar,iSample,iEvent);
}

double samplePDFDUNEBeamND::ReturnKinematicParameter(std::string KinematicParameter, int iSample, int iEvent) {
  KinematicTypes KinPar = static_cast<KinematicTypes>(ReturnKinematicParameterFromString(KinematicParameter));
  return ReturnKinematicParameter(KinPar,iSample,iEvent);
}

double samplePDFDUNEBeamND::ReturnKinematicParameter(KinematicTypes KinPar, int iSample, int iEvent) {
  // HH: Special cases for when we are dealing with ints
  switch(KinPar){
    case kIsCC:
      return static_cast<double>(dunendmcSamples[iSample].rw_isCC[iEvent]);
    case kRecoNumu:
      return static_cast<double>(dunendmcSamples[iSample].rw_reco_numu[iEvent]);
    case kRecoNue:
      return static_cast<double>(dunendmcSamples[iSample].rw_reco_nue[iEvent]);
    case kNuPDG:
      return static_cast<double>(dunendmcSamples[iSample].rw_nuPDG[iEvent]);
    case kCCNumu:
      return static_cast<double>(dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].rw_nuPDG[iEvent])==14);
    case kCCNue:
      return static_cast<double>(dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].rw_nuPDG[iEvent])==12);
    case kNotCCNumu:
      // HH: Adapted from TDR definition
      // For details see:
      // https://github.com/DUNE/lblpwgtools/blob/3d475f50a998fbfa6266df9a0c4eb3056c0cdfe5/CAFAna/Systs/EnergySysts.h#L39
      // Not (CC Numu)
      return static_cast<double>(!(dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].rw_nuPDG[iEvent])==14));
    // HH: Otherwise use the old function
    default:
      return *GetPointerToKinematicParameter(KinPar, iSample, iEvent);
  }
}

void samplePDFDUNEBeamND::setupFDMC(int iSample) {
  dunemc_base *duneobj = &(dunendmcSamples[iSample]);
  FarDetectorCoreInfo *fdobj = &(MCSamples[iSample]);

  fdobj->SampleDetID = SampleDetID;
  for(int iEvent = 0 ;iEvent < fdobj->nEvents ; ++iEvent){
    fdobj->rw_etru[iEvent] = &(duneobj->rw_etru[iEvent]);
    fdobj->mode[iEvent] = &(duneobj->mode[iEvent]);
    fdobj->Target[iEvent] = &(duneobj->Target[iEvent]); 
    fdobj->isNC[iEvent] = !(duneobj->rw_isCC[iEvent]);
    fdobj->nupdgUnosc[iEvent] = &(duneobj->nupdgUnosc[iEvent]);
    fdobj->nupdg[iEvent] = &(duneobj->nupdg[iEvent]);
  }
}

/*
void samplePDFDUNEBeamND::applyShifts(int iSample, int iEvent) {
  // reset erec back to original value
  
  dunendmcSamples[iSample].rw_erec_shifted[iEvent] = dunendmcSamples[iSample].rw_erec[iEvent];

  //Calculate values needed
  double sqrtErecHad =  sqrt(dunendmcSamples[iSample].rw_erec_had[iEvent]);
  double sqrtErecLep =  sqrt(dunendmcSamples[iSample].rw_erec_lep[iEvent]);
  double sqrteRecoPi0 = sqrt(dunendmcSamples[iSample].rw_eRecoPi0[iEvent]);
  double sqrteRecoN = sqrt(dunendmcSamples[iSample].rw_eRecoN[iEvent]);
  double sumEhad = dunendmcSamples[iSample].rw_eRecoP[iEvent] + dunendmcSamples[iSample].rw_eRecoPip[iEvent] + dunendmcSamples[iSample].rw_eRecoPim[iEvent];
  double sqrtSumEhad = sqrt(sumEhad);

  double invSqrtErecHad =  1/(sqrtErecHad+0.1);
  double invSqrtErecLep =  1/(sqrtErecLep+0.1);
  double invSqrteRecoPi0 =  1/(sqrteRecoPi0+0.1);
  double invSqrteRecoN =  1/(sqrteRecoN+0.1);
  double invSqrtSumEhad =  1/(sqrtSumEhad+0.1);

  bool CCnumu {dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].rw_nuPDG[iEvent])==NuPDG::kNumu && dunendmcSamples[iSample].rw_nuPDGunosc[iEvent]==NuPDG::kNumu};
  bool CCnue {dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].rw_nuPDG[iEvent])==NuPDG::kNue && dunendmcSamples[iSample].rw_nuPDGunosc[iEvent]==NuPDG::kNue};
  bool NotCCnumu {!(dunendmcSamples[iSample].rw_isCC[iEvent]==1 && abs(dunendmcSamples[iSample].rw_nuPDG[iEvent])==14) && dunendmcSamples[iSample].rw_nuPDGunosc[iEvent]==NuPDG::kNumu};

/*
  TotalEScaleND(NDDetectorSystPointers[0], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_erec_had[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], NotCCnumu);

  TotalEScaleSqrtND(NDDetectorSystPointers[1], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_erec_had[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], sqrtErecHad, sqrtErecLep, NotCCnumu);

  TotalEScaleInvSqrtND(NDDetectorSystPointers[2], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_erec_had[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecHad, invSqrtErecLep, NotCCnumu);

  HadEScaleND(NDDetectorSystPointers[3], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], sumEhad);

  HadEScaleSqrtND(NDDetectorSystPointers[4], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], sumEhad, sqrtSumEhad);

  HadEScaleInvSqrtND(NDDetectorSystPointers[5], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], sumEhad, invSqrtSumEhad);

  MuEScaleND(NDDetectorSystPointers[6], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], CCnumu);

  MuEScaleSqrtND(NDDetectorSystPointers[7], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], sqrtErecLep, CCnumu);

  MuEScaleInvSqrtND(NDDetectorSystPointers[8], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecLep, CCnumu);

  NEScaleND(NDDetectorSystPointers[9], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoN[iEvent]);

  NEScaleSqrtND(NDDetectorSystPointers[10], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoN[iEvent], sqrteRecoN);

  NEScaleInvSqrtND(NDDetectorSystPointers[11], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoN[iEvent], invSqrteRecoN);

  EMEScaleND(NDDetectorSystPointers[12], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoPi0[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], CCnue);

  EMEScaleSqrtND(NDDetectorSystPointers[13], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoPi0[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], sqrtErecLep, sqrteRecoPi0, CCnue);

  EMEScaleInvSqrtND(NDDetectorSystPointers[14], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoPi0[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], invSqrtErecLep, invSqrteRecoPi0, CCnue);

  HadResND(NDDetectorSystPointers[15], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoP[iEvent], dunendmcSamples[iSample].rw_eRecoPip[iEvent], dunendmcSamples[iSample].rw_eRecoPim[iEvent], dunendmcSamples[iSample].rw_eP[iEvent], dunendmcSamples[iSample].rw_ePip[iEvent], dunendmcSamples[iSample].rw_ePim[iEvent]);

  MuResND(NDDetectorSystPointers[16], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], dunendmcSamples[iSample].rw_LepE[iEvent], CCnumu);

  NResND(NDDetectorSystPointers[17], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoN[iEvent], dunendmcSamples[iSample].rw_eN[iEvent]);

  EMResND(NDDetectorSystPointers[18], &dunendmcSamples[iSample].rw_erec_shifted[iEvent], dunendmcSamples[iSample].rw_eRecoPi0[iEvent], dunendmcSamples[iSample].rw_ePi0[iEvent], dunendmcSamples[iSample].rw_erec_lep[iEvent], dunendmcSamples[iSample].rw_LepE[iEvent], CCnue);
}
*/

std::vector<double> samplePDFDUNEBeamND::ReturnKinematicParameterBinning(std::string KinematicParameterStr) 
{
  std::vector<double> binningVector;
  return binningVector;
}

int samplePDFDUNEBeamND::ReturnKinematicParameterFromString(std::string KinematicParameterStr) {
  if(KinematicParameterStr == "TrueNeutrinoEnergy") return kTrueNeutrinoEnergy;
  if(KinematicParameterStr == "RecoQ") return kRecoQ;
  if(KinematicParameterStr == "RecoNeutrinoEnergy") return kRecoNeutrinoEnergy;
  if(KinematicParameterStr == "IsFHC") return kIsFHC;
  if(KinematicParameterStr == "RecoY") return kRecoY;
  if(KinematicParameterStr == "IsCC") return kIsCC;
  if(KinematicParameterStr == "RecoNumu") return kRecoNumu;
  if(KinematicParameterStr == "RecoNue") return kRecoNue;
  if(KinematicParameterStr == "NuPDG") return kNuPDG;
  if(KinematicParameterStr == "NotCCNumu") return kNotCCNumu;
  if(KinematicParameterStr == "CCNumu") return kCCNumu;
  if(KinematicParameterStr == "CCNue") return kCCNue;
  return -1;
}

std::string samplePDFDUNEBeamND::ReturnStringFromKinematicParameter(int KinematicParameter) {
  switch(KinematicParameter){
    case kTrueNeutrinoEnergy: return "TrueNeutrinoEnergy";
    case kRecoQ: return "RecoQ";
    case kRecoNeutrinoEnergy: return "RecoNeutrinoEnergy";
    case kIsFHC: return "IsFHC";
    case kRecoY: return "RecoY";
    case kIsCC: return "IsCC";
    case kRecoNumu: return "RecoNumu";
    case kRecoNue: return "RecoNue";
    case kNuPDG: return "NuPDG";
    case kNotCCNumu: return "NotCCnumu";
    default: return "";
  }
}
