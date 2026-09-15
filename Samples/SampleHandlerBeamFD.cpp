#include "SampleHandlerBeamFD.h"

SampleHandlerBeamFD::SampleHandlerBeamFD(std::string mc_version_, ParameterHandlerGeneric* ParHandler_,const std::shared_ptr<OscillationHandler>&  Oscillator_) : SampleHandlerBase(mc_version_, ParHandler_, Oscillator_) {
  KinematicParameters = &KinematicParametersDUNE;
  ReversedKinematicParameters = &ReversedKinematicParametersDUNE;

  Initialise();

  downsamplingStep = 1;
}

SampleHandlerBeamFD::~SampleHandlerBeamFD() {
}

void SampleHandlerBeamFD::Init() {
  beamFDSampleDetails.resize(GetNSamples());
  // dunemcSamples.resize(nSamples,dunemc_beamfd());

  auto EnabledSamples = Get<std::vector<std::string>>(SampleManager->raw()["Samples"], __FILE__ , __LINE__);

  for(int i = 0; i < GetNSamples(); i++){
    const auto TempTitle = EnabledSamples[i];
    beamFDSampleDetails[i].isFHC = Get<double>(SampleManager->raw()[TempTitle]["DUNESampleBools"]["isFHC"], __FILE__ , __LINE__);
    beamFDSampleDetails[i].iselike = Get<bool>(SampleManager->raw()[TempTitle]["DUNESampleBools"]["iselike"], __FILE__ , __LINE__);
    beamFDSampleDetails[i].pot = Get<double>(SampleManager->raw()[TempTitle]["POT"], __FILE__ , __LINE__);

    MACH3LOG_INFO("Setting up beam sample {}", GetSampleTitle(i));
    MACH3LOG_INFO("- isFHC: {}", beamFDSampleDetails[i].isFHC);
    MACH3LOG_INFO("- iselike: {}", beamFDSampleDetails[i].iselike);
  }

  downsamplingStep = GetFromManager<unsigned int>(SampleManager->raw()["DownsamplingStep"], 1, __FILE__ , __LINE__);
  if (downsamplingStep == 0) {
    throw MaCh3Exception(__FILE__, __LINE__,
      "Downsampling step cannot be zero. Please set it to a positive integer in the Beam FD sample config file."
    );
  }
  MACH3LOG_INFO("Beam FD downsampling step: {}", downsamplingStep);

  MACH3LOG_INFO("-------------------------------------------------------------------");
}

// ************************************************
void SampleHandlerBeamFD::InititialiseData()
{
  // ************************************************
  // Reweight MC to match
  Reweight();
  // set asimov data
  for (int iSample = 0; iSample < GetNSamples(); iSample++)
  {
    AddData(iSample, GetMCArray(iSample));
  }
}


void SampleHandlerBeamFD::SetupSplines() {

  ///@todo move all of the spline setup into core
  if(ParHandler->GetNumParamsFromSampleName(SampleHandlerName, kSpline) > 0){
    MACH3LOG_INFO("Found {} splines for this sample so I will create a spline object", ParHandler->GetNumParamsFromSampleName(SampleHandlerName, kSpline));
    SplineHandler = std::unique_ptr<BinnedSplineHandler>(new BinnedSplineHandlerDUNE(ParHandler,Modes.get()));
    InitialiseSplineObject();
  }
  else{
    MACH3LOG_INFO("Found {} splines for this sample so I will not load or evaluate splines", ParHandler->GetNumParamsFromSampleName(SampleHandlerName, kSpline));
    SplineHandler = nullptr;
  }

  return;
}


// === HH: Functional parameters ===
void TotalEScale(double const & par_val, dunemc_beamfd &ev) {
  // Total energy scale uncertainties for anything but CC Numu, see:
  // https://github.com/DUNE/lblpwgtools/blob/3d475f50a998fbfa6266df9a0c4eb3056c0cdfe5/CAFAna/Systs/EnergySysts.h#L39

  ev.rw_erec_shifted += par_val * ev.rw_erec_had;
}

void TotalEScaleNotCCNumu(double const & par_val, dunemc_beamfd &ev) {
  // A special case for Not (CC Numu), where we also scale Erec by lepton energy
  // Since we reconstruct muon energy in a different way, see:
  // https://github.com/DUNE/lblpwgtools/blob/3d475f50a998fbfa6266df9a0c4eb3056c0cdfe5/CAFAna/Systs/EnergySysts.h#L39
  ev.rw_erec_shifted += par_val * ev.rw_erec_lep;
}

void TotalEScaleSqrt(double const & par_val, dunemc_beamfd &ev) {
  ev.rw_erec_shifted += par_val * ev.rw_erec_had * ev.rw_erec_had_sqrt;
}

void TotalEScaleSqrtNotCCNumu(double const & par_val, dunemc_beamfd &ev) {
  // See comments in TotalEScaleNotCCNumu
  ev.rw_erec_shifted += par_val * ev.rw_erec_lep * ev.rw_erec_lep_sqrt;
}

void TotalEScaleInvSqrt(double const & par_val, dunemc_beamfd &ev) {
  // Erec/sqrt(Erec) = sqrt(Erec)
  ev.rw_erec_shifted += par_val * ev.rw_erec_had_sqrt;
}

void TotalEScaleInvSqrtNotCCNumu(double const & par_val, dunemc_beamfd &ev) {
  // See comments in TotalEScaleNotCCNumu
  ev.rw_erec_shifted += par_val * ev.rw_erec_lep_sqrt;
}

void HadEScale(double const & par_val, dunemc_beamfd &ev) {
  ev.rw_erec_shifted += par_val * ev.rw_sum_ehad;
}

void HadEScaleSqrt(double const & par_val, dunemc_beamfd &ev) {
  ev.rw_erec_shifted += par_val * ev.rw_sum_ehad * ev.rw_sum_ehad_sqrt;
}

void HadEScaleInvSqrt(double const & par_val, dunemc_beamfd &ev) {
  ev.rw_erec_shifted += par_val * ev.rw_sum_ehad_sqrt;
}

void MuEScale(double const & par_val, dunemc_beamfd &ev) {
  // HH TODO: Functionally this is the same as TotalEScaleNotCCNumu, not sure if this function is even needed
  TotalEScaleNotCCNumu(par_val, ev);
}

void MuEScaleSqrt(double const & par_val, dunemc_beamfd &ev) {
  // See comments in MuEScale
  TotalEScaleSqrtNotCCNumu(par_val, ev);
}

void MuEScaleInvSqrt(double const & par_val, dunemc_beamfd &ev) {
  // See comments in MuEScale
  TotalEScaleInvSqrtNotCCNumu(par_val, ev);
}

void NEScale(double const & par_val, dunemc_beamfd &ev) {
  ev.rw_erec_shifted += par_val * ev.rw_eRecoN;
}

void NEScaleSqrt(double const & par_val, dunemc_beamfd &ev) {
  ev.rw_erec_shifted += par_val * ev.rw_eRecoN * ev.rw_eRecoN_sqrt;
}

void NEScaleInvSqrt(double const & par_val, dunemc_beamfd &ev) {
  ev.rw_erec_shifted += par_val * ev.rw_eRecoN_sqrt;
}

void EMEScale(double const & par_val, dunemc_beamfd &ev) {
  ev.rw_erec_shifted += par_val * ev.rw_eRecoPi0;
}

void EMEScaleCCNue(double const & par_val, dunemc_beamfd &ev) {
  // Again this is the same as TotalEScaleNotCCNumu, not sure if this function is needed
  TotalEScaleNotCCNumu(par_val, ev);
}

void EMEScaleSqrt(double const & par_val, dunemc_beamfd &ev) {
  ev.rw_erec_shifted += par_val * ev.rw_eRecoPi0 * ev.rw_eRecoPi0_sqrt;
}

void EMEScaleSqrtCCNue(double const & par_val, dunemc_beamfd &ev) {
  // See comments in EMEScaleCCNue
  TotalEScaleSqrtNotCCNumu(par_val, ev);
}

void EMEScaleInvSqrt(double const & par_val, dunemc_beamfd &ev) {
  ev.rw_erec_shifted += par_val * ev.rw_eRecoPi0_sqrt;
}

void EMEScaleInvSqrtCCNue(double const & par_val, dunemc_beamfd &ev) {
  // See comments in EMEScaleCCNue
  TotalEScaleInvSqrtNotCCNumu(par_val, ev);
}

void HadRes(double const & par_val, dunemc_beamfd &ev) {
  // True sum - reco sum
  ev.rw_erec_shifted += par_val * (ev.rw_eP
    + ev.rw_ePip
    + ev.rw_ePim
    - ev.rw_sum_ehad);
}

void MuRes(double const & par_val, dunemc_beamfd &ev) {
  // True muon energy - reco muon energy
  ev.rw_erec_shifted += par_val * (ev.rw_LepE - ev.rw_erec_lep);
}

void NRes(double const & par_val, dunemc_beamfd &ev) {
  // True neutron energy - reco neutron energy
  ev.rw_erec_shifted += par_val * (ev.rw_eN - ev.rw_eRecoN);
}

void EMRes(double const & par_val, dunemc_beamfd &ev) {
  // True pi0 energy - reco pi0 energy
  ev.rw_erec_shifted += par_val * (ev.rw_ePi0 - ev.rw_eRecoPi0);
}

void EMResCCNue(double const & par_val, dunemc_beamfd &ev) {
  // This is the same as MuRes, again not sure if this function is needed
  MuRes(par_val, ev);
}

void RecoCVNNumu(double const & par_val, dunemc_beamfd &ev) {
  // CVN numu uncertainty
  ev.rw_cvnnumu_shifted += par_val;
}

void RecoCVNNue(double const & par_val, dunemc_beamfd &ev) {
  // CVN nue uncertainty
  ev.rw_cvnnue_shifted += par_val;
}

void SampleHandlerBeamFD::RegisterFunctionalParameters() {
  MACH3LOG_INFO("Registering functional parameters");
  // This function manually populates the map of functional parameters
  // Maps the name of the functional parameter to the pointer of the function

  RegisterIndividualFunctionalParameter(dunemcSamples, "TotalEScaleFD",
                                        TotalEScale);

  RegisterIndividualFunctionalParameter(dunemcSamples, "TotalEScaleNotCCNumuFD",
                                        TotalEScaleNotCCNumu);

  RegisterIndividualFunctionalParameter(dunemcSamples, "TotalEScaleSqrtFD",
                                        TotalEScaleSqrt);

  RegisterIndividualFunctionalParameter(
      dunemcSamples, "TotalEScaleSqrtNotCCNumuFD", TotalEScaleSqrtNotCCNumu);

  RegisterIndividualFunctionalParameter(dunemcSamples, "TotalEScaleInvSqrtFD",
                                        TotalEScaleInvSqrt);

  RegisterIndividualFunctionalParameter(dunemcSamples,
                                        "TotalEScaleInvSqrtNotCCNumuFD",
                                        TotalEScaleInvSqrtNotCCNumu);

  RegisterIndividualFunctionalParameter(dunemcSamples, "HadEScaleFD",
                                        HadEScale);

  RegisterIndividualFunctionalParameter(dunemcSamples, "HadEScaleSqrtFD",
                                        HadEScaleSqrt);

  RegisterIndividualFunctionalParameter(dunemcSamples, "HadEScaleInvSqrtFD",
                                        HadEScaleInvSqrt);

  RegisterIndividualFunctionalParameter(dunemcSamples, "MuEScaleFD", MuEScale);

  RegisterIndividualFunctionalParameter(dunemcSamples, "MuEScaleSqrtFD",
                                        MuEScaleSqrt);

  RegisterIndividualFunctionalParameter(dunemcSamples, "MuEScaleInvSqrtFD",
                                        MuEScaleInvSqrt);

  RegisterIndividualFunctionalParameter(dunemcSamples, "NEScaleFD", NEScale);

  RegisterIndividualFunctionalParameter(dunemcSamples, "NEScaleSqrtFD",
                                        NEScaleSqrt);

  RegisterIndividualFunctionalParameter(dunemcSamples, "NEScaleInvSqrtFD",
                                        NEScaleInvSqrt);

  RegisterIndividualFunctionalParameter(dunemcSamples, "EMEScaleFD", EMEScale);

  RegisterIndividualFunctionalParameter(dunemcSamples, "EMEScaleCCNueFD",
                                        EMEScaleCCNue);

  RegisterIndividualFunctionalParameter(dunemcSamples, "EMEScaleSqrtFD",
                                        EMEScaleSqrt);

  RegisterIndividualFunctionalParameter(dunemcSamples, "EMEScaleSqrtCCNueFD",
                                        EMEScaleSqrtCCNue);

  RegisterIndividualFunctionalParameter(dunemcSamples, "EMEScaleInvSqrtFD",
                                        EMEScaleInvSqrt);

  RegisterIndividualFunctionalParameter(dunemcSamples, "EMEScaleInvSqrtCCNueFD",
                                        EMEScaleInvSqrtCCNue);

  RegisterIndividualFunctionalParameter(dunemcSamples, "HadResFD", HadRes);

  RegisterIndividualFunctionalParameter(dunemcSamples, "MuResFD", MuRes);

  RegisterIndividualFunctionalParameter(dunemcSamples, "NResFD", NRes);

  RegisterIndividualFunctionalParameter(dunemcSamples, "EMResFD", EMRes);

  RegisterIndividualFunctionalParameter(dunemcSamples, "EMResCCNueFD",
                                        EMResCCNue);

  RegisterIndividualFunctionalParameter(dunemcSamples, "RecoCVNNumuFD",
                                        RecoCVNNumu);

  RegisterIndividualFunctionalParameter(dunemcSamples, "RecoCVNNueFD",
                                        RecoCVNNue);

  MACH3LOG_INFO("Finished registering functional parameters");
}

// HH: Reset the shifted values to the original values
void SampleHandlerBeamFD::ResetShifts(int iEvent) {
  dunemcSamples[iEvent].rw_erec_shifted = dunemcSamples[iEvent].rw_erec;
  dunemcSamples[iEvent].rw_cvnnumu_shifted = dunemcSamples[iEvent].rw_cvnnumu;
  dunemcSamples[iEvent].rw_cvnnue_shifted = dunemcSamples[iEvent].rw_cvnnue;
}
// =================================

void SampleHandlerBeamFD::AddAdditionalWeightPointers() {
  for (size_t i = 0; i < dunemcSamples.size(); ++i) {
    MCEvents[i].total_weight_pointers.push_back(&(dunemcSamples[i].pot_s));
    MCEvents[i].total_weight_pointers.push_back( &(dunemcSamples[i].norm_s));
    MCEvents[i].total_weight_pointers.push_back( &(dunemcSamples[i].rw_berpaacvwgt));
    MCEvents[i].total_weight_pointers.push_back( &(dunemcSamples[i].flux_w));
  }
}


int SampleHandlerBeamFD::SetupExperimentMC() {

  // dunemc_base *duneobj = &(dunemcSamples[iSample]);

  MACH3LOG_INFO("-------------------------------------------------------------------");
  TChain* _data = new TChain("caf");
  // Maps the file index within the TChain (GetTreeNumber()) to its sample index and
  // per-file norm values.
  std::vector<size_t> fileIndexToSample;
  std::vector<std::array<double, 2>> fileIndexToNorm; // [norm_s, pot_s]
  for (size_t iSample=0; iSample<SampleDetails.size(); iSample++) {
    for (const std::vector<std::string>& files : SampleDetails[iSample].mc_files) {
      for (const std::string& filename : files){

        MACH3LOG_INFO("Adding file to TChain: {}", filename);

        TFile* _sampleFile = TFile::Open(filename.c_str(), "READ");
        // HH: still have the read the individual ROOT file to get the norm histograms
        TH1D* norm = _sampleFile->Get<TH1D>("norm");
        if(!norm){
          MACH3LOG_ERROR("Add a norm KEY to the root file using MakeNormHists.cxx");
          throw MaCh3Exception(__FILE__, __LINE__);
        }
        fileIndexToSample.push_back(iSample);
        fileIndexToNorm.push_back({norm->GetBinContent(1), beamFDSampleDetails[iSample].pot / norm->GetBinContent(2)});
        _sampleFile->Close();
        // HH: Check whether the file exists, see https://root.cern/doc/master/classTChain.html#a78a896924ac6c7d3691b7e013bcbfb1c
        int _add_rtn = _data->Add(filename.c_str(), -1);
        if(_add_rtn == 0){
          MACH3LOG_ERROR("Could not add file {} to TChain, please check the file exists and is readable", filename);
          throw MaCh3Exception(__FILE__, __LINE__);
        }
      }
    }
  }


  if(_data){
    MACH3LOG_INFO("Number of entries in TChain: {}", _data->GetEntries());
  }
  else{
    MACH3LOG_ERROR("Failed to create the TChain.");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  //Reco Variables
  double _erec;
  double _erec_nue;
  double _erec_had;
  double _erec_had_nue;
  double _erec_lep;
  double _erec_lep_nue;

  double _eRecoP;
  double _eRecoPip;
  double _eRecoPim;
  double _eRecoPi0;
  double _eRecoN;

  double _cvnnumu;
  double _cvnnue;
  double _vtx_x;
  double _vtx_y;
  double _vtx_z;

  //Truth Variables
  int _mode;
  double _ev;
  double _LepE;
  double _eP;
  double _ePip;
  double _ePim;
  double _ePi0;
  double _eN;
  double _BeRPA_cvwgt;
  int _isCC;
  int _nuPDGunosc;
  int _nuPDG;

  _data->SetBranchStatus("*", 0);
  _data->SetBranchStatus("Ev", 1);
  _data->SetBranchAddress("Ev", &_ev);
  _data->SetBranchStatus("Ev_reco_numu", 1);
  _data->SetBranchAddress("Ev_reco_numu", &_erec);
  _data->SetBranchStatus("Ev_reco_nue", 1);
  _data->SetBranchAddress("Ev_reco_nue", &_erec_nue);
  _data->SetBranchStatus("RecoHadEnNumu", 1);
  _data->SetBranchAddress("RecoHadEnNumu", &_erec_had);
  _data->SetBranchStatus("RecoHadEnNue", 1);
  _data->SetBranchAddress("RecoHadEnNue", &_erec_had_nue);
  _data->SetBranchStatus("RecoLepEnNumu", 1);
  _data->SetBranchAddress("RecoLepEnNumu", &_erec_lep);
  _data->SetBranchStatus("RecoLepEnNue", 1);
  _data->SetBranchAddress("RecoLepEnNue", &_erec_lep_nue);

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

  _data->SetBranchStatus("mode",1);
  _data->SetBranchAddress("mode",&_mode);
  _data->SetBranchStatus("cvnnumu",1);
  _data->SetBranchAddress("cvnnumu", &_cvnnumu);
  _data->SetBranchStatus("cvnnue",1);
  _data->SetBranchAddress("cvnnue", &_cvnnue);
  _data->SetBranchStatus("isCC", 1);
  _data->SetBranchAddress("isCC", &_isCC);
  _data->SetBranchStatus("nuPDGunosc", 1);
  _data->SetBranchAddress("nuPDGunosc", &_nuPDGunosc);
  _data->SetBranchStatus("nuPDG", 1);
  _data->SetBranchAddress("nuPDG", &_nuPDG);
  _data->SetBranchStatus("BeRPA_A_cvwgt", 1);
  _data->SetBranchAddress("BeRPA_A_cvwgt", &_BeRPA_cvwgt);
  _data->SetBranchStatus("vtx_x", 1);
  _data->SetBranchAddress("vtx_x", &_vtx_x);
  _data->SetBranchStatus("vtx_y", 1);
  _data->SetBranchAddress("vtx_y", &_vtx_y);
  _data->SetBranchStatus("vtx_z", 1);
  _data->SetBranchAddress("vtx_z", &_vtx_z);

  size_t nEntries = static_cast<size_t>(_data->GetEntries());
  size_t nDownsampledEntries = nEntries / downsamplingStep;
  if (nDownsampledEntries == 0) {
    throw MaCh3Exception(__FILE__, __LINE__,
      "Downsampling step is too large, resulting in zero entries. Please set it to a smaller positive integer in the Beam FD sample config file."
    );
  }
  size_t countwidth = nDownsampledEntries / 5;
  dunemcSamples.resize(nDownsampledEntries);
  _data->GetEntry(0);

  // HH: A map to keep track of negative energies
  std::unordered_map<std::string, int> negative_counts;
  // Initialize the negative counts for each energy variable
  negative_counts["rw_erec_had"] = 0;
  negative_counts["rw_erec_lep"] = 0;
  negative_counts["rw_eRecoN"] = 0;
  negative_counts["rw_eRecoPi0"] = 0;
  negative_counts["rw_sum_ehad"] = 0;

  //FILL DUNE STRUCT
  for (unsigned int i = 0; i < nDownsampledEntries; ++i) { // Loop through tree
    _data->GetEntry(i * downsamplingStep);

    if (i % countwidth == 0) {
      M3::Utils::PrintProgressBar(i, static_cast<Long64_t>(nDownsampledEntries));
    }

    const size_t sample_index = fileIndexToSample[static_cast<size_t>(_data->GetTreeNumber())];
    const bool iselike_temp = beamFDSampleDetails[sample_index].iselike;
    dunemcSamples[i].SampleIndex = static_cast<int>(sample_index);
    dunemcSamples[i].nupdgUnosc = _nuPDGunosc;
    dunemcSamples[i].nupdg = _nuPDG;
    dunemcSamples[i].OscChannelIndex = static_cast<double>(GetOscChannel(SampleDetails[sample_index].OscChannels, dunemcSamples[i].nupdgUnosc, dunemcSamples[i].nupdg));

    // POT stuff
    dunemcSamples[i].norm_s = fileIndexToNorm[static_cast<size_t>(_data->GetTreeNumber())][0]; // Norm in sample
    dunemcSamples[i].pot_s = fileIndexToNorm[static_cast<size_t>(_data->GetTreeNumber())][1] * downsamplingStep; // POT in sample

    dunemcSamples[i].rw_cvnnumu = (_cvnnumu);
    dunemcSamples[i].rw_cvnnue = (_cvnnue);
    dunemcSamples[i].rw_cvnnumu_shifted = (_cvnnumu);
    dunemcSamples[i].rw_cvnnue_shifted = (_cvnnue);
    if (iselike_temp) {
      dunemcSamples[i].rw_erec = (_erec_nue);
      dunemcSamples[i].rw_erec_shifted = (_erec_nue);
      dunemcSamples[i].rw_erec_had = (_erec_had_nue);
      dunemcSamples[i].rw_erec_lep = (_erec_lep_nue);
    } else {
      dunemcSamples[i].rw_erec = (_erec);
      dunemcSamples[i].rw_erec_shifted = (_erec);
      dunemcSamples[i].rw_erec_had = (_erec_had);
      dunemcSamples[i].rw_erec_lep = (_erec_lep);
    }

    dunemcSamples[i].rw_eRecoP = (_eRecoP);
    dunemcSamples[i].rw_eRecoPip = (_eRecoPip);
    dunemcSamples[i].rw_eRecoPim = (_eRecoPim);
    dunemcSamples[i].rw_eRecoPi0 = (_eRecoPi0);
    dunemcSamples[i].rw_eRecoN = (_eRecoN);
    dunemcSamples[i].rw_LepE = (_LepE);
    dunemcSamples[i].rw_eP = (_eP);
    dunemcSamples[i].rw_ePip = (_ePip);
    dunemcSamples[i].rw_ePim = (_ePim);
    dunemcSamples[i].rw_ePi0 = (_ePi0);
    dunemcSamples[i].rw_eN = (_eN);

    // HH: Add checks to make sure the energies are not negative
    if (dunemcSamples[i].rw_erec_had < 0) {
      dunemcSamples[i].rw_erec_had = 0;
      negative_counts["rw_erec_had"]++;
    }
    if (dunemcSamples[i].rw_erec_lep < 0) {
      dunemcSamples[i].rw_erec_lep = 0;
      negative_counts["rw_erec_lep"]++;
    }
    if (dunemcSamples[i].rw_eRecoN < 0) {
      dunemcSamples[i].rw_eRecoN = 0;
      negative_counts["rw_eRecoN"]++;
    }
    if (dunemcSamples[i].rw_eRecoPi0 < 0) {
      dunemcSamples[i].rw_eRecoPi0 = 0;
      negative_counts["rw_eRecoPi0"]++;
    }

    dunemcSamples[i].rw_erec_had_sqrt = sqrt(dunemcSamples[i].rw_erec_had);
    dunemcSamples[i].rw_erec_lep_sqrt = sqrt(dunemcSamples[i].rw_erec_lep);
    dunemcSamples[i].rw_eRecoN_sqrt = sqrt(dunemcSamples[i].rw_eRecoN);
    dunemcSamples[i].rw_eRecoPi0_sqrt = sqrt(dunemcSamples[i].rw_eRecoPi0);

    dunemcSamples[i].rw_sum_ehad = dunemcSamples[i].rw_eRecoP + dunemcSamples[i].rw_eRecoPip + dunemcSamples[i].rw_eRecoPim;
    if (dunemcSamples[i].rw_sum_ehad < 0) {
      dunemcSamples[i].rw_sum_ehad = 0;
      negative_counts["rw_sum_ehad"]++;
    }
    dunemcSamples[i].rw_sum_ehad_sqrt = sqrt(dunemcSamples[i].rw_sum_ehad);

    dunemcSamples[i].enu_true = (_ev);
    dunemcSamples[i].rw_isCC = _isCC;
    // dunemcSamples[i].rw_nuPDGunosc = _nuPDGunosc;
    // dunemcSamples[i].rw_nuPDG = _nuPDG;
    dunemcSamples[i].rw_berpaacvwgt = (_BeRPA_cvwgt);
    dunemcSamples[i].rw_vtx_x = (_vtx_x);
    dunemcSamples[i].rw_vtx_y = (_vtx_y);
    dunemcSamples[i].rw_vtx_z = (_vtx_z);

    dunemcSamples[i].rw_trueccnumu = static_cast<double>(dunemcSamples[i].rw_isCC==1 && abs(dunemcSamples[i].nupdg)==14);
    dunemcSamples[i].rw_trueccnue = static_cast<double>(dunemcSamples[i].rw_isCC==1 && abs(dunemcSamples[i].nupdg)==12);

    //Assume everything is on Argon40 for now....
    dunemcSamples[i].Target = kTarget_Ar;

    int M3Mode = Modes->GetModeFromGenerator(std::abs(_mode));
    if (!_isCC) M3Mode += 14; //Account for no ability to distinguish CC/NC
    if (M3Mode > 15) M3Mode -= 1; //Account for no NCSingleKaon
    dunemcSamples[i].mode = M3Mode;

    dunemcSamples[i].flux_w = 1.0;
  }

  // HH: Give a warning if any negative energies were found
  for (const auto& pair : negative_counts) {
    if (pair.second > 0) {
      MACH3LOG_WARN("Found {} negative values for {}.", pair.second, pair.first);
    }
  }

  delete _data;
  return static_cast<int>(nDownsampledEntries);
}

const double* SampleHandlerBeamFD::GetPointerToKinematicParameter(const int KinPar, const int iEvent) const {
  switch(KinPar){
  case kTrueNeutrinoEnergy:
    return &(dunemcSamples[iEvent].enu_true);
  case kRecoNeutrinoEnergy:
    return &(dunemcSamples[iEvent].rw_erec_shifted);
    break;
  case kTrueXPos:
    return &(dunemcSamples[iEvent].rw_vtx_x);
  case kTrueYPos:
    return &(dunemcSamples[iEvent].rw_vtx_y);
  case kTrueZPos:
    return &(dunemcSamples[iEvent].rw_vtx_z);
  case kCVNNumu:
    return &(dunemcSamples[iEvent].rw_cvnnumu_shifted);
  case kCVNNue:
    return &(dunemcSamples[iEvent].rw_cvnnue_shifted);
  case kM3Mode:
    return &(dunemcSamples[iEvent].mode);
  case kOscChannel:
    return &(dunemcSamples[iEvent].OscChannelIndex);
  case kIsFHC:
    return &(beamFDSampleDetails[MCEvents[iEvent].NominalSample].isFHC);
  case kTrueCCnue:
	return &(dunemcSamples[iEvent].rw_trueccnue);
  case kTrueCCnumu:
	  return &(dunemcSamples[iEvent].rw_trueccnumu);
  case kTargetNucleus:
    return &(dunemcSamples[iEvent].Target);
  default:
    MACH3LOG_ERROR("Did not recognise Kinematic Parameter type {}...", KinPar);
    throw MaCh3Exception(__FILE__, __LINE__);
  }
}


double SampleHandlerBeamFD::ReturnKinematicParameter(const int KinematicVariable, const int iEvent) const{
  KinematicTypes KinPar = static_cast<KinematicTypes>(KinematicVariable);
  return *GetPointerToKinematicParameter(KinPar, iEvent);
}

void SampleHandlerBeamFD::SetupMC() {
  // dunemc_base *duneobj = &(dunemcSamples[iSample]);
  // FarDetectorCoreInfo *fdobj = &(MCEvents[iSample]);

  for (unsigned int iEvent = 0; iEvent < GetNEvents(); ++iEvent) {
    MCEvents[iEvent].enu_true = dunemcSamples[iEvent].enu_true;
    MCEvents[iEvent].isNC = !(dunemcSamples[iEvent].rw_isCC);
    MCEvents[iEvent].nupdg = dunemcSamples[iEvent].nupdg;
    MCEvents[iEvent].nupdgUnosc = dunemcSamples[iEvent].nupdgUnosc;
    MCEvents[iEvent].NominalSample = dunemcSamples[iEvent].SampleIndex;
  }

}

