#include <iostream>
#include <chrono>
#include <iomanip>
#include <vector>

#include <TFile.h>

#include "Samples/MaCh3DUNEFactory.h"
#include "Samples/StructsDUNE.h"
#include "Fitters/MaCh3Factory.h"

using MaCh3Inst = std::pair<std::unique_ptr<ParameterHandlerGeneric>, std::vector<SampleHandlerBase*>>;

//###############################################################################################################################
// HW: DO NOT CHANGE ME!
MaCh3Inst InitMaCh3(int argc, char * argv[]){
//###############################################################################################################################
  auto FitManager = MaCh3ManagerFactory(argc, argv);

  //###############################################################################################################################
  //Create sample handler + parameter_handler objects
  auto [param_handler, samples] = MaCh3DuneFactory(FitManager);

  //###############################################################################################################################
  //Perform reweight, print total integral, and set data

  std::vector<std::unique_ptr<TH1>> DUNEHists;
  for(auto handler : samples){
    for (unsigned iSample = 0; iSample < handler->GetNSamples(); ++iSample) {
      handler->Reweight();
      DUNEHists.push_back(M3::Clone(handler->GetMCHist(iSample)));
      MACH3LOG_INFO("Event rate for {} : {:<5.2f}", handler->GetSampleTitle(iSample), handler->GetMCHist(iSample)->Integral());

      if (handler->GetNDim(iSample) == 1) {
        handler->AddData(iSample, static_cast<TH1D*>(DUNEHists.back().get()));
      } else if (handler->GetNDim(iSample) == 2) {
        handler->AddData(iSample, static_cast<TH2D*>(DUNEHists.back().get()));
      }
    }
  }

  return {std::move(param_handler), samples};
}

// Helper Function
std::vector<double> Simulate(MaCh3Inst& mach3_handler, std::vector<double> param_values){
    /// Runs a single simulation, returns x_obs
    mach3_handler.first->SetParameters(param_values);

    /// We can check they're in bounds
    if(mach3_handler.first->GetLikelihood()>=M3::_LARGE_LOGL_){
        return {};
    }

    std::vector<double> mc_vec;
    for(auto s : mach3_handler.second){
        s->Reweight();
        auto mc = s->GetMCArray();
        mc_vec.insert(mc_vec.end(), mc.begin(), mc.end());
    }
    return mc_vec;
}

void FillTree(TTree* tree, std::vector<std::string>& branch_names, const std::vector<std::vector<double>>& table){

    std::vector<double> values(branch_names.size());

    for(size_t i =0; i<branch_names.size(); ++i){
        tree->Branch(branch_names[i].c_str(), &values[i]);
    }

    for(const auto& row : table){
        std::copy(row.begin(), row.end(), values.begin());
        tree->Fill();
    }
    tree->Write();

}

void SaveToOutput(MaCh3Inst& mach3_handler, const std::string& outfile, const std::vector<std::vector<double>>& theta, const std::vector<std::vector<double>>& x_obs){

    // Open TFile
    TFile* out = TFile::Open(outfile.c_str(), "RECREATE");
    out->cd();

    if(theta.size()!=x_obs.size()){
        MACH3LOG_ERROR("Provided {} parameter values != {} observations", theta.size(), x_obs.size());
        throw MaCh3Exception(__FILE__, __LINE__);
    }

    // Create TTrees
    auto theta_tree = new TTree("theta", "theta");
    
    // Now we need some parameter names
    int n_pars = mach3_handler.first->GetNumParams();
    auto theta_pars_names = std::vector<std::string>(n_pars);
    
    
    // Let make the theta names
    for(int i=0; i<n_pars; ++i){
        theta_pars_names[i] = mach3_handler.first->GetParName(i);
    }
    
    FillTree(theta_tree, theta_pars_names, theta);
    
    auto obs_tree = new TTree("x_obs", "x_obs");
    // For xobs we can make something a bit smarter!
    int n_mc_bins = x_obs[0].size();
    auto x_obs_par_names  = std::vector<std::string>(n_mc_bins);
    for(int i=0; i<n_mc_bins; i++){
        x_obs_par_names[i] = "x_obs_" + std::to_string(i);
    }
    FillTree(obs_tree, x_obs_par_names, x_obs);

    out->Close();

}

std::vector<double> Linspace(double start, double end, int n_bins){
    double increment = (end-start)/n_bins;

    std::vector<double> ls(n_bins);
    for(int i=0; i<n_bins; i++){
        ls[i]=increment*i;
    }
    return ls;

}

//// =========================================
//// THIS IS THE METHOD YOU'LL NEED TO EDIT!!!!
//// =========================================
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> CreateObservations(MaCh3Inst& mach3_handler){
    /// THIS IS THE CODE YOU NEED TO GENERATE!
    auto& par_handler = mach3_handler.first;

    // Let's try doing this with delta_cp
    int delta_index = par_handler->GetParIndex("delta_cp");

    // We can now get lower/upper bounds with
    double lower_bound = par_handler->GetLowerBound(delta_index);
    double upper_bound = par_handler->GetUpperBound(delta_index);
    
    int n_sims = 100;

    /// We can now generate the range
    std::vector<double> dcp_vals = Linspace(lower_bound, upper_bound, n_sims);

    /// Now we make the vectors we want to fill
    std::vector<std::vector<double>> theta(n_sims);
    std::vector<std::vector<double>> x_obs(n_sims);

    /// We can also get the nominal
    auto nominal = par_handler->GetPreFitValues();

    for(int i=0; i<n_sims; i++){
        nominal[delta_index] = dcp_vals[i];
        /// Store the parameter value
        theta[i] = nominal;
        x_obs[i] = Simulate(mach3_handler, nominal);
    }

    return {theta, x_obs};
}

int main(int argc, char * argv[]) {
    // Firstly we initialise MaCh3
    auto mach3_handlers = InitMaCh3(argc, argv);

    /// YOUR METHOD GOES HERE
    auto [theta, x_obs] = CreateObservations(mach3_handlers);

    std::string output_file = "observations.root";
    SaveToOutput(mach3_handlers, output_file, theta, x_obs);
}
