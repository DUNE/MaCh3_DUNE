#pragma once

// Include the input manager
#include "manager/manager.h"

// Include the samplePDFs
#include "samplePDF/samplePDFFDBase.h"

/// @brief Factory function that generates MaCh3 DUNE instance including configured samples
/// @param fitMan Configuration Manager 
/// @param sample_vec Vector of samplePDF objects
/// @param xsec Cross-section covariance matrix
/// @param osc Oscillation covariance matrix
void MakeMaCh3DuneInstance(manager *fitMan, std::vector<samplePDFFDBase*> &sample_vec,  covarianceXsec *&xsec, covarianceOsc *&osc);

/// @brief Gets MaCh3 DUNE samplePDF instance
/// @param SampleType Sample type (BeamFD, BeamND, Atm)
/// @param SampleConfig Configuration file 
/// @param xsec Cross-section covariance matrix
/// @return samplePDF instance
samplePDFFDBase* GetMaCh3DuneInstance(std::string SampleType, std::string SampleConfig, covarianceXsec* &xsec);
