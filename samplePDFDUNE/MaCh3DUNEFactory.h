#pragma once
#include "samplePDFDUNE/samplePDFDUNEBeamFD_actuallyFD.h"
#ifdef BUILD_NDGAR
#include "samplePDFDUNE/samplePDFDUNEBeamNDGAr.h"
#else
#include "samplePDFDUNE/samplePDFDUNEBeamFD.h"
#include "samplePDFDUNE/samplePDFDUNEBeamND.h"
#include "samplePDFDUNE/samplePDFDUNEAtm.h"
#endif

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
