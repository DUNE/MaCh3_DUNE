#pragma once

// Include the input manager
#include "Manager/Manager.h"

// Include the SampleHandlers
#include "Samples/SampleHandlerFD.h"

/// @brief Factory function that generates MaCh3 DUNE instance including configured samples
/// @param fitMan Configuration Manager 
/// @param sample_vec Vector of SampleHandler objects
/// @param xsec Cross-section covariance matrix
/// @param osc Oscillation covariance matrix
void MakeMaCh3DuneInstance(std::unique_ptr<manager>& fitMan, std::vector<SampleHandlerFD*> &sample_vec,  ParameterHandlerGeneric *&xsec);
