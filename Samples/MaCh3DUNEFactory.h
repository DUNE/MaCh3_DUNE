#pragma once

// Std includes
#include "memory"

// Include the input manager
#include "Manager/Manager.h"

// Include the SampleHandlers
#include "Samples/SampleHandlerBase.h"
#include "Fitters/MaCh3Factory.h"

// DUNE Handlers
#include "Samples/SampleHandlerBeamFD.h"
#include "Samples/SampleHandlerBeamND.h"
#include "Samples/SampleHandlerBeamNDGAr.h"
#include "Samples/SampleHandlerAtm.h"

/// @brief Factory function that generates MaCh3 DUNE instance including configured samples
/// @param fitMan Configuration Manager
/// @param param_handler Cross-section covariance matrix
std::vector<SampleHandlerBase*> MaCh3DuneSampleFactory(std::unique_ptr<Manager> &FitManager, std::unique_ptr<ParameterHandlerGeneric> &param_handler);
