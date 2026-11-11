#pragma once

// Std includes
#include "memory"
#include <string>
#include <utility>
#include <vector>

// Include the input manager
#include "Manager/Manager.h"

// Include the SampleHandlers
#include "Samples/SampleHandlerBase.h"
#include "Fitters/Algorithms/MaCh3Factory.h"

// BeamNDCov; SampleHandlerBase.h above pulls in OscillationHandler
#include "Samples/StructsDUNE.h"

/// @brief Construct a single DUNE sample handler of the requested type
/// @param SampleType One of "BeamFD", "BeamND", "BeamNDGAr", "Atm", "PDSP"
/// @param SampleConfig Path to the sample's YAML config
/// @param param_handler Parameter handler; ownership stays with the caller
/// @param BeamOscillator_ Oscillation handler for beam samples, may be null
/// @param AtmOscillator_ Oscillation handler for atmospheric samples, may be null
/// @param beamNDCov ND covariance matrices, only required for "BeamND"
/// @returns Newly allocated sample handler, owned by the caller
SampleHandlerBase* GetMaCh3DuneInstance(std::string SampleType,
                                        std::string SampleConfig,
                                        std::unique_ptr<ParameterHandlerGeneric>& param_handler,
                                        const std::shared_ptr<OscillationHandler>& BeamOscillator_,
                                        const std::shared_ptr<OscillationHandler>& AtmOscillator_,
                                        BeamNDCov beamNDCov);

/// @brief Build every sample handler listed in the configuration
/// @param FitManager Configuration Manager
/// @param param_handler Parameter handler the samples are attached to
/// @returns Vector of sample handlers, each owned by the caller
std::vector<SampleHandlerBase*> MaCh3DuneSampleFactory(std::unique_ptr<Manager>& FitManager,
                                                       std::unique_ptr<ParameterHandlerGeneric>& param_handler);

/// @brief Was real data configured, i.e. has MaCh3DuneFactory already set every sample's data?
/// @param FitManager Configuration Manager
/// @returns true when General.DataFile is set
/// @details Apps which would otherwise make Asimov data out of the nominal prediction have to ask
/// this first. AddData overwrites whatever data a sample already holds, so calling it
/// unconditionally silently throws away the histograms MaCh3DuneFactory loaded from the data file
/// and fits the nominal MC instead
bool UsingDataFromFile(std::unique_ptr<Manager>& FitManager);

/// @brief Factory function that generates MaCh3 DUNE instance including configured samples
/// @param fitMan Configuration Manager
/// @returns parameter handler and vector of sample handler
/// @note When General.DataFile is set, each sample's data is loaded from it here. See
/// UsingDataFromFile() before adding Asimov data on top
std::pair<std::unique_ptr<ParameterHandlerGeneric>, std::vector<SampleHandlerBase*>> MaCh3DuneFactory(std::unique_ptr<Manager> &FitManager);
