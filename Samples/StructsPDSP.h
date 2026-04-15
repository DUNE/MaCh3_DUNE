#pragma once

/// @brief class holding basic information about tutorial MC
struct PDSPMCInfo {
  /// @brief Constructor
  PDSPMCInfo() {
  }
  /// @brief Destructor
  ~PDSPMCInfo() {
  }
  /// interaction mode
  double Mode = M3::_BAD_INT_;
  /// True Interacting Kinetic Energy
  double TrueKEInt = M3::_BAD_DOUBLE_;
  /// Reconstructed Interacting Kinetic Energy
  double RecoKEInt = M3::_BAD_DOUBLE_;

};

struct PDSPMCPlottingInfo {
  /// True Interacting Kinetic Energy
  double TrueKEInt = M3::_BAD_DOUBLE_;
  /// Reconstructed Interacting Kinetic Energy
  double RecoKEInt = M3::_BAD_DOUBLE_;
};
