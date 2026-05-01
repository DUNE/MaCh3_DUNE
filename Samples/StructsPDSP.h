#pragma once

/// @brief class holding basic information about MC
struct PDSPMCInfo {
  /// @brief Constructor
  PDSPMCInfo() {
  }
  /// @brief Destructor
  ~PDSPMCInfo() {
  }
  /// interaction mode, might relate to neutrino interaction mode, so also leave this unused?
  double Mode = M3::_BAD_INT_;
  /// Apparently is a MaCh3 required variable (but will be unused.)
  double OscillationChannel = M3::_BAD_INT_;

  /// True Interacting Kinetic Energy
  double TrueKEInt = M3::_BAD_DOUBLE_;
  /// Reconstructed Interacting Kinetic Energy
  double RecoKEInt = M3::_BAD_DOUBLE_;
  /// Reconstructed track length (cm), from track_length_reco
  double RecoEndZ = M3::_BAD_DOUBLE_;

};

struct MetaData {
  /// not entirely clear what this needs to be (file index, MC sample number, process type etc.)
  int SampleIndex = M3::_BAD_INT_;

};

struct PDSPMCPlottingInfo {
  /// True Interacting Kinetic Energy
  double TrueKEInt = M3::_BAD_DOUBLE_;
  /// Reconstructed Interacting Kinetic Energy
  double RecoKEInt = M3::_BAD_DOUBLE_;
};
