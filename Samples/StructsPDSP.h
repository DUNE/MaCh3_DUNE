#pragma once

#include "Samples/SampleStructs.h"

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
  /// ProtoDUNE-SP interactions are on argon.
  double Target = kTarget_Ar;

  /// True Initial Kinetic Energy
  double TrueKEIni = M3::_BAD_DOUBLE_;
  /// True Interacting Kinetic Energy
  double TrueKEInt = M3::_BAD_DOUBLE_;
  /// True track end psoition, from track_length_reco
  double TrueEndZ = M3::_BAD_DOUBLE_;

  /// Reconstructed kinetic energy at the TPC front face.
  double RecoKEFF = M3::_BAD_DOUBLE_;
  double RecoKEFFShifted = M3::_BAD_DOUBLE_;
  /// Reconstructed kinetic energy at the start of the fiducial volume.
  double RecoKEIni = M3::_BAD_DOUBLE_;
  /// Reconstructed Interacting Kinetic Energy
  double RecoKEInt = M3::_BAD_DOUBLE_;
  /// Reconstructed track length (cm), from track_length_reco
  double RecoEndZ = M3::_BAD_DOUBLE_;

  /// Instrumental reconstructed beam momentum (MeV/c).
  double RecoPinst = M3::_BAD_DOUBLE_;
  double RecoPinstShifted = M3::_BAD_DOUBLE_;
  /// Reconstructed trajectory length (cm), distinct from the endpoint z.
  double RecoTrackLength = M3::_BAD_DOUBLE_;
  double RecoTrackLengthShifted = M3::_BAD_DOUBLE_;
  /// Fixed fractional N(0, 0.026) track-length smear drawn once per event.
  double TrackLengthSmearFraction = 0.0;
  /// Shifted fit observables used by functional systematics.
  double RecoKEIniShifted = M3::_BAD_DOUBLE_;
  double RecoKEIntShifted = M3::_BAD_DOUBLE_;
  /// Multiplicative beam-spectrum correction.
  M3::float_t BeamMomentumWeight = 1.0;

};

struct MetaData {
  /// not entirely clear what this needs to be (file index, MC sample number, process type etc.)
  int SampleIndex = M3::_BAD_INT_;

};

struct PDSPMCPlottingInfo {
  /// True Interacting Kinetic Energy
  double TrueKEIni = M3::_BAD_DOUBLE_;
  /// True Interacting Kinetic Energy
  double TrueKEInt = M3::_BAD_DOUBLE_;
  /// Reconstructed Interacting Kinetic Energy
  double RecoKEIni = M3::_BAD_DOUBLE_;
  /// Reconstructed Interacting Kinetic Energy
  double RecoKEInt = M3::_BAD_DOUBLE_;
};
