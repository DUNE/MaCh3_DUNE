add_library(FDTDRFakeDataStudies SHARED WeightToNuWro.cxx)

target_include_directories(FDTDRFakeDataStudies PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

install(TARGETS FDTDRFakeDataStudies
  EXPORT mach3dune-targets
  LIBRARY DESTINATION lib/)
