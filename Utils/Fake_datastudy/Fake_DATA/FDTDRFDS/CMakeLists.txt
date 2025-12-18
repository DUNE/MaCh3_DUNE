add_library(FakeDataStudies SHARED WeightToNuWro.cxx)

target_include_directories(FakeDataStudies PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

install(TARGETS FakeDataStudies
  EXPORT mach3dune-targets
  LIBRARY DESTINATION lib/)