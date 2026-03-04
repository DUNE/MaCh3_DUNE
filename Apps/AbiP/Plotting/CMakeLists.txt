cmake_minimum_required(VERSION 3.14)

# Collect all source files

set(PLOTTING_SOURCES
  Compare_onaxis_vs_offaxis.cpp
  Compare_runplans.cpp
  PosteriorPredictive_forFD.cpp
  PosteriorPredictive_forND.cpp
  Post_xsec_diagnostics_withcorr.cpp
  Post_xsec_2D.cpp
)

# Add executables

foreach(src ${PLOTTING_SOURCES})
  get_filename_component(exeName ${src} NAME_WE)
  add_executable(${exeName} ${src})

  # Include directories
  target_include_directories(${exeName} PRIVATE
      ${CMAKE_SOURCE_DIR}/Samples
      ${CMAKE_SOURCE_DIR}/Apps/AbiP/Plotting
      ${YAML_CPP_INCLUDE_DIR}   # If yaml-cpp is found via find_package
  )

  # Link libraries
  target_link_libraries(${exeName} PRIVATE
      MaCh3DUNE::All
  )
endforeach()
