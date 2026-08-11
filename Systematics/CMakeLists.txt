set(HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/Beam.h
)
set(SOURCE_FILES
  Beam.cpp
)

add_library(SystematicsDUNE SHARED
  ${SOURCE_FILES}
)

set_target_properties(SystematicsDUNE PROPERTIES
  PUBLIC_HEADER "${HEADERS}"
  EXPORT_NAME SystematicsDUNE)

target_link_libraries(SystematicsDUNE PUBLIC MaCh3DUNECompilerOptions)
target_link_libraries(SystematicsDUNE PRIVATE duneanafluxtools::all DUNEMaCh3Warnings)

target_include_directories(SystematicsDUNE PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS SystematicsDUNE
  EXPORT mach3dune-targets
  LIBRARY DESTINATION lib/
  PUBLIC_HEADER DESTINATION include/Samples)

add_library(MACH3DUNE::SystematicsDUNE ALIAS SystematicsDUNE)

add_executable(DumpBeamSystematicRatios DumpBeamSystematicRatios.cpp)
target_link_libraries(DumpBeamSystematicRatios SystematicsDUNE ROOT::Hist)

install(TARGETS DumpBeamSystematicRatios DESTINATION bin)
