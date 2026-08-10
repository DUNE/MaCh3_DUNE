set(SOURCE_FILES
    ReadEvents.cpp
    SampleHandlerBeamFDStandardRecord.cpp
)
set(HEADERS
  EventInfo.h
  Projections.h
  ReadEvents.h
  SampleHandlerBeamFDStandardRecord.h
)

add_library(SampleHandlerBeamFDStandardRecord SHARED
  ${SOURCE_FILES}
)

set_target_properties(SampleHandlerBeamFDStandardRecord PROPERTIES
  PUBLIC_HEADER "${HEADERS}"
  EXPORT_NAME SampleHandlerBeamFDStandardRecord)

target_include_directories(SampleHandlerBeamFDStandardRecord PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS SampleHandlerBeamFDStandardRecord
  EXPORT mach3dune-targets
  LIBRARY DESTINATION lib/
  PUBLIC_HEADER DESTINATION include/Samples)

target_link_libraries(SampleHandlerBeamFDStandardRecord PUBLIC SplinesDUNE SystematicsDUNE MaCh3::All MaCh3DUNECompilerOptions)
target_link_libraries(SampleHandlerBeamFDStandardRecord PRIVATE duneanaobj::StandardRecord duneanasel::all DUNEMaCh3Warnings ROOT::TreePlayer)

