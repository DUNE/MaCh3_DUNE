set(SOURCE_FILES
    ReadEvents.cpp
    Systematics.cpp
    Utility.cpp
)
set(HEADERS
  EventInfo.h
  Projections.h
  ReadEvents.h
  Systematics.h
  Utility.h
)

add_library(SamplesDUNEBeamOffAxis SHARED
  ${SOURCE_FILES}
)

set_target_properties(SamplesDUNEBeamOffAxis PROPERTIES
  PUBLIC_HEADER "${HEADERS}"
  EXPORT_NAME SamplesDUNEBeamOffAxis)

target_include_directories(SamplesDUNEBeamOffAxis PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS SamplesDUNEBeamOffAxis
  EXPORT mach3dune-targets
  LIBRARY DESTINATION lib/
  PUBLIC_HEADER DESTINATION include/Samples)

target_link_libraries(SamplesDUNEBeamOffAxis PUBLIC SplinesDUNE SystematicsDUNE duneanaobj::all MaCh3::All MaCh3DUNECompilerOptions)
target_link_libraries(SamplesDUNEBeamOffAxis PRIVATE FDTDRFakeDataStudies DUNEMaCh3Warnings ROOT::TreePlayer)
